"""MD pose-stability check on the top leads.

For each lead ligand:
  1. Build receptor+ligand complex from the docked best pose
  2. Parameterize with amber14 + SMIRNOFF (ligand) + OBC2 implicit solvent
  3. Minimize → 100 ps equilibration → production MD (default 5 ns)
  4. Track ligand heavy-atom RMSD-from-start and center-of-mass drift
  5. Verdict: ligand stable in pocket?  PASS if RMSD < 3 Å and COM drift < 4 Å

Output: per-ligand DCD + CSV + plot, plus a stability.csv summary.

Run:
    PATH=$CONDA_PREFIX/bin:$PATH python md_stability.py
    PATH=$CONDA_PREFIX/bin:$PATH python md_stability.py --ns 10 --ligands EGFRm-009 EGFRm-002
"""
from __future__ import annotations
import argparse
import csv
import sys
import time
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT))

from lib.refine import _ligand_from_smiles_and_pose   # noqa: E402

MUT_RECEPTOR = ROOT / "projects" / "egfr_demo" / "receptor" / "protein_prepared.pdb"
ROUNDS_DIR   = HERE.parent / "rounds"
FINAL_UNIQUE = HERE.parent / "final_ranked_unique.csv"
OUT_DIR      = HERE.parent / "_md_stability"
SUMMARY_CSV  = HERE.parent / "stability.csv"


def find_pose_sdf(name: str) -> Path | None:
    for m in sorted(ROUNDS_DIR.glob(f"round_*/docking/poses/{name}_poses.sdf"),
                    reverse=True):
        return m
    return None


def first_record(src_sdf: Path, dst_sdf: Path) -> Path:
    text = src_sdf.read_text()
    end = text.find("$$$$")
    dst_sdf.write_text(text[:end+4] + "\n" if end >= 0 else text)
    return dst_sdf


def load_lead_rows(names: list[str]) -> list[dict]:
    """Resolve ligand names → CSV rows. final_ranked_unique.csv is dedup'd by
    canonical SMILES but still contains name collisions across rounds (Claude
    re-used names). The CSV is rank-sorted, so the FIRST occurrence per name
    is the best-scored one — that's the one we want."""
    rows = list(csv.DictReader(FINAL_UNIQUE.open()))
    by_name: dict[str, dict] = {}
    for r in rows:
        if r["name"] not in by_name:
            by_name[r["name"]] = r
    return [by_name[n] for n in names if n in by_name]


def run_one_ligand(name: str, smiles: str, work: Path, ns: float,
                   report_ps: float = 50.0, platform_name: str = "OpenCL") -> dict:
    """Run implicit-solvent MD on the complex, return stability metrics."""
    from openmm.app import (PDBFile, ForceField, Modeller, NoCutoff, HBonds,
                            Simulation, DCDReporter, StateDataReporter)
    from openmm import LangevinMiddleIntegrator, Platform, unit
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    import numpy as np

    work.mkdir(parents=True, exist_ok=True)
    src = find_pose_sdf(name)
    if not src:
        raise FileNotFoundError(f"no pose SDF for {name}")
    pose_sdf = first_record(src, work / "pose.sdf")

    # 1. Ligand: bonds from SMILES, coords from docked pose
    lig = _ligand_from_smiles_and_pose(smiles, pose_sdf)
    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
    gen = SMIRNOFFTemplateGenerator(molecules=lig, forcefield="openff-2.1.0")
    ff.registerTemplateGenerator(gen.generator)

    # 2. Complex topology
    rec = PDBFile(str(MUT_RECEPTOR))
    lig_top = lig.to_topology().to_openmm()
    lig_pos = lig.conformers[0].to_openmm()
    modeller = Modeller(rec.topology, rec.positions)
    modeller.add(lig_top, lig_pos)
    cpx_top = modeller.topology
    n_rec = len(list(rec.topology.atoms()))
    n_tot = len(list(cpx_top.atoms()))
    lig_atom_idx = list(range(n_rec, n_tot))

    # Heavy-atom indices for ligand RMSD
    lig_heavy = [n_rec + i for i, atom in enumerate(lig_top.atoms())
                 if atom.element.symbol != "H"]

    # 3. System with implicit solvent, no PME, no barostat
    system = ff.createSystem(cpx_top, nonbondedMethod=NoCutoff,
                             constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin,
                                          1.0 / unit.picosecond,
                                          2.0 * unit.femtoseconds)
    try:
        plat = Platform.getPlatformByName(platform_name)
        sim = Simulation(cpx_top, system, integrator, plat)
    except Exception:
        sim = Simulation(cpx_top, system, integrator)
    sim.context.setPositions(modeller.positions)

    print(f"  [{name}] minimizing …", flush=True)
    sim.minimizeEnergy(maxIterations=500)

    eq_steps = int(100 * 1000 / 2)   # 100 ps at 2 fs
    print(f"  [{name}] equilibrating 100 ps ({eq_steps} steps) …", flush=True)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
    sim.step(eq_steps)

    # Capture the START-of-production structure as the RMSD reference
    ref_pdb = work / "start.pdb"
    with open(ref_pdb, "w") as fh:
        PDBFile.writeFile(sim.topology,
                          sim.context.getState(getPositions=True).getPositions(),
                          fh)

    # 4. Production MD
    prod_steps = int(ns * 1000 * 1000 / 2)
    report_steps = int(report_ps * 1000 / 2)
    n_frames = prod_steps // report_steps

    traj_dcd = work / "traj.dcd"
    log_csv  = work / "log.csv"
    sim.reporters.append(DCDReporter(str(traj_dcd), report_steps))
    sim.reporters.append(StateDataReporter(str(log_csv), report_steps,
                                            step=True, time=True,
                                            potentialEnergy=True,
                                            temperature=True, speed=True))

    print(f"  [{name}] production {ns} ns ({prod_steps} steps, {n_frames} frames) …",
          flush=True)
    t0 = time.time()
    sim.step(prod_steps)
    wall = time.time() - t0
    print(f"  [{name}] MD wall time: {wall:.0f}s ({prod_steps/wall:.0f} steps/s)",
          flush=True)

    # 5. Analyze stability with MDAnalysis
    import MDAnalysis as mda
    from MDAnalysis.analysis import align as mda_align

    # Load trajectory with the topology
    u = mda.Universe(str(ref_pdb), str(traj_dcd))
    ref = mda.Universe(str(ref_pdb))

    # Align trajectory to receptor backbone, then compute ligand RMSD
    prot_sel = "protein and name CA"
    lig_sel  = f"bynum {lig_heavy[0]+1}:{lig_heavy[-1]+1} and not name H*"

    aligner = mda_align.AlignTraj(u, ref, select=prot_sel, in_memory=True).run()

    lig_rmsd = []
    com_drift = []
    ref_lig_atoms = ref.select_atoms(f"bynum {lig_heavy[0]+1}:{lig_heavy[-1]+1}")
    ref_com = ref_lig_atoms.center_of_mass()
    ref_lig_pos = ref_lig_atoms.positions.copy()

    u_lig = u.select_atoms(f"bynum {lig_heavy[0]+1}:{lig_heavy[-1]+1}")
    for ts in u.trajectory:
        cur = u_lig.positions
        # heavy-atom RMSD after alignment to receptor
        diff = cur - ref_lig_pos
        lig_rmsd.append(float(np.sqrt((diff*diff).sum() / len(cur))))
        com_drift.append(float(np.linalg.norm(u_lig.center_of_mass() - ref_com)))

    lig_rmsd = np.array(lig_rmsd)
    com_drift = np.array(com_drift)
    times_ps = np.arange(len(lig_rmsd)) * report_ps

    rmsd_csv = work / "rmsd.csv"
    with rmsd_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_ps", "lig_rmsd_A", "com_drift_A"])
        for t, r, c in zip(times_ps, lig_rmsd, com_drift):
            w.writerow([f"{t:.1f}", f"{r:.3f}", f"{c:.3f}"])

    return {
        "name": name,
        "ns": ns,
        "n_frames": len(lig_rmsd),
        "rmsd_mean": float(lig_rmsd.mean()),
        "rmsd_final": float(lig_rmsd[-1]),
        "rmsd_max": float(lig_rmsd.max()),
        "com_drift_final": float(com_drift[-1]),
        "com_drift_max": float(com_drift.max()),
        "stable": bool(lig_rmsd.max() < 3.0 and com_drift.max() < 4.0),
        "wall_s": round(wall, 1),
        "rmsd_csv": str(rmsd_csv),
        "traj_dcd": str(traj_dcd),
    }


def plot_stability(results: list[dict], out_png: Path):
    """Plot ligand RMSD vs time for all ligands on one chart."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    for r in results:
        if not Path(r["rmsd_csv"]).exists():
            continue
        rows = list(csv.DictReader(open(r["rmsd_csv"])))
        t = [float(x["time_ps"]) for x in rows]
        rmsd = [float(x["lig_rmsd_A"]) for x in rows]
        com  = [float(x["com_drift_A"]) for x in rows]
        ax1.plot(t, rmsd, label=r["name"], linewidth=1.8)
        ax2.plot(t, com,  label=r["name"], linewidth=1.8)
    for ax, ylabel, threshold in [(ax1, "Ligand heavy-atom RMSD (Å)", 3.0),
                                  (ax2, "COM drift from start (Å)",   4.0)]:
        ax.axhline(threshold, color="red", linestyle="--", alpha=0.4, label=f"pass < {threshold} Å")
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel.split(" (")[0])
        ax.legend(loc="upper left", fontsize=9)
        ax.grid(alpha=0.3)
    fig.suptitle("MD stability of EGFR T790M/L858R docked poses (implicit solvent, OBC2)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    print(f"  Plot → {out_png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligands", nargs="+",
                    default=["EGFRm-009", "EGFRm-002", "EGFR-TM-014"],
                    help="Names from final_ranked_unique.csv (default: top 2 leads + 1 negative control)")
    ap.add_argument("--ns", type=float, default=5.0,
                    help="Production MD length per ligand (ns)")
    ap.add_argument("--report-ps", type=float, default=50.0,
                    help="Frame interval (ps)")
    ap.add_argument("--platform", default="OpenCL")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    leads = load_lead_rows(args.ligands)
    if not leads:
        sys.exit("None of the requested ligand names were found in final_ranked_unique.csv")

    print(f"MD stability check ({args.ns} ns implicit-OBC2, platform={args.platform})")
    print(f"  Ligands: {', '.join(r['name'] for r in leads)}")
    print(f"  Output : {OUT_DIR}\n")

    results = []
    for r in leads:
        print(f"=== {r['name']}  (Vina mut={r['score']}) ===")
        try:
            out = run_one_ligand(r["name"], r["smiles"],
                                 OUT_DIR / r["name"], args.ns,
                                 report_ps=args.report_ps,
                                 platform_name=args.platform)
            results.append(out)
            print(f"  RMSD mean={out['rmsd_mean']:.2f} Å  "
                  f"final={out['rmsd_final']:.2f}  max={out['rmsd_max']:.2f}")
            print(f"  COM  final={out['com_drift_final']:.2f} Å  "
                  f"max={out['com_drift_max']:.2f}")
            verdict = "✓ STABLE" if out["stable"] else "✗ UNSTABLE"
            print(f"  Verdict: {verdict}\n")
        except Exception as e:
            print(f"  ! FAILED: {e}\n")
            import traceback; traceback.print_exc()
            results.append({"name": r["name"], "error": str(e), "stable": False})

    # Summary CSV
    fields = ["name","ns","n_frames","rmsd_mean","rmsd_final","rmsd_max",
              "com_drift_final","com_drift_max","stable","wall_s"]
    with SUMMARY_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    print(f"\n✓ stability.csv → {SUMMARY_CSV}")

    plot_stability(results, OUT_DIR / "rmsd_plot.png")

    # Final summary
    print("\n" + "=" * 70)
    print(f"{'Ligand':<14} {'rmsd_mean':>10} {'rmsd_max':>10} {'com_max':>10}  verdict")
    print("-" * 70)
    for r in results:
        if "error" in r:
            print(f"{r['name']:<14}  FAILED: {r['error']}")
            continue
        v = "✓ STABLE" if r["stable"] else "✗ UNSTABLE"
        print(f"{r['name']:<14} {r['rmsd_mean']:>10.2f} {r['rmsd_max']:>10.2f} "
              f"{r['com_drift_max']:>10.2f}  {v}")
    print("=" * 70)


if __name__ == "__main__":
    main()
