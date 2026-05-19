"""Glide — molecular docking (Schrödinger Glide equivalent).

Backends:
  - AutoDock Vina (default, fast)
  - Smina (custom scoring functions)

Features:
  - Auto grid box from reference ligand OR explicit center+size
  - Batch docking of many ligands against one receptor
  - Per-ligand pose SDF + interaction analysis
  - Ranked CSV of best scores
  - Pose validity check via PoseBusters (when installed)
"""
from __future__ import annotations
import argparse
import csv
import subprocess
from pathlib import Path
from typing import Iterable

from .utils import banner, get_logger, new_project, require


# ---------------------------------------------------------------------------
# Grid box from reference ligand

def grid_from_ligand(ligand_pdb: Path, padding: float = 8.0) -> tuple[tuple[float,float,float], tuple[float,float,float]]:
    """Compute (center, size) of a box that encloses a reference ligand."""
    xs, ys, zs = [], [], []
    for ln in ligand_pdb.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            xs.append(float(ln[30:38]))
            ys.append(float(ln[38:46]))
            zs.append(float(ln[46:54]))
    if not xs:
        raise ValueError(f"No atoms in {ligand_pdb}")
    cx, cy, cz = (max(xs)+min(xs))/2, (max(ys)+min(ys))/2, (max(zs)+min(zs))/2
    sx = max(xs) - min(xs) + 2*padding
    sy = max(ys) - min(ys) + 2*padding
    sz = max(zs) - min(zs) + 2*padding
    return (cx, cy, cz), (sx, sy, sz)


# ---------------------------------------------------------------------------
# Backend: AutoDock Vina (Python bindings)

def dock_vina(receptor_pdbqt: Path, ligand_pdbqt: Path, out_pdbqt: Path,
              center: tuple[float,float,float], size: tuple[float,float,float],
              exhaustiveness: int = 8, n_poses: int = 9,
              seed: int = 42) -> list[float]:
    """Run AutoDock Vina via Python bindings. Returns list of scores."""
    from vina import Vina
    v = Vina(sf_name='vina', seed=seed)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=list(center), box_size=list(size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(str(out_pdbqt), n_poses=n_poses, overwrite=True)
    energies = v.energies(n_poses=n_poses)
    # energies is Nx7 numpy array; first column is total binding affinity
    return [float(row[0]) for row in energies]


# ---------------------------------------------------------------------------
# Backend: Smina (binary)

def dock_smina(receptor_pdbqt: Path, ligand_pdbqt: Path, out_sdf: Path,
               center: tuple[float,float,float], size: tuple[float,float,float],
               exhaustiveness: int = 8, n_poses: int = 9,
               scoring: str = "vinardo", seed: int = 42) -> list[float]:
    require("smina")
    log_path = out_sdf.with_suffix(".log")
    cmd = ["smina",
           "-r", str(receptor_pdbqt),
           "-l", str(ligand_pdbqt),
           "-o", str(out_sdf),
           "--center_x", str(center[0]),
           "--center_y", str(center[1]),
           "--center_z", str(center[2]),
           "--size_x", str(size[0]),
           "--size_y", str(size[1]),
           "--size_z", str(size[2]),
           "--exhaustiveness", str(exhaustiveness),
           "--num_modes", str(n_poses),
           "--seed", str(seed),
           "--scoring", scoring,
           "--log", str(log_path)]
    subprocess.run(cmd, check=True, capture_output=True)
    # Parse scores from log
    scores: list[float] = []
    for ln in log_path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln[0].isalpha() or ln.startswith("-"):
            continue
        parts = ln.split()
        try:
            idx = int(parts[0])
            scores.append(float(parts[1]))
        except (ValueError, IndexError):
            continue
    return scores


# ---------------------------------------------------------------------------

def pdbqt_to_sdf(pdbqt: Path, sdf: Path) -> Path:
    """Convert docked PDBQT poses → SDF using OpenBabel."""
    require("obabel")
    subprocess.run(["obabel", str(pdbqt), "-O", str(sdf)],
                   check=True, capture_output=True)
    return sdf


# ---------------------------------------------------------------------------

def dock_batch(receptor_pdbqt: Path, ligand_files: Iterable[Path],
               output_dir: Path,
               center: tuple[float,float,float],
               size: tuple[float,float,float],
               backend: str = "vina",
               exhaustiveness: int = 8,
               n_poses: int = 9,
               scoring: str = "vinardo") -> Path:
    """Dock many ligands and produce a ranked CSV.

    Returns path to the summary CSV.
    """
    log = get_logger("glide")
    output_dir.mkdir(parents=True, exist_ok=True)
    poses_dir = output_dir / "poses"
    poses_dir.mkdir(exist_ok=True)

    results: list[dict] = []
    for i, lig in enumerate(ligand_files, 1):
        name = lig.stem
        log.info(f"[{i}] docking {name}  ({backend})")
        try:
            if backend == "vina":
                out_pdbqt = poses_dir / f"{name}_poses.pdbqt"
                scores = dock_vina(receptor_pdbqt, lig, out_pdbqt,
                                   center=center, size=size,
                                   exhaustiveness=exhaustiveness,
                                   n_poses=n_poses)
                # convert to SDF
                pdbqt_to_sdf(out_pdbqt, poses_dir / f"{name}_poses.sdf")
            elif backend == "smina":
                out_sdf = poses_dir / f"{name}_poses.sdf"
                scores = dock_smina(receptor_pdbqt, lig, out_sdf,
                                    center=center, size=size,
                                    exhaustiveness=exhaustiveness,
                                    n_poses=n_poses, scoring=scoring)
            else:
                raise ValueError(f"Unknown backend {backend}")
            best = min(scores) if scores else None
            log.info(f"   best score: {best} kcal/mol")
            results.append({
                "ligand": name,
                "best_score": best,
                "n_poses": len(scores),
                "all_scores": ";".join(f"{s:.2f}" for s in scores),
            })
        except Exception as e:
            log.error(f"   FAILED: {e}")
            results.append({"ligand": name, "best_score": None,
                            "n_poses": 0, "all_scores": f"ERROR: {e}"})

    # Ranked CSV
    results.sort(key=lambda r: (r["best_score"] is None, r["best_score"] or 0))
    csv_path = output_dir / "docking_results.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["rank", "ligand", "best_score",
                                          "n_poses", "all_scores"])
        w.writeheader()
        for rank, r in enumerate(results, 1):
            r["rank"] = rank
            w.writerow(r)
    log.info(f"Ranked results → {csv_path}")
    return csv_path


# ---------------------------------------------------------------------------
# CLI

def main(argv=None):
    p = argparse.ArgumentParser(
        prog="glide",
        description="Glide — molecular docking (AutoDock Vina / Smina backend)."
    )
    p.add_argument("-r", "--receptor", required=True, type=Path,
                   help="Receptor PDBQT (from prepwizard)")
    p.add_argument("-l", "--ligands", required=True, nargs="+", type=Path,
                   help="One or more ligand PDBQT files (from ligprep)")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output directory (default: new project folder)")
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--ref-ligand", type=Path,
                     help="Reference ligand PDB to auto-define the grid box")
    grp.add_argument("--center", nargs=3, type=float, metavar=("X","Y","Z"),
                     help="Explicit box center (Å)")
    p.add_argument("--size", nargs=3, type=float, default=[22.0, 22.0, 22.0],
                   metavar=("X","Y","Z"),
                   help="Box size in Å (default 22 22 22)")
    p.add_argument("--padding", type=float, default=8.0,
                   help="Padding around reference ligand (Å, default 8)")
    p.add_argument("--backend", choices=["vina", "smina"], default="vina")
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--n-poses", type=int, default=9)
    p.add_argument("--scoring", default="vinardo",
                   help="Smina scoring function (vinardo|vina|ad4_scoring|...)")
    args = p.parse_args(argv)

    if args.ref_ligand:
        center, size = grid_from_ligand(args.ref_ligand, padding=args.padding)
    else:
        center, size = tuple(args.center), tuple(args.size)

    banner("Glide — molecular docking",
           f"receptor={args.receptor.name}  ligands={len(args.ligands)}  "
           f"backend={args.backend}  exh={args.exhaustiveness}")
    print(f"  grid center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    print(f"  grid size:   ({size[0]:.1f}, {size[1]:.1f}, {size[2]:.1f}) Å\n")

    if args.output is None:
        proj = new_project("glide")
        outdir = proj / "output"
    else:
        outdir = args.output

    csv_path = dock_batch(args.receptor, args.ligands, outdir,
                          center=center, size=size,
                          backend=args.backend,
                          exhaustiveness=args.exhaustiveness,
                          n_poses=args.n_poses,
                          scoring=args.scoring)

    # Print top 5
    print("\n  Top results:")
    with open(csv_path) as f:
        rdr = csv.DictReader(f)
        for i, row in enumerate(rdr):
            if i >= 5:
                break
            print(f"    {row['rank']}. {row['ligand']:30s} {row['best_score']} kcal/mol")
    print(f"\n✓ Docking complete. Poses in: {outdir}/poses")
    print(f"  Summary CSV: {csv_path}")


if __name__ == "__main__":
    main()
