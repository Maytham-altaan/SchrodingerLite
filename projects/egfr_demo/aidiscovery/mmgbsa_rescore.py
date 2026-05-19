"""MM-GBSA rescoring of the top EGFR hits in both mutant and wild-type.

For each of the top-N unique hits:
  1. Take the best Vina pose against the mutant receptor (4ZAU)
  2. Run MM-GBSA-style ΔG_bind via amber14 + OBC2 + SMIRNOFF on the ligand
  3. Repeat against the wild-type receptor (1M17) using the WT-docked pose
  4. ΔΔG_mmgbsa = wt_dG − mut_dG;  positive = mutant-selective
  5. Compare to Vina ΔΔG to validate or refute the selectivity ordering

Output: mmgbsa.csv with per-ligand mut_dG, wt_dG, ΔΔG, plus an HTML report.

Run (env's PATH must include obabel/vina from schrodinger-lite):
    PATH=$CONDA_PREFIX/bin:$PATH python mmgbsa_rescore.py --top 5
"""
from __future__ import annotations
import argparse
import csv
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT))

from lib.refine import mmgbsa_score_v2   # noqa: E402

MUT_RECEPTOR    = ROOT / "projects" / "egfr_demo" / "receptor" / "protein_prepared.pdb"
WT_RECEPTOR     = HERE.parent / "_wt_receptor" / "protein_prepared.pdb"

ROUNDS_DIR      = HERE.parent / "rounds"
WT_POSES_DIR    = HERE.parent / "_selectivity" / "poses"

FINAL_UNIQUE    = HERE.parent / "final_ranked_unique.csv"
SELECTIVITY_CSV = HERE.parent / "selectivity.csv"
MMGBSA_CSV      = HERE.parent / "mmgbsa.csv"
WORKDIR         = HERE.parent / "_mmgbsa"


def find_mut_pose_sdf(name: str) -> Path | None:
    for m in sorted(ROUNDS_DIR.glob(f"round_*/docking/poses/{name}_poses.sdf"),
                    reverse=True):
        return m
    return None


def extract_first_record(src_sdf: Path, dst_sdf: Path) -> Path:
    """Pull MODEL 1 (first $$$$-terminated record) out of a Vina multi-pose SDF."""
    text = src_sdf.read_text()
    end = text.find("$$$$")
    if end < 0:
        dst_sdf.write_text(text)
    else:
        dst_sdf.write_text(text[: end + 4] + "\n")
    return dst_sdf


def load_top(csv_path: Path, n: int) -> list[dict]:
    rows = list(csv.DictReader(csv_path.open()))
    # Already deduplicated and ranked
    return rows[:n]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--iters", type=int, default=800,
                    help="Minimization iterations per component (default 800)")
    args = ap.parse_args()

    if not FINAL_UNIQUE.exists():
        sys.exit(f"missing {FINAL_UNIQUE} — run view_top_hits.py first")
    if not WT_RECEPTOR.exists():
        sys.exit(f"missing {WT_RECEPTOR} — run selectivity_check.py first")

    WORKDIR.mkdir(parents=True, exist_ok=True)
    top = load_top(FINAL_UNIQUE, args.top)
    print(f"MM-GBSA rescore of top {len(top)} unique hits ({args.iters} iters/component)")
    print(f"  mut receptor: {MUT_RECEPTOR.relative_to(ROOT)}")
    print(f"  wt  receptor: {WT_RECEPTOR.relative_to(ROOT)}\n")

    results: list[dict] = []
    for i, r in enumerate(top, 1):
        name  = r["name"]
        smi   = r["smiles"]
        vina_mut = float(r["score"])

        # Build SDFs of best mut & wt poses
        mut_src = find_mut_pose_sdf(name)
        wt_src  = WT_POSES_DIR / f"{name}.pdbqt"   # we get pose from selectivity
        wt_sdf_src = WT_POSES_DIR / f"{name}_poses.sdf"
        if not mut_src or not mut_src.exists():
            print(f"[{i}/{len(top)}] {name}  ! missing mut pose")
            continue
        if not wt_sdf_src.exists():
            print(f"[{i}/{len(top)}] {name}  ! missing wt pose")
            continue

        mut_pose = extract_first_record(mut_src,    WORKDIR / f"{name}_mut.sdf")
        wt_pose  = extract_first_record(wt_sdf_src, WORKDIR / f"{name}_wt.sdf")

        print(f"[{i}/{len(top)}] {name}   Vina mut={vina_mut:+.2f}")

        t0 = time.time()
        try:
            mut_res = mmgbsa_score_v2(MUT_RECEPTOR, mut_pose, smi,
                                     WORKDIR / f"{name}_mut", minimize_iter=args.iters)
            print(f"   mut MM-GBSA: dG={mut_res['dG_bind_kcal_per_mol']:+.2f}  "
                  f"({time.time()-t0:.0f}s)")
        except Exception as e:
            print(f"   mut MM-GBSA FAILED: {e}")
            mut_res = None

        t0 = time.time()
        try:
            wt_res = mmgbsa_score_v2(WT_RECEPTOR, wt_pose, smi,
                                     WORKDIR / f"{name}_wt", minimize_iter=args.iters)
            print(f"   wt  MM-GBSA: dG={wt_res['dG_bind_kcal_per_mol']:+.2f}  "
                  f"({time.time()-t0:.0f}s)")
        except Exception as e:
            print(f"   wt MM-GBSA FAILED: {e}")
            wt_res = None

        ddG = None
        if mut_res and wt_res:
            ddG = wt_res["dG_bind_kcal_per_mol"] - mut_res["dG_bind_kcal_per_mol"]

        results.append({
            "name": name,
            "smiles": smi,
            "vina_mut": round(vina_mut, 3),
            "mmgbsa_mut": mut_res["dG_bind_kcal_per_mol"] if mut_res else None,
            "mmgbsa_wt":  wt_res["dG_bind_kcal_per_mol"]  if wt_res  else None,
            "ddG_mmgbsa": None if ddG is None else round(ddG, 2),
        })

    # Sort by mmgbsa_mut (best/most negative first)
    results.sort(key=lambda r: (r["mmgbsa_mut"] is None, r["mmgbsa_mut"] or 0))

    fieldnames = ["name","smiles","vina_mut","mmgbsa_mut","mmgbsa_wt","ddG_mmgbsa"]
    with MMGBSA_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results)

    # --- Combined report (Vina ΔΔG vs MM-GBSA ΔΔG) ---
    print("\n" + "=" * 78)
    print(f"{'Name':<13} {'Vina mut':>9} {'MMG mut':>9} {'MMG wt':>9} "
          f"{'Vina ΔΔG':>10} {'MMG ΔΔG':>10}  agreement")
    print("-" * 78)

    vina_ddG = {}
    if SELECTIVITY_CSV.exists():
        for row in csv.DictReader(SELECTIVITY_CSV.open()):
            try:
                vina_ddG[row["name"]] = float(row["ddG_wt_minus_mut"])
            except (ValueError, TypeError, KeyError):
                pass

    for r in results:
        vd = vina_ddG.get(r["name"])
        md = r["ddG_mmgbsa"]
        if vd is None or md is None:
            agree = "—"
        else:
            agree = ("✓ both selective" if vd > 0.5 and md > 1.0
                     else "≈ both ~tie" if abs(vd) <= 0.5 and abs(md) <= 1.0
                     else "✗ disagree" if (vd > 0) != (md > 0)
                     else "weak")
        print(f"{r['name']:<13} {r['vina_mut']:>9.2f} "
              f"{(r['mmgbsa_mut'] if r['mmgbsa_mut'] is not None else 0):>9.2f} "
              f"{(r['mmgbsa_wt']  if r['mmgbsa_wt']  is not None else 0):>9.2f} "
              f"{(vd if vd is not None else 0):>10.2f} "
              f"{(md if md is not None else 0):>10.2f}  {agree}")

    print("=" * 78)
    print(f"\n✓ mmgbsa.csv → {MMGBSA_CSV}")


if __name__ == "__main__":
    main()
