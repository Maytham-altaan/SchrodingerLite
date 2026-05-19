"""Open the top-ranked EGFR hits docked into the prepared receptor in Chrome.

Reads final_ranked.csv, picks the best pose of each of the top-N ligands
(MODEL 1 from each *_poses.pdbqt), strips Vina headers, and feeds receptor +
poses to the Maestro HTML viewer.

Usage:
    python view_top_hits.py            # top 5
    python view_top_hits.py --top 10
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
from pathlib import Path

# Make the SchrodingerLite library importable when running this file directly
HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]                       # .../SchrodingerLite
sys.path.insert(0, str(ROOT))
from lib.viewer import html_view             # noqa: E402

RECEPTOR_PDB = ROOT / "projects" / "egfr_demo" / "receptor" / "protein_prepared.pdb"
ROUNDS_DIR   = HERE.parent / "rounds"
FINAL_CSV_UNIQUE = HERE.parent / "final_ranked_unique.csv"  # preferred
FINAL_CSV    = HERE.parent / "final_ranked.csv"
OUT_HTML     = HERE.parent / "top_hits_view.html"


def _canonical(smiles: str) -> str:
    """Canonical SMILES via RDKit if available; raw lowercase fallback."""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        m = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(m, canonical=True) if m else smiles
    except ImportError:
        return smiles


def find_best_pose(name: str) -> Path | None:
    """Locate <round>/docking/poses/<name>_poses.pdbqt across all rounds.
    Prefer the highest-numbered round (most recent docking)."""
    matches = sorted(ROUNDS_DIR.glob(f"round_*/docking/poses/{name}_poses.pdbqt"),
                     reverse=True)
    return matches[0] if matches else None


_MODEL_BLOCK = re.compile(
    r"^MODEL\s+1\s*\n(.*?)^ENDMDL",
    flags=re.MULTILINE | re.DOTALL,
)


def extract_first_pose(poses_pdbqt: Path, out_path: Path, lig_name: str) -> Path:
    """Write the MODEL 1 block of a Vina poses.pdbqt to a single-model .pdb.

    Translates PDBQT ATOM lines to plain PDB by trimming the autodock-specific
    columns past col 66, so 3Dmol.js can render the molecule as a normal model.
    """
    text = poses_pdbqt.read_text(errors="replace")
    m = _MODEL_BLOCK.search(text)
    body = m.group(1) if m else text  # fallback: whole file

    pdb_lines = []
    for ln in body.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            # Keep the first 66 chars, force HETATM, rename residue → LIG
            row = ln.ljust(80)
            row = "HETATM" + row[6:17] + "LIG" + row[20:]
            pdb_lines.append(row[:66].rstrip())
    pdb_lines.append("END")

    out_path.write_text("\n".join(pdb_lines) + "\n", encoding="utf-8")
    return out_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=5, help="How many top hits to overlay (default 5)")
    ap.add_argument("--no-open", action="store_true")
    args = ap.parse_args()

    csv_path = FINAL_CSV_UNIQUE if FINAL_CSV_UNIQUE.exists() else FINAL_CSV
    if not csv_path.exists():
        sys.exit(f"No ranked CSV found (looked for {FINAL_CSV_UNIQUE} and {FINAL_CSV})")
    if not RECEPTOR_PDB.exists():
        sys.exit(f"Receptor not found at {RECEPTOR_PDB}")
    print(f"Reading: {csv_path.name}")

    # Dedup by canonical SMILES, keep best (lowest) score per molecule
    best: dict[str, dict] = {}
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            c = row.get("canonical_smiles") or _canonical(row["smiles"])
            score = float(row["score"])
            if c not in best or score < float(best[c]["score"]):
                best[c] = row
    top_rows = sorted(best.values(), key=lambda r: float(r["score"]))[:args.top]
    top: list[tuple[str, str]] = [(r["name"], r["score"]) for r in top_rows]

    workdir = HERE.parent / "_top_pose_cache"
    workdir.mkdir(exist_ok=True)

    entries: list[tuple[Path, str]] = [(RECEPTOR_PDB, "receptor")]
    print(f"Top {len(top)} hits:")
    for name, score in top:
        poses = find_best_pose(name)
        if poses is None:
            print(f"  · {name}  (score {score})   ! no pose file found")
            continue
        single = workdir / f"{name}_best.pdb"
        extract_first_pose(poses, single, name)
        entries.append((single, "ligand"))
        print(f"  · {name}  (score {score})   {poses.relative_to(ROUNDS_DIR.parent)}")

    title = f"EGFR T790M/L858R — top {len(entries)-1} docked hits"
    out = html_view(entries, out_html=OUT_HTML, title=title,
                    open_browser=not args.no_open)
    print(f"\n✓ Viewer HTML written → {out}")


if __name__ == "__main__":
    main()
