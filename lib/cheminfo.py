"""Canvas — cheminformatics (Schrödinger Canvas equivalent).

Descriptors, fingerprints, similarity search, clustering, Lipinski/Veber
filters, pharmacophore-style queries, scaffold analysis.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Draw
from rdkit.Chem.Scaffolds import MurckoScaffold

from .utils import banner, new_project


# ---- descriptors -----------------------------------------------------------

DESCRIPTORS = {
    "MW":      Descriptors.MolWt,
    "LogP":    Descriptors.MolLogP,
    "HBD":     rdMolDescriptors.CalcNumHBD,
    "HBA":     rdMolDescriptors.CalcNumHBA,
    "TPSA":    Descriptors.TPSA,
    "RotB":    rdMolDescriptors.CalcNumRotatableBonds,
    "Rings":   rdMolDescriptors.CalcNumRings,
    "AromaticRings": rdMolDescriptors.CalcNumAromaticRings,
    "HeavyAtoms": Chem.Mol.GetNumHeavyAtoms,
    "FCsp3":   rdMolDescriptors.CalcFractionCSP3,
    "QED":     None,  # filled below to avoid recompute
}

def _compute(mol: Chem.Mol) -> dict:
    from rdkit.Chem import QED
    out = {}
    for k, fn in DESCRIPTORS.items():
        try:
            if k == "QED":
                out[k] = round(QED.qed(mol), 3)
            else:
                out[k] = round(float(fn(mol)), 3)
        except Exception:
            out[k] = None
    out["Lipinski_pass"] = (out["MW"] <= 500 and out["LogP"] <= 5
                            and out["HBD"] <= 5 and out["HBA"] <= 10)
    out["Veber_pass"] = (out["RotB"] <= 10 and out["TPSA"] <= 140)
    out["Scaffold"] = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    return out


def descriptors_table(input_file: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    ext = input_file.suffix.lower()
    mols = []
    if ext == ".sdf":
        for m in Chem.SDMolSupplier(str(input_file)):
            if m: mols.append(m)
    else:
        for ln in input_file.read_text().splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("#"): continue
            parts = ln.split(None, 1)
            m = Chem.MolFromSmiles(parts[0])
            if m:
                m.SetProp("_Name", parts[1] if len(parts)>1 else f"mol_{len(mols)+1}")
                mols.append(m)

    rows = []
    for m in mols:
        d = {"Name": m.GetProp("_Name") if m.HasProp("_Name") else "",
             "SMILES": Chem.MolToSmiles(m)}
        d.update(_compute(m))
        rows.append(d)

    out_csv = output_dir / "descriptors.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    return out_csv


# ---- similarity ------------------------------------------------------------

def similarity_search(query_smiles: str, library: Path,
                      output_dir: Path, top: int = 20,
                      radius: int = 2, nbits: int = 2048) -> Path:
    qmol = Chem.MolFromSmiles(query_smiles)
    qfp = AllChem.GetMorganFingerprintAsBitVect(qmol, radius, nBits=nbits)

    hits = []
    if library.suffix.lower() == ".sdf":
        for m in Chem.SDMolSupplier(str(library)):
            if not m: continue
            fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nbits)
            tan = DataStructs.TanimotoSimilarity(qfp, fp)
            hits.append((m.GetProp("_Name") if m.HasProp("_Name") else Chem.MolToSmiles(m),
                         Chem.MolToSmiles(m), tan))
    else:
        for ln in library.read_text().splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("#"): continue
            parts = ln.split(None, 1)
            m = Chem.MolFromSmiles(parts[0])
            if not m: continue
            fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nbits)
            tan = DataStructs.TanimotoSimilarity(qfp, fp)
            hits.append((parts[1] if len(parts)>1 else parts[0], parts[0], tan))

    hits.sort(key=lambda x: -x[2])
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / "similarity.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["rank","name","smiles","tanimoto"])
        for r,(n,s,t) in enumerate(hits[:top], 1):
            w.writerow([r,n,s,f"{t:.3f}"])
    return out


# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(prog="canvas",
        description="Canvas — cheminformatics, descriptors, similarity.")
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("descriptors", help="Compute descriptor table from SMI/SDF")
    a.add_argument("input", type=Path)
    a.add_argument("-o","--output", type=Path, default=None)

    b = sub.add_parser("similarity", help="Tanimoto similarity search")
    b.add_argument("query", help="Query SMILES")
    b.add_argument("library", type=Path)
    b.add_argument("--top", type=int, default=20)
    b.add_argument("-o","--output", type=Path, default=None)

    args = p.parse_args(argv)
    banner(f"Canvas — {args.cmd}")
    out = (new_project(f"canvas_{args.cmd}") / "output") if args.output is None else args.output

    if args.cmd == "descriptors":
        path = descriptors_table(args.input, out)
        print(f"\n✓ Descriptors → {path}")
    else:
        path = similarity_search(args.query, args.library, out, top=args.top)
        print(f"\n✓ Similarity hits → {path}")


if __name__ == "__main__":
    main()
