"""Jaguar — quantum mechanics (Schrödinger Jaguar equivalent).

Supports:
  - Single-point energies, geometry optimization, frequencies (Psi4)
  - Fast semi-empirical optimization (xtb)
  - Common DFT functionals + basis sets
"""
from __future__ import annotations
import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from .utils import banner, get_logger, new_project


def _xyz_block(mol: Chem.Mol) -> str:
    """RDKit Mol → XYZ string (first conformer)."""
    if mol.GetNumConformers() == 0:
        m = Chem.AddHs(mol)
        AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(m)
        mol = m
    conf = mol.GetConformer()
    lines = []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():2s}  {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(lines)


def psi4_calc(mol: Chem.Mol, output_dir: Path,
              method: str = "b3lyp", basis: str = "6-31g*",
              task: str = "energy",
              charge: int = 0, multiplicity: int = 1) -> dict:
    """Run a Psi4 calculation."""
    import psi4
    output_dir.mkdir(parents=True, exist_ok=True)
    psi4.set_output_file(str(output_dir / "psi4.out"), append=False)
    psi4.set_memory("4 GB")

    geom = f"{charge} {multiplicity}\n{_xyz_block(mol)}\n"
    m = psi4.geometry(geom)
    spec = f"{method}/{basis}"

    if task == "energy":
        e = psi4.energy(spec, molecule=m)
        return {"energy_Eh": float(e), "method": spec}
    elif task == "optimize":
        e = psi4.optimize(spec, molecule=m)
        return {"energy_Eh": float(e), "method": spec, "optimized": True}
    elif task == "frequency":
        e, wfn = psi4.frequency(spec, molecule=m, return_wfn=True)
        freqs = list(wfn.frequencies().to_array())
        return {"energy_Eh": float(e), "method": spec, "frequencies_cm-1": freqs}
    else:
        raise ValueError(f"Unknown task {task}")


def xtb_optimize(mol: Chem.Mol, output_dir: Path) -> dict:
    """Fast semi-empirical optimization with xtb."""
    import subprocess
    output_dir.mkdir(parents=True, exist_ok=True)
    xyz = output_dir / "input.xyz"
    natoms = mol.GetNumAtoms() if mol.GetNumConformers() else Chem.AddHs(mol).GetNumAtoms()
    xyz.write_text(f"{natoms}\n\n{_xyz_block(mol)}\n")
    result = subprocess.run(["xtb", str(xyz), "--opt", "tight"],
                            cwd=output_dir, capture_output=True, text=True)
    (output_dir / "xtb.log").write_text(result.stdout)
    energy = None
    for ln in result.stdout.splitlines():
        if "TOTAL ENERGY" in ln:
            energy = float(ln.split()[-3])
    return {"energy_Eh": energy, "method": "GFN2-xTB"}


# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(prog="jaguar",
        description="Jaguar — quantum mechanics (Psi4 / xtb).")
    p.add_argument("input", help="SMILES string, .smi, or .sdf file")
    p.add_argument("-o", "--output", type=Path, default=None)
    p.add_argument("--engine", choices=["psi4", "xtb"], default="xtb",
                   help="QM engine. xtb (semi-empirical, default) is installed by default; "
                        "psi4 (ab initio DFT) requires manual install: "
                        "`conda install -n schrodinger-lite -c conda-forge psi4`")
    p.add_argument("--method", default="b3lyp")
    p.add_argument("--basis", default="6-31g*")
    p.add_argument("--task", choices=["energy", "optimize", "frequency"],
                   default="energy")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--mult", type=int, default=1)
    args = p.parse_args(argv)

    # Load molecule
    inp = Path(args.input)
    if inp.exists() and inp.suffix.lower() == ".sdf":
        mol = next(Chem.SDMolSupplier(str(inp), removeHs=False))
    elif inp.exists():
        mol = Chem.MolFromSmiles(inp.read_text().strip().split()[0])
    else:
        mol = Chem.MolFromSmiles(args.input)
    if mol is None:
        raise SystemExit(f"Could not parse molecule from: {args.input}")

    banner("Jaguar — quantum mechanics",
           f"engine={args.engine}  task={args.task}  "
           f"{args.method}/{args.basis if args.engine=='psi4' else ''}")

    out = new_project("jaguar") / "output" if args.output is None else args.output
    out.mkdir(parents=True, exist_ok=True)

    if args.engine == "psi4":
        res = psi4_calc(mol, out, method=args.method, basis=args.basis,
                        task=args.task, charge=args.charge,
                        multiplicity=args.mult)
    else:
        res = xtb_optimize(mol, out)

    import json
    (out / "result.json").write_text(json.dumps(res, indent=2))
    print(f"\n✓ {res}")
    print(f"  Output: {out}")


if __name__ == "__main__":
    main()
