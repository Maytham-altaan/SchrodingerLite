"""MacroModel — conformational search (Schrödinger MacroModel equivalent).

Uses RDKit ETKDGv3 + clustering and MMFF94s minimization for low-energy
conformer ensembles.
"""
from __future__ import annotations
import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from .utils import banner, get_logger, new_project


def conformer_search(mol: Chem.Mol, n_confs: int = 200,
                     prune_rms: float = 0.5,
                     ff: str = "MMFF94s",
                     keep_top: int | None = None) -> Chem.Mol:
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    params.pruneRmsThresh = prune_rms
    AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    energies = []
    if ff.upper().startswith("MMFF"):
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=ff)
        for cid in range(mol.GetNumConformers()):
            f = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
            f.Minimize()
            energies.append((cid, f.CalcEnergy()))
    else:
        for cid in range(mol.GetNumConformers()):
            f = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            f.Minimize()
            energies.append((cid, f.CalcEnergy()))

    energies.sort(key=lambda x: x[1])

    if keep_top:
        kept_ids = {cid for cid, _ in energies[:keep_top]}
        # Build a new mol containing only kept conformers
        new = Chem.Mol(mol)
        new.RemoveAllConformers()
        for cid, e in energies[:keep_top]:
            conf = mol.GetConformer(cid)
            new_cid = new.AddConformer(conf, assignId=True)
            new.GetConformer(new_cid).SetIntProp("OrigID", cid) if False else None
        # store energies as prop on each conf via molblock - simpler: keep mol but track e
        return new, energies[:keep_top]
    return mol, energies


def main(argv=None):
    p = argparse.ArgumentParser(prog="macromodel",
        description="MacroModel — conformational search.")
    p.add_argument("input", help="SMILES string or .sdf file")
    p.add_argument("-o", "--output", type=Path, default=None)
    p.add_argument("-n", "--n-confs", type=int, default=200)
    p.add_argument("--prune-rms", type=float, default=0.5,
                   help="RMS pruning threshold during embedding")
    p.add_argument("--ff", default="MMFF94s", choices=["MMFF94", "MMFF94s", "UFF"])
    p.add_argument("--keep-top", type=int, default=50,
                   help="Keep N lowest-energy conformers")
    args = p.parse_args(argv)

    inp = Path(args.input)
    if inp.exists() and inp.suffix.lower() == ".sdf":
        mol = next(Chem.SDMolSupplier(str(inp)))
    else:
        mol = Chem.MolFromSmiles(args.input if not inp.exists() else inp.read_text().strip().split()[0])
    if mol is None:
        raise SystemExit("Could not parse molecule")

    banner("MacroModel — conformational search",
           f"n_confs={args.n_confs}  ff={args.ff}  keep={args.keep_top}")

    out = new_project("macromodel") / "output" if args.output is None else args.output
    out.mkdir(parents=True, exist_ok=True)

    new_mol, energies = conformer_search(mol, n_confs=args.n_confs,
                                          prune_rms=args.prune_rms,
                                          ff=args.ff, keep_top=args.keep_top)
    sdf = out / "conformers.sdf"
    with Chem.SDWriter(str(sdf)) as W:
        for i, (_, e) in enumerate(energies):
            new_mol.SetProp("Energy_kcal_per_mol", f"{e:.3f}")
            new_mol.SetProp("ConfRank", str(i+1))
            W.write(new_mol, confId=i)
    print(f"\n✓ Wrote {len(energies)} conformers → {sdf}")
    print(f"  Lowest E: {energies[0][1]:.3f} kcal/mol")
    print(f"  Highest E: {energies[-1][1]:.3f} kcal/mol")


if __name__ == "__main__":
    main()
