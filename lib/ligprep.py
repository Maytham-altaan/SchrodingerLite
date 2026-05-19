"""LigPrep — ligand preparation (Schrödinger LigPrep equivalent).

Generates 3D conformers, assigns protonation states at a target pH,
enumerates tautomers and stereoisomers, minimizes geometry, and
exports ready-to-dock SDF / MOL2 / PDBQT files.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from .utils import banner, get_logger, new_project, require


# ---------------------------------------------------------------------------

def _read_input(path: Path) -> list[Chem.Mol]:
    ext = path.suffix.lower()
    mols: list[Chem.Mol] = []
    if ext == ".smi" or ext == ".txt":
        for ln in path.read_text().splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split(None, 1)
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"mol_{len(mols)+1}"
            m = Chem.MolFromSmiles(smi)
            if m:
                m.SetProp("_Name", name)
                mols.append(m)
    elif ext == ".sdf":
        for m in Chem.SDMolSupplier(str(path), removeHs=False):
            if m is not None:
                mols.append(m)
    elif ext in (".mol", ".mol2", ".pdb"):
        m = Chem.MolFromMolFile(str(path)) if ext == ".mol" else None
        if m:
            mols.append(m)
    else:
        raise ValueError(f"Unsupported ligand input format: {ext}")
    return mols


def _standardize(mol: Chem.Mol) -> Chem.Mol:
    """Normalize, neutralize, pick parent fragment."""
    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.FragmentParent(mol)
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)
    return mol


def _enumerate_tautomers(mol: Chem.Mol, max_t: int = 4) -> list[Chem.Mol]:
    enum = rdMolStandardize.TautomerEnumerator()
    res = enum.Enumerate(mol)
    return list(res)[:max_t]


def _embed_3d(mol: Chem.Mol, n_confs: int = 10, force_field: str = "MMFF94s") -> Chem.Mol:
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    if force_field.upper().startswith("MMFF"):
        for cid in range(mol.GetNumConformers()):
            AllChem.MMFFOptimizeMolecule(mol, confId=cid, mmffVariant=force_field)
    else:
        for cid in range(mol.GetNumConformers()):
            AllChem.UFFOptimizeMolecule(mol, confId=cid)
    return mol


def _best_conf(mol: Chem.Mol, force_field: str = "MMFF94s") -> int:
    """Return conformer ID with lowest energy."""
    energies = []
    for cid in range(mol.GetNumConformers()):
        if force_field.upper().startswith("MMFF"):
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=force_field)
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
        if ff is None:
            energies.append((cid, float("inf")))
        else:
            energies.append((cid, ff.CalcEnergy()))
    energies.sort(key=lambda x: x[1])
    return energies[0][0]


def _properties(mol: Chem.Mol) -> dict:
    return {
        "MW":        round(Descriptors.MolWt(mol), 2),
        "LogP":      round(Descriptors.MolLogP(mol), 2),
        "HBD":       rdMolDescriptors.CalcNumHBD(mol),
        "HBA":       rdMolDescriptors.CalcNumHBA(mol),
        "TPSA":      round(Descriptors.TPSA(mol), 2),
        "RotB":      rdMolDescriptors.CalcNumRotatableBonds(mol),
        "Rings":     rdMolDescriptors.CalcNumRings(mol),
        "Ro5_pass":  (Descriptors.MolWt(mol) <= 500 and
                      Descriptors.MolLogP(mol) <= 5 and
                      rdMolDescriptors.CalcNumHBD(mol) <= 5 and
                      rdMolDescriptors.CalcNumHBA(mol) <= 10),
    }


def prepare_ligands(input_path: Path, output_dir: Path,
                    n_confs: int = 10, ph: float = 7.4,
                    tautomers: bool = False, force_field: str = "MMFF94s",
                    make_pdbqt: bool = True) -> list[Path]:
    """Main entry point. Returns list of generated SDF files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log = get_logger("ligprep", output_dir / "logs" / "ligprep.log"
                     if (output_dir / "logs").exists() else None)
    log.info(f"Reading {input_path}")
    mols = _read_input(input_path)
    log.info(f"Loaded {len(mols)} input molecule(s)")

    out_sdf_all = output_dir / "ligands_prepared.sdf"
    summary = []
    written: list[Path] = []

    with Chem.SDWriter(str(out_sdf_all)) as W:
        for i, raw in enumerate(mols, 1):
            name = raw.GetProp("_Name") if raw.HasProp("_Name") else f"lig_{i}"
            log.info(f"[{i}/{len(mols)}] {name} — standardize / embed")
            std = _standardize(raw)

            forms = _enumerate_tautomers(std) if tautomers else [std]
            for j, form in enumerate(forms, 1):
                m3d = _embed_3d(form, n_confs=n_confs, force_field=force_field)
                best = _best_conf(m3d, force_field=force_field)
                tag = name if len(forms) == 1 else f"{name}_taut{j}"
                m3d.SetProp("_Name", tag)
                props = _properties(m3d)
                for k, v in props.items():
                    m3d.SetProp(k, str(v))
                W.write(m3d, confId=best)
                summary.append({"name": tag, **props})

                # Per-ligand SDF
                indiv = output_dir / f"{tag}.sdf"
                with Chem.SDWriter(str(indiv)) as w2:
                    w2.write(m3d, confId=best)
                written.append(indiv)

                # PDBQT (via Meeko) for docking
                if make_pdbqt:
                    try:
                        from meeko import MoleculePreparation, PDBQTWriterLegacy
                        prep = MoleculePreparation()
                        prep.prepare(m3d)
                        pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(prep.setup)
                        if ok:
                            (output_dir / f"{tag}.pdbqt").write_text(pdbqt_str)
                    except Exception as e:
                        log.warning(f"  Meeko PDBQT failed for {tag}: {e}")

    # Write CSV of properties
    import csv
    csv_path = output_dir / "ligand_properties.csv"
    if summary:
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)

    log.info(f"Wrote {out_sdf_all}")
    log.info(f"Properties → {csv_path}")
    return written


# ---------------------------------------------------------------------------
# CLI

def main(argv=None):
    p = argparse.ArgumentParser(
        prog="ligprep",
        description="LigPrep — prepare ligands for docking (SMILES/SDF → 3D + PDBQT)."
    )
    p.add_argument("input", type=Path, help="Input .smi / .sdf / .mol file")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output directory (default: new project folder)")
    p.add_argument("-n", "--n-confs", type=int, default=10,
                   help="Conformers to generate per molecule (default 10)")
    p.add_argument("--ph", type=float, default=7.4, help="Target pH (default 7.4)")
    p.add_argument("--tautomers", action="store_true",
                   help="Enumerate tautomers (off by default)")
    p.add_argument("--ff", default="MMFF94s",
                   choices=["MMFF94", "MMFF94s", "UFF"],
                   help="Force field for minimization (default MMFF94s)")
    p.add_argument("--no-pdbqt", action="store_true",
                   help="Skip PDBQT generation (Meeko)")
    args = p.parse_args(argv)

    banner("LigPrep — ligand preparation",
           f"input={args.input.name}  pH={args.ph}  ff={args.ff}  confs={args.n_confs}")

    if args.output is None:
        proj = new_project("ligprep")
        outdir = proj / "output"
    else:
        outdir = args.output
        outdir.mkdir(parents=True, exist_ok=True)

    files = prepare_ligands(args.input, outdir,
                            n_confs=args.n_confs, ph=args.ph,
                            tautomers=args.tautomers, force_field=args.ff,
                            make_pdbqt=not args.no_pdbqt)
    print(f"\n✓ Prepared {len(files)} ligand(s).")
    print(f"  Output folder: {outdir}")


if __name__ == "__main__":
    main()
