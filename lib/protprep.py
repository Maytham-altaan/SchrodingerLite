"""PrepWizard — protein preparation (Schrödinger Protein Preparation Wizard).

- Download from PDB or accept local file
- Remove waters / heteroatoms (optional, can keep cofactors)
- Add missing residues / atoms with PDBFixer
- Add hydrogens at target pH
- Restrained energy minimization (optional, OpenMM)
- Write cleaned PDB and Vina-ready PDBQT
"""
from __future__ import annotations
import argparse
import subprocess
import urllib.request
from pathlib import Path

from .utils import banner, get_logger, new_project, require


def fetch_pdb(pdb_id: str, dest: Path) -> Path:
    """Download a PDB file by 4-letter ID."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    out = dest / f"{pdb_id.upper()}.pdb"
    urllib.request.urlretrieve(url, out)
    return out


def clean_with_pdbfixer(pdb_in: Path, pdb_out: Path,
                        ph: float = 7.4,
                        keep_water: bool = False,
                        keep_hetero: bool = True) -> Path:
    """Use PDBFixer to add missing atoms/residues and hydrogens."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=str(pdb_in))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    if not keep_hetero:
        fixer.removeHeterogens(keepWater=keep_water)
    elif not keep_water:
        # Keep ligands but drop waters
        fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(pdb_out, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    return pdb_out


def minimize_protein(pdb_in: Path, pdb_out: Path,
                     max_iter: int = 200) -> Path:
    """Brief restrained minimization with OpenMM (Amber14)."""
    from openmm.app import PDBFile, ForceField, Modeller, Simulation, HBonds, NoCutoff
    from openmm import LangevinIntegrator, unit, Platform

    pdb = PDBFile(str(pdb_in))
    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff)
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff,
                             constraints=HBonds)
    integrator = LangevinIntegrator(300 * unit.kelvin,
                                    1.0 / unit.picosecond,
                                    0.002 * unit.picoseconds)
    sim = Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=max_iter)
    state = sim.context.getState(getPositions=True)
    with open(pdb_out, "w") as fh:
        PDBFile.writeFile(modeller.topology, state.getPositions(), fh, keepIds=True)
    return pdb_out


def pdb_to_pdbqt(pdb_in: Path, pdbqt_out: Path) -> Path:
    """Convert a prepared PDB into Vina-style PDBQT using OpenBabel."""
    require("obabel")
    cmd = ["obabel", str(pdb_in), "-O", str(pdbqt_out),
           "-xr",   # rigid receptor
           "--partialcharge", "gasteiger"]
    subprocess.run(cmd, check=True, capture_output=True)
    return pdbqt_out


_LIGAND_BLACKLIST = {
    # waters
    "HOH","WAT","H2O","DOD","OH","SOL",
    # ions
    "NA","CL","K","MG","CA","ZN","FE","MN","CU","NI","CO","CD","HG",
    "BR","IOD","F","RB","CS","BA","SR","LI",
    # buffers / cryo / crystallization additives
    "SO4","PO4","NO3","ACT","EDO","GOL","PEG","PG4","PGE","P6G",
    "DMS","TRS","MES","HEPES","FMT","ACE","BME","DTT","MPD","BTB",
    "EPE","CIT","TLA","MLA","MAL","ICT","FLC","BCT",
    # sugars / common cofactors that aren't usually "the" ligand
    "NAG","BMA","MAN","GAL","GLC","FUC","XYS","BGC","FRU",
}

def auto_detect_ligand(pdb_in: Path) -> str | None:
    """Find the largest non-water/non-ion HETATM group in a PDB.
    Returns the residue name or None."""
    counts: dict[str, int] = {}
    for ln in pdb_in.read_text().splitlines():
        if not ln.startswith("HETATM"):
            continue
        rn = ln[17:20].strip()
        if rn in _LIGAND_BLACKLIST or len(rn) == 0:
            continue
        counts[rn] = counts.get(rn, 0) + 1
    if not counts:
        return None
    return max(counts, key=counts.get)


def extract_ligand(pdb_in: Path, resname: str, out_file: Path) -> Path:
    """Pull a co-crystallized ligand out of a PDB by residue name.

    Pass 'auto' to auto-detect the largest non-water HETATM."""
    if resname.lower() == "auto":
        auto = auto_detect_ligand(pdb_in)
        if auto is None:
            raise ValueError(f"No co-crystallized ligand found in {pdb_in.name}")
        resname = auto
    lines = pdb_in.read_text().splitlines()
    kept = [ln for ln in lines
            if ln.startswith("HETATM") and ln[17:20].strip() == resname.upper()]
    if not kept:
        # Helpful error: tell them what ligands ARE in the file
        seen = sorted({ln[17:20].strip() for ln in lines
                       if ln.startswith("HETATM")
                       and ln[17:20].strip() not in _LIGAND_BLACKLIST})
        raise ValueError(
            f"No HETATM with resname {resname} in {pdb_in.name}. "
            f"Candidate ligand resnames found: {seen or '(none)'}\n"
            f"Tip: use --extract-ligand auto to pick the largest automatically."
        )
    out_file.write_text("\n".join(kept) + "\nEND\n")
    # Stamp the chosen resname onto the output filename for traceability
    return out_file


# ---------------------------------------------------------------------------

def prepare_protein(source: str, output_dir: Path,
                    ph: float = 7.4,
                    keep_water: bool = False,
                    keep_hetero: bool = True,
                    minimize: bool = False,
                    extract_ligand_resname: str | None = None) -> dict:
    """Main entry point.

    Args:
        source: PDB ID (4 chars) OR path to a local .pdb file.
    """
    log = get_logger("prepwizard")
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Fetch or copy
    src = Path(source) if Path(source).exists() else None
    if src is None and len(source) == 4 and source.isalnum():
        log.info(f"Downloading PDB {source.upper()} …")
        src = fetch_pdb(source, output_dir)
    elif src is None:
        raise FileNotFoundError(f"Not a PDB ID or existing file: {source}")
    else:
        # copy into project
        dest = output_dir / src.name
        dest.write_bytes(src.read_bytes())
        src = dest
    log.info(f"Source: {src}")

    # 2. Extract ligand if requested (before cleaning waters/hetero)
    lig_file = None
    if extract_ligand_resname:
        # Resolve 'auto' before naming the file
        resolved = extract_ligand_resname
        if resolved.lower() == "auto":
            detected = auto_detect_ligand(src)
            if detected is None:
                log.warning("auto-ligand: no co-crystal ligand found, continuing without")
                resolved = None
            else:
                log.info(f"auto-ligand: detected resname '{detected}'")
                resolved = detected
        if resolved:
            lig_file = output_dir / f"ref_ligand_{resolved.upper()}.pdb"
            extract_ligand(src, resolved, lig_file)
            log.info(f"Extracted reference ligand → {lig_file.name}")

    # 3. PDBFixer pass
    fixed = output_dir / "protein_prepared.pdb"
    log.info("Running PDBFixer (residues, atoms, hydrogens) …")
    clean_with_pdbfixer(src, fixed, ph=ph,
                        keep_water=keep_water, keep_hetero=False)
    # we drop hetero in main; ligand already extracted above
    log.info(f"Cleaned PDB → {fixed.name}")

    # 4. Optional minimization
    if minimize:
        log.info("Restrained minimization (OpenMM, Amber14) …")
        mini = output_dir / "protein_minimized.pdb"
        try:
            minimize_protein(fixed, mini)
            fixed = mini
            log.info(f"Minimized → {mini.name}")
        except Exception as e:
            log.warning(f"Minimization failed (continuing without): {e}")

    # 5. PDBQT for docking
    pdbqt = output_dir / "receptor.pdbqt"
    log.info("Converting to PDBQT (OpenBabel) …")
    pdb_to_pdbqt(fixed, pdbqt)
    log.info(f"Receptor PDBQT → {pdbqt.name}")

    return {
        "pdb":     fixed,
        "pdbqt":   pdbqt,
        "ligand":  lig_file,
        "source":  src,
    }


# ---------------------------------------------------------------------------
# CLI

def main(argv=None):
    p = argparse.ArgumentParser(
        prog="prepwizard",
        description="Protein Preparation Wizard — clean & ready a PDB for docking."
    )
    p.add_argument("source", help="4-letter PDB ID (e.g. 1HSG) or local .pdb file")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output directory (default: new project folder)")
    p.add_argument("--ph", type=float, default=7.4, help="Target pH (default 7.4)")
    p.add_argument("--keep-water", action="store_true", help="Keep crystallographic waters")
    p.add_argument("--keep-hetero", action="store_true", help="Keep cofactors / ligands")
    p.add_argument("--minimize", action="store_true",
                   help="Run brief OpenMM minimization")
    p.add_argument("--extract-ligand", default=None, metavar="RESNAME",
                   help="Extract co-crystallized ligand. Pass 'auto' to use the "
                        "largest non-water/non-ion HETATM, or a 3-letter resname (e.g. MK1).")
    args = p.parse_args(argv)

    banner("PrepWizard — protein preparation",
           f"source={args.source}  pH={args.ph}  minimize={args.minimize}")

    if args.output is None:
        proj = new_project("prepwizard")
        outdir = proj / "output"
    else:
        outdir = args.output

    result = prepare_protein(args.source, outdir,
                             ph=args.ph,
                             keep_water=args.keep_water,
                             keep_hetero=args.keep_hetero,
                             minimize=args.minimize,
                             extract_ligand_resname=args.extract_ligand)

    print("\n✓ Protein prepared.")
    print(f"  PDB:     {result['pdb']}")
    print(f"  PDBQT:   {result['pdbqt']}")
    if result['ligand']:
        print(f"  Ref lig: {result['ligand']}")


if __name__ == "__main__":
    main()
