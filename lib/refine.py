"""Prime — structure refinement & MM-GBSA scoring (Schrödinger Prime equivalent).

Provides:
  - Side-chain / loop optimization via OpenMM relaxation
  - Implicit-solvent (GBSA) rescoring of docked complexes
"""
from __future__ import annotations
import argparse
from pathlib import Path

from .utils import banner, get_logger, new_project


def refine(pdb_in: Path, pdb_out: Path, max_iter: int = 1000) -> Path:
    """Energy-minimize with GB-OBC implicit solvent (Amber14)."""
    from openmm.app import PDBFile, ForceField, Modeller, Simulation, HBonds, NoCutoff
    from openmm import LangevinIntegrator, unit

    pdb = PDBFile(str(pdb_in))
    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff)
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff,
                             constraints=HBonds)
    integrator = LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond,
                                    0.002*unit.picoseconds)
    sim = Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=max_iter)
    state = sim.context.getState(getPositions=True, getEnergy=True)
    with open(pdb_out, "w") as fh:
        PDBFile.writeFile(modeller.topology, state.getPositions(), fh, keepIds=True)
    return pdb_out


def _ligand_from_smiles_and_pose(smiles: str, pose_sdf: Path):
    """Build a hydrogen-bearing RDKit Mol whose bonds + aromaticity come from
    a canonical SMILES and whose heavy-atom 3D coordinates come from a docked
    SDF pose. PDBQT→SDF via OpenBabel doesn't preserve aromaticity, so we
    treat the SDF as a heavy-atom-coord source only."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from openff.toolkit import Molecule

    # 1. Clean reference from SMILES (correct topology, no Hs yet)
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # 2. Heavy-atom positions from the docked pose
    supplier = Chem.SDMolSupplier(str(pose_sdf), removeHs=True, sanitize=False)
    pose = next((m for m in supplier if m is not None), None)
    if pose is None:
        raise ValueError(f"No pose in {pose_sdf}")
    pose_conf = pose.GetConformer()

    n_heavy_ref = ref.GetNumAtoms()
    n_heavy_pose = pose.GetNumAtoms()
    if n_heavy_pose < n_heavy_ref:
        raise ValueError(f"Pose has {n_heavy_pose} heavy atoms but SMILES has {n_heavy_ref}")

    # Match by index first (likely same atom order since pose came from this SMILES);
    # verify by element. Fall back to substructure match.
    match = list(range(n_heavy_ref))
    elem_ok = all(ref.GetAtomWithIdx(i).GetAtomicNum() ==
                  pose.GetAtomWithIdx(i).GetAtomicNum()
                  for i in range(n_heavy_ref))
    if not elem_ok:
        substr = pose.GetSubstructMatch(ref)
        if substr:
            match = list(substr)

    # 3. Plant the docked heavy-atom coords on ref (heavies only, no Hs)
    heavy_conf = Chem.Conformer(n_heavy_ref)
    for i in range(n_heavy_ref):
        heavy_conf.SetAtomPosition(i, pose_conf.GetAtomPosition(match[i]))
    ref.AddConformer(heavy_conf, assignId=True)

    # 4. Add Hs with positions inferred from the docked heavy-atom geometry
    mol_h = Chem.AddHs(ref, addCoords=True)

    # 5. Brief H-only MMFF relax so Hs sit at sensible bond angles
    try:
        mp = AllChem.MMFFGetMoleculeProperties(mol_h)
        ff = AllChem.MMFFGetMoleculeForceField(mol_h, mp)
        if ff:
            for i in range(mol_h.GetNumAtoms()):
                if mol_h.GetAtomWithIdx(i).GetAtomicNum() != 1:
                    ff.AddFixedPoint(i)
            ff.Minimize(maxIts=300)
    except Exception:
        pass

    return Molecule.from_rdkit(mol_h, allow_undefined_stereo=True)


def mmgbsa_score_v2(receptor_pdb: Path,
                    ligand_pose_sdf: Path,
                    ligand_smiles: str,
                    output_dir: Path,
                    smirnoff_ff: str = "openff-2.1.0",
                    minimize_iter: int = 600) -> dict:
    """Single-trajectory MM-GBSA-style ΔG for arbitrary small molecules.

        ΔG_bind ≈ E(complex) − E(receptor_at_cpx_coords) − E(ligand_at_cpx_coords)

    The complex is minimized with heavy-atom positional restraints (keeps the
    docked pose intact while relaxing hydrogens). The receptor-only and
    ligand-only energies are then single-point evaluations at the SAME
    coordinates extracted from the minimized complex.

    Why single-trajectory: with 3-trajectory endpoint MM-GBSA, the protein
    converges to different conformations in apo vs holo, leaving thousands of
    kcal/mol of uncancelled protein-internal energy. Single-trajectory cancels
    that exactly — only the protein-ligand interaction (vdW + Coulomb + GBSA
    desolvation) survives.

    Args:
        receptor_pdb:    clean receptor PDB (no ligand HETATM).
        ligand_pose_sdf: SDF with the docked pose (heavy-atom coords only).
        ligand_smiles:   canonical SMILES — rebuilds correct bond orders since
                         OpenBabel's PDBQT→SDF strips aromaticity.
    """
    from openmm.app import (PDBFile, ForceField, Modeller, NoCutoff, HBonds,
                            Simulation)
    from openmm import LangevinIntegrator, CustomExternalForce, unit
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    import warnings
    warnings.filterwarnings("ignore")

    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Ligand: bonds from SMILES, 3D coords from docked pose
    lig = _ligand_from_smiles_and_pose(ligand_smiles, ligand_pose_sdf)
    if lig.n_conformers == 0:
        raise ValueError(f"No 3D conformer recovered from {ligand_pose_sdf}")

    # 2. Combined force field
    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
    gen = SMIRNOFFTemplateGenerator(molecules=lig, forcefield=smirnoff_ff)
    ff.registerTemplateGenerator(gen.generator)

    # 3. Build complex topology
    rec = PDBFile(str(receptor_pdb))
    lig_top = lig.to_topology().to_openmm()
    lig_pos = lig.conformers[0].to_openmm()
    modeller = Modeller(rec.topology, rec.positions)
    modeller.add(lig_top, lig_pos)
    cpx_top = modeller.topology
    cpx_pos = modeller.positions
    n_rec = len(list(rec.topology.atoms()))
    n_tot = len(list(cpx_top.atoms()))

    # 4. Heavy-atom-restrained minimization of the complex
    cpx_system = ff.createSystem(cpx_top, nonbondedMethod=NoCutoff,
                                 constraints=HBonds)
    restraint = CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    # Strong restraint so heavy atoms barely move; Hs are free
    restraint.addGlobalParameter(
        "k", 40 * 4.184 * unit.kilojoule_per_mole / unit.angstrom**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    for i, atom in enumerate(cpx_top.atoms()):
        if atom.element.symbol == "H":
            continue
        p = cpx_pos[i].value_in_unit(unit.nanometer)
        restraint.addParticle(i, [p[0], p[1], p[2]])
    cpx_system.addForce(restraint)

    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                               0.002 * unit.picoseconds)
    sim = Simulation(cpx_top, cpx_system, integ)
    sim.context.setPositions(cpx_pos)
    sim.minimizeEnergy(maxIterations=minimize_iter)

    # Recompute energy WITHOUT the restraint contribution
    state = sim.context.getState(getPositions=True)
    cpx_min_pos = state.getPositions()
    sim_no_restraint_sys = ff.createSystem(cpx_top, nonbondedMethod=NoCutoff,
                                           constraints=HBonds)
    sim2 = Simulation(cpx_top, sim_no_restraint_sys, LangevinIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds))
    sim2.context.setPositions(cpx_min_pos)
    E_complex = float(sim2.context.getState(getEnergy=True)
                      .getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))

    # 5. Receptor single-point at the complex's protein coords
    rec_pos = cpx_min_pos[:n_rec]
    rec_system = ff.createSystem(rec.topology, nonbondedMethod=NoCutoff,
                                 constraints=HBonds)
    rec_sim = Simulation(rec.topology, rec_system, LangevinIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds))
    rec_sim.context.setPositions(rec_pos)
    E_receptor = float(rec_sim.context.getState(getEnergy=True)
                       .getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))

    # 6. Ligand single-point at the complex's ligand coords
    lig_pos_min = cpx_min_pos[n_rec:n_tot]
    lig_system = ff.createSystem(lig_top, nonbondedMethod=NoCutoff,
                                 constraints=HBonds)
    lig_sim = Simulation(lig_top, lig_system, LangevinIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds))
    lig_sim.context.setPositions(lig_pos_min)
    E_ligand = float(lig_sim.context.getState(getEnergy=True)
                     .getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))

    dG = E_complex - E_receptor - E_ligand
    return {
        "E_complex":  round(E_complex, 2),
        "E_receptor": round(E_receptor, 2),
        "E_ligand":   round(E_ligand, 2),
        "dG_bind_kcal_per_mol": round(dG, 2),
    }


def mmgbsa_score(complex_pdb: Path, output_dir: Path) -> dict:
    """Legacy single-PDB MM-GBSA — only works for ligands that are valid
    amber14 residues. Prefer `mmgbsa_score_v2` for arbitrary small molecules.
    """
    from openmm.app import PDBFile, ForceField, Modeller, Simulation, HBonds, NoCutoff
    from openmm import LangevinIntegrator, unit

    output_dir.mkdir(parents=True, exist_ok=True)
    lines = complex_pdb.read_text().splitlines()
    prot_lines = [ln for ln in lines if ln.startswith("ATOM")]
    lig_lines  = [ln for ln in lines if ln.startswith("HETATM")]
    if not prot_lines or not lig_lines:
        raise ValueError("Need both ATOM (protein) and HETATM (ligand) records.")

    prot = output_dir / "rec.pdb"
    lig  = output_dir / "lig.pdb"
    cpx  = output_dir / "cpx.pdb"
    prot.write_text("\n".join(prot_lines) + "\nEND\n")
    lig.write_text("\n".join(lig_lines)  + "\nEND\n")
    cpx.write_text("\n".join(prot_lines + lig_lines) + "\nEND\n")

    def energy(pdb_path: Path) -> float:
        pdb = PDBFile(str(pdb_path))
        ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
        modeller = Modeller(pdb.topology, pdb.positions)
        try:
            modeller.addHydrogens(ff)
        except Exception:
            pass
        sys_ = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff,
                                constraints=HBonds)
        ig = LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
        sim = Simulation(modeller.topology, sys_, ig)
        sim.context.setPositions(modeller.positions)
        sim.minimizeEnergy(maxIterations=50)
        s = sim.context.getState(getEnergy=True)
        return float(s.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))

    E_complex = energy(cpx)
    E_receptor = energy(prot)
    try:
        E_ligand = energy(lig)
    except Exception:
        E_ligand = 0.0
    dG = E_complex - E_receptor - E_ligand
    return {"E_complex": E_complex, "E_receptor": E_receptor,
            "E_ligand": E_ligand, "dG_bind_kcal_per_mol": round(dG, 2)}


def main(argv=None):
    p = argparse.ArgumentParser(prog="prime",
        description="Prime — refinement + MM-GBSA rescoring.")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("refine", help="Energy-minimize a protein structure")
    r.add_argument("pdb", type=Path)
    r.add_argument("-o","--output", type=Path, default=None)
    r.add_argument("--max-iter", type=int, default=1000)

    g = sub.add_parser("mmgbsa", help="Rescore a protein-ligand complex (MM-GBSA-like)")
    g.add_argument("complex", type=Path, help="PDB with ATOM (protein) + HETATM (ligand)")
    g.add_argument("-o","--output", type=Path, default=None)

    args = p.parse_args(argv)
    banner(f"Prime — {args.cmd}")
    out = (new_project(f"prime_{args.cmd}") / "output") if args.output is None else args.output

    if args.cmd == "refine":
        out.mkdir(parents=True, exist_ok=True)
        result = refine(args.pdb, out / "refined.pdb", max_iter=args.max_iter)
        print(f"\n✓ Refined → {result}")
    else:
        scores = mmgbsa_score(args.complex, out)
        import json
        (out / "mmgbsa.json").write_text(json.dumps(scores, indent=2))
        print(f"\n✓ ΔG(bind) ≈ {scores['dG_bind_kcal_per_mol']} kcal/mol")
        print(f"  Details: {out / 'mmgbsa.json'}")


if __name__ == "__main__":
    main()
