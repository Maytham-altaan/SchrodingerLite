"""Desmond — molecular dynamics (Schrödinger Desmond equivalent, OpenMM-powered)."""
from __future__ import annotations
import argparse
from pathlib import Path

from .utils import banner, get_logger, new_project


def run_md(pdb_in: Path, output_dir: Path,
           ns: float = 1.0, temperature: float = 300.0,
           timestep_fs: float = 2.0, report_every_ps: float = 10.0,
           solvate: bool = True, box_padding_nm: float = 1.0,
           platform: str = "auto") -> dict:
    """Solvate + equilibrate + production MD with OpenMM (Amber14)."""
    from openmm.app import (PDBFile, ForceField, Modeller, Simulation,
                            PME, HBonds, DCDReporter, StateDataReporter,
                            CheckpointReporter, PDBReporter)
    from openmm import (LangevinMiddleIntegrator, MonteCarloBarostat,
                        unit, Platform)

    output_dir.mkdir(parents=True, exist_ok=True)
    log = get_logger("desmond")

    log.info(f"Reading {pdb_in}")
    pdb = PDBFile(str(pdb_in))
    ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff)
    if solvate:
        log.info(f"Solvating (TIP3P, {box_padding_nm} nm padding) …")
        modeller.addSolvent(ff, model='tip3p',
                            padding=box_padding_nm * unit.nanometer,
                            ionicStrength=0.15 * unit.molar)

    system = ff.createSystem(modeller.topology,
                             nonbondedMethod=PME,
                             nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=HBonds)
    system.addForce(MonteCarloBarostat(1 * unit.bar,
                                        temperature * unit.kelvin, 25))
    integrator = LangevinMiddleIntegrator(temperature * unit.kelvin,
                                          1.0 / unit.picosecond,
                                          timestep_fs * unit.femtoseconds)

    if platform == "auto":
        plat = None
    else:
        plat = Platform.getPlatformByName(platform)
    sim = Simulation(modeller.topology, system, integrator, plat) if plat else \
          Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)

    log.info("Minimizing …")
    sim.minimizeEnergy(maxIterations=500)

    eq_steps = int(100 * unit.picoseconds / (timestep_fs * unit.femtoseconds))
    log.info(f"Equilibrating ({eq_steps} steps) …")
    sim.context.setVelocitiesToTemperature(temperature * unit.kelvin)
    sim.step(eq_steps)

    prod_steps = int(ns * 1000 * unit.picoseconds / (timestep_fs * unit.femtoseconds))
    report_steps = int(report_every_ps * unit.picoseconds / (timestep_fs * unit.femtoseconds))

    traj = output_dir / "trajectory.dcd"
    csv  = output_dir / "report.csv"
    fin  = output_dir / "final.pdb"

    sim.reporters.append(DCDReporter(str(traj), report_steps))
    sim.reporters.append(StateDataReporter(str(csv), report_steps,
                                            step=True, time=True,
                                            potentialEnergy=True,
                                            kineticEnergy=True,
                                            totalEnergy=True,
                                            temperature=True,
                                            volume=True, density=True,
                                            speed=True))
    sim.reporters.append(StateDataReporter(__import__("sys").stdout, report_steps*10,
                                            step=True, time=True,
                                            temperature=True, speed=True))

    log.info(f"Production: {ns} ns ({prod_steps} steps) …")
    sim.step(prod_steps)

    state = sim.context.getState(getPositions=True)
    with open(fin, "w") as fh:
        PDBFile.writeFile(sim.topology, state.getPositions(), fh, keepIds=True)

    log.info(f"Final structure → {fin.name}")
    return {"trajectory": traj, "report": csv, "final": fin}


def main(argv=None):
    p = argparse.ArgumentParser(prog="desmond",
        description="Desmond — molecular dynamics (OpenMM, Amber14).")
    p.add_argument("pdb", type=Path, help="Prepared protein PDB")
    p.add_argument("-o", "--output", type=Path, default=None)
    p.add_argument("--ns", type=float, default=1.0, help="Production length (ns)")
    p.add_argument("--temperature", type=float, default=300.0)
    p.add_argument("--no-solvate", action="store_true")
    p.add_argument("--platform", default="auto",
                   choices=["auto", "CPU", "OpenCL", "CUDA", "Metal"])
    args = p.parse_args(argv)

    banner("Desmond — molecular dynamics",
           f"input={args.pdb.name}  {args.ns} ns @ {args.temperature} K")

    out = new_project("desmond") / "output" if args.output is None else args.output
    run_md(args.pdb, out, ns=args.ns, temperature=args.temperature,
           solvate=not args.no_solvate, platform=args.platform)
    print(f"\n✓ MD complete. Output: {out}")


if __name__ == "__main__":
    main()
