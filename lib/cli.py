"""Master 'schrodinger' dispatcher.

Usage:
    schrodinger <module> [args...]
    schrodinger workflow run <yaml>

Modules:
    maestro      view structures
    prepwizard   protein preparation
    ligprep      ligand preparation
    glide        molecular docking
    prime        refinement + MM-GBSA rescoring
    desmond      molecular dynamics
    jaguar       quantum mechanics
    macromodel   conformational search
    canvas       cheminformatics
    strike       QSAR
    aidiscover   AI de novo drug discovery (Claude + Vina)
    workflow     run a YAML pipeline (prep → dock → score)
    version
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

MODULES = {
    "maestro":    ("viewer",      "main"),
    "prepwizard": ("protprep",    "main"),
    "ligprep":    ("ligprep",     "main"),
    "glide":      ("docking",     "main"),
    "prime":      ("refine",      "main"),
    "desmond":    ("md",          "main"),
    "jaguar":     ("qm",          "main"),
    "macromodel": ("confsearch",  "main"),
    "canvas":     ("cheminfo",    "main"),
    "strike":     ("qsar",        "main"),
    "aidiscover": ("ai_pipeline", "main"),
}


def show_help():
    print(__doc__)
    print("Examples:")
    print("  schrodinger prepwizard 1HSG --extract-ligand MK1 --keep-hetero")
    print("  schrodinger ligprep ligands.smi -n 20")
    print("  schrodinger glide -r receptor.pdbqt -l lig.pdbqt --ref-ligand ref.pdb")
    print("  schrodinger maestro receptor.pdb poses.sdf")
    print("  schrodinger workflow run pipeline.yml")


def run_workflow(yaml_path: Path):
    """Tiny YAML-driven pipeline runner: prepwizard → ligprep → glide."""
    import yaml
    cfg = yaml.safe_load(yaml_path.read_text())

    from . import protprep, ligprep, docking
    from .utils import new_project, banner

    banner("Workflow", f"file={yaml_path.name}")
    proj = new_project(cfg.get("name", "workflow"))
    out  = proj / "output"; out.mkdir(parents=True, exist_ok=True)

    # 1. Protein
    rec_cfg = cfg["receptor"]
    rec_out = out / "receptor"; rec_out.mkdir(exist_ok=True)
    rec = protprep.prepare_protein(
        rec_cfg["source"], rec_out,
        ph=rec_cfg.get("ph", 7.4),
        keep_water=rec_cfg.get("keep_water", False),
        keep_hetero=rec_cfg.get("keep_hetero", True),
        minimize=rec_cfg.get("minimize", False),
        extract_ligand_resname=rec_cfg.get("extract_ligand"),
    )

    # 2. Ligands
    lig_cfg = cfg["ligands"]
    lig_out = out / "ligands"; lig_out.mkdir(exist_ok=True)
    lig_files = ligprep.prepare_ligands(
        Path(lig_cfg["input"]), lig_out,
        n_confs=lig_cfg.get("n_confs", 10),
        ph=lig_cfg.get("ph", 7.4),
        tautomers=lig_cfg.get("tautomers", False),
        force_field=lig_cfg.get("ff", "MMFF94s"),
        make_pdbqt=True,
    )
    pdbqts = [f.with_suffix(".pdbqt") for f in lig_files
              if f.with_suffix(".pdbqt").exists()]

    # 3. Dock
    dock_cfg = cfg.get("docking", {})
    if rec.get("ligand"):
        from .docking import grid_from_ligand
        center, size = grid_from_ligand(rec["ligand"],
                                        padding=dock_cfg.get("padding", 8))
    else:
        center = tuple(dock_cfg["center"])
        size   = tuple(dock_cfg.get("size", [22,22,22]))

    dock_out = out / "docking"; dock_out.mkdir(exist_ok=True)
    csv = docking.dock_batch(rec["pdbqt"], pdbqts, dock_out,
                              center=center, size=size,
                              backend=dock_cfg.get("backend", "vina"),
                              exhaustiveness=dock_cfg.get("exhaustiveness", 8),
                              n_poses=dock_cfg.get("n_poses", 9))
    print(f"\n✓ Workflow complete. See {csv}")


def main(argv=None):
    argv = argv or sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help", "help"):
        show_help(); return 0
    if argv[0] == "version":
        from . import __version__
        print(f"SchrodingerLite v{__version__}")
        return 0
    if argv[0] == "workflow":
        if len(argv) < 3 or argv[1] != "run":
            print("Usage: schrodinger workflow run <pipeline.yml>"); return 2
        run_workflow(Path(argv[2])); return 0

    mod = argv[0]
    if mod not in MODULES:
        print(f"Unknown module: {mod}")
        show_help(); return 2
    modname, fn = MODULES[mod]
    import importlib
    m = importlib.import_module(f".{modname}", package="schrodingerlite")
    getattr(m, fn)(argv[1:])
    return 0


if __name__ == "__main__":
    sys.exit(main())
