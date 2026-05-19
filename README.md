# SchrodingerLite

**An open-source local replica of the Schrödinger Suite, runnable from any terminal on your Mac.**

Built for pharmacists, medicinal chemists, and computational biologists who need
Maestro/Glide/LigPrep capabilities without a Schrödinger license.

> 📄 **Companion preprint:** *Iterative LLM-driven de novo drug design with multi-method consensus scoring: a case study on EGFR T790M/L858R.*
> Manuscript and figures in [`projects/egfr_demo/aidiscovery/`](projects/egfr_demo/aidiscovery/) ([`preprint.md`](projects/egfr_demo/aidiscovery/preprint.md), [`preprint.docx`](projects/egfr_demo/aidiscovery/preprint.docx)).
> Triple-confirmed lead (Vina + MM-GBSA + 2 ns MD): **EGFRm-009**, a 3-(2-trifluoromethyl-4-fluorobenzyl) pyridopyrimidinone. Reproduce in ~3 hours on a laptop with the EGFR demo below.

---

## Modules (and their Schrödinger counterparts)

| Command       | Schrödinger module        | Powered by                |
|---------------|---------------------------|---------------------------|
| `maestro`     | Maestro (3D viewer)       | Py3Dmol + PyMOL           |
| `prepwizard`  | Protein Preparation Wizard| PDBFixer + OpenBabel      |
| `ligprep`     | LigPrep                   | RDKit + Meeko + OpenBabel |
| `glide`       | Glide docking             | AutoDock Vina + Smina     |
| `prime`       | Prime / MM-GBSA           | OpenMM (Amber14 + GBSA)   |
| `desmond`     | Desmond MD                | OpenMM (PME, TIP3P)       |
| `jaguar`      | Jaguar QM                 | Psi4 + xtb                |
| `macromodel`  | MacroModel confsearch     | RDKit ETKDGv3 + MMFF94s   |
| `canvas`      | Canvas cheminformatics    | RDKit + scikit-learn      |
| `strike`      | Strike QSAR               | scikit-learn + XGBoost    |
| `aidiscover`  | **(new)** AI de novo drug discovery | **Claude API** + RDKit + Vina |

All commands are first-class CLI tools available **anywhere** in any terminal:

```bash
schrodinger glide -r rec.pdbqt -l lig.pdbqt --ref-ligand ref.pdb
# equivalently:
glide -r rec.pdbqt -l lig.pdbqt --ref-ligand ref.pdb
```

---

## Install

```bash
git clone https://github.com/Maytham-altaan/SchrodingerLite.git
cd SchrodingerLite
bash install.sh
```

What it does:
1. Checks for **conda** (installs Miniforge via Homebrew if missing).
2. Creates a conda env named `schrodinger-lite` from `environment.yml`.
3. Installs the `schrodingerlite` Python package.
4. Symlinks all 11 launchers into `/usr/local/bin` (or `~/.local/bin`),
   so you can call them from anywhere.

Time: 5–15 minutes (most of it is conda resolving dependencies).
Disk: ~3.5 GB.

After installation, open a **new terminal** and verify:

```bash
schrodinger version
schrodinger --help
```

---

## Quick start (5 lines, full docking pipeline)

```bash
# 1. Prepare the protein (downloads 1HSG, extracts the bound ligand)
prepwizard 1HSG --extract-ligand MK1 --keep-hetero

# 2. Prepare your library (SMILES file → 3D + PDBQT)
ligprep examples/ligands.smi -n 20 -o ./prepared_ligs

# 3. Dock everything
glide -r <receptor.pdbqt> -l ./prepared_ligs/*.pdbqt \
      --ref-ligand <ref_ligand_MK1.pdb> --exhaustiveness 16

# 4. View the top poses in your browser
maestro <receptor.pdb> ./poses/best_pose.sdf

# 5. (Optional) Rescore the top hit with MM-GBSA
prime mmgbsa <complex.pdb>
```

Or run **the whole pipeline in one shot** from a YAML file:

```bash
schrodinger workflow run examples/hiv_protease_demo.yml
```

---

## Finding NEW drug molecules with Claude AI

The `aidiscover` module is a **closed-loop AI drug-discovery pipeline**:
Claude designs novel molecules → RDKit filters them → Vina docks them →
the top hits go back to Claude → Claude designs better analogues. Repeat.

```bash
# Set your Anthropic API key once
export ANTHROPIC_API_KEY='sk-ant-...'

# Run the EGFR T790M demo (3 rounds × 15 molecules, ~15 minutes, < $1 API)
bash examples/aidiscover_quickrun.sh
```

Full details and tuning knobs:
**[docs/AI_DRUG_DISCOVERY.md](docs/AI_DRUG_DISCOVERY.md)**

---

## Where do my files go?

Every command creates a **timestamped project folder** in:

```
$SCHRODINGER_LITE_HOME/projects/   # defaults to the cloned repo's projects/ dir
```

Each project has:
```
2026-05-18_142312_glide/
├── input/        ← copies of your inputs
├── output/       ← results (SDF poses, CSV ranks, PDBQT, etc.)
└── logs/         ← log files
```

Override the output dir at any time with `-o /path/of/your/choice`.

---

## Tips for pharmacists

- **Pose validity:** the `posebusters` Python package (installed with this) flags
  geometrically implausible docking poses. Use it to filter Glide output.
- **Reference ligand grid:** if your PDB has a co-crystal ligand, passing
  `--ref-ligand` is almost always better than guessing a center.
- **Lipinski / Veber:** `canvas descriptors` produces a CSV with pass/fail
  for both rule sets, plus QED drug-likeness.
- **MM-GBSA caveat:** Schrödinger's Prime MM-GBSA uses OPLS + VSGB2.1; this
  replica uses Amber14 + OBC2. Ranks correlate; absolute ΔG values do not.
- **Free GPU acceleration:** Desmond MD will automatically use your Mac's GPU
  (Metal/OpenCL) via OpenMM.

---

## Uninstall

```bash
conda env remove -n schrodinger-lite
rm /usr/local/bin/{schrodinger,maestro,prepwizard,ligprep,glide,prime,desmond,jaguar,macromodel,canvas,strike} 2>/dev/null
# or, for ~/.local/bin
rm ~/.local/bin/{schrodinger,maestro,prepwizard,ligprep,glide,prime,desmond,jaguar,macromodel,canvas,strike} 2>/dev/null
```

The SchrodingerLite folder itself can be deleted freely once removed.

---

## License & differences from real Schrödinger

This is **not** affiliated with Schrödinger, LLC. It is an independent open-source
toolkit using free academic and community packages that solve the same problems:

| Capability                  | Schrödinger | SchrodingerLite |
|----------------------------|-------------|-----------------|
| OPLS force field            | ✅          | ❌ (Amber14)    |
| GPU-accelerated Desmond MD  | ✅          | ✅ (OpenMM)     |
| Validated MM-GBSA           | ✅          | ⚠️ approximation |
| Glide SP/XP scoring         | ✅          | ❌ (Vina)       |
| WaterMap, FEP+              | ✅          | ❌              |
| Free, local, runs anywhere  | ❌          | ✅              |

For production drug discovery work, the commercial Schrödinger Suite remains
the gold standard. SchrodingerLite is excellent for teaching, methods
prototyping, and most academic virtual screening campaigns.
