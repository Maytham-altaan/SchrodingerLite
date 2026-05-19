# SchrodingerLite — QUICKSTART

## One-time install

Open Terminal and run **this single line**:

```bash
bash "/Users/maythamaltaan/Desktop/briefcase/AI DOCKING/SchrodingerLite/install.sh"
```

Takes 5–15 minutes. After it finishes, **close Terminal and open a new one.**

## Verify

```bash
schrodinger version          # should print: SchrodingerLite v1.0.0
schrodinger --help           # lists all 10 modules
```

If `schrodinger: command not found`, run:
```bash
source ~/.zshrc
```

## Your first docking run (HIV protease demo)

```bash
schrodinger workflow run "/Users/maythamaltaan/Desktop/briefcase/AI DOCKING/SchrodingerLite/examples/hiv_protease_demo.yml"
```

This will:
1. Download PDB 1HSG (HIV-1 protease + indinavir).
2. Strip waters, add hydrogens at pH 7.4, build the receptor PDBQT.
3. Extract the co-crystal ligand (MK1) to define the docking box.
4. Prepare 5 HIV-protease inhibitor analogues (amprenavir, ritonavir, etc.).
5. Dock all 5 into the active site with AutoDock Vina.
6. Save ranked CSV + pose SDF files inside `projects/`.

When done, you'll see something like:
```
✓ Workflow complete. See projects/2026-05-18_…/output/docking/docking_results.csv

  Top results:
    1. saquinavir                     -10.4 kcal/mol
    2. indinavir                       -9.8 kcal/mol
    3. ritonavir                       -9.2 kcal/mol
    ...
```

## View poses in 3D

```bash
maestro projects/<your-run>/output/docking/poses/saquinavir_poses.sdf \
        projects/<your-run>/output/receptor/protein_prepared.pdb
```
Opens an interactive 3D viewer in your default browser.

## Where everything lives

- **Install:** `/Users/maythamaltaan/Desktop/briefcase/AI DOCKING/SchrodingerLite/`
- **All run outputs:** `…/SchrodingerLite/projects/<timestamp>_<module>/`
- **Examples + demo SMILES:** `…/SchrodingerLite/examples/`
- **Full docs:** `…/SchrodingerLite/README.md`
