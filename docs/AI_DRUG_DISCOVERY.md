# Finding New Drug Molecules with Claude AI

The `aidiscover` module turns SchrodingerLite into a closed-loop AI drug discovery
platform. Claude proposes novel molecules → local code validates and docks them →
top hits go back to Claude → Claude designs better analogues. Repeat.

## What it actually does (under the hood)

```
┌────────────────────────────────────────────────────────────────────────┐
│  ROUND N                                                               │
│                                                                        │
│   target description ─────┐                                            │
│   prior top hits ─────────┤                                            │
│                           ▼                                            │
│                      ┌─────────┐   "design 20 novel molecules"         │
│                      │ Claude  │                                       │
│                      │  API    │                                       │
│                      └────┬────┘                                       │
│                           │  JSON: [{smiles, rationale, ...}, ...]     │
│                           ▼                                            │
│                  ┌─────────────────┐  RDKit                            │
│                  │  ai_filter.py   │  • valid SMILES                   │
│                  │                 │  • Lipinski / Veber               │
│                  │                 │  • SA score (synthesizable?)      │
│                  │                 │  • PAINS / Brenk / NIH alerts     │
│                  │                 │  • novelty vs reference library   │
│                  └────────┬────────┘                                   │
│                           │  ~50% survive                              │
│                           ▼                                            │
│                  ┌─────────────────┐                                   │
│                  │   ligprep.py    │  RDKit + Meeko                    │
│                  │                 │  SMILES → 3D conformer → PDBQT    │
│                  └────────┬────────┘                                   │
│                           ▼                                            │
│                  ┌─────────────────┐                                   │
│                  │  docking.py     │  AutoDock Vina                    │
│                  │                 │  → poses + ΔG (kcal/mol)          │
│                  └────────┬────────┘                                   │
│                           ▼                                            │
│                  top 5 hits  ─────────────┐                            │
│                                            │ feedback to Claude        │
└────────────────────────────────────────────┼───────────────────────────┘
                                              ▼
                                        ROUND N+1
```

After K rounds, you get:
- `final_ranked.csv` — every molecule, sorted by predicted binding affinity.
- `summary.txt` — Claude's medicinal-chemistry interpretation of the top 10.
- `rounds/round_NN/docking/poses/` — 3D poses you can open in Maestro.

## Prerequisites

1. **SchrodingerLite installed** (`bash install.sh`).
2. **Anthropic API key** for Claude:

   ```bash
   # Get a key at https://console.anthropic.com/settings/keys
   export ANTHROPIC_API_KEY='sk-ant-...'
   # Persist it:
   echo "export ANTHROPIC_API_KEY='sk-ant-...'" >> ~/.zshrc
   ```

   Typical cost: a 3-round × 20-molecule run uses ~30k input + 15k output tokens
   on Sonnet — well under $1.

3. A **prepared receptor PDBQT** (use `prepwizard <PDBID>` to make one).

## Three ways to run

### A. The one-line shortcut (EGFR demo)

```bash
bash "/Users/maythamaltaan/Desktop/briefcase/AI DOCKING/SchrodingerLite/examples/aidiscover_quickrun.sh"
```

Downloads PDB 4ZAU (EGFR T790M/L858R), preps it, runs 3 rounds of AI design
against the ATP-binding site, opens the top poses in your browser. ~15 minutes.

### B. Your own target, free-form

```bash
# 1. Prep the target (any PDB ID)
prepwizard 6LU7 --extract-ligand 02J --keep-hetero -o ./mpro_prep

# 2. Run AI discovery — describe the target in plain English
aidiscover \
  --target "SARS-CoV-2 main protease (Mpro). Covalent-irreversible inhibitors
            are common (warheads vs Cys145), but we want REVERSIBLE binders.
            Exploit the S1/S2 hydrophobic pockets and the catalytic dyad
            (Cys145-His41). MW 300-500, logP 1-4." \
  --receptor ./mpro_prep/receptor.pdbqt \
  --ref-ligand ./mpro_prep/ref_ligand_02J.pdb \
  --rounds 3 --per-round 20 --top-k 5
```

### C. As a Python script (full control)

```python
from schrodingerlite.ai_pipeline import run_discovery
from pathlib import Path

run_discovery(
    target_description=Path("my_target.txt").read_text(),
    receptor_pdbqt=Path("./prep/receptor.pdbqt"),
    ref_ligand=Path("./prep/ref_ligand.pdb"),
    output_dir=Path("./my_run"),
    rounds=5,
    per_round=25,
    top_k_feedback=5,
    reference_library=Path("known_actives.smi"),  # for novelty scoring
    sa_max=5.5,                                    # tighter synth-accessibility
    novelty_min=0.5,                               # stricter novelty
)
```

## Writing a good target description

This is the most important input. The better your description, the better
Claude's designs. Include:

1. **Target identity** — UniProt ID or PDB code if you have one.
2. **Disease context** — gives Claude SAR priors from training data.
3. **Binding-site features** — key residues, H-bond donors/acceptors,
   hydrophobic pockets, gatekeeper/selectivity residues.
4. **Hard requirements** — physchem ranges, scaffolds to avoid, warheads OK?
5. **Negative examples** — drugs whose scaffolds you DON'T want copied.

See [examples/egfr_aidiscover.txt](../examples/egfr_aidiscover.txt) for a
worked example.

## Tuning knobs

| Flag | Default | When to change |
|------|---------|----------------|
| `--rounds` | 3 | More rounds = more convergence but more API cost. 5 is a sweet spot for hard targets. |
| `--per-round` | 20 | Bigger = more diversity, more cost. 10-30 is sensible. |
| `--top-k` | 5 | How many top hits Claude sees as feedback. Increase to 8-10 for fragment-like targets where many similar binders are useful. |
| `--exhaustiveness` | 8 | Vina docking thoroughness. 16-32 for final rounds; 4 for fast triage. |
| `--sa-max` | 6.0 | Synthetic accessibility cutoff (1 easy, 10 hard). 5.5 for medchem-friendly. |
| `--novelty-min` | 0.4 | 1−Tanimoto vs reference library. 0.5+ for serious novelty. |
| `--reference-library` | (none) | Path to SMILES file of known actives. Hits too similar are filtered out. |

## What the output looks like

```
final_ranked.csv
┌──────┬──────────────────┬─────────────────────────────────────┬───────┐
│ rank │ name             │ smiles                              │ score │
├──────┼──────────────────┼─────────────────────────────────────┼───────┤
│  1   │ NM_R3_07         │ Cc1nn(C)c2c1ncn2-c3ccc(N4CCOCC4)cc3 │ -11.2 │
│  2   │ NM_R3_02         │ Fc1cc(Nc2ncc3ccccc3n2)ccc1OC        │ -10.8 │
│  3   │ NM_R2_11         │ ...                                 │ -10.5 │
└──────┴──────────────────┴─────────────────────────────────────┴───────┘

summary.txt
"The top 3 hits share a pyrazolo[3,4-d]pyrimidine hinge binder, with the
 morpholine in NM_R3_07 reaching into the solvent-exposed front pocket.
 The fluorinated phenyl in NM_R3_02 occupies the gatekeeper region as
 designed for T790M selectivity. Watch for CYP3A4 inhibition risk in
 NM_R3_07; the morpholine-aniline motif is a known soft spot. Suggested
 next steps: synthesize NM_R3_07 and a Boc-protected piperidine analogue,
 measure WT vs T790M EGFR IC50 ratio in a kinase assay…"
```

## Honest limitations

- **Vina scoring is approximate.** A −10 kcal/mol Vina score does not mean
  10 nM affinity. Use it for ranking, not absolute prediction.
- **Claude is not a magic oracle.** It draws on med-chem patterns from
  training data. Wildly novel targets (orphan GPCRs, undruggable interfaces)
  will produce more speculative output. Always sanity-check before synthesis.
- **No PK/PD or toxicity guarantees.** Lipinski + PAINS + SA filters catch
  the obvious garbage, not subtle metabolic liabilities.
- **Synthesis isn't tested.** The SA score is a heuristic. A retrosynthesis
  tool (e.g., AiZynthFinder) is a good follow-up before committing to a hit.

## Cost & time guide

| Setting | API cost | Time |
|---------|----------|------|
| 3 rounds × 15 mols, exhaustiveness 8 | ~$0.30 | 8-15 min |
| 5 rounds × 30 mols, exhaustiveness 16 | ~$1.50 | 30-60 min |
| 10 rounds × 50 mols, exhaustiveness 32 | ~$5 | 2-4 h |

Most of the wall-clock time is docking, not the API. If you want to go faster,
turn down `--exhaustiveness` for early rounds and crank it up for the final one.

## See also

- [Main README](../README.md) — overview of all modules.
- [QUICKSTART](../QUICKSTART.md) — basic install + first docking run.
- [examples/egfr_aidiscover.txt](../examples/egfr_aidiscover.txt) — sample target description.
- [examples/aidiscover_quickrun.sh](../examples/aidiscover_quickrun.sh) — end-to-end demo script.
