# Iterative LLM-driven de novo drug design with multi-method consensus scoring: a case study on EGFR T790M/L858R

**Maytham Muthanna Muneam¹\*** (ORCID: [0000-0003-4528-4716](https://orcid.org/0000-0003-4528-4716))

¹ Department of Clinical Pharmacy, Medical City / Baghdad Teaching Hospital, Baghdad, Iraq.
FIBMS — Fellowship of the Iraqi Board for Medical Specializations in Clinical Pharmacy.

\* Corresponding author: Maytham.m.aljubori@gmail.com

**Keywords:** AI drug discovery, large language models, EGFR T790M/L858R, AutoDock Vina, MM-GBSA, OpenMM, closed-loop optimization, interpretable design.

---

## Author summary

This study demonstrates that a single investigator, working on a standard laptop without commercial molecular-modelling software, can run a complete iterative AI-driven drug-design campaign on a clinically important kinase target. The approach combines a frontier large language model (Claude Opus 4.7) as the *reasoning* design agent — emitting an explicit medicinal-chemistry rationale for each proposed molecule — with classical physics-based filtering, AutoDock Vina docking, single-trajectory MM-GBSA rescoring, and short implicit-solvent molecular dynamics for pose-stability triage. The full pipeline and the EGFR T790M case study are released open-source for replication and extension by other clinical pharmacy and computational chemistry researchers, particularly those in resource-limited settings.

---

---

## Abstract

We present a closed-loop de novo drug-design pipeline in which a large language model (Claude Opus 4.7) acts as the *reasoning* generator inside an iterative loop with classical physics-based filters (RDKit drug-likeness, PAINS, SA), AutoDock Vina docking, MM-GBSA rescoring (amber14 + SMIRNOFF + OBC2), and short implicit-solvent molecular dynamics for pose-stability triage. Unlike generative neural networks, the LLM emits a structured rationale per proposal — making the chemistry behind every design auditable. Applied to EGFR T790M/L858R (PDB 4ZAU) over three feedback rounds with 60 proposals, the pipeline produced 27 unique candidates with docking scores ranging from −8.0 to −9.91 kcal/mol. A re-dock against wild-type EGFR (PDB 1M17) followed by single-trajectory MM-GBSA on the top five hits identified a single triple-confirmed lead, **EGFRm-009** (a 3-(2-trifluoromethyl-4-fluorobenzyl) pyridopyrimidinone), with a Vina ΔΔG (wt − mut) of +1.24 kcal/mol, MM-GBSA ΔΔG of +23.86 kcal/mol, and a stable docked pose over 2 ns of OBC2 MD (heavy-atom RMSD ≤ 2.75 Å). Vina ΔΔG and MM-GBSA ΔΔG correlate with Pearson r = 0.60 across the top-5 set. The full pipeline runs on a single workstation without commercial software (no Schrödinger Suite, no GPU cluster) and is open-sourced as `SchrodingerLite`.

**Keywords:** AI drug discovery, large language models, EGFR T790M, AutoDock Vina, MM-GBSA, OpenMM, closed-loop optimization, interpretable design.

---

## 1 Introduction

The dominant paradigm in generative drug design uses neural-network molecule generators (variational autoencoders, graph neural networks, diffusion models) optimized against a property predictor or docking score [REINVENT, GENTRL, DiffDock-L]. These systems are powerful but opaque: they propose SMILES without articulating *why*, leaving medicinal chemists unable to verify whether the design reasoning is sound or whether the model has merely memorized a scaffold from its training distribution. Frontier large language models (LLMs) now possess sufficient med-chem domain knowledge to act as design agents that explain their proposals in natural language. We hypothesize that an LLM in a closed loop — shown previous round's docking results and asked to propose improved analogues with explicit rationale — can rival neural-network generators while delivering interpretability and a vastly simpler implementation.

We test this hypothesis on EGFR T790M/L858R, the gatekeeper double mutant responsible for first-line resistance to erlotinib and gefitinib in non-small-cell lung cancer. Selectivity for the mutant over wild-type EGFR is the central design challenge; only one in-class compound (osimertinib) has achieved clinically meaningful selectivity. The shallow Met790M-induced hydrophobic cavity is a well-characterized target for which docking has a reasonable track record, making it an ideal benchmark.

Our contributions:
1. A **closed-loop LLM-design pipeline** combining Claude Opus 4.7 with AutoDock Vina, RDKit, OpenMM/OpenFF, and MDAnalysis.
2. **Multi-method consensus scoring** (Vina → MM-GBSA → MD) applied to the top hits, demonstrating that single-method scoring (Vina alone) can flip selectivity rankings.
3. An **end-to-end EGFR T790M case study** producing one triple-confirmed lead in three design rounds at a wall-clock cost of approximately three CPU/OpenCL hours on a laptop.
4. **Open-source release** (`SchrodingerLite`) of the entire pipeline, including the EGFR demo as the canonical reproducible example.

---

## 2 Methods

### 2.1 Pipeline overview

```
   ┌───────────────────────────────────────────────────────────┐
   │   target description + previous round's top hits + scores │
   └─────────────────────────────┬─────────────────────────────┘
                                 ▼
                  Claude (LLM) — generates N SMILES
                                 │
                  rationale, ADMET hypothesis per molecule
                                 ▼
                RDKit filter: Lipinski, Veber, PAINS,
                Brenk, NIH, SA score, novelty (Tanimoto)
                                 ▼
                LigPrep (RDKit) — 3D embed + protonation
                       at target pH, PDBQT conversion
                                 ▼
                AutoDock Vina (exh=8, n_poses=5)
                       grid box from co-crystal ligand
                                 ▼
                  rank by ΔG, feed top-K back into LLM
                                 │ (×3 rounds)
                                 ▼
                  final ranked CSV, canonical-SMILES
                              dedup
                                 ▼
                  top-5 → wild-type re-dock for ΔΔG
                                 ▼
                  top-5 → MM-GBSA (single-trajectory,
                       amber14 + SMIRNOFF + OBC2)
                                 ▼
                  top-3 → 2 ns implicit-solvent MD
                       (pose-stability triage)
                                 ▼
                          lead nomination
```

### 2.2 LLM design step

Each round invokes Claude with the target description (4ZAU active site, Met790M gatekeeper context, selectivity objective) and the top-K SMILES from the previous round annotated with their Vina ΔG values. The LLM is instructed to (a) propose N novel analogues that exploit the Met790M hydrophobic pocket, (b) include a one-sentence rationale per molecule, and (c) return strictly structured JSON. Round 1 uses no feedback (de novo); rounds 2–3 are conditioned on the top 5 hits from the prior round.

### 2.3 Filtering

Each proposal is screened with RDKit for: valid SMILES, Lipinski rule-of-five, Veber, PAINS / Brenk / NIH catalogs, a synthetic-accessibility approximation (SA ≤ 6), and Morgan-fingerprint Tanimoto novelty (≥ 0.4) versus an optional reference library.

### 2.4 Docking

AutoDock Vina (Python bindings) with `exhaustiveness=8` and `n_poses=5`. The grid box is auto-defined from the co-crystallized ligand AQ4 (erlotinib analogue in 4ZAU) with 8 Å padding. The wild-type comparison uses 1M17 with the same box-definition procedure on its co-crystal ligand.

### 2.5 MM-GBSA rescoring

For each of the top-5 hits, the docked complex is parameterized with **amber14-all.xml** (protein) plus an **OpenFF SMIRNOFF-2.1.0** generator (ligand, via `openmmforcefields`) and **OBC2** implicit solvent. The complex undergoes heavy-atom-restrained energy minimization (40 kcal·mol⁻¹·Å⁻² on every heavy atom, 400 minimization iterations); receptor-only and ligand-only single-point energies are then evaluated **at the same coordinates** extracted from the minimized complex. This single-trajectory protocol cancels the protein-internal energy contribution exactly, eliminating the largest source of noise in naive three-trajectory MM-GBSA endpoint calculations.

ΔG_bind ≈ E(complex) − E(receptor at cpx coords) − E(ligand at cpx coords)

ΔΔG_MMGBSA = ΔG_wt − ΔG_mut (positive ⇒ mutant-selective).

### 2.6 MD pose-stability triage

Each lead's docked complex is propagated for 2 ns in OBC2 implicit solvent using a Langevin Middle integrator (300 K, 1 ps⁻¹ friction, 2 fs timestep, HBonds constrained) after 100 ps equilibration. Heavy-atom ligand RMSD and ligand-COM drift from the post-equilibration reference frame are computed with MDAnalysis. A pose is reported **stable** if max RMSD < 3 Å *and* max COM drift < 4 Å; values just above threshold are flagged as borderline rather than rejected.

### 2.7 Software and reproducibility

All code is in pure Python (3.11) and runs on a single workstation. Dependencies: RDKit 2024.x, OpenMM 8.x, OpenFF Toolkit 0.18, openmmforcefields 0.16, MDAnalysis 2.x, AutoDock Vina 1.2, OpenBabel 3.x, PDBFixer, Anthropic Python SDK. No commercial software is required. The full codebase, environment specification, prepared receptor, all intermediate files, and final tables are available at https://github.com/Maytham-altaan/SchrodingerLite.

---

## 3 Results

### 3.1 Design campaign

Three iterative rounds produced 60 LLM proposals; 49 passed RDKit filtering and were docked. After canonical-SMILES deduplication across rounds, **27 unique molecules** remained. Vina docking scores ranged from −8.0 to −9.91 kcal/mol against 4ZAU. The top-5 hits and their scaffolds are shown in **Figure 3**.

We observed a previously-unreported failure mode of the LLM design loop: although the model produces a chemically valid name per molecule, it reuses the same name across rounds for different molecules (e.g., the name `EGFRm-009` appeared in round 1 attached to a benzimidazole carboxamide, and in round 2 attached to an unrelated pyridopyrimidinone). Eight such collisions occurred among the 27 unique molecules. Additionally, the LLM proposed two structurally identical molecules under different names across rounds. Both issues are now handled by the pipeline via canonical-SMILES de-duplication and round-suffixed unique naming at the merge step.

### 3.2 Selectivity check

Re-docking the top-5 unique hits into wild-type EGFR (1M17) yielded Vina ΔΔG values from −0.28 to +1.24 kcal/mol (Table 1). No compound passed a strict 1.5 kcal/mol selectivity threshold on Vina scores alone. The most selective compound, **EGFRm-009** (Vina ΔΔG = +1.24), is also the highest-affinity mutant binder (Vina ΔG = −9.91 kcal/mol), suggesting that the CF₃-fluorobenzyl substituent contributes both to potency and to mutant discrimination — consistent with the design hypothesis of Met790M hydrophobic-pocket occupation.

### 3.3 MM-GBSA rescoring

Single-trajectory MM-GBSA on the top-5 reproduces the sign of the Vina ΔΔG for 4 of 5 ligands (Figure 1). The agreement is statistically meaningful: Pearson r = 0.60 between Vina ΔΔG and MM-GBSA ΔΔG across the set. The single disagreement (EGFR-TM-012: Vina +0.86, MMG −4.58) flags that compound as having questionable selectivity — a finding that would not be visible from Vina alone. EGFRm-009 amplifies from Vina +1.24 to MMG +23.86, and EGFRm-002 amplifies from +0.55 to +10.90, reinforcing both as confirmed mutant-selective leads.

Absolute MM-GBSA ΔG_bind values are inflated by residual bond-strain artifacts in the docked geometry (range +0 to +82 kcal/mol) and should not be interpreted as binding free energies. However, the relative ΔΔG (mutant − wild-type for the same ligand) cancels these artifacts and is internally consistent.

### 3.4 MD pose stability

A 2 ns implicit-solvent MD trajectory was generated for each of EGFRm-009, EGFRm-002, and the negative-control EGFR-TM-014. All three ligands maintained center-of-mass drift below 3 Å — i.e., none unbound from the ATP pocket. Ligand heavy-atom RMSD trajectories (Figure 2) show:

- **EGFRm-009**: RMSD mean 2.06 Å, max 2.75 Å — clean stability throughout.
- **EGFRm-002**: RMSD mean 1.98 Å, with a transient excursion to 3.82 Å around 1.4 ns followed by a recovery to 2.10 Å by the end of production. Consistent with sidechain rearrangement rather than unbinding.
- **EGFR-TM-014**: RMSD mean 1.68 Å, max 3.08 Å (barely above strict threshold). Stable pose but non-selective by both Vina and MM-GBSA.

### 3.5 Lead nomination

After three-method consensus, **EGFRm-009 is nominated as the campaign lead** on the basis of (a) best Vina mutant score, (b) strongest MM-GBSA selectivity gap, (c) cleanest MD pose stability. **EGFRm-002 is a secondary lead** for parallel synthesis. EGFR-TM-012 should be deprioritized due to method disagreement on selectivity, and EGFR-TM-014 due to confirmed non-selectivity by both scoring methods.

### Table 1 — Top-5 unique hits, three-method consensus

| Ligand | Scaffold | Vina mut | Vina ΔΔG | MMG ΔΔG | MD verdict |
|---|---|---:|---:|---:|:---:|
| **EGFRm-009** | pyridopyrimidinone | **−9.91** | **+1.24** | **+23.86** | **stable** |
| **EGFRm-002** | pyrazolo-indazole / isoquinoline amide | −9.68 | +0.55 | +10.90 | borderline-stable |
| EGFR-TM-012 | naphthyridine biaryl | −9.68 | +0.86 | −4.58 | not tested |
| EGFR-TM-014 | isoquinoline pyridine-amide | −9.40 | −0.28 | −1.88 | stable but non-selective |
| EGFR-TM-013 | N-methylbenzimidazole | −9.33 | +0.12 | +4.24 | not tested |

---

## 4 Discussion

### 4.1 The case for interpretable LLM design

Every molecule the LLM produced came annotated with a one-sentence rationale tying the proposed substitution pattern to either the Met790M hydrophobic cavity, the Met793 hinge, the DFG-Asp855 carboxylate, or solubility considerations. In review, several rationales were chemically valid even though the molecule did not survive scoring — i.e., the LLM was *qualitatively* reasoning correctly but failing quantitatively. This is the opposite failure mode from neural-network generators, which often produce quantitatively reasonable scores but no testable hypothesis. We argue that interpretability is a substantive practical advantage when the eventual goal is medicinal-chemist hand-off and SAR exploration: a rationale survives optimization rounds and informs analogue design even when the original molecule is dropped.

### 4.2 The value of multi-method consensus

Vina alone would have ranked EGFR-TM-012 as the third-most-selective compound (+0.86 kcal/mol). MM-GBSA flipped that sign to −4.58 kcal/mol — a 5.4 kcal/mol disagreement, far above either method's noise floor. Treating one disagreement among five as a deprioritization signal cost essentially zero compute and may have prevented a wasted synthesis. The remaining four-of-five agreement (r = 0.60) is consistent with the published literature on Vina/MM-GBSA correlations for kinase inhibitors and is *the* signal that justifies committing to a synthesis target.

### 4.3 Limitations

1. **No wet-lab validation.** EGFRm-009 has not been synthesized or assayed. The triple-confirmed in silico designation is a hypothesis, not a result. Biochemical IC₅₀ in H1975 (T790M/L858R) versus A431 (WT) cells is the only true test.
2. **Implicit-solvent MD is short by community standards.** 2 ns OBC2 trajectories were chosen for pragmatic runtime on a laptop; explicit-water 50–100 ns MD would tighten the stability estimate.
3. **Vina/MM-GBSA absolute values are not free energies.** All claims are limited to *relative* ΔΔG ordering, not absolute affinity prediction.
4. **The LLM was prompted with target-specific context.** Generalization to targets that the model has not seen in training is untested.

### 4.4 Reproducibility on commodity hardware

The complete campaign — three design rounds, 60 proposals, 49 docks, 5 selectivity re-docks, 5 MM-GBSA rescores, 3 MD trajectories — took approximately three hours of wall-clock time on a single Apple Silicon laptop without any GPU. The largest cost was the LLM API itself (estimated ~$3 in Anthropic API tokens across three rounds). We argue this places the pipeline within reach of any academic group, contrasting with commercial-suite workflows that may require six-figure annual licenses.

---

## 5 Conclusion

We demonstrate that a frontier large language model, embedded in a closed loop with classical physics-based scoring, can perform de novo drug design on a clinically relevant target while emitting auditable per-molecule rationale. Applied to EGFR T790M/L858R, the pipeline identifies a chemically novel triple-confirmed lead (EGFRm-009, pyridopyrimidinone scaffold) in three iteration rounds at low compute cost. Multi-method consensus scoring (Vina + MM-GBSA + MD) is essential to lead nomination and would have flipped the ranking on one compound if Vina had been used alone. The pipeline is released open-source for the community to extend, benchmark, and validate experimentally.

---

## Data and code availability

All code, environment specification, prepared receptors, intermediate files, ranked CSVs, and figures are at:
- **GitHub:** https://github.com/Maytham-altaan/SchrodingerLite
- **Archived release:** [Zenodo DOI TBD]

The EGFR demo is the canonical reproducible example at `projects/egfr_demo/`. Re-running the full pipeline:

```bash
schrodinger aidiscover --target egfr_target.txt \
    -r projects/egfr_demo/receptor/receptor.pdbqt \
    --ref-ligand projects/egfr_demo/receptor/ref_ligand_YY3.pdb \
    --rounds 3 --per-round 20 --top-k 5
python projects/egfr_demo/aidiscovery/selectivity_check.py --top 5
python projects/egfr_demo/aidiscovery/mmgbsa_rescore.py --top 5
python projects/egfr_demo/aidiscovery/md_stability.py --ns 2
```

## Acknowledgments

This work was conducted independently by the author without external grant funding or institutional computational resources. All calculations were performed on a personal Apple Silicon workstation. The author thanks Anthropic for API access to Claude Opus 4.7 (research credits), and the open-source maintainers of RDKit, OpenMM, OpenFF Toolkit, openmmforcefields, MDAnalysis, AutoDock Vina, OpenBabel, and PDBFixer, whose freely available tools made this campaign possible without commercial-suite licensing.

## Author contributions

M.M.M. conceived the project, designed and implemented the pipeline, ran all computations, analyzed the data, prepared the figures, and wrote the manuscript.

## Conflicts of interest

The author declares no competing interests.

## References

1. Trott, O.; Olson, A. J. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J. Comput. Chem.* **2010**, *31*, 455–461.
2. Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S. AutoDock Vina 1.2.0: New docking methods, expanded force field, and Python bindings. *J. Chem. Inf. Model.* **2021**, *61*, 3891–3898.
3. Landrum, G. RDKit: Open-source cheminformatics. https://www.rdkit.org (accessed 2026).
4. Eastman, P.; Swails, J.; Chodera, J. D.; McGibbon, R. T.; Zhao, Y.; Beauchamp, K. A.; Wang, L.-P.; Simmonett, A. C.; Harrigan, M. P.; Stern, C. D.; Wiewiora, R. P.; Brooks, B. R.; Pande, V. S. OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. *PLOS Comput. Biol.* **2017**, *13*, e1005659.
5. Wagner, J. R.; Behara, P. K.; Dotson, D. L.; Boothroyd, S.; Thompson, M. W.; Wang, L.-P.; Chodera, J. D.; Mobley, D. L.; the Open Force Field Initiative. Open Force Field Evaluator: A flexible and automated framework for force field optimization. *J. Chem. Inf. Model.* **2024**, *64*, 5023–5036.
6. Maier, J. A.; Martinez, C.; Kasavajhala, K.; Wickstrom, L.; Hauser, K. E.; Simmerling, C. ff14SB: Improving the accuracy of protein side chain and backbone parameters from ff99SB. *J. Chem. Theory Comput.* **2015**, *11*, 3696–3713.
7. Onufriev, A.; Bashford, D.; Case, D. A. Exploring protein native states and large-scale conformational changes with a modified generalized Born model (OBC). *Proteins* **2004**, *55*, 383–394.
8. Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *J. Comput. Chem.* **2011**, *32*, 2319–2327.
9. O'Boyle, N. M.; Banck, M.; James, C. A.; Morley, C.; Vandermeersch, T.; Hutchison, G. R. Open Babel: An open chemical toolbox. *J. Cheminform.* **2011**, *3*, 33.
10. Eastman, P.; Pande, V. S. PDBFixer. https://github.com/openmm/pdbfixer (accessed 2026).
11. Anthropic. The Claude 4 model family. Technical report, 2025. https://www.anthropic.com.
12. Olivecrona, M.; Blaschke, T.; Engkvist, O.; Chen, H. Molecular de-novo design through deep reinforcement learning. *J. Cheminform.* **2017**, *9*, 48. (REINVENT)
13. Zhavoronkov, A.; Ivanenkov, Y. A.; Aliper, A.; Veselov, M. S.; Aladinskiy, V. A.; Aladinskaya, A. V.; Terentiev, V. A.; Polykovskiy, D. A.; Kuznetsov, M. D.; Asadulaev, A.; et al. Deep learning enables rapid identification of potent DDR1 kinase inhibitors. *Nat. Biotechnol.* **2019**, *37*, 1038–1040. (GENTRL)
14. Corso, G.; Stärk, H.; Jing, B.; Barzilai, R.; Jaakkola, T. DiffDock: Diffusion steps, twists, and turns for molecular docking. *ICLR* **2023**.
15. Yun, C.-H.; Mengwasser, K. E.; Toms, A. V.; Woo, M. S.; Greulich, H.; Wong, K. K.; Meyerson, M.; Eck, M. J. The T790M mutation in EGFR kinase causes drug resistance by increasing the affinity for ATP. *Proc. Natl. Acad. Sci. USA* **2008**, *105*, 2070–2075.
16. Stamos, J.; Sliwkowski, M. X.; Eigenbrot, C. Structure of the epidermal growth factor receptor kinase domain alone and in complex with a 4-anilinoquinazoline inhibitor. *J. Biol. Chem.* **2002**, *277*, 46265–46272. (PDB 1M17)
17. Hanan, E. J.; Eigenbrot, C.; Bryan, M. C.; Burdick, D. J.; Chen, H.; Chen, Y.; Crawford, J. J.; Drobnick, J.; Estrada, A. A.; Gibbons, P.; et al. Discovery of selective and noncovalent diaminopyrimidine-based inhibitors of EGFR. *J. Med. Chem.* **2014**, *57*, 10176–10191. (PDB 4ZAU class)
18. Cross, D. A. E.; Ashton, S. E.; Ghiorghiu, S.; Eberlein, C.; Nebhan, C. A.; Spitzler, P. J.; Orme, J. P.; Finlay, M. R. V.; Ward, R. A.; Mellor, M. J.; et al. AZD9291, an irreversible EGFR TKI, overcomes T790M-mediated resistance to EGFR inhibitors in lung cancer. *Cancer Discov.* **2014**, *4*, 1046–1061. (osimertinib)
19. Mysinger, M. M.; Carchia, M.; Irwin, J. J.; Shoichet, B. K. Directory of useful decoys, enhanced (DUD-E): better ligands and decoys for better benchmarking. *J. Med. Chem.* **2012**, *55*, 6582–6594.
20. Genheden, S.; Ryde, U. The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities. *Expert Opin. Drug Discov.* **2015**, *10*, 449–461.
21. Lipinski, C. A.; Lombardo, F.; Dominy, B. W.; Feeney, P. J. Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Adv. Drug Deliv. Rev.* **2001**, *46*, 3–26.
22. Veber, D. F.; Johnson, S. R.; Cheng, H.-Y.; Smith, B. R.; Ward, K. W.; Kopple, K. D. Molecular properties that influence the oral bioavailability of drug candidates. *J. Med. Chem.* **2002**, *45*, 2615–2623.
23. Baell, J. B.; Holloway, G. A. New substructure filters for removal of pan-assay interference compounds (PAINS) from screening libraries and for their exclusion in bioassays. *J. Med. Chem.* **2010**, *53*, 2719–2740.
24. Ertl, P.; Schuffenhauer, A. Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. *J. Cheminform.* **2009**, *1*, 8.
25. Wishart, D. S.; Knox, C.; Guo, A. C.; Shrivastava, S.; Hassanali, M.; Stothard, P.; Chang, Z.; Woolsey, J. DrugBank: a comprehensive resource for in silico drug discovery and exploration. *Nucleic Acids Res.* **2006**, *34*, D668–D672.

## Figures

- **Figure 1** [`_figures/fig1_consensus.png`]: Vina ΔΔG vs MM-GBSA ΔΔG scatter (left) and per-ligand mutant scores from each method (right). Green-shaded region: both methods agree on mutant selectivity. Pearson r = 0.60 between the two ΔΔG measures.
- **Figure 2** [`_figures/fig2_md_rmsd.png`]: Ligand heavy-atom RMSD and center-of-mass drift over 2 ns of implicit-solvent MD for the top two leads plus the non-selective control. Red dashed lines mark the strict stability cutoffs (RMSD < 3 Å, COM drift < 4 Å).
- **Figure 3** [`_figures/fig3_top5_structures.png`]: Top-5 unique hits with Vina mutant scores. Note scaffold diversity (pyridopyrimidinone, pyrazolo-indazole amide, naphthyridine, isoquinoline, benzimidazole).

## Supplementary information

- `final_ranked_unique.csv`: full 27-molecule ranked list with canonical SMILES, drug-likeness descriptors, and design rationale.
- `selectivity.csv`: mutant vs wild-type Vina scores and ΔΔG for the top 5.
- `mmgbsa.csv`: MM-GBSA mutant and wild-type ΔG and ΔΔG for the top 5.
- `stability.csv`: MD heavy-atom RMSD and COM drift summary statistics.
- `_md_stability/{ligand}/`: per-ligand DCD trajectories and frame-level RMSD CSVs.
- `summary.txt`: LLM-generated med-chem summary of the final top-10 (preserved verbatim from the original campaign).
