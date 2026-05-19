# Cover letter — manuscript submission

**Maytham Muthanna Muneam, FIBMS Clinical Pharmacy**
Department of Clinical Pharmacy
Medical City / Baghdad Teaching Hospital
Baghdad, Iraq
Maytham.m.aljubori@gmail.com

[Date]

[Editor-in-Chief Name]
[Journal Name]
[Journal Address]

**RE: Manuscript submission — "Iterative LLM-driven de novo drug design with multi-method consensus scoring: a case study on EGFR T790M/L858R"**

Dear Editor,

I am pleased to submit the enclosed manuscript for consideration. The work describes a complete closed-loop pipeline for de novo small-molecule design in which a frontier large language model (Claude Opus 4.7) serves as the reasoning generator inside an iterative loop with classical physics-based scoring (AutoDock Vina, MM-GBSA via OpenMM + OpenFF SMIRNOFF, implicit-solvent molecular dynamics). The pipeline is applied to EGFR T790M/L858R, the clinically important gatekeeper mutant responsible for first-line resistance to gefitinib and erlotinib in non-small-cell lung cancer.

**Novelty and significance.** Most published AI-driven design approaches use opaque neural-network generators that produce SMILES strings without articulating *why*. The present work demonstrates that a large language model can perform the design step while emitting an explicit medicinal-chemistry rationale for every proposed molecule, making the chemistry behind each design auditable by the medicinal chemist. Combined with multi-method consensus scoring, the pipeline identifies a single triple-confirmed lead (EGFRm-009, a pyridopyrimidinone bearing a 2-trifluoromethyl-4-fluorobenzyl substituent) from 27 unique designed molecules across three feedback rounds.

**Reproducibility on commodity hardware.** The entire campaign was performed on a single Apple Silicon laptop without commercial software (no Schrödinger Suite, no GPU cluster), in approximately three hours of wall-clock time and at an LLM API cost of approximately $3. The complete pipeline is open-sourced as `SchrodingerLite` with the EGFR T790M demo as the canonical reproducible example. This places sophisticated computational drug-design workflows within reach of academic groups in resource-limited settings, an issue of particular relevance to clinical pharmacy and pharmaceutical sciences research in low- and middle-income countries.

**Relevance to your journal.** This work fits your journal's scope by combining methodology innovation (interpretable LLM-driven design, single-trajectory MM-GBSA with SMIRNOFF parameterization) with a concrete therapeutic case study on a clinically important kinase target. The honest treatment of limitations — particularly the absence of wet-lab validation and the noise floor of implicit-solvent MM-GBSA — is intended to give readers a clear picture of what the pipeline can and cannot demonstrate in silico.

**Author and competing interests.** This is a single-author submission. I declare no competing financial or personal interests. The work received no external funding and was performed independently. I confirm the manuscript has not been published elsewhere and is not under consideration by any other journal.

**Suggested reviewers** (the author has no professional relationship with any of these, listed alphabetically):
- [Name 1], [Affiliation], [email] — expertise in computational drug discovery, kinase inhibitors
- [Name 2], [Affiliation], [email] — expertise in MM-GBSA / free-energy methods
- [Name 3], [Affiliation], [email] — expertise in AI for drug design

Thank you for your consideration. I look forward to your editorial decision.

Sincerely,

Maytham Muthanna Muneam, FIBMS Clinical Pharmacy
