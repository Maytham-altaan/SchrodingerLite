# Submission guide — where and how to publish

## Realistic venues for this manuscript

Ranked by best fit for a single-author, no-wet-lab, methodology+case-study paper from a clinical-pharmacy author.

### Tier A — Submit first (these are likely accepts)

| Venue | Type | APC (USD) | Review time | Why it fits |
|---|---|---|---|---|
| **ChemRxiv** | Preprint (no peer review) | $0 | < 48 h | Locks priority date; citeable; standard practice. Submit *first*, in parallel with a journal submission. |
| **F1000Research** | Open peer review, post-publication | ~$1,150 (waivers available for LMICs) | ~1–2 weeks to first publication, refereeing in public | Accepts methodology + case studies; transparent review; LMIC-friendly. **Strong fit.** |
| **PLOS ONE** | Peer-reviewed, "scientific soundness" | ~$1,805 (waivers for LMICs) | 2–4 months | Reviews on technical soundness, not novelty. Methodology papers fit. |

### Tier B — Worth trying after one revision cycle

| Venue | Type | APC (USD) | Review time | Notes |
|---|---|---|---|---|
| **Pharmaceuticals (MDPI)** | Peer-reviewed open access | ~$2,900 (often waived) | ~1 month | Fast turnaround, accepts in silico drug-design case studies. Caveat: MDPI's reputation is mixed; cite responsibly. |
| **Molecules (MDPI)** | Peer-reviewed open access | ~$2,900 | ~1 month | Same caveat as above; broader scope. |
| **Heliyon (Elsevier)** | Peer-reviewed broad-scope | ~$1,950 | 1–3 months | Accepts methodology with case study. Reasonable impact. |
| **Iraqi Journal of Pharmaceutical Sciences** | National journal | $0–low | varies | Free / very low cost. Useful if international APC is a barrier. |
| **Journal of Computer-Aided Molecular Design** (Springer) | Peer-reviewed specialty | ~$3,890 | 3–6 months | Right audience but tougher reviewers; will likely demand DUD-E benchmark and baselines (see below). |

### Tier C — Aspirational, will demand significant additional work

| Venue | What they'll demand |
|---|---|
| **J. Chem. Inf. Model. (ACS)** | DUD-E benchmark, REINVENT baseline, larger sample, longer MD |
| **RSC Digital Discovery** | Same as JCIM plus a comparison to at least one published LLM-design pipeline |
| **J. Med. Chem.** | Wet-lab IC50; not feasible without a synthesis collaborator |

## Recommended submission strategy

**Step 1 — This week:**
1. **Submit to ChemRxiv** as a preprint. No review; standard upload. Use the `preprint.md` file (export to PDF first via pandoc or Word).
2. Get the DOI; cite it in your CV immediately.

**Step 2 — Within 2 weeks:**
3. **Submit to F1000Research** in parallel. They run open peer review *after* publication, so your paper is live and citeable within ~1–2 weeks of submission. Strongest fit for the work as it currently stands.

**Step 3 — If F1000 reviewers are constructive:**
4. Respond to F1000 reviewers (their comments are public, which is good — you can show responsiveness). If they ask for the DUD-E benchmark or REINVENT comparison, that's a reasonable revision request.

**Backup path** if F1000 doesn't work or isn't preferred:
- PLOS ONE → Pharmaceuticals (MDPI) → Heliyon

## Practical checklist before any submission

- [x] Author info filled in (`preprint.md`)
- [x] References list populated (25 refs)
- [x] Figures saved at 300 dpi (current is 150 dpi — re-export Figs 1 and 2 at 300 dpi before final submission)
- [ ] Convert `preprint.md` to PDF via pandoc: `pandoc preprint.md -o preprint.pdf --pdf-engine=xelatex`
- [ ] Add ORCID iD (register free at orcid.org if you don't have one)
- [ ] Push the codebase to GitHub and add the URL to the manuscript (replace `[URL TBD]` placeholders)
- [ ] Tag a release on GitHub and archive on Zenodo for a permanent DOI (free; Zenodo auto-mints DOIs from GitHub releases)
- [ ] Write the supplementary information ZIP containing: `final_ranked_unique.csv`, `selectivity.csv`, `mmgbsa.csv`, `stability.csv`, the `_figures/` folder, and the `summary.txt`
- [ ] Cover letter customized per venue (template in `cover_letter.md`)
- [ ] If submitting to a journal with APC: apply for the LMIC waiver. Iraq is on the World Bank low-to-middle income list and qualifies for most publisher waivers (F1000, PLOS, Wellcome-Trust-funded venues, Plan-S signatories).

## What to expect from reviewers (honest preview)

Likely objections, in decreasing order of frequency:

1. **"Where is the experimental validation?"** — Address head-on in the cover letter and limitations section. Frame the paper as methodology, not drug discovery. The honest answer is in section 4.3 of the preprint already.
2. **"How does this compare to REINVENT/GENTRL?"** — You can either run REINVENT yourself (1–2 days work, free) and add the comparison in revision, or argue in the response that the interpretability dimension makes direct numerical comparison difficult (weaker but defensible).
3. **"DUD-E benchmark?"** — Run it in revision if asked. Takes ~6 hours of compute. Strengthens the paper substantially.
4. **"MM-GBSA absolute values look unphysical."** — Already documented honestly. Point to section 4.3 and the explanation in 3.3.
5. **"2 ns MD is too short."** — Re-run at 50 ns explicit solvent on EGFRm-009 only (overnight job). Add as a revision figure.
6. **"Statistical significance with n=5?"** — Either expand the consensus chain to 20+ ligands (1–2 days of compute, no new tools needed) or rewrite the relevant claims to remove implicit statistical assertions.

## Time investment per venue path

- **ChemRxiv only**: 1–2 days (format conversion + figure tweaks + upload)
- **ChemRxiv + F1000Research**: 3–4 days
- **ChemRxiv + PLOS ONE / Pharmaceuticals**: 1 week (more careful formatting; longer cover letter)
- **ChemRxiv + JCIM**: 3–4 weeks (must add DUD-E + REINVENT comparison first)

## Estimated probability of acceptance

Best-effort honest estimates for the manuscript *as currently written*:

| Venue | First-submission acceptance probability |
|---|---|
| ChemRxiv | ~100% (no peer review) |
| F1000Research | ~80% (their bar is technical soundness, which is met) |
| PLOS ONE | ~50% (some reviewers strict on lack of wet-lab) |
| Pharmaceuticals (MDPI) | ~70% (fast review, less rigorous) |
| Heliyon | ~50% |
| JCIM (without DUD-E benchmark) | ~10% |
| J. Chem. Inf. Model. (after revisions) | ~40% |
| Iraqi Journal of Pharmaceutical Sciences | ~85% (national journal, supportive) |
