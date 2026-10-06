# Response to external reviews (draft 2)

Two informal reviews of the study pitch and draft 1. For each point: what we did and where it is. "A10" and "A11" are the amendments frozen on 2026-10-06 before the corresponding analyses (Supplementary Note 1).

## Review 2: numbered points

| # | Point | Response | Where |
|---|---|---|---|
| 1 | Conditioning finer than PP2/PM1 | Added the direct flag-based test: PP3 within PM1-like, PM1-strict and PP2-like strata | A10.1; Results Q4; Table 2 |
| 2 | No direct overlap test | **Added:** R = LR(PP3 ∧ F) / [LR(PP3) · LR(F)]. DD primary, autism replication, public cohort descriptive. Implied PP3 points given F. | A10.1; Results Q4; Figure 4 |
| 3 | Estimator not validated | Planted-signal simulation: bias and coverage at real sample size | A11.5; Supplementary Note 3 |
| 4 | Mutation-rate sensitivity | CpG-transition rates × 0.919 and × 1.088, applied to all sites | A11.1; Supplementary Note 2 |
| 5 | PM1 strength and o/e threshold | ClinVar results at PM1 +1/+2/+4 and at o/e < 0.2 (Chao et al.'s Moderate level) | A10.2; Supplementary Table 3 |
| 6 | ClinVar enrichment is weak evidence | Agreed. No longer used as evidence of overlap; the joint-LR test replaces it. | Results Q4 |
| 7 | RMC coverage artefact | Strata: genes with ≥ 2 RMC regions vs one whole-gene region | A11.3 |
| 8 | Circularity and leakage | ClinVar P/LP sites excluded (DECIPHER and HGMD are not publicly available). Leakage comparison by training labels reported. | A11.2; Results Q3 |
| 9 | Retention rests on the top-2% tail | Thresholds at 1/2/10% (stable); per-predictor intervals; predictor-class table | Results Q1; Supplementary Table 1 |
| 10 | Ceiling effect | Strata by gene constraint (s_het) and mechanism | A11.3 |
| 11 | Window size | 25/50/100-residue windows, each with a same-size random partition | A11.4; Supplementary Figure 2 |
| 12 | Transfer to all ClinVar | Labelled as extrapolation; DDG2P-monoallelic restriction reported | Methods; Limitation 5; A10.2 |
| 13 | Aim 4 overreach | "Ceiling" claim removed. Two one-sided equivalence tests added; conclusion stated as "comparable; superiority untested". | A10.3; Results Q5 |
| 14 | Vacuous AUC check | Removed | A10.2 |
| 15 | Wording on the primary test and registration | The failed DD primary is stated first. "Pre-specified and timestamped", not "registered". | Methods; Results Q3 |
| 16 | "Label-free" qualifier | "No variant-level clinical labels"; curated gene and mechanism labels are acknowledged | Methods; Limitation 2 |
| 17 | Lead with retention | Abstract and Results lead with retention; the 14,343-VUS counterfactual is secondary | Abstract; Q4 |
| 18 | Fixed-prior point | Stated, with gene-specific priors (Chen et al.) as the alternative | Discussion |
| 19 | Prior work | Chen, Harrizi, Pathak; the AlphaMissense and popEVE de novo benchmarks; the ENIGMA VCEP. **The CanVIG claim could not be verified and is not cited.** | Discussion |
| 20 | Title | It now states the finding | Title |
| 21 | Numbers | Every value states its site set and denominator. AlphaMissense 13.4 (common set) vs 15.4 (all sites) are both labelled. | Throughout |
| 22 | Structure | Results organised by question; amendment codes and hashes moved to Supplementary Note 1 | Throughout |
| 23 | Tone | Informal phrasing removed | Throughout |
| 24 | Tables and figures | Table 1 (pre-specified tests), Table 2 (overlap test), Figure 4 (overlap), Supplementary Table 1 (per-predictor) | Tables and figures |
| 25 | References | Completed and verified | References |
| 26 | Mapping | GENCODE v32 canonical (17,788/18,397 pass QC); MANE concordance not assessed; ClinVar 87.3% mapped | Methods |
| 27 | Versions | Stated for every resource | Methods |
| 28–31 | Process | Preprint, competing-preprint re-search, ClinGen-adjacent reader, containerised pipeline and Zenodo, journal rules | **Author actions** |

## Review 2: scorecard gaps

| Gap | Response |
|---|---|
| **External validation (1.5, "none")** | Autism cohort (12,316, independent) and a public multi-phenotype cohort (4,256 trios, five phenotypes; A9). The pre-registered independent test of a model awaits a large NDD cohort (SPARK). |
| **Practical impact ("overlap exists but not what to do")** | Direct overlap test with implied PP3 points given PM1/PP2 (A10.1), plus per-predictor points (Supplementary Table 1) |
| **Validity: gene discovery reuses the same data** | Stated (Limitation 2). The all-genes analyses do not depend on DDG2P. |
| **Pre-registration integrity** | Stated plainly. The autism D2 primary was declared after the DD result, before the autism analysis. |
| **Novelty of the conclusion** | Reframed: the contribution is quantifying the overlap per predictor, without labels, and testing additivity directly. |

## Review 1 (pitch)

| Point | Response |
|---|---|
| Headline vs table (2/36 at E1) | Fixed. Every count reports its level. |
| "Three independent tests" | Reworded: one replication cohort, one new-predictor extension, plus a third (public) cohort |
| Percentages vs counts | Denominators are informative predictors (19/24, 16/20, 13/14) and are now shown |
| 54 vs 55 predictors | Corrected to 55 |
| "Validated to 1e-11" | Reworded as an implementation consistency check |
| Dilution argument | **Does not apply:** the LR is excess-based. Draft 1's "diluted" wording was wrong and has been corrected. |
| "Fixed effects shrink any signal" | No fixed effects are fitted; the random-partition control shows small units remove little evidence |
| Sibling control too narrow | Extended to every level and predictor (0/37) |
| Null control / CIs | Random partition, shifted windows, InterPro domains; paired bootstrap CIs on retention |
| ClinGen cap wording | The primary source caps PP3 + PM1 at Strong (Pejaver et al. 2022) |
| Zhou trio count | Explained: the supplementary tables list only individuals with a coding de novo (13,103 affected); 12,316 after filters |

## Review 4 (draft 2)

| # | Point | Response |
|---|---|---|
| 1 | The impact analysis is null; the 14,343 counterfactual is a straw man | Agreed. The headline is now the **joint evidence**: Strong-tier PP3 + PM1-like gives 2.5–3.9 points against the cap's 4 and the uncapped 6, so the data support the cap and contradict uncapped co-use (CanVIG). Added a classification-change analysis with PS2 + PM2 (A12.4), where P vs LP can actually change. The ClinGen Evidence Repository recomputation is pre-registered (A12.5). |
| 2 | Headline rests on one cell | The abstract now reports every primary-tool R, including the public cohort (0.67 [0.42–0.95]) and the weaker REVEL/VEST4 cells (0.57–0.82). Table 2 shows all cells. |
| 3 | Overlap may be by construction; R from the null side | R decomposition into excess-side and null-side P(PP3 \| F) (A12.1). R for all predictors related to constraint and frequency inputs (A12.2, pre-registered). |
| 4 | Two adjustments; R by tier | **R by PP3 tier is reported.** Strong tier: R 0.33 [0.30–0.36] DD, 0.44 autism, 0.48 public. The joint-points estimate is now the practical adjustment; retention is descriptive. |
| 5 | Cliff edge on the Supporting claim | All points are reported on both scales (2.406 and 2.08), with intervals (A10 addendum) |
| 6 | Open credibility gaps | Leave-top-genes-out (A12.3) and cross-cohort selection, predictor resampling and headroom (A10 addendum) are pre-registered. The degenerate 13/14 interval is explained in a footnote. External validation: autism plus the public cohort. |
| 7 | Sloppiness | Fixed: public-cohort power wording ("passed the check; wide intervals"); reference header; 0.51 vs 0.56 (both all sites); v3 0.96 (cold genes, 95%) vs 0.93 (all genes, 90%) labelled; the Table 2 placeholder removed |
| — | Chao et al.'s ACMG calibration of constraint; Pejaver as the origin of the cap | Both added |
| — | Kayaalp et al. 2026 | **Not verified** in several searches; not cited until the DOI is confirmed |
| — | Move Q5 to the supplement | Done (Supplementary Note 4) |
