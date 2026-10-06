# Pre-registration A16: a label-free gene → region → PP3 model against flat-prior rules, on temporally held-out ClinVar truth

*For deposit on OSF before any validation data are analysed. Written 7 October 2026, 00:30–01:15 IST. Study: "PP3 and PM1 evidence overlap for de novo missense variants in NDD" (repository `tanukusaitejesh-prog/ndd-pp3-calibration`, amendments A1–A15).*

## 1. What has and has not been seen

**Seen before this document:** every result in amendments A1–A15. That includes the overlap estimates R on the de novo scale, on ClinVar 2026-10 (P/LP vs B/LB, the estimand bridge) and on ClinGen expert-panel truth (eRepo; the cap was the best-calibrated rule there). It also includes the ProteinGym MAVE result for 7 DDG2P genes.

**Feasibility checks done (counts only):**
- MaveDB has score sets for 64 DDG2P monoallelic genes, mostly 50–90-residue domain-stability fragments.
- ClinVar archive files exist for 2023-01 and 2026-10.

**Not seen:**
- No analysis of the 2023-01 → 2026-10 reclassified variants (validation set V1). Their number is unknown at writing.
- The hierarchical model has not been fitted.
- No comparison of the four rules below has been computed on any dataset.

## 2. Hypotheses

- **H1 (primary).** On V1 truth, a label-free hierarchical model (gene prior → region prior → AlphaMissense LR conditional on region), fitted on DD de novos only, gives better-calibrated probabilities than the ClinGen cap rule. It is applied to the computational, gene-level and position-level evidence only.
- **H2.** PP3 and PM1-like evidence overlap on V1 truth (R < 1).

## 3. The hierarchical model (fitted on DD de novos only; no clinical labels)

**Data:** the same site universe, mutation rates and DD de novo missense variants (all three centres) as the main study. Expected counts are E = scaled Roulette rates; O are observed de novos.

**Units:**
- gene *g* (analysed transcript);
- 50-residue window *w* within the gene;
- region stratum *s*: PM1-strict (RMC o/e < 0.2), PM1-like only (o/e 0.2–0.4 or CCR ≥ 95), or other;
- AlphaMissense bin *b*: < 0.34, 0.34–0.792 (reference), 0.792–0.906, 0.906–0.972, ≥ 0.972. The 0.34 cut-off is AlphaMissense's published likely-benign boundary (Cheng et al. 2023); the others are the ClinGen PP3 Supporting, Moderate and Strong thresholds used throughout the study (A10).

**Rate model:** O_v ~ Poisson(E_v · ρ_v) with log ρ_v = log θ_g + log ψ_s + log φ_w + log π_{b|s}.

1. **Gene level:** θ_g ~ Gamma(a, a/μ). a and μ are set by maximum negative-binomial marginal likelihood across transcripts. θ̂_g is the posterior mean (O_g + a)/(E_g + a/μ).
2. **Region stratum and predictor bins:** ψ_s and π_{b|s} are fixed effects, estimated by Poisson regression with offset log(E · θ̂_g · φ̂_w), pooled over genes.
3. **Window level:** φ_w ~ Gamma(c, c), with c set by maximum marginal likelihood. φ̂_w = (O_w + c)/(Ê_w + c), where Ê_w is the expected count in the window under steps 1–2.
4. Steps 2–3 alternate three times. Model fitting uses no validation data.

**Output LR (de novo-excess scale, as in the main study):** LR_hier(v) = max(ρ̂_v − 1, 0.01 · (Θ − 1)) / (Θ − 1), where Θ = ΣO/ΣE. The calibration table of π_{b|s}, the PP3 LR by score bin and region stratum, is a deliverable.

## 4. Comparator rules (computational, gene and position evidence only)

| Rule | Points |
|---|---|
| 1 Uncapped | AlphaMissense PP3 points (1/2/4 at 0.792/0.906/0.972) + PM1 2 if PM1-like + PP2 1 if gene missense Z ≥ 3.09. BP4 is not applied in any rule. |
| 2 Cap (ClinGen) | As rule 1, with PP3 + PM1 limited to 4 |
| 3 Conditional (retention) | As rule 1, with PP3 points replaced by the DD conditional points of Supplementary Table 4 when PM1-like (else PP2-like) applies |
| 4 Hierarchical | ln LR_hier |

**Recalibration, identical for every rule.** A logistic model, P(P/LP) = expit(α + β · score), is fitted per rule on training set T, where score = points · ln 2.406 for rules 1–3 and ln LR_hier for rule 4.

**T:** ClinVar 2023-01 missense SNVs at analysed sites, labelled P/LP or B/LB by all of their submissions, in the primary gene stratum. All V1 variants and all DD de novo sites are excluded. Because ClinVar classifications partly use PP3 and PM1, T favours rules 1–3; this bias works against H1.

## 5. Validation set V1 (temporal ClinVar, primary)

1. Submission records (SCV) present in both submission_summary 2023-01 and 2026-10, classified *Uncertain significance* in 2023-01 and P, LP, B or LB in 2026-10 (germline).
2. Autosomal single-nucleotide missense variants mapped to analysed sites with an AlphaMissense score; DD de novo sites excluded.
3. **Truth, independent of computational and position evidence.** ACMG/AMP codes are parsed from the 2026-10 SCV comment with the parser in `make_erepo_truth.py`. Points are summed over all codes except PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5. Truth is P if points ≥ 6 and B if points ≤ −1; anything else, and any SCV without parsable codes, is excluded.
4. A variant with several qualifying SCVs takes the agreed truth; variants with conflicting truths are excluded.
5. **Primary stratum:** DDG2P monoallelic genes, if V1 has ≥ 100 P and ≥ 50 B there; otherwise all genes. The rule is fixed now; the count is the first thing run after deposit.
6. **Futility:** if V1 has < 30 P or < 15 B in all genes, H1 is not tested and is reported as infeasible.

## 6. Endpoints

**Primary:** Δ = mean log-loss(rule 4) − mean log-loss(rule 2) on V1, recalibrated as in section 4. Uncertainty comes from a gene-clustered bootstrap (2,000 replicates, seed 20261007). **H1 is supported if the 95% upper bound of Δ is < 0.** It is refuted if the lower bound is > 0; otherwise the result is inconclusive.

**Secondary (reported regardless of outcome):**
- **S1.** Δ for rule 4 against rules 1 and 3.
- **S1b. PP3 side only.** The primary contrast is repeated with rule 4's < 0.34 bin merged into the reference. This checks that any advantage does not come only from benign-side evidence, which rules 1–3 do not use.
- **S2 (H2).** R for AlphaMissense (Supporting and Strong) × PM1-like on V1, with a gene-clustered bootstrap. Supported if the 95% upper bound is < 1.
- **S3. Does R predict over-calling?** In six cells (AlphaMissense Supporting/Strong × PM1-like/PM1-strict/PP2-like), the overstatement by uncapped co-use observed on V1 (uncapped points − observed ln LR(PP3 ∧ F)/ln 2.406) is compared with the overstatement predicted from de novo R (−ln R/ln 2.406). Reported as a Spearman correlation and a table; exploratory (n = 6).
- **S4.** Pooled within-gene AUROC of each rule on V1 (genes with ≥ 3 P and ≥ 3 B).
- **S5. Label-free versus ClinVar gene-specific calibration (Task 7).** In label-rich DDG2P genes (≥ 10 P/LP and ≥ 10 B/LB missense in ClinVar 2026-10), compare:
  - ClinVar gene-specific ln LR(AlphaMissense ≥ 0.792), P/LP vs B/LB;
  - the label-free equivalent, ln of the ratio of mean LR_hier at AlphaMissense ≥ 0.792 to that at < 0.792 within the gene, times β from T.

  Reported as Spearman ρ and the mean absolute difference in points. Exploratory.
- **S6.** H1 and H2 on eRepo independent truth. The cap's performance there is already known, so this is **not confirmatory**.
- **S7.** H2 on MaveDB/ProteinGym assays with within-gene PM1-like variation. Exploratory.

## 7. Reporting

- Every endpoint is reported whatever its outcome.
- Any change after deposit goes in the deviation log with its reason.
- If H1 is refuted or inconclusive, the paper reports the overlap measurement with the cap as the best available rule, and does not claim the hierarchical rule.
