# Pre-registration A17: measurement, correction and prospective test of PP3 + PM1 double counting

*For deposit on OSF **before any analysis below is run**. Written 7 October 2026. Study repository `tanukusaitejesh-prog/ndd-pp3-calibration`; amendments A1–A16 and their results (draft 5) are already public.*

**Spine question.** Does PP3 double-count PM1, by how much, and what is the right correction?

**Inclusion rule.** An analysis enters the main text only if it passes its kill rule. Otherwise it goes to the supplement as a short negative result. Every endpoint is reported, whatever its outcome.

**Shared conventions:**
- Intervals come from a gene-clustered bootstrap: 1,000 replicates (2,000 for primaries), seed 20261008.
- Points are ln LR / ln 2.406.
- "Expert truth" means points from applied codes excluding PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5, with P ≥ 6 and B ≤ −1 (as in `make_erepo_truth.py`). The population-frequency sensitivity analysis also excludes BA1, BS1, BS2 and PM2.

## Phase 2. Measurement and truth

**T5. Expanded expert truth.**
- *Data:* the ClinGen Evidence Repository plus ClinVar 3-star (expert-panel) submissions whose comments carry codes, parsed with `a16_points`. Variants are de-duplicated on GRCh38 SNV and pooled across AlphaMissense, REVEL and VEST4 at each predictor's Supporting and Strong thresholds.
- *Endpoint:* the PP3 points earned once panel-applied PM1 applies, pooled over predictors and stratified by tier.
- *Kill rule:* main text only if the 95% interval at Strong excludes 4 points.

**T6. Explaining the clinical-scale gap.**
- *Analysis:* decompose R into the PP3 rate among B/LB variants and among P/LP variants, inside and outside PM1, on the de novo, ClinVar and expert scales. The analogue on the de novo scale is null-side versus excess-side.
- *Endpoint:* the share of the log R difference between scales attributable to the benign side. Descriptive; no kill rule; supplement unless ≥ 50% is explained by one side.

**T7. Does de novo retention predict clinical performance?**
- *Endpoints:* Spearman ρ, across informative predictors, between de novo within-window retention and (a) pooled within-gene ClinVar AUROC in DDG2P genes and (b) mean within-gene Spearman with ProteinGym assays in DDG2P genes.
- *Kill rule:* main text only if ρ > 0 with a 95% interval excluding 0, using a predictor-level bootstrap.

**T8. Replication extensions.**
- (a) Leave-one-centre-out for DD (DDD, GeneDx, Radboudumc): R and within-DDG2P R, AlphaMissense × PM1-like.
- (b) A re-scan of MaveDB and ProteinGym for DDG2P assays in which the PM1-like fraction is 5–95%; R is pooled with gene-stratified LRs.
- (c) The three-way overlap R3 = LR(PP3 ∧ PM1 ∧ PP2) / [LR(PP3) · LR(PM1) · LR(PP2)].
- *Kill rule:* (a) enters the main text if every leave-one-out upper bound is < 1. (b) and (c) are descriptive.

**T9. Double-counting index.**
- *Endpoint:* Spearman ρ across the 40 informative predictors between de novo retention r(D2) and R × PM1-like (A12.2).
- *Prediction:* ρ > 0, meaning predictors that keep less evidence within windows overlap more with PM1.
- *Kill rule:* main text only if the 95% interval excludes 0.

## Phase 3. The correction

**T10. Residual PP3, the primary correction.**
- *Score:* each variant's predictor score is converted to its percentile within its own 50-residue window (primary) or gene (secondary). The percentile is taken over the mutation-rate-weighted distribution of all possible missense SNVs in that unit, so no labels are used.
- *Thresholds:* Residual PP3 applies at the percentile thresholds that reproduce the published Supporting and Strong opportunity shares genome-wide.
- *Tests, all three needed for main text:*
  1. R (residual PP3 × PM1-like) on DD de novos has a 95% interval that includes 1, with upper bound ≥ 0.85.
  2. On expert truth, uncapped residual PP3 + PM1 achieves lower mean log-loss than the cap with original PP3, with the same logistic recalibration as A16. The recalibration is fitted on ClinVar 2023-01 unanimous labels, excluding expert-truth variants. The 95% upper bound of Δ must be < 0.
  3. Within-gene AUROC of residual PP3 on expert truth is ≥ that of the original score minus 0.02.
- *Before running:* a scoop check (literature search for "within-gene score normalization" or "percentile within gene" for PP3) is recorded in the deviation log.

**T11. Region-conditional PP3 table.**
- *Estimate:* PP3 LR by score bin × region stratum (other / PM1-like / PM1-strict), from de novo within-window LRs, with monotone constraints on bins within a stratum (isotonic). Five-fold cross-validation with genes held out.
- *Comparator:* the cap, on expert truth, using the T10 recalibration protocol.
- *Kill rule:* main text only if the 95% upper bound of Δ log-loss versus the cap is < 0.

**T12. Main table.** PP3 points once PM1 applies, by predictor (AlphaMissense, REVEL, VEST4) and tier, with gene-clustered intervals, on de novo (DD, autism) and expert truth. Always reported; it is the paper's practical output.

## Phase 4. Can de novo data calibrate within-gene effect?

**T13. Within-unit conditional-likelihood model.**
- *Model:* a Poisson model with offset log(E · window rate), so only within-window variation is learned. Features: AlphaMissense, pLDDT, Cα neighbours, amino-acid change class and conservation (phyloP). Gradient boosting with ≤ 50 trees and depth ≤ 3; five-fold cross-validation with genes held out, on DD only.
- *Tests:*
  - held-out within-window Poisson deviance improves on AlphaMissense alone (95% interval of Δ excludes 0);
  - within-gene ClinVar AUROC in DDG2P genes improves;
  - within-gene assay Spearman improves.
- *Kill rule:* main text only if all three pass. A null is reported as support for the claim that de novo data cannot calibrate within-gene effect.

## Phase 5. Prospective test

**T14. Prospective ClinVar.**
- *Validation set:* VUS in ClinVar 2026-10 reclassified P/LP or B/LB by the same submitter in any release up to 2027-04, with truth from submitter codes as in A16 (free-text parser, negations skipped).
- *Rules, frozen now:* uncapped, cap, conditional (Supplementary Table 4) and residual PP3 (T10, frozen at its specification above).
- *Primary:* Δ mean log-loss, residual PP3 − cap; recalibration on ClinVar 2026-10 unanimous labels.
- *Feasibility:* at least 30 P and 15 B are needed, or the test is reported as infeasible. Results are added in revision.

## Ordering

T9 and T12 use existing outputs. T5–T8 and T10–T13 are run only after this document is deposited on OSF. T14 is analysed only after 2027-04.
