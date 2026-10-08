# A20 pre-registration: rules for PP3 (and PP2) beside PM1, evaluated on the de novo scale, out of cohort

Written and committed before any rule was evaluated (8 October 2026). Already seen: the cell-level de novo estimates
themselves (results_fullstudy/a10_overlap.tsv, triple_summary.txt, pm1pp2_summary.txt; reported in drafts 13-17, Tables
1-2). Not yet computed: any rule's prediction error, any fitted discount, any comparison between rules on these data.
Motivation: the A19 rule test (H1) used classification-derived laboratory truth, which A19 H4 negative controls showed is
not specific to dependence. De novo truth does not derive from the codes.

## Quantity and why it is scale-free
For a cell (predictor band b x location flag F), the observed joint evidence is ln LR(PP3_b AND F). A rule keeps a
fraction kappa of PP3's own evidence beside PM1, predicting ln LR(F) + kappa x ln LR(PP3_b). Because kappa is a fraction of
PP3's own evidence on the same scale, it is unaffected by a uniform compression of the de novo scale relative to the
clinical scale. Prediction error e = observed - predicted = ln R + (1 - kappa) ln LR(PP3_b), in points (/ ln 2.406).

## Rules (kappa by nominal PP3 strength s in {1 Supporting, 2 Moderate, 4 Strong}; PM1 nominal Moderate = 2)
R0 uncapped: 1. R1 ClinGen cap (PP3 + PM1 <= 4): min(1, 2/s) -> 1, 1, 0.5. R2 PP3 at most Supporting: 1/s -> 1, 0.5, 0.25.
R3 PP3 one strength level lower beside PM1: 0, 0.5, 0.5. R4 learned constant kappa (fitted in DD). R5 learned kappa per
band (fitted in DD). R6 no PP3 credit: 0.

## Cells, training and test
Predictors AlphaMissense (primary) and REVEL; bands "Supporting only", "Moderate only", ">= Strong"; flags PM1-like and
PM1-strict. Cells used if >= 50 de novos in the joint cell and the R interval is strictly positive. Training: DD (fits R4,
R5). Test: autism and the public cohort (independent cohorts). Weights 1/SE^2 with SE = [ln R_hi - ln R_lo]/3.92.
Metric: weighted mean squared error of e (points^2) over test cells. Uncertainty: parametric bootstrap (10,000) drawing
ln R per cell from N(ln R, SE^2), refitting R4/R5 on drawn DD values each time; cells treated as independent
(approximation; cells sharing sites are positively correlated, so intervals are optimistic).
Primary contrast H5a: best of R2, R3, R4, R5 by training-cohort (DD) MSE vs R1 (cap) on the test cohorts. A rule beats the
cap if the upper bound of [MSE(rule) - MSE(cap)] < 0. Secondary: every rule vs the cap; AlphaMissense and REVEL separately.

## PP2 beside PM1 (H5b) and the triple (H5c)
H5b: PM1-like x PP2-like, kappa_PP2 in {1 (add PP2), 0 (no PP2 beside PM1), learned in DD}; test on autism; same error.
H5c: PP3 x PM1 x PP2 joint cells (triple_summary), predicted ln LR(PM1) + kappa_PP3 ln LR(PP3) + kappa_PP2 ln LR(PP2);
combined rules: cap + PP2 added; cap + no PP2; Supporting limit + no PP2; learned (DD) PP3 and PP2 kappas; test on autism.
Descriptive (2 cohorts; no public-cohort triple estimates).

## Failure
If no rule beats the cap on the de novo scale, the cap is confirmed on truth independent of the codes, and the paper says
so. Results reported whatever they are.
