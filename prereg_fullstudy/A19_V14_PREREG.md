# A19 pre-registration: rule comparison, PP2 confirmation, calibrate-vs-cap and multi-vintage replication (draft 14)

Written and committed before any of the analyses below were run (8 October 2026). Inputs already inspected before
this document: the v13 laboratory and panel results (lab_stacking_summary.txt, stacking_cv_summary.txt,
erepo_truth_summary.txt, a12_classify_summary.txt), and a re-run of make_lab_stacking.py on the freshly downloaded
ClinVar 2026-10 files that reproduced v13 exactly. No vintage other than 2026-10 has been tabulated; no rule has been
evaluated out of sample on laboratory truth; no multilevel model, PM1-rationale classification or gene-calibrated PP3
has been run. OSF deposit: not possible from the analysis environment; the git commit hash is the timestamp. If the
authors deposit later, it is a later copy.

## Data
- Laboratory truth: germline, non-expert-panel ClinVar SCVs on missense SNVs whose comments list ACMG/AMP codes
  (make_lab_dataset.py; parser of make_lab_stacking.py unchanged). Truth = points from all codes except PP3, BP4, PM1,
  PP2, BP1, BP3, PM4, PM5: P/LP >= 6, B/LB <= -1; others excluded. Vintages 2026-10 (primary), 2024-01, 2022-01, 2020-01.
- Expert truth: ClinGen Evidence Repository (erepo_truth_sites.tsv), as in v13.

## H1 (primary): which rule for PP3 beside PM1 is best calibrated out of sample?
Population: laboratory SCVs with PM1 applied and with truth (2026-10).
Prediction for each SCV: logit p = a + (PP3 points under the rule) x ln 2.406, where a = logit of the truth rate among
PM1-applied, PP3-absent SCVs in the training genes. Rules for PP3 points when PM1 applies:
 R0 uncapped (nominal 1/2/4); R1 ClinGen cap (PP3 + PM1 <= 4 points); R2 PP3 at most Supporting (1 point);
 R3 discount (nominal x k, k = earned/applied in training genes); R4 learned credit per PP3 strength (training genes);
 R5 no PP3 credit (0 points; reference).
Evaluation: 5-fold cross-validation by gene, repeated 20 times; out-of-sample mean log-loss; 95% intervals by
submitter-clustered bootstrap (1,000) of the held-out log-loss differences.
Primary contrast H1a: R4 vs R1 (learned credit vs ClinGen cap). Secondary: R3 vs R1, R2 vs R1, R0 vs R1.
Decision: R4 better than the cap if the upper 95% bound of [log-loss(R4) - log-loss(R1)] < 0. If not, the cap is
retained as the recommended rule ("no evidence a learned rule beats the cap").
Also reported (descriptive): classifications changed per rule; for ClinVar VUS (supp_data1, DDG2P monoallelic genes),
upgrades to LP allowed or withheld per rule and the disease-associated fraction implied by the held-out credit.

## H2 (primary): PP2 credit beside PM1 is below its nominal 1 point
Population: laboratory SCVs with PM1 applied and with truth (2026-10).
Model: Bayesian binomial mixed GLM (statsmodels BinomialBayesMixedGLM), truth ~ PP2_applied + PP3_applied, random
intercepts for submitter and gene. PP2 credit = coefficient / ln 2.406 points; one-sided test that credit < 1.
Decision: confirmed if the upper bound of the 95% credible interval < 1 point.
H2b (mechanism): PM1 rationale classified from the comment text by pre-specified keyword rules into constraint-based
("constraint", "missense depleted", "regional", "o/e", "z-score", "tolerat", "intolerant", "lack of benign"),
hotspot/domain-based ("hotspot", "domain", "functional domain", "critical", "residue", "region", "motif", "binding",
"active site") and other/unclassified; if both constraint and hotspot terms occur, constraint-based. Test: PP2 credit
lower under constraint-based than hotspot/domain PM1 (interaction term; descriptive if fewer than 20 B/LB in a class).
H2c (contrast): within-submitter PP2 - PP3 credit difference beside PM1 (submitter-clustered bootstrap); descriptive.
Baselines (descriptive): PP2 and PP3 credit with and without PM1.

## H3 (secondary): calibrate, don't cap
AlphaMissense PP3 flags recalibrated per gene: within each gene, the score threshold at which the ClinVar 2025-05
P/LP vs B/LB likelihood ratio reaches 2.406 (Supporting), estimated leaving out every variant in the laboratory truth set;
genes with < 10 P/LP and < 10 B/LB fall back to the published global thresholds. Test on laboratory truth: R(PP3, PM1)
with gene-calibrated flags vs published flags (gene-clustered bootstrap). Interpretation: R moves to ~1 (interval
includes 1) -> "calibrate, don't cap"; R stays < 1 -> dependence is not an artefact of global thresholds.

## H4 (secondary): replication across ClinVar vintages and evaluation eras
SCVs pooled across vintages 2020-01, 2022-01, 2024-01, 2026-10, each SCV counted once at its earliest vintage carrying
codes (accession without version). Eras by DateLastEvaluated: <= 2019, 2020-2022, 2023-2024, 2025-2026. In each era:
R(PP3, PM1) and PP2 credit beside PM1 (gene-clustered bootstrap). Replication = R < 1 (upper bound < 1) in every era
with >= 20 B/LB in the PP3 & PM1 cell; PP2 credit point estimate < 1 in every era with >= 20 B/LB in the PM1 & PP2 cell.
Negative control: Likely Pathogenic classifications whose total reaches 6-9 with the last point from PP1 or PP4
instead of PP2, matched on number of codes and gene; their credit is reported beside PP2's.

## Multiplicity
Primary family H1a and H2: Holm at alpha 0.05. Everything else secondary or descriptive, labelled as such.

## What failure looks like
If H1a fails, the paper recommends the cap and reports that learned rules do not beat it. If H2 fails, the PP2 rule
is dropped from the recommendations. Both are reported whatever the result.
