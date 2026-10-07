# Supplementary material: PP3 and PM1 evidence overlap for de novo missense variants in NDD

*Draft 10, 7 October 2026. Accompanies `MANUSCRIPT_v10.md`. Every value comes from a committed output.*

---

## Supplementary Note 1. Pre-specification, amendments and deviations

The plan (`FULL_STUDY_PLAN.md`) and amendments A1–A17 are in the study repository, each committed before the analysis it governs. Snapshots of the plan at each amendment are in `prereg_fullstudy/`. The deviation log at the end of the plan lists every change, its date and its reason.

**What each amendment changed:**

- **A1** nested conditioning.
- **A2** autism replication.
- **A3** dbNSFP 5.4a predictors.
- **A4** clinical impact and position-unit controls.
- **A5, A7, A8, A15** de novo-trained models.
- **A6** intervals and the sibling control.
- **A9** the public cohort.
- **A10** the direct overlap test.
- **A11** robustness and simulation.
- **A12** decomposition, all predictors, leave-top-genes-out, classification change and Evidence Repository.
- **A13** overlap within disease genes.
- **A14** ClinVar validation of the models.
- **A16** the hierarchical model and its temporal validation.

**A16 registration.** The A16 specification (`prereg_fullstudy/A16_HIERARCHICAL_PREREG.md`) was committed and pushed (`34bac58`) before the model was fitted and before the validation set was analysed. The study lead could not deposit it on OSF before the analysis ran, so any OSF deposit is a later copy and is not claimed to predate the analysis.

**Deviations recorded before the A16 validation ran:**

- a free-text code parser, because ClinVar comments are prose;
- within-gene renormalisation of site rates;
- a joint estimate of the stratum × bin effects.

The pre-specified count rule then selected all genes as the primary stratum.

---

## Supplementary Note 2. Robustness of retention


**Negative control.** In 4,051 unaffected siblings, no predictor showed calibrated evidence at any level (0 of 37: the 36 dbNSFP 4.1a predictors plus AlphaMissense).

**Leakage.** Retention did not differ between predictors trained on clinical labels and label-free predictors (median r(D2) 0.29 against 0.37; Mann–Whitney p = 0.48). Leakage would in any case inflate, not reduce, within-gene evidence.

**Robustness** (Supplementary Note 2; median retention after gene and window, DD, 24 informative predictors):

| Analysis | Median r(D2) | AlphaMissense within-window LR |
|---|---|---|
| Baseline | 0.33 | 2.19 |
| CpG-transition rates × 0.919 | 0.33 | 2.17 |
| CpG-transition rates × 1.088 | 0.33 | 2.17 |
| Excluding 64,573 ClinVar P/LP sites (any review status) | 0.34 | 2.07 |

- **Neither mutation-rate model error nor known pathogenic variants drive the result:** excluding every ClinVar P/LP site leaves retention unchanged.
- **It is not an artefact of constraint-region coverage.** Retention was similar in genes with ≥ 2 constraint regions (0.32) and in genes without a regional split (0.37), where the within-gene to within-window drop cannot come from constraint regions.
- **It did not depend on gene constraint** (s_het high 0.35, low 0.36).

**Where PP3 adds the most.** Within DDG2P monoallelic genes, predictors retained more evidence beyond position: median r(D2) 0.47 in LoF genes and 0.57 in non-LoF (dominant-negative and gain-of-function) genes, where fewer than a third of predictors lost most of their evidence. Part of this is mechanical: within a curated gene set the marginal LR is smaller, so less of it is gene identity. It nonetheless indicates that within established NDD genes, and particularly for non-LoF mechanisms, predictor evidence is least redundant with position.



---

## Supplementary Note 3. Estimator simulation


**Estimator behaviour** (Supplementary Note 3; 100 simulated datasets per scenario at the real sample size): the within-unit estimator recovered planted within-gene signals with relative bias ≤ 1% at λ = 1, 2 and 4, for gene, window and random units. 95% coverage ranged from 0.86 to 0.99. A purely variant-level planted signal (λ = 4: true within-gene LR 3.11) still fell to 2.48 within windows, because AlphaMissense's top scores cluster in protein regions. That is a retention of about 0.80 relative to within-gene. The observed within-window loss in real data is much larger. Across unit sizes, real windows always retained far less than random partitions with the same number of units: median r 0.27 against 0.48 at 25 residues, 0.33 against 0.52 at 50, and 0.39 against 0.53 at 100 (common site set). Retention therefore falls with finer positional resolution, while random partitioning of any size costs little.


---

## Supplementary Note 4. Models trained on the de novo excess

**Status.** v1 and v2 were pre-specified. v3 and v4 were designed after the v1/v2 results and specified in amendments A7 and A8 before training. v5 and the ClinVar validation of all three were specified in A14 and A15 before they ran.


**Methods.**
**Training.** Gradient-boosted classifiers separate real de novo missense variants from sites drawn in proportion to mutation rate. Their inputs:

- 53 predictor scores;
- amino acids;
- regional constraint;
- AlphaFold structure (v2 onward).

Gene identity and mutation-rate features are excluded.

**Versions:**

| Model | Null sites drawn from | Training de novos | Status |
|---|---|---|---|
| v1, v2 | All genes | GeneDx | Pre-specified |
| v3 | The same transcript as each real de novo | GeneDx | Post hoc (A7) |
| v4 | The same transcript as each real de novo | All three DD centres | Post hoc (A8) |
| v5 | As v4, label-free inputs only (102 features) | All three DD centres | Specified in A15 before training |

v4 was evaluated only outside DD.

**Comparison with AlphaMissense:**

- Within-gene top-2% LR, paired bootstrap.
- Equivalence assessed with two one-sided tests, margin [0.80, 1.25].

**Summary.** Gradient-boosted models trained to separate real de novo missense variants from mutation-rate-matched sites had larger marginal LRs than any published predictor, but did not exceed AlphaMissense within gene. With a gene-matched null, they were statistically equivalent to AlphaMissense within gene in DD and autism. In the public cohort the point estimates were higher (v4: 1.30, 90% CI 1.03–1.80); this was the fifth, post hoc comparison and is exploratory.

**Genome-wide null (v1, v2).** Models trained against a genome-wide mutational null had the largest marginal LRs in the study (v2: 27.3), but less within-gene evidence than AlphaMissense:

| Comparison (v2 / AlphaMissense, within gene) | Ratio [95% CI] |
|---|---|
| DD cold genes | 0.69 [0.57, 0.84] |
| Autism | 0.81 [0.68, 0.95] |

The objective rewarded gene-level information.

**Gene-matched null (v3, v4).** Matching the null on gene raised within-gene evidence by 34–51% in every evaluation set. Within gene relative to AlphaMissense:

| Model | Cohort | Ratio [95% CI] |
|---|---|---|
| v3 | DD cold genes | 0.96 [0.82, 1.14] |
| v3 | Autism | 1.03 [0.89, 1.19] |
| v3 | Public cohort | 1.24 [0.91, 1.89] |
| v4 (all DD centres; evaluated outside DD only) | Autism | 1.09 [0.94, 1.29] |
| v4 | Public cohort | 1.30 [0.98, 1.91] |

- **Equivalence** (two one-sided tests, 90% CI within [0.80, 1.25]) was shown for v3 in DD, all genes (0.93 [0.89, 0.99]; the 0.96 above is for cold genes only, with a 95% interval) and in autism (1.03 [0.91, 1.17]), and for v4 in autism (1.09 [0.96, 1.25]).
- **In the public cohort,** v4's 90% interval lay above 1 (1.30 [1.03, 1.80]). As the fifth, post hoc model comparison, this is exploratory.
- The point estimates in the two independent cohorts favour the models, but no comparison was significant. These are the fourth and fifth model comparisons, two of them post hoc.
- **External validation on ClinVar** (A14, A15; DD de novo sites removed; pooled within-gene AUROC, P/LP vs B/LB, genes with ≥ 3 of each; gene bootstrap):

| Stratum | Model | AUROC | AlphaMissense | Difference [95% CI] |
|---|---|---|---|---|
| All genes (1,458) | v3 | 0.881 | 0.965 | −0.085 [−0.106, −0.068] |
| All genes | v4 | 0.845 | 0.965 | −0.121 [−0.148, −0.098] |
| All genes | v5 | 0.755 | 0.965 | −0.211 [−0.254, −0.174] |
| DDG2P (331) | v3 | 0.923 | 0.967 | −0.044 [−0.062, −0.031] |
| DDG2P | v4 (A14 primary) | 0.905 | 0.967 | −0.062 [−0.083, −0.048] |
| DDG2P | v5 (A15 primary) | 0.879 | 0.967 | −0.088 [−0.118, −0.067] |

AlphaMissense was also the best published predictor within gene: REVEL 0.949, VEST4 0.945, ESM1b 0.944, MisFit_S 0.944, popEVE 0.941 (DDG2P). ClinVar classifications partly use computational evidence, so these AUROCs are optimistic for every predictor.

- **v5** (A15): v4 restricted to inputs built without variant-level clinical labels (102 features). Autism within-gene top-2% LR ratio against AlphaMissense: 0.96 [0.83, 1.12].
- **Conclusion:** models trained on the de novo excess are comparable to AlphaMissense on de novos and worse within gene on ClinVar. The de novo excess is a calibration target, not a training target, for within-gene variant classification.


---

## Supplementary Note 5. Overlap: full results

**Supplementary Table 5. Direct overlap test, all cells, gene-clustered 95% intervals** (`results_fullstudy/a12_gclust2.tsv`). R: marginal. f: fraction of summed log-evidence kept. Within-gene R: all genes and DDG2P monoallelic genes.

| Cohort | Predictor | Tier | Flag | R | f | Within gene, all | Within gene, DDG2P |
|---|---|---|---|---|---|---|---|
| DD | AlphaMissense | Supporting | PM1-like | 0.44 [0.41, 0.47] | 0.77 [0.76, 0.78] | 0.75 [0.71, 0.78] | 0.69 [0.65, 0.73] |
| DD | AlphaMissense | Strong | PM1-like | 0.33 [0.31, 0.35] | 0.73 [0.72, 0.75] | 0.70 [0.66, 0.74] | 0.65 [0.60, 0.69] |
| DD | AlphaMissense | Supporting | PM1-strict | 0.37 [0.35, 0.40] | 0.77 [0.75, 0.78] | 0.67 [0.63, 0.71] | 0.61 [0.57, 0.65] |
| DD | AlphaMissense | Strong | PM1-strict | 0.25 [0.23, 0.28] | 0.72 [0.70, 0.74] | 0.58 [0.54, 0.63] | 0.53 [0.48, 0.58] |
| DD | AlphaMissense | Supporting | PP2-like | 0.61 [0.57, 0.65] | 0.85 [0.83, 0.86] | 0.96 [0.93, 1.00] | 0.95 [0.92, 0.97] |
| DD | AlphaMissense | Strong | PP2-like | 0.52 [0.48, 0.56] | 0.83 [0.81, 0.84] | 0.97 [0.93, 1.02] | 0.96 [0.92, 0.99] |
| DD | REVEL | Supporting | PM1-like | 0.57 [0.53, 0.61] | 0.84 [0.82, 0.86] | 0.78 [0.73, 0.83] | 0.72 [0.68, 0.77] |
| DD | REVEL | Strong | PM1-like | 0.53 [0.48, 0.58] | 0.84 [0.82, 0.86] | 0.77 [0.71, 0.84] | 0.68 [0.63, 0.75] |
| DD | REVEL | Supporting | PM1-strict | 0.48 [0.43, 0.54] | 0.83 [0.80, 0.85] | 0.70 [0.64, 0.76] | 0.64 [0.58, 0.71] |
| DD | REVEL | Strong | PM1-strict | 0.45 [0.37, 0.54] | 0.83 [0.79, 0.87] | 0.69 [0.60, 0.79] | 0.62 [0.52, 0.74] |
| DD | REVEL | Supporting | PP2-like | 0.73 [0.67, 0.79] | 0.90 [0.88, 0.93] | 0.98 [0.93, 1.04] | 0.97 [0.92, 1.01] |
| DD | REVEL | Strong | PP2-like | 0.75 [0.67, 0.84] | 0.92 [0.90, 0.95] | 1.03 [0.96, 1.12] | 0.99 [0.93, 1.05] |
| DD | VEST4 | Supporting | PM1-like | 0.67 [0.63, 0.70] | 0.87 [0.86, 0.89] | 0.88 [0.83, 0.93] | 0.80 [0.75, 0.84] |
| DD | VEST4 | Strong | PM1-like | 0.72 [0.66, 0.78] | 0.91 [0.88, 0.93] | 0.90 [0.82, 0.97] | 0.77 [0.71, 0.84] |
| DD | VEST4 | Supporting | PM1-strict | 0.60 [0.55, 0.64] | 0.87 [0.84, 0.88] | 0.81 [0.75, 0.87] | 0.73 [0.67, 0.80] |
| DD | VEST4 | Strong | PM1-strict | 0.59 [0.50, 0.67] | 0.87 [0.84, 0.91] | 0.75 [0.66, 0.86] | 0.63 [0.53, 0.74] |
| DD | VEST4 | Supporting | PP2-like | 0.77 [0.72, 0.82] | 0.91 [0.89, 0.93] | 1.01 [0.97, 1.07] | 0.98 [0.94, 1.01] |
| DD | VEST4 | Strong | PP2-like | 0.93 [0.84, 1.01] | 0.98 [0.95, 1.00] | 1.12 [1.04, 1.20] | 1.01 [0.94, 1.07] |
| autism | AlphaMissense | Supporting | PM1-like | 0.56 [0.47, 0.66] | 0.79 [0.75, 0.84] | 0.96 [0.84, 1.10] | 0.74 [0.67, 0.82] |
| autism | AlphaMissense | Strong | PM1-like | 0.44 [0.36, 0.54] | 0.75 [0.71, 0.79] | 0.92 [0.80, 1.07] | 0.73 [0.64, 0.82] |
| autism | AlphaMissense | Supporting | PM1-strict | 0.50 [0.42, 0.61] | 0.79 [0.74, 0.84] | 0.88 [0.75, 1.03] | 0.68 [0.60, 0.77] |
| autism | AlphaMissense | Strong | PM1-strict | 0.39 [0.30, 0.50] | 0.75 [0.69, 0.80] | 0.83 [0.67, 1.02] | 0.68 [0.55, 0.83] |
| autism | AlphaMissense | Supporting | PP2-like | 0.65 [0.55, 0.76] | 0.84 [0.79, 0.89] | 1.05 [0.93, 1.18] | 0.90 [0.84, 0.97] |
| autism | AlphaMissense | Strong | PP2-like | 0.57 [0.47, 0.68] | 0.82 [0.78, 0.87] | 1.07 [0.93, 1.22] | 0.90 [0.81, 0.98] |
| autism | REVEL | Supporting | PM1-like | 0.82 [0.69, 0.99] | 0.92 [0.85, 1.00] | 1.08 [0.92, 1.30] | 0.74 [0.65, 0.84] |
| autism | REVEL | Strong | PM1-like | 0.74 [0.57, 0.97] | 0.90 [0.82, 0.99] | 1.05 [0.83, 1.45] | 0.66 [0.56, 0.80] |
| autism | REVEL | Supporting | PM1-strict | 0.73 [0.56, 0.96] | 0.89 [0.82, 0.98] | 0.95 [0.74, 1.22] | 0.66 [0.53, 0.80] |
| autism | REVEL | Strong | PM1-strict | 0.68 [0.42, 1.02] | 0.89 [0.77, 1.01] | 0.96 [0.63, 1.40] | 0.59 [0.38, 0.84] |
| autism | REVEL | Supporting | PP2-like | 0.93 [0.78, 1.13] | 0.97 [0.90, 1.06] | 1.19 [1.02, 1.42] | 0.90 [0.80, 0.99] |
| autism | REVEL | Strong | PP2-like | 0.94 [0.73, 1.21] | 0.98 [0.90, 1.07] | 1.19 [0.95, 1.53] | 0.90 [0.79, 1.02] |
| autism | VEST4 | Supporting | PM1-like | 0.77 [0.64, 0.91] | 0.89 [0.83, 0.96] | 1.06 [0.92, 1.25] | 0.86 [0.77, 0.97] |
| autism | VEST4 | Strong | PM1-like | 0.84 [0.67, 1.05] | 0.94 [0.86, 1.02] | 1.13 [0.92, 1.43] | 0.83 [0.71, 0.98] |
| autism | VEST4 | Supporting | PM1-strict | 0.79 [0.64, 0.97] | 0.92 [0.85, 0.99] | 1.13 [0.92, 1.37] | 0.89 [0.73, 1.06] |
| autism | VEST4 | Strong | PM1-strict | 0.76 [0.51, 1.08] | 0.92 [0.80, 1.03] | 1.04 [0.72, 1.43] | 0.78 [0.53, 1.05] |
| autism | VEST4 | Supporting | PP2-like | 0.80 [0.69, 0.95] | 0.91 [0.85, 0.97] | 1.09 [0.96, 1.25] | 0.97 [0.89, 1.06] |
| autism | VEST4 | Strong | PP2-like | 1.02 [0.82, 1.29] | 1.01 [0.93, 1.11] | 1.24 [1.03, 1.53] | 0.95 [0.83, 1.08] |
| public | AlphaMissense | Supporting | PM1-like | 0.67 [0.48, 0.96] | 0.83 [0.73, 0.98] | 0.96 [0.74, 1.28] | 0.75 [0.62, 0.89] |
| public | AlphaMissense | Strong | PM1-like | 0.48 [0.34, 0.70] | 0.75 [0.66, 0.85] | 0.83 [0.60, 1.13] | 0.64 [0.48, 0.79] |
| public | AlphaMissense | Supporting | PM1-strict | 0.55 [0.38, 0.83] | 0.81 [0.69, 0.92] | 0.84 [0.61, 1.20] | 0.73 [0.61, 0.95] |
| public | AlphaMissense | Strong | PM1-strict | 0.34 [0.20, 0.54] | 0.70 [0.57, 0.81] | 0.63 [0.40, 0.99] | 0.57 [0.40, 0.87] |
| public | AlphaMissense | Supporting | PP2-like | 0.78 [0.57, 1.07] | 0.89 [0.78, 1.04] | 1.00 [0.80, 1.29] | 0.93 [0.83, 1.05] |
| public | AlphaMissense | Strong | PP2-like | 0.69 [0.49, 0.97] | 0.87 [0.78, 0.99] | 1.01 [0.80, 1.30] | 0.93 [0.82, 1.06] |
| public | REVEL | Supporting | PM1-like | 0.78 [0.54, 1.15] | 0.90 [0.77, 1.07] | 0.89 [0.66, 1.25] | 0.78 [0.64, 0.94] |
| public | REVEL | Strong | PM1-like | 0.71 [0.45, 1.09] | 0.89 [0.74, 1.04] | 0.75 [0.50, 1.13] | 0.77 [0.57, 1.00] |
| public | REVEL | Supporting | PM1-strict | 0.54 [0.26, 0.88] | 0.80 [0.56, 0.95] | 0.66 [0.35, 0.99] | 0.65 [0.40, 0.90] |
| public | REVEL | Strong | PM1-strict | 0.48 [0.15, 0.94] | 0.81 [0.51, 0.98] | 0.55 [0.21, 1.06] | 0.62 [0.28, 1.04] |
| public | REVEL | Supporting | PP2-like | 1.02 [0.73, 1.43] | 1.01 [0.88, 1.19] | 1.18 [0.93, 1.59] | 0.96 [0.86, 1.07] |
| public | REVEL | Strong | PP2-like | 0.81 [0.48, 1.19] | 0.93 [0.76, 1.07] | 0.85 [0.57, 1.16] | 0.85 [0.64, 1.01] |
| public | VEST4 | Supporting | PM1-like | 0.75 [0.54, 1.04] | 0.87 [0.76, 1.02] | 1.02 [0.77, 1.36] | 0.91 [0.75, 1.11] |
| public | VEST4 | Strong | PM1-like | 0.84 [0.56, 1.26] | 0.94 [0.80, 1.11] | 1.09 [0.76, 1.58] | 0.85 [0.66, 1.10] |
| public | VEST4 | Supporting | PM1-strict | 0.63 [0.43, 0.95] | 0.85 [0.72, 0.98] | 0.96 [0.67, 1.45] | 0.86 [0.66, 1.26] |
| public | VEST4 | Strong | PM1-strict | 0.80 [0.45, 1.37] | 0.93 [0.77, 1.11] | 1.16 [0.68, 2.00] | 0.96 [0.52, 1.73] |
| public | VEST4 | Supporting | PP2-like | 0.74 [0.54, 0.97] | 0.87 [0.73, 0.99] | 0.95 [0.72, 1.22] | 0.92 [0.79, 1.05] |
| public | VEST4 | Strong | PP2-like | 0.97 [0.67, 1.41] | 0.99 [0.86, 1.15] | 1.12 [0.82, 1.55] | 0.96 [0.80, 1.15] |

**Supplementary Table 5b. The same cells with the original individual-level bootstrap and implied points** (draft 3):


**Table 2. Direct overlap test** (PP3 at the ClinGen Supporting threshold; R = 1 means additive).

| Predictor | Flag | Cohort | R [95% CI] | PP3 LR given flag | Points (2.406 / 2.08 scale) | PP3 LR without flag |
|---|---|---|---|---|---|---|
| AlphaMissense | PM1-like | DD | 0.44 [0.41, 0.48] | 1.85 | 0.70 / 0.84 | 3.67 |
| AlphaMissense | PM1-like | autism | 0.56 [0.44, 0.67] | 1.68 | 0.59 / 0.71 | 2.50 |
| AlphaMissense | PM1-like | public3 | 0.67 [0.42, 0.95] | 1.76 | 0.64 / 0.77 | 2.08 |
| AlphaMissense | PM1-strict | DD | 0.37 [0.34, 0.41] | 1.57 | 0.51 / 0.61 | 4.07 |
| AlphaMissense | PM1-strict | autism | 0.50 [0.39, 0.62] | 1.51 | 0.47 / 0.57 | 2.88 |
| AlphaMissense | PM1-strict | public3 | 0.55 [0.33, 0.83] | 1.45 | 0.42 / 0.50 | 2.45 |
| AlphaMissense | PP2-like | DD | 0.61 [0.56, 0.66] | 2.56 | 1.07 / 1.28 | 3.96 |
| AlphaMissense | PP2-like | autism | 0.65 [0.52, 0.77] | 1.96 | 0.77 / 0.92 | 2.81 |
| AlphaMissense | PP2-like | public3 | 0.78 [0.48, 1.09] | 2.03 | 0.81 / 0.96 | 2.17 |
| REVEL | PM1-like | DD | 0.57 [0.52, 0.61] | 2.43 | 1.01 / 1.21 | 3.68 |
| REVEL | PM1-like | autism | 0.82 [0.65, 1.00] | 2.00 | 0.79 / 0.95 | 1.75 |
| REVEL | PM1-like | public3 | 0.78 [0.51, 1.13] | 2.25 | 0.92 / 1.11 | 2.46 |
| REVEL | PM1-strict | DD | 0.48 [0.44, 0.52] | 2.05 | 0.82 / 0.98 | 4.15 |
| REVEL | PM1-strict | autism | 0.73 [0.55, 0.96] | 1.78 | 0.66 / 0.79 | 2.27 |
| REVEL | PM1-strict | public3 | 0.54 [0.28, 0.87] | 1.56 | 0.51 / 0.61 | 2.83 |
| REVEL | PP2-like | DD | 0.73 [0.67, 0.79] | 3.10 | 1.29 / 1.55 | 3.85 |
| REVEL | PP2-like | autism | 0.93 [0.74, 1.13] | 2.27 | 0.93 / 1.12 | 1.76 |
| REVEL | PP2-like | public3 | 1.02 [0.67, 1.48] | 2.93 | 1.23 / 1.47 | 2.09 |
| VEST4 | PM1-like | DD | 0.67 [0.62, 0.71] | 2.00 | 0.79 / 0.95 | 2.72 |
| VEST4 | PM1-like | autism | 0.77 [0.62, 0.91] | 1.69 | 0.60 / 0.71 | 1.98 |
| VEST4 | PM1-like | public3 | 0.75 [0.48, 1.02] | 1.85 | 0.70 / 0.84 | 2.35 |
| VEST4 | PM1-strict | DD | 0.60 [0.55, 0.65] | 1.78 | 0.66 / 0.79 | 2.96 |
| VEST4 | PM1-strict | autism | 0.79 [0.62, 0.98] | 1.75 | 0.64 / 0.76 | 2.12 |
| VEST4 | PM1-strict | public3 | 0.63 [0.38, 0.96] | 1.57 | 0.51 / 0.61 | 2.45 |
| VEST4 | PP2-like | DD | 0.77 [0.71, 0.83] | 2.30 | 0.95 / 1.14 | 2.86 |
| VEST4 | PP2-like | autism | 0.80 [0.66, 0.95] | 1.77 | 0.65 / 0.78 | 2.08 |
| VEST4 | PP2-like | public3 | 0.74 [0.48, 0.98] | 1.83 | 0.69 / 0.82 | 2.53 |

**Overlap within disease genes** (A13). Where PP3 is used clinically, within known disease genes, the overlap persisted:

- **Within DDG2P monoallelic genes, comparing variants within the same gene:** R for AlphaMissense × PM1-like was 0.69 [0.68, 0.70] in DD and 0.74 [0.68, 0.80] in autism.
- **Across DDG2P genes:** R was 0.53 in DD and 0.59 in autism.
- **Within gene across all genes:** R was 0.75 [0.71, 0.78] in DD but 0.96 [0.84, 1.08] in autism.
- **Reading:** part of the overlap is between genes. Across all genes in autism, the within-gene overlap cannot be distinguished from additivity. Within DDG2P disease genes, where PP3 is applied clinically, a substantial within-gene component remains in both cohorts.

**Where the overlap comes from** (A12.1). R < 1 arises mainly on the null side. AlphaMissense flags 51% of mutation-expected sites at PM1-like positions, against 19% overall, a 2.7-fold rise. Among excess de novos the rise is only from 79% to 94% (1.2-fold) in DD; autism shows the same pattern. Inside constrained regions the predictor calls most sites damaging whether or not they are disease-causing, so PP3 separates less there. REVEL and VEST4 show the same pattern.

**Every predictor, and without the strongest genes.**

- **All predictors** (A12.2; PP3 flag = each predictor's top 2% of mutational opportunity):
  - In DD, R × PM1-like had an upper bound below 1 for all 40 informative predictors (median 0.51, range 0.18–0.80). In autism the median was 0.68.
  - Overlap was no larger for predictors whose inputs or training labels use population constraint or allele frequency: median 0.56 vs 0.51, Mann–Whitney p = 0.87.
  - The overlap is therefore not an artefact of constraint-based inputs.
- **Leave-top-genes-out** (A12.3): excluding the 50 genes with the largest de novo missense excess (1,571 of 23,851 de novos) left R at 0.48 [0.42, 0.53]. Retention within windows stayed below the random partition (median 0.35 vs 0.56).

**Ordered and disordered regions** (exploratory; AlphaFold pLDDT ≥ 70 vs < 70, PP3 at Supporting):

- **The overlap is not confined to structured domains.** In DD, R for AlphaMissense × PM1-like was 0.53 [0.50, 0.57] at pLDDT ≥ 70 and 0.39 [0.28, 0.51] at pLDDT < 70. REVEL gave 0.68 and 0.51; VEST4 0.78 and 0.75.
- **Autism:** 0.54 [0.43, 0.64] in ordered regions; in disordered regions 0.71 [0.38, 1.41] (74 de novos with PP3 and the flag), too few to separate from additivity.
- **Reading:** in disordered regions PP3 is rarer but carries a larger LR (AlphaMissense 9.0 against 3.1 in DD), but at PM1-like sites its LR halves (3.47, 1.42 points, against 6.86, 2.19 points, without the flag).

**Strong tier.** R with PP3 at its Strong threshold and a PM1-like flag:

| Predictor | DD | Autism | Public cohort |
|---|---|---|---|
| AlphaMissense | 0.33 [0.30, 0.36] | 0.44 | 0.48 |
| REVEL | 0.53 | 0.74 | 0.71 |
| VEST4 | 0.72 | 0.84 | 0.84 |

Because Strong-tier LRs are large, these R values still keep most of the log-evidence. The fraction kept, f, and the absolute joint points:

| Predictor | f, DD | f, autism | f, public | Joint points, DD / replication (2.406 scale) |
|---|---|---|---|---|
| AlphaMissense | 0.73 [0.71, 0.76] | 0.75 [0.67, 0.81] | 0.75 [0.57, 0.89] | 3.5 / 2.5–2.8 |
| REVEL | 0.84 [0.82, 0.87] | 0.90 [0.80, 0.99] | 0.89 [0.71, 1.04] | 3.9 / 3.0–3.2 |
| VEST4 | 0.91 [0.88, 0.93] | 0.94 [0.85, 1.02] | 0.94 [0.73, 1.09] | 3.6 / 2.9 |

f = 1 is additivity; ClinGen's cap corresponds to f = 0.67. In relative terms the overlap at the Strong tier is smaller than the reduction the cap imposes.

**Classification change with de novo evidence** (A12.4, DD; PS2 + PM2 + PM1-like + AlphaMissense PP3, Tavtigian points):

- **Uncapped:** 2,438 de novo missense variants reach Pathogenic (≥ 10 points), but the de novo excess implies only 84% of variants in that class are disease-associated, below what a Pathogenic classification requires.
- **Capped or conditional:** both rules move all 2,438 to Likely Pathogenic.
- **Caveat:** this is a class-level estimate on the de novo scale, not a per-variant posterior.

**The same overlap on the clinical scale** (estimand bridge; exploratory):

- **Method:** R recomputed with LRs from ClinVar P/LP against B/LB, with DD de novo sites removed and a gene-clustered bootstrap.
- **Result:** overlap was stronger than on the de novo scale.
  - DDG2P monoallelic genes: R 0.32 [0.28, 0.37] for AlphaMissense Supporting × PM1-like, 0.23 [0.18, 0.28] at Strong.
  - All genes: R 0.20 and 0.15.
- **Bias works against the finding:** curators who apply both PP3 and PM1 enrich P/LP variants in the joint category, which pushes R towards 1.
- **Reading:** the overlap is not specific to the de novo estimand.

**Independent truth: expert-panel evidence without computational or position criteria** (exploratory):

- **Truth.** ClinGen expert-panel variants were reclassified from their applied criteria with PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5 removed. P/LP means ≥ 6 remaining points, B/LB means ≤ −1. This gave 2,516 variants (1,627 P/LP, 889 B/LB) in 108 genes. The truth then depends on neither the predictor nor position.
- **Overlap with the panel's own PM1.** R for AlphaMissense was 0.47 [0.22, 1.18] at Supporting and 0.17 [0.08, 0.61] at Strong (gene-clustered bootstrap).
- **Points PP3 actually earns once PM1 applies:** 1.68 at Supporting and 1.83 at Strong. Without PM1, Strong-tier PP3 had LR 30.8.
- **Which rule fits.** The cap leaves at most 2 points for Strong PP3 next to Moderate PM1, which matches the observed 1.83. Uncapped co-use assigns 4, and the conditional rule derived from de novos assigns 1.1.
- **Reading:** on truth independent of computational and position evidence, the overlap replicates, and the PP3 points kept beside PM1 are consistent with the cap and exclude uncapped co-use.
- **Limitations:** with the o/e-based PM1-like flag the estimates are imprecise (Supporting R 0.90, Strong 0.22, upper bounds unbounded), because few B/LB variants fall in constrained regions. Panel variants are concentrated in a few genes.

**Functional truth from multiplexed assays** (exploratory; ProteinGym v1.3, DDG2P genes; `make_mave_truth.py`):

- **Truth.** A variant counted as damaging if its assay score fell below the ProteinGym cutoff (`DMS_score_bin` = 0). Only single substitutions reachable by an SNV were used, mapped to the analysed transcript. Eight assays in seven genes gave 8,907 variants.
- **Only two genes can test overlap.** A within-gene R needs the PM1-like flag to vary inside the gene. It did so only in PPM1D (73% of variants flagged) and SHOC2 (19%). In DLG4, PITX2 and PTEN nearly every assayed residue is PM1-like; in RAF1 and YAP1 none is.
- **Result.** PPM1D: R 0.72 [0.65, 0.78] at Supporting (0.71 at Strong), close to the de novo within-gene estimate in DDG2P genes (0.69). SHOC2: 0.86 [0.71, 1.04]. Pooled over both: 0.91 [0.82, 1.02] at Supporting and 0.85 [0.74, 0.98] at Strong (position-clustered bootstrap).
- **Reading:** on functional truth the overlap is in the same direction but smaller than on the de novo or clinical scale, and it is clear in one gene of two. One possible reason: the de novo overlap is mainly null-side (A12.1), weighted by mutation rate, whereas an assay weights all substitutions equally.
- **Limitations:** two informative genes; most cutoffs are median splits, which define assay-relative damage rather than pathogenicity; 41% of PTEN residues are not in the analysed site set.

**Expert-panel classifications** (A12.5; ClinGen Evidence Repository, 13,293 classifications):

- PP3 and PM1 were both applied in 1,018, mostly at PP3 Supporting (HNF1A, RUNX1, RYR1, HNF4A and FBN1 lead).
- **The cap** changed 24 classifications (12 P→LP, 12 LP→VUS), because most panels already keep PP3 + PM1 within 4 points.
- **The conditional rule** changed 238 (98 P→LP, 139 LP→VUS).
- **Caveat:** this rule extrapolates de novo NDD estimates to other diseases, where panel-specific PM1 definitions apply, so it is an upper bound on impact.

**Per-predictor guidance** (Supplementary Table 1). Multiplying 4 published points by each predictor's retention gives the points a predictor can justify when combined with position-based criteria:

| Predictor | Points |
|---|---|
| VEST4, MutPred, PROVEAN | About 2 |
| AlphaMissense, REVEL, and most others | About 1 |
| MPC, which is built from regional constraint | About 0 |

**ClinVar.** In 1,987,084 ClinVar missense variants:

- 30% of PP3-positive (AlphaMissense) P/LP variants lie at PM1-like sites.
- This enrichment is expected, partly because classifications use PP3 and PM1 themselves, so it is not used as evidence of overlap.
- If PP3 and PM1 were combined without ClinGen's cap, PP3 + PM1 alone would reach Likely Pathogenic for 14,343 VUS (AlphaMissense). The cap prevents this.
- **The overlap is concentrated in NDD genes.** In DDG2P monoallelic genes, 37.3% of AlphaMissense PP3-positive VUS lie at PM1-like sites, against 13.2% across all genes. Without the cap, PP3 + PM1 at Moderate alone would reach Likely Pathogenic for 7,101 VUS in these genes (13,245 with PM1 at Strong).
- **Stricter PM1 definition** (o/e < 0.2, Chao et al.'s Moderate level): the corresponding numbers are 9.7% and 2,148.
- **The cap and the conditional rule** both prevent PP3 + PM1 alone from reaching Likely Pathogenic in every setting (Supplementary Table 3).

**Buried and exposed sites** (exploratory; AlphaFold Cα neighbours within 10 Å, split at the median of 9; PP3 at Supporting):

| Predictor × flag | DD, buried | DD, exposed | Autism, buried | Autism, exposed |
|---|---|---|---|---|
| AlphaMissense × PM1-like | 0.55 [0.52, 0.59] | 0.36 [0.29, 0.43] | 0.57 [0.46, 0.68] | 0.58 [0.36, 0.86] |
| AlphaMissense × PP2-like | 0.64 [0.60, 0.67] | 0.61 [0.51, 0.73] | 0.64 [0.51, 0.76] | 0.69 [0.43, 1.02] |
| REVEL × PM1-like | 0.69 [0.65, 0.74] | 0.48 [0.38, 0.57] | 0.90 [0.72, 1.10] | 0.69 [0.40, 1.13] |
| VEST4 × PM1-like | 0.79 [0.74, 0.83] | 0.65 [0.54, 0.76] | 0.78 [0.64, 0.91] | 0.86 [0.53, 1.44] |

At exposed sites PP3 is rarer but carries a larger LR (AlphaMissense 7.07 against 2.88 at buried sites in DD); at PM1-like exposed sites it keeps 2.54 (1.06 points), against 5.90 elsewhere.

**Gene-clustered intervals** (AlphaMissense; transcripts resampled 1,000 times; synonymous scale factor fixed at its point estimate):

| Flag | DD Supporting | DD Strong | Autism Supporting | Autism Strong |
|---|---|---|---|---|
| PM1-like | 0.44 [0.41, 0.47] | 0.33 [0.31, 0.35] | 0.56 [0.47, 0.66] | 0.44 [0.36, 0.54] |
| PM1-strict | 0.37 [0.35, 0.40] | 0.25 [0.23, 0.28] | 0.50 [0.42, 0.61] | 0.39 [0.30, 0.50] |
| PP2-like | 0.61 [0.57, 0.65] | 0.52 [0.48, 0.56] | 0.65 [0.55, 0.76] | 0.57 [0.48, 0.68] |

---

## Supplementary Note 6. Pre-specified hierarchical model (A16): full results

**Specification:** `prereg_fullstudy/A16_HIERARCHICAL_PREREG.md` (commit `34bac58`). **Code:** `fullstudy.py hier` and `fullstudy.py a16val`. **Outputs:** `results_fullstudy/a16_*`.

**Fit (DD de novo missense only).**

- 23,851 de novos in 17,436 transcripts and 202,870 windows; Θ = 1.249.
- Gene gamma shape 2.47 (mean 1.24). Window gamma shape 2.18 after three iterations.

**Supplementary Table 6. Within-window rate multipliers by region stratum and AlphaMissense bin** (relative to the 0.34–0.792 bin in the same stratum).

| AlphaMissense bin | Other | PM1-like only | PM1-strict |
|---|---|---|---|
| < 0.34 | 0.98 | 0.79 | 0.63 |
| 0.34–0.792 (reference) | 1.00 | 1.00 | 1.00 |
| 0.792–0.906 (PP3 Supporting) | 1.07 | 1.13 | 1.54 |
| 0.906–0.972 (PP3 Moderate) | 1.11 | 1.59 | 2.37 |
| ≥ 0.972 (PP3 Strong) | 1.58 | 3.33 | 4.41 |

**Validation set V1.**

- 2,321 submissions went from VUS (2023-01) to P/LP/B/LB (2026-10) with parsable codes. 1,734 had independent truth, and 1,695 variants had agreed truth.
- At analysed missense sites: 87 P and 887 B in 546 genes. DDG2P monoallelic genes held 17 P and 108 B in 83 genes, below the pre-specified 100/50, so the primary stratum was all genes.
- Training set T (ClinVar 2023-01, unanimous): 29,623 P/LP and 46,155 B/LB.

**Results.**

| Endpoint | Result |
|---|---|
| **Primary:** Δ log-loss, hierarchical − cap | +0.147 [+0.106, +0.190]: **H1 not supported** |
| S1: hierarchical − uncapped | +0.148 [+0.107, +0.190] |
| S1: hierarchical − conditional | +0.144 [+0.101, +0.187] |
| S1b: hierarchical (PP3 side only) − cap | +0.148 [+0.107, +0.189] |
| Mean log-loss: uncapped / cap / conditional / hierarchical | 0.364 / 0.365 / 0.368 / 0.512 |
| Recalibration slope on T: uncapped / cap / conditional / hierarchical | 1.32 / 1.37 / 1.57 / 0.27 |
| S2: R, AlphaMissense Supporting × PM1-like | 0.47 [0.22, unbounded] (13 P, 3 B in the joint cell) |
| S2: R, AlphaMissense Strong × PM1-like | 0.42 [0.11, unbounded] (11 P, 1 B) |
| S2: R, AlphaMissense Supporting × PP2-like | 0.69 [0.38, 1.78] |
| S2: PM1-strict cells | No B/LB variant in the joint cell; R undefined |
| S3: overstatement observed vs predicted from de novo R (6 cells) | Spearman −0.43; uninformative because of the empty cells |
| S4: within-gene AUROC (83 pairs) | Rules 1–3: 0.825; hierarchical: 0.843; hierarchical PP3 side only: 0.735 |
| S5: label-free vs ClinVar gene-specific PP3 points, 142 label-rich DDG2P genes | Spearman 0.01; mean absolute difference 2.39 points. This is a failure in the model's home domain (NDD genes). |
| S6: eRepo truth (not confirmatory), hierarchical − cap | +0.140 [+0.047, +0.232] |
| S7: assay truth | Supplementary Note 5 (PPM1D, SHOC2) |

**Reading.** The model encodes region-dependent PP3 evidence (Supplementary Table 6), but its evidence scale does not transfer to clinical labels across genes. With a recalibration slope of 0.27 it discriminates poorly on ClinVar, and its gene-specific PP3 evidence does not track ClinVar's. The flat rules hardly differ on V1, because only 16 V1 variants carry both PP3 and a PM1-like flag. As pre-specified, the paper does not claim the hierarchical rule.

---

## Supplementary Table 7. Gene- and position-adjusted predictor benchmark

The file is `supp_table7_benchmark.tsv` (built by `make_benchmark.py` from committed outputs). It covers the 39 informative predictors in DD: 25 from dbNSFP 4.1a plus AlphaMissense, and 14 from dbNSFP 5.4a.

**Columns:**

- marginal, within-gene and within-50-residue-window top-2% LRs, with 95% intervals on the within-window LR from the A6 bootstrap;
- retention shares;
- points on the 2.406 scale;
- the inflation factor (marginal / within-window LR);
- ranks and the change in rank;
- DDG2P-stratum LRs;
- autism and public-cohort marginal and within-window LRs.

**Summary** (`results_fullstudy_backup/benchmark_summary.txt`):

- Median LR: 6.08 marginal, 2.99 within gene, 1.78 within window. Median retention 0.32, so 0.68 is gene and position.
- Median inflation 3.58 (range 1.47–10.62).
- Spearman between marginal and within-window LR: −0.18 (p = 0.28).
- Largest drops in rank: MPC −35, PrimateAI −25, MutFormer −25, PHACTboost −23, DEOGEN2 −22.
- Largest gains: MutPred +29, ClinPred +28, PROVEAN +28, MutationAssessor +27, VEST4 +24.
- Replication of within-window LR across predictors: autism Spearman 0.85 (n 25); public 0.49 (n 25). Marginal LR, DD against autism: 0.77.
- popEVE (DDD only): marginal 17.86; within gene 3.13; within 50-residue window 1.76.

---

## Supplementary Note 7. Pre-specified within-window model (A17, task 13)

**Specification.** A17 task 13, committed in `8b500f8` before the run. **Code:** `fullstudy.py t13`. **Output:** `results_fullstudy/t13_summary.txt`.

**Model.**
- Poisson histogram gradient boosting: 50 iterations, depth 3, monotone in AlphaMissense.
- Offset log(E · window rate), re-normalised once, so that only within-window variation is learned.
- Inputs: AlphaMissense, pLDDT, Cα neighbours within 10 Å, phyloP100way, reference and alternate amino acid.
- Data: 6,335,481 sites in the 20,595 windows with ≥ 1 DD de novo (23,851 de novos), with five-fold cross-validation and genes held out.
- Comparator: the same model with AlphaMissense as its only input.

**Results.**

| Test (kill rule: all three must pass) | Result | Verdict |
|---|---|---|
| Held-out within-window log-likelihood per de novo: AlphaMissense-only − mutation rate only | +0.0249 [+0.0203, +0.0300] | (reference) |
| Full model − AlphaMissense-only | +0.0008 [−0.0006, +0.0023] | fail |
| ClinVar within-gene AUROC, DDG2P (331 genes), model − AlphaMissense | 0.894 vs 0.967; −0.073 [−0.100, −0.055] | fail |
| Assays: Spearman with damage, model vs AlphaMissense | better in 1 of 8; mean 0.417 vs 0.502 | fail |

As pre-specified, the result is reported in the supplement as a negative result.

---

## Supplementary Note 8. Oracle derivation

**Setting.** Sites are pathogenic or benign; pathogenic sites carry de novo rate k × mutation rate. A window's excess X_w = E_w(ρ_w − 1)⁺ is taken from the fitted A16 rate model (DD only), so the pathogenic opportunity in w is X_w/(k − 1).

**The oracle predictor.** Its top bin holds exactly the pathogenic sites; when these are less than 2% of opportunity, the bin is padded with randomly chosen sites in proportion to opportunity.

**The estimators.** In the excess estimators used throughout, the bin holds all excess (X_b = ΣX), so:

- marginal LR = (X_b/ΣX)/(E_b/ΣE);
- within-unit LR = ΣX / Σ_u X_u E_ub/E_u.

The within-unit LR depends on how the bin's opportunity is distributed relative to each unit's excess. For a bin holding only pathogenic sites, it equals the reciprocal of the excess-weighted pathogenic share of the units.

**What k means.** k is unknown. k ≈ 30 corresponds to about 1% of all missense opportunity being DD-pathogenic, which is roughly 1,000 dominant DD genes, each with about a fifth of its missense opportunity pathogenic.

**Output** (`fullstudy.py oracle`, `results_fullstudy/oracle_summary.txt`):

```
Oracle ceiling (world defined at 50-residue-window resolution from the fitted A16 rates): windows 202870; total excess / opportunity 0.317
k    10 natural : pathogenic share 0.0352; marginal LR 28.4; within gene 5.7; within window 2.8; retention gene 0.52, window 0.30
k    10 q=2%    : pathogenic share 0.0352; marginal LR 28.4; within gene 5.7; within window 2.8; retention gene 0.52, window 0.30
k    30 natural : pathogenic share 0.0109; marginal LR 91.5; within gene 18.2; within window 8.9; retention gene 0.64, window 0.48
k    30 q=2%    : pathogenic share 0.0109; marginal LR 50.0; within gene 15.6; within window 8.2; retention gene 0.70, window 0.54
k   100 natural : pathogenic share 0.0032; marginal LR 312.3; within gene 62.2; within window 30.4; retention gene 0.72, window 0.59
k   100 q=2%    : pathogenic share 0.0032; marginal LR 50.0; within gene 30.4; within window 20.1; retention gene 0.87, window 0.77
k   300 natural : pathogenic share 0.0011; marginal LR 943.3; within gene 187.9; within window 91.8; retention gene 0.76, window 0.66
k   300 q=2%    : pathogenic share 0.0011; marginal LR 50.0; within gene 41.2; within window 33.5; retention gene 0.95, window 0.90
k  1000 natural : pathogenic share 0.0003; marginal LR 3151.8; within gene 627.8; within window 306.6; retention gene 0.80, window 0.71
k  1000 q=2%    : pathogenic share 0.0003; marginal LR 50.0; within gene 47.0; within window 43.6; retention gene 0.98, window 0.96
```

---

## Supplementary Note 9. Benchmark robustness and validation

**Sources.** `fullstudy.py bench2` (`results_fullstudy/bench2.tsv`) is the full benchmark with transcript-clustered intervals, rank intervals, 1/2/10% thresholds, the A1 and common site sets, the DDG2P stratum and a hash-based null predictor. It supersedes the individual-level Supplementary Table 7 values. `fullstudy.py benchval` (`results_fullstudy/benchval.tsv`) holds the validation.

**Benchmark summary:**

```
all_genes a1     top  1%: informative 39; median r(window) 0.33; Spearman(marginal, window) -0.11 [-0.17, +0.17]; share r<0.5 0.87 [0.69, 0.97]; NULL marginal 0.98, window 1.01
all_genes a1     top  2%: informative 39; median r(window) 0.32; Spearman(marginal, window) -0.17 [-0.26, +0.06]; share r<0.5 0.87 [0.69, 0.97]; NULL marginal 1.11, window 1.06
all_genes a1     top 10%: informative 39; median r(window) 0.29; Spearman(marginal, window) -0.22 [-0.31, -0.06]; share r<0.5 0.95 [0.77, 1.00]; NULL marginal 0.84, window 0.82
all_genes common top  1%: informative 39; median r(window) 0.34; Spearman(marginal, window) -0.02 [-0.11, +0.27]; share r<0.5 0.85 [0.67, 0.97]; NULL marginal 0.90, window 0.97
all_genes common top  2%: informative 39; median r(window) 0.32; Spearman(marginal, window) -0.15 [-0.25, +0.11]; share r<0.5 0.90 [0.69, 0.97]; NULL marginal 1.10, window 1.07
all_genes common top 10%: informative 39; median r(window) 0.30; Spearman(marginal, window) -0.17 [-0.27, -0.01]; share r<0.5 0.92 [0.79, 1.00]; NULL marginal 0.87, window 0.84
ddg2p     a1     top  1%: informative 38; median r(window) 0.44; Spearman(marginal, window) +0.55 [+0.34, +0.67]; share r<0.5 0.79 [0.50, 0.92]; NULL marginal 0.94, window 0.99
ddg2p     a1     top  2%: informative 39; median r(window) 0.43; Spearman(marginal, window) +0.68 [+0.53, +0.76]; share r<0.5 0.82 [0.62, 0.92]; NULL marginal 1.36, window 1.30
ddg2p     a1     top 10%: informative 34; median r(window) 0.39; Spearman(marginal, window) +0.68 [+0.57, +0.76]; share r<0.5 0.91 [0.71, 1.00]; NULL marginal 1.00, window 1.00
ddg2p     common top  1%: informative 38; median r(window) 0.43; Spearman(marginal, window) +0.56 [+0.33, +0.70]; share r<0.5 0.79 [0.47, 0.92]; NULL marginal 0.61, window 0.71
ddg2p     common top  2%: informative 38; median r(window) 0.43; Spearman(marginal, window) +0.67 [+0.45, +0.73]; share r<0.5 0.79 [0.50, 0.92]; NULL marginal 1.31, window 1.27
ddg2p     common top 10%: informative 34; median r(window) 0.41; Spearman(marginal, window) +0.67 [+0.58, +0.78]; share r<0.5 0.91 [0.68, 1.00]; NULL marginal 1.02, window 1.00
```

**Validation summary:**

```
Benchmark validity across 39 informative predictors (DD within-window and marginal top-2% LR, gene-clustered bench2)
Spearman(LR_marginal, clinvar_auc_pooled) = +0.28 [-0.08, +0.62] (n 39)
Spearman(LR_window, clinvar_auc_pooled) = +0.44 [+0.10, +0.68] (n 39)
Spearman(LR_marginal, clinvar_auc_within_gene) = +0.37 [+0.06, +0.63] (n 39)
Spearman(LR_window, clinvar_auc_within_gene) = +0.55 [+0.22, +0.77] (n 39)
Spearman(LR_marginal, clinvar_auc_within_gene_ddg2p) = +0.37 [+0.04, +0.64] (n 39)
Spearman(LR_window, clinvar_auc_within_gene_ddg2p) = +0.56 [+0.25, +0.77] (n 39)
Spearman(LR_marginal, assay_spearman_mean) = +0.37 [+0.04, +0.62] (n 39)
Spearman(LR_window, assay_spearman_mean) = +0.63 [+0.39, +0.79] (n 39)
window - marginal, for clinvar_auc_within_gene: +0.18 [-0.32, +0.63]
window - marginal, for clinvar_auc_within_gene_ddg2p: +0.19 [-0.29, +0.62]
window - marginal, for assay_spearman_mean: +0.26 [-0.11, +0.65]
```

ClinVar AUROCs are optimistic for predictors trained on ClinVar labels, which may inflate both correlations. The assay comparison is free of that circularity.

---

## Supplementary Note 10. PM1 × PP2 dependence and the triple stack

**Output** (`fullstudy.py pm1pp2`; transcript-clustered bootstrap, 1,000 replicates):

```
DD     PM1-like   x PP2-like: R 0.31 [0.28, 0.33]
DD     PM1-strict x PP2-like: R 0.27 [0.24, 0.31]
autism PM1-like   x PP2-like: R 0.36 [0.30, 0.45]
autism PM1-strict x PP2-like: R 0.34 [0.27, 0.44]
```

**Triple stack** (`fullstudy.py triple`; AlphaMissense PP3 × PM1 × PP2-like; transcript-clustered bootstrap, 1,000 replicates):

```
DD     AlphaMissense Supporting x PM1-like   x PP2-like: R3 0.13 [0.11, 0.14]; joint LR 24.6 (3.65 points [3.50, 3.79]) vs product 196.7 (6.02 points); n de novos in joint cell 2144
DD     AlphaMissense Strong     x PM1-like   x PP2-like: R3 0.09 [0.08, 0.11]; joint LR 35.1 (4.05 points [3.89, 4.20]) vs product 373.7 (6.75 points); n de novos in joint cell 1747
DD     AlphaMissense Supporting x PM1-strict x PP2-like: R3 0.10 [0.08, 0.12]; joint LR 38.7 (4.16 points [3.89, 4.40]) vs product 394.7 (6.81 points); n de novos in joint cell 817
DD     AlphaMissense Strong     x PM1-strict x PP2-like: R3 0.07 [0.06, 0.08]; joint LR 49.8 (4.45 points [4.17, 4.70]) vs product 749.9 (7.54 points); n de novos in joint cell 671
autism AlphaMissense Supporting x PM1-like   x PP2-like: R3 0.18 [0.13, 0.24]; joint LR 13.4 (2.96 points [2.67, 3.22]) vs product 75.5 (4.92 points); n de novos in joint cell 445
autism AlphaMissense Strong     x PM1-like   x PP2-like: R3 0.14 [0.10, 0.20]; joint LR 17.6 (3.27 points [2.96, 3.53]) vs product 124.5 (5.49 points); n de novos in joint cell 318
autism AlphaMissense Supporting x PM1-strict x PP2-like: R3 0.16 [0.11, 0.23]; joint LR 19.5 (3.38 points [2.98, 3.76]) vs product 122.6 (5.48 points); n de novos in joint cell 146
autism AlphaMissense Strong     x PM1-strict x PP2-like: R3 0.12 [0.08, 0.17]; joint LR 24.0 (3.62 points [3.19, 4.01]) vs product 202.2 (6.05 points); n de novos in joint cell 111
```

---

## Supplementary Note 11. Stacking in expert practice

**Code and data.** `make_stacking_audit.py` runs on the ClinGen Evidence Repository export of 2026-10-07 (13,293 classifications; retracted entries removed).

**Truth.** Points from the applied criteria, excluding PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5. P/LP at ≥ 6, B/LB at ≤ −1. Gene-clustered bootstrap, 1,000 replicates.

**Rules.** Class thresholds follow the Tavtigian points scale: P ≥ 10, LP 6–9, VUS 0–5, LB −1 to −6, B ≤ −7.

**Missense classifications:**

```
Classifications 6177 (retracted removed); panels 46
applied: PP3 2990, PM1 1092, PP2 592; PP3+PM1 787, PM1+PP2 174, PP3+PP2 449, all three 138
PM1+PP2 by panel: RASopathy VCEP 62; FBN1 VCEP 61; Monogenic Diabetes VCEP 22; Brain Malformations VCEP 11; PTEN VCEP 11; Pulmonary Hypertension VCEP 5; ACADVL VCEP 2
PM1+PP2 classifications with points from PP3+PM1+PP2 > 4: 35 of 174
rule [as applied]: classifications changed 0 (none)
rule [ClinGen cap only (PP3+PM1 <= 4)]: classifications changed 24 (LP->VUS 12; P->LP 12)
rule [PM1+PP2 <= PM1 (PP2 dropped when PM1 applies)]: classifications changed 39 (LP->VUS 11; P->LP 28)
rule [location cap: PP3+PM1+PP2 <= 4]: classifications changed 35 (LP->VUS 18; P->LP 17)
rule [PP3 limited to Supporting (1 point) when PM1 applies]: classifications changed 43 (LP->VUS 21; P->LP 22)
rule [PP3 credited at measured 0.6 points when PM1 applies]: classifications changed 187 (LP->VUS 102; P->LP 84; VUS->LB 1)
rule [PP3 credited at 1.10 (upper 95% bound) when PM1 applies]: classifications changed 43 (LP->VUS 21; P->LP 22)
rule [PP3 credited at 0.11 (lower 95% bound) when PM1 applies]: classifications changed 187 (LP->VUS 102; P->LP 84; VUS->LB 1)
independent expert truth: 2942 classifications (P/LP 1721, B/LB 1221, genes 130)
expert truth R(PM1,PP2) 1.07 [0.23, 1391124.13]; LR PM1 7.08, PP2 8.77, joint 66.69; n joint 95 (B/LB in joint 1)
expert truth R(PP3,PM1) 0.21 [0.12, 0.31]; LR PP3 8.24, PM1 7.08, joint 11.98; n joint 304 (B/LB in joint 17)
expert truth R(PP3,PP2) 1.18 [0.46, 2258935.42]; LR PP3 8.24, PP2 8.77, joint 85.14; n joint 242 (B/LB in joint 2)
points earned by panel-applied PP3 (any strength) once panel PM1 applies: 0.60 [0.11, 1.10]; mean PP3 points applied in that cell 1.24
expert truth R3(PP3,PM1,PP2) 79548.409 [19242.524, 177977.275]; joint LR 40674026.73 vs product 511.31; summed points 7.10 vs joint points 19.96
```

**Strength-specific credit, sensitivity analyses, gene-split cross-validation and temporal split** (`make_stacking_cv.py`):

```
Missense classifications 6177; independent truth 2942 (P/LP 1721, B/LB 1221, genes 130)
PP3 any strength           beside panel PM1: earns 0.60 points [0.14, 1.05]; applied 1.24; n 304 (B/LB 17)
PP3 Supporting (1 point)   beside panel PM1: earns 0.81 points [0.24, 1.68]; applied 1.00; n 235 (B/LB 11)
PP3 Moderate (2 points)    beside panel PM1: earns 0.02 points [-1.69, 0.89]; applied 2.00; n 67 (B/LB 6)
Gene-split CV (200 splits): training-half credit median 0.61 (IQR 0.43-0.72); held-out credit median 0.59 (IQR 0.48-0.73); held-out credit below points applied in 98% of splits
Held-out classifications changed when PP3 is credited at the training-half estimate: median 89 of 3103 (2.5-97.5%: 13-140); scaled to all classifications: median 182
Sensitivity of PP3 credit beside panel PM1 (points [95% CI]; n joint; B/LB in joint):
  truth >=6 / <=-1 (primary)     PP3 any       : 0.60 [0.08, 1.08]; n 304; B/LB 17
  truth >=6 / <=-1 (primary)     PP3 Supporting: 0.81 [0.28, 1.66]; n 235; B/LB 11
  truth >=6 / <=-1 (primary)     PP3 Moderate  : 0.02 [-1.82, 0.83]; n 67; B/LB 6
  truth >=4 / <=-1               PP3 any       : 0.61 [0.15, 1.07]; n 426; B/LB 17
  truth >=4 / <=-1               PP3 Supporting: 0.82 [0.28, 1.77]; n 329; B/LB 11
  truth >=4 / <=-1               PP3 Moderate  : 0.05 [-1.13, 0.80]; n 94; B/LB 6
  truth >=8 / <=-1               PP3 any       : 0.61 [0.13, 1.04]; n 201; B/LB 17
  truth >=8 / <=-1               PP3 Supporting: 0.89 [0.30, 1.84]; n 163; B/LB 11
  truth >=8 / <=-1               PP3 Moderate  : -0.23 [-2.61, 0.59]; n 37; B/LB 6
  truth >=6 / <=-3               PP3 any       : 0.65 [0.16, 1.30]; n 299; B/LB 12
  truth >=6 / <=-3               PP3 Supporting: 0.83 [0.30, 2.06]; n 232; B/LB 8
  truth >=6 / <=-3               PP3 Moderate  : 0.14 [-2.08, 2.34]; n 65; B/LB 4
  frequency codes also removed   PP3 any       : 0.49 [-0.06, 1.31]; n 245; B/LB 6
  frequency codes also removed   PP3 Supporting: 0.67 [-0.05, 17.93]; n 190; B/LB 4
  frequency codes also removed   PP3 Moderate  : 0.01 [-0.92, 15.34]; n 54; B/LB 2
Temporal split (approved before vs from 2025-01-01): training credit 0.70 [-0.14, 1.26] (n truth 1772); held-out credit 0.40 [-0.25, 1.07] (n truth 1170; joint 67, B/LB 9); held-out applied 1.40
  2025-2026 classifications (2507): rule [PP3 at most Supporting when PM1 applies] changes 15
  2025-2026 classifications (2507): rule [PP3 at training-period credit when PM1 applies] changes 75
  2025-2026 classifications (2507): rule [ClinGen cap] changes 16
```

**All variant types** (for reference):

```
Classifications 13286 (retracted removed); panels 46
applied: PP3 4353, PM1 1439, PP2 757; PP3+PM1 1018, PM1+PP2 208, PP3+PP2 588, all three 163
PM1+PP2 by panel: RASopathy VCEP 65; FBN1 VCEP 61; Monogenic Diabetes VCEP 47; Brain Malformations VCEP 17; PTEN VCEP 11; Pulmonary Hypertension VCEP 5; ACADVL VCEP 2
PM1+PP2 classifications with points from PP3+PM1+PP2 > 4: 35 of 208
rule [as applied]: classifications changed 0 (none)
rule [ClinGen cap only (PP3+PM1 <= 4)]: classifications changed 24 (LP->VUS 12; P->LP 12)
rule [PM1+PP2 <= PM1 (PP2 dropped when PM1 applies)]: classifications changed 45 (LB->B 1; LP->VUS 12; P->LP 32)
rule [location cap: PP3+PM1+PP2 <= 4]: classifications changed 35 (LP->VUS 18; P->LP 17)
rule [PP3 credited at measured 0.6 points when PM1 applies]: classifications changed 239 (LP->VUS 139; P->LP 99; VUS->LB 1)
rule [PP3 credited at 1.10 (upper 95% bound) when PM1 applies]: classifications changed 51 (LP->VUS 29; P->LP 22)
rule [PP3 credited at 0.11 (lower 95% bound) when PM1 applies]: classifications changed 239 (LP->VUS 139; P->LP 99; VUS->LB 1)
independent expert truth: 7563 classifications (P/LP 4929, B/LB 2634, genes 162)
expert truth R(PM1,PP2) 0.59 [0.21, 413337.22]; LR PM1 5.55, PP2 5.74, joint 18.70; n joint 108 (B/LB in joint 3)
expert truth R(PP3,PM1) 0.25 [0.16, 0.38]; LR PP3 6.62, PM1 5.55, joint 9.37; n joint 352 (B/LB in joint 19)
expert truth R(PP3,PP2) 1.92 [0.85, 1547637.84]; LR PP3 6.62, PP2 5.74, joint 73.21; n joint 276 (B/LB in joint 2)
points earned by panel-applied PP3 (any strength) once panel PM1 applies: 0.60 [0.15, 1.06]; mean PP3 points applied in that cell 1.22
expert truth R3(PP3,PM1,PP2) 76874.580 [22388.077, 186117.164]; joint LR 16230472.71 vs product 211.13; summed points 6.10 vs joint points 18.91
```

---

## Supplementary Note 12. Replication in laboratory practice

**Code.** `make_lab_stacking.py`, run on ClinVar submission_summary 2026-10 and variant_summary 2026-10.

**Data.** Germline SCVs, excluding expert-panel and practice-guideline review status, whose comments list any ACMG/AMP code, for missense SNVs.

**Parsing.** Codes and strengths were parsed from free text, with negations skipped.

**Truth.** As for the panels (PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5 excluded; P/LP ≥ 6, B/LB ≤ −1). Gene-clustered bootstrap, 1,000 replicates, seed 20261015.

**An earlier version** restricted the data to SCVs that mention PP3 or PM1. That selects on the exposure, and it biased the marginal LRs (R 1.53); the version reported here applies no such selection.

```
Lab submissions (non-expert-panel, germline, missense SNV) listing any ACMG code: 87873; submitters 472
PP3 applied 20760; PM1 applied 6795; both 4477; PP3 strength when both: 1 4212, 2 185, 4 80
independent truth: 39075 (P/LP 9487, B/LB 29588, genes 5976, submitters 291)
R(PP3, PM1) 0.18 [0.14, 0.25]
PP3 any strength beside PM1: earns 0.85 points [0.53, 1.20]; applied 1.09; n 1775; B/LB 75
PP3 Supporting   beside PM1: earns 0.83 points [0.53, 1.24]; applied 1.00; n 1681; B/LB 72
PP3 Moderate     beside PM1: earns 0.75 points [-0.37, 13.82]; applied 2.00; n 65; B/LB 3
PP3 Strong       beside PM1: earns 13.01 points [12.33, 13.54]; applied 4.00; n 29; B/LB 0
```

---

## Supplementary Note 13. Functional assays

**Code.** `make_mave_all.py`.

**Data.** All human ProteinGym v1.3 single-substitution assays. UniProt identifiers were mapped to gene symbols through the UniProt REST API (`work_pg_uniprot_genes.tsv`), and assay sequences were aligned to the analysed transcripts.

**Measures.**
- Truth: `DMS_score_bin` = 0 (damaging).
- PP3: AlphaMissense at Supporting or Strong.
- PM1-like: as in the main text.
- Within each assay, R with a position-clustered bootstrap (500 replicates). An assay is informative if 5–95% of its variants fall at PM1-like sites and it has ≥ 10 damaging and ≥ 10 non-damaging variants.
- Pooling: inverse-variance (fixed) and DerSimonian–Laird (random) effects on log R.

**Output:**

```
Human ProteinGym single-substitution assays mapped: 82 (70 genes); informative (PM1-like 5-95%, >= 10 damaging and >= 10 not): 16 assays, 13 genes; outside DDG2P: 13 assays, 10 genes
outside DDG2P (primary)  AlphaMissense Supporting x PM1-like: assays 13; fixed 0.98 [0.96, 1.00]; random 0.94 [0.88, 0.99]; share of assays with R < 1 0.77; tau2 0.007
outside DDG2P (primary)  AlphaMissense Strong     x PM1-like: assays 13; fixed 0.99 [0.96, 1.03]; random 0.98 [0.92, 1.03]; share of assays with R < 1 0.69; tau2 0.003
all genes                AlphaMissense Supporting x PM1-like: assays 16; fixed 0.96 [0.94, 0.98]; random 0.89 [0.83, 0.95]; share of assays with R < 1 0.81; tau2 0.013
all genes                AlphaMissense Strong     x PM1-like: assays 16; fixed 0.95 [0.92, 0.98]; random 0.90 [0.83, 0.99]; share of assays with R < 1 0.75; tau2 0.018
  EPHB2_HUMAN_Tsuboyama_2023_1F0M               other n   376 F 0.21 R_sup 0.40 [0.24, 0.62] R_str 0.33
  HXK4_HUMAN_Gersing_2022_activity              other n  2413 F 0.13 R_sup 0.84 [0.59, 1.50] R_str 0.76
  HXK4_HUMAN_Gersing_2023_abundance             other n  2388 F 0.13 R_sup 0.79 [0.59, 1.06] R_str 0.95
  HMDH_HUMAN_Jiang_2019                         other n  4994 F 0.46 R_sup 1.40 [1.16, 1.66] R_str 1.56
  RASH_HUMAN_Bandaru_2017                       DDG2P n   939 F 0.45 R_sup 0.74 [0.66, 0.82] R_str 0.69
  ILF3_HUMAN_Tsuboyama_2023_2L33                other n   399 F 0.57 R_sup 0.92 [0.82, 1.03] R_str 0.86
  MK01_HUMAN_Brenan_2016                        other n  2070 F 0.81 R_sup 0.97 [0.93, 1.02] R_str 0.99
  MET_HUMAN_Estevam_2023                        other n  1600 F 0.58 R_sup 0.83 [0.76, 0.92] R_str 0.85
  PPARG_HUMAN_Majithia_2016                     other n  2943 F 0.12 R_sup 0.65 [0.52, 0.81] R_str 0.69
  PPM1D_HUMAN_Miller_2022                       DDG2P n  2386 F 0.73 R_sup 0.72 [0.65, 0.79] R_str 0.71
  RBP1_HUMAN_Tsuboyama_2023_2KWH                other n   278 F 0.23 R_sup 0.60 [0.38, 0.89] R_str 0.71
  SHOC2_HUMAN_Kwon_2022                         DDG2P n  3294 F 0.19 R_sup 0.86 [0.71, 1.02] R_str 0.82
  SRC_HUMAN_Ahler_2019                          other n  1038 F 0.89 R_sup 1.00 [0.96, 1.04] R_str 1.01
  SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM   other n  1107 F 0.89 R_sup 1.00 [0.96, 1.06] R_str 1.02
  SRC_HUMAN_Nguyen_2022                         other n  1036 F 0.89 R_sup 1.02 [0.98, 1.07] R_str 1.02
  UBR5_HUMAN_Tsuboyama_2023_1I2T                other n   325 F 0.70 R_sup 0.95 [0.88, 1.03] R_str 0.91
```

---

## Supplementary Note 14. Selection into PP3, the Moderate–Supporting test and power

**Selection into PP3** (`make_erepo_eligibility.py`).
- Panel classifications were mapped to analysed sites. Truth excludes PP3, BP4, PM1, PP2, BP1, BP3, PM4, PM5, PS2 and PM6.
- PP3 eligibility is defined by score: REVEL ≥ 0.644 or AlphaMissense ≥ 0.792.

```
Mapped missense classifications 5990; with truth (PS2/PM6 also excluded) 2381 (P/LP 1492, B/LB 889, genes 103)
Within panel-PM1 classifications (n 345; P/LP 309, B/LB 36):
  PP3 applied                     : n 254 (P/LP 241, B/LB 13); LR 2.16
  eligible, PP3 not applied       : n 70 (P/LP 63, B/LB 7); LR 1.05
  not eligible, PP3 not applied   : n 21 (P/LP 5, B/LB 16); LR 0.04
credit beside panel PM1, panel-applied PP3                                   : 0.88 points [0.35, 1.49]; n 254; B/LB 13
credit beside panel PM1, score-defined PP3 (REVEL >= 0.644 or AM >= 0.792)   : 0.71 points [0.35, 1.20]; n 322; B/LB 19
credit beside panel PM1, score-defined PP3 (REVEL >= 0.644)                  : 0.68 points [0.34, 1.19]; n 315; B/LB 19
credit beside panel PM1, score-defined PP3 (AM >= 0.792)                     : 1.66 points [0.98, 2.58]; n 264; B/LB 7
```

**Credit by strength, the Moderate − Supporting difference, power, threshold and frequency sensitivity, gene-split cross-validation and the temporal split** (`make_stacking_cv.py`):

```
Missense classifications 6177; independent truth 2942 (P/LP 1721, B/LB 1221, genes 130)
PP3 any strength           beside panel PM1: earns 0.60 points [0.14, 1.05]; applied 1.24; n 304 (B/LB 17)
PP3 Supporting (1 point)   beside panel PM1: earns 0.81 points [0.24, 1.68]; applied 1.00; n 235 (B/LB 11)
PP3 Moderate (2 points)    beside panel PM1: earns 0.02 points [-1.69, 0.89]; applied 2.00; n 67 (B/LB 6)
Gene-split CV (200 splits): training-half credit median 0.61 (IQR 0.43-0.72); held-out credit median 0.59 (IQR 0.48-0.73); held-out credit below points applied in 98% of splits
Held-out classifications changed when PP3 is credited at the training-half estimate: median 89 of 3103 (2.5-97.5%: 13-140); scaled to all classifications: median 182
Sensitivity of PP3 credit beside panel PM1 (points [95% CI]; n joint; B/LB in joint):
  truth >=6 / <=-1 (primary)     PP3 any       : 0.60 [0.08, 1.08]; n 304; B/LB 17
  truth >=6 / <=-1 (primary)     PP3 Supporting: 0.81 [0.28, 1.66]; n 235; B/LB 11
  truth >=6 / <=-1 (primary)     PP3 Moderate  : 0.02 [-1.82, 0.83]; n 67; B/LB 6
  truth >=4 / <=-1               PP3 any       : 0.61 [0.15, 1.07]; n 426; B/LB 17
  truth >=4 / <=-1               PP3 Supporting: 0.82 [0.28, 1.77]; n 329; B/LB 11
  truth >=4 / <=-1               PP3 Moderate  : 0.05 [-1.13, 0.80]; n 94; B/LB 6
  truth >=8 / <=-1               PP3 any       : 0.61 [0.13, 1.04]; n 201; B/LB 17
  truth >=8 / <=-1               PP3 Supporting: 0.89 [0.30, 1.84]; n 163; B/LB 11
  truth >=8 / <=-1               PP3 Moderate  : -0.23 [-2.61, 0.59]; n 37; B/LB 6
  truth >=6 / <=-3               PP3 any       : 0.65 [0.16, 1.30]; n 299; B/LB 12
  truth >=6 / <=-3               PP3 Supporting: 0.83 [0.30, 2.06]; n 232; B/LB 8
  truth >=6 / <=-3               PP3 Moderate  : 0.14 [-2.08, 2.34]; n 65; B/LB 4
  frequency codes also removed   PP3 any       : 0.49 [-0.06, 1.31]; n 245; B/LB 6
  frequency codes also removed   PP3 Supporting: 0.67 [-0.05, 17.93]; n 190; B/LB 4
  frequency codes also removed   PP3 Moderate  : 0.01 [-0.92, 15.34]; n 54; B/LB 2
Temporal split (approved before vs from 2025-01-01): training credit 0.70 [-0.14, 1.26] (n truth 1772); held-out credit 0.40 [-0.25, 1.07] (n truth 1170; joint 67, B/LB 9); held-out applied 1.40
  2025-2026 classifications (2507): rule [PP3 at most Supporting when PM1 applies] changes 15
  2025-2026 classifications (2507): rule [PP3 at training-period credit when PM1 applies] changes 75
  2025-2026 classifications (2507): rule [ClinGen cap] changes 16
Moderate - Supporting difference in PP3 credit beside panel PM1 (gene-clustered bootstrap, 2,000):
  primary truth                         : Supporting 0.81, Moderate 0.02; difference -0.79 [-2.81, +0.27], bootstrap p 0.139; n Supporting 235 (B/LB 11), Moderate 67 (B/LB 6); SE 1.55; minimum detectable difference (80% power, two-sided 5%) 4.35 points; overall credit 0.60; truth n 2942
  PS2/PM6 also excluded                 : Supporting 0.75, Moderate 0.15; difference -0.60 [-3.03, +0.45], bootstrap p 0.240; n Supporting 191 (B/LB 11), Moderate 64 (B/LB 6); SE 1.61; minimum detectable difference (80% power, two-sided 5%) 4.52 points; overall credit 0.58; truth n 2768
  PS2/PM6 and BA1/BS1/BS2/PM2 excluded  : Supporting 0.58, Moderate 0.29; difference -0.29 [-17.96, +15.10], bootstrap p 0.658; n Supporting 128 (B/LB 4), Moderate 50 (B/LB 2); SE 5.92; minimum detectable difference (80% power, two-sided 5%) 16.58 points; overall credit 0.50; truth n 1531
PP3 credit beside PM1 with PS2/PM6 excluded from truth (gene-clustered 95% CI):
  PP3 any       : 0.58 [0.08, 1.04]; applied 1.27; n 257
  PP3 Supporting: 0.75 [0.17, 1.56]; applied 1.00; n 191
  PP3 Moderate  : 0.15 [-2.10, 0.97]; applied 2.00; n 64
```

