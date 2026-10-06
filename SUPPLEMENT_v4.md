# Supplementary material: PP3 and PM1 evidence overlap for de novo missense variants in NDD

*Draft 4, 7 October 2026. Accompanies `MANUSCRIPT_v4.md`. Every value comes from a committed output.*

---

## Supplementary Note 1. Pre-specification, amendments and deviations

The plan (`FULL_STUDY_PLAN.md`) and amendments A1–A16 are in the study repository, each committed before the analysis it governs. Snapshots of the plan at each amendment are in `prereg_fullstudy/`. The deviation log at the end of the plan lists every change, its date and its reason.

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

**Supplementary Table 5. Direct overlap test, all predictors, flags and cohorts.**


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
- **Reading:** on truth independent of computational and position evidence, the overlap replicates and the cap is the best-calibrated of the three rules.
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

[BURIAL_SUPP]

[GCLUST_SUPP]

---

## Supplementary Note 6. Pre-registered hierarchical model (A16): full results

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
| **Primary:** Δ log-loss, hierarchical − cap | +0.147 [+0.106, +0.190]: **H1 refuted** |
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
| S5: label-free vs ClinVar gene-specific PP3 points, 142 label-rich DDG2P genes | Spearman 0.01; mean absolute difference 2.39 points |
| S6: eRepo truth (not confirmatory), hierarchical − cap | +0.140 [+0.047, +0.232] |
| S7: assay truth | Supplementary Note 5 (PPM1D, SHOC2) |

**Reading.** The model encodes region-dependent PP3 evidence (Supplementary Table 6), but its evidence scale does not transfer to clinical labels across genes. With a recalibration slope of 0.27 it discriminates poorly on ClinVar, and its gene-specific PP3 evidence does not track ClinVar's. The flat rules hardly differ on V1, because only 16 V1 variants carry both PP3 and a PM1-like flag. As pre-specified, the paper does not claim the hierarchical rule.
