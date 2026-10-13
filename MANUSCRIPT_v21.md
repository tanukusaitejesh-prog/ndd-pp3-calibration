# Computational and location evidence are not independent in missense variant classification: label-free and functional measurement, cost, and tests of rules

*Draft 21 (8–13 October 2026), formatted for Genetics in Medicine (Article: structured abstract ≤ 200 words; main text ≤ 4,000 words; ≤ 5 display items; numbered references). Supplementary material: `SUPPLEMENT_v21.md`. Builds on draft 20 (`MANUSCRIPT_v20.md`) with the same analyses and numbers.*

**Title page**

- **Authors:** [names, degrees, affiliations]
- **Corresponding author:** [name, address, email]
- **Running title:** Dependence of PP3, PM1 and PP2
- **Keywords:** ACMG/AMP; variant classification; computational evidence; PM1; de novo variants
- **Word count:** abstract 195; main text about 3,300 (excluding abstract, tables, figure legends, references and end matter)
- **Display items:** 3 figures, 2 tables

---

## Abstract

**Purpose:** ACMG/AMP classification sums computational (PP3), hotspot/domain (PM1) and gene-level (PP2) evidence as though independent; ClinGen caps PP3 with PM1 only at Strong. We measured their dependence on truth not derived from these codes and tested which rules correct it.

**Methods:** Likelihood ratios came from de novo missense excess over mutation-rate expectation in 47,630 trios and probands, and from functional assays with expert-curated PM1 residues. Candidate rules were trained in one de novo cohort and tested in two others. Practice was audited in 94,050 classifications.

**Results:** Joint PP3 and PM1-like evidence was 0.44 (95% CI 0.41–0.47) of the product of their separate evidence, replicated in two cohorts and concentrated in known disease genes. On functional assays in five genes it was 0.59 (0.44–0.80), whereas random residues gave 1.00. Stacking PP3, PM1 and PP2 over-counted 2.3 points. In independent cohorts, no rule bettered ClinGen's cap for PP3, and the cap without PP2 beside PM1 matched the joint evidence within 0.5 points, whereas adding PP2 over-counted by 1.3–2.1. The cap bound in 4.7% of laboratory PP3-with-PM1 combinations.

**Conclusion:** Computational and location evidence overlap by up to a classification step. Keeping the cap and not adding PP2 beside PM1 corrects the overlap.

---

## Introduction

Missense variants are the largest class of variants of uncertain significance. Under the ACMG/AMP guidelines,[1] their classification draws on a computational prediction (PP3), location in a mutational hotspot or critical functional domain (PM1), and membership of a gene with a low rate of benign missense variation (PP2). The Bayesian points framework adds these in log-likelihood-ratio units, which assumes that they carry independent information.[2,3]

There are reasons to doubt the assumption. Predictors such as AlphaMissense, REVEL and VEST4[4–6] learn which genes and protein regions tolerate variation, which is the information PM1 and PP2 encode. Guidance has responded inconsistently. ClinGen caps PP3 with PM1 at Strong;[7,8] UK guidance advises against counting constraint evidence in both PM1 and PP2;[9] CanVIG-UK permits PP3 and PM1 together;[10] Sherloc builds dependencies into its rules;[11] and one expert panel withholds PP2 beside PM1.[12] These positions rest on expert judgement rather than measurement.

Measuring the dependence requires truth that the classifications themselves did not produce. ClinVar labels partly reflect PP3 and PM1, and several predictors were trained on them.[13,14] Truth re-scored from curators' own classifications has a subtler problem, which we show with negative-control criteria: it records which combinations of codes curators use as well as what the codes mean. Two kinds of truth avoid both problems. The excess of de novo variants over mutation-rate expectation in affected children defines a likelihood ratio without variant-level labels, for any set of sites and within any gene. Functional assays measure molecular damage directly.

We used both to measure the dependence between PP3, PM1 and PP2, estimated what it costs, tested candidate rules in cohorts independent of those used to choose them, and audited how panels and laboratories combine the criteria.

---

## Materials and Methods

### Analysis log and pre-registration
Analyses followed a public, version-controlled analysis log, each amendment committed before the analysis it governs (Supplementary Note 1). The original primary criterion, how much predictor evidence survives conditioning on gene and protein position, was not met under its first position control. The PP3 × PM1 dependence test was specified before its results were seen. Four later sets of tests were pre-registered in documents committed before they were run: temporal reclassification (A18), laboratory tests (A19), rule tests on de novo truth (A20) and functional assays from MaveDB (A21). None was deposited in an independent registry. Other practice analyses are exploratory. Supplementary Table S0 dates every hypothesis relative to data access.

### Data
**De novo cohorts.** We used 31,058 developmental-disorder (DD) trios,[15] 12,316 autism probands with 4,051 unaffected siblings,[16] and 4,256 public trios with congenital heart disease, schizophrenia, Tourette syndrome, epileptic encephalopathy or intellectual disability.[17,18] Overlapping studies were excluded, and individuals sharing two or more identical de novo variants with another cohort were removed (64 autism probands, 8 public trios).

**Sites and annotations.** Sites were coding positions on GENCODE v32 canonical transcripts, excluding those within 3 bp of an exon boundary. Annotations were predictor scores from dbNSFP 4.1a and 5.4a[19] and AlphaMissense,[4] gene and regional missense constraint,[20,21] constrained coding regions,[22] and DDG2P monoallelic genes.[23]

**Functional assays.** We used human single-substitution assays from ProteinGym v1.3[24] (82 assays in 70 genes) and MaveDB[25] score sets for genes with curated PM1 residues (81 examined). Curated PM1 residues are the protein positions at which any ClinGen expert panel applied PM1 to a missense variant in the Evidence Repository (export of 7 October 2026): 841 residues in 67 genes.

**Classifications.** We used 6,177 missense classifications by 46 expert panels from the Evidence Repository, and 87,873 germline ClinVar laboratory submissions (October 2026) for missense SNVs from 472 submitters whose comments list ACMG/AMP codes. A parser read codes and strengths from the free text and skipped negated codes. Against the structured codes of 13,091 Evidence Repository classifications, it had precision 0.92 and recall 0.98 (≥ 0.95 and ≥ 0.99 for PP3, PM1 and PP2).

### Dependence and cost on the de novo scale
Expected de novo counts E came from Roulette mutation rates,[26] scaled to each cohort by its synonymous de novo count. The excess is X = O − E, and for a set of sites b, LR = (X_b/ΣX)/(E_b/ΣE).

PP3 was flagged at the ClinGen thresholds for AlphaMissense (primary), REVEL and VEST4.[7,27] A site was PM1-like if its regional missense observed/expected ratio was below 0.4 or its constrained-coding-region percentile at least 95. It was PM1-strict if observed/expected was below 0.2. A gene was PP2-like if its missense Z was at least 3.09. These flags capture the constraint component of PM1 and PP2, not their full definitions.

Dependence was measured as R = LR(A ∧ B)/[LR(A) · LR(B)]: R = 1 means independence, and R < 1 means summing overstates the joint evidence by −ln R/ln 2.406 points. PP3's own evidence was compared alone and beside PM1, as the share of its log-evidence kept. The de novo scale is compressed relative to the clinical scale, so comparisons stay within it. Intervals resample genes (1,000 replicates).

To estimate cost, DD de novo variants were classified with PS2 + PM2 + PM1-like + AlphaMissense PP3, with and without the cap. Each class's disease-associated fraction, 1 − E/O, was compared with its nominal posterior.

### Dependence on functional truth
R was computed within each assay with a position-clustered bootstrap and pooled by random effects. In ProteinGym, damage was the assay's binary label. In MaveDB (pre-registered A21), damage meant a score closer to the median of the score set's nonsense than its synonymous variants. In MaveDB, R was combined within genes and pooled across genes, and random residue sets of the same size as the curated set (200 draws per score set) served as a negative control.

### Rules tested on de novo truth (A20)
For each cell of PP3 band and PM1 flag, a rule keeps a fraction κ of PP3's own evidence beside PM1. The predicted joint evidence is ln LR(PM1) + κ ln LR(PP3), a comparison on one scale that is unaffected by its compression. The rules for PP3 were the uncapped sum, ClinGen's cap, PP3 at most Supporting, PP3 one strength lower, a learned constant κ, a learned κ per band, and no credit. For PP2 they were adding it, not adding it, or a learned fraction, with the triple evaluated in combination. Learned rules were fitted in DD; every rule was tested in autism and the public cohort, with errors weighted by inverse variance and a parametric bootstrap. The pre-registered primary contrast was the rule with the lowest DD error, other than the cap, against the cap.

### Practice
**Co-application.** Co-application and strength were tabulated in panel and laboratory classifications and by year of last evaluation. Classifications were re-scored on the points scale (P ≥ 10, LP 6–9, VUS 0–5) under alternative rules.

**Classification-derived truth.** For laboratory and panel truth, classifications were re-scored without the computational and location criteria (PP3, BP4, PM1, PP2, BP1, BP3, PM4, PM5), with P/LP at ≥ 6 remaining points and B/LB at ≤ −1. On this truth the credit a criterion earns beside PM1 is ln[LR(C ∧ PM1)/LR(PM1)]/ln 2.406.

**Pre-registered tests on laboratory data (A19).** Submissions from four ClinVar releases were used for four tests:
- **H1:** six rules for PP3 beside PM1, compared by cross-validated log-loss;
- **H2:** PP2 credit beside PM1, in a mixed model with submitter and gene random effects;
- **H3:** per-gene PP3 thresholds;
- **H4:** replication by era, with PP1 and PP4 as negative-control criteria.

Holm correction covered H1 and H2.

**Temporal test (A18).** ClinVar 2019 Likely Pathogenic variants were compared on downgrade by 2026, by whether they had both PP3-like and PM1-like evidence or one.

---

## Results

### PP3, PM1 and PP2 are not independent on the de novo scale
For AlphaMissense PP3 at Supporting with a PM1-like flag, R was 0.44 [0.41, 0.47] in DD and replicated in autism (0.56) and the public cohort (0.67) (Table 1, Figure 2a). Within DDG2P monoallelic genes it was 0.69–0.75. Across all genes the within-gene ratio was 0.75 in DD and 0.96 in the other cohorts, so much of the dependence lies between genes, and within genes it is concentrated in known disease genes. Dependence was weaker for REVEL (0.57) and VEST4 (0.67). Every one of 40 informative predictors in DD had an upper bound below 1, and the dependence held in ordered and disordered regions, at buried and exposed sites, and without the 50 genes with the largest excess (Supplementary Note 3).

Beside PM1, AlphaMissense kept 43–49% of its own evidence and REVEL 61–78%, similarly at Supporting and Strong (Supplementary Table S19). The overlap arises because predictors largely carry location. Across 39 predictors, a median 0.32 of the log-evidence of their top 2% was incremental to gene and 50-residue window, against 0.51 for a random partition of genes. Inside constrained regions AlphaMissense flags 51% of mutation-expected sites, against 19% overall (Supplementary Note 2).

### The dependence on functional truth
With PM1 defined as expert panels apply it, functional assays showed the same dependence (Table 1, Figure 2b).
- **ProteinGym.** In 7 assays in 3 genes, R was 0.49 [0.40, 0.60] and below 1 in every assay. At curated residues AlphaMissense PP3 kept a median 40% of its log-evidence for damage, close to the de novo 43–49%, because the predictor already scores most of those variants as damaging (77% PP3-positive, against 34% elsewhere).
- **MaveDB (pre-registered).** In 5 genes and 12 score sets, R was 0.59 [0.44, 0.80] at Supporting and 0.48 [0.30, 0.76] at Strong. It was below 1 in GCK, LDLR, PTEN and TP53, and F9 was degenerate. Leaving out any gene gave 0.53–0.71.
- **Negative control.** Random residue sets of the same size gave 1.00 (0.98–1.08 by gene), so the dependence is specific to the residues curators flag.
- **Sensitivity.** Excluding a dominant-negative TP53 readout with inverted truth, and three F9 score sets without non-damaging variants in the joint cell, gave 0.52 [0.47, 0.57]. This analysis was added after the pre-registered run (Supplementary Note 13).

The constraint proxy, by contrast, barely predicted damage within assays (median LR 1.10) and showed little dependence (0.94).

On classification-derived truth, R for PP3 with PM1 was also below 1, in every independent half of the laboratory data and every era from before 2020 to 2026 (0.10–0.27; H4). Pre-registered negative controls behaved the same way. PP1 and PP4, which carry no information about location, had R with PM1 of 0.04 and 0.03 (Table 1). Low R on such truth therefore reflects which combinations benign classifications contain, so we treat it as a description of practice.

### The cost of stacking
Supporting AlphaMissense PP3 beside a PM1-like flag over-counts 0.97 points in DD. The triple PP3 × PM1-like × PP2-like (R 0.13 [0.11, 0.14]) over-counts 2.32 points [2.24, 2.51] in DD and 1.95 [1.63, 2.32] in autism, the width of a classification step. With PS2 + PM2 + PM1-like + AlphaMissense PP3 and no cap, 2,438 DD de novo variants were Pathogenic, a class 84% [84, 85] disease-associated against a nominal 99% (Figure 3). Under the cap they are Likely Pathogenic, a class that is 48% disease-associated against a nominal 90%. Part of that shortfall reflects the compressed estimand: the class comparison shows direction, and the points show size. In ClinVar, PP3 with Moderate PM1 alone would make 7,101 PP3-positive VUS in DDG2P genes Likely Pathogenic without the cap (Supplementary Note 7).

### Which rules correct the overlap (A20)
Rules were trained in DD and tested in autism and the public cohort (Table 2; Supplementary Note 12).
- **PP3 beside PM1.** The uncapped sum over-counted (error relative to the cap +0.80 points², 95% interval 0.51–1.09), and no credit under-counted. The pre-registered primary contrast did not beat the cap. The rule with the lowest DD error, κ learned per band, overfitted and failed in the test cohorts (+0.42 [0.05, 0.90]). Counting PP3 one level lower beside PM1 beat the cap for AlphaMissense (−0.06 [−0.11, −0.02]) but not REVEL (+0.02 [−0.11, +0.15]). These were secondary contrasts, with intervals that treat cells as independent.
- **PP2 beside PM1.** In autism, adding PP2 over-counted the joint evidence by 1.16 points (SE 0.12), and not adding it under-counted by 0.62.
- **The triple (pre-specified as descriptive).** The cap with PP2 added, as laboratories combine them, over-counted in all four autism cells (1.31–2.09 points). The cap without PP2 beside PM1 was within 0.5 points in all four.

### How panels and laboratories combine the criteria
Panels applied PP3 beside PM1 in 787 classifications and laboratories in 4,471. PP3 was at Supporting or Moderate in 99% and 98% of these, and PM1 was Moderate in 90% of laboratory combinations. The cap, which binds only when the two exceed 4 points, therefore bound in 211 of 4,471 (4.7%). Laboratories added PP2 beside PM1 in 1,481 submissions, 1,035 with PP3 as well (Supplementary Table S20). PP3 at Moderate or Strong beside PM1 rose from 0–0.9% of laboratory PM1 submissions in each year from 2018 to 2023 to 8.5% in 2025, after calibrated PP3 strengths were introduced, and PP2 beside PM1 rose from 18–20% to 31% (both trends p < 10⁻⁸; Supplementary Note 8).

Not adding PP2 beside PM1 would change 352 laboratory classifications (0.4% of 87,873), which is 24% of those combining the two. Of their Likely Pathogenic calls, 36% would become VUS. The changes span 160 genes and 66 laboratories, led by NOTCH3, FBN1, COL4A3, CASR, POLG, CACNA1A, KCNH2 and SCN1A (Table 2; Supplementary Table S21). The cap changes 40. Where it would change a class, laboratories had reported the uncapped class in 32 of 38 (Supplementary Note 10).

### Pre-registered tests on laboratory data and later reclassification
- **H1, rule comparison.** No PP3 credit beside PM1 was worse calibrated than the cap (Δ log-loss +0.038 [0.024, 0.050]), but no rule was better calibrated: a learned credit gave −0.004 [−0.010, +0.011], and PP3 at most Supporting +0.000 [−0.001, +0.001]. The verdict was unchanged under the ACMG convention for an absent code (a deviation).
- **H3, per-gene thresholds.** Recalibrating AlphaMissense thresholds per gene left the dependence (R 0.41 [0.30, 0.61]).
- **H2, PP2 credit.** PP2 earned 0.35 points [credible interval −0.12, 0.81] beside PM1 against 1.12 alone, which passes Holm correction. On the same truth, however, PP1 and PP4 also lost most of their credit beside PM1 (5.40 to 1.34 and 5.00 to 0.76 points), so this loss is not by itself evidence of redundancy; the de novo test is the one that measures it.
- **Credit on classification-derived truth (descriptive).** Supporting PP3 beside PM1 earned 0.81 points in panels and 0.84 in laboratories. Without 18 heavily templated submitters it earned 0.28 [0.02, 0.66], so the estimate is not robust (Supplementary Table S22).
- **A18, later reclassification (null).** Of 5,547 Likely Pathogenic variants in ClinVar in 2019, those with both PP3-like and PM1-like evidence were downgraded by 2026 less often than those with one (2.5% against 6.2%; gene-stratified odds ratio 0.53 [0.22, 1.07]; minimum detectable about 3). A third had left the 2026 table, about half of them to conflicting status (Supplementary Note 14).

---

## Discussion

Computational evidence is not independent of the location evidence that ACMG/AMP credits separately. We measured this on two kinds of truth that the codes did not produce. On the de novo scale the dependence held in three cohorts and within known disease genes. On functional assays, with PM1 as curators apply it, it held in five of six genes tested, and random residues showed none. Predictors have learned where pathogenic variation falls, which is what PM1 and PP2 describe, and beside PM1 PP3 keeps about 40–50% of its evidence on both scales. The dependence is concentrated between genes and in known disease genes.

Classification-derived truth could not measure this, because negative-control criteria with no location content showed the same low R. The lesson is general. Estimates of dependence or credit on truth re-scored from curators' own classifications, including our earlier ones, need negative-control criteria, because low R can reflect which combinations benign classifications contain rather than shared information.

The overlap matters. Stacked PP3, PM1 and PP2 over-count about two points, the width of a classification step, and an uncapped Pathogenic class built on them was 84% rather than 99% disease-associated. ClinGen's cap handles most of what matters for PP3. No pre-registered rule bettered it, on laboratory data or on independent de novo truth. It leaves PP2 uncorrected. On independent truth, the cap without PP2 beside PM1 reproduced the joint evidence of all three criteria within half a point, whereas adding PP2 over-counted by one to two points. UK guidance and the ClinGen cardiomyopathy panel already advise this;[9,12] our data supply the measurement they lacked.

In practice the cap rarely acts, because PP3 is almost always Supporting or Moderate and PM1 Moderate. Stronger PP3 beside PM1 has spread since calibrated strengths were introduced. Most such combinations are Moderate beside Moderate, which sum exactly to the cap's limit and are allowed in full, although PP3 keeps only about half of its evidence there. Counting Moderate PP3 one level lower beside PM1 beat the cap for AlphaMissense but not REVEL. It would apply to 264 laboratory classifications, so it should be tested prospectively before it is adopted.

Later reclassification did not show over-crediting in practice: the pre-registered temporal test was null, with a point estimate in the other direction. Variants with both kinds of evidence may truly be more often pathogenic, and the test had little power and heavy drop-out, but it argues against reading the de novo over-count as a large burden of wrong calls already in ClinVar.

**Implications.**
1. **Keep ClinGen's cap.**
2. **Do not add PP2 when PM1 applies**, unless a specification documents that the two rest on independent evidence.
3. **Test, rather than adopt, counting Moderate PP3 one level lower beside PM1.**
4. **Report conditional evidence in calibrations.** Calibrations of new predictors and location criteria should report evidence conditional on the other criteria. The de novo method needs no variant labels and runs on any predictor or region definition with the released code.
5. **Include negative-control criteria** in analyses on classification-derived truth.

The forthcoming ACMG/AMP/CAP/ClinGen standard reportedly groups related evidence to stop it being counted twice;[28] our measurements support placing PP3, PM1 and PP2 in a shared, capped group [verify against the version 4 specification].

**Relation to previous work.** That these criteria can double-count is not new.[7–12] Calibration has been done one criterion at a time,[7,21,27] and gene- and domain-aware calibration addresses heterogeneity between genes rather than dependence between criteria.[29–31] Laboratories are known to apply the criteria inconsistently.[32,33] Our contribution is the measurement: a label-free and functional estimate of the dependence, its cost, and rules tested on independent truth (Supplementary Table S17).

**Limitations.**
- **Pre-registration.** It was in the repository, not an independent registry, and the emphasis on dependence followed the failure of the original primary criterion.
- **The de novo scale.** It reflects severe early-onset disorders and is compressed, so evidence is compared within it.
- **Proxy flags.** Our de novo PM1 and PP2 flags are constraint-based proxies.
- **Rule tests.** They are cell-level, rest on two test cohorts, and their primary contrast failed through overfitting.
- **Functional assays.** Curated PM1 was testable in six genes, its residues are where pathogenic variants had already been classified, and the MaveDB damage rule fails for dominant-negative readouts.
- **Classification-derived analyses.** They are descriptive, sensitive to templated submissions, and rest on frequency-based benign truth and a parser not yet validated by double annotation.
- **Ancestry.** The cohorts are mostly of European ancestry.

**Conclusion.** PP3, PM1 and PP2 are not independent. On de novo and functional truth, PP3 keeps about half of its evidence beside PM1, and the three stacked over-count by up to a classification step. ClinGen's cap should be kept and PP2 should not be added beside PM1. Stronger PP3 beside PM1, which is spreading below the cap, should be tested for a further limit.

---

## Tables

**Table 1. Dependence between PP3 and PM1 on each kind of truth.** R = LR(A ∧ B)/[LR(A)·LR(B)]; 1 = independent. AlphaMissense PP3 at Supporting unless stated.

| Truth | Criteria | Scope | R [95% CI] |
|---|---|---|---|
| **De novo excess (label-free)** | PP3 × PM1-like | DD, 31,058 trios | 0.44 [0.41, 0.47] |
| | | Autism, 12,316 probands | 0.56 [0.47, 0.66] |
| | | Public, 4,256 trios | 0.67 [0.48, 0.96] |
| | — within DDG2P genes | DD / autism / public | 0.69 / 0.74 / 0.75 |
| | — within genes, all genes | DD / autism / public | 0.75 / 0.96 / 0.96 |
| | PP3 × PP2-like | DD | 0.61 [0.57, 0.65] |
| | PM1-like × PP2-like | DD | 0.31 [0.28, 0.33] |
| | PP3 × PM1-like × PP2-like | DD | 0.13 [0.11, 0.14] |
| **Functional assays** | PP3 × curated PM1 | ProteinGym, 7 assays, 3 genes | 0.49 [0.40, 0.60] |
| | PP3 × curated PM1 (pre-registered) | MaveDB, 12 score sets, 5 genes | 0.59 [0.44, 0.80] |
| | PP3 × random residues (negative control) | same | 1.00 (0.98–1.08 by gene) |
| | PP3 × PM1-like proxy | ProteinGym, 13 assays, 10 genes | 0.94 [0.88, 0.99] |
| **Classification-derived (descriptive)** | Laboratory PP3 × PM1 | 39,015 submissions | 0.18 [0.14, 0.25] |
| | PP1 × PM1; PP4 × PM1 (negative controls) | 38,007 submissions | 0.04; 0.03 |
| | Panel PP3 × panel PM1 | 2,942 classifications | 0.21 [0.12, 0.31] |

Intervals: de novo, gene-clustered; assays, position-clustered and pooled by random effects. Within-gene intervals and further strata are in Supplementary Table S3. Negative-control intervals: PP1 0.02–0.22; PP4 0.02–0.12.

**Table 2. Rules for PP3 and PP2 beside PM1: accuracy on independent de novo truth and reach in practice.**

| Rule | De novo test cohorts: error relative to the cap (points², 95% interval)* | | Laboratory classifications changed (of 87,873) |
|---|---|---|---|
| | AlphaMissense | REVEL | |
| **PP3 beside PM1** | | | |
| Uncapped sum | +0.80 [0.51, 1.09] | −0.02 [−0.25, 0.23] | — (as reported) |
| ClinGen cap (PP3 + PM1 ≤ 4 points) | 0 (reference) | 0 (reference) | 40 |
| PP3 one strength level lower | −0.06 [−0.11, −0.02] | +0.02 [−0.11, 0.15] | applies to 264 |
| PP3 at most Supporting | +0.07 [−0.09, 0.22] | +0.27 [0.12, 0.42] | 90 |
| Learned constant κ | −0.06 [−0.14, 0.03] | −0.08 [−0.18, 0.02] | — |
| Learned κ per band (pre-registered primary) | +0.42 [0.05, 0.90] | −0.09 [−0.22, 0.05] | — |
| No PP3 credit | +0.51 [0.21, 0.82] | +0.86 [0.55, 1.18] | — |
| **PP2 beside PM1 (autism; error in points)†** | | | |
| PP2 added | −1.16 (SE 0.12) | | — |
| PP2 not added | +0.62 (SE 0.12) | | 352 |
| **PP3 × PM1 × PP2 (autism, 4 cells; points)†** | | | |
| Cap + PP2 added | −1.31 to −2.09 | | — |
| Cap + PP2 not added | −0.31 to +0.47 | | — |
| PP3 + PM1 + PP2 ≤ 4 points | — | | 63 |

\*Weighted mean squared error in points² relative to the cap; negative favours the rule. Intervals treat cells as independent. †Negative error = over-count. Trained in DD; tested in autism and the public cohort (the triple in autism). Full results: Supplementary Note 12; laboratory re-scoring by class: Supplementary Table S21.

---

## Figure legends

**Figure 1. Study design.** (a) Label-free dependence: de novo excess over mutation-rate expectation in three cohorts gives likelihood ratios for PP3, PM1-like and PP2-like sites, from which dependence, PP3's evidence beside PM1 and points over-counted are estimated. (b) Functional truth: dependence on ProteinGym and MaveDB assays with PM1 as curated by expert panels, against random residue sets and the constraint proxy. (c) Rules and practice: rules trained in one de novo cohort and tested in two others; co-application in expert-panel and laboratory classifications; pre-registered tests on laboratory data with negative-control criteria and on later reclassification.

**Figure 2. PP3 and PM1 are dependent on truth that does not come from the codes.** (a) De novo excess: R for AlphaMissense PP3 (Supporting) with a PM1-like flag in three cohorts, across all sites, within DDG2P monoallelic genes and within genes; gene-clustered 95% intervals. (b) Functional assays: R with PM1 residues as curated by ClinGen expert panels, per gene and pooled (MaveDB, pre-registered; ProteinGym), against random residue sets of the same size (open circles); position-clustered 95% intervals pooled by random effects. F9 is omitted (degenerate estimate). R = 1 is independence.

**Figure 3. The cost of stacking.** Disease-associated fraction (de novo excess, DD) of the classes produced by PS2 + PM2 + PM1-like + AlphaMissense PP3 without and with ClinGen's cap, with exact 95% intervals, against each class's nominal posterior.

---

## Data Availability
All data are public: published de novo tables, ClinVar, the ClinGen Evidence Repository, ProteinGym, MaveDB, dbNSFP (academic licence), AlphaMissense, gnomAD, DDG2P and Roulette. Code, the analysis log, pre-registrations and all outputs are in the study repository (https://github.com/tanukusaitejesh-prog/ndd-pp3-calibration; to be made public with a Zenodo DOI at submission). `REPRODUCE.md` maps every number, table and figure to its script.

## Acknowledgments
[To be completed.]

## Funding Statement
[To be completed.]

## Author Contributions
[CRediT roles to be completed.]

## Ethics Declaration
This study used only publicly available, de-identified, aggregate or summary-level data (published de novo variant tables, ClinVar, the ClinGen Evidence Repository, functional-assay databases and population constraint resources). No human participants were recruited and no individual-level identifiable data were accessed, so institutional review board approval and informed consent were not required. [Confirm against the authors' institutional policy.]

## Conflict of Interest
[To be completed.]

---

## References
1. Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med*. 2015;17:405–424.
2. Tavtigian SV, Greenblatt MS, Harrison SM, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med*. 2018;20:1054–1060.
3. Tavtigian SV, Harrison SM, Boucher KM, Biesecker LG. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Hum Mutat*. 2020;41:1734–1737.
4. Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*. 2023;381:eadg7492.
5. Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *Am J Hum Genet*. 2016;99:877–885.
6. Carter H, Douville C, Stenson PD, Cooper DN, Karchin R. Identifying Mendelian disease genes with the variant effect scoring tool. *BMC Genomics*. 2013;14(Suppl 3):S3.
7. Pejaver V, Byrne AB, Feng BJ, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet*. 2022;109:2163–2177.
8. Stenton SL, Pejaver V, Bergquist T, et al. Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med*. 2024;26:101213.
9. Durkie M, Cassidy E-J, Berry I, et al. ACGS best practice guidelines for variant classification in rare disease 2024. Association for Clinical Genomic Science; 2024.
10. Garrett A, Durkie M, Callaway A, et al. Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. *J Med Genet*. 2021;58:297–304.
11. Nykamp K, Anderson M, Powers M, et al. Sherloc: a comprehensive refinement of the ACMG-AMP variant classification criteria. *Genet Med*. 2017;19:1105–1117.
12. Kelly MA, Caleshu C, Morales A, et al. Adaptation and validation of the ACMG/AMP variant classification framework for MYH7-associated inherited cardiomyopathies: recommendations by ClinGen's Inherited Cardiomyopathy Expert Panel. *Genet Med*. 2018;20:351–359.
13. Grimm DG, Azencott CA, Aicheler F, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat*. 2015;36:513–523.
14. Harrizi S, Nait Irahal I, Mostafa K, Arnoult D. Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. *bioRxiv*. 2026. doi:10.64898/2026.08.21.746180.
15. Kaplanis J, Samocha KE, Wiel L, et al. Evidence for 28 genetic disorders discovered by combining healthcare and research data. *Nature*. 2020;586:757–762.
16. Zhou X, Feliciano P, Shu C, et al. Integrating de novo and inherited variants in 42,607 autism cases identifies mutations in new moderate-risk genes. *Nat Genet*. 2022;54:1305–1319.
17. Jin SC, Homsy J, Zaidi S, et al. Contribution of rare inherited and de novo variants in 2,871 congenital heart disease probands. *Nat Genet*. 2017;49:1593–1601.
18. Turner TN, Yi Q, Krumm N, et al. denovo-db: a compendium of human de novo variants. *Nucleic Acids Res*. 2017;45:D804–D811.
19. Liu X, Li C, Mou C, Dong Y, Tu Y. dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs. *Genome Med*. 2020;12:103.
20. Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature*. 2020;581:434–443.
21. Chao KR, Wang L, Panchal R, et al. The landscape of regional missense mutational intolerance quantified from 125,748 exomes. *bioRxiv*. 2024. doi:10.1101/2024.04.11.588920.
22. Havrilla JM, Pedersen BS, Layer RM, Quinlan AR. A map of constrained coding regions in the human genome. *Nat Genet*. 2019;51:88–95.
23. Thormann A, Halachev M, McLaren W, et al. Flexible and scalable diagnostic filtering of genomic variants using G2P with Ensembl VEP. *Nat Commun*. 2019;10:2373.
24. Notin P, Kollasch AW, Ritter D, et al. ProteinGym: large-scale benchmarks for protein fitness prediction and design. *Adv Neural Inf Process Syst*. 2023.
25. Rubin AF, Stone J, Bianchi AH, Capodanno BJ, Da EY, Dias M, et al. MaveDB 2024: a curated community database with over seven million variant effects from multiplexed functional assays. *Genome Biol*. 2025;26:13.
26. Seplyarskiy V, Koch EM, Lee DJ, et al. A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription. *Nat Genet*. 2023;55:2235–2242.
27. Bergquist T, Stenton SL, Nadeau EAW, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med*. 2025;27:101402.
28. ACMG/AMP/CAP/ClinGen. Standards for sequence variant classification, version 4 (draft and pilot materials). [Cite the publication or pilot report once available.]
29. Tejura M, Fayer S, McEwen AE, et al. Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes. *Am J Hum Genet*. 2024;111:2031–2043.
30. Badonyi M, Marsh JA. acmgscaler: an R package and Colab for standardized gene-level variant effect score calibration within the ACMG/AMP framework. *Bioinformatics*. 2025;41(10):btaf503.
31. Chen Y, Fayer S, Jain S, et al. Gene- and domain-aware calibration increases the clinical utility of variant effect predictors. *bioRxiv*. 2026. doi:10.1101/2026.02.17.706269.
32. Amendola LM, Jarvik GP, Leo MC, et al. Performance of ACMG-AMP variant-interpretation guidelines among nine laboratories in the Clinical Sequencing Exploratory Research Consortium. *Am J Hum Genet*. 2016;98:1067–1076.
33. Harrison SM, Dolinsky JS, Knight Johnson AE, et al. Clinical laboratories collaborate to resolve differences in variant interpretations submitted to ClinVar. *Genet Med*. 2017;19:1096–1104.
