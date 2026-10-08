# Computational and location evidence are not independent in missense variant classification: label-free and functional measurement, cost, and pre-registered tests of rules

*Draft 20, 8 October 2026. Supplementary material: `SUPPLEMENT_v20.md`. Pre-registrations (each committed before its analysis; not deposited in an independent registry): `prereg_fullstudy/A18_TEMPORAL_RECLASSIFICATION_PREREG.md` (commit 8a57850), `A19_V14_PREREG.md` (ceb1a11), `A20_DENOVO_RULES_PREREG.md` (2e9692a) and `A21_MAVEDB_CURATED_PM1_PREREG.md` (8cdd9c2); deviations in `A19_DEVIATIONS.md`.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** ACMG/AMP classification sums computational (PP3), hotspot/domain (PM1) and gene-level (PP2) evidence as though independent; ClinGen limits PP3 with PM1 only at Strong. We measured their dependence on truth that does not come from the codes, its cost, and which rules correct it.

**Methods.** Likelihood ratios came from de novo missense excess over mutation-rate expectation in 47,630 trios and probands, and from functional assays with expert-curated PM1 residues. Candidate rules were trained in one de novo cohort and tested in two others. Practice was audited in 87,873 laboratory and 6,177 expert-panel classifications. The assay expansion and the rule tests were pre-registered.

**Results.** Joint PP3 and PM1-like evidence was 0.44 (95% CI 0.41–0.47) of the product of their separate evidence in developmental disorders, replicated in two cohorts; within genes the dependence was confined largely to known disease genes. On functional assays in five genes, PP3 and curated PM1 were dependent (0.59, 0.44–0.80), whereas random residues were not (1.00). Stacking PP3, PM1 and PP2 over-counted 2.3 points. On independent cohorts, no rule bettered ClinGen's cap for PP3, and the cap without PP2 beside PM1 reproduced the joint evidence of all three criteria within 0.5 points, against 1.3–2.1 points over-counted with PP2. In practice the cap bound in 4.7% of laboratory PP3-with-PM1 combinations, and Moderate or Strong PP3 beside PM1 rose to 8.5% of laboratory PM1 submissions in 2025.

**Conclusion.** Computational and location evidence overlap by up to a classification step. Keeping ClinGen's cap and not adding PP2 beside PM1 corrects the overlap on independent truth.

---

## Introduction

Missense variants are the largest class of variants of uncertain significance. Under the ACMG/AMP guidelines (Richards et al. 2015), their classification draws on three kinds of evidence: a computational prediction (PP3), location in a mutational hotspot or critical functional domain (PM1), and membership of a gene with a low rate of benign missense variation (PP2). The Bayesian points framework adds these in log-likelihood-ratio units, which assumes that they carry independent information (Tavtigian et al. 2018; Tavtigian et al. 2020).

There are reasons to doubt the assumption. Predictors such as AlphaMissense, REVEL and VEST4 (Cheng et al. 2023; Ioannidis et al. 2016; Carter et al. 2013) learn which genes and protein regions tolerate variation, and that is the information PM1 and PP2 encode. Guidance has responded inconsistently. ClinGen limits the combination of PP3 and PM1 to Strong (Pejaver et al. 2022; Stenton et al. 2024); UK guidance advises against counting constraint evidence in both PM1 and PP2 (Durkie et al. 2024); and CanVIG-UK permits PP3 and PM1 together (Garrett et al. 2021). These positions rest on expert judgement rather than on a measurement of the dependence.

Measuring the dependence requires truth that the classifications themselves did not produce. ClinVar labels partly reflect PP3 and PM1, and several predictors were trained on them (Grimm et al. 2015; Harrizi et al. 2026). Truth re-scored from curators' own classifications has a subtler problem, which we show with negative-control criteria: it records which combinations of codes curators use, not only what the codes mean. Two kinds of truth avoid both problems. The excess of de novo variants over mutation-rate expectation in affected children defines a likelihood ratio without variant-level labels, for any set of sites and within any gene. Functional assays measure molecular damage directly.

We used both. We measured the dependence between PP3, PM1 and PP2 on the de novo scale and on functional assays with PM1 defined as expert panels apply it. We estimated what the overlap costs in points and in class calibration. We tested candidate rules on de novo truth in cohorts independent of those used to choose them. Finally, we audited how often, and at what strength, panels and laboratories combine the criteria.

---

## Methods

### Analysis log and pre-registration
Analyses followed a public, version-controlled analysis log of seventeen amendments, each committed before the analysis it governs (Supplementary Note 1). The primary criterion of the original plan, how much predictor evidence survives conditioning on gene and protein position, was not met under its original position control, and a fixed-window control was adopted afterwards. The PP3 × PM1 dependence test was specified before its results were seen. Four later sets of tests were pre-registered in documents committed to the repository before any of them was run: a temporal reclassification test (A18), four tests on laboratory data (A19), rule tests on de novo truth (A20) and the expanded functional-assay test (A21). None was deposited in an independent registry. All other analyses of expert-panel and laboratory practice are exploratory. Supplementary Table S0 lists every hypothesis and when it was specified relative to data access.

### Data
**De novo cohorts.** We used three cohorts (Figure 1): 31,058 developmental-disorder (DD) trios (Kaplanis et al. 2020); 12,316 autism probands, with 4,051 unaffected siblings as a negative control (Zhou et al. 2022); and a public cohort of 4,256 trios with congenital heart disease, schizophrenia, Tourette syndrome, epileptic encephalopathy or intellectual disability (Jin et al. 2017; denovo-db v1.6.1, Turner et al. 2017). Studies known to overlap the DD or autism cohorts were excluded from the public cohort. Individuals sharing two or more identical de novo variants with a DD proband were removed (64 autism probands and 8 public trios), and the public cohort was also deduplicated against the autism cohort. Unrecognised overlap cannot be excluded entirely.

**Sites and annotations.** Sites were coding positions on GENCODE v32 canonical transcripts, excluding those within 3 bp of an exon boundary. Predictor scores came from dbNSFP 4.1a and 5.4a (Liu et al. 2020) and AlphaMissense; gene and regional missense constraint from gnomAD (Karczewski et al. 2020; Chao et al. 2024); constrained coding regions from Havrilla et al. (2019); disease genes from the DDG2P monoallelic list (Thormann et al. 2019); and structure from AlphaFold DB (Varadi et al. 2022).

**Functional assays.** We used human single-substitution assays from ProteinGym v1.3 (Notin et al. 2023; 82 assays in 70 genes mapped to analysed transcripts; damaging if `DMS_score_bin` = 0) and score sets from MaveDB for genes with curated PM1 residues (81 examined).

**Curated PM1.** Curated PM1 residues are the protein positions at which any ClinGen expert panel applied PM1 to a missense variant in the Evidence Repository (export of 7 October 2026): 841 residues in 67 genes.

**Classifications.** We used 6,177 missense classifications by 46 expert panels from the Evidence Repository, and 87,873 germline ClinVar laboratory submissions (October 2026) for missense SNVs, from 472 submitters, whose comments list ACMG/AMP codes; expert-panel and practice-guideline submissions were excluded. A parser read codes and strengths from the free text and skipped negated codes. On the free-text summaries of 13,091 Evidence Repository classifications, checked against their structured codes, its precision was 0.92 and recall 0.98, and for PP3, PM1 and PP2 precision was 0.95, 0.98 and 0.98, with recall and agreement on strength of at least 0.99.

### Dependence on the de novo scale
Expected de novo counts E were derived from Roulette mutation rates (Seplyarskiy et al. 2023) and scaled to each cohort by its synonymous de novo count. The excess is X = O − E, and for a set of sites b the likelihood ratio is LR = (X_b/ΣX)/(E_b/ΣE).

PP3 was flagged at the ClinGen thresholds for AlphaMissense (primary), REVEL and VEST4 (Pejaver et al. 2022; Bergquist et al. 2025). A site was PM1-like if its regional missense observed/expected ratio was below 0.4 or its constrained-coding-region percentile was at least 95, and PM1-strict if observed/expected was below 0.2, the level Chao et al. (2024) calibrated to Moderate. A gene was PP2-like if its missense Z score was at least 3.09. This flag captures only the constraint component of PP2, not the requirement that missense variants are a common disease mechanism.

Dependence was measured as R = LR(A ∧ B)/[LR(A) · LR(B)], the joint evidence of criteria A and B relative to the product of their separate evidence: R = 1 means independence, and R < 1 means that summing them overstates the joint evidence. The points over-counted by summing are −ln R/ln 2.406. To compare PP3's evidence on one scale, we report its points alone, ln LR(PP3)/ln 2.406, and beside PM1, ln LR(PP3 | PM1)/ln 2.406, with the share kept. The de novo scale is compressed relative to the clinical scale, so we compare evidence within it and never against nominal ACMG weights. Intervals resample genes (1,000 replicates).

### Dependence on functional truth
R was computed within each assay with a position-clustered bootstrap and pooled by random effects. With the PM1-like proxy, assays were used if 5–95% of their variants lay at flagged sites. With curated PM1 in ProteinGym, genes were used if at least 90% of curated residues matched the transcript's reference amino acid, and assays if at least 10 assayed variants lay at curated residues. The MaveDB analysis (A21) was pre-registered: damage was defined within each score set as a score closer to the median of its nonsense than of its synonymous variants (at least 10 of each required), R was combined within genes and pooled across genes (gene as the unit), and random residue sets of the same size as the curated set (200 draws per score set) served as a negative control.

### The cost of stacking
DD de novo missense variants were classified with PS2 + PM2 + PM1-like + AlphaMissense PP3 with and without ClinGen's cap. Each class's disease-associated fraction is 1 − E/O, with an exact Poisson interval on O; it is a class-level estimate of association with any severe developmental disorder, compared with the class's nominal posterior (Pathogenic ≥ 0.99, Likely Pathogenic ≥ 0.90; prior 0.0441).

### Rules tested on de novo truth (A20)
For each cell of PP3 band and PM1 flag, a rule keeps a fraction κ of PP3's own evidence beside PM1 and predicts ln LR(PM1) + κ ln LR(PP3). Because κ is a fraction of evidence on the same scale, the comparison is unaffected by the compression of the de novo scale. The rules were the uncapped sum, ClinGen's cap, PP3 at most Supporting, PP3 one strength level lower, a learned constant κ, a learned κ per band, and no PP3 credit. For PP2 beside PM1 they were adding it, not adding it, or a learned fraction, and the triple was evaluated in combination. Learned rules were fitted in DD and every rule was tested in autism and the public cohort, with error in points weighted by the inverse variance of each cell and a parametric bootstrap (10,000 replicates). The pre-registered primary contrast was the rule with the lowest DD error, other than the cap, against the cap.

### Practice
**Co-application.** Co-application of PP3, PM1 and PP2 and PP3 strength beside PM1 were tabulated in panel and laboratory classifications, and among laboratory PM1 submissions by year of last evaluation. Classifications were re-scored on the points scale (P ≥ 10, LP 6–9, VUS 0–5, LB −1 to −6, B ≤ −7) under four rules: the ClinGen cap; a limit on PP3 + PM1 + PP2; PP2 not added when PM1 applies; and PP3 limited to Supporting when PM1 applies.

**Classification-derived truth.** Classifications were re-scored after removing the computational and location criteria (PP3, BP4, PM1, PP2, BP1, BP3, PM4, PM5). Truth was P/LP at 6 or more remaining points and B/LB at −1 or fewer. On this truth LR(C) = P(C applied | P/LP)/P(C applied | B/LB), and the credit C earns once PM1 applies is ln[LR(C ∧ PM1)/LR(PM1)]/ln 2.406.

**Pre-registered tests on laboratory data (A19).** Submissions from four ClinVar releases (2020–2026) were parsed identically.
- **H1 (primary):** six rules for PP3 beside PM1, compared by held-out log-loss in gene-grouped cross-validation.
- **H2 (primary):** PP2 credit beside PM1 in a mixed logistic model with submitter and gene random intercepts.
- **H3:** PP3 thresholds recalibrated per gene.
- **H4:** replication across eras of last evaluation, with PP1 and PP4 as negative-control criteria, on truth that also excluded them.

Holm correction was applied across H1 and H2.

**Temporal test (A18).** ClinVar 2019 Likely Pathogenic missense variants were compared by whether they had both PP3-like (REVEL) and PM1-like evidence or only one, with downgrade by 2026 as the outcome. The odds ratio was Mantel–Haenszel stratified by gene.

**Exploratory.** Exploratory analyses of credit (within-submitter, frequency-handled, sensitivity and replication analyses) and a held-out test of the PP2 rule are described in Supplementary Notes 5, 6 and 8.

---

## Results

### PP3, PM1 and PP2 are not independent on the de novo scale
For AlphaMissense PP3 at Supporting with a PM1-like flag, R was 0.44 [0.41, 0.47] in DD, and the dependence replicated in autism (0.56) and in the public cohort (0.67) (Table 1, Figure 2). It persisted when variants were compared within DDG2P monoallelic genes (0.69–0.75). Across all genes, however, the within-gene ratio was 0.75 in DD and 0.96 in autism and the public cohort. Much of the dependence therefore lies between genes, and within genes it is concentrated in known disease genes.

Dependence was weaker for REVEL (0.57) and VEST4 (0.67). Across 40 informative predictors in DD, every R with PM1-like had an upper bound below 1. It held in ordered and disordered regions, at buried and exposed sites, and after removing the 50 genes with the largest excess (Supplementary Note 3). PM1 and PP2 were the most dependent pair, largely by construction, because both of our flags derive from missense constraint.

**Table 1. Dependence on the de novo scale** (R; gene-clustered 95% intervals).

| Criteria | DD | Autism | Public |
|---|---|---|---|
| PP3 (AlphaMissense) × PM1-like | 0.44 [0.41, 0.47] | 0.56 [0.47, 0.66] | 0.67 [0.48, 0.96] |
| — within DDG2P genes | 0.69 [0.65, 0.73] | 0.74 [0.67, 0.82] | 0.75 [0.62, 0.89] |
| — within genes, all genes | 0.75 [0.71, 0.78] | 0.96 [0.84, 1.10] | 0.96 [0.74, 1.28] |
| PP3 (AlphaMissense) × PP2-like | 0.61 [0.57, 0.65] | 0.65 [0.55, 0.76] | 0.78 [0.57, 1.07] |
| PM1-like × PP2-like | 0.31 [0.28, 0.33] | 0.36 [0.30, 0.45] | — |

Beside PM1, AlphaMissense kept 43–49% of its own evidence and REVEL 61–78%, with similar shares at Supporting and at Strong (Table 2). The overlap arises because missense predictors largely carry location information. Across 39 informative predictors, a median 0.32 of the log-evidence of their highest-scoring 2% was incremental to gene and 50-residue window, against 0.51 for a random partition of genes into the same number of units. Inside constrained regions AlphaMissense flags 51% of mutation-expected sites, against 19% overall. A perfect predictor would also lose evidence within windows, by an amount the data cannot identify, so this share measures evidence incremental to location rather than predictor skill (Supplementary Note 2).

**Table 2. PP3 evidence alone and beside PM1 on the de novo scale** (points on the 2.406 scale; R with gene-clustered 95% interval).

| Predictor, cohort | PP3 at | Alone | Beside PM1-like | Share kept | R |
|---|---|---|---|---|---|
| AlphaMissense, DD | ≥ Supporting | 1.63 | 0.70 | 43% | 0.44 [0.41, 0.47] |
| AlphaMissense, DD | Strong | 2.36 | 1.10 | 47% | 0.33 [0.31, 0.35] |
| REVEL, DD | ≥ Supporting | 1.65 | 1.01 | 61% | 0.57 [0.53, 0.61] |
| REVEL, DD | Strong | 2.30 | 1.58 | 68% | 0.53 [0.48, 0.58] |
| AlphaMissense, autism | ≥ Supporting | 1.26 | 0.59 | 47% | 0.56 [0.47, 0.66] |
| AlphaMissense, autism | Strong | 1.83 | 0.90 | 49% | 0.44 [0.36, 0.54] |
| REVEL, autism | ≥ Supporting | 1.01 | 0.79 | 78% | 0.82 [0.69, 0.99] |
| REVEL, autism | Strong | 1.59 | 1.25 | 78% | 0.74 [0.57, 0.97] |

### The dependence on functional truth
Functional assays measure damage without reference to the codes. With PM1 defined as expert panels apply it, the dependence appeared there too (Table 3).
- **ProteinGym.** In 7 assays in 3 genes (TP53, GCK, HRAS), R was 0.49 [0.40, 0.60] and below 1 in every assay. At curated PM1 residues AlphaMissense PP3 kept a median 40% of its log-evidence for damage, close to the 43–49% it keeps beside PM1-like sites on the de novo scale. This happens because the predictor already scores most variants at those residues as damaging: 77% were PP3-positive there, against 34% elsewhere.
- **MaveDB (pre-registered).** Across 5 genes and 12 score sets, R was 0.59 [0.44, 0.80] at Supporting and 0.48 [0.30, 0.76] at Strong. R was below 1 in four genes (GCK 0.49, LDLR 0.52, PTEN 0.57, TP53 0.61), and the fifth, F9, was degenerate. Leaving out any gene gave 0.53–0.71.
- **Negative control.** Random residue sets of the same size gave R of 1.00 (0.98–1.08 by gene), so the dependence is specific to the residues curators flag.
- **Sensitivity.** Two problems were found after the pre-registered run and are reported as sensitivity analyses. One TP53 score set with a dominant-negative readout had inverted truth, and three F9 score sets had a joint cell with no non-damaging variants. Excluding them gave 0.52 [0.47, 0.57] (Supplementary Note 13).

The constraint proxy behaved differently. It barely predicted damage within assays (median LR 1.10 outside DDG2P genes), and its dependence was correspondingly weak (0.94). It flags none of the assayed TP53 residues, 3% in BRCA1 and nearly all in PTEN, so none of these genes could be tested with it.

**Table 3. Dependence on truth independent of the de novo data.**

| Truth | Criteria | Scope | R [95% CI] |
|---|---|---|---|
| Functional assays (ProteinGym) | AlphaMissense PP3 × curated PM1 | 7 assays, 3 genes | 0.49 [0.40, 0.60] |
| Functional assays (MaveDB; pre-registered) | AlphaMissense PP3 × curated PM1 | 12 score sets, 5 genes | 0.59 [0.44, 0.80] |
| — negative control | AlphaMissense PP3 × random residues, same number | same | 1.00 (0.98–1.08 by gene) |
| Functional assays (ProteinGym) | AlphaMissense PP3 × PM1-like proxy | 13 assays, 10 genes outside DDG2P | 0.94 [0.88, 0.99] |
| ClinVar P/LP vs B/LB* | AlphaMissense PP3 × PM1-like | DDG2P genes | 0.32 [0.28, 0.37] |
| Laboratory truth† | Laboratory PP3 × PM1 | 39,015 submissions, 5,976 genes | 0.18 [0.14, 0.25] |
| — negative controls† | PP1 × PM1; PP4 × PM1 | 38,007 submissions | 0.04 [0.02, 0.22]; 0.03 [0.02, 0.12] |
| Expert-panel truth† | Panel PP3 × panel PM1 | 2,942 classifications, 130 genes | 0.21 [0.12, 0.31] |

\*Partly circular, because ClinVar classifications use PP3 and PM1 themselves, which biases R towards 1. †Classification-derived truth. The negative controls show that low R on such truth is not specific to shared information, so these rows describe practice rather than measure dependence.

On classification-derived truth, R was also below 1 for PP3 with PM1, in every independent half of the laboratory data and in every era of last evaluation from before 2020 to 2026 (0.10–0.27; pre-registered H4). The pre-registered negative controls, however, behaved the same way. PP1 and PP4, which carry no information about location, had R with PM1 of 0.04 and 0.03, lower than PP3's. Low R on truth re-scored from curators' classifications therefore arises largely from which combinations benign classifications contain. We treat it as a description of practice, and rest the measurement of dependence on the de novo and functional truth.

### The cost of stacking
In points, Supporting AlphaMissense PP3 beside a PM1-like flag over-counts 0.97 in DD. The triple PP3 × PM1-like × PP2-like (R3 0.13 [0.11, 0.14]) over-counts 2.32 points [2.24, 2.51] in DD and 1.95 [1.63, 2.32] in autism, about the width of a classification step. Classified with PS2 + PM2 + PM1-like + AlphaMissense PP3 and no cap, 2,438 DD de novo variants were Pathogenic, a class the de novo excess implies is 84% [84, 85] disease-associated, against a nominal 99% (Figure 3). Under ClinGen's cap the same variants are Likely Pathogenic. The cap does not make that class calibrated either (48% against a nominal 90%). Part of that shortfall reflects the estimand, which is compressed relative to the clinical scale, so the class comparison shows direction and the points over-counted show size. In ClinVar, PP3 with Moderate PM1 alone would make 7,101 PP3-positive VUS in DDG2P monoallelic genes Likely Pathogenic without the cap (Supplementary Note 7).

### Which rules correct the overlap: tests on independent de novo truth (A20)
Rules were trained in DD and tested in autism and the public cohort (Supplementary Note 12).
- **PP3 beside PM1.** The uncapped sum over-counted (error relative to the cap +0.80 points², 95% interval 0.51–1.09), and no PP3 credit under-counted (+0.51). The pre-registered primary contrast did not beat the cap. The rule with the lowest DD error, a learned κ per band, overfitted and failed in the test cohorts (+0.42 [0.05, 0.90]), because AlphaMissense PP3 in its Supporting band carries almost no evidence alone on this scale (LR 1.05). Among secondary contrasts, counting PP3 one strength level lower beside PM1 was better calibrated than the cap for AlphaMissense (−0.06 [−0.11, −0.02]) but not for REVEL (+0.02 [−0.11, +0.15]). These intervals treat cells as independent.
- **PP2 beside PM1.** In autism, adding PP2 over-counted the joint evidence by 1.16 points (SE 0.12), not adding it under-counted by 0.62, and a fraction learned in DD (0.33) was within 0.02 points.
- **The triple (descriptive).** The combination laboratories use, the cap with PP2 added, over-counted in all four autism cells by 1.31–2.09 points. The cap without PP2 beside PM1 was within 0.5 points in all four (−0.31 to +0.47).

On truth that does not derive from the codes, then, the cap is not bettered for PP3, and not adding PP2 beside PM1 is the change that corrects the remaining overlap.

### How panels and laboratories combine the criteria
Panels applied PP3 beside PM1 in 787 classifications and laboratories in 4,471 (Table 4). In both, PP3 was almost always at Supporting or Moderate strength (99% and 98%), and in laboratories PM1 was Moderate in 90%, so the cap, which binds only when the two sum to more than 4 points, bound in 211 of 4,471 laboratory combinations (4.7%). Panels applied PP2 beside PM1 in 174 classifications, 145 of them from three panels whose specifications may permit it (RASopathy, FBN1 and Monogenic Diabetes); laboratories did so in 1,481 submissions, 1,035 of which also applied PP3.

**Table 4. Co-application in practice.**

| | Expert panels | Laboratories |
|---|---|---|
| Missense classifications | 6,177 | 87,873 |
| PP3 beside PM1 | 787 | 4,471 |
| — PP3 Supporting / Moderate / Strong | 605 / 178 / 4 | 4,207 / 184 / 80 |
| Share with PP3 at Strong | 0.5% | 1.8% |
| PP2 beside PM1 | 174 | 1,481 |
| All three | 138 | 1,035 |

Stronger combinations are becoming more common. Among laboratory PM1 submissions, PP3 at Moderate or Strong was applied in 0–0.9% each year from 2018 to 2023, then in 4.3% in 2024 and 8.5% in 2025, after calibrated PP3 strengths were introduced (Pejaver et al. 2022). PP2 beside PM1 rose from 18–20% in 2022–2023 to 24% in 2025 and 31% in 2026 (both trends p < 10⁻⁸; Supplementary Note 8).

Rules that remove these overlaps would change few classifications overall but many of those they apply to (Table 5). Not adding PP2 beside PM1 would change 352 of 87,873 laboratory classifications (0.4%). That is 24% of the 1,481 that combine the two, and 223 of their 624 Likely Pathogenic calls (36%) would become VUS. The changes span 160 genes and 66 laboratories, led by NOTCH3, FBN1, COL4A3, CASR, POLG, CACNA1A, KCNH2 and SCN1A. ClinGen's cap changes 40. Of the 38 submissions applying both PP3 and PM1 whose class the cap would change, laboratories reported the uncapped class in 32 (Supplementary Note 10).

**Table 5. Classifications that change class under alternative rules.**

| Rule | Panels: changed | P→LP | LP→VUS | Other | Laboratories: changed | P→LP | LP→VUS | Other |
|---|---|---|---|---|---|---|---|---|
| ClinGen cap (PP3 + PM1 ≤ 4) | 24 | 12 | 12 | 0 | 40 | 30 | 10 | 0 |
| PP3 + PM1 + PP2 ≤ 4 | 35 | 17 | 18 | 0 | 63 | 39 | 22 | 2 |
| PP2 not added when PM1 applies | 39 | 28 | 11 | 0 | 352 | 124 | 223 | 5 |
| PP3 limited to Supporting when PM1 applies | 43 | 22 | 21 | 0 | 90 | 39 | 51 | 0 |

"Other" counts any other change of class.

### Pre-registered tests on laboratory data (A19, A18)
In the pre-registered rule comparison (H1; 2,812 PM1-applied submissions with truth, 701 genes), giving PP3 no credit beside PM1 was worse calibrated than the cap (Δ log-loss +0.038 [0.024, 0.050]), so PP3 is informative there. No rule was better calibrated than the cap: a learned credit per strength gave −0.004 [−0.010, +0.011] and PP3 limited to Supporting +0.000 [−0.001, +0.001]. The verdict was unchanged under the ACMG convention that an absent code carries no evidence (deviation 1). Recalibrating AlphaMissense thresholds per gene (H3; 341 genes) reduced but did not remove the dependence (R 0.41 [0.30, 0.61], against 0.31 with published thresholds), so it is not an artefact of global thresholds.

In the pre-registered mixed model (H2), PP2 earned 0.35 points [95% credible interval −0.12, 0.81] beside PM1, below its nominal point. The posterior probability of at least 1 point was 0.003, which passes Holm correction. Alone, PP2 was calibrated (1.12 [0.90, 1.31]). The negative controls qualify this result as they qualify R: on the same truth, PP1 and PP4 also lost most of their credit beside PM1 (5.40 to 1.34 and 5.00 to 0.76 points). On classification-derived truth, therefore, a loss of credit beside PM1 is not by itself evidence of redundancy; the de novo test above is the one that measures it.

The pre-registered temporal test (A18) was null. Among 5,547 ClinVar 2019 Likely Pathogenic missense variants, those with both REVEL PP3 and PM1-like evidence were downgraded by 2026 less often, not more, than those with one (2.5% against 6.2%). The gene-stratified odds ratio was 0.53 [0.22, 1.07], and the minimum detectable odds ratio was about 3. A third of 2019 Likely Pathogenic variants had dropped out of the 2026 table, about half of them to conflicting status, so the test is weak (Supplementary Note 14).

### Credit on classification-derived truth (descriptive)
Table 6 gives the points PP2 and PP3 earned once PM1 applied on classification-derived truth. Given the negative controls, these describe how criteria co-occur in classifications rather than what they are worth. Supporting PP3 beside PM1 earned 0.81 points [0.24, 1.68] in panels and 0.84 [0.52, 1.26] in laboratories. This was not robust: removing 18 heavily templated submitters (1,681 submissions) lowered laboratory PP3 credit beside PM1 to 0.28 [0.02, 0.66] and PP2 credit to −0.80 [−1.24, −0.13] (Supplementary Note 10). Moderate PP3 beside PM1 earned 0.02 points [−1.69, 0.89] against 2 in panels, on 6 B/LB variants, with a minimum detectable difference from Supporting of about 4 points. Exploratory within-submitter, frequency-handled and held-out analyses of PP2 agreed in direction but were not robust (Supplementary Note 8).

**Table 6. Points earned once PM1 applies on classification-derived truth, against points applied** (descriptive; gene-clustered 95% intervals; n and B/LB in the joint cell).

| Criterion | Expert panels | Laboratories | Laboratories, mixed model (H2) | Applied |
|---|---|---|---|---|
| PP2 | not estimable (1 B/LB) | −0.22 [−0.67, 0.40] (565; 57) | 0.35 [−0.12, 0.81] | 1.00 |
| PP3, any strength | 0.60 [0.14, 1.05] (304; 17) | 0.85 [0.54, 1.22] (1,768; 74) | — | 1.24 / 1.09 |
| PP3 Supporting | 0.81 [0.24, 1.68] (235; 11) | 0.84 [0.52, 1.26] (1,674; 71) | — | 1.00 |
| PP3 Moderate | 0.02 [−1.69, 0.89] (67; 6) | not estimable (3 B/LB) | — | 2.00 |

---

## Discussion

Computational evidence is not independent of the location evidence that ACMG/AMP credits separately. We measured this on two kinds of truth that the codes did not produce. On the de novo scale the dependence held in three cohorts and within known disease genes. On functional assays, with PM1 defined as curators apply it, it held in five of the six genes tested (the sixth, F9, was degenerate), and random residues showed none. It arises because predictors have learned where pathogenic variation falls, which is what PM1 and PP2 describe. Beside PM1, PP3 keeps about 40–50% of its evidence on both scales. The dependence is concentrated between genes and in known disease genes; within genes elsewhere it is small.

Classification-derived truth could not measure this. R was low on it, but it was low for negative-control criteria with no location content too. That finding applies beyond this study: estimates of dependence or credit on truth re-scored from curators' own classifications, including our earlier ones, need negative-control criteria, because without them low R can reflect which combinations of codes benign classifications contain rather than shared information.

The overlap is large enough to matter. Stacked PP3, PM1 and PP2 over-count about two points, the width of a classification step, and an uncapped Pathogenic class built on them was 84% rather than 99% disease-associated. ClinGen's cap handles most of what matters for PP3. On independent de novo truth no rule bettered it for PP3, and on laboratory data no rule was better calibrated. The cap leaves PP2 uncorrected. On independent truth, not adding PP2 beside PM1, with the cap, reproduced the joint evidence of all three criteria within half a point, whereas adding PP2 over-counted by one to two points. This is what UK guidance and the ClinGen cardiomyopathy panel already advise (Durkie et al. 2024; Kelly et al. 2018). Our data supply the measurement those positions lacked.

In practice the cap rarely acts: it bound in 4.7% of laboratory PP3-with-PM1 combinations, because PP3 is almost always Supporting or Moderate and PM1 Moderate. Stronger PP3 beside PM1 has spread since calibrated PP3 strengths were introduced. Most such combinations (155 of 264 in laboratories) are Moderate PP3 beside Moderate PM1, which sums exactly to the cap's limit, so the cap allows them in full, although on de novo and functional truth PP3 keeps only about half of its evidence beside PM1. Whether Moderate PP3 should count one level lower beside PM1 is unresolved. That rule was better calibrated than the cap for AlphaMissense on de novo truth but not for REVEL, and clinical truth is too thin to test it. It would apply to the 264 laboratory classifications with Moderate or Strong PP3 beside PM1, 6% of PP3-with-PM1 co-applications, so it should be tested prospectively before it is adopted.

Later reclassification gave no support for over-crediting in practice. The pre-registered temporal test was null, and its point estimate ran the other way: stacked Likely Pathogenic calls from 2019 were downgraded less often than single-evidence ones. Variants with both kinds of evidence may truly be more often pathogenic, and the test had little power and heavy drop-out. Even so, it argues against reading the de novo over-count as a large burden of wrong calls already in ClinVar.

**Implications for curation.**
1. **Keep ClinGen's cap.** No pre-registered rule bettered it for PP3, on laboratory truth or on independent de novo truth.
2. **Do not add PP2 when PM1 applies,** unless a specification documents that the two rest on independent evidence. This rule corrects the overlap on independent truth. In laboratory practice it would change 24% of the classifications that combine the two, mostly from Likely Pathogenic to VUS.
3. **Test, rather than adopt, counting Moderate PP3 one level lower beside PM1.**
4. **Report conditional evidence in calibrations.** Calibrations of new predictors and location criteria should report evidence conditional on the other criteria, not only marginal evidence. The de novo method needs no variant labels and can be applied to any predictor or region definition with the released code.
5. **Include negative-control criteria** in analyses on classification-derived truth.

**Version 4.** The forthcoming ACMG/AMP/CAP/ClinGen standard reorganises evidence hierarchically, in part to stop related codes being counted twice, and reportedly groups predictive and functional evidence apart from observational evidence [verify against the version 4 specification]. Our measurements support placing PP3, PM1 and PP2 in a shared, capped group. If version 4 retires PP2 or folds it into location evidence, the PP2 findings support that design rather than patch the 2015 rules.

**Relation to previous work.** That PP3, PM1 and PP2 can double-count the same attribute is not new. It motivates ClinGen's cap, UK guidance, Sherloc's rule dependencies (Nykamp et al. 2017) and at least one expert-panel specification (Supplementary Table S17), all of which rest on expert judgement. PP3 thresholds (Pejaver et al. 2022; Bergquist et al. 2025) and regional constraint for PM1 (Chao et al. 2024) are calibrated one criterion at a time, and gene- and domain-aware calibration addresses heterogeneity between genes rather than dependence between criteria (Tejura et al. 2024; Badonyi & Marsh 2025; Chen et al. 2026). Laboratories are known to apply the criteria inconsistently (Amendola et al. 2016; Harrison et al. 2017). Our contribution is the measurement: a label-free and functional estimate of the dependence, its cost, and rules tested on independent truth. We found no previous study doing any of these.

**Limitations.**
- **Pre-registration.** The analysis log and the pre-registrations were committed to the repository, not an independent registry. The original primary criterion was not met, and the emphasis on dependence followed that failure, although the dependence test itself was specified before its results.
- **The de novo scale.** It reflects high-penetrance variants in severe early-onset disorders and is compressed relative to the clinical scale, so we compare evidence within it.
- **Proxy flags.** Our PM1 and PP2 flags on the de novo scale are constraint-based proxies, and panel definitions of PM1 vary.
- **Rule tests (A20).** They are cell-level and use intervals that treat cells as independent. They rest on two test cohorts (one for the triple), their primary contrast failed through overfitting, and the triple comparison was descriptive.
- **Functional assays.** Curated PM1 could be tested in six genes, three in ProteinGym and five in MaveDB, two of them in both. Curated residues are where pathogenic variants had already been classified. Damage thresholds from nonsense and synonymous variants fail for dominant-negative readouts.
- **Classification-derived analyses.** They are descriptive. Benign truth there is mostly frequency-based, estimates vary between submitters and are sensitive to templated submissions, and the parser was validated against structured codes rather than double-annotated free text (a 400-comment annotation sheet is provided).
- **Coverage.** Transcripts match MANE Select in 74.4% of genes, and the cohorts are predominantly of European ancestry.

**Conclusion.** PP3, PM1 and PP2 are not independent. On de novo and functional truth, PP3 keeps about half of its evidence beside PM1, and the three criteria stacked over-count by up to a classification step. ClinGen's cap should be kept, PP2 should not be added beside PM1, and stronger PP3 beside PM1, which is spreading below the cap, should be tested for a further limit.

---

## Figure legends
**Figure 1. Study design.** (a) Label-free dependence: de novo excess over mutation-rate expectation in three cohorts gives likelihood ratios for PP3, PM1-like and PP2-like sites, from which dependence R, PP3's evidence beside PM1 and points over-counted are estimated. (b) Functional truth: R on ProteinGym and MaveDB assays with PM1 as curated by expert panels, against random residue sets and the constraint proxy. (c) Rules and practice: rules trained in one de novo cohort and tested in two others (A20); co-application of the criteria in expert-panel and laboratory classifications, with trends and re-scoring; pre-registered tests on laboratory data with negative-control criteria (A19) and on later reclassification (A18) (`fig1_design.png`).

**Figure 2. Dependence on the de novo scale.** R for PP3 with PM1-like and PP2-like flags by predictor, flag and cohort, with gene-clustered 95% intervals; R = 1 is independence (`fig2_dependence.png`).

**Figure 3. The cost of stacking.** Disease-associated fraction (de novo excess, DD) of the classes produced by PS2 + PM2 + PM1-like + AlphaMissense PP3 without and with ClinGen's cap, with exact 95% intervals, against each class's nominal posterior band (`figures/fig_cost_calibration.png`).

## Data and code availability
All data are public: published de novo tables, ClinVar, the ClinGen Evidence Repository, ProteinGym, MaveDB, dbNSFP (academic licence), AlphaMissense, AlphaFold DB, gnomAD, DDG2P and Roulette. Code, the analysis log, pre-registrations and all outputs are in the study repository (https://github.com/tanukusaitejesh-prog/ndd-pp3-calibration; to be made public with a Zenodo DOI at submission). `REPRODUCE.md` maps every number, table and figure to its script. `reproduce_v14.sh` re-runs the pre-registered laboratory analyses from public ClinVar files, and `make_denovo_rules.py` re-runs the de novo rule tests from the saved de novo estimates.

## Competing interests
[To be completed.]

## Funding
[To be completed.]

## Author contributions
[To be completed.]

---

## References
- ACMG/AMP/CAP/ClinGen. Standards for sequence variant classification, version 4 (draft and pilot materials). [Cite the publication or pilot report once available.]
- Amendola LM, Jarvik GP, Leo MC, et al. Performance of ACMG-AMP variant-interpretation guidelines among nine laboratories in the Clinical Sequencing Exploratory Research Consortium. *Am J Hum Genet* 2016;98:1067–1076.
- Badonyi M, Marsh JA. acmgscaler: an R package and Colab for standardized gene-level variant effect score calibration within the ACMG/AMP framework. 2025.
- Bergquist T, Stenton SL, Nadeau EAW, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 2025;27:101402.
- Carter H, Douville C, Stenson PD, Cooper DN, Karchin R. Identifying Mendelian disease genes with the variant effect scoring tool. *BMC Genomics* 2013;14(Suppl 3):S3.
- Chao KR, Wang L, Panchal R, et al. The landscape of regional missense mutational intolerance quantified from 125,748 exomes. bioRxiv 2024. doi:10.1101/2024.04.11.588920.
- Chen Y, Fayer S, Jain S, Benazouz M, Sverchkov Y, Stone J, Sharma H, Bergquist T, Stewart R, Mooney SD, Craven M, Radivojac P, Starita LM, Fowler DM, Pejaver V. Gene- and domain-aware calibration increases the clinical utility of variant effect predictors. bioRxiv 2026. doi:10.1101/2026.02.17.706269.
- Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 2023;381:eadg7492.
- Durkie M, Cassidy E-J, Berry I, et al. ACGS best practice guidelines for variant classification in rare disease 2024. Association for Clinical Genomic Science, 2024.
- Garrett A, Durkie M, Callaway A, et al. Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. *J Med Genet* 2021;58:297–304.
- Grimm DG, Azencott CA, Aicheler F, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat* 2015;36:513–523.
- Harrizi S, Nait Irahal I, Mostafa K, Arnoult D. Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. bioRxiv 2026. doi:10.64898/2026.08.21.746180.
- Harrison SM, Dolinsky JS, Knight Johnson AE, et al. Clinical laboratories collaborate to resolve differences in variant interpretations submitted to ClinVar. *Genet Med* 2017;19:1096–1104.
- Havrilla JM, Pedersen BS, Layer RM, Quinlan AR. A map of constrained coding regions in the human genome. *Nat Genet* 2019;51:88–95.
- Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *Am J Hum Genet* 2016;99:877–885.
- Jin SC, Homsy J, Zaidi S, et al. Contribution of rare inherited and de novo variants in 2,871 congenital heart disease probands. *Nat Genet* 2017;49:1593–1601.
- Kaplanis J, Samocha KE, Wiel L, et al. Evidence for 28 genetic disorders discovered by combining healthcare and research data. *Nature* 2020;586:757–762.
- Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 2020;581:434–443.
- Kelly MA, Caleshu C, Morales A, et al. Adaptation and validation of the ACMG/AMP variant classification framework for MYH7-associated inherited cardiomyopathies: recommendations by ClinGen's Inherited Cardiomyopathy Expert Panel. *Genet Med* 2018;20:351–359.
- Liu X, Li C, Mou C, Dong Y, Tu Y. dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs. *Genome Med* 2020;12:103.
- Notin P, Kollasch AW, Ritter D, et al. ProteinGym: large-scale benchmarks for protein fitness prediction and design. *Adv Neural Inf Process Syst* 2023.
- Nykamp K, Anderson M, Powers M, et al. Sherloc: a comprehensive refinement of the ACMG-AMP variant classification criteria. *Genet Med* 2017;19:1105–1117.
- Pejaver V, Byrne AB, Feng BJ, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet* 2022;109:2163–2177.
- Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med* 2015;17:405–424.
- Seplyarskiy V, Koch EM, Lee DJ, et al. A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription. *Nat Genet* 2023;55:2235–2242.
- Stenton SL, Pejaver V, Bergquist T, et al. Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med* 2024;26:101213.
- Tavtigian SV, Greenblatt MS, Harrison SM, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med* 2018;20:1054–1060.
- Tavtigian SV, Harrison SM, Boucher KM, Biesecker LG. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Hum Mutat* 2020;41:1734–1737.
- Tejura M, Fayer S, McEwen AE, et al. Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes. *Am J Hum Genet* 2024;111:2031–2043.
- Thormann A, Halachev M, McLaren W, et al. Flexible and scalable diagnostic filtering of genomic variants using G2P with Ensembl VEP. *Nat Commun* 2019;10:2373.
- Turner TN, Yi Q, Krumm N, et al. denovo-db: a compendium of human de novo variants. *Nucleic Acids Res* 2017;45:D804–D811.
- Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Res* 2022;50:D439–D444.
- Zhou X, Feliciano P, Shu C, et al. Integrating de novo and inherited variants in 42,607 autism cases identifies mutations in new moderate-risk genes. *Nat Genet* 2022;54:1305–1319.
