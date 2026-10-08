# Double counting of computational and location evidence in missense classification: its cost, and a pre-registered test of the rules that contain it

*Draft 14, 8 October 2026. Supplementary material: `SUPPLEMENT_v14.md`. Pre-registration of the new tests:
`prereg_fullstudy/A19_V14_PREREG.md` (committed before analysis, commit ceb1a11); deviations:
`prereg_fullstudy/A19_DEVIATIONS.md`.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** ACMG/AMP classification adds computational (PP3), hotspot/domain (PM1) and gene-level (PP2) evidence as
though independent. We measured what this costs, and tested the rules that limit it.

**Methods.** Dependence was estimated without variant-level labels from the excess of de novo missense variants over
mutation-rate expectation (47,630 trios and probands). Pre-registered tests used 87,873 laboratory ClinVar
classifications re-scored on evidence independent of these criteria, across four ClinVar releases (2020–2026).

**Results.** Joint PP3 and PM1-like evidence was 0.44 (95% CI 0.41–0.47) of the product of their separate evidence in
developmental disorders; adding PP2 over-counted 2.3 points (2.2–2.5). Uncapped stacking placed 2,438 de novo variants
at Pathogenic in a class 84% disease-associated, against a nominal 99%. Out of sample, no rule calibrated PP3 beside PM1
better than ClinGen's cap. PP2 earned 0.35 points (−0.12 to 0.81) beside PM1 against 1 nominal, although alone it was
calibrated (1.12). On classification-derived truth, negative-control criteria showed the same apparent loss of credit,
so laboratory estimates describe practice rather than measure dependence.

**Conclusion.** Computational and location evidence overlap enough to move variants a class when stacked. ClinGen's cap
should stay and should extend to PP2 beside PM1; the forthcoming version 4 framework's grouping of related evidence is
supported.

---

## Introduction

Missense variants are the largest class of variants of uncertain significance. Under the ACMG/AMP guidelines (Richards
et al. 2015) their classification draws on a computational prediction (PP3), location in a mutational hotspot or
critical domain (PM1), and membership of a gene with little benign missense variation (PP2). The Bayesian points
framework adds these in log-likelihood-ratio units, which assumes they carry independent information (Tavtigian et al.
2018; 2020).

Predictors such as AlphaMissense, REVEL and VEST4 (Cheng et al. 2023; Ioannidis et al. 2016; Carter et al. 2013) learn
which genes and protein regions tolerate variation — the information PM1 and PP2 encode. Guidance has responded
inconsistently: ClinGen limits PP3 + PM1 to Strong (Pejaver et al. 2022; Stenton et al. 2024); UK guidance advises
against counting constraint in both PM1 and PP2 (Durkie et al. 2024); CanVIG-UK permits PP3 and PM1 together (Garrett
et al. 2021). The forthcoming ACMG/AMP/CAP/ClinGen version 4 standard reorganises evidence hierarchically, in part to
stop related codes being counted twice [verify against the v4 specification]. None of these positions rests on a
measurement of the dependence or of its consequence.

Measuring dependence needs evidence the classifications did not produce. De novo variants provide it: their excess over
mutation-rate expectation in affected children defines a likelihood ratio without variant-level labels. We used it to
measure the dependence and its cost, then ran pre-registered tests of the rules that would contain it on laboratory
practice, with negative controls, out-of-sample comparison and replication across ClinVar releases.

---

## Methods

### Pre-registration and analysis log
- **De novo analyses** followed a version-controlled log of seventeen amendments, each committed before the analysis it
  governs (Supplementary Note 1). It was not deposited in an independent registry; the pre-specified primary criterion
  on predictor evidence was not met under its original position control.
- **Draft-14 tests** (H1–H4) were pre-registered in a document committed to the public repository before any of them
  was run (A19; commit ceb1a11). Supplementary Table S0 lists every hypothesis with the date it was specified relative
  to data access. One deviation was made after running H1 and is reported beside the pre-registered result.
- **Everything else** in the practice analyses is labelled exploratory.

### De novo dependence
Data: developmental disorders (DD; 31,058 trios; Kaplanis et al. 2020), autism (12,316 probands, 4,051 unaffected
siblings; Zhou et al. 2022) and a public cohort (4,256 trios; Jin et al. 2017; Turner et al. 2017), deduplicated
against each other. Sites: coding positions on GENCODE v32 canonical transcripts. Expected counts E from Roulette rates
(Seplyarskiy et al. 2023) scaled by synonymous de novos; excess X = O − E; for sites b, LR = (X_b/ΣX)/(E_b/ΣE).
Flags: PP3 at the ClinGen thresholds for AlphaMissense (primary), REVEL and VEST4 (Pejaver et al. 2022; Bergquist et
al. 2025); PM1-like = regional missense o/e < 0.4 or CCR percentile ≥ 95 (Chao et al. 2024; Havrilla et al. 2019);
PP2-like = gene missense Z ≥ 3.09 (constraint component only). Dependence: R = LR(A ∧ B)/[LR(A)·LR(B)]; points
over-counted = −ln R / ln 2.406 (ClinGen scale; Pejaver et al. 2022). Intervals resample genes (1,000).

### The cost of stacking
DD de novo missense variants were classified with PS2 + PM2 + PM1-like + AlphaMissense PP3 under uncapped, capped and
conditional rules. Each class's disease-associated fraction is 1 − E/O (de novo excess), with an exact Poisson interval
on O. It is a class-level estimate of association with any severe developmental disorder, compared with the class's
nominal posterior (Pathogenic ≥ 0.99; Likely Pathogenic ≥ 0.90; prior 0.0441).

### Laboratory classifications
87,873 germline ClinVar submissions (October 2026; 472 submitters) for missense SNVs whose comments list ACMG/AMP codes;
expert-panel submissions excluded; parser validated on 13,091 Evidence Repository summaries (precision 0.92, recall
0.98). Truth = points from all codes except PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5: P/LP ≥ 6, B/LB ≤ −1. Credit of a
criterion C beside PM1 = ln[odds(truth | C, PM1)/odds(truth | PM1)]/ln 2.406 (ACMG convention: an absent code carries
no evidence). Earlier releases (January 2020, 2022, 2024) were processed identically; each submission counted once.

### Pre-registered tests (A19)
- **H1, rules for PP3 beside PM1 (primary).** Six rules (uncapped; ClinGen cap; PP3 ≤ Supporting; discount;
  learned credit per strength; no credit) predicted truth among PM1-applied submissions; gene-grouped 5-fold
  cross-validation × 20; held-out log-loss; submitter-clustered bootstrap.
- **H2, PP2 beside PM1 (primary).** Bayesian mixed logistic model with submitter and gene random intercepts; PP2 credit
  < 1 point if the 95% credible interval excludes 1. Mechanism (H2b): PM1 rationale from the comment text.
- **H3, calibrate rather than cap (secondary).** AlphaMissense PP3 thresholds recalibrated per gene on ClinVar,
  leaving out the evaluation variants.
- **H4, replication (secondary).** Four evaluation eras; negative-control criteria PP1 and PP4 (with truth also
  excluding them).
Holm correction across the two primary tests.

---

## Results

### 1. PP3, PM1 and PP2 are not independent, and stacking them over-counts by up to three points
For AlphaMissense PP3 at Supporting with a PM1-like flag, R was 0.44 [0.41, 0.47] in DD, replicating in autism (0.56)
and the public cohort (0.67; Table 1, Figure 2). Every one of 40 informative predictors had R × PM1-like below 1 in DD.
Within DDG2P genes the dependence persisted (0.69–0.75); across all genes much of it lay between genes (within-gene
0.75 in DD, 0.96 in autism). PP3 kept 43–49% of its own evidence beside PM1 (AlphaMissense) and 61–78% (REVEL).

In points, a Supporting PP3 beside a PM1-like flag over-counts 0.97 in DD; with PP2 added, the triple over-counts
2.32 points [2.24, 2.51] in DD and 1.95 [1.63, 2.32] in autism, and 2.6–3.0 with strict PM1 (Table 2). This is the
dominant real-world combination: 1,035 of the 1,481 laboratory submissions combining PM1 and PP2 also apply PP3.

**Table 1. Dependence on the de novo scale** (R; gene-clustered 95% intervals).

| Criteria | DD | Autism | Public |
|---|---|---|---|
| PP3 (AlphaMissense) × PM1-like | 0.44 [0.41, 0.47] | 0.56 [0.47, 0.66] | 0.67 [0.48, 0.96] |
| — within DDG2P genes | 0.69 [0.65, 0.73] | 0.74 [0.67, 0.82] | 0.75 [0.62, 0.89] |
| — within genes, all genes | 0.75 [0.71, 0.78] | 0.96 [0.84, 1.10] | 0.96 [0.74, 1.28] |
| PP3 (AlphaMissense) × PP2-like | 0.61 [0.57, 0.65] | 0.65 [0.55, 0.76] | 0.78 [0.57, 1.07] |
| PM1-like × PP2-like | 0.31 [0.28, 0.33] | 0.36 [0.30, 0.45] | — |
| PP3 × PM1-like × PP2-like (R3) | 0.13 [0.11, 0.14] | 0.18 [0.13, 0.24] | — |

### 2. The cost: stacked evidence moves variants a class
Uncapped PS2 + PM2 + PM1-like + AlphaMissense PP3 placed 2,438 DD de novo variants at Pathogenic. The de novo excess
implies the class is 84% [84, 85] disease-associated, against a nominal ≥ 99% (Table 2, Figure 3). Under ClinGen's cap
the same variants are Likely Pathogenic. The cap does not make the Likely Pathogenic class calibrated (48% against a
nominal ≥ 90%); part of that shortfall reflects the estimand (association with any severe developmental disorder,
compressed relative to the clinical scale), so we use the class-level comparison to show direction and the
points over-counted to show size.

In ClinVar, 37.3% of AlphaMissense PP3-positive VUS in DDG2P monoallelic genes lie at PM1-like sites (13.2% across all
genes). Without the cap, PP3 + PM1 at Moderate alone would make 7,101 of them Likely Pathogenic; under every limiting rule
tested (cap, Supporting limit, learned credit), none.

**Table 2. Cost of stacking.**

| Quantity | Value |
|---|---|
| Points over-counted, PP3 Supporting × PM1-like (DD / autism) | 0.97 / 0.23 |
| Points over-counted, PP3 × PM1-like × PP2-like (DD / autism) | 2.32 [2.24, 2.51] / 1.95 [1.63, 2.32] |
| Uncapped Pathogenic class, DD (n; disease-associated vs nominal) | 2,438; 84% vs ≥ 99% |
| ClinVar VUS reaching LP from PP3 + PM1 alone, DDG2P genes: uncapped / any limiting rule | 7,101 / 0 |

### 3. Pre-registered test of rules: no rule beats ClinGen's cap (H1)
Among 2,812 PM1-applied laboratory submissions with truth (238 B/LB; 701 genes; 150 submitters), giving PP3 no credit
was clearly worse calibrated (Δ log-loss vs cap +0.038 [0.024, 0.050]): PP3 beside PM1 is informative. No rule was better
calibrated than the cap: learned credit −0.004 [−0.010, +0.011], discount −0.004 [−0.009, +0.006], Supporting limit
+0.000 [−0.001, +0.001]. H1a was therefore not met and the cap is retained. With the ACMG convention for an absent code
(deviation 1), learned credit for Supporting PP3 beside PM1 was 0.83 points, the same verdict followed, and a calibrated
fractional credit would move 1,032 laboratory classifications that sit exactly on a class boundary (Table 3).

### 4. PP2 beside PM1 earns less than its point (H2, confirmed)
In the mixed model, PP2 earned 0.35 points [95% CrI −0.12, 0.81] beside PM1 (one-sided posterior probability of ≥ 1
point 0.003; passes Holm). PP2 without PM1 earned 1.12 [0.90, 1.31] — calibrated alone — and −0.23 [−0.67, 0.37]
beside PM1 (ACMG convention). Within submitters PP2 earned 1.08 points [0.36, 1.69] less than PP3 beside PM1 (H2c,
descriptive). The pre-registered mechanism test (H2b) could not be run as intended: only 17 of 2,812 laboratory PM1
calls documented a constraint-based rationale, so whether PP2 is redundant specifically with constraint-based PM1 is
untested.

### 5. Negative controls: classification-derived truth does not measure dependence (H4)
R(PP3, PM1) was below 1 in every evaluation era (≤ 2019 0.27; 2020–2022 0.10; 2023–2024 0.24; 2025–2026 0.17; all upper
bounds < 1; only the last era had ≥ 20 B/LB in the joint cell). But the negative-control criteria, which should not
depend on PM1, behaved the same way: on truth that also excluded them, R(PP1, PM1) was 0.04 [0.02, 0.22] and R(PP4, PM1)
0.03 [0.02, 0.12], and both lost most of their credit beside PM1 (PP1 5.40 → 1.34; PP4 5.00 → 0.76 points). Low R on
classification-derived truth is therefore not specific to dependence: it arises when benign classifications rarely
carry PM1 together with a strong criterion. We consequently treat laboratory and expert-panel estimates of R and credit
as descriptions of practice, not as measurements of the dependence; the de novo scale, where truth is independent of the
codes, carries the measurement. The PP2 result (section 4) is a comparison of credit with and without PM1 and within
submitters; it is not free of this concern and is presented with it.

### 6. Calibration does not remove the dependence (H3)
With AlphaMissense thresholds recalibrated per gene (341 genes), R(PP3, PM1) on laboratory truth was 0.41 [0.30, 0.61]
(published thresholds 0.31). Recalibration reduced but did not remove the dependence, so it is not an artefact of global
thresholds. (Both estimates share the truth construction of section 5; their ratio is the informative quantity, 1.33
[0.65, 1.96].)

### 7. How laboratories combine the criteria, and whether the analysis reflects their practice
Laboratories applied PP3 beside PM1 in 4,471 submissions, 98% at Supporting or Moderate, where the cap does not act;
PM1 was Moderate in 90%; the cap binds in 211 (Table 3). Laboratories reported the uncapped class in 32 of the 38
submissions whose class the cap would change, so the re-scoring baseline reflects reported practice. Eighteen templated
submitters (1,681 submissions) did not drive any estimate (without them, R 0.11, PP2 credit beside PM1 −0.80).

**Table 3. Practice and the classifications each rule would change** (laboratories; 87,873 submissions).

| | n |
|---|---|
| PP3 beside PM1 (Supporting / Moderate / Strong) | 4,471 (4,207 / 184 / 80) |
| PP2 beside PM1; all three | 1,481; 1,035 |
| ClinGen cap binds (PP3 + PM1 > 4) | 211 |
| Changed by: ClinGen cap | 40 |
| Changed by: PP3 + PM1 + PP2 ≤ 4 | 63 |
| Changed by: PP2 not added when PM1 applies | 352 |
| Changed by: PP3 limited to Supporting beside PM1 | 90 |
| Boundary-sitting classifications moved by calibrated fractional PP3 credit | 1,032 |

---

## Discussion

**Stacking related evidence over-counts it, by enough to change a class.** On a label-free scale in three cohorts,
computational evidence was not independent of the location and gene evidence ACMG/AMP credits separately, and the
combination laboratories use most (PP3 + PM1 + PP2) over-counts by two to three points — the width of a classification
step. Uncapped, it puts variants at Pathogenic in a class whose disease association falls well short of 99%.

**What to do.** ClinGen's cap should stay: no rule we tested, including rules fitted to laboratory data, calibrated PP3
beside PM1 better out of sample, and gene-level recalibration did not remove the dependence. The gap is PP2: it is
calibrated alone but adds little beside PM1, in a pre-registered model accounting for submitter and gene. Not adding
PP2 when PM1 applies, as UK guidance already advises (Durkie et al. 2024), would change about 0.4% of laboratory
classifications, mostly from Likely Pathogenic to VUS. Because laboratories rarely document why PM1 applies, we cannot
yet say whether the redundancy is confined to constraint-based PM1; specifications should record the PM1 rationale.

**Version 4.** The forthcoming standard groups related evidence hierarchically and caps it, reportedly separating
predictive from observational evidence [verify against the v4 specification]. Our measurements support placing PP3,
PM1 and PP2 within a shared, capped group: their dependence is strongest exactly where version 2015 adds them freely.
If version 4 retires PP2 or folds it into location evidence, the PP2 finding becomes evidence for that design rather
than a patch to the 2015 rules.

**What classification-derived truth can and cannot show.** Negative-control criteria showed the same apparent loss of
credit beside PM1 as PP3 and PP2. Studies that measure dependence or credit on truth re-scored from curators' own
classifications — including our earlier drafts — should include such controls; without them, low R can reflect which
combinations benign classifications contain rather than overlapping information.

**Limitations.** (1) The de novo scale reflects high-penetrance variants in severe early-onset disorders and is
compressed relative to the clinical scale; class-level comparisons show direction, not exact posteriors. (2) PM1-like
and PP2-like are constraint proxies; panel definitions of PM1 vary. (3) Laboratory and panel truth are derived from the
same classifications (section 5). (4) Several estimates rest on few benign variants and vary between submitters. (5) The
parser was validated against structured codes, not double-annotated free text; a stratified 400-comment sheet is
provided for independent annotation (Supplementary Note 9). (6) The de novo analyses were not independently registered.
(7) Transcripts match MANE Select in 74.4% of genes; cohorts are predominantly of European ancestry.

**Conclusion.** PP3, PM1 and PP2 overlap; stacked, they move variants a class. The cap on PP3 + PM1 is the best rule we
could find and should be kept; PP2 should not be added beside PM1 without documented independence.

---

## Figures and tables
- **Figure 1.** Study design (`fig1_design.png`).
- **Figure 2.** Dependence on the de novo scale by predictor, flag and cohort (`fig2_dependence.png`).
- **Figure 3.** Class-level calibration of stacked evidence (`figures/fig_cost_calibration.png`).
- **Tables 1–3** above; pre-registered results in Supplementary Table S0 and Notes 8–9.

## Data and code availability
All data are public (de novo tables, ClinVar releases 2020-01 to 2026-10, ClinGen Evidence Repository, ProteinGym,
dbNSFP [academic licence; scores not redistributed], AlphaMissense, AlphaFold DB, gnomAD, DDG2P, Roulette). Code,
analysis log, pre-registrations and outputs: study repository [URL; Zenodo DOI at submission]. `REPRODUCE.md` maps every
number, table and figure to the script that produces it.

## Competing interests
[To be completed.]

## Funding
[To be completed.]

## Author contributions
[To be completed.]

---

## References
- ACMG/AMP/CAP/ClinGen. Standards for sequence variant classification, version 4 (draft and pilot materials). [Cite publication or pilot report once available; verify claims marked in the text.]
- Badonyi M, Marsh JA. acmgscaler: an R package and Colab for standardized gene-level variant effect score calibration within the ACMG/AMP framework. 2025.
- Bergquist T, Stenton SL, Nadeau EAW, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 2025;27:101402.
- Carter H, Douville C, Stenson PD, Cooper DN, Karchin R. Identifying Mendelian disease genes with the variant effect scoring tool. *BMC Genomics* 2013;14(Suppl 3):S3.
- Chao KR, Wang L, Panchal R, et al. The landscape of regional missense mutational intolerance quantified from 125,748 exomes. bioRxiv 2024. doi:10.1101/2024.04.11.588920.
- Chen Y, et al. Gene- and domain-aware calibration increases the clinical utility of variant effect predictors. bioRxiv 2026. doi:10.1101/2026.02.17.706269.
- Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 2023;381:eadg7492.
- Durkie M, Cassidy E-J, Berry I, et al. ACGS best practice guidelines for variant classification in rare disease 2024. Association for Clinical Genomic Science, 2024.
- Garrett A, Durkie M, Callaway A, et al. Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. *J Med Genet* 2021;58:297–304.
- Grimm DG, Azencott CA, Aicheler F, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat* 2015;36:513–523.
- Harrizi S, Nait Irahal I, Mostafa K, Arnoult D. Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. bioRxiv 2026. doi:10.64898/2026.08.21.746180.
- Havrilla JM, Pedersen BS, Layer RM, Quinlan AR. A map of constrained coding regions in the human genome. *Nat Genet* 2019;51:88–95.
- Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *Am J Hum Genet* 2016;99:877–885.
- Jin SC, Homsy J, Zaidi S, et al. Contribution of rare inherited and de novo variants in 2,871 congenital heart disease probands. *Nat Genet* 2017;49:1593–1601.
- Kaplanis J, Samocha KE, Wiel L, et al. Evidence for 28 genetic disorders discovered by combining healthcare and research data. *Nature* 2020;586:757–762.
- Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 2020;581:434–443.
- Liu X, Li C, Mou C, Dong Y, Tu Y. dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs. *Genome Med* 2020;12:103.
- Notin P, Kollasch AW, Ritter D, et al. ProteinGym: large-scale benchmarks for protein fitness prediction and design. *Adv Neural Inf Process Syst* 2023.
- Pejaver V, Byrne AB, Feng BJ, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet* 2022;109:2163–2177.
- Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med* 2015;17:405–424.
- Seplyarskiy V, Koch EM, Lee DJ, et al. A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription. *Nat Genet* 2023.
- Stenton SL, Pejaver V, Bergquist T, et al. Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med* 2024;26:101213.
- Tavtigian SV, Greenblatt MS, Harrison SM, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med* 2018;20:1054–1060.
- Tavtigian SV, Harrison SM, Boucher KM, Biesecker LG. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Hum Mutat* 2020;41:1734–1737.
- Tejura M, Fayer S, McEwen AE, et al. Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes. *Am J Hum Genet* 2024.
- Thormann A, Halachev M, McLaren W, et al. Flexible and scalable diagnostic filtering of genomic variants using G2P with Ensembl VEP. *Nat Commun* 2019;10:2373.
- Turner TN, Yi Q, Krumm N, et al. denovo-db: a compendium of human de novo variants. *Nucleic Acids Res* 2017;45:D804–D811.
- Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Res* 2022;50:D439–D444.
- Zhou X, Feliciano P, Shu C, et al. Integrating de novo and inherited variants in 42,607 autism cases identifies mutations in new moderate-risk genes. *Nat Genet* 2022;54:1305–1319.
