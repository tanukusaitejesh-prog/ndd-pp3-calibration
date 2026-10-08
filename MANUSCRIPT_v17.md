# Dependence between computational and location evidence in missense variant classification: a label-free measurement, its cost, and a pre-registered test of practice

*Draft 17, 8 October 2026. Supplementary material: `SUPPLEMENT_v17.md`. Pre-registration of the draft-17 tests: `prereg_fullstudy/A19_V14_PREREG.md` (committed before analysis, commit ceb1a11); deviations: `prereg_fullstudy/A19_DEVIATIONS.md`.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** ACMG/AMP sums computational (PP3), hotspot/domain (PM1) and gene-level (PP2) evidence as though independent; ClinGen limits PP3 with PM1 only at Strong. We measured their dependence without variant-level labels, its consequence, and how practice combines them.

**Methods.** Likelihood ratios were estimated from de novo missense excess over mutation-rate expectation in 47,630 trios and probands, and on functional assays. Practice was examined on 87,873 laboratory and 6,177 expert-panel classifications; four tests on laboratory data across ClinVar releases 2020–2026 were pre-registered, with negative-control criteria. Other practice analyses are exploratory.

**Results.** Joint PP3 and PM1-like evidence was 0.44 (95% CI 0.41–0.47) of the product of their separate evidence in developmental disorders, replicated in two cohorts; beside PM1, PP3 kept 43–49% of its evidence, and 40% on functional assays. Adding PP2 over-counted 2.3 points; uncapped, 2,438 de novo variants were Pathogenic in a class 84% disease-associated against a nominal 99%. In pre-registered tests, no rule calibrated PP3 beside PM1 better than ClinGen's cap, and PP2 earned 0.35 points (−0.12 to 0.81) beside PM1 although calibrated alone. Negative-control criteria showed the same apparent dependence on classification-derived truth. Moderate or Strong PP3 beside PM1 rose to 8.5% of laboratory PM1 submissions in 2025.

**Conclusion.** Computational and location evidence overlap enough to move variants a class. The cap should stay, PP2 beside PM1 over-counts, and stronger PP3 beside PM1 is spreading outside the limit.

---

## Introduction

Missense variants are the largest class of variants of uncertain significance. Under the ACMG/AMP guidelines (Richards et al. 2015), their classification draws on three kinds of evidence: a computational prediction (PP3), location in a mutational hotspot or critical functional domain (PM1), and membership of a gene with a low rate of benign missense variation (PP2). The Bayesian points framework adds these in log-likelihood-ratio units, which assumes that they carry independent information (Tavtigian et al. 2018; Tavtigian et al. 2020).

There are reasons to doubt the assumption. Predictors such as AlphaMissense, REVEL and VEST4 (Cheng et al. 2023; Ioannidis et al. 2016; Carter et al. 2013) learn which genes and protein regions tolerate variation, and that is the information PM1 and PP2 encode. Guidance has responded inconsistently. ClinGen limits the combination of PP3 and PM1 to Strong (Pejaver et al. 2022; Stenton et al. 2024); UK guidance advises against counting constraint evidence in both PM1 and PP2 (Durkie et al. 2024); and CanVIG-UK permits PP3 and PM1 together (Garrett et al. 2021). None of these positions rests on a measurement of the dependence, and it is not known how often, or at what strength, curators combine the criteria.

Measuring the dependence requires evidence that the classifications themselves did not produce. ClinVar labels partly reflect PP3 and PM1, and several predictors were trained on them (Grimm et al. 2015; Harrizi et al. 2026). De novo variants avoid this circularity: their excess over mutation-rate expectation in affected children defines a likelihood ratio without variant-level labels, for any set of sites and within any gene.

We used this property in three steps. We first measured the dependence between PP3, PM1 and PP2 on the de novo scale. We then examined it on clinical and functional truth independent of the de novo data. Finally, we measured how much PP2 and PP3 add once PM1 applies in expert-panel and laboratory practice, and how many classifications would change if the overlaps were removed.

---

## Methods

### Analysis log
Analyses followed a public, version-controlled analysis log of seventeen amendments, each committed before the analysis it governs (Supplementary Note 1); the log was not deposited in an independent registry beforehand. The primary pre-specified criterion concerned how much predictor evidence survives conditioning on gene and protein position. It was not met under its original position control, which did not control position in genes without a regional-constraint split, and a fixed-window control was adopted afterwards. The PP3 × PM1 dependence test was specified before its results were seen. All analyses of expert-panel and laboratory practice were added later and are hypothesis-generating, except four tests on laboratory data (H1–H4 below), which were pre-registered in a document committed to the repository before any of them was run (A19). Supplementary Table S0 lists every hypothesis and when it was specified relative to data access; one deviation, made after running H1, is reported beside the pre-registered result.

### Data
**De novo cohorts.** We used three cohorts (Figure 1): 31,058 developmental-disorder (DD) trios (Kaplanis et al. 2020); 12,316 autism probands, with 4,051 unaffected siblings as a negative control (Zhou et al. 2022); and a public cohort of 4,256 trios with congenital heart disease, schizophrenia, Tourette syndrome, epileptic encephalopathy or intellectual disability (Jin et al. 2017; denovo-db v1.6.1, Turner et al. 2017). Studies known to overlap the DD or autism cohorts were excluded from the public cohort. Individuals sharing two or more identical de novo variants with a DD proband were removed (64 autism probands and 8 public trios), and the public cohort was also deduplicated against the autism cohort. Unrecognised overlap cannot be excluded entirely.

**Sites and annotations.** Sites were coding positions on GENCODE v32 canonical transcripts, excluding those within 3 bp of an exon boundary. Predictor scores came from dbNSFP 4.1a and 5.4a (Liu et al. 2020) and AlphaMissense; gene and regional missense constraint from gnomAD (Karczewski et al. 2020; Chao et al. 2024); constrained coding regions from Havrilla et al. (2019); disease genes from the DDG2P monoallelic list (Thormann et al. 2019); and structure from AlphaFold DB (Varadi et al. 2022).

**Expert-panel classifications.** We used 6,177 missense classifications by 46 expert panels from the ClinGen Evidence Repository (export of 7 October 2026).

**Laboratory classifications.** We used 87,873 germline ClinVar submissions (October 2026) for missense SNVs, from 472 submitters, whose comments list ACMG/AMP codes; expert-panel and practice-guideline submissions were excluded. A parser read codes and strengths from the free text and skipped negated codes. We validated it on the free-text summaries of 13,091 Evidence Repository classifications against their structured codes: overall precision was 0.92 and recall 0.98, and for PP3, PM1 and PP2 precision was 0.95, 0.98 and 0.98, with recall and agreement on strength of at least 0.99.

**Functional assays.** We used human single-substitution assays in ProteinGym v1.3 (Notin et al. 2023) that mapped to analysed transcripts: 82 assays in 70 genes. Damaging variants were those with `DMS_score_bin` = 0.

### Dependence on the de novo scale
Expected de novo counts E were derived from Roulette mutation rates (Seplyarskiy et al. 2023) and scaled to each cohort by its synonymous de novo count. The excess is X = O − E, and for a set of sites b the likelihood ratio is LR = (X_b/ΣX)/(E_b/ΣE).

PP3 was flagged at the ClinGen thresholds for AlphaMissense (primary), REVEL and VEST4 (Pejaver et al. 2022; Bergquist et al. 2025). A site was PM1-like if its regional missense observed/expected ratio was below 0.4 or its constrained-coding-region percentile was at least 95, and PM1-strict if observed/expected was below 0.2, the level Chao et al. (2024) calibrated to Moderate. A gene was PP2-like if its missense Z score was at least 3.09. This flag captures only the constraint component of PP2, not the requirement that missense variants are a common disease mechanism.

We used two measures of dependence. R = LR(A ∧ B)/[LR(A) · LR(B)] is the joint evidence of criteria A and B relative to the product of their separate evidence: R = 1 means independence, and R < 1 means that summing them overstates the joint evidence. f = ln LR(A ∧ B)/[ln LR(A) + ln LR(B)] is the fraction of their summed log-evidence that is kept. The two describe the same dependence on different scales; f is milder because it includes the evidence of the stronger criterion, which is largely kept.

To compare PP3's evidence on a single scale, we report its points alone, ln LR(PP3)/ln 2.406, and beside PM1, ln LR(PP3 | PM1)/ln 2.406, together with the share kept. We do not compare de novo points with nominal ACMG weights, because the de novo scale is compressed relative to the clinical scale. Intervals resample genes (1,000 replicates).

### Dependence and credit on clinical truth
Expert-panel and laboratory classifications were re-scored after removing the computational and location criteria (PP3, BP4, PM1, PP2, BP1, BP3, PM4, PM5); a sensitivity analysis also removed the de novo criteria (PS2, PM6) and PS1, PP5 and BP6. Truth was pathogenic or likely pathogenic (P/LP) at 6 or more remaining points, and benign or likely benign (B/LB) at −1 or fewer. On this truth, the likelihood ratio of a criterion C is LR(C) = P(C applied | P/LP)/P(C applied | B/LB), and R is defined as above.

The credit that C earns once PM1 applies is ln[LR(C ∧ PM1)/LR(PM1)]/ln 2.406, on the ClinGen points scale, and is compared with the points curators applied. We report credit as well as R because classification uses credit: R can be well below 1 while credit equals the nominal points, if C alone carries more than its nominal points.

Laboratories differ in how liberally they apply criteria and in how they classify, so a pooled estimate can reflect differences between laboratories. We therefore also estimated credit within submitters: among PM1-applied submissions, the likelihood ratio of C for P/LP against B/LB was pooled across submitters with the Mantel–Haenszel risk-ratio estimator, so that only contrasts within a laboratory contribute, with a submitter-clustered bootstrap. Because laboratory benign truth rests largely on frequency criteria, PP2 credit was also estimated in three further ways: in submissions without BA1, BS1 or BS2; in submissions with PM2; and with benign truth requiring at least one benign criterion other than frequency. Sensitivity analyses varied the truth thresholds (P/LP ≥ 4, 6 or 8; B/LB ≤ −1 or −3), removed the frequency criteria (BA1, BS1, BS2, PM2), removed the de novo criteria with PS1, PP5 and BP6, and clustered by submitter as well as by gene. For replication, laboratory data were split into independent halves by date of last evaluation and by submitter, and into 100 random submitter halves. For Supporting PP3, we tested equivalence by asking whether the 90% interval of its credit lay within ±0.5 points of 1. In expert-panel data we also used a temporal split (classifications approved before 2025 against 2025–2026), defined PP3 by predictor score rather than by curator choice, and compared Moderate with Supporting PP3 directly, with its minimum detectable difference.

### ClinVar, assays and classification rules
In ClinVar we compared P/LP with B/LB variants after removing DD de novo sites. In functional assays, R was computed within each assay with a position-clustered bootstrap and pooled by random effects, for assays in which 5–95% of variants lay at PM1-like sites. Because the PM1-like flag is a constraint proxy, we repeated the assay analysis with curated PM1: the protein residues at which any expert panel applied PM1 to a missense variant in the Evidence Repository (841 residues in 67 genes). Genes were used if at least 90% of these residues matched the reference amino acid of the analysed transcript, and assays if at least 10 assayed variants lay at curated residues. Finally, classifications were re-scored on the points scale (P ≥ 10, LP 6–9, VUS 0–5, LB −1 to −6, B ≤ −7) under four rules: the ClinGen limit; a limit on PP3 + PM1 + PP2; PP2 not added when PM1 applies; and PP3 limited to Supporting when PM1 applies. To describe trends, co-application among laboratory PM1 submissions was tabulated by year of last evaluation.

### Held-out test of the PP2 rule
To test the rule "PP2 not added when PM1 applies" on outcomes not used to derive it, we parsed ClinVar laboratory missense submissions from January 2023 as above. Rule-affected submissions were Likely Pathogenic submissions whose class would fall to VUS if PP2 were removed beside PM1. There were two control sets. The first was Likely Pathogenic submissions on the same total points (6) without PP2 beside PM1. The second was the subset of these with PM1, which differ from rule-affected submissions only in whether their sixth point came from PP2 or from other evidence. Variants used to estimate PP2 credit (2026 PM1 submissions with truth) were excluded. The primary outcome, read from ClinVar in October 2026, was whether most of the other submitters' classifications of the variant were VUS, likely benign or benign. The secondary outcome was whether the same submission was no longer P/LP. Groups were compared with the exact conditional odds ratio, its exact central interval and the matching exact p value, and with Mantel–Haenszel odds ratios stratified by gene and by submitter.

### The cost of stacking
DD de novo missense variants were classified with PS2 + PM2 + PM1-like + AlphaMissense PP3 under uncapped, capped and conditional rules. Each class's disease-associated fraction is 1 − E/O (de novo excess), with an exact Poisson interval on O; it is a class-level estimate of association with any severe developmental disorder, compared with the class's nominal posterior (Pathogenic ≥ 0.99, Likely Pathogenic ≥ 0.90; prior 0.0441). Points over-counted by dependence are −ln R / ln 2.406.

### Pre-registered tests (A19)
Laboratory submissions from four ClinVar releases (January 2020, 2022 and 2024; October 2026) were parsed identically, each submission counted once. Credit followed the ACMG convention that an absent code carries no evidence. **H1 (primary)** compared six rules for PP3 beside PM1 (uncapped; ClinGen cap; PP3 at most Supporting; a discount; a learned credit per strength; no credit) by held-out log-loss of laboratory truth among PM1-applied submissions, with gene-grouped five-fold cross-validation repeated 20 times and submitter-clustered intervals. **H2 (primary)** estimated PP2 credit beside PM1 in a Bayesian mixed logistic model with submitter and gene random intercepts; a deficit was declared if the 95% credible interval excluded 1 point, and H2b classified the documented PM1 rationale. **H3** recalibrated AlphaMissense PP3 thresholds per gene on ClinVar, leaving out the evaluation variants. **H4** repeated R(PP3, PM1) and PP2 credit in four eras of last evaluation and, as negative controls, estimated R and credit beside PM1 for PP1 and PP4 on truth that also excluded them. Holm correction was applied across H1 and H2.

---

## Results

### PP3, PM1 and PP2 are not independent on the de novo scale
For AlphaMissense PP3 at Supporting with a PM1-like flag, R was 0.44 [0.41, 0.47] in DD, and the dependence replicated in autism (0.56) and in the public cohort (0.67) (Table 1, Figure 2). It persisted when variants were compared within DDG2P monoallelic genes (0.69–0.75). Across all genes, however, the within-gene ratio was 0.75 in DD and 0.96 in autism and the public cohort. Much of the dependence therefore lies between genes, and within genes it is concentrated in known disease genes.

Dependence was weaker for REVEL (0.57) and VEST4 (0.67), and across 40 informative predictors in DD every R with PM1-like had an upper bound below 1. It held in ordered and disordered regions, at buried and exposed sites, and after removing the 50 genes with the largest excess (Supplementary Note 3). PM1 and PP2 were the most dependent pair, largely by construction, because both of our flags derive from missense constraint.

**Table 1. Dependence on the de novo scale** (R; gene-clustered 95% intervals).

| Criteria | DD | Autism | Public |
|---|---|---|---|
| PP3 (AlphaMissense) × PM1-like | 0.44 [0.41, 0.47] | 0.56 [0.47, 0.66] | 0.67 [0.48, 0.96] |
| — within DDG2P genes | 0.69 [0.65, 0.73] | 0.74 [0.67, 0.82] | 0.75 [0.62, 0.89] |
| — within genes, all genes | 0.75 [0.71, 0.78] | 0.96 [0.84, 1.10] | 0.96 [0.74, 1.28] |
| PP3 (AlphaMissense) × PP2-like | 0.61 [0.57, 0.65] | 0.65 [0.55, 0.76] | 0.78 [0.57, 1.07] |
| PM1-like × PP2-like | 0.31 [0.28, 0.33] | 0.36 [0.30, 0.45] | — |

Beside PM1, AlphaMissense kept 43–49% of its own evidence and REVEL 61–78%, with similar shares at Supporting and at Strong (Table 2). The fraction of combined log-evidence kept, f, was 0.73–0.77 for AlphaMissense in DD, and 0.72 [0.70, 0.74] at Strong with PM1-strict. Because the de novo scale is compressed relative to the clinical scale, we do not compare these values with the nominal weights behind ClinGen's limit.

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

The overlap arises because missense predictors largely carry location information. Across 39 informative predictors, a median 0.32 of the log-evidence of their highest-scoring 2% was incremental to gene and 50-residue window, against 0.51 for a random partition of genes into the same number of units. Inside constrained regions AlphaMissense flags 51% of mutation-expected sites, against 19% overall. A perfect predictor would also lose evidence within windows, by an amount the data cannot identify, so this share measures evidence incremental to location rather than predictor skill (Supplementary Note 2).

### The cost of stacking
In points, Supporting PP3 beside a PM1-like flag over-counts 0.97 in DD, and the triple PP3 × PM1-like × PP2-like (R3 0.13 [0.11, 0.14]) over-counts 2.32 points [2.24, 2.51] in DD and 1.95 [1.63, 2.32] in autism, about the width of a classification step; 1,035 of the 1,481 laboratory submissions that combine PM1 and PP2 also apply PP3. Classified with PS2 + PM2 + PM1-like + AlphaMissense PP3 and no limit, 2,438 DD de novo variants were Pathogenic, a class the de novo excess implies is 84% [84, 85] disease-associated against a nominal 99% (Figure 3). Under ClinGen's limit the same variants are Likely Pathogenic. The limit does not make that class calibrated (48% against a nominal 90%), and part of the shortfall reflects the estimand, which is compressed relative to the clinical scale; the class comparison shows direction, the points over-counted show size. In ClinVar, PP3 with Moderate PM1 alone would make 7,101 PP3-positive VUS in DDG2P monoallelic genes Likely Pathogenic without the limit, and none under any limiting rule we tested (Supplementary Note 9).

### The dependence on independent truth
Every estimate of R on clinical truth was below 1 (Table 3). On laboratory truth, clustering by submitter widened the intervals without changing any conclusion (PP3 × PM1, 0.18 [0.08, 0.43]), and PP3 × PM1 and PM1 × PP2 were below 1, with intervals excluding 1, in every independent half of the laboratory data (Supplementary Note 6) and, in the pre-registered replication (H4), in every era of last evaluation from before 2020 to 2026 (0.10–0.27). Pre-registered negative controls, however, showed that low R on classification-derived truth is not specific to dependence. PP1 and PP4, which carry no information about location, had R with PM1 of 0.04 [0.02, 0.22] and 0.03 [0.02, 0.12] on truth that also excluded them, lower than PP3's. Low R on such truth therefore arises largely from which combinations benign classifications contain, and we treat the expert-panel and laboratory values in Table 3 as descriptions of practice rather than measurements of dependence. The measurements are the de novo scale and the functional assays, on which truth does not derive from the codes.

**Table 3. Dependence on truth independent of the de novo data.**

| Truth | Criteria | Scope | R [95% CI] |
|---|---|---|---|
| ClinVar P/LP vs B/LB* | AlphaMissense PP3 × PM1-like | DDG2P genes | 0.32 [0.28, 0.37] |
| Expert-panel truth | Panel PP3 × panel PM1 | 2,942 classifications, 130 genes | 0.21 [0.12, 0.31] |
| Laboratory truth | Laboratory PP3 × PM1 | 39,015 submissions, 5,976 genes | 0.18 [0.14, 0.25] |
| Laboratory truth | Laboratory PM1 × PP2 | same | 0.27 [0.19, 0.44] |
| Laboratory truth | Laboratory PP3 × PP2 | same | 0.42 [0.35, 0.51] |
| Functional assays | AlphaMissense PP3 × curated PM1 | 7 assays, 3 genes (TP53, GCK, HRAS) | 0.49 [0.40, 0.60] |
| Functional assays | AlphaMissense PP3 × PM1-like | 13 assays, 10 genes outside DDG2P | 0.94 [0.88, 0.99] |

\*Partly circular, because ClinVar classifications use PP3 and PM1 themselves. This biases R towards 1.

Functional assays are the one truth set that curators did not produce. On them, the quantity that matters for classification is what PP3 adds at residues where PM1 applies. At curated PM1 residues, AlphaMissense PP3 had a median likelihood ratio for damage of 1.78, against 3.41 across each assay, so it kept a median 40% of its log-evidence (Supplementary Note 8). This is close to the 43–49% it keeps beside PM1-like sites on the de novo scale. It happens because the predictor already scores most variants at these residues as damaging: 77% were PP3-positive there, against 34% elsewhere. R was correspondingly below 1 in all seven assays (0.49 [0.40, 0.60] at Supporting; 0.43 [0.36, 0.51] at Strong). This included two assays in which curated PM1 itself carried no information about damage (GCK abundance and HRAS; LR 1.00 and 0.85). There, R < 1 reflects PP3 adding nothing at those residues rather than PM1 duplicating PP3, which for classification is the same loss of credit. The constraint proxy for PM1, by contrast, barely predicted damage within assays (median LR 1.10 outside DDG2P genes), and its dependence was weak (0.94). Coverage also differed. The proxy flags none of the assayed TP53 residues, 3% in BRCA1 and nearly all in PTEN, so none of these genes could be tested with it. Curated PM1 made TP53 testable; PTEN failed the numbering check, and BRCA1 has no curated PM1 residues.

### How panels and laboratories combine the criteria
Panels applied PP3 beside PM1 in 787 classifications and laboratories in 4,471 (Table 4). In both, PP3 was almost always at Supporting or Moderate strength, where ClinGen's limit does not act. Panels applied PP2 beside PM1 in 174 classifications, 145 of them from three panels whose specifications may permit it (RASopathy, FBN1 and Monogenic Diabetes); laboratories did so in 1,481 submissions.

**Table 4. Co-application in practice.**

| | Expert panels | Laboratories |
|---|---|---|
| Missense classifications | 6,177 | 87,873 |
| PP3 beside PM1 | 787 | 4,471 |
| — PP3 Supporting / Moderate / Strong | 605 / 178 / 4 | 4,207 / 184 / 80 |
| Share at Strong, where the limit acts | 0.5% | 1.8% |
| PP2 beside PM1 | 174 | 1,481 |
| All three | 138 | 1,035 |

Combinations at higher strength are increasing. Among laboratory PM1 submissions, PP3 at Moderate or Strong was applied in 0–0.9% each year from 2018 to 2023, then in 4.3% in 2024 and 8.5% in 2025. This followed the introduction of calibrated PP3 strengths (Pejaver et al. 2022). PP2 beside PM1 rose from 18–20% in 2022–2023 to 24% in 2025 and 31% in 2026 (both trends p < 10⁻⁸; Supplementary Note 8).

### The evidence PP2 and PP3 add once PM1 applies
Table 5 gives the points earned once PM1 applies, against the points applied.

In laboratory practice, PP2 beside PM1 earned −0.22 points [−0.67, 0.40] against the 1 point applied (565 submissions; 57 B/LB in the joint cell), and −0.41 [−0.89, 0.24] with stricter truth. The point estimate was below the nominal point in every independent half, from −1.23 to 0.14 by date of evaluation and by submitter; the interval excluded 1 in two of the four halves and in 62% of 100 random submitter halves. Part of this deficit reflects differences between laboratories: clustering by submitter widens the interval to [−1.23, 0.79], and within submitters PP2 earned 0.37 points [−0.36, 0.95] (35 laboratories with both P/LP and B/LB PM1 submissions; 432 with PP2, 44 of them B/LB). The within-submitter estimate is the less confounded of the two. It is higher than the pooled estimate, and its interval excludes the nominal point only narrowly. Laboratory benign truth rests mostly on frequency criteria: only 43 of 238 B/LB PM1 submissions lacked BA1, BS1 and BS2. Handling frequency gave unstable estimates whose intervals included 1. Within submitters, PP2 earned 0.90 [−0.49, 2.18] without BA1/BS1/BS2, −0.42 [−0.90, 0.63] among PM2-applied submissions, and 0.16 [−0.71, 1.14] with benign truth requiring a non-frequency criterion. A deficit for PP2 beside PM1 is therefore likely but not established independently of frequency-based benign truth. Expert-panel truth could not test this pair (1 B/LB). In the pre-registered mixed model with submitter and gene random effects (H2), PP2 earned 0.35 points [95% credible interval −0.12, 0.81] beside PM1, below its nominal point (posterior probability of at least 1 point 0.003; passes Holm). Alone, PP2 was calibrated (1.12 [0.90, 1.31]), and within submitters it earned 1.08 points [0.36, 1.69] less than PP3 beside PM1. The negative controls qualify this too: on the same truth PP1 and PP4 also lost most of their credit beside PM1 (5.40 to 1.34 and 5.00 to 0.76 points), so a loss of credit beside PM1 on classification-derived truth is not by itself evidence of redundancy. The mechanism test (H2b) could not be run: only 17 of 2,812 laboratory PM1 calls documented a constraint-based rationale.

Supporting PP3 beside PM1 earned 0.81 points [0.24, 1.68] in panels and 0.84 [0.52, 1.26] in laboratories. Within submitters it earned 0.93 [−0.03, 1.62]. It was equivalent to its nominal point within ±0.5 points with gene clustering (90% interval 0.58–1.17) but not with submitter clustering (0.10–1.47), and varied from 0.24 to 1.21 across independent halves. Two cautions apply. Panels apply PP3 preferentially to more convincing variants: within panel-PM1 classifications, the likelihood ratio for pathogenicity was 2.16 where PP3 was applied and 1.05 where the variant was eligible by score but PP3 was not applied. And defining PP3 by score rather than by panel choice gave 0.68 points [0.34, 1.19] for REVEL and 1.66 [0.98, 2.58] for AlphaMissense. We find no consistent evidence of over-credit at Supporting (a difference of 0.5 points is not excluded), but cannot establish that it is calibrated.

Moderate PP3 beside PM1 earned 0.02 points [−1.69, 0.89] against 2 in expert panels, on 6 B/LB variants. It could not be distinguished from Supporting (difference −0.79 [−2.81, 0.27], p = 0.14; minimum detectable difference about 4 points), and laboratories had too few B/LB variants to test it. We regard possible over-credit at Moderate as a hypothesis.

Over all strengths, panel PP3 beside PM1 earned 0.60 points [0.14, 1.05] against 1.24 applied (Figure 3). This was stable in held-out genes (0.59), with PS1, PP5, BP6, PS2 and PM6 removed from truth (0.57 [0.08, 1.04]) and across truth thresholds (0.60–0.65), and lower in classifications approved in 2025–2026 (0.40 [−0.25, 1.07]).

**Table 5. Points earned once PM1 applies, against points applied** (clinical truth; gene-clustered 95% intervals; n and B/LB in the joint cell).

| Criterion | Expert panels | Laboratories | Laboratories, within submitter† | Applied |
|---|---|---|---|---|
| PP2 | not estimable (1 B/LB) | −0.22 [−0.67, 0.40] (565; 57) | 0.37 [−0.36, 0.95] | 1.00 |
| PP3, any strength | 0.60 [0.14, 1.05] (304; 17) | 0.85 [0.54, 1.22] (1,768; 74) | 1.00 [0.05, 1.68] | 1.24 / 1.09 |
| PP3 Supporting | 0.81 [0.24, 1.68] (235; 11) | 0.84 [0.52, 1.26] (1,674; 71) | 0.93 [−0.03, 1.62] | 1.00 |
| PP3 Moderate | 0.02 [−1.69, 0.89] (67; 6) | not estimable (3 B/LB) | not estimable | 2.00 |

†Mantel–Haenszel across 35 laboratories with both P/LP and B/LB PM1 submissions; submitter-clustered 95% intervals.

### Consequences for classification
Rules that remove these overlaps would change few classifications overall (Table 4, continued), but many of those they apply to. Not adding PP2 when PM1 applies would change 352 of 87,873 laboratory classifications (0.4%). That is 24% of the 1,481 classifications that combine the two, and 36% of their Likely Pathogenic calls (223 of 624) would become VUS. The changed classifications span 160 genes and 66 laboratories, led by NOTCH3, FBN1, COL4A3, CASR, POLG, CACNA1A, KCNH2 and SCN1A. Limiting PP3 to Supporting when PM1 applies would change 90, and ClinGen's limit changes 40.

In the pre-registered comparison (H1; 2,812 PM1-applied submissions with truth, 701 genes), giving PP3 no credit beside PM1 was worse calibrated than the cap (Δ log-loss +0.038 [0.024, 0.050]), so PP3 is informative there, but no rule was better calibrated than the cap: a learned credit per strength −0.004 [−0.010, +0.011], a discount −0.004 [−0.009, +0.006], PP3 limited to Supporting +0.000 [−0.001, +0.001]. The cap is retained. Under the ACMG convention for an absent code (deviation 1) the learned credit for Supporting PP3 beside PM1 was 0.83 points, the verdict was the same, and a fractional credit would move 1,032 laboratory classifications that sit exactly on a class boundary. Recalibrating AlphaMissense thresholds per gene (H3; 341 genes) reduced but did not remove the dependence (R 0.41 [0.30, 0.61], against 0.31 with published thresholds). PM1 was Moderate in 90% of laboratory PP3-and-PM1 submissions and the cap binds in 211. Of the 38 submissions applying both whose class the cap would change, laboratories reported the uncapped class in 32, so re-scoring reflects reported practice; 18 heavily templated submitters (1,681 submissions) did not drive any estimate (Supplementary Note 10).

**Table 4 (continued). Classifications that change class under alternative rules.**

| Rule | Panels: changed | P→LP | LP→VUS | Other | Laboratories: changed | P→LP | LP→VUS | Other |
|---|---|---|---|---|---|---|---|---|
| ClinGen limit (PP3 + PM1 ≤ 4) | 24 | 12 | 12 | 0 | 40 | 30 | 10 | 0 |
| PP3 + PM1 + PP2 ≤ 4 | 35 | 17 | 18 | 0 | 63 | 39 | 22 | 2 |
| PP2 not added when PM1 applies | 39 | 28 | 11 | 0 | 352 | 124 | 223 | 5 |
| PP3 limited to Supporting when PM1 applies | 43 | 22 | 21 | 0 | 90 | 39 | 51 | 0 |

"Other" counts any other change of class.

### A held-out test of the PP2 rule
Of 8,510 laboratory missense submissions with ACMG codes in January 2023, 1,748 were Likely Pathogenic. After excluding variants used to estimate PP2 credit, 41 reached that class only through PP2 beside PM1. There were 450 controls on the same points without that combination, 125 of them with PM1. By October 2026, other laboratories classified rule-affected variants as VUS or benign by majority more often than controls: 9 of 34 (26%) against 48 of 411 (12%) (exact odds ratio 2.71 [1.05, 6.46], p = 0.039; stratified by submitter, 4.69 [1.39, 15.85]; by gene, 7.43 [0.98, 56.5]). Against controls that also had PM1, the difference was in the same direction but not significant (26% against 15%; 2.02 [0.71, 5.51], p = 0.21). The submitting laboratory withdrew the pathogenic call for 2 of 41 rule-affected submissions against 2 of 442 controls (11.2 [0.79, 158], p = 0.076). The test supports the direction of the PP2 estimates, but it is small and does not separate PP2 from location evidence in general (Supplementary Note 8).

---

## Discussion

Computational evidence is not independent of the location evidence that ACMG/AMP credits separately. The dependence appears on three kinds of evidence. On a label-free de novo scale it holds in three cohorts and within known disease genes. On expert-panel and laboratory truth it holds in every independent half of the laboratory data. And on functional assays, the one truth that curators did not produce, it holds once PM1 is defined as curators apply it. On clinical truth R was also below 1, but negative-control criteria with no location content behaved the same way, so those values describe practice rather than measure dependence; the functional assays, which curators did not produce, carry that check. The dependence arises because predictors have learned where pathogenic variation falls, which is what PM1 and PP2 describe. Its reach is limited: much of it lies between genes, so that within genes outside the known disease genes it is small; and a constraint proxy for PM1 captures it only where constraint itself predicts damage.

Dependence does not imply a deficit at every strength. Beside PM1, PP3 keeps about half of its own evidence (AlphaMissense) or more (REVEL), but whether that falls short of the points curators give depends on how much PP3 carries on its own. On clinical truth PP3 alone carries more than its nominal point at Supporting, so Supporting PP3 beside PM1 showed no consistent over-credit. ClinGen's limit acts only at Strong, which is 1–2% of practice; whether it is adequate there could not be tested on clinical truth, where too few Strong combinations had benign outcomes. The overlap matters more as stronger PP3 spreads. Beside PM1, PP3 keeps 40–50% of its evidence on two scales that curators did not produce, the de novo and the functional. At Supporting, the clinical evidence PP3 carries alone covers that loss. At Moderate, nominally two points, it may not: on panel truth, Moderate PP3 beside PM1 earned 0.02 points [−1.69, 0.89]. Moderate and Strong PP3 beside PM1 rose from under 1% of laboratory PM1 submissions before 2024 to 8.5% in 2025, most of it at Moderate, where the limit does not act.

Adjustment may be needed in two places. Moderate PP3 beside PM1 may be over-credited. The clinical evidence is thin, but the label-free and functional evidence both show that PP3 keeps only about half of its evidence there, and such combinations are increasing. PP2 beside PM1 earned less than its nominal point in most laboratory analyses, and Likely Pathogenic calls that depended on it were less often confirmed by other laboratories. This is consistent with UK guidance (Durkie et al. 2024), but the deficit could not be separated from frequency-based benign truth, and expert-panel truth could not test it.

**Implications for curation.** Our results support four practical positions. First, Supporting PP3 beside PM1 showed no detectable over-credit; a difference of half a point is not excluded. Second, we propose testing an extension of ClinGen's limit to Moderate PP3: when PM1 applies, PP3 is counted at Supporting. This rule has the most support from evidence that curators did not produce. It would change 90 laboratory classifications (Table 4, continued), a number that will grow as Moderate PP3 spreads. Third, PP2 should not be added when PM1 rests on the same constraint evidence, as UK guidance and the ClinGen cardiomyopathy panel already do. Our data support that advice without establishing a deficit independent of it. In laboratory practice it would change 24% of the classifications that combine the two. Fourth, calibrations of new predictors and new location criteria should report evidence conditional on the other criteria, not only marginal evidence; the de novo method used here requires no variant labels and can be applied to any predictor or region definition with the released code. Studies that estimate dependence or credit on truth re-scored from curators' own classifications, as our earlier analyses did, should include negative-control criteria: without them, low R can reflect which combinations benign classifications contain rather than shared information. ClinGen's limit was not bettered by any rule we tested out of sample, including rules fitted to laboratory data.

**Version 4.** The forthcoming ACMG/AMP/CAP/ClinGen standard reorganises evidence hierarchically, in part to stop related codes being counted twice, and reportedly groups predictive and functional evidence apart from observational evidence [verify against the version 4 specification]. Our measurements support placing PP3, PM1 and PP2 in a shared, capped group. If version 4 retires PP2 or folds it into location evidence, the PP2 findings are evidence for that design rather than a patch to the 2015 rules.

**Relation to previous work.** That PP3, PM1 and PP2 can double-count the same attribute is not new: it motivates ClinGen's limit, UK guidance, Sherloc's rule dependencies and at least one expert-panel specification (Supplementary Table S17). Those positions rest on expert judgement. PP3 thresholds (Pejaver et al. 2022; Bergquist et al. 2025) and regional constraint for PM1 (Chao et al. 2024) are calibrated one criterion at a time, and gene- and domain-aware calibration addresses heterogeneity between genes rather than dependence between criteria (Tejura et al. 2024; Badonyi & Marsh 2025; Chen et al. 2026). Our contribution is empirical: a label-free measurement of the dependence, its examination on independent truth, and an audit of how the criteria are combined and credited in practice. We found no previous study doing any of these (Supplementary Table S17).

**Limitations.** The analysis log was not independently registered and the primary pre-specified criterion was not met. Of the practice analyses, only H1–H4 were pre-registered (in the repository, not an independent registry), and H1's baseline was changed after running it (both versions reported). The dependence test was specified before its results, but the paper's emphasis on it followed the failure of the primary criterion. The analyses added afterwards (curated PM1 on assays, within-submitter estimates and the held-out reclassification test) used data not involved in that choice. The de novo scale reflects high-penetrance variants in severe early-onset disorders, so we compare evidence within that scale rather than against nominal weights. Panel and laboratory truth come from the same classifications whose PP3 and PM1 we evaluate, after those criteria are removed: truth thresholds select on the remaining evidence, curators choose PP3 strength knowing that evidence, and benign truth relies heavily on frequency criteria. Our PM1 and PP2 flags on the de novo scale are constraint-based proxies, and panel definitions of PM1 vary. Curated PM1 could be tested on assays in only three genes, and the residues at which panels applied PM1 are those where pathogenic variants had already been classified. Pre-registered negative controls showed that R and credit on classification-derived truth are not specific to dependence. The parser was validated against structured codes, not double-annotated free text; a stratified 400-comment sheet is provided for independent annotation. Several estimates rest on few benign variants, laboratory estimates vary between submitters, and laboratory benign truth is mostly frequency-based. The held-out test is small, and its comparison with PM1-matched controls was not significant. Functional assays covered ten genes with the proxy and three with curated PM1. Transcripts match MANE Select in 74.4% of genes, and the cohorts are predominantly of European ancestry.

**Conclusion.** PP3, PM1 and PP2 are not independent on de novo, clinical or functional evidence: beside PM1, PP3 keeps about half of its own evidence. In practice, the overlap is concentrated where the existing limit does not act, and stronger PP3 beside PM1 is spreading. Stacked, they over-count by up to a classification step. ClinGen's limit was not bettered by any rule tested and should stay; Supporting PP3 beside PM1 showed no detectable over-credit; Moderate PP3 beside PM1 should be tested for a limit; and PP2 beside PM1 over-counts in a pre-registered model, which supports existing UK advice.

---

## Figure legends
**Figure 1. Study design.** (a) De novo excess over mutation-rate expectation in three cohorts defines label-free likelihood ratios for PP3, PM1-like and PP2-like sites, from which dependence R is estimated marginally and within genes. (b) Dependence is re-examined on expert-panel and laboratory classifications re-scored without computational and location criteria, on ClinVar, and on functional assays with proxy and panel-curated PM1. (c) Practice is audited: how often and at what strength the criteria are combined, the points PP2 and PP3 earn once PM1 applies (pooled and within laboratory), the classifications that change under alternative rules, and a held-out test of one rule against 2026 reclassification (`fig1_design.png`).

**Figure 2. Dependence on the de novo scale.** R for PP3 with PM1-like and PP2-like flags by predictor, flag and cohort, with gene-clustered 95% intervals; R = 1 is independence (`fig2_dependence.png`).

**Figure 3. The cost of stacking.** Disease-associated fraction (de novo excess, DD) of the classes produced by PS2 + PM2 + PM1-like + AlphaMissense PP3 without and with ClinGen's limit, with exact 95% intervals, against each class's nominal posterior band (`figures/fig_cost_calibration.png`). The former Figure 3 (PP3 credit beside PM1) is Supplementary Figure S1; its values are in Table 5.

## Data and code availability
All data are public: published de novo tables, ClinVar, the ClinGen Evidence Repository, ProteinGym, dbNSFP (academic licence), AlphaMissense, AlphaFold DB, gnomAD, DDG2P and Roulette. Code, the analysis log, pre-registrations and all outputs are in the study repository; `REPRODUCE.md` maps every number, table and figure to its script, and `reproduce_v14.sh` re-runs the pre-registered laboratory analyses from public ClinVar files (https://github.com/tanukusaitejesh-prog/ndd-pp3-calibration; to be made public with a Zenodo DOI at submission).

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
