# Dependence between computational and location evidence in missense variant classification, and its consequences in expert and laboratory practice

*Draft 12, 7 October 2026. Supplementary material: `SUPPLEMENT_v12.md`.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** ACMG/AMP classification sums computational evidence (PP3) with hotspot/domain (PM1) and gene-level (PP2) evidence as if they were independent; ClinGen caps PP3 + PM1 only at Strong. We measured their dependence without variant labels, and the evidence each adds beside PM1 in practice.

**Methods.** Likelihood ratios came from the excess of de novo missense variants over mutation-rate expectation in 47,630 trios and probands. Results were checked on ClinVar, functional assays, and 6,177 expert-panel and 87,873 laboratory missense classifications re-scored without computational, location or de novo criteria. Practice analyses were not pre-specified.

**Results.**
- **No pair of criteria was independent.** On the de novo scale, PP3 kept about half its own evidence beside PM1 (AlphaMissense 43%; R 0.44 [0.41–0.47]). On panel and laboratory truth, R was 0.21 and 0.18.
- **PP2 beside PM1 showed no detectable benefit in laboratory practice:** −0.22 points [−0.67–0.40] against 1 applied. Laboratories co-applied the two in 1,481 submissions; not adding PP2 would change 352 classifications.
- **Supporting PP3 beside PM1 was close to its nominal point** on panel and laboratory truth (0.81, 0.84).
- **Over-credit of Moderate PP3 was suggestive** (panels 0.02 [−1.69–0.89] against 2).
- **The cap**, acting only at Strong, applied to 0.5–2% of co-applications.

**Conclusion.** Location-type criteria overlap. PP2 should not be added beside PM1 without independent support. Supporting PP3 beside PM1 appears sound; Moderate PP3 beside PM1 warrants caution.

---

## Introduction

ACMG/AMP classification of missense variants combines three kinds of evidence (Richards et al. 2015):
- computational prediction (PP3);
- location in a mutational hotspot or critical domain (PM1);
- membership of a gene with little benign missense variation (PP2).

The Bayesian points framework adds these in log-likelihood-ratio units, which assumes they are independent (Tavtigian et al. 2018, 2020).

The assumption is questionable. Predictors such as AlphaMissense, REVEL and VEST4 (Cheng et al. 2023; Ioannidis et al. 2016; Carter et al. 2013) learn which genes and protein regions are intolerant of variation, and that is the information PM1 and PP2 encode.
- **ClinGen** capped PP3 + PM1 at Strong (Pejaver et al. 2022; Stenton et al. 2024).
- **UK guidance** advises against counting constraint evidence in both PM1 and PP2 (Durkie et al. 2024).
- **CanVIG-UK** permits PP3 and PM1 together (Garrett et al. 2021).

None of these rules rests on a measurement of the dependence. It is also unknown how often, and at what strength, curators combine the criteria.

Measuring the dependence needs evidence that the classifications themselves did not produce. ClinVar labels partly reflect PP3 and PM1, and several predictors were trained on them (Grimm et al. 2015). De novo variants avoid this. Their excess over mutation-rate expectation in affected children defines a likelihood ratio (LR) without variant-level labels.

**What we did:**
1. We measured the dependence between PP3, PM1 and PP2 on the de novo scale.
2. We checked it on independent truth: ClinVar, functional assays and re-scored expert-panel and laboratory classifications.
3. We measured the evidence PP2 and PP3 add once PM1 applies in expert-panel and laboratory practice, against the points curators give them.

---

## Methods

### Analysis log
Analyses followed a public, version-controlled analysis log (seventeen amendments, each committed before the analysis it governs; Supplementary Note 1). It was not deposited in an independent registry beforehand.

- **The primary pre-specified criterion** concerned how much predictor evidence survives conditioning on gene and position. It was not met under its original position control.
- **The PP3 × PM1 dependence test** was specified before its results.
- **Not pre-specified, and descriptive:** the PM1 × PP2 and triple tests, and all analyses of expert-panel and laboratory practice.

### Data
**De novo cohorts.**
- **Developmental disorders (DD):** 31,058 trios (Kaplanis et al. 2020).
- **Autism:** 12,316 probands, with 4,051 unaffected siblings as a negative control (Zhou et al. 2022).
- **Public:** 4,256 trios with five other phenotypes (Jin et al. 2017; denovo-db v1.6.1, Turner et al. 2017).

Sites are coding positions on GENCODE v32 canonical transcripts, excluding those within 3 bp of an exon boundary.

**Annotations.**
- **Predictor scores:** dbNSFP 4.1a and 5.4a (Liu et al. 2020), and AlphaMissense.
- **Constraint:** gene and regional missense constraint from gnomAD (Karczewski et al. 2020; Chao et al. 2024), and constrained coding regions (Havrilla et al. 2019).
- **Disease genes:** DDG2P monoallelic genes (Thormann et al. 2019).
- **Structure:** AlphaFold DB (Varadi et al. 2022).

**Classifications.**
- **Expert panels.** 6,177 missense classifications by 46 expert panels in the ClinGen Evidence Repository (export of 2026-10-07).
- **Laboratories.** 87,873 ClinVar 2026-10 germline submissions for missense SNVs from 472 submitters, excluding expert-panel and practice-guideline submissions. Each lists ACMG/AMP codes in free text.
- **Code parser.** A parser read codes and strengths from the free text and skipped negated codes. It was validated against the structured codes of 13,091 Evidence Repository classifications, using their free-text summaries:
  - overall precision 0.92 and recall 0.98;
  - for PP3, PM1 and PP2: precision 0.95, 0.98 and 0.98, recall ≥ 0.99, and strength agreement ≥ 0.99.

**Functional assays.** Human single-substitution assays in ProteinGym v1.3 (Notin et al. 2023), mapped to analysed transcripts: 82 assays in 70 genes. Damaging was defined as `DMS_score_bin` = 0.

### Dependence on the de novo scale
**Likelihood ratios.**
- **Expected counts** E come from Roulette mutation rates (Seplyarskiy et al. 2023), scaled by synonymous de novos; the excess is X = O − E.
- **For a set of sites b**, LR = (X_b/ΣX)/(E_b/ΣE).

**Flags.**
- **PP3:** the ClinGen thresholds for AlphaMissense (primary), REVEL and VEST4 (Pejaver et al. 2022; Bergquist et al. 2025).
- **PM1-like:** regional missense o/e < 0.4 or constrained-coding-region percentile ≥ 95.
- **PM1-strict:** o/e < 0.2, the level Chao et al. calibrated to Moderate.
- **PP2-like:** gene missense Z ≥ 3.09.

**Measures.**
- **Dependence:** R = LR(A ∧ B)/[LR(A) · LR(B)]. R = 1 is independence; R < 1 means summing A and B overstates their joint evidence. R3 is the three-criterion analogue.
- **PP3's own evidence alone and beside PM1, on the same scale:**
  - points alone, ln LR(PP3)/ln 2.406;
  - points beside PM1, ln LR(PP3 | PM1)/ln 2.406;
  - their ratio, the share of PP3's own evidence kept beside PM1.

  We do not compare de novo points with nominal ACMG weights, because the de novo scale is compressed relative to the clinical scale.
- **At Strong**, the fraction of summed log-evidence kept, f, can be compared with ClinGen's cap (f = 0.67 for Strong PP3 with Moderate PM1).

**Uncertainty.** Intervals resample genes (1,000 replicates).

### Independent truth and the evidence each criterion adds
**Re-scored truth.** Expert-panel and laboratory classifications were re-scored after removing:
- computational and location criteria: PP3, BP4, PM1, PP2, BP1, BP3, PM4, PM5;
- de novo criteria: PS2, PM6.

Truth was P/LP at ≥ 6 remaining points and B/LB at ≤ −1.

**Credit.** On this truth, the evidence criterion C adds once PM1 applies is ln[LR(C ∧ PM1)/LR(PM1)]/ln 2.406, on the ClinGen points scale. It is compared with the points curators applied. For PP3 it was computed by applied strength. We report credit rather than R because classification uses credit: R can be well below 1 while the credit equals the nominal points, when C alone is worth more than its nominal points.

**Sensitivity and out-of-sample analyses:**
- truth thresholds (P/LP ≥ 4, 6 or 8; B/LB ≤ −1 or −3);
- frequency criteria (BA1, BS1, BS2, PM2) also removed;
- PS1, PP5 and BP6 also removed;
- 200 random gene splits;
- for panels, a temporal split (approved before 2025 against 2025–2026);
- for panels, PP3 defined by predictor score instead of curator choice;
- for laboratories, clustering by submitter as well as by gene.

**ClinVar and assays.**
- **ClinVar:** P/LP versus B/LB, with DD de novo sites removed.
- **Assays:** R within each assay with a position-clustered bootstrap, pooled by random effects. An assay was included where 5–95% of its variants lie at PM1-like sites.

### Classification rules
Classifications were re-scored on the points scale (P ≥ 10, LP 6–9, VUS 0–5, LB −1 to −6, B ≤ −7) under four rules:
1. the ClinGen cap;
2. PP3 + PM1 + PP2 ≤ 4;
3. PP2 not added when PM1 applies;
4. PP3 limited to Supporting when PM1 applies.

---

## Results

### The criteria are not independent on the de novo scale
For AlphaMissense PP3 at Supporting with a PM1-like flag, R was 0.44 [0.41, 0.47] in DD. The dependence replicated in autism and in the public cohort, and persisted within DDG2P disease genes (Table 1, Figure 3).

**Table 1. Dependence on the de novo scale** (R; gene-clustered 95% intervals).

| Criteria | DD | Autism | Public |
|---|---|---|---|
| PP3 (AlphaMissense) × PM1-like | 0.44 [0.41, 0.47] | 0.56 [0.47, 0.66] | 0.67 [0.48, 0.96] |
| — within DDG2P genes, same-gene comparisons | 0.69 [0.65, 0.73] | 0.74 [0.67, 0.82] | 0.75 [0.62, 0.89] |
| PP3 (AlphaMissense) × PP2-like | 0.61 [0.57, 0.65] | 0.65 [0.55, 0.76] | 0.78 [0.57, 1.07] |
| PM1-like × PP2-like | 0.31 [0.28, 0.33] | 0.36 [0.30, 0.45] | — |
| PP3 × PM1-like × PP2-like (R3) | 0.13 [0.11, 0.14] | 0.18 [0.13, 0.24] | — |

**Robustness.**
- **Other predictors.** Dependence was weaker for REVEL (0.57) and VEST4 (0.67). Across 40 informative predictors in DD, every R × PM1-like had an upper bound below 1.
- **Within genes.** Within genes across all genes, autism and the public cohort were consistent with independence (0.96). Part of the dependence therefore lies between genes.
- **Other checks.** It held in ordered and disordered regions, at buried and exposed sites, and after removing the 50 genes with the largest excess (Supplementary Note 4).

**PM1 × PP2.** This was the most dependent pair. It is partly dependent by construction, because both of our flags are built from missense constraint.

**PP3 keeps about half its own evidence beside PM1** (Table 2):
- It does so at every strength, for both predictors and in both cohorts. The dependence is therefore not specific to one PP3 tier.
- At Strong with PM1-strict, f was 0.72 [0.70, 0.74] in DD. That is milder than the cap (0.67), so on this scale the cap removes enough at Strong.

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

**Why the criteria depend on each other.** Missense predictors carry largely location information.
- Across 39 informative predictors, a median 0.32 of top-2% log-evidence was incremental to gene and 50-residue window.
- A random partition of genes into the same number of units kept 0.51.
- Inside constrained regions, AlphaMissense flags 51% of mutation-expected sites, against 19% overall.

A perfect predictor would also lose evidence within windows, by an amount the data cannot pin down. This share therefore measures evidence incremental to location, not skill (Supplementary Note 3).

### The dependence on independent truth
**Table 3. Dependence on truth independent of the de novo data.**

| Truth | Criteria | Scope | R [95% CI] |
|---|---|---|---|
| ClinVar P/LP vs B/LB* | AlphaMissense PP3 × PM1-like | DDG2P genes | 0.32 [0.28, 0.37] |
| Expert-panel truth | Panel PP3 × panel PM1 | 2,942 classifications, 130 genes | 0.21 [0.12, 0.31] |
| Laboratory truth | Laboratory PP3 × PM1 | 39,015 submissions, 5,976 genes | 0.18 [0.14, 0.25] |
| Laboratory truth | Laboratory PM1 × PP2 | same | 0.27 [0.19, 0.44] |
| Laboratory truth | Laboratory PP3 × PP2 | same | 0.42 [0.35, 0.51] |
| Functional assays | AlphaMissense PP3 × PM1-like | 13 assays, 10 genes outside DDG2P | 0.94 [0.88, 0.99] |

\*Partly circular: ClinVar classifications use PP3 and PM1 themselves. This biases R towards 1.

- **Clinical truth.** On expert-panel and laboratory truth the dependence was strong and consistent. Clustering laboratory intervals by submitter widened them but left every conclusion unchanged; for example, PP3 × PM1 was 0.18 [0.08, 0.43].
- **Functional truth.** On assays it was weak and varied between genes. Within the assayed regions, the PM1-like flag barely predicted damage (median LR 1.10 outside DDG2P genes), so little overlap is possible there. BRCA1, TP53 and PTEN could not be tested, because nearly all their assayed residues lie at PM1-like sites.

### How panels and laboratories combine the criteria
**Table 4. Co-application in practice.**

| | Expert panels | Laboratories |
|---|---|---|
| Missense classifications | 6,177 | 87,873 |
| PP3 beside PM1 | 787 | 4,471 |
| — PP3 Supporting / Moderate / Strong | 605 / 178 / 4 | 4,207 / 184 / 80 |
| Share at Strong, where the cap acts | 0.5% | 1.8% |
| PP2 beside PM1 | 174 | 1,481 |
| All three | 138 | 1,035 |

Three panels whose specifications may permit co-application of PM1 and PP2 account for 145 of the 174 panel cases: RASopathy, FBN1 and Monogenic Diabetes. Some classifications predate current specifications (ClinGen RASopathy Expert Panels 2025).

### The evidence PP2 and PP3 add once PM1 applies
**Table 5. Points earned once PM1 applies, against points applied** (clinical truth; gene-clustered 95% intervals; B/LB = benign variants in the joint cell).

| Criterion | Expert panels | Laboratories | Applied |
|---|---|---|---|
| PP2 | not estimable (1 B/LB) | −0.22 [−0.67, 0.40] (565; 57 B/LB) | 1.00 |
| PP3, any strength | 0.60 [0.14, 1.05] (304; 17 B/LB) | 0.85 [0.54, 1.22] (1,768; 74 B/LB) | 1.24 / 1.09 |
| PP3 Supporting | 0.81 [0.24, 1.68] (235; 11 B/LB) | 0.84 [0.52, 1.26] (1,674; 71 B/LB) | 1.00 |
| PP3 Moderate | 0.02 [−1.69, 0.89] (67; 6 B/LB) | not estimable (3 B/LB) | 2.00 |

**PP2 beside PM1** showed no detectable benefit in laboratory practice. The finding held:
- clustering by submitter (−0.22 [−1.23, 0.79]);
- with PS1, PP5, BP6, PS2 and PM6 also removed from truth (−0.41 [−0.89, 0.24]).

Expert-panel truth could not test this pair.

**Supporting PP3 beside PM1** was close to its nominal point in both panels and laboratories. That fits Table 2: PP3 keeps about half its own evidence beside PM1, and on clinical truth PP3 alone is worth more than its nominal point (LR 8.2 on panel truth, against 2.4 for one point).

**Moderate PP3 beside PM1** earned at most 0.89 points (upper 95% bound) against 2 in expert panels. This is suggestive rather than established:
- it rests on 6 B/LB variants;
- Moderate could not be distinguished from Supporting (difference −0.79 [−2.81, 0.27]; p = 0.14; minimum detectable difference about 4 points);
- laboratory data had too few B/LB to test it.

**Panel-level robustness of the overall PP3 credit** (Figure 2):
- held-out genes: 0.59, below the points applied in 98% of 200 splits;
- classifications approved in 2025–2026: 0.40 [−0.25, 1.07] against 1.40;
- PS1, PP5, BP6, PS2 and PM6 removed from truth: 0.57 [0.08, 1.04];
- truth thresholds: 0.60–0.65;
- frequency criteria removed: 0.49 [−0.06, 1.31].

**Defining PP3 by score rather than by panel choice** gave mixed results (Supplementary Note 7):
- REVEL ≥ 0.644 earned 0.68 [0.34, 1.19];
- AlphaMissense ≥ 0.792 earned 1.66 [0.98, 2.58];
- by score band, only REVEL at Moderate (0.58 [0.08, 1.46]) and AlphaMissense at Strong (1.75 [0.76, 3.74]) had bounded intervals below their nominal points.

### Consequences for classification
**Table 6. Classifications that change class under alternative rules.**

| Rule | Expert panels (6,177) | | | | Laboratories (87,873) | | | |
|---|---|---|---|---|---|---|---|---|
| | Changed | P→LP | LP→VUS | Other | Changed | P→LP | LP→VUS | Other |
| ClinGen cap (PP3 + PM1 ≤ 4) | 24 | 12 | 12 | 0 | 40 | 30 | 10 | 0 |
| PP3 + PM1 + PP2 ≤ 4 | 35 | 17 | 18 | 0 | 63 | 39 | 22 | 2 |
| PP2 not added when PM1 applies | 39 | 28 | 11 | 0 | 352 | 124 | 223 | 5 |
| PP3 limited to Supporting when PM1 applies | 43 | 22 | 21 | 0 | 90 | 39 | 51 | 0 |

"Other" counts any other class transition.

---

## Discussion

**What we found.**
- **No independence.** Computational evidence is not independent of the location evidence that ACMG/AMP credits separately. The dependence appears on a label-free de novo scale in three cohorts and within disease genes, and on independent expert-panel and laboratory truth.
- **Why.** Predictors have learned where pathogenic variation falls, which is what PM1 and PP2 describe.
- **Where it is weak.** On functional assays, because there PM1 carries little information about molecular damage.

**PP2 beside PM1.** This is the clearest practical finding.
- **The guidance.** UK guidance advises against counting constraint evidence in both criteria (Durkie et al. 2024).
- **Practice.** Laboratories co-applied them in 1,481 submissions.
- **The evidence.** On laboratory truth PP2 beside PM1 showed no detectable benefit, robust to submitter clustering and to stricter truth.
- **Recommendation.** Do not add PP2 when PM1 applies, unless the two rest on independent evidence. That would change 352 of 87,873 laboratory classifications.
- **Panels.** Expert panels co-apply the pair mainly where their specifications appear to permit it. Our data cannot test those panel-specific definitions.

**PP3 beside PM1.** PP3 keeps about half its own evidence beside PM1. Whether that is less than the points curators give depends on how much PP3 is worth alone:
- **At Supporting,** where 77% of panel and 94% of laboratory co-application occurs, PP3 beside PM1 earned about its nominal point. Practice there appears sound.
- **At Moderate,** the panel data suggest over-credit, but the estimate is imprecise.
- **At Strong,** the cap removes as much as the dependence does, on the de novo scale.
- **Our suggestion.** Limiting PP3 to Supporting when PM1 applies changes few classifications (43 panel, 90 laboratory) and is consistent with the data. Confirming it needs more classifications with Moderate PP3 and benign truth.

**Relation to previous work.**
- **Calibration** of PP3 (Pejaver et al. 2022; Bergquist et al. 2025) and of regional constraint for PM1 (Chao et al. 2024) treats each criterion alone. Gene- and domain-aware calibration addresses heterogeneity, not dependence between criteria (Tejura et al. 2024; Badonyi & Marsh 2025; Chen et al. 2026).
- **Evaluation.** Gene identity dominates ClinVar evaluation of predictors (Grimm et al. 2015; Harrizi et al. 2026).
- **What is new.** We are not aware of a previous measurement of the dependence between PP3, PM1 and PP2, or of how they are combined in expert and laboratory practice.

**Limitations.**
1. **Pre-specification.** The analysis log was not independently registered, and the practice analyses are descriptive.
2. **Scale.** The de novo scale reflects high-penetrance variants in severe early-onset disorders; we compare evidence within that scale, not against nominal points.
3. **Truth.** Panel and laboratory truth come from the same classifications whose PP3 and PM1 we evaluate, after those criteria are removed.
   - Truth thresholds select on the remaining evidence.
   - Curators choose PP3 strength knowing that evidence.
   - Benign truth relies heavily on frequency criteria.

   The sensitivity analyses reduce, but do not remove, these concerns.
4. **Proxies.** Our PM1 and PP2 flags on the de novo scale are constraint-based proxies, and panel definitions of PM1 vary.
5. **Thin cells.** Several estimates rest on few benign variants, notably panel Moderate PP3 (6) and PP3 overall (17).
6. **Assays** tested ten genes outside NDD.
7. **Coverage.** Transcripts match MANE Select in 74.4% of genes, and the cohorts are predominantly of European ancestry.

**Conclusion.** PP3, PM1 and PP2 are not independent. In practice:
- **PP2 beside PM1** showed no detectable benefit and is commonly added by laboratories.
- **Supporting PP3 beside PM1** appears sound.
- **Moderate PP3 beside PM1** warrants caution.

Not adding PP2 when PM1 applies is the change best supported by these data.

---

## Figures and tables
- **Figure 1.** Design: de novo excess and conditioning (`fig1_design.png`).
- **Figure 2.** PP3 credit beside PM1, earned against applied, by data source and check (`fig2_credit.png`).
- **Figure 3.** Dependence on the de novo scale by predictor, flag and cohort (`fig4_overlap.png`).
- **Tables 1–6** in the text.
- **Supplementary Notes 1–10 and Supplementary Tables 1–7** in `SUPPLEMENT_v12.md`.

## Data and code availability
All data are public: published de novo tables, ClinVar, the ClinGen Evidence Repository, ProteinGym, dbNSFP (academic licence), AlphaMissense, AlphaFold DB, gnomAD, DDG2P and Roulette. Code, the analysis log and all outputs are in the study repository [URL; Zenodo DOI at submission].

## Competing interests
[To be completed.]

## Funding
[To be completed.]

---

## References
- Badonyi M, Marsh JA. acmgscaler: standardized gene-level variant effect score calibration within the ACMG/AMP framework. 2025.
- Bergquist T, Stenton SL, Nadeau EAW, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 2025.
- Carter H, Douville C, Stenson PD, Cooper DN, Karchin R. Identifying Mendelian disease genes with the variant effect scoring tool. *BMC Genomics* 2013;14(Suppl 3):S3.
- Chao KR, Wang L, Panchal R, et al. The landscape of regional missense mutational intolerance quantified from 125,748 exomes. bioRxiv 2024. doi:10.1101/2024.04.11.588920.
- Chen et al. Gene- and domain-aware calibration increases the clinical utility of variant effect predictors. Research Square preprint, 2026.
- Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 2023;381:eadg7492.
- ClinGen RASopathy Expert Panels. Updated ACMG/AMP specifications for variant interpretation and gene curations from the ClinGen RASopathy expert panels. *Genet Med Open* 2025.
- Durkie M, Cassidy E-J, Berry I, et al. ACGS best practice guidelines for variant classification in rare disease 2024. Association for Clinical Genomic Science, 2024.
- Garrett A, Durkie M, Callaway A, et al. Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. *J Med Genet* 2021;58:297–304.
- Grimm DG, Azencott CA, Aicheler F, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat* 2015;36:513–523.
- Harrizi et al. Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. bioRxiv 2026. doi:10.64898/2026.08.21.746180.
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
- Stenton SL, Pejaver V, Bergquist T, et al. Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med* 2024.
- Tavtigian SV, Greenblatt MS, Harrison SM, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med* 2018;20:1054–1060.
- Tavtigian SV, Harrison SM, Boucher KM, Biesecker LG. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Hum Mutat* 2020;41:1734–1737.
- Tejura M, Fayer S, McEwen AE, et al. Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes. *Am J Hum Genet* 2024.
- Thormann A, Halachev M, McLaren W, et al. Flexible and scalable diagnostic filtering of genomic variants using G2P with Ensembl VEP. *Nat Commun* 2019;10:2373.
- Turner TN, Yi Q, Krumm N, et al. denovo-db: a compendium of human de novo variants. *Nucleic Acids Res* 2017;45:D804–D811.
- Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Res* 2022;50:D439–D444.
- Zhou X, Feliciano P, Shu C, et al. Integrating de novo and inherited variants in 42,607 autism cases identifies mutations in new moderate-risk genes. *Nat Genet* 2022;54:1305–1319.
