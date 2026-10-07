# Computational and location evidence are not independent in missense variant classification, and current practice over-credits their combination

*Draft 10, 7 October 2026. Supplementary material: `SUPPLEMENT_v10.md`.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** ACMG/AMP classification sums computational (PP3) and location (PM1, PP2) evidence as if independent. ClinGen caps PP3 + PM1 only at Strong. We tested their independence without variant labels and audited how expert panels and laboratories combine them.

**Methods.** Likelihood ratios came from the excess of de novo missense variants over mutation-rate expectation in 47,630 trios and probands across three cohorts. Dependence was checked on ClinVar, on functional assays, and on 6,177 expert-panel and 87,873 laboratory classifications re-scored without computational, location or de novo criteria.

**Results.**
- PP3 and PM1 were not independent. Their joint evidence was 0.44 [0.41–0.47] of the product of their separate likelihood ratios in developmental disorders, 0.69 within disease genes, and 0.18–0.21 on panel and laboratory truth. On functional assays the dependence was weak (0.94).
- Panels and laboratories applied PP3 beside PM1 almost always below Strong (99% and 98%), where the cap does not act.
- There, PP3 earned 0.60 points [0.14–1.05] against 1.24 applied. This held in held-out genes, in later classifications and with PP3 defined by score rather than curator choice. Moderate PP3 earned less than its 2 points (upper bound 0.89).

**Conclusion.** Computational and location evidence overlap, and practice over-credits PP3 beside PM1 at strengths the cap does not reach. Limiting PP3 to Supporting when PM1 applies is consistent with these data.

---

## Introduction

Missense variants are the largest class of variants of uncertain significance. ACMG/AMP classification combines three kinds of evidence for them (Richards et al. 2015):
- a computational prediction (PP3);
- location in a hotspot or critical domain (PM1);
- membership of a gene with little benign missense variation (PP2).

The Bayesian points framework adds these in log-likelihood-ratio units, which assumes they are independent (Tavtigian et al. 2018, 2020).

They may not be. Modern predictors learn which genes and protein regions are intolerant of variation, and that is also what PM1 and PP2 encode.
- **ClinGen** recognised the overlap and capped PP3 + PM1 at Strong (Pejaver et al. 2022; Stenton et al. 2024).
- **UK guidance** advises against counting the same constraint evidence in both PM1 and PP2 (Durkie et al. 2024).
- **CanVIG-UK** permits PP3 and PM1 together at full weight (Garrett et al. 2021).

None of these rules rests on a measurement of how much the criteria overlap. Nor is it known how often curators actually stack them, or at what strength.

Measuring the overlap needs evidence that does not come from the classifications themselves. ClinVar labels partly reflect PP3 and PM1, and several predictors were trained on them (Grimm et al. 2015). De novo variants provide such evidence. Their excess over mutation-rate expectation in affected children defines a likelihood ratio (LR) without variant-level labels, for any set of sites and within any gene.

We used the de novo excess to measure the dependence between PP3, PM1 and PP2, checked it on four independent kinds of truth, and audited how expert panels and clinical laboratories combine these criteria.

---

## Methods

### Analysis log
The analysis plan and seventeen amendments were committed to a public version-controlled repository before the analyses they govern. None was deposited in an independent registry beforehand, so we refer to an analysis log rather than a pre-registration.

- **What was pre-specified.** The primary pre-specified criterion concerned how much predictor evidence survives conditioning on gene and position. It was not met under its original position control, for reasons given in Supplementary Note 1. The PP3 × PM1 dependence test was specified (amendment A10) before its results.
- **What was not.** The PM1 × PP2 and triple tests, the practice audit and all analyses of panel and laboratory classifications were added afterwards. They are descriptive.
- **Deviation log:** Supplementary Note 1.

### Data
**De novo cohorts.**
- **Developmental disorders (DD):** 31,058 trios (Kaplanis et al. 2020).
- **Autism:** 12,316 probands, with 4,051 unaffected siblings as a negative control (Zhou et al. 2022).
- **Public:** 4,256 trios with five other phenotypes (Jin et al. 2017; denovo-db v1.6.1).

Sites are coding positions on GENCODE v32 canonical transcripts, excluding those within 3 bp of an exon boundary.

**Classifications.**
- **ClinGen Evidence Repository** (export of 2026-10-07): 6,177 missense classifications by 46 expert panels.
- **ClinVar 2026-10:** germline submissions from laboratories (expert-panel and practice-guideline submissions excluded) whose comments list ACMG/AMP codes for missense SNVs. This gave 87,873 submissions from 472 submitters; codes were parsed from free text, with negated codes skipped.

**Functional assays.** All human single-substitution assays in ProteinGym v1.3 that could be mapped to analysed transcripts (82 assays, 70 genes), with damaging defined as `DMS_score_bin` = 0.

### Measuring dependence on the de novo scale
**Likelihood ratios.**
- **Expected counts** E come from Roulette mutation rates (Seplyarskiy et al. 2023), scaled to each cohort by its synonymous de novos. The excess is X = O − E.
- **For a set of sites b**, LR = (X_b/ΣX)/(E_b/ΣE).
- **Within-gene LRs** compare each de novo with the mutational-opportunity profile of its own gene.

**Flags.**
- **PP3:** the ClinGen Supporting or Strong threshold for AlphaMissense (primary), REVEL or VEST4 (Pejaver et al. 2022; Bergquist et al. 2025).
- **PM1-like:** regional missense o/e < 0.4 or constrained-coding-region percentile ≥ 95 (Chao et al. 2024; Havrilla et al. 2019).
- **PM1-strict:** o/e < 0.2, the level Chao et al. calibrated to Moderate.
- **PP2-like:** gene missense Z ≥ 3.09.

**Dependence.** For criteria A and B, R = LR(A ∧ B) / [LR(A) · LR(B)]. R = 1 is independence, as the points framework assumes; R < 1 means summing the two overstates their joint evidence. For three criteria, R3 is defined in the same way.

**Uncertainty.** Intervals come from resampling genes (1,000 replicates).

### Independent truth
Four sources of truth were used, none derived from the de novo data:

1. **ClinVar** P/LP versus B/LB, with DD de novo sites removed.
2. **Expert-panel truth.** Classifications were re-scored from their applied criteria after removing the computational and location criteria (PP3, BP4, PM1, PP2, BP1, BP3, PM4, PM5). P/LP is ≥ 6 remaining points and B/LB is ≤ −1. The panels' own PM1 and PP3 are then compared against evidence they did not contribute to. A sensitivity analysis also removed the de novo criteria (PS2, PM6).
3. **Laboratory truth**, constructed in the same way from laboratory submissions.
4. **Functional assays.** Within each assay, R was computed with a position-clustered bootstrap, in assays where 5–95% of variants fall at PM1-like sites. Results were pooled by random-effects meta-analysis, primarily in genes outside DDG2P.

### Practice audit
**Co-application.** We counted how often panels and laboratories apply PP3 beside PM1 (and PP2), and at what strength.

**Credit.** The credit PP3 earns once PM1 applies is ln[LR(PP3 ∧ PM1)/LR(PM1)]/ln 2.406, on the ClinGen points scale. It was computed by PP3 strength, with gene-clustered intervals, and Moderate and Supporting PP3 were compared directly (difference test and its minimum detectable difference).

**Checks:**
- truth thresholds (P/LP ≥ 4, 6 or 8; B/LB ≤ −1 or −3);
- frequency criteria (BA1, BS1, BS2, PM2) also removed from truth;
- 200 random gene splits;
- a temporal split: estimated on classifications approved before 2025, tested on those approved in 2025–2026;
- **selection into PP3:** among panel classifications mapped to analysed sites, PP3 was also defined by score (REVEL ≥ 0.644 or AlphaMissense ≥ 0.792), independent of curator choice.

**Rules.** Classifications were re-scored on the points scale (P ≥ 10, LP 6–9, VUS 0–5, LB −1 to −6, B ≤ −7) under:
- the ClinGen cap;
- a cap on PP3 + PM1 + PP2;
- PP2 not added when PM1 applies;
- PP3 limited to Supporting when PM1 applies.

---

## Results

### PP3, PM1 and PP2 are not independent on the de novo scale
For AlphaMissense PP3 at Supporting with a PM1-like flag, R was 0.44 [0.41, 0.47] in DD (Table 1, Figure 3).

**Table 1. Dependence of ACMG/AMP criteria on the de novo scale** (R; gene-clustered 95% intervals).

| Criteria | DD | Autism | Public |
|---|---|---|---|
| PP3 (AlphaMissense) × PM1-like | 0.44 [0.41, 0.47] | 0.56 [0.47, 0.66] | 0.67 [0.48, 0.96] |
| — within DDG2P genes, same-gene comparisons | 0.69 [0.65, 0.73] | 0.74 [0.67, 0.82] | 0.75 [0.62, 0.89] |
| PP3 (AlphaMissense) × PM1-strict | 0.37 [0.35, 0.40] | 0.50 [0.42, 0.61] | 0.55 [0.38, 0.83] |
| PP3 (AlphaMissense) × PP2-like | 0.61 [0.57, 0.65] | 0.65 [0.55, 0.76] | 0.78 [0.57, 1.07] |
| PM1-like × PP2-like | 0.31 [0.28, 0.33] | 0.36 [0.30, 0.45] | — |
| PP3 × PM1-like × PP2-like (R3) | 0.13 [0.11, 0.14] | 0.18 [0.13, 0.24] | — |

**Robustness.**
- **Other predictors.** Dependence was weaker for REVEL (0.57) and VEST4 (0.67). Across 40 informative predictors, every R × PM1-like had an upper bound below 1 in DD (median 0.51).
- **Within genes.** The within-gene dependence persisted in disease genes. Across all genes in autism and the public cohort it was indistinguishable from independence (0.96), so part of the dependence lies between genes.
- **Other checks.** It held in ordered and disordered regions, at buried and exposed sites, and after excluding the 50 genes with the largest excess (Supplementary Note 5).

**PM1 and PP2** were the most dependent pair. Part of this is by construction, because both of our flags are built from missense constraint, which UK guidance assigns to PP2 rather than PM1 (Durkie et al. 2024).

**Against the nominal ACMG weights.** When all three criteria apply:

| PP3 tier | Nominal sum | Under the cap | De novo evidence |
|---|---|---|---|
| Supporting PP3 | 4 points | 4 points | 3.65 [3.50, 3.79] |
| Strong PP3 | 7 points | 5 points | 4.05 [3.89, 4.20] |

The de novo scale is not the clinical scale, so this compares direction and size only.

**Why.** Missense predictors carry mostly location information.
- Across 39 informative predictors, a median 0.32 of top-2% log-evidence was incremental to the gene and 50-residue window.
- A random partition of genes into the same number of units kept 0.51, so the loss reflects position, not the act of conditioning.
- Inside constrained regions, AlphaMissense flags 51% of mutation-expected sites, against 19% overall. Among excess de novos the rise is only from 79% to 94%, so within those regions PP3 separates little.
- How far a perfect predictor would also lose evidence within windows depends on an enrichment parameter the data do not identify. We therefore treat this share as evidence incremental to location, not as skill (Supplementary Notes 8–9, including a gene- and position-adjusted predictor benchmark).

### The dependence on independent truth
**Table 2. Dependence of PP3 and PM1 across truth sets.**

| Truth | PP3 and PM1 definitions | Scope | R [95% CI] |
|---|---|---|---|
| De novo excess, DD | AlphaMissense; PM1-like | All genes | 0.44 [0.41, 0.47] |
| ClinVar P/LP vs B/LB | AlphaMissense; PM1-like | DDG2P genes | 0.32 [0.28, 0.37] |
| Expert-panel truth | Panel-applied PP3 and PM1 | 2,942 classifications, 130 genes | 0.21 [0.12, 0.31] |
| Laboratory truth | Laboratory-applied PP3 and PM1 | 39,075 submissions, 5,976 genes | 0.18 [0.14, 0.25] |
| Functional assays | AlphaMissense; PM1-like | 13 assays, 10 genes outside DDG2P | 0.94 [0.88, 0.99] |

**On classification truth the dependence was strong and consistent.** Panels and laboratories agree (0.21 and 0.18), on evidence they did not contribute to.

**On functional assays it was weak, and it varied.**
- Some assays showed strong dependence: EPHB2 0.40, PPARG 0.65, MET 0.83.
- Others showed none: SRC and MAPK1 about 1.00.
- One showed the reverse: HMGCR 1.40.

Including DDG2P genes, the pooled R was 0.89 [0.83, 0.95] across 16 assays. Assays for BRCA1, TP53 and PTEN could not test the question, because nearly all their residues fall at PM1-like sites. Assay truth measures molecular damage rather than disease, and the dependence measured here is largest where truth reflects disease.

### How panels and laboratories stack the criteria
**Table 3. Co-application of PP3, PM1 and PP2.**

| | Expert panels | Laboratories |
|---|---|---|
| Missense classifications | 6,177 | 87,873 |
| PP3 beside PM1 | 787 | 4,477 |
| — PP3 Supporting / Moderate / Strong | 605 / 178 / 4 | 4,212 / 185 / 80 |
| Share of PP3 + PM1 below Strong | 99% | 98% |
| PM1 beside PP2 | 174 (7 panels) | — |
| All three | 138 | — |

ClinGen's cap acts only when Strong PP3 meets Moderate PM1, so it constrains almost none of this stacking. PM1 + PP2 co-application was concentrated in seven panels: RASopathy, FBN1, Monogenic Diabetes, Brain Malformations, PTEN, Pulmonary Hypertension and ACADVL. Some of these classifications predate current specifications (ClinGen RASopathy Expert Panels 2025).

### What PP3 earns beside PM1
**Figure 2** shows, for each data source and check, the points PP3 earned once PM1 applied, against the points curators gave it.

**Expert panels.**
- **Overall.** PP3 earned 0.60 points [0.14, 1.05] against 1.24 applied. With the de novo criteria also removed from truth: 0.58 [0.08, 1.04] against 1.27.
- **By strength.**
  - **Supporting PP3** earned 0.81 [0.24, 1.68] against 1.00 (235 classifications; 11 B/LB).
  - **Moderate PP3** earned 0.02 [−1.69, 0.89] against 2.00 (67 classifications; 6 B/LB).
- **Moderate against Supporting.** The difference was −0.79 [−2.81, 0.27] (p = 0.14). With these numbers the minimum detectable difference is about 4 points. The data therefore show that Moderate PP3 earns less than it is given, but cannot show that it earns less than Supporting PP3.

**Checks.**
- **Truth thresholds:** the overall credit stayed at 0.60–0.65 across the four combinations tested.
- **Frequency criteria also removed:** 0.49 [−0.06, 1.31], with only 6 B/LB in the joint cell.
- **Held-out genes:** across 200 random splits, the credit was 0.61 on the estimating half and 0.59 on the held-out half. The held-out credit fell below the points applied in 98% of splits.
- **Later classifications:** for those approved in 2025–2026, PP3 earned 0.40 [−0.25, 1.07] against 1.40 applied.
- **Curator choice.** Panels applied PP3 preferentially. Within panel-PM1 classifications, variants where PP3 was applied had LR 2.16 for pathogenicity, against 1.05 for score-eligible variants where it was not. Defined by score alone, PP3 beside panel PM1 earned 0.71 [0.35, 1.20] against the 1 point a Supporting threshold implies.

**Laboratories.**
- **Overall.** PP3 beside PM1 earned 0.85 [0.53, 1.20] against 1.09 applied.
- **Supporting PP3:** 0.83 [0.53, 1.24].
- **Moderate PP3** earned 0.75 against 2.00, but with 3 B/LB variants its interval was unbounded.

**Reading.** Supporting PP3 beside PM1 is close to its nominal weight in both panels and laboratories. The excess credit comes from PP3 applied above Supporting, which the cap does not reach.

### Consequences for classification
**Table 4. Expert-panel missense classifications that change class under alternative rules** (6,177 classifications).

| Rule | Changed | P → LP | LP → VUS |
|---|---|---|---|
| ClinGen cap (PP3 + PM1 ≤ 4) | 24 | 12 | 12 |
| PP3 + PM1 + PP2 ≤ 4 | 35 | 17 | 18 |
| PP2 not added when PM1 applies | 39 | 28 | 11 |
| PP3 limited to Supporting when PM1 applies | 43 | 22 | 21 |

**Later classifications.** Among the 2,507 classifications approved in 2025–2026, limiting PP3 to Supporting changed 15 and the cap 16.

**Developmental disorders.** In DD, uncapped co-use with de novo evidence placed 2,438 variants at Pathogenic, in a class the de novo excess implies is only 84% disease-associated.

**ClinVar.** In DDG2P genes, 37.3% of AlphaMissense PP3-positive VUS lie at PM1-like sites. Without the cap, 7,101 would reach Likely Pathogenic from PP3 + PM1 alone.

---

## Discussion

**What we found.** Computational evidence is not independent of the location evidence that ACMG/AMP credits separately. The dependence appears:
- on a label-free de novo scale, in three cohorts;
- within disease genes;
- on independent expert-panel and laboratory truth.

It is strongest between PM1 and PP2, and between PP3 and PM1. The reason is that predictors have learned where pathogenic variation falls, and that is what PM1 and PP2 describe.

**Practice.**
- **The existing safeguard rarely acts.** It caps PP3 + PM1 only at Strong, while panels and laboratories apply PP3 beside PM1 almost always at Supporting or Moderate.
- **At Supporting, the double counting is small.** PP3 earned close to its nominal point in both panels and laboratories.
- **Above Supporting it is not.** On expert-panel truth, Moderate PP3 beside PM1 earned at most 0.89 points (upper 95% bound) against the two points given, and at Strong the nominal sum exceeds the de novo evidence by several points.
- **A rule that fits.** Limiting PP3 to Supporting when PM1 applies removes the credit the data cannot support. It needs no estimated parameter, and it would change 43 of 6,177 expert classifications.
- **Its limits.** The evidence against Moderate PP3 rests on 67 panel classifications, and we could not show that Moderate earns less than Supporting. We therefore present the rule as consistent with the data, for panels to check against their own classifications, rather than as an established standard.

**PM1 and PP2.** The data quantify what UK guidance already advises against: their dependence is the strongest of the three pairs. Seven ClinGen panels co-apply them. Our flags share constraint by construction, and expert truth had too few benign variants to test panel-defined PM1 + PP2. Panels that co-apply them should check that each rests on independent evidence, as ACGS recommends.

**Functional truth.** Dependence was weak on assays. That fits the explanation above: PM1 and predictors share information about where disease-causing variation falls, more than about molecular damage as such. It also cautions against generalising the size of the dependence beyond clinical classification.

**Relation to previous work.**
- **Calibration.** Predictor calibration for PP3 (Pejaver et al. 2022; Bergquist et al. 2025) and region calibration for PM1 (Chao et al. 2024) treat each criterion alone. Gene- and domain-aware calibration (Tejura et al. 2024; Badonyi & Marsh 2025; Chen et al. 2026) recognises heterogeneity but not dependence between criteria.
- **Evaluation.** Gene identity dominates ClinVar evaluation of predictors (Grimm et al. 2015; Harrizi et al. 2026).
- **What is new.** To our knowledge, this is the first measurement of the dependence between PP3, PM1 and PP2, and the first audit of how they are combined in expert and laboratory practice.

**Limitations.**
1. **The analysis log was not independently registered**, and the practice analyses were added after the main analyses.
2. **Estimand.** The de novo scale reflects high-penetrance variants in severe early-onset disorders; the clinical comparisons bridge this only partly.
3. **Truth.** Panel and laboratory truth are built from the same classifications whose PP3 and PM1 we evaluate, after removing those criteria. Truth thresholds select on other evidence, and laboratory codes were parsed from free text. The checks above reduce, but do not remove, these concerns.
4. **Proxies.** Our PM1 and PP2 flags are constraint-based proxies, and the panels' own PM1 definitions vary.
5. **Assays** tested dependence in ten genes outside NDD; key clinical genes such as BRCA1 and TP53 could not be tested.
6. **Coverage.** Transcripts match MANE Select in 74.4% of genes, and the cohorts are predominantly of European ancestry.

**Conclusion.** PP3, PM1 and PP2 are not independent, and the dependence carries into practice. Curators stack PP3 beside PM1 at strengths below the reach of the current cap, and PP3 above Supporting in that position adds little. Counting PP3 at most as Supporting when PM1 applies is a simple rule consistent with these data. A cap on location-type evidence as a whole deserves consideration.

---

## Figures and tables
- **Figure 1.** Design: de novo excess and conditioning (`fig1_design.png`).
- **Figure 2.** PP3 credit beside PM1, earned versus applied, by data source and check (`fig2_credit.png`).
- **Figure 3.** Dependence on the de novo scale, by predictor, flag and cohort (`fig4_overlap.png`).
- **Tables:**
  - Table 1: dependence on the de novo scale;
  - Table 2: dependence across truth sets;
  - Table 3: co-application in practice;
  - Table 4: classification changes.

**Supplementary:**
- **Notes 1–14:**
  - Note 1: analysis log and deviations;
  - Note 2: robustness of retention;
  - Note 3: simulation;
  - Note 4: de novo-trained classifiers (negative);
  - Note 5: full dependence results;
  - Note 6: pre-specified hierarchical model (negative);
  - Note 7: within-window model (negative);
  - Note 8: oracle;
  - Note 9: predictor benchmark and validation;
  - Note 10: PM1 × PP2 and the triple stack;
  - Note 11: expert-panel audit;
  - Note 12: laboratory replication;
  - Note 13: functional assays;
  - Note 14: selection into PP3 and the Moderate–Supporting test.
- **Tables 1–7.**

## Data and code availability
All data are public: published de novo tables, ClinVar, the ClinGen Evidence Repository, ProteinGym, dbNSFP (academic licence), AlphaMissense, AlphaFold DB, gnomAD constraint and Roulette. Code, the analysis log and all outputs are in the study repository [URL; Zenodo DOI at submission].

## Competing interests
[To be completed.]

## Funding
[To be completed.]

---

## References
- Badonyi M, Marsh JA. acmgscaler: standardized gene-level variant effect score calibration within the ACMG/AMP framework. 2025.
- Bergquist T, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 2025.
- Chao KR, et al. The landscape of regional missense mutational intolerance quantified from 125,748 exomes. bioRxiv 2024. doi:10.1101/2024.04.11.588920.
- Chen et al. Gene- and domain-aware calibration increases the clinical utility of variant effect predictors. Research Square preprint, 2026.
- ClinGen RASopathy Expert Panels. Updated ACMG/AMP specifications for variant interpretation and gene curations from the ClinGen RASopathy expert panels. *Genet Med Open* 2025.
- Durkie M, Cassidy E-J, Berry I, et al. ACGS best practice guidelines for variant classification in rare disease 2024. Association for Clinical Genomic Science, 2024.
- Garrett A, et al. Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. *J Med Genet* 2021;58:297.
- Grimm DG, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat* 2015.
- Harrizi et al. Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. bioRxiv 2026. doi:10.64898/2026.08.21.746180.
- Havrilla JM, et al. A map of constrained coding regions in the human genome. *Nat Genet* 2019.
- Jin SC, et al. Contribution of rare inherited and de novo variants in 2,871 congenital heart disease probands. *Nat Genet* 2017;49:1593–1601.
- Kaplanis J, et al. Evidence for 28 genetic disorders discovered by combining healthcare and research data. *Nature* 2020.
- Notin P, et al. ProteinGym: large-scale benchmarks for protein fitness prediction and design. *NeurIPS* 2023.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet* 2022.
- Richards S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med* 2015.
- Seplyarskiy V, et al. A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription (Roulette). *Nat Genet* 2023.
- Stenton SL, et al. Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med* 2024.
- Tavtigian SV, et al. *Genet Med* 2018; Tavtigian SV, et al. *Hum Mutat* 2020.
- Tejura M, et al. Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes. *Am J Hum Genet* 2024.
- Turner TN, et al. denovo-db. *Nucleic Acids Res* 2017.
- Zhou X, et al. Integrating de novo and inherited variants in 42,607 autism cases identifies mutations in new moderate-risk genes. *Nat Genet* 2022.
