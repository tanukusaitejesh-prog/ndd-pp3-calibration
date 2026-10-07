# Location-type evidence is double counted in missense variant classification: PP2 adds nothing beside PM1, and PP3 above Supporting adds little

*Draft 11, 7 October 2026. Supplementary material: `SUPPLEMENT_v11.md`.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** ACMG/AMP classification sums computational evidence (PP3) with location evidence: hotspot or domain (PM1) and gene-level (PP2). It treats them as independent, and ClinGen caps PP3 + PM1 only at Strong. We measured how much each criterion adds once another applies.

**Methods.** Likelihood ratios came from the excess of de novo missense variants over mutation-rate expectation in 47,630 trios and probands. We estimated the points each criterion earns beside PM1 on 6,177 expert-panel and 87,873 laboratory missense classifications, re-scored without computational or location criteria.

**Results.**
- **No pair of criteria was independent** on the de novo scale.
- **PP2 adds nothing beside PM1 in laboratory practice:** it earned −0.20 points [−0.66–0.42] against 1 applied. Laboratories stacked the two in 1,511 submissions, and not adding PP2 would change 357 classifications.
- **Supporting PP3 beside PM1 was close to its nominal point** on panel and laboratory truth (0.81, 0.83).
- **Above Supporting, PP3 beside PM1 earned well below its nominal points.** On the de novo scale, Moderate AlphaMissense earned −0.12 [−0.38–0.10] against 2 and Strong 1.10 [1.04–1.16] against 4. In expert panels, Moderate PP3 earned 0.02 [−1.69–0.89].
- **The cap acts only at Strong,** which accounts for 0.5–2% of PP3 + PM1 co-application.

**Conclusion.** Location-type criteria are double counted. Adding PP2 beside PM1 is not supported, and PP3 above Supporting adds little beside PM1, where the existing cap rarely acts.

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

We used the de novo excess to measure the dependence between PP3, PM1 and PP2, and checked it on four independent kinds of truth. We then measured how many points PP2 and PP3 earn once PM1 applies, in expert-panel and laboratory practice, and compared these with the points curators give them.

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

**Credit.**
- **Definition.** The points a criterion C earns once PM1 applies are ln[LR(C ∧ PM1)/LR(PM1)]/ln 2.406, on the ClinGen points scale. We report this credit rather than R because it is what classification uses. R can be well below 1 while the credit equals the nominal points, when C alone is worth more than its nominal points on that truth.
- **Where it was computed:**
  - PP2 and PP3 on panel and laboratory truth, with PP3 by applied strength and, for panels, by predictor score band;
  - PP3 on the de novo scale, by AlphaMissense and REVEL score band (Supporting, Moderate, Strong), with gene-clustered intervals.
- **Moderate against Supporting.** Panel-applied Moderate and Supporting PP3 were compared directly, with a difference test and its minimum detectable difference.

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
| ClinVar P/LP vs B/LB* | AlphaMissense; PM1-like | DDG2P genes | 0.32 [0.28, 0.37] |
| Expert-panel truth | Panel-applied PP3 and PM1 | 2,942 classifications, 130 genes | 0.21 [0.12, 0.31] |
| Laboratory truth | Laboratory-applied PP3 and PM1 | 39,075 submissions, 5,976 genes | 0.18 [0.14, 0.25] |
| Functional assays | AlphaMissense; PM1-like | 13 assays, 10 genes outside DDG2P | 0.94 [0.88, 0.99] |
| Laboratory truth: PM1 × PP2 | Laboratory-applied | 39,075 submissions | 0.28 [0.19, 0.44] |
| Laboratory truth: PP3 × PP2 | Laboratory-applied | 39,075 submissions | 0.42 [0.35, 0.51] |

\*ClinVar classifications partly use PP3 and PM1 themselves, so this row is partly circular. The bias runs towards R = 1, because curators who apply both enrich P/LP variants in the joint cell.

**On classification truth the dependence was strong and consistent.** Panels and laboratories agree (0.21 and 0.18), on evidence they did not contribute to.

**On functional assays it was weak, and it varied.**
- Some assays showed strong dependence: EPHB2 0.40, PPARG 0.65, MET 0.83.
- Others showed none: SRC and MAPK1 about 1.00.
- One showed the reverse: HMGCR 1.40.

Including DDG2P genes, the pooled R was 0.89 [0.83, 0.95] across 16 assays. Assays for BRCA1, TP53 and PTEN could not test the question, because nearly all their residues fall at PM1-like sites. Within the assayed regions, the PM1-like flag barely predicted molecular damage: median LR 1.10 outside DDG2P genes (IQR 1.04–1.42), against about 5–8 on clinical and de novo truth. Two criteria cannot overlap on information that one of them does not carry, so R near 1 on assays is expected. The dependence matters where PM1 is informative, which is in clinical classification.

### How panels and laboratories stack the criteria
**Table 3. Co-application of PP3, PM1 and PP2.**

| | Expert panels | Laboratories |
|---|---|---|
| Missense classifications | 6,177 | 87,873 |
| PP3 beside PM1 | 787 | 4,477 |
| — PP3 Supporting / Moderate / Strong | 605 / 178 / 4 | 4,212 / 185 / 80 |
| Share of PP3 + PM1 at Strong (where the cap acts) | 0.5% | 1.8% |
| PP2 beside PM1 | 174 (7 panels) | 1,511 |
| All three | 138 | 1,059 |

### PP2 adds nothing beside PM1
- **Laboratory truth.** PM1 and PP2 were strongly dependent: R 0.28 [0.19, 0.44]. Once PM1 applied, PP2 earned −0.20 points [−0.66, 0.42] against the 1 point applied (571 submissions; 57 B/LB).
- **De novo scale.** The constraint-based flags were similarly dependent (R 0.31 [0.28, 0.33]); part of this is by construction.
- **Expert-panel truth** had only one B/LB variant in the joint cell, so it could not test the pair.

UK guidance already advises against counting the same constraint evidence in both criteria (Durkie et al. 2024). These data show that in laboratory practice PP2 beside PM1 carries no detectable independent evidence.

### PP3 beside PM1, by strength
**Table 4. Points PP3 earns once PM1 applies** (95% intervals, gene-clustered; nominal points: Supporting 1, Moderate 2, Strong 4).

| Source of truth and definition of PP3 | Supporting | Moderate | Strong |
|---|---|---|---|
| **De novo excess, DD** — AlphaMissense band | −0.92 [−1.37, −0.59] | −0.12 [−0.38, 0.10] | 1.10 [1.04, 1.16] |
| **De novo excess, DD** — REVEL band | 0.47 [0.27, 0.64] | 1.04 [0.91, 1.16] | 1.82 [1.61, 2.01] |
| **De novo excess, autism** — AlphaMissense band | −0.04 [−0.72, 0.45] | −0.13 [−0.94, 0.35] | 0.90 [0.75, 1.04] |
| **De novo excess, autism** — REVEL band | 0.39 [−0.09, 0.71] | 0.86 [0.58, 1.08] | 1.39 [0.86, 1.76] |
| **Expert panels** — strength applied by the panel | 0.81 [0.24, 1.68] (n 235) | 0.02 [−1.69, 0.89] (n 67) | n 4 |
| **Expert panels** — REVEL band | −0.52 [−1.46, 1.17] | 0.58 [0.08, 1.46] | 1.84 [0.81, —] |
| **Expert panels** — AlphaMissense band | 1.63 [0.09, —] | 1.44 [0.34, —] | 1.75 [0.76, 3.74] |
| **Laboratories** — strength applied | 0.83 [0.53, 1.24] (n 1,681) | 0.75 [—] (3 B/LB) | 0 B/LB |

"—" marks an unbounded interval (≤ 3 B/LB in the joint cell). De novo bands are relative to all other PM1-like sites, and the de novo scale is not the clinical scale.

**Two independent lines of evidence.** The de novo excess and the expert panels are independent of each other; the panel rows share classifications.
- **On the de novo scale, every estimate above Supporting is below its nominal points, with intervals excluding them,** for both predictors and in both cohorts. Moderate AlphaMissense beside PM1 earns nothing (DD −0.12 [−0.38, 0.10]). Strong earns about 1 point (1.10 [1.04, 1.16]) where 4 are given.
- **In expert panels, Moderate PP3 applied beside PM1 earned at most 0.89 points** against 2. Where the score-band estimates are bounded, they also fall below their nominal points: REVEL Moderate below 2 and AlphaMissense Strong below 4.
- **Not every definition agrees.** AlphaMissense bands on panel truth were imprecise and their point estimates were at or above nominal at Supporting and Moderate. Defining PP3 by AlphaMissense eligibility alone, it earned 1.66 points [0.98, 2.58] beside panel PM1. REVEL eligibility gave 0.68 [0.34, 1.19], and either predictor 0.71 [0.35, 1.20] (Supplementary Note 14).
- **Moderate against Supporting.** In panel-applied data, Moderate could not be distinguished from Supporting (difference −0.79 [−2.81, 0.27]; p = 0.14; minimum detectable difference about 4 points).

**Supporting PP3.**
- **On clinical truth it is calibrated:** 0.81 in panels and 0.83 in laboratories, both intervals including 1.
- **On the de novo scale it earned less than 1** (REVEL 0.47; AlphaMissense below 0). Because the scales differ, we do not conclude that Supporting PP3 is over-credited in practice.
- **Laboratories overall** did not over-credit PP3 beside PM1 (0.85 [0.53, 1.20] against 1.09 applied), because 94% of their PP3 + PM1 co-applications are at Supporting.

**Why R < 1 does not imply over-credit at Supporting.** On panel truth PP3 alone has LR 8.2, more than its nominal Supporting LR of 2.4. Beside PM1 that falls to about the nominal value. The dependence removes the excess, not the nominal point.

**Expert-panel robustness of the overall credit** (0.60 [0.14, 1.05] against 1.24 applied; 17 B/LB in the joint cell; Figure 2):
- held-out genes: 0.59 on held-out halves, below the points applied in 98% of 200 splits;
- classifications approved in 2025–2026: 0.40 [−0.25, 1.07] against 1.40;
- de novo criteria removed from truth: 0.58 [0.08, 1.04];
- truth thresholds: 0.60–0.65;
- frequency criteria removed: 0.49 [−0.06, 1.31] (Supplementary Notes 11 and 14).

### Consequences for classification
**Table 5. Classifications that change class under alternative rules.**

| Rule | Expert panels (6,177) | Laboratories (87,873) |
|---|---|---|
| ClinGen cap (PP3 + PM1 ≤ 4) | 24 (12 P→LP, 12 LP→VUS) | 40 (30 P→LP, 10 LP→VUS) |
| PP3 + PM1 + PP2 ≤ 4 | 35 (17, 18) | 64 (40, 22) |
| **PP2 not added when PM1 applies** | 39 (28, 11) | **357 (125, 227)** |
| **PP3 limited to Supporting when PM1 applies** | 43 (22, 21) | 91 (40, 51) |

- **Later panel classifications.** Among the 2,507 approved in 2025–2026, limiting PP3 to Supporting changed 15 and the cap 16.
- **Developmental disorders.** In DD, uncapped co-use with de novo evidence placed 2,438 variants at Pathogenic, in a class the de novo excess implies is only 84% disease-associated.

---

## Discussion

**What we found.** Location-type evidence is double counted in missense classification.
- **The dependence.** On a label-free de novo scale, computational evidence (PP3) is not independent of hotspot or domain evidence (PM1) or gene-level evidence (PP2), and PM1 is not independent of PP2.
- **Why.** Predictors have learned where pathogenic variation falls, which is what PM1 and PP2 describe.
- **Where it holds.** On expert-panel and laboratory truth the dependence holds. On functional assays it is weak, because within assayed regions PM1 carries little information about molecular damage.

**PP2 beside PM1.** This is the clearest practical finding.
- **The guidance.** UK guidance advises against counting constraint evidence in both criteria.
- **Practice.** Laboratories stack them in 1,511 submissions and seven ClinGen panels in 174 classifications.
- **The evidence.** On laboratory truth, PP2 beside PM1 earned no detectable evidence.
- **Recommendation.** Do not add PP2 when PM1 applies, unless PP2 rests on evidence independent of PM1. This would change 357 laboratory and 39 expert classifications.

**PP3 above Supporting beside PM1.**
- **The evidence.** On the de novo scale, Moderate and Strong PP3 beside PM1 earned well under their nominal points, with intervals excluding them, for both predictors and in two cohorts. Panel-applied Moderate PP3 points the same way.
- **The cap misses it.** ClinGen's cap acts only when Strong PP3 meets Moderate PM1, which accounts for 0.5% of panel and 1.8% of laboratory PP3 + PM1 co-application.
- **Recommendation.** Counting PP3 at most as Supporting when PM1 applies is consistent with the data and needs no estimated parameter. It would change 43 expert and 91 laboratory classifications.
- **Its limits.** In panel-applied data, Moderate could not be distinguished from Supporting, and AlphaMissense-defined PP3 on panel truth did not show the deficit. We therefore present it as a rule for panels to check against their own classifications, not as an established standard.

**Supporting PP3 beside PM1** was calibrated on clinical truth. It earned less on the de novo scale, which measures association with severe developmental disorders rather than pathogenicity in general. Practice at Supporting therefore looks sound. Supporting accounts for 77% of panel and 94% of laboratory PP3 + PM1 co-application.

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
6. **Thin cells.** Several estimates rest on few benign variants: 17 for the overall panel PP3 credit and 6 for panel Moderate PP3. The de novo estimates above Supporting carry the claim.
7. **Coverage.** Transcripts match MANE Select in 74.4% of genes, and the cohorts are predominantly of European ancestry.

**Conclusion.** PP3, PM1 and PP2 are not independent, and the dependence matters in practice in two places:
- **PP2 beside PM1** adds no detectable evidence, yet laboratories and some expert panels add it routinely.
- **PP3 above Supporting beside PM1** adds little, below the reach of the current cap.

Two simple rules follow, both consistent with these data: do not add PP2 when PM1 applies, and count PP3 at most as Supporting when PM1 applies.

---

## Figures and tables
- **Figure 1.** Design: de novo excess and conditioning (`fig1_design.png`).
- **Figure 2.** PP3 credit beside PM1, earned versus applied, by data source and check (`fig2_credit.png`).
- **Figure 3.** Dependence on the de novo scale, by predictor, flag and cohort (`fig4_overlap.png`).
- **Tables:**
  - Table 1: dependence on the de novo scale;
  - Table 2: dependence across truth sets;
  - Table 3: co-application in practice;
  - Table 4: PP3 credit beside PM1 by strength, de novo and clinical;
  - Table 5: classification changes in panels and laboratories.

**Supplementary:**
- **Notes 1–15:**
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
  - Note 12: laboratory practice (replication, PM1 × PP2, re-scoring);
  - Note 13: functional assays;
  - Note 14: selection into PP3 and the Moderate–Supporting test;
  - Note 15: PP3 credit beside PM1 on the de novo scale, by score band.
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
- Bergquist T, Stenton SL, Nadeau EAW, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 2025.
- Chao KR, Wang L, Panchal R, et al. The landscape of regional missense mutational intolerance quantified from 125,748 exomes. bioRxiv 2024. doi:10.1101/2024.04.11.588920.
- Chen et al. Gene- and domain-aware calibration increases the clinical utility of variant effect predictors. Research Square preprint, 2026.
- ClinGen RASopathy Expert Panels. Updated ACMG/AMP specifications for variant interpretation and gene curations from the ClinGen RASopathy expert panels. *Genet Med Open* 2025.
- Durkie M, Cassidy E-J, Berry I, et al. ACGS best practice guidelines for variant classification in rare disease 2024. Association for Clinical Genomic Science, 2024.
- Garrett A, Durkie M, Callaway A, et al. Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. *J Med Genet* 2021;58:297–304.
- Grimm DG, Azencott CA, Aicheler F, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat* 2015;36:513–523.
- Harrizi et al. Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. bioRxiv 2026. doi:10.64898/2026.08.21.746180.
- Havrilla JM, Pedersen BS, Layer RM, Quinlan AR. A map of constrained coding regions in the human genome. *Nat Genet* 2019;51:88–95.
- Jin SC, Homsy J, Zaidi S, et al. Contribution of rare inherited and de novo variants in 2,871 congenital heart disease probands. *Nat Genet* 2017;49:1593–1601.
- Kaplanis J, Samocha KE, Wiel L, et al. Evidence for 28 genetic disorders discovered by combining healthcare and research data. *Nature* 2020;586:757–762.
- Notin P, Kollasch AW, Ritter D, et al. ProteinGym: large-scale benchmarks for protein fitness prediction and design. *Adv Neural Inf Process Syst* 2023.
- Pejaver V, Byrne AB, Feng BJ, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet* 2022;109:2163–2177.
- Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med* 2015;17:405–424.
- Seplyarskiy V, Koch EM, Lee DJ, et al. A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription. *Nat Genet* 2023.
- Stenton SL, Pejaver V, Bergquist T, et al. Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med* 2024.
- Tavtigian SV, Greenblatt MS, Harrison SM, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med* 2018;20:1054–1060.
- Tavtigian SV, Harrison SM, Boucher KM, Biesecker LG. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Hum Mutat* 2020;41:1734–1737.
- Tejura M, Fayer S, McEwen AE, et al. Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes. *Am J Hum Genet* 2024.
- Turner TN, Yi Q, Krumm N, et al. denovo-db: a compendium of human de novo variants. *Nucleic Acids Res* 2017;45:D804–D811.
- Zhou X, Feliciano P, Shu C, et al. Integrating de novo and inherited variants in 42,607 autism cases identifies mutations in new moderate-risk genes. *Nat Genet* 2022;54:1305–1319.
