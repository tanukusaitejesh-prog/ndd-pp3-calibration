# Computational (PP3) and hotspot/domain (PM1) evidence overlap for de novo missense variants in neurodevelopmental disorders: a label-free test of additivity

*Draft 4, 7 October 2026. Supersedes draft 3 (`MANUSCRIPT_v2.md`). Supplementary material: `SUPPLEMENT_v4.md`. Every value comes from a committed output.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** ACMG/AMP classification adds computational evidence (PP3) to hotspot/domain evidence (PM1) as if the two were independent. ClinGen caps their combination at Strong; other guidance allows uncapped co-use. We tested additivity without variant-level clinical labels.

**Methods.** Likelihood ratios (LRs) came from the excess of de novo missense variants over a mutation-rate expectation in 31,058 developmental-disorder (DD) trios, 12,316 autism probands and 4,256 public trios. Additivity was measured as R, the joint LR of PP3 and a PM1-like flag divided by the product of their separate LRs (R = 1 is additive). We checked the result on three truth sets independent of the de novo data: ClinVar, ClinGen expert-panel classifications with computational and position criteria removed, and multiplexed assays. We then pre-registered a label-free gene → region → PP3 model and tested it against the cap on ClinVar variants reclassified between 2023 and 2026.

**Results.** For AlphaMissense, R was 0.44 [95% CI 0.41–0.48] in DD, 0.56 [0.44–0.67] in autism and 0.67 [0.42–0.95] in the public cohort. Within DDG2P disease genes, comparing variants inside the same gene, R was 0.69 [0.68–0.70] and 0.74 [0.68–0.80]. Across all genes in autism, the within-gene R (0.96 [0.84–1.08]) was consistent with additivity. Overlap was larger on the clinical scale: ClinVar R 0.23 at Strong in DDG2P genes, expert-panel R 0.17 [0.08–0.61]. On expert-panel truth, Strong PP3 earned 1.83 points once PM1 applied. The cap allows 2 and uncapped co-use 4. The pre-registered hierarchical model calibrated worse than the cap on reclassified ClinVar variants (log-loss difference +0.147 [+0.106, +0.190]).

**Conclusion.** PP3 and PM1 evidence overlap substantially for de novo missense variants, including within known disease genes. ClinGen's cap is the best-calibrated of the rules tested. A label-free model built from de novo data measures the overlap but does not replace the cap.

---

## Introduction

Missense variants are the largest class of variants of uncertain significance. ACMG/AMP classification combines three kinds of evidence for them: computational predictions (PP3), location in a mutational hotspot or critical domain (PM1), and membership of a gene with little benign missense variation (PP2) (Richards et al. 2015).

In the Bayesian points framework these criteria add in log-likelihood-ratio units, with one prior for every gene (Tavtigian et al. 2018, 2020). ClinGen calibrated predictor thresholds against ClinVar (Pejaver et al. 2022; Bergquist et al. 2025). Recognising that PP3 and PM1 can describe the same property, Pejaver et al. capped their sum at Strong, a rule Stenton et al. (2024) restated as no double counting. CanVIG-UK instead permits both at full weight, regarding regional enrichment and in silico prediction as largely different evidence (Garrett et al. 2021).

How much the two overlap has not been measured. A predictor that has learned which genes and domains are intolerant separates pathogenic from benign variants pooled across genes, yet may add little once the gene and domain are known. ClinVar is a poor place to measure this:

- its labels are concentrated in a few genes;
- several predictors were trained on those labels;
- gene identity dominates its benchmarks (Grimm et al. 2015; Tejura et al. 2024; Harrizi et al. 2026).

De novo variants avoid these problems. Their excess over mutation-rate expectation in affected children defines a likelihood ratio without variant-level labels, and that ratio can be computed within a gene, a window or a domain.

We used the de novo excess to answer three questions:

1. How much predictor evidence survives conditioning on gene and protein position?
2. Is PP3 additive with PM1-like evidence, and does the answer hold on independent truth?
3. Does a label-free alternative calibrate better than the cap?

---

## Methods

### Design and pre-specification
The analysis plan was committed to a version-controlled repository before the full analysis and extended by sixteen amendments (A1–A16). Each amendment was committed before the analysis it governs and states what had been seen when it was written.

- **Not independently registered.** No amendment was deposited in an independent registry before its analysis ran. The timestamps fix when the text was written, not that no results had been seen.
- **The primary DD criterion was not met** (Table 1). The fixed-window position control used for all later tests was adopted after that result, for a reason recorded before those tests.
- **Post hoc analyses are labelled as such.**

The amendment timeline and deviation log are in Supplementary Note 1.

### Cohorts and sites
**Cohorts:**

- **DD:** 31,058 parent–offspring trios (Kaplanis et al. 2020).
- **Autism:** 12,316 probands, with 4,051 unaffected siblings as a negative control (Zhou et al. 2022). The source tables list 13,103 affected individuals with a coding de novo, 12,501 of whom carry an SNV. After restriction to autosomes, deduplication, liftover and site matching, 12,316 carry an analysable SNV. Individuals sharing ≥ 2 identical de novo SNVs with a DD proband were removed (n = 64).
- **Public:** 4,256 trios with congenital heart disease, schizophrenia, Tourette syndrome, epileptic encephalopathy or intellectual disability (Jin et al. 2017; denovo-db v1.6.1). Studies overlapping the other cohorts, and studies reporting fewer than 0.15 synonymous de novos per proband, were excluded; eight duplicates of DD probands were removed.

**Sites.** Coding sites on GENCODE v32 canonical transcripts (17,788 passing QC), excluding sites within 3 bp of an exon boundary. In DD this gave 23,851 missense and 8,222 synonymous de novos.

**Scores and annotations:**

- predictors from dbNSFP 4.1a and 5.4a, and AlphaMissense (Cheng et al. 2023);
- regional missense constraint (Chao et al. 2024) and constrained coding regions (Havrilla et al. 2019);
- AlphaFold DB v6;
- ClinVar 2026-10.

### Likelihood ratios from the de novo excess
Expected counts E come from Roulette mutation rates (Seplyarskiy et al. 2023), scaled to each cohort by its synonymous de novos. The excess is X = O − E.

- **Marginal LR** for a set of sites *b*: LR = (X_b / ΣX) / (E_b / ΣE).
- **Within-unit LR:** each observed de novo is compared with the opportunity profile of its own unit (gene, window or domain).
- **Uncertainty:** a Poisson bootstrap over individuals (1,000 replicates).

**Estimand.** This LR is the evidence a score provides for a de novo missense variant in an affected child, against mutational background. It is not ClinVar's pathogenic-versus-benign LR, so we compare it with the clinical scale only through the bridges described below.

### Retention
Retention r(level) = ln LR(level) / ln LR(marginal), for each predictor's top 2% of scores by mutational opportunity. Because ACMG points are proportional to ln LR, r is the share of a predictor's points that survives conditioning.

- **Levels:** gene; gene × regional-constraint region; gene × 50-residue window (D2); window × constraint (E2).
- **Control:** a random partition of each gene into as many groups as it has windows separates "conditioning on small units" from "conditioning on position".
- **Pre-specified readout:** the share of informative predictors (marginal LR ≥ 2.406) with r(D2) < 0.5, with ≥ 75% called systematic.

### Direct overlap test
**Flags at each site:**

| Flag | Definition |
|---|---|
| PP3 | Published ClinGen Supporting (or Strong) threshold for AlphaMissense, REVEL and VEST4 |
| PM1-like | Regional missense o/e < 0.4 or CCR ≥ 95 |
| PM1-strict | o/e < 0.2 |
| PP2-like | Gene missense Z ≥ 3.09 |

**The test.** R = LR(PP3 ∧ F) / [LR(PP3) · LR(F)] = LR(PP3 | F) / LR(PP3). The PP3 points earned once F applies are ln LR(PP3 | F) / ln 2.406.

**Within-gene R** uses the within-unit estimator, comparing variants only within the same gene, in all genes and in DDG2P monoallelic genes. The difference from the marginal R is the between-gene component.

**At the Strong tier** we also report f, the fraction of the summed log-evidence kept in the joint evidence. ClinGen's cap corresponds to f = 0.67: Strong PP3 (4 points) plus Moderate PM1 (2 points), capped at 4.

### Truth sets independent of the de novo data (exploratory)
1. **ClinVar scale (estimand bridge).** R from ClinVar P/LP vs B/LB LRs, with DD de novo sites removed and a gene-clustered bootstrap.
2. **Expert-panel truth.** ClinGen Evidence Repository variants were reclassified from their applied criteria after removing PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5. P/LP is ≥ 6 remaining points and B/LB is ≤ −1. The truth then depends on neither the predictor nor position.
3. **Multiplexed assays.** ProteinGym v1.3 assays in DDG2P genes, with damaging defined as `DMS_score_bin` = 0.

### Pre-registered label-free model (A16)
**The model** was fitted to DD de novo missense variants only, with no clinical labels.

- **Structure:** a Poisson rate model with a gamma-shrunk gene multiplier, a gamma-shrunk 50-residue-window multiplier, and a fixed effect for each combination of region stratum (other, PM1-like, PM1-strict) and AlphaMissense bin.
- **Output:** each variant's evidence is its de novo-excess LR.

**Comparators.** Four rules were compared, each using computational, gene and position evidence only:

1. uncapped PP3 + PM1 + PP2;
2. the ClinGen cap;
3. the conditional rule (PP3 points given the flag, Supplementary Table 4);
4. the hierarchical model.

Every rule received the same two-parameter logistic recalibration, fitted on ClinVar 2023-01 variants with unanimous P/LP or B/LB classifications. Because those classifications partly use PP3 and PM1, this favours rules 1–3.

**Validation set V1.** V1 contains submissions classified as VUS in ClinVar 2023-01 and reclassified P/LP or B/LB by the same submitter by 2026-10. Truth comes from the ACMG/AMP codes in the submitter's comment, with the same eight criteria removed and negated codes ignored.

**Primary endpoint:** the difference in mean log-loss between the hierarchical rule and the cap on V1, with a gene-clustered bootstrap (2,000 replicates).

**Registration.** The specification (`prereg_fullstudy/A16_HIERARCHICAL_PREREG.md`) was committed and pushed before the model was fitted or V1 was opened. It was not deposited on OSF beforehand.

---

## Results

### Most predictor evidence is information about gene and position
**Retention.** Of 36 dbNSFP 4.1a predictors, 24 were informative in DD. Their median retention was:

| Conditioning | Median retention |
|---|---|
| Within gene | 0.56 |
| Within constraint region | 0.38 |
| Within 50-residue window | 0.33 |
| Window × constraint | 0.32 |

(Figures 2–3; Supplementary Table 1.)

**AlphaMissense.** Its top-2% LR fell from 13.4 (marginal) to 4.10 within gene and 2.19 within window, with r(window) 0.30 [0.28, 0.32].

**The loss reflects position, not the act of conditioning.**

- A random partition of each gene into as many groups as it has windows retained 0.51, close to gene alone (0.56). Real windows on the same sites retained 0.31.
- The paired difference was 0.21 [0.19, 0.23], positive for all 24 predictors.
- InterPro domains retained 0.33, like windows.
- A simulation recovered planted within-gene signals with relative bias ≤ 1% (Supplementary Note 3).

**Pre-specified tests** (Table 1):

- **Primary, not met.** The primary DD test used a position control that, in genes without a constraint split, does not control position. It was not met: 62% of predictors lost most of their evidence, against the 75% criterion.
- **Fixed windows.** With fixed windows the share was 79% in DD, 80% in autism and 93% for 17 newer predictors. The DD and autism intervals include 75%.
- **The continuous estimate was stable:** median retention 0.33 (DD), 0.37 (autism), 0.31 (newer predictors) and 0.28 (public cohort). Per-predictor retention correlated across cohorts (Spearman ρ 0.76).
- **Negative control.** No predictor showed evidence in unaffected siblings.

**Table 1. Pre-specified tests.**

| Cohort / predictors | Informative | r(D2) < 0.5 | 95% CI |
|---|---|---|---|
| DD, 36 (primary; D1 position control) | 24 | 15 (62%): **not met** | — |
| DD, D2 position control | 24 | 19 (79%) | 71–88% |
| Autism, 36 (D2 pre-specified) | 20 | 16 (80%) | 60–95% |
| DD, 17 newer predictors (D2 pre-specified) | 14 | 13 (93%) | 93–100% |
| Public, 36 (descriptive) | 15 | 11 (73%) | — |

Retention was unchanged by ±9% CpG mutation-rate perturbation and by removing all 64,573 ClinVar P/LP sites. It did not differ between predictors trained on clinical labels and label-free ones (Supplementary Note 2).

### PP3 and PM1-like evidence are not additive
**Direct overlap test.** For AlphaMissense PP3 with a PM1-like flag in DD, R was 0.44 [0.41, 0.48] (Table 2, Figure 4):

- the joint LR was 15.0, against 33.9 if the two were multiplied;
- at PM1-like sites PP3 corresponded to LR 1.85 (0.70 points), against LR 3.66 elsewhere.

**Replication.** R was 0.56 [0.44, 0.67] in autism and 0.67 [0.42, 0.95] in the public cohort.

**Dose-response.** Overlap was stronger for the strict PM1 definition (0.37) and weaker for the gene-level PP2-like flag (0.61).

**Table 2. Direct overlap test** (PP3 at Supporting; the full 27-cell table is Supplementary Table 5).

| Predictor × flag | DD | Autism | Public | PP3 points given flag, DD |
|---|---|---|---|---|
| AlphaMissense × PM1-like | 0.44 [0.41, 0.48] | 0.56 [0.44, 0.67] | 0.67 [0.42, 0.95] | 0.70 |
| AlphaMissense × PM1-strict | 0.37 [0.34, 0.41] | 0.50 [0.39, 0.62] | 0.55 [0.33, 0.83] | 0.51 |
| AlphaMissense × PP2-like | 0.61 [0.56, 0.66] | 0.65 [0.52, 0.77] | 0.78 [0.48, 1.09] | 1.07 |
| REVEL × PM1-like | 0.57 [0.52, 0.61] | 0.82 [0.65, 1.00] | 0.78 [0.51, 1.13] | 1.01 |
| VEST4 × PM1-like | 0.67 [0.62, 0.71] | 0.77 [0.62, 0.91] | 0.75 [0.48, 1.02] | 0.79 |

[GCLUST]

**How much of the overlap is between genes.**

| Comparison | DD | Autism |
|---|---|---|
| Within DDG2P monoallelic genes, same-gene comparisons | 0.69 [0.68, 0.70] | 0.74 [0.68, 0.80] |
| Within gene, across all genes | 0.75 [0.71, 0.78] | 0.96 [0.84, 1.08] |

Part of the overlap is between genes. Across all genes in autism, the within-gene overlap is indistinguishable from additivity. In the disease genes where PP3 is applied clinically, a within-gene component remains in both cohorts.

**Where the overlap comes from.**

- **The null side.** AlphaMissense calls 51% of mutation-expected sites at PM1-like positions damaging, against 19% overall. Among excess de novos the rise is only from 79% to 94%. Inside constrained regions the predictor calls most sites damaging whether or not they cause disease.
- **Robustness.** R × PM1-like was below 1 for all 40 informative predictors in DD (median 0.51, range 0.18–0.80; autism median 0.68). Predictors using population constraint or allele frequency overlapped no more than others (p = 0.87). Excluding the 50 genes with the largest de novo excess left R at 0.48 [0.42, 0.53].
- **Structure.** The overlap held in ordered (pLDDT ≥ 70; R 0.53 [0.50, 0.57]) and disordered regions (0.39 [0.28, 0.51]). [BURIAL]

**Strong tier.** With PP3 at its Strong threshold, R for AlphaMissense × PM1-like was 0.33 [0.30, 0.36] in DD. Because Strong-tier LRs are large, the joint evidence still kept most of the summed log-evidence:

| Predictor | f, DD | f, autism | f, public |
|---|---|---|---|
| AlphaMissense | 0.73 [0.71, 0.76] | 0.75 [0.67, 0.81] | 0.75 [0.57, 0.89] |
| REVEL | 0.84 [0.82, 0.87] | 0.90 [0.80, 0.99] | 0.89 [0.71, 1.04] |
| VEST4 | 0.91 [0.88, 0.93] | 0.94 [0.85, 1.02] | 0.94 [0.73, 1.09] |

On the de novo scale, then, the cap (f = 0.67) removes more than the overlap does (6–27%). The next section shows that this does not hold on the clinical scale.

### The overlap is larger on clinical truth
**Table 3** compares R for AlphaMissense × PM1 across truth sets.

| Truth set | Flag | Supporting | Strong |
|---|---|---|---|
| De novo excess, DD (all genes) | PM1-like | 0.44 [0.41, 0.48] | 0.33 [0.30, 0.36] |
| De novo excess, DD, within DDG2P genes | PM1-like | 0.69 [0.68, 0.70] | — |
| ClinVar P/LP vs B/LB, DDG2P genes | PM1-like | 0.32 [0.28, 0.37] | 0.23 [0.18, 0.28] |
| ClinVar P/LP vs B/LB, all genes | PM1-like | 0.20 | 0.15 |
| Expert panel, computational and position criteria removed (2,516 variants, 108 genes) | Panel-applied PM1 | 0.47 [0.22, 1.18] | 0.17 [0.08, 0.61] |
| Multiplexed assays: PPM1D | PM1-like | 0.72 [0.65, 0.78] | 0.71 |
| Multiplexed assays: SHOC2 | PM1-like | 0.86 [0.71, 1.04] | 0.82 |
| Reclassified ClinVar, V1 (87 P, 887 B) | PM1-like | 0.47 [0.22, —] | 0.42 [0.11, —] |

**ClinVar.** Curators who apply both PP3 and PM1 enrich P/LP variants in the joint category, which biases ClinVar R towards 1, against the finding.

**Expert panel.** On truth that excludes computational and position evidence, Strong PP3 earned 1.83 points once the panel applied PM1. Without PM1 its LR was 30.8. The cap leaves Strong PP3 at most 2 points beside Moderate PM1, which matches. Uncapped co-use assigns 4, and the conditional rule derived from de novos assigns 1.1.

**Assays and V1.** Within-gene assay truth could test the overlap in only two genes. In the others the PM1-like flag was constant across the assayed residues. V1 had too few benign variants in the joint category to bound R (13 P and 3 B at Supporting).

### A label-free hierarchical model does not beat the cap (pre-registered)
**The fit.** Within windows, AlphaMissense Strong multiplied the de novo rate 1.6-fold in unconstrained regions, 3.3-fold in PM1-like regions and 4.4-fold in PM1-strict regions (Supplementary Table 6). The same score therefore carries different evidence in different regions, which a flat PP3 + PM1 rule cannot express.

**The validation set.** V1 held 87 P and 887 B variants in 546 genes. DDG2P genes alone held 17 P, below the pre-specified threshold, so the primary analysis used all genes.

**The result (Table 4).** The hierarchical model calibrated worse than the cap: the log-loss difference was +0.147 [+0.106, +0.190], so H1 is refuted. It was also worse than the uncapped and conditional rules, with the PP3-side-only version, and on expert-panel truth (+0.140 [+0.047, +0.232]; not confirmatory).

**Table 4. Pre-registered comparison on reclassified ClinVar variants (V1).**

| Rule | Recalibration slope | Mean log-loss | Δ vs cap [95% CI] |
|---|---|---|---|
| Uncapped | 1.32 | 0.364 | — |
| **ClinGen cap** | 1.37 | **0.365** | reference |
| Conditional | 1.57 | 0.368 | — |
| Hierarchical (label-free) | 0.27 | 0.512 | +0.147 [+0.106, +0.190] |

**Why it failed.**

- **The training slope.** On ClinVar training data the hierarchical score's slope was 0.27. It discriminates ClinVar pathogenicity poorly outside NDD genes, where most V1 variants lie (463 of 546 genes).
- **Gene-specific calibration.** In 142 label-rich DDG2P genes, its gene-specific PP3 evidence did not track ClinVar gene-specific calibration (Spearman 0.01).
- **Rules 1–3 barely differ on V1.** Only 16 V1 variants carried both PP3 and PM1-like, so the choice among flat rules hardly matters at this sample size.

This is the same pattern as the models trained on the de novo excess (Supplementary Note 4). They matched AlphaMissense on de novos but ranked ClinVar variants worse within gene (AUROC 0.88–0.92 vs 0.97). The de novo excess measures overlap well, but it does not provide a classifier that transfers to clinical labels.

### Clinical consequences
**Classification with de novo evidence** (DD; PS2 + PM2 + PM1-like + AlphaMissense PP3):

- Uncapped co-use called 2,438 variants Pathogenic. The de novo excess implies only 84% of that class is disease-associated, below what Pathogenic requires.
- The cap and the conditional rule both moved all 2,438 to Likely Pathogenic.

**Expert-panel classifications** (13,293 in the ClinGen Evidence Repository):

- 1,018 applied both PP3 and PM1.
- The cap changed 24 classifications, because most panels already keep the pair within 4 points. The conditional rule changed 238.

**ClinVar.**

- In DDG2P monoallelic genes, 37.3% of AlphaMissense PP3-positive VUS lie at PM1-like sites, against 13.2% across all genes.
- Without the cap, PP3 + PM1 at Moderate would alone reach Likely Pathogenic for 7,101 VUS in these genes. The cap prevents this.

**Conditional PP3 points** are tabulated by predictor, tier and cohort in Supplementary Table 4. In DD, with a PM1-like flag, AlphaMissense earns 0 points at Supporting and Moderate and 1.1 at Strong.

---

## Discussion

**Main finding.** For de novo missense variants in NDD, most of what missense predictors contribute is information about which gene and which part of the protein a variant falls in. A predictor that scores a hotspot highly is not wrong; it restates what location already says. When PP3 is combined with PM1-like evidence, the joint evidence is substantially less than the sum. The overlap:

- replicates in three cohorts;
- holds within known disease genes;
- appears for all 40 informative predictors;
- survives removal of the strongest genes.

**It is larger on clinical truth.** On ClinVar and on expert-panel truth stripped of computational and position criteria, R at the Strong tier was 0.17–0.23.

**The cap is the best-calibrated rule we tested.** ClinGen's cap (Pejaver et al. 2022; Stenton et al. 2024) fits the data better than the alternatives.

- **Uncapped co-use**, as permitted by CanVIG-UK (Garrett et al. 2021), overstated Strong PP3 by about 2 points on expert-panel truth. With de novo evidence it produced Pathogenic calls in a class only 84% disease-associated.
- **On the de novo scale the cap over-corrects:** it removes 33% of the log-evidence where the overlap removes 6–27%. On expert-panel truth the evidence PP3 keeps beside PM1 (1.83 points) matches the cap's allowance (2 points). We do not claim the cap's exact size is optimal, only that it is closer than the alternatives.
- **The conditional rule** from de novo data under-credits PP3 at the Strong tier on clinical truth (1.1 against 1.83 points).

**The pre-registered alternative failed.** We expected a label-free prior that conditions PP3 on gene and region to beat the cap, and pre-registered that test. It did not.

- **What it captures.** The model shows that PP3 evidence is region-dependent: at the Strong tier the within-window rate ratio is two to three times larger in constrained regions than outside them.
- **Why it doesn't transfer.** It inherits the de novo estimand: high-penetrance, dominant, NDD-specific. Outside constrained regions de novo data give AlphaMissense little within-window signal (1.6-fold at Strong). ClinVar truth across all genes still rewards that signal, and the model discounts it.
- **The same pattern as before.** Three de novo-trained classifiers (Supplementary Note 4) behaved the same way.
- **The lesson** is about using the data, not about the overlap. De novo data are well suited to measuring how criteria interact, and poorly suited to supplying the calibrated evidence itself outside the disorders that generated them.

**Relation to previous work.**

- Chao et al. (2024) calibrated regional missense constraint to ACMG/AMP strengths and showed that constraint plus AlphaMissense enriches case de novos more than either alone. Our results agree: the joint LR exceeds each marginal LR. What fails is the product assumption; it is not that combination adds nothing.
- Gene identity dominates ClinVar benchmarks (Harrizi et al. 2026), and genome-wide calibration masks gene heterogeneity (Tejura et al. 2024). This has motivated gene-level calibration (Badonyi & Marsh 2025; Chen et al. 2026). Our negative A16 result cautions that a label-free route to gene- and region-specific calibration does not yet transfer to clinical classification.
- **What is new here:**
  - a label-free measurement of PP3 + PM1 additivity on the points scale;
  - its decomposition into between-gene and within-gene parts;
  - its replication on independent clinical truth;
  - a pre-registered negative test of the obvious alternative.

**Limitations.**

1. **Pre-specification.** The primary pre-specified DD criterion was not met. The fixed-window control was adopted afterwards, and sixteen amendments followed earlier results. None was deposited in an independent registry before its analysis; A16 was fixed by a public repository push only. We present continuous estimates and the direct test as the findings.
2. **Labels are not entirely absent.** DDG2P strata use curated gene lists, some discovered in these cohorts. The all-genes analyses do not depend on them.
3. **Mutation-rate model.** Residual error is about 9% at CpG transitions; ±9% perturbation leaves the results unchanged.
4. **Estimand.** The de novo LR reflects high-penetrance variants in severe early-onset disorders. The clinical-scale analyses bridge this only partly.
5. **Independent truth sets are small or concentrated.**
   - Expert-panel truth covers 108 genes, and its o/e-flag estimates are imprecise.
   - Assay truth tests overlap in two genes.
   - V1 has few benign variants in the joint category, and its primary stratum fell back to all genes.
6. **Transcripts.** GENCODE v32 canonical transcripts match MANE Select v1.5 for 74.4% of genes (66.5% in DDG2P). Window assignments may differ for the rest.
7. **Ancestry.** The cohorts are predominantly of European ancestry.

**Conclusion.** Computational and hotspot evidence overlap substantially for de novo missense variants, including within known disease genes. The overlap is larger on clinical truth. ClinGen's cap on PP3 + PM1 is the best-calibrated of the rules we tested, and uncapped co-use is not supported. A pre-registered label-free alternative did not improve on the cap.

---

## Tables and figures

- **Table 1.** Pre-specified retention tests.
- **Table 2.** Direct overlap test (summary; full table in Supplementary Table 5).
- **Table 3.** Overlap across truth sets.
- **Table 4.** Pre-registered comparison of rules on reclassified ClinVar variants.
- **Figure 1.** Design: de novo excess, conditioning levels and retention (`fig1_design.png`).
- **Figure 2.** Evidence ladder: top-2% LR by conditioning level (`fig2_evidence_ladder.png`).
- **Figure 3.** Per-predictor retention (`fig3_retention.png`).
- **Figure 4.** Direct overlap test, R with 95% intervals by predictor, flag and cohort (`fig4_overlap.png`).

---

## References

- Badonyi M, Marsh JA. acmgscaler: standardized gene-level variant effect score calibration within the ACMG/AMP framework. 2025 (PMC12496131).
- Bergquist T, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 2025.
- Chao KR, et al. The landscape of regional missense mutational intolerance quantified from 125,748 exomes. bioRxiv 2024. doi:10.1101/2024.04.11.588920.
- Chen et al. Gene- and domain-aware calibration increases the clinical utility of variant effect predictors. Research Square preprint, 2026.
- Cheng J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 2023.
- Garrett A, et al. Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. *J Med Genet* 2021;58:297.
- Grimm DG, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat* 2015.
- Harrizi et al. Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. bioRxiv 2026. doi:10.64898/2026.08.21.746180.
- Havrilla JM, et al. A map of constrained coding regions in the human genome. *Nat Genet* 2019.
- Jin SC, et al. Contribution of rare inherited and de novo variants in 2,871 congenital heart disease probands. *Nat Genet* 2017;49:1593–1601.
- Kaplanis J, et al. Evidence for 28 genetic disorders discovered by combining healthcare and research data. *Nature* 2020.
- Liu X, et al. dbNSFP v4. *Genome Med* 2020;12:103.
- Notin P, et al. ProteinGym: large-scale benchmarks for protein fitness prediction and design. *NeurIPS* 2023.
- Orenbuch R, et al. popEVE. *Nat Genet* 2025.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet* 2022.
- Richards S, et al. *Genet Med* 2015.
- Seplyarskiy V, et al. A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription (Roulette). *Nat Genet* 2023.
- Stenton SL, et al. Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med* 2024.
- Tavtigian SV, et al. *Genet Med* 2018; Tavtigian SV, et al. *Hum Mutat* 2020.
- Tejura M, et al. Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes. *Am J Hum Genet* 2024.
- Turner TN, et al. denovo-db. *Nucleic Acids Res* 2017.
- Zhou X, et al. Integrating de novo and inherited variants in 42,607 autism cases identifies mutations in new moderate-risk genes. *Nat Genet* 2022.
