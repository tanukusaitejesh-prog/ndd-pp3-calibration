# De novo benchmarks overstate missense predictor skill: a gene- and position-adjusted evaluation and the double counting of computational and hotspot evidence (PP3 + PM1)

*Draft 6, 7 October 2026. Supersedes draft 5 (`MANUSCRIPT_v5.md`). Supplementary material: `SUPPLEMENT_v6.md`. Every value comes from a committed output.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** Missense predictors are benchmarked on data pooled across genes. We asked how much of their evidence concerns variant effect within a gene, and what this implies for combining computational (PP3) with hotspot/domain (PM1) evidence.

**Methods.** Likelihood ratios (LRs) were estimated without variant-level labels. They came from the excess of de novo missense variants over mutation-rate expectation in 31,058 developmental-disorder (DD) trios, 12,316 autism probands and 4,256 public trios. LRs were computed marginally, within genes and within 50-residue windows. Additivity of PP3 and PM1-like evidence was measured as R, the joint LR over the product of the separate LRs, and checked on ClinVar and expert-panel truth.

**Results.** Across 39 informative predictors, a median 0.32 of log-evidence survived conditioning on gene and window, and pooling inflated LRs 3.6-fold. Marginal and within-window rankings were unrelated (Spearman −0.18). The four strongest marginal predictors fell to ranks 17–39. The within-window ranking replicated in autism (ρ 0.85). Consequently PP3 and PM1-like evidence were not additive: for AlphaMissense, R was 0.44 [0.41–0.47] in DD, 0.56 [0.47–0.66] in autism and 0.67 [0.48–0.96] in the public cohort, and 0.69–0.75 within DDG2P genes. On expert-panel truth without computational or position criteria, Strong PP3 earned 1.83 points [0.84–3.54] once PM1 applied, against 4 under uncapped co-use. A pre-registered label-free correction did not beat ClinGen's cap.

**Conclusion.** De novo case–control benchmarks overstate predictors' variant-level skill. Gene- and position-adjusted evaluation re-ranks predictors and explains why PP3 double-counts PM1. The evidence is consistent with ClinGen's cap and excludes uncapped co-use.

---

## Introduction

Missense variants are the largest class of variants of uncertain significance. ACMG/AMP classification combines three kinds of evidence for them: computational predictions (PP3), location in a mutational hotspot or critical domain (PM1), and membership of a gene with little benign missense variation (PP2) (Richards et al. 2015).

In the Bayesian points framework these criteria add in log-likelihood-ratio units, with one prior for every gene (Tavtigian et al. 2018, 2020). ClinGen calibrated predictor thresholds against ClinVar (Pejaver et al. 2022; Bergquist et al. 2025). Recognising that PP3 and PM1 can describe the same property, Pejaver et al. capped their sum at Strong, a rule Stenton et al. (2024) restated as no double counting. CanVIG-UK instead permits both at full weight, regarding regional enrichment and in silico prediction as largely different evidence (Garrett et al. 2021).

Predictor evaluations share a blind spot. A predictor that has learned which genes and domains are intolerant separates pathogenic from benign variants pooled across genes, yet may add little once the gene and domain are known. This matters twice:

- **Comparing predictors.** Pooled evaluations reward gene and region knowledge.
- **Classification.** ACMG/AMP credits gene and region separately through PP2 and PM1, so how much PP3 overlaps them has not been measured.

ClinVar is a poor place to measure this:

- its labels are concentrated in a few genes;
- several predictors were trained on those labels;
- gene identity dominates its benchmarks (Grimm et al. 2015; Tejura et al. 2024; Harrizi et al. 2026).

De novo variants avoid these problems, but published de novo evaluations do not exploit this. AlphaMissense, PrimateAI-3D, popEVE and MVP were all assessed by separating case from control de novos pooled across genes (Qi et al. 2021; Cheng et al. 2023; Gao et al. 2023; Orenbuch et al. 2025). Such comparisons reward a predictor for knowing which genes and regions are intolerant, which is information ACMG/AMP already credits through PP2 and PM1. The de novo excess over mutation-rate expectation defines a likelihood ratio without variant-level labels, and that ratio can be computed within a gene or a window.

We used it to answer four questions:

1. How much predictor evidence survives conditioning on gene and protein position?
2. Does a gene- and position-adjusted benchmark rank predictors differently?
3. Is PP3 additive with PM1-like evidence, and does the answer hold on independent truth?
4. Does a label-free alternative calibrate better than ClinGen's cap?

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
- **Uncertainty:** retention uses a Poisson bootstrap over individuals (1,000 replicates). Every R and f interval resamples transcripts instead (gene-clustered, 1,000 replicates; synonymous scale factor fixed), because simulated coverage of individual-level intervals dipped to 0.86 and genes are the natural cluster.

**Estimand.** This LR is the evidence a score provides for a de novo missense variant in an affected child, against mutational background. It is not ClinVar's pathogenic-versus-benign LR, so we compare it with the clinical scale only through the bridges described below.

### Retention
Retention r(level) = ln LR(level) / ln LR(marginal), for each predictor's top 2% of scores by mutational opportunity. Because ACMG points are proportional to ln LR, r is the share of a predictor's points that survives conditioning.

- **Levels:** gene; gene × regional-constraint region; gene × 50-residue window (D2); window × constraint (E2).
- **Control:** a random partition of each gene into as many groups as it has windows separates "conditioning on small units" from "conditioning on position".
- **Pre-specified readout:** the share of informative predictors (marginal LR ≥ 2.406) with r(D2) < 0.5, with ≥ 75% called systematic.

### Gene- and position-adjusted benchmark
Each informative predictor was ranked by its within-window LR for its top 2% of scores. This is the evidence it adds once gene and 50-residue window are fixed. The ranking was compared with the marginal ranking, and within-window LRs were compared across cohorts. Predictors from dbNSFP 5.4a were evaluated on their own site set, and intervals are from the individual-level bootstrap. The full resource is Supplementary Table 7 (`make_benchmark.py`).

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

### Clinical impact
We applied three rules for combining PP3 with PM1 (uncapped; ClinGen cap; conditional, with PP3 points given the flag) in three settings:

- **ClinVar missense VUS**, with PM1 at PM1-like sites;
- **ClinGen Evidence Repository classifications**, using the panels' own applied criteria;
- **DD de novo missense variants**, with PS2 + PM2 added. Here the de novo excess gives the fraction of each resulting class that is disease-associated.

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

- A random partition of each gene into as many groups as it has windows retained 0.51, close to gene alone (0.56). Real windows on the same sites retained 0.31. This control used all analysed sites; on the common site set above, windows retained 0.33.
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

### A gene- and position-adjusted benchmark re-ranks predictors
**Table 2** ranks the 39 informative predictors (25 from dbNSFP 4.1a plus AlphaMissense, and 14 newer) by within-window LR. The ten best are shown, with the four strongest marginal predictors appended (all 39, with replications, in Supplementary Table 7).

**Table 2. Predictors ranked by evidence within 50-residue windows** (DD; top 2% of scores; ranks among 39).

| Predictor | Marginal LR (rank) | Within gene | Within window [95% CI] (rank) | Share retained |
|---|---|---|---|---|
| VEST4 | 5.6 (25) | 4.38 | 2.79 [2.60, 2.98] (1) | 0.59 |
| MutPred | 5.1 (31) | 4.05 | 2.52 [2.33, 2.74] (2) | 0.56 |
| MetaRNN | 6.1 (20) | 4.10 | 2.49 [2.33, 2.65] (3) | 0.51 |
| ClinPred | 5.1 (32) | 3.48 | 2.26 [2.11, 2.42] (4) | 0.50 |
| AlphaMissense | 13.4 (5) | 4.10 | 2.19 [2.09, 2.29] (5) | 0.30 |
| MutationAssessor | 4.7 (33) | 3.52 | 2.16 [1.99, 2.34] (6) | 0.50 |
| Eigen | 5.8 (23) | 3.40 | 2.16 [2.02, 2.30] (7) | 0.44 |
| MisFit_D | 7.5 (13) | 4.04 | 2.12 [1.99, 2.24] (8) | 0.37 |
| PROVEAN | 3.1 (37) | 2.88 | 2.11 [1.93, 2.27] (9) | 0.66 |
| VARITY_ER_LOO | 5.6 (26) | 3.36 | 2.10 [1.97, 2.24] (10) | 0.43 |
| MisFit_S | 16.2 (1) | 3.08 | 1.90 [1.83, 1.97] (17) | 0.23 |
| MutScore | 15.1 (2) | 3.88 | 1.78 [1.72, 1.85] (21) | 0.21 |
| PHACTboost | 15.0 (3) | 2.81 | 1.70 [1.63, 1.76] (26) | 0.20 |
| MPC | 13.6 (4) | 1.95 | 1.28 [1.24, 1.32] (39) | 0.10 |

**Inflation.**
- Pooling across genes inflated the top-2% LR a median 3.6-fold (range 1.5–10.6).
- The median LR fell from 6.08 (marginal) to 2.99 within gene and 1.78 within window.

**The rankings were unrelated** (Spearman −0.18, p = 0.28).
- The four strongest marginal predictors, MisFit_S, MutScore, PHACTboost and MPC (marginal LR 13.6–16.2), ranked 17th, 21st, 26th and 39th within windows.
- VEST4, MutPred, MetaRNN and ClinPred, ranked 20th–32nd marginally, led within windows.
- popEVE, evaluated in DDD only, fell from 17.9 to 1.76.
- AlphaMissense was the only predictor in the top five on both scales.
- Only VEST4's within-window LR had a lower 95% bound above the ClinGen Supporting LR (2.406).

**Replication.** Within-window LRs correlated across predictors between DD and autism (Spearman 0.85, 25 predictors) and, more weakly, with the smaller public cohort (0.49).

The predictors that look strongest in pooled de novo evaluations are largely those that encode gene and regional constraint. That is the information that also underlies PM1.

**De novo data cannot teach within-window effect** (pre-registered in A17, task 13; Supplementary Note 7).

- **The model.** A Poisson gradient-boosted model was fitted with an offset that fixes each window's rate, so it could learn only within-window variation. Its inputs were AlphaMissense, AlphaFold confidence and burial, phyloP conservation and the amino-acid change; it was trained on DD with genes held out.
- **Within windows it added nothing to AlphaMissense.** The held-out log-likelihood gain was +0.0008 [−0.0006, +0.0023] nats per de novo, against +0.025 for AlphaMissense over mutation rate alone.
- **It was worse on clinical and functional truth.** Within-gene ClinVar AUROC in DDG2P genes was 0.894 against 0.967. Assay correlation was lower in 7 of 8 assays.
- **Reading.** The within-window signal in de novo data is small enough to rank existing predictors but too small to train a better one.

### PP3 and PM1-like evidence are not additive
**Direct overlap test.** For AlphaMissense PP3 with a PM1-like flag in DD, R was 0.44 [0.41, 0.47] (Table 3, Figure 4):

- the joint LR was 15.0, against 33.9 if the two were multiplied;
- at PM1-like sites PP3 corresponded to LR 1.85 (0.70 points), against LR 3.66 elsewhere.

**Replication.** R was 0.56 [0.47, 0.66] in autism and 0.67 [0.48, 0.96] in the public cohort.

**Dose-response.** Overlap was stronger for the strict PM1 definition (0.37) and weaker for the gene-level PP2-like flag (0.61).

**Table 3. Direct overlap test** (PP3 at Supporting; gene-clustered 95% intervals; all cells, both tiers, in Supplementary Table 5).

| Predictor × flag | DD | Autism | Public | PP3 points given flag, DD |
|---|---|---|---|---|
| AlphaMissense × PM1-like | 0.44 [0.41, 0.47] | 0.56 [0.47, 0.66] | 0.67 [0.48, 0.96] | 0.70 |
| AlphaMissense × PM1-strict | 0.37 [0.35, 0.40] | 0.50 [0.42, 0.61] | 0.55 [0.38, 0.83] | 0.51 |
| AlphaMissense × PP2-like | 0.61 [0.57, 0.65] | 0.65 [0.55, 0.76] | 0.78 [0.57, 1.07] | 1.07 |
| REVEL × PM1-like | 0.57 [0.53, 0.61] | 0.82 [0.69, 0.99] | 0.78 [0.54, 1.15] | 1.01 |
| VEST4 × PM1-like | 0.67 [0.63, 0.70] | 0.77 [0.64, 0.91] | 0.75 [0.54, 1.04] | 0.79 |

Gene-clustered intervals were as narrow as individual-level ones, so the overlap is not carried by a few genes.

**How much of the overlap is between genes.**

| Comparison (AlphaMissense × PM1-like, Supporting) | DD | Autism | Public |
|---|---|---|---|
| Within DDG2P monoallelic genes, same-gene comparisons | 0.69 [0.65, 0.73] | 0.74 [0.67, 0.82] | 0.75 [0.62, 0.89] |
| Within gene, across all genes | 0.75 [0.71, 0.78] | 0.96 [0.84, 1.10] | 0.96 [0.74, 1.28] |

Part of the overlap is between genes. Across all genes in autism, the within-gene overlap is indistinguishable from additivity. In the disease genes where PP3 is applied clinically, a within-gene component remains in all three cohorts.

**Where the overlap comes from.**

- **The null side.** AlphaMissense calls 51% of mutation-expected sites at PM1-like positions damaging, against 19% overall. Among excess de novos the rise is only from 79% to 94%. Inside constrained regions the predictor calls most sites damaging whether or not they cause disease.
- **Robustness.** R × PM1-like had an upper 95% bound below 1 for all 40 informative predictors in DD (median 0.51, range 0.18–0.80; autism median 0.68). Predictors using population constraint or allele frequency overlapped no more than others (p = 0.87). Excluding the 50 genes with the largest de novo excess left R at 0.48 [0.42, 0.53].
- **Structure.** The overlap held in ordered (pLDDT ≥ 70; R 0.53 [0.50, 0.57]) and disordered regions (0.39 [0.28, 0.51]). It also held at buried and exposed sites (AlphaFold Cα neighbours within 10 Å ≥ 9 vs < 9; R 0.55 [0.52, 0.59] and 0.36 [0.29, 0.43]), and in autism (0.57 and 0.58). PP3 is not independent of PM1-like evidence in any structural context we examined (Supplementary Note 5).

**Strong tier.** With PP3 at its Strong threshold, R for AlphaMissense × PM1-like was 0.33 [0.31, 0.35] in DD (within DDG2P genes 0.65 [0.60, 0.69]). Because Strong-tier LRs are large, the joint evidence still kept most of the summed log-evidence:

| Predictor | f, DD | f, autism | f, public |
|---|---|---|---|
| AlphaMissense | 0.73 [0.72, 0.75] | 0.75 [0.71, 0.79] | 0.75 [0.66, 0.85] |
| REVEL | 0.84 [0.82, 0.86] | 0.90 [0.82, 0.99] | 0.89 [0.74, 1.04] |
| VEST4 | 0.91 [0.88, 0.93] | 0.94 [0.86, 1.02] | 0.94 [0.80, 1.11] |

On the de novo scale, then, the cap (f = 0.67) removes more than the overlap does (6–27%). The next section shows that this does not hold on the clinical scale.

### The overlap is larger on clinical truth
**Table 4** compares R for AlphaMissense × PM1 across truth sets.

| Truth set | Flag | Supporting | Strong |
|---|---|---|---|
| De novo excess, DD (all genes) | PM1-like | 0.44 [0.41, 0.47] | 0.33 [0.31, 0.35] |
| De novo excess, DD, within DDG2P genes | PM1-like | 0.69 [0.65, 0.73] | 0.65 [0.60, 0.69] |
| ClinVar P/LP vs B/LB, DDG2P genes | PM1-like | 0.32 [0.28, 0.37] | 0.23 [0.18, 0.28] |
| ClinVar P/LP vs B/LB, all genes | PM1-like | 0.20 | 0.15 |
| Expert panel, computational and position criteria removed (2,516 variants, 108 genes) | Panel-applied PM1 | 0.47 [0.22, 1.18] | 0.17 [0.08, 0.61] |
| Multiplexed assays: PPM1D | PM1-like | 0.72 [0.65, 0.78] | 0.71 |
| Multiplexed assays: SHOC2 | PM1-like | 0.86 [0.71, 1.04] | 0.82 |
| Reclassified ClinVar, V1 (87 P, 887 B) | PM1-like | Not estimable (3 B/LB in the joint cell) | Not estimable (1 B/LB) |

**ClinVar.** Curators who apply both PP3 and PM1 enrich P/LP variants in the joint category, which biases ClinVar R towards 1, against the finding.

**Expert panel.** On truth that excludes computational and position evidence, Strong PP3 earned 1.83 points [0.84, 3.54] once the panel applied PM1 (gene-clustered bootstrap). Without PM1 its LR was 30.8. The interval includes the cap's allowance (at most 2 points beside Moderate PM1) and excludes the 4 points of uncapped co-use. It also includes the 1.1 points of the conditional rule derived from de novos, so these data cannot separate the cap from the conditional rule.

**Population-frequency sensitivity.** BA1, BS1, BS2 and PM2 were also removed from the truth points, so that frequency evidence cannot drive the truth. This left 1,540 variants (1,201 P/LP, 339 B/LB). The overlap was, if anything, stronger: R with the panel's PM1 was 0.34 [0.17, 1.55] at Supporting and 0.12 at Strong. With so few benign variants the Strong-tier interval was unbounded.

**Assays and V1.** Within-gene assay truth could test the overlap in only two genes. In the others the PM1-like flag was constant across the assayed residues. V1 had too few benign variants in the joint category to estimate R (13 P and 3 B at Supporting), so no point estimate is reported.

### A label-free hierarchical model does not beat the cap (pre-registered)
**The fit.** Within windows, AlphaMissense Strong multiplied the de novo rate 1.6-fold in unconstrained regions, 3.3-fold in PM1-like regions and 4.4-fold in PM1-strict regions (Supplementary Table 6). The same score therefore carries different evidence in different regions, which a flat PP3 + PM1 rule cannot express.

**The validation set.** V1 held 87 P and 887 B variants in 546 genes. DDG2P genes alone held 17 P, below the pre-specified threshold, so the primary analysis used all genes.

**The result (Table 5).** The hierarchical model calibrated worse than the cap: the log-loss difference was +0.147 [+0.106, +0.190], so H1 was not supported. It was also worse than the uncapped and conditional rules, with the PP3-side-only version, and on expert-panel truth (+0.140 [+0.047, +0.232]; not confirmatory).

**Table 5. Pre-registered comparison on reclassified ClinVar variants (V1).**

| Rule | Recalibration slope | Mean log-loss | Δ vs cap [95% CI] |
|---|---|---|---|
| Uncapped | 1.32 | 0.364 | — |
| **ClinGen cap** | 1.37 | **0.365** | reference |
| Conditional | 1.57 | 0.368 | — |
| Hierarchical (label-free) | 0.27 | 0.512 | +0.147 [+0.106, +0.190] |

**Why it failed.**

- **The training slope.** On ClinVar training data the hierarchical score's slope was 0.27. It discriminates ClinVar pathogenicity poorly outside NDD genes, where most V1 variants lie (463 of 546 genes).
- **It failed even in its home domain.** In 142 label-rich DDG2P genes, the NDD genes the model was fitted for, its gene-specific PP3 evidence did not track ClinVar gene-specific calibration (Spearman 0.01). The failure is therefore not only a transfer problem to non-NDD genes.
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

**Main finding.** De novo case–control benchmarks reward missense predictors for knowing which genes and regions matter. Once gene and position are held fixed:

- about two-thirds of the evidence disappears;
- the ranking of predictors changes completely;
- the strongest predictors in pooled evaluations, several of them built partly from constraint, add among the least within a window.

A predictor that scores a hotspot highly is not wrong; it restates what location already says. Pooled de novo benchmarks therefore overstate variant-level skill, and the gene- and position-adjusted ranking is the relevant one wherever PP3 is applied alongside gene- or region-level criteria. The same redundancy is the mechanism of PP3 + PM1 double counting. When PP3 is combined with PM1-like evidence, the joint evidence is substantially less than the sum. The overlap:

- replicates in three cohorts;
- holds within known disease genes;
- appears for all 40 informative predictors;
- survives removal of the strongest genes.

**It is larger on clinical truth.** On ClinVar and on expert-panel truth stripped of computational and position criteria, R at the Strong tier was 0.15–0.23.

**The evidence is consistent with the cap and excludes uncapped co-use.** ClinGen's cap (Pejaver et al. 2022; Stenton et al. 2024) is consistent with the data; uncapped co-use is not.

- **Uncapped co-use**, as permitted by CanVIG-UK (Garrett et al. 2021), overstated Strong PP3 by about 2 points on expert-panel truth. With de novo evidence it produced Pathogenic calls in a class only 84% disease-associated.
- **On the de novo scale the cap over-corrects:** it removes 33% of the log-evidence where the overlap removes 6–27%. On expert-panel truth the evidence PP3 keeps beside PM1 (1.83 points [0.84, 3.54]) is consistent with the cap's allowance (2 points). We do not claim the cap's exact size is optimal: the interval also includes the conditional rule.
- **The conditional rule** from de novo data assigns 1.1 points at the Strong tier, inside the expert-panel interval. The data cannot distinguish it from the cap.

**The pre-registered alternative failed.** We expected a label-free prior that conditions PP3 on gene and region to beat the cap, and pre-registered that test. It did not.

- **What it captures.** The model shows that PP3 evidence is region-dependent: at the Strong tier the within-window rate ratio is two to three times larger in constrained regions than outside them.
- **Why it fails.** It inherits the de novo estimand: high-penetrance, dominant, NDD-specific. Outside constrained regions de novo data give AlphaMissense little within-window signal (1.6-fold at Strong). ClinVar truth still rewards that signal, and the model discounts it. The failure is not only one of transfer: even in NDD genes, its gene-specific PP3 evidence did not track ClinVar's (Spearman 0.01).
- **The same pattern as before.** Three de novo-trained classifiers (Supplementary Note 4) and a pre-registered within-window model (task 13) behaved the same way: de novo data rank predictors but cannot train a better one.
- **The lesson** is about using the data, not about the overlap. De novo data are well suited to measuring how criteria interact, and poorly suited to supplying the calibrated evidence itself outside the disorders that generated them.

**Relation to previous work.**

- Chao et al. (2024) calibrated regional missense constraint to ACMG/AMP strengths and report that constraint and predictor scores together enrich case de novos beyond either alone [to verify against the source before submission]. Our results agree: the joint LR exceeds each marginal LR (15.0 against 4.2 and 8.1 for AlphaMissense × PM1-like in DD). What fails is the product assumption, not the value of combining the two.
- Gene identity dominates ClinVar benchmarks (Harrizi et al. 2026), and genome-wide calibration masks gene heterogeneity (Tejura et al. 2024). This has motivated gene-level calibration (Badonyi & Marsh 2025; Chen et al. 2026). Our negative A16 result cautions that a label-free route to gene- and region-specific calibration does not yet transfer to clinical classification.
- **Predictor benchmarks.** Published de novo evaluations of AlphaMissense, PrimateAI-3D, popEVE and MVP compare cases with controls pooled across genes. Gene-level heterogeneity and circularity of predictor performance are documented on ClinVar (Grimm et al. 2015; Tejura et al. 2024; Harrizi et al. 2026). To our knowledge, no previous benchmark evaluates predictors within gene and window on label-free de novo data.
- **What is new here:**
  - a label-free, gene- and position-adjusted benchmark of 39 predictors, released as a resource;
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
6. **Transcripts.** Of GENCODE v32 canonical transcripts in MANE genes, 74.4% match MANE Select v1.5 (66.5% in DDG2P genes). Window assignments may differ for the rest.
7. **Ancestry.** The cohorts are predominantly of European ancestry.
8. **Benchmark scope.** The benchmark scores each predictor's top 2% of scores against de novo excess in severe developmental disorders. Predictors that are informative mainly at lower thresholds, or for other phenotypes, may rank differently.

**Conclusion.** Computational and hotspot evidence overlap substantially for de novo missense variants, including within known disease genes. The overlap is larger on clinical truth. On expert-panel truth the evidence PP3 keeps beside PM1 is consistent with ClinGen's cap and excludes uncapped co-use. A pre-registered label-free alternative did not improve on the cap.

---

## Tables and figures

- **Table 1.** Pre-specified retention tests.
- **Table 2.** Predictors ranked by within-window evidence.
- **Table 3.** Direct overlap test (summary; full table in Supplementary Table 5).
- **Table 4.** Overlap across truth sets.
- **Table 5.** Pre-registered comparison of rules on reclassified ClinVar variants.
- **Figure 1.** Design: de novo excess, conditioning levels and retention (`fig1_design.png`).
- **Figure 2.** Evidence ladder: top-2% LR by conditioning level (`fig2_evidence_ladder.png`).
- **Figure 3.** Per-predictor retention (`fig3_retention.png`).
- **Figure 4.** Direct overlap test, R with 95% intervals by predictor, flag and cohort (`fig4_overlap.png`).

**Supplementary (`SUPPLEMENT_v4.md`):**

- **Notes:**
  - Note 1: pre-specification and deviations;
  - Note 2: robustness;
  - Note 3: estimator simulation;
  - Note 4: models trained on the de novo excess;
  - Note 5: full overlap results;
  - Note 6: pre-registered hierarchical model;
  - Note 7: pre-registered within-window model (task 13).
- **Tables:**
  - Table 1: per-predictor retention;
  - Table 2: training-data ledger (`aim3_training_ledger.tsv`);
  - Table 3: ClinVar sensitivity;
  - Table 4: conditional PP3 points (`supp_table4_conditional_pp3.tsv`);
  - Table 5: full direct overlap test;
  - Table 6: hierarchical model calibration table (`results_fullstudy/a16_calibration_table.tsv`);
  - Table 7: gene- and position-adjusted predictor benchmark, all 39 predictors (`supp_table7_benchmark.tsv`).
- **Data 1:** per-variant ClinVar PP3 + PM1 points.

---

## References

- Badonyi M, Marsh JA. acmgscaler: standardized gene-level variant effect score calibration within the ACMG/AMP framework. 2025 (PMC12496131).
- Bergquist T, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 2025.
- Chao KR, et al. The landscape of regional missense mutational intolerance quantified from 125,748 exomes. bioRxiv 2024. doi:10.1101/2024.04.11.588920.
- Chen et al. Gene- and domain-aware calibration increases the clinical utility of variant effect predictors. Research Square preprint, 2026.
- Cheng J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 2023.
- Gao H, et al. The landscape of tolerated genetic variation in humans and primates (PrimateAI-3D). *Science* 2023.
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
- Qi H, et al. MVP predicts the pathogenicity of missense variants by deep learning. *Nat Commun* 2021;12:510.
- Richards S, et al. *Genet Med* 2015.
- Seplyarskiy V, et al. A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription (Roulette). *Nat Genet* 2023.
- Stenton SL, et al. Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med* 2024.
- Tavtigian SV, et al. *Genet Med* 2018; Tavtigian SV, et al. *Hum Mutat* 2020.
- Tejura M, et al. Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes. *Am J Hum Genet* 2024.
- Turner TN, et al. denovo-db. *Nucleic Acids Res* 2017.
- Zhou X, et al. Integrating de novo and inherited variants in 42,607 autism cases identifies mutations in new moderate-risk genes. *Nat Genet* 2022.
