# Computational, regional and gene-level evidence are stacked but not independent: a label-free audit of missense criteria and of expert classification practice

*Draft 9.1, 7 October 2026. Supplementary material: `SUPPLEMENT_v9.md`.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** ACMG/AMP classification sums computational (PP3), hotspot/domain (PM1) and gene-level (PP2) evidence as if independent. ClinGen caps PP3 + PM1 at Strong. We measured, without variant-level labels, whether these criteria are independent and how expert panels combine them.

**Methods.** Likelihood ratios came from the excess of de novo missense variants over mutation-rate expectation in 31,058 developmental-disorder trios, 12,316 autism probands and 4,256 public trios. Dependence was also measured on 6,177 expert-panel and 87,873 laboratory missense classifications re-scored without computational or location criteria. The primary pre-specified criterion was not met; the practice analyses were not pre-specified.

**Results.** Most predictor evidence concerned gene and position: a median 0.32 of log-evidence survived conditioning on 50-residue windows. No pair of criteria was independent (e.g. PP3 × PM1-like R 0.44 [0.41–0.47]). With Strong PP3, the three criteria nominally sum to 7 points (5 when capped) against 4.05 points of de novo evidence. At Supporting the sum, 4, approximately matched the evidence (3.65). Expert panels applied PP3 beside PM1 almost always below Strong (783 of 787), where the cap does not act. There, Moderate PP3 added no detectable evidence (0.02 points [−1.69–0.89]), whereas Supporting PP3 was approximately calibrated (0.81 [0.24–1.68]). This held in held-out genes and later classifications. Laboratory submissions replicated the dependence (R 0.18 [0.14–0.25]) and the Supporting calibration.

**Conclusion.** Location-type criteria overlap, and the overlap that matters in expert practice lies above Supporting PP3, below the reach of the current cap. Limiting PP3 to Supporting when PM1 applies is consistent with these data.

---

## Introduction

Missense variants are the largest class of variants of uncertain significance. Under ACMG/AMP, their classification draws on three kinds of evidence (Richards et al. 2015):

- computational prediction (PP3);
- location in a hotspot or critical domain (PM1);
- membership of a gene with little benign missense variation (PP2).

The Bayesian points framework adds these in log-likelihood-ratio units, assuming they are independent (Tavtigian et al. 2018, 2020). ClinGen calibrated predictor thresholds against ClinVar (Pejaver et al. 2022; Bergquist et al. 2025). Recognising that PP3 and PM1 can describe the same property, it capped their sum at Strong (Pejaver et al. 2022; Stenton et al. 2024). CanVIG-UK permits both at full weight (Garrett et al. 2021). UK guidance cautions against counting constraint in both PM1 and PP2 (Durkie et al. 2024), and some expert panels restrict PP2 for this reason. The size of that dependence has not been measured, and several ClinGen panels still apply both.

The question is how much a predictor adds once gene and position are known. A predictor that has learned which genes and domains are intolerant separates pathogenic from benign variants pooled across genes, yet may add little within them. Pooled evaluation rewards that knowledge twice:

- **in benchmarks**, which rank predictors on data pooled across genes;
- **in classification**, which credits gene and region separately through PP2 and PM1.

ClinVar is a poor place to measure this. Its labels are concentrated in few genes, several predictors were trained on them, and gene identity dominates its benchmarks (Grimm et al. 2015; Tejura et al. 2024; Harrizi et al. 2026). Published de novo evaluations avoid the labels but pool across genes (Qi et al. 2021; Cheng et al. 2023; Gao et al. 2023; Orenbuch et al. 2025).

The excess of de novo missense variants over mutation-rate expectation defines a likelihood ratio without variant-level labels, and that ratio can be computed within any unit: a gene, a window or a domain. We used it to measure the evidence predictors add beyond gene and position, to bound what is attainable, to check the resulting ranking against independent truth, to test the independence of PP3, PM1 and PP2, and to measure how expert panels stack them in practice.

---

## Methods

### Pre-specification
The analysis plan and seventeen amendments (A1–A17) were committed to a public version-controlled repository, each before the analysis it governs. None was deposited in an independent registry before its analysis ran. Each amendment records what had been seen when it was written.

- **The primary pre-specified criterion was not met** (Table 1). It required that ≥ 75% of predictors lose most of their evidence under a position control (D1). In genes without a regional-constraint split, D1 does not control position. A fixed-window control (D2) was adopted after that result and pre-specified for all later tests.
- **Not pre-specified:** the benchmark, the oracle, the PM1 × PP2 and triple tests, the predictor-validation analyses and the expert-practice audit were added after the main analyses. They are descriptive.
- **The amendment timeline and deviation log** are in Supplementary Note 1.

### Cohorts and sites
- **DD:** 31,058 parent–offspring trios (Kaplanis et al. 2020).
- **Autism:** 12,316 probands with an analysable de novo SNV, and 4,051 unaffected siblings as a negative control (Zhou et al. 2022). The source tables list 13,103 affected individuals with a coding de novo, 12,501 of them with an SNV. The rest were removed by autosome restriction, liftover and site matching. Sixty-four individuals sharing ≥ 2 de novo SNVs with a DD proband were also removed.
- **Public:** 4,256 trios with congenital heart disease, schizophrenia, Tourette syndrome, epileptic encephalopathy or intellectual disability (Jin et al. 2017; denovo-db v1.6.1), excluding studies that overlap the other cohorts.
- **Sites:** coding sites on GENCODE v32 canonical transcripts (17,788 passing QC), excluding sites within 3 bp of an exon boundary. In DD this gave 23,851 missense and 8,222 synonymous de novos.
- **Predictors:** dbNSFP 4.1a and 5.4a, and AlphaMissense (Cheng et al. 2023).
- **Annotations:** regional missense constraint (Chao et al. 2024), constrained coding regions (Havrilla et al. 2019) and AlphaFold DB v6.

### Likelihood ratios from the de novo excess
**Expected counts and excess.** Expected counts E come from Roulette mutation rates (Seplyarskiy et al. 2023), scaled to each cohort by its synonymous de novos. The excess is X = O − E.

**Estimators.**
- **Marginal LR** for a set of sites *b*: LR = (X_b/ΣX)/(E_b/ΣE).
- **Within-unit LR:** each observed de novo is compared with the mutational-opportunity profile of its own unit: LR_b = X_b / (Σ_u O_u E_ub/E_u − E_b).
- **Retention** r = ln LR(within) / ln LR(marginal), evaluated at each predictor's top 2% of scores by mutational opportunity (also 1% and 10%). Because ACMG points are proportional to ln LR, r is the share of a predictor's points that is incremental to the unit.

**Uncertainty.** Intervals resample transcripts (gene-clustered bootstrap, 1,000 replicates), with the synonymous scale factor fixed. Shares of predictors additionally resample predictors.

**Estimand.** This is the evidence a score gives for a de novo missense variant in an affected child against mutational background. It is not ClinVar's pathogenic-versus-benign LR.

**Controls:**
- a random partition of each gene into as many groups as it has windows;
- 25-, 50- and 100-residue windows;
- InterPro domains;
- unaffected siblings;
- a null predictor (a hash of the variant key);
- a simulation with planted signals (Supplementary Note 3).

### Oracle
A predictor can retain little within windows either because it lacks variant-level information or because little within-window evidence is attainable. To separate these, we computed the evidence available to a perfect predictor, under the excess estimators used throughout:

- **The world.** Sites are pathogenic or not; pathogenic sites carry de novo rate k times their mutation rate. Each window's excess is taken from a fitted empirical-Bayes rate model (gene, window, region stratum and AlphaMissense bin; DD only), so its pathogenic opportunity is X_w/(k − 1).
- **The oracle predictor.** It flags exactly the pathogenic sites, padded at random to 2% of opportunity.
- **k is unknown**, and the de novo data do not identify it. Results are reported for k = 10–300. The oracle is used only to show that retention is not a measure of skill, not to place predictors on a skill scale.

### Benchmark and its validation
**Ranking.** Informative predictors (marginal top-2% LR ≥ 2.406) were ranked by within-window LR. Rank intervals come from the gene-clustered bootstrap. Robustness was checked at 1%, 2% and 10% thresholds, on a common site set scored by every predictor, and within DDG2P monoallelic genes.

**Validity.** Across predictors, within-window and marginal LR were correlated with:
- pooled and within-gene ClinVar AUROC, for P/LP versus B/LB missense variants with DD de novo sites removed, in genes with ≥ 3 of each;
- mean within-assay Spearman correlation with damage across eight ProteinGym assays in DDG2P genes.

Intervals come from a predictor-level bootstrap.

### Independence of PP3, PM1 and PP2
**Flags:**

| Flag | Definition |
|---|---|
| PP3 | ClinGen Supporting or Strong threshold (AlphaMissense, REVEL, VEST4) |
| PM1-like | Regional missense o/e < 0.4 or CCR ≥ 95 |
| PM1-strict | o/e < 0.2, the level Chao et al. calibrated to Moderate |
| PP2-like | Gene missense Z ≥ 3.09 |

**Measures.**
- **Dependence:** for criteria A and B, R = LR(A ∧ B) / [LR(A) · LR(B)]. R = 1 is independence, as the points framework assumes; R < 1 means summing them overstates the joint evidence.
- **Within-gene R** uses the within-unit estimator, in all genes and in DDG2P monoallelic genes.
- **At the Strong tier**, f = ln LR(joint) / [ln LR(PP3) + ln LR(F)] is the fraction of summed log-evidence kept. ClinGen's cap on Strong PP3 (4 points) plus Moderate PM1 (2) corresponds to f = 0.67, so PM1-strict is the primary flag for that comparison.

### Truth sets independent of the de novo data
1. **ClinVar P/LP versus B/LB**, with DD de novo sites removed.
2. **Expert-panel truth.** ClinGen Evidence Repository classifications were re-scored from their applied criteria after removing PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5. Truth was P/LP at ≥ 6 remaining points and B/LB at ≤ −1. A sensitivity analysis also removed the frequency criteria BA1, BS1, BS2 and PM2.
3. **ProteinGym assays** in DDG2P genes.
4. **ClinVar submissions reclassified** from VUS (2023-01) to P/LP or B/LB (2026-10), with truth from the submitter's codes.

### Stacking in expert practice
**Data.** All missense classifications in the ClinGen Evidence Repository (export of 2026-10-07; 6,177 after removing retracted entries and non-missense variants).

**How the criteria were measured.**
- We counted co-application of PP3, PM1 and PP2, and the strength at which each was applied.
- On truth re-scored without computational and location criteria, we estimated the dependence of the panels' own PP3 and PM1, and the PP3 points earned once PM1 applies (gene-clustered bootstrap).
- On the de novo scale, the triple dependence was R3 = LR(PP3 ∧ PM1 ∧ PP2) / [LR(PP3) · LR(PM1) · LR(PP2)].

**Alternative rules.** Every classification was re-scored on the points scale under five rules:
1. the ClinGen cap;
2. a cap on all three location-type criteria;
3. PP2 not added when PM1 applies;
4. PP3 limited to Supporting when PM1 applies;
5. PP3 credited at its measured value beside PM1.

**Robustness and out-of-sample tests.**
- **Truth thresholds:** P/LP at ≥ 4, 6 or 8 points; B/LB at ≤ −1 or −3.
- **Frequency criteria** (BA1, BS1, BS2, PM2) also removed from the truth points.
- **Gene-split cross-validation:** 200 random halves of the genes.
- **Temporal split:** credit estimated on classifications approved before 2025 and tested on those approved in 2025–2026.

---

## Results

### Most predictor evidence is information about gene and position
**Retention** (design in Figure 1). Across the 39 informative predictors in DD, the median top-2% LR was 6.06 marginally, 3.00 within gene and 1.79 within 50-residue windows. The median retention was 0.32, and pooling inflated LRs a median 3.6-fold. Retention was similar at the 1% and 10% thresholds (0.33, 0.29) and on a common site set (0.34, 0.32, 0.30). The null predictor gave LRs of 0.84–1.11 throughout.

**The loss reflects position, not the act of conditioning.**
- A random partition of each gene into as many groups as it has windows retained 0.51, close to gene alone (0.56), whereas real windows retained 0.31 (24 dbNSFP 4.1a predictors, all sites).
- Retention fell with finer windows: 0.39 at 100 residues, 0.33 at 50 and 0.27 at 25. Equal-sized random partitions kept 0.53, 0.52 and 0.48.
- InterPro domains retained 0.33.
- No predictor showed evidence in unaffected siblings.
- Removing all 64,573 ClinVar P/LP sites and perturbing CpG mutation rates by ±9% left retention unchanged (Supplementary Note 2).

**Table 1. Pre-specified tests** (share of informative predictors with within-window retention < 0.5; exact binomial intervals, which ignore correlation between predictors).

| Cohort and predictors | Position control | Informative | r < 0.5 | 95% CI |
|---|---|---|---|---|
| DD, dbNSFP 4.1a (primary) | D1 | 24 | 15 (62%): **not met** | 41–81% |
| DD, dbNSFP 4.1a | D2 | 24 | 19 (79%) | 58–93% |
| Autism, dbNSFP 4.1a (pre-specified D2) | D2 | 20 | 16 (80%) | 56–94% |
| DD, 17 newer predictors (pre-specified D2) | D2 | 14 | 13 (93%) | 66–100% |
| Public, dbNSFP 4.1a (descriptive) | D2 | 15 | 11 (73%) | 45–92% |

Across all 39 predictors, with genes and predictors both resampled, the share with r < 0.5 was 0.87 [0.69, 0.97].

### How much could a perfect predictor retain?
A perfect predictor also loses evidence within windows, because the de novo excess is concentrated in some windows: each window with excess holds pathogenic and benign sites in proportion to its excess. How much it loses depends on how strongly pathogenic sites are enriched (Table 2).

**Table 2. Evidence available to a perfect predictor** (top 2% of opportunity; DD).

| Enrichment of pathogenic sites, k | Pathogenic share of missense opportunity | Marginal LR | Within-gene LR | Within-window LR | Retention within window |
|---|---|---|---|---|---|
| 10 | 3.5% | 28.4 | 5.7 | 2.8 | 0.30 |
| 30 | 1.1% | 50.0 | 15.6 | 8.2 | 0.54 |
| 100 | 0.3% | 50.0 | 30.4 | 20.1 | 0.77 |
| 300 | 0.1% | 50.0 | 41.2 | 33.5 | 0.90 |

Across plausible k, a perfect predictor retains 0.30–0.90, a range that contains the observed median (0.32; AlphaMissense 0.30 [0.27, 0.33]; VEST4 0.59). Because k is not identified, retention cannot be read as a fraction of attainable skill. It measures evidence incremental to location, which is the quantity that matters when PP3 is combined with PM1 or PP2.

### A gene- and position-adjusted benchmark (resource)
**Figure 2** and **Table 3** compare pooled and within-window rankings. We present the benchmark as a mechanism and a resource; the re-ranking is not a headline finding.

**Table 3. Predictors ranked by within-window evidence** (DD, top 2%; gene-clustered 95% intervals; the ten best and the four strongest pooled predictors; all 39 in Supplementary Table 7).

| Predictor | Pooled LR (rank) | Within-window LR [95% CI] | Rank [95% interval] | Retention |
|---|---|---|---|---|
| VEST4 | 5.75 (25) | 2.83 [2.44, 3.28] | 1 [1, 2] | 0.59 |
| MutPred | 5.15 (31) | 2.53 [2.16, 2.91] | 2 [1, 7] | 0.57 |
| MetaRNN | 6.06 (20) | 2.45 [2.14, 2.77] | 3 [2, 7] | 0.50 |
| ClinPred | 5.07 (32) | 2.25 [1.97, 2.53] | 4 [2, 13] | 0.50 |
| AlphaMissense | 13.25 (5) | 2.18 [2.00, 2.38] | 5 [3, 13] | 0.30 |
| MutationAssessor | 4.65 (33) | 2.16 [1.84, 2.47] | 6 [3, 17] | 0.50 |
| Eigen | 5.85 (23) | 2.16 [1.93, 2.40] | 7 [3, 14] | 0.44 |
| MisFit_D | 7.53 (13) | 2.12 [1.87, 2.39] | 8 [4, 16] | 0.37 |
| PROVEAN | 3.12 (37) | 2.12 [1.77, 2.50] | 9 [3, 22] | 0.66 |
| VARITY_ER_LOO | 5.56 (26) | 2.10 [1.86, 2.34] | 10 [4, 15] | 0.43 |
| MisFit_S | 16.25 (1) | 1.90 [1.77, 2.05] | 17 [11, 21] | 0.23 |
| MutScore | 14.85 (2) | 1.77 [1.65, 1.92] | 21 [16, 27] | 0.21 |
| PHACTboost | 14.76 (3) | 1.68 [1.59, 1.79] | 27 [21, 31] | 0.19 |
| MPC | 13.90 (4) | 1.31 [1.21, 1.42] | 39 [35, 39] | 0.10 |

**Where the rankings differ, and where they do not.**
- Across all genes, the two rankings were not positively associated at any threshold or site set: Spearman ρ between −0.02 and −0.22. At the 2% threshold it was −0.17 [−0.26, 0.06].
- The four strongest pooled predictors ranked 17th to 39th within windows. MPC, built from regional constraint, ranked last.
- Only VEST4's within-window LR had a lower bound above the ClinGen Supporting LR (2.406). Below the top two, rank intervals span 5–19 places.

**Within DDG2P genes the rankings agreed** (ρ 0.68 [0.53, 0.76]), with the same five predictors leading. The divergence across all genes therefore comes from between-gene information.

**Reproducibility.** The within-window ranking replicated in autism (ρ 0.85 across 25 predictors) and, more weakly, in the public cohort (0.49).

**Within-window evidence tracks independent variant-level performance, but not demonstrably better than pooled evidence.**
- Across the 39 predictors, within-window LR correlated with within-gene ClinVar AUROC (ρ 0.55 [0.22, 0.77]; 0.56 in DDG2P genes) and with assay performance (ρ 0.63 [0.39, 0.79]).
- Marginal LR also correlated with them (ρ 0.37 for each).
- The differences were not significant: +0.18 [−0.32, 0.63] for ClinVar and +0.26 [−0.11, 0.65] for assays.

Where PP3 is applied (within disease genes), pooled and within-window rankings largely agree.

### PP3, PM1 and PP2 are not independent
**PP3 and PM1** (Figure 3). For AlphaMissense PP3 at Supporting with a PM1-like flag, R was:

| Comparison | DD | Autism | Public |
|---|---|---|---|
| Marginal | 0.44 [0.41, 0.47] | 0.56 [0.47, 0.66] | 0.67 [0.48, 0.96] |
| Within gene, DDG2P monoallelic genes | 0.69 [0.65, 0.73] | 0.74 [0.67, 0.82] | 0.75 [0.62, 0.89] |
| Within gene, all genes | 0.75 [0.71, 0.78] | 0.96 [0.84, 1.10] | 0.96 [0.74, 1.28] |

- Part of the dependence is between genes. Across all genes in autism and the public cohort, the within-gene dependence was indistinguishable from independence.
- Within DDG2P disease genes, where PP3 is applied clinically, a within-gene component remained in all three cohorts.
- Dependence was stronger with PM1-strict (DD 0.37 [0.35, 0.40]) and weaker for REVEL (0.57) and VEST4 (0.67) (Supplementary Table 5).
- Across the 40 predictors informative on all analysed sites, R × PM1-like had an upper bound below 1 in DD (median 0.51). This set includes popEVE, phyloP470way and fathmm-MKL; the benchmark's site set (constraint annotated) instead admits SIFT and PolyPhen2_HVAR, and 36 predictors are in both. Predictors using population constraint or allele frequency were no more dependent than others (p = 0.87).

**PP3 and PP2, PM1 and PP2.**
- PP3 was less dependent on the gene-level flag (median R across predictors 0.70 in DD, 0.83 in autism).
- PM1 and PP2 were the most dependent pair: R 0.31 [0.28, 0.33] in DD and 0.36 [0.30, 0.45] in autism for PM1-like, and 0.27 and 0.34 for PM1-strict.
- Part of this is by construction. Both of our flags are built from missense constraint, which ACGS assigns to PP2 rather than PM1 (Durkie et al. 2024). Panel-defined PM1 rests also on hotspots and functional domains, so its dependence on PP2 may be smaller.

**Where the dependence comes from.** It arises mainly on the null side. AlphaMissense flags 51% of mutation-expected sites at PM1-like positions, against 19% overall. Among excess de novos the rise is only from 79% to 94%.

**Robustness.** The dependence held:
- in ordered (pLDDT ≥ 70) and disordered regions;
- at buried and exposed sites;
- after excluding the 50 genes with the largest excess (R 0.48 [0.42, 0.53]) (Supplementary Note 5).

**Strong tier, against PM1-strict.** f for AlphaMissense was 0.72 [0.70, 0.74] in DD, 0.75 [0.69, 0.80] in autism and 0.70 [0.57, 0.81] in the public cohort. For REVEL it was 0.83 and for VEST4 0.87 in DD. On the de novo scale, the dependence removes 25–30% of AlphaMissense's summed log-evidence, close to the cap's 33%, and less for REVEL and VEST4.

### Independent truth
**Table 4. Dependence of AlphaMissense PP3 and PM1 across truth sets.**

| Truth set | PM1 flag | R, Supporting | R, Strong |
|---|---|---|---|
| De novo excess, DD | PM1-like | 0.44 [0.41, 0.47] | 0.33 [0.31, 0.35] |
| De novo excess, DD, within DDG2P genes | PM1-like | 0.69 [0.65, 0.73] | 0.65 [0.60, 0.69] |
| ClinVar, DDG2P genes | PM1-like | 0.32 [0.28, 0.37] | 0.23 [0.18, 0.28] |
| Expert panel (2,516 variants, 108 genes) | Panel-applied PM1 | 0.47 [0.22, 1.18] | 0.17 [0.08, 0.61] |
| Expert panel, frequency criteria also removed (1,540) | Panel-applied PM1 | 0.34 [0.17, 1.55] | 0.12 [0.05, unbounded] |
| Assays: PPM1D / SHOC2 | PM1-like | 0.72 [0.65, 0.78] / 0.86 [0.71, 1.04] | 0.71 / 0.82 |
| Reclassified ClinVar (87 P, 887 B) | PM1-like | Not estimable (3 B/LB in the joint cell) | Not estimable |
| Expert panel, panel-applied PP3 × PM1 (2,942 missense, 130 genes) | Panel-applied PM1 | 0.21 [0.12, 0.31] (any PP3 strength) | — |

**Every estimate is below 1, but the size varies.**
- The dependence is larger on ClinVar and expert-panel truth than on the de novo scale.
- It is smallest in the two assays: the most independent truth, but only two genes.
- ClinVar estimates are biased towards 1, because curators applying both criteria enrich P/LP variants in the joint cell. They may also be biased away from 1 through PM1-region ascertainment.
- At Supporting, the expert-panel interval includes independence.

**How many points PP3 earns beside the panel's own PM1.**
- At Strong: 1.83 points [0.84, 3.54]. This excludes the 4 points of uncapped co-use, and includes both the cap's allowance (2) and the conditional points derived from de novos (1.1).
- At Supporting: 1.68 [1.00, 2.65].

### Stacking in expert practice
**How panels combine the criteria.** Among 6,177 missense classifications in the ClinGen Evidence Repository (46 panels):

| Criteria applied | Classifications |
|---|---|
| PP3 | 2,990 |
| PM1 | 1,092 |
| PP2 | 592 |
| PP3 + PM1 | 787 |
| PM1 + PP2 | 174 |
| All three | 138 |

- **The cap rarely applies.** Where PP3 accompanied PM1, it was applied at Supporting in 605 classifications, at Moderate in 178 and at Strong in 4. The ClinGen cap binds only at Strong, so it constrains almost none of this stacking.
- **PM1 + PP2 is concentrated in seven panels:** RASopathy 62, FBN1 61, Monogenic Diabetes 22, Brain Malformations 11, PTEN 11, Pulmonary Hypertension 5 and ACADVL 2.
- **Some classifications predate current specifications.** For example, the RASopathy panels updated theirs in 2025 (ClinGen RASopathy Expert Panels 2025). The repository records each classification as approved under the specification then in force.

**Dependence of the panels' own criteria on independent truth.**
- **The truth set.** Classifications were re-scored after removing PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5: 2,942 missense classifications in 130 genes reached P/LP (≥ 6 points) or B/LB (≤ −1). Unlike the Table 4 analysis, this uses the panels' own PP3 rather than AlphaMissense, so no site mapping is needed.
- **PP3 × PM1.** Panel-applied PP3 and PM1 were strongly dependent: R 0.21 [0.12, 0.31] (gene-clustered).
- **What PP3 earns beside PM1.** Once a panel applied PM1, its PP3 earned 0.60 points [0.14, 1.05] of independent evidence, against a mean 1.24 points applied. By strength:
  - **Supporting PP3** (235 classifications) earned 0.81 points [0.24, 1.68] against 1.00 applied.
  - **Moderate PP3** (67 classifications; 6 B/LB in the joint cell) earned 0.02 points [−1.69, 0.89] against 2.00 applied. The interval excludes the points applied, though it is wide. The ClinGen cap allows Moderate PP3 with Moderate PM1 (4 points), so it never intervenes here.
- **Robustness.** The credit was stable across truth thresholds: 0.60–0.65 overall, 0.81–0.89 for Supporting PP3 and −0.23 to 0.14 for Moderate. With frequency criteria also removed it was 0.49 [−0.06, 1.31], on only 6 B/LB variants in the joint cell (Supplementary Note 11).
- **Truth that conditions on other evidence.** Truth requires the remaining evidence to reach ≥ 6 or ≤ −1 points, and panels choose PP3 strength knowing that evidence. Both can bias the credit; neither the threshold analysis nor the out-of-sample tests suggest a large effect.
- **Out of sample.**
  - **Held-out genes.** Across 200 random splits of the genes, the credit estimated on one half was 0.61 (median) and on the held-out half 0.59. The held-out credit fell below the points applied in 98% of splits.
  - **Later classifications.** Estimated on classifications approved before 2025, the credit was 0.70 [−0.14, 1.26]. On classifications approved in 2025–2026 it was 0.40 [−0.25, 1.07], against 1.40 points applied.
- **PM1 × PP2 and PP3 × PP2** could not be estimated on this truth: there were 1 and 2 B/LB variants in the joint cells.

**The triple stack on the de novo scale.** When AlphaMissense PP3, PM1-like and PP2-like all apply:

- **Dependence.** Their joint de novo evidence was far below the product of their separate de novo LRs (R3 0.13 [0.11, 0.14] at Supporting, 0.09 at Strong; autism 0.18 and 0.14).
- **Against the nominal ACMG weights** (PP3 1 or 4, PM1 2, PP2 1):

| PP3 tier | Nominal sum | Under the cap | De novo evidence, DD | De novo evidence, autism |
|---|---|---|---|---|
| Supporting | 4 | 4 | 3.65 [3.50, 3.79] | 2.96 [2.67, 3.22] |
| Strong | 7 | 5 | 4.05 [3.89, 4.20] | 3.27 [2.96, 3.53] |

At Supporting, the nominal weights approximately match the joint evidence. With Strong PP3 they exceed it by about 3 points uncapped and 1 point capped. The de novo scale is not the clinical scale, so these are comparisons of direction and size, not calibrated points.

**Effect on expert classifications.** Re-scoring the 6,177 missense classifications under alternative rules changed:

| Rule | Classifications changed | P → LP | LP → VUS |
|---|---|---|---|
| ClinGen cap (PP3 + PM1 ≤ 4) | 24 | 12 | 12 |
| Cap on all three location-type criteria (PP3 + PM1 + PP2 ≤ 4) | 35 | 17 | 18 |
| PP2 not added when PM1 applies | 39 | 28 | 11 |
| PP3 limited to Supporting when PM1 applies | 43 | 22 | 21 |
| PP3 credited at its measured value beside PM1 (0.60 points; in brackets, at the upper 95% bound) | 187 (43) | 84 (22) | 102 (21) |

**Out of sample.**
- **Held-out genes.** With the credit estimated on the other half of the genes, the measured-credit rule changed a median of 89 of 3,103 classifications per split (2.9%; 95% range 13–140), equivalent to about 177 of 6,177.
- **Later classifications.** Among the 2,507 missense classifications approved in 2025–2026, limiting PP3 to Supporting would change 15 and the cap 16.

Limiting PP3 to Supporting needs no estimated parameter.

**Replication in laboratory practice.**
- **Data.** We repeated the analysis on ClinVar 2026-10 germline submissions from diagnostic laboratories, excluding expert-panel and practice-guideline submissions. We kept submissions whose comments list ACMG/AMP codes for missense SNVs: 87,873 submissions from 472 submitters.
- **How labs stack the criteria.** Laboratories applied PP3 beside PM1 in 4,477 submissions, almost always below Strong: Supporting 4,212, Moderate 185, Strong 80.
- **Truth.** Points from the remaining codes, with the same exclusions and thresholds as for the panels: 39,075 submissions in 5,976 genes.
- **Dependence.** PP3 and PM1 were dependent, R 0.18 [0.14, 0.25], close to the panels' 0.21.
- **Supporting PP3** beside PM1 earned 0.83 points [0.53, 1.24] against 1.00 applied, matching the panels' 0.81.
- **Moderate PP3** earned 0.75 points against 2.00 applied, but with only 3 B/LB variants in the joint cell its interval is unbounded. Strong PP3 had none.

Laboratory data therefore replicate the dependence and the calibration of Supporting PP3, but cannot independently confirm the result for Moderate PP3.

### Clinical consequences
**DD de novo classification** (PS2 + PM2 + PM1-like + AlphaMissense PP3):
- Uncapped co-use placed 2,438 variants at Pathogenic. The de novo excess implies that 84% of that class is disease-associated.
- Taken literally, 84% falls below the 90% floor for Likely Pathogenic as well as the 99% floor for Pathogenic. The estimate is class-level, and refers to association with any severe developmental disorder rather than a per-variant posterior, so we do not use it to choose between rules. It shows that the uncapped sum is over-confident for this class.

**Expert-panel classifications.** See *Stacking in expert practice*. The cap changes 24 classifications whether all 13,293 classifications or only the 6,177 missense ones are considered. Of the 72 classifications in which PP3 + PM1 exceed 4 points, 70 are missense, and the 2 that are not do not change class.

**ClinVar.** In DDG2P monoallelic genes, 37.3% of AlphaMissense PP3-positive VUS lie at PM1-like sites (13.2% across all genes). Without the cap, PP3 + PM1 at Moderate alone would reach Likely Pathogenic for 7,101 of them.

**Conditional points.** Supplementary Table 4 gives the PP3 points that remain once PM1 or PP2 applies, by predictor and tier.

### Attempts to build a better rule from de novo data
Two pre-specified attempts failed (Supplementary Notes 6 and 7).

- **A label-free gene → window → AlphaMissense model** calibrated worse than the cap on reclassified ClinVar variants (log-loss difference +0.147 [+0.106, +0.190]). Most of those variants lay outside NDD genes. Even within NDD genes, its gene-specific PP3 evidence did not track ClinVar's (ρ 0.01).
- **A model trained only on within-window variation** added nothing to AlphaMissense on held-out de novos (+0.0008 [−0.0006, +0.0023] nats per de novo). It was worse on ClinVar (AUROC 0.894 vs 0.967) and in 7 of 8 assays.

Earlier de novo-trained classifiers behaved the same way (Supplementary Note 4). The de novo signal within windows is enough to rank predictors but not to train a better one.

---

## Discussion

**What we measured.** We measured the evidence missense predictors add beyond the gene and the 50-residue window a variant falls in, without variant-level labels:
- **About a third of the log-evidence is incremental to gene and window.** That is the quantity that matters where location is credited separately, as ACMG/AMP does through PP2 and PM1.
- **It is not a statement about skill.** A perfect predictor would also lose evidence within windows, by an amount that depends on an enrichment parameter the data do not identify.

**Benchmarks.** Pooled de novo comparisons favour predictors that encode gene and regional constraint.
- **Across all genes the rankings differ**, consistent with Harrizi et al. (2026) on ClinVar. The predictors leading pooled de novo comparisons tend to encode constraint; MPC, built from regional constraint, falls to last.
- **Within known disease genes the rankings agree**, and within-window evidence tracks within-gene clinical and functional performance no better than pooled evidence does, as far as these data can tell.
- **We therefore offer the benchmark as a mechanism and a resource**, not as evidence that a different predictor should be used for PP3.

**Independence of criteria.** No pair of PP3, PM1 and PP2 was independent on the de novo scale.
- **PM1 × PP2** was the most dependent pair. Applying both counts gene-level and regional constraint largely twice. UK guidance already advises against counting the same constraint evidence in both (Durkie et al. 2024). Our data quantify the dependence for constraint-based flags, and show that seven ClinGen panels co-apply the two.
- **PP3 × PM1** was dependent in three cohorts, within disease genes and on clinical truth, though with a size that varies across truth sets.
- **Expert practice.** Panels stack these criteria, mostly at strengths the existing cap does not reach.
  - On the panels' own classifications, Supporting PP3 beside PM1 was approximately calibrated, while Moderate PP3 added no detectable evidence.
  - On the de novo scale the nominal weights match the joint evidence when PP3 is Supporting, and exceed it when PP3 is stronger.
  - The two scales agree that the overlap that matters lies above Supporting PP3.
- **What it implies for rules.** On expert-panel truth the evidence PP3 keeps beside PM1 excludes uncapped co-use at the Strong tier. It cannot distinguish ClinGen's cap from a conditional rule. On the de novo scale the cap removes slightly more than the measured dependence. We therefore support the cap, or a conditional rule, over uncapped co-use, and do not claim the cap's size is optimal.
- **Recommendations.**
  - **Limit PP3 to Supporting when PM1 applies.** This is consistent with the data and needs no estimated parameter. Supporting PP3 beside PM1 was calibrated in both expert-panel and laboratory practice. The evidence against Moderate PP3 rests on 67 panel classifications, which laboratory data point the same way but cannot confirm, so we suggest panels check it on their own data.
  - **Check PM1 + PP2 for shared evidence**, as ACGS advises. Our flags quantify a strong dependence, partly by construction, and expert truth was too sparse in benign variants to test panel-defined PM1 + PP2.
  - **Calibrate co-application jointly.** Panels that keep co-applying location-type criteria should calibrate the combination against independent truth, as done here.

**Relation to previous work.**
- Chao et al. (2024) showed that regions with less than 20% of expected missense variation reach Moderate support, and that MPC, which combines regional constraint with variant-level information, separates DD cases from controls. Our results agree that combining the two adds evidence: the joint LR exceeds each marginal LR (15.0 against 4.2 and 8.1). What fails is the product assumption.
- Gene-level heterogeneity and circularity of predictor performance are documented on ClinVar (Grimm et al. 2015; Tejura et al. 2024; Harrizi et al. 2026). They have motivated gene-level calibration (Badonyi & Marsh 2025; Chen et al. 2026). To our knowledge, no previous evaluation has measured predictors within gene and window on label-free de novo data, or tested the independence of PP3, PM1 and PP2 directly.

**Limitations.**
1. **Pre-specification.** The primary criterion failed because its position control did not control position in genes without a constraint split, and the replacement control was adopted afterwards. The amendments were public but not independently registered. We present continuous estimates rather than the binary criterion.
2. **Estimand.** The de novo LR reflects high-penetrance variants in severe early-onset disorders, not pathogenicity for a specific condition. The oracle depends on an unknown enrichment k.
3. **The PM1 and PP2 flags are proxies.** Panel-specific PM1 regions differ; expert-panel truth used the panel's own PM1.
4. **Truth sets.** Laboratory truth is self-reported in free-text comments, and its codes were parsed automatically. Expert-panel truth covers 108–130 genes, is ascertained through panels, and was used both to estimate PP3's credit beside PM1 and to show the effect of applying it. Gene-split cross-validation reproduced the credit and the size of the effect in held-out genes. The panels are, however, a single source. Assay truth tests dependence in two genes. Reclassified ClinVar had too few benign variants in the joint cell.
5. **Predictors are correlated.** Shares and rank correlations across predictors are less certain than their nominal counts suggest; we resampled predictors where it matters.
6. **Transcripts.** 74.4% of GENCODE v32 canonical transcripts in MANE genes match MANE Select v1.5 (66.5% in DDG2P genes). Window assignments may differ for the rest.
7. **Ancestry.** The cohorts are predominantly of European ancestry.

**Conclusion.** Most of what missense predictors contribute concerns where a variant lies, and the criteria that credit location separately are not independent. In expert practice, PP3 is stacked beside PM1 mainly at strengths the existing cap does not reach. There, Supporting PP3 appears calibrated and Moderate PP3 does not. Limiting PP3 to Supporting when PM1 applies is a simple rule consistent with these data.

---

## Figures and tables

- **Figure 1.** Design: de novo excess, conditioning levels and retention (`fig1_design.png`).
- **Figure 2.** Pooled versus within-window predictor ranks (`fig2_benchmark.png`).
- **Figure 3.** Dependence of PP3 and PM1/PP2 flags, R with 95% intervals by predictor, flag and cohort (`fig4_overlap.png`).
- **Tables:**
  - Table 1: pre-specified tests;
  - Table 2: oracle;
  - Table 3: benchmark;
  - Table 4: dependence across truth sets.

**Supplementary Notes:**
1. Pre-specification and deviations.
2. Robustness of retention.
3. Estimator simulation.
4. De novo-trained classifiers.
5. Full dependence results, including structural strata and gene-clustered intervals.
6. Hierarchical model.
7. Within-window model.
8. Oracle derivation.
9. Benchmark robustness and validation.
10. PM1 × PP2 dependence and the triple stack.
11. Stacking in expert practice.

**Supplementary Tables:**
1. Per-predictor retention.
2. Training-data ledger.
3. ClinVar sensitivity.
4. Conditional PP3 points.
5. Full dependence test.
6. Hierarchical-model calibration table.
7. Benchmark, all predictors, thresholds and site sets.

## Data and code availability
All analyses use public data: published de novo tables, ClinVar, the ClinGen Evidence Repository, ProteinGym, MaveDB, dbNSFP (academic licence), AlphaMissense, AlphaFold DB, gnomAD constraint and Roulette mutation rates. Code, the pre-specification history and all outputs are in the study repository [URL; Zenodo DOI to be minted at submission].

## Competing interests
[To be completed by the authors.]

## Funding
[To be completed by the authors.]

---

## References

- Badonyi M, Marsh JA. acmgscaler: standardized gene-level variant effect score calibration within the ACMG/AMP framework. 2025.
- Bergquist T, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 2025.
- Chao KR, et al. The landscape of regional missense mutational intolerance quantified from 125,748 exomes. bioRxiv 2024. doi:10.1101/2024.04.11.588920.
- Chen et al. Gene- and domain-aware calibration increases the clinical utility of variant effect predictors. Research Square preprint, 2026.
- Cheng J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 2023.
- Durkie M, Cassidy E-J, Berry I, et al. ACGS best practice guidelines for variant classification in rare disease 2024. Association for Clinical Genomic Science, 2024.
- Gao H, et al. The landscape of tolerated genetic variation in humans and primates (PrimateAI-3D). *Science* 2023.
- Garrett A, et al. Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. *J Med Genet* 2021;58:297.
- Grimm DG, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat* 2015.
- Harrizi et al. Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. bioRxiv 2026. doi:10.64898/2026.08.21.746180.
- Havrilla JM, et al. A map of constrained coding regions in the human genome. *Nat Genet* 2019.
- Jin SC, et al. Contribution of rare inherited and de novo variants in 2,871 congenital heart disease probands. *Nat Genet* 2017;49:1593–1601.
- Kaplanis J, et al. Evidence for 28 genetic disorders discovered by combining healthcare and research data. *Nature* 2020.
- Liu X, et al. dbNSFP v4. *Genome Med* 2020;12:103.
- Notin P, et al. ProteinGym: large-scale benchmarks for protein fitness prediction and design. *NeurIPS* 2023.
- Orenbuch R, et al. Proteome-wide model for human disease genetics (popEVE). *Nat Genet* 2025.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet* 2022.
- Qi H, et al. MVP predicts the pathogenicity of missense variants by deep learning. *Nat Commun* 2021;12:510.
- Richards S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med* 2015.
- Seplyarskiy V, et al. A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription (Roulette). *Nat Genet* 2023.
- Stenton SL, et al. Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med* 2024.
- Tavtigian SV, et al. *Genet Med* 2018; Tavtigian SV, et al. *Hum Mutat* 2020.
- Tejura M, et al. Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes. *Am J Hum Genet* 2024.
- Turner TN, et al. denovo-db. *Nucleic Acids Res* 2017.
- ClinGen RASopathy Expert Panels. Updated ACMG/AMP specifications for variant interpretation and gene curations from the ClinGen RASopathy expert panels. *Genet Med Open* 2025.
- Zhou X, et al. Integrating de novo and inherited variants in 42,607 autism cases identifies mutations in new moderate-risk genes. *Nat Genet* 2022.
