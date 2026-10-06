# Computational (PP3) and hotspot/domain (PM1) evidence overlap for de novo missense variants in neurodevelopmental disorders: a label-free test of additivity

*Draft 3, 6 October 2026. Restructured in response to four external reviews (point-by-point responses in `REVIEW_RESPONSE.md`). All values come from committed outputs.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** ACMG/AMP classification adds computational evidence (PP3) to hotspot/domain evidence (PM1) and gene-level evidence (PP2) as if they were independent. ClinGen caps PP3 + PM1 at Strong; other guidance permits uncapped co-use. We tested, without variant-level clinical labels, whether these criteria are additive.

**Methods.** Likelihood ratios (LRs) were estimated from the excess of de novo missense variants over a mutation-rate expectation in three cohorts: 31,058 developmental-disorder (DD) trios, 12,316 autism probands and 4,256 public trios with other phenotypes. Additivity was measured as R, the joint LR of PP3 and a PM1-like or PP2-like flag divided by the product of their separate LRs (R = 1: additive). R was decomposed into between-gene and within-gene components, and we measured how much predictor evidence survives conditioning on gene and protein position. A planted-signal simulation validated the estimator.

**Results.** PP3 and PM1-like evidence were not additive. For AlphaMissense, R was 0.44 [95% CI 0.41–0.48] in DD, 0.56 [0.44–0.67] in autism and 0.67 [0.42–0.95] in the public cohort. Across REVEL and VEST4 and all flags, R ranged from 0.48 to 1.02, with overlap weakest for PP2-like flags outside DD. Part of the overlap is between genes. Comparing variants within the same gene across all genes, R was 0.75 [0.71–0.78] in DD and 0.96 [0.84–1.08] in autism, which is consistent with additivity; within DDG2P disease genes it was 0.69 [0.68–0.70] and 0.74 [0.68–0.80]. At the Strong PP3 tier, overlap removed 6–27% of the summed log-evidence, against the 33% removed by ClinGen's cap. In a classification analysis, uncapped combination with de novo evidence called 2,438 DD variants Pathogenic, a class in which only 84% were disease-associated according to the de novo excess; capping moved all of them to Likely Pathogenic. Predictors retained a median 0.33 of their log-evidence within 50-residue windows, against 0.51 under a same-size random partition.

**Conclusion.** For de novo missense variants, PP3 and PM1-like evidence overlap substantially, including within known disease genes. The overlap is consistent with ClinGen's cap, which is large enough to absorb it, and not with uncapped co-use.

---

## Introduction

**The setting.** Missense variants are the largest class of variants of uncertain significance. Their classification under ACMG/AMP combines several kinds of evidence:

- computational predictions (PP3);
- location in a mutational hotspot or critical functional domain (PM1);
- membership of a gene with a low rate of benign missense variation (PP2).

**How the pieces are combined.** In the Bayesian points framework the evidence is additive in log-likelihood-ratio units, with a single prior probability for every gene. The ClinGen Sequence Variant Interpretation working group calibrated predictor thresholds to LRs using ClinVar (Pejaver et al. 2022; Bergquist et al. 2025) and, recognising that PP3 and PM1 can describe the same property, capped their combination at Strong (Stenton et al. 2024). Other guidance differs: CanVIG-UK permits PM1 and PP3 together, regarding regional enrichment and in silico prediction as largely different evidence types (Garrett et al. 2021). Genome-wide calibration also masks heterogeneous performance across genes (Tejura et al. 2024), which has motivated gene-level calibration (Badonyi & Marsh 2025; Chen et al. 2026).

**What is not known: how much the evidence overlaps.** A predictor that has learned which genes are intolerant, or which domains are hotspots, separates pathogenic from benign variants pooled across genes. Yet it may add little once the gene and the domain are known. ClinVar-based calibration cannot easily measure this:

- its labels are concentrated in a few genes;
- several predictors were trained on them;
- the evaluation is dominated by gene identity (Harrizi et al. 2026).

**Our approach.** Here we use de novo variants, whose excess over mutation-rate expectation in affected children defines a likelihood ratio without variant-level labels, and which can be estimated *within* a gene, window or domain. We ask four questions:

1. How much of each predictor's evidence survives conditioning on gene, protein position and regional constraint?
2. Is the loss due to position, or to the act of conditioning?
3. Does it replicate across cohorts, phenotypes and predictor generations?
4. What does it imply for combining PP3 with PM1 and PP2?


---

## Methods

### Study design and governance
**Pre-specification.** The analysis plan was pre-specified before the full analysis and extended by eleven amendments, each written before the analyses it governs. Each document was committed to a version-controlled repository and timestamped with OpenTimestamps.

- **Not an independent registry.** The documents were not deposited in an independent registry; the timestamps establish when the text was fixed, not that no results were seen. Each amendment states what had been seen when it was written.
- **The primary pre-specified test in the DD cohort did not meet its criterion** (Results, Q3). The position control used as primary in later tests was adopted after that result, for a reason stated in advance of those tests.
- **Post hoc analyses are labelled as such.** These include two of the four models and the threshold sensitivity.

The amendment timeline, the deviation log and the commit identifiers are in Supplementary Note 1.

### Cohorts
**DD cohort.** 31,058 parent–offspring trios with severe developmental disorders (Kaplanis et al. 2020), from three centres: DDD, GeneDx and Radboudumc.

**Autism cohort.** De novo SNVs in 12,316 autism probands, with 4,051 unaffected siblings as a negative control (Zhou et al. 2022; Supplementary Data 1–2).

- The supplementary tables list only individuals who carry a coding de novo: 13,103 affected, of whom 12,501 carry an SNV.
- Individuals sharing ≥ 2 identical de novo SNVs with a DD proband were removed (n = 64).

**Public multi-phenotype cohort.** 4,256 trios assembled from public sources:

| Source | Phenotype | Trios |
|---|---|---|
| Jin et al. 2017 | Congenital heart disease | 2,645 |
| denovo-db v1.6.1 | Schizophrenia | 785 |
| denovo-db v1.6.1 | Tourette syndrome | 511 |
| denovo-db v1.6.1 | Epileptic encephalopathy | 264 |
| denovo-db v1.6.1 | Intellectual disability | 51 |

- Excluded: studies overlapping the DD or autism cohorts, and studies reporting fewer than 0.15 synonymous de novos per proband (incomplete calls).
- Eight individuals were removed as duplicates of DD probands.

### Sites, annotation and versions
**Universe.**

- Coding sites on canonical GENCODE v32 transcripts: 18,397 transcripts, of which 17,788 passed quality control. Concordance with MANE Select was not assessed.
- Sites within 3 bp of an exon boundary were excluded on both the observed and the expected side.
- After filtering, the DD cohort contributed 23,851 missense and 8,222 synonymous de novos.

**Data versions:**

- **Predictor scores:** dbNSFP 4.1a (36 tools) and 5.4a (17 tools plus popEVE).
- **AlphaMissense:** the published hg38 release.
- **Regional missense constraint:** gnomAD v2.1.1, lifted to GRCh38 (Chao et al. 2024). Chao et al. also calibrated regional missense constraint to ACMG/AMP evidence strengths and showed that MPC stratifies case and control de novo missense variants; our analysis differs in conditioning predictor evidence on such regions rather than calibrating the regions themselves.
- **Constrained coding regions:** Havrilla et al. 2019.
- **Gene constraint:** gnomAD v2.1.1.
- **InterPro domains:** dbNSFP 5.4a.
- **AlphaFold DB:** human proteome v6.
- **ClinVar:** release 2026-10-04.
- **Mapping:** of 2,276,459 qualifying ClinVar missense SNVs, 1,987,084 (87.3%) mapped to analysed sites; the rest fell on non-canonical transcripts or filtered sites.

### Likelihood ratios from the de novo excess
**Expected counts.** Each site's mutation rate is from the Roulette model (Seplyarskiy et al. 2023). A cohort-specific scale factor, estimated from synonymous de novos, converts rates to expected counts E. The excess is X = O − E.

**Marginal LR** for a set of sites *b*:

LR = (X_b / ΣX) / (E_b / ΣE)

Because it is built on the excess, the LR is not diluted by the background de novos.

**Within-unit LR.** Each observed de novo is compared with the opportunity profile of its own unit (gene, window, domain):

LR_b = X_b / (Σ_u O_u · E_ub / E_u − E_b)

**The estimand.** This is the evidence a score provides for a de novo missense variant in an affected individual relative to mutational background. It is not the ClinVar pathogenic-versus-benign LR. Its comparison class includes damaging variants in genes unrelated to the phenotype, and it reflects variants of high penetrance in probands. Retention shares, which compare a predictor with itself, are therefore the primary quantity; absolute comparisons with ClinGen cut-offs are secondary.

**Uncertainty.**

- 1,000 Poisson bootstrap replicates over individuals (seed 20261004).
- Replicates share individual weights across analyses of the same cohort, which gives paired intervals.
- A pre-specified margin δ = 0.091, the largest calibration error of the mutation-rate model, defines one-sided lower bounds.

### Conditioning levels and retention
**Levels:**

| Level | Held fixed |
|---|---|
| A | Nothing (marginal) |
| B | Gene |
| C | Gene × gnomAD regional missense constraint region |
| D1 | Gene × evenly spaced chunks, as many as the gene's constraint regions |
| D2 | Gene × fixed 50-residue window |
| E1 | D1 chunk × constraint region |
| E2 | Window × constraint region |

**Two position controls.** The original plan specified D1/E1 as the position control.

- **D1's limitation:** in the ~13,500 genes without a regional-constraint split, D1 has a single unit. It therefore does not control position in those genes.
- **D2/E2:** this limitation was recorded after the DD analysis and before the autism and newer-predictor analyses, and D2/E2 (fixed windows) was pre-specified as primary for those.
- **Both are reported throughout.**

**Retention.** For each predictor, r(level) = ln LR(level) / ln LR(A), computed for the top 2% of the predictor's scores by mutational opportunity. Points in the ACMG framework are proportional to ln LR, so r is the share of a predictor's points that survives conditioning.

**Informative predictors** have a marginal top-2% LR ≥ 2.406.

**Pre-specified readout.** The share of informative predictors with r(D2) < 0.5, called "systematic" if ≥ 75%.

**Controls:**

- a random partition of each gene into the same number of non-contiguous groups as D2;
- windows shifted by 25 residues;
- gene × InterPro domain units;
- unaffected siblings.

### Direct test of PP3 + PM1/PP2 overlap
**Flags at each site:**

- **PP3:** the published ClinGen Supporting threshold for AlphaMissense, REVEL and VEST4.
- **PM1-like:** regional missense o/e < 0.4 or CCR ≥ 95.
- **PM1-strict:** o/e < 0.2.
- **PP2-like:** gene missense Z ≥ 3.09.

**Test.** For each predictor and flag F:

R = LR(PP3 ∧ F) / [LR(PP3) · LR(F)] = LR(PP3 | F) / LR(PP3)

- R = 1 means the evidence is additive, as the points framework assumes.
- R < 1 means applying both criteria overstates the joint evidence.
- The PP3 points implied when F is applied are ln LR(PP3 | F) / ln 2.406.

**Analysis.** DD is the primary cohort; autism is the replication; the public cohort is descriptive.

**Within-gene and between-gene overlap** (A13). R was re-estimated with the within-unit estimator, comparing variants only within the same gene. This was done across all genes and within DDG2P monoallelic genes. The difference from the marginal R is the between-gene component.

**Strong tier.** The test was repeated with PP3 at each predictor's published Strong threshold. Two summaries are reported:

- **Relative (primary):** the fraction of the summed log-evidence kept in the joint evidence, f = ln LR(PP3 ∧ F) / [ln LR(PP3) + ln LR(F)].
  - ClinGen's cap corresponds to f = 0.67: Strong PP3 (4 points) plus Moderate PM1 (2 points), capped at 4.
  - This comparison does not depend on the scale of the LR.
  - It was derived after the A10 outputs were committed and is labelled exploratory.
- **Absolute (secondary):** joint points, ln LR(PP3 ∧ F) / ln 2.406. The de novo-excess LR measures association with high-penetrance disease across all genes, not pathogenicity for a specific disorder. Absolute points are therefore not directly comparable with ClinGen's scale.

### Clinical-impact analysis
**Points.** ClinVar missense variants were scored with published PP3 points, with PM1 at +1, +2 or +4 at PM1-like sites, under three rules:

1. the ClinGen cap;
2. a conditional rule (published PP3 points × the predictor's retention);
3. this study's within-unit thresholds.

**Scope.** Retention was measured in dominant de novo DD/NDD variants. Applying it to ClinVar variants of other inheritance modes and phenotypes is an extrapolation, so results are also reported for DDG2P monoallelic genes only.

---

## Results

### Q1. How much predictor evidence survives conditioning on gene and position?
Of the 36 dbNSFP 4.1a predictors, 24 were informative in DD. Median retention:

| Level | Median r |
|---|---|
| B: within gene | 0.56 |
| C: within constraint region | 0.38 |
| D2: within 50-residue window | 0.33 |
| E2: window × constraint | 0.32 |

(Figure 2, Figure 3, Supplementary Table 1.)

**AlphaMissense** (top 2%, common site set): LR 13.4 marginal, 4.10 within gene, 2.19 within window, 2.08 after constraint. Its retention was r(B) 0.54 [0.53, 0.56] and r(D2) 0.30 [0.28, 0.32].

**Constraint adds little once position is fixed:** median r(D2) − r(E2) was 0.01–0.02.

**The strongest marginal predictors retained the least.** Examples:

| Predictor | Marginal LR | Within gene |
|---|---|---|
| MPC (built from regional constraint) | 13.6 | 1.95 |
| popEVE (DDD only) | 17.9 | 3.13 |

popEVE ended at 1.69 after window and constraint.

**Threshold sensitivity.** Retention was stable across score thresholds: median r(D2) 0.33, 0.32 and 0.28 at the top 1%, 2% and 10% (exploratory).

### Q2. Is the loss due to position or to conditioning on small units?
**Random partitions retained nearly all the within-gene evidence.** This control used all analysed sites; the main analyses used the common site set, which is why its values differ slightly from Q1 (0.31 against 0.33). On this site set, a random partition of each gene's residues into the same number of groups as the 50-residue windows retained a median 0.51, against 0.56 for gene alone (both all sites). Real windows retained 0.31.

- The paired difference was 0.21 [0.19, 0.23], positive for all 24 predictors.
- Shifting window boundaries by 25 residues changed nothing (−0.02).
- Curated InterPro domains retained 0.33, close to windows.

So the within-gene loss reflects where variants fall in the protein, and the information lost is of the kind PM1 encodes.

**Estimator behaviour** (Supplementary Note 3; 100 simulated datasets per scenario at the real sample size): the within-unit estimator recovered planted within-gene signals with relative bias ≤ 1% at λ = 1, 2 and 4, for gene, window and random units. 95% coverage ranged from 0.86 to 0.99. A purely variant-level planted signal (λ = 4: true within-gene LR 3.11) still fell to 2.48 within windows, because AlphaMissense's top scores cluster in protein regions. That is a retention of about 0.80 relative to within-gene. The observed within-window loss in real data is much larger. Across unit sizes, real windows always retained far less than random partitions with the same number of units: median r 0.27 against 0.48 at 25 residues, 0.33 against 0.52 at 50, and 0.39 against 0.53 at 100 (common site set). Retention therefore falls with finer positional resolution, while random partitioning of any size costs little.

### Q3. Does the result replicate?
**Pre-specified tests:**

| Cohort / predictors | Informative | r(D2) < 0.5 | 95% CI | Supporting after full adjustment |
|---|---|---|---|---|
| DD, 36 predictors (primary; D1 position control) | 24 | 15 (62%): criterion not met | — | 2/36 (E1) |
| DD, D2 position control | 24 | 19 (79%) | 71–88% | 0/36 (E2) |
| Autism, 36 predictors (D2 pre-specified) | 20 | 16 (80%) | 60–95% | 0/36 |
| DD, 17 newer predictors (D2 pre-specified) | 14 | 13 (93%) | 93–100%\* | 0/17 |
| Public cohort, 36 predictors (descriptive; adequately powered) | 15 | 11 (73%) | — | 0/36 (E2 and E1) |

\*With 13 of 14 predictors below 0.5, the bootstrap share can only take a few values (13/14 or 14/14 in most replicates), so the percentile interval sits at the point estimate.

**Reading the binary criterion honestly.**

- The primary DD criterion failed under its pre-specified position control (D1). In genes without a constraint region, D1 does not control position.
- Under D2 the point estimates exceed 75% in every cohort, but the intervals for DD and autism include 75%.
- **The continuous measure is the more informative result, and it is stable.** Median retention after gene and position was 0.33 (DD), 0.37 (autism), 0.31 (newer predictors) and 0.28 (public cohort). The same four predictors (ClinPred, VEST4, MutationAssessor, MutPred) retained the most evidence in the DD, autism and public cohorts. Per-predictor retention correlated across cohorts (Spearman ρ 0.76).

**Negative control.** In 4,051 unaffected siblings, no predictor showed calibrated evidence at any level (0 of 37: the 36 dbNSFP 4.1a predictors plus AlphaMissense).

**Leakage.** Retention did not differ between predictors trained on clinical labels and label-free predictors (median r(D2) 0.29 against 0.37; Mann–Whitney p = 0.48). Leakage would in any case inflate, not reduce, within-gene evidence.

**Robustness** (Supplementary Note 2; median retention after gene and window, DD, 24 informative predictors):

| Analysis | Median r(D2) | AlphaMissense within-window LR |
|---|---|---|
| Baseline | 0.33 | 2.19 |
| CpG-transition rates × 0.919 | 0.33 | 2.17 |
| CpG-transition rates × 1.088 | 0.33 | 2.17 |
| Excluding 64,573 ClinVar P/LP sites (any review status) | 0.34 | 2.07 |

- **Neither mutation-rate model error nor known pathogenic variants drive the result:** excluding every ClinVar P/LP site leaves retention unchanged.
- **It is not an artefact of constraint-region coverage.** Retention was similar in genes with ≥ 2 constraint regions (0.32) and in genes without a regional split (0.37), where the within-gene to within-window drop cannot come from constraint regions.
- **It did not depend on gene constraint** (s_het high 0.35, low 0.36).

**Where PP3 adds the most.** Within DDG2P monoallelic genes, predictors retained more evidence beyond position: median r(D2) 0.47 in LoF genes and 0.57 in non-LoF (dominant-negative and gain-of-function) genes, where fewer than a third of predictors lost most of their evidence. Part of this is mechanical: within a curated gene set the marginal LR is smaller, so less of it is gene identity. It nonetheless indicates that within established NDD genes, and particularly for non-LoF mechanisms, predictor evidence is least redundant with position.

### Q4. What does this imply for combining PP3 with PM1 and PP2?
**Direct overlap test** (DD, AlphaMissense, PM1-like flag): R = 0.44 [0.41, 0.48], in the pre-committed "substantial overlap" tier (upper bound < 0.85). The joint LR was 15.0, against 33.9 if the two were multiplied.

- **At PM1-like sites,** AlphaMissense PP3 corresponds to LR 1.85: 0.70 points on the ClinGen scale (2.406 per point) or 0.84 on the Tavtigian scale (2.08 per point), against 1–4 points from the published intervals. Away from PM1-like sites it keeps LR 3.66.
- **Dose-response:** overlap was stronger for the strict PM1 definition (o/e < 0.2; R 0.37 [0.34, 0.41]) and weaker for PP2-like genes (0.61 [0.56, 0.66]).
- **Replication:**
  - The autism cohort replicated the primary result in the same tier: R 0.56 [0.44, 0.67]; PP3 at PM1-like sites 0.59 points.
  - In the public cohort (descriptive), R was 0.67 [0.42, 0.95].
  - The overlap was tool-dependent. For REVEL × PM1-like, R was 0.57 [0.52, 0.61] in DD but 0.82 [0.65, 1.00] in autism.
- Results for REVEL, VEST4, the strict PM1 definition and the PP2-like flag are in Table 2.

**Table 2. Direct overlap test** (PP3 at the ClinGen Supporting threshold; R = 1 means additive).

| Predictor | Flag | Cohort | R [95% CI] | PP3 LR given flag | Points (2.406 / 2.08 scale) | PP3 LR without flag |
|---|---|---|---|---|---|---|
| AlphaMissense | PM1-like | DD | 0.44 [0.41, 0.48] | 1.85 | 0.70 / 0.84 | 3.67 |
| AlphaMissense | PM1-like | autism | 0.56 [0.44, 0.67] | 1.68 | 0.59 / 0.71 | 2.50 |
| AlphaMissense | PM1-like | public3 | 0.67 [0.42, 0.95] | 1.76 | 0.64 / 0.77 | 2.08 |
| AlphaMissense | PM1-strict | DD | 0.37 [0.34, 0.41] | 1.57 | 0.51 / 0.61 | 4.07 |
| AlphaMissense | PM1-strict | autism | 0.50 [0.39, 0.62] | 1.51 | 0.47 / 0.57 | 2.88 |
| AlphaMissense | PM1-strict | public3 | 0.55 [0.33, 0.83] | 1.45 | 0.42 / 0.50 | 2.45 |
| AlphaMissense | PP2-like | DD | 0.61 [0.56, 0.66] | 2.56 | 1.07 / 1.28 | 3.96 |
| AlphaMissense | PP2-like | autism | 0.65 [0.52, 0.77] | 1.96 | 0.77 / 0.92 | 2.81 |
| AlphaMissense | PP2-like | public3 | 0.78 [0.48, 1.09] | 2.03 | 0.81 / 0.96 | 2.17 |
| REVEL | PM1-like | DD | 0.57 [0.52, 0.61] | 2.43 | 1.01 / 1.21 | 3.68 |
| REVEL | PM1-like | autism | 0.82 [0.65, 1.00] | 2.00 | 0.79 / 0.95 | 1.75 |
| REVEL | PM1-like | public3 | 0.78 [0.51, 1.13] | 2.25 | 0.92 / 1.11 | 2.46 |
| REVEL | PM1-strict | DD | 0.48 [0.44, 0.52] | 2.05 | 0.82 / 0.98 | 4.15 |
| REVEL | PM1-strict | autism | 0.73 [0.55, 0.96] | 1.78 | 0.66 / 0.79 | 2.27 |
| REVEL | PM1-strict | public3 | 0.54 [0.28, 0.87] | 1.56 | 0.51 / 0.61 | 2.83 |
| REVEL | PP2-like | DD | 0.73 [0.67, 0.79] | 3.10 | 1.29 / 1.55 | 3.85 |
| REVEL | PP2-like | autism | 0.93 [0.74, 1.13] | 2.27 | 0.93 / 1.12 | 1.76 |
| REVEL | PP2-like | public3 | 1.02 [0.67, 1.48] | 2.93 | 1.23 / 1.47 | 2.09 |
| VEST4 | PM1-like | DD | 0.67 [0.62, 0.71] | 2.00 | 0.79 / 0.95 | 2.72 |
| VEST4 | PM1-like | autism | 0.77 [0.62, 0.91] | 1.69 | 0.60 / 0.71 | 1.98 |
| VEST4 | PM1-like | public3 | 0.75 [0.48, 1.02] | 1.85 | 0.70 / 0.84 | 2.35 |
| VEST4 | PM1-strict | DD | 0.60 [0.55, 0.65] | 1.78 | 0.66 / 0.79 | 2.96 |
| VEST4 | PM1-strict | autism | 0.79 [0.62, 0.98] | 1.75 | 0.64 / 0.76 | 2.12 |
| VEST4 | PM1-strict | public3 | 0.63 [0.38, 0.96] | 1.57 | 0.51 / 0.61 | 2.45 |
| VEST4 | PP2-like | DD | 0.77 [0.71, 0.83] | 2.30 | 0.95 / 1.14 | 2.86 |
| VEST4 | PP2-like | autism | 0.80 [0.66, 0.95] | 1.77 | 0.65 / 0.78 | 2.08 |
| VEST4 | PP2-like | public3 | 0.74 [0.48, 0.98] | 1.83 | 0.69 / 0.82 | 2.53 |

**Overlap within disease genes** (A13). Where PP3 is used clinically, within known disease genes, the overlap persisted:

- **Within DDG2P monoallelic genes, comparing variants within the same gene:** R for AlphaMissense × PM1-like was 0.69 [0.68, 0.70] in DD and 0.74 [0.68, 0.80] in autism.
- **Across DDG2P genes:** R was 0.53 in DD and 0.59 in autism.
- **Within gene across all genes:** R was 0.75 [0.71, 0.78] in DD but 0.96 [0.84, 1.08] in autism.
- **Reading:** part of the overlap is between genes. Across all genes in autism, the within-gene overlap cannot be distinguished from additivity. Within DDG2P disease genes, where PP3 is applied clinically, a substantial within-gene component remains in both cohorts.

**Strong tier.** R with PP3 at its Strong threshold and a PM1-like flag:

| Predictor | DD | Autism | Public cohort |
|---|---|---|---|
| AlphaMissense | 0.33 [0.30, 0.36] | 0.44 | 0.48 |
| REVEL | 0.53 | 0.74 | 0.71 |
| VEST4 | 0.72 | 0.84 | 0.84 |

Because Strong-tier LRs are large, these R values still keep most of the log-evidence. The fraction kept, f, and the absolute joint points:

| Predictor | f, DD | f, autism | f, public | Joint points, DD / replication (2.406 scale) |
|---|---|---|---|---|
| AlphaMissense | 0.73 [0.71, 0.76] | 0.75 [0.67, 0.81] | 0.75 [0.57, 0.89] | 3.5 / 2.5–2.8 |
| REVEL | 0.84 [0.82, 0.87] | 0.90 [0.80, 0.99] | 0.89 [0.71, 1.04] | 3.9 / 3.0–3.2 |
| VEST4 | 0.91 [0.88, 0.93] | 0.94 [0.85, 1.02] | 0.94 [0.73, 1.09] | 3.6 / 2.9 |

f = 1 is additivity; ClinGen's cap corresponds to f = 0.67. In relative terms the overlap at the Strong tier is smaller than the reduction the cap imposes.

**Classification change with de novo evidence** (A12.4, DD; PS2 + PM2 + PM1-like + AlphaMissense PP3, Tavtigian points):

- **Uncapped:** 2,438 de novo missense variants reach Pathogenic (≥ 10 points), but the de novo excess implies only 84% of variants in that class are disease-associated, below what a Pathogenic classification requires.
- **Capped or conditional:** both rules move all 2,438 to Likely Pathogenic.
- **Caveat:** this is a class-level estimate on the de novo scale, not a per-variant posterior.

**Per-predictor guidance** (Supplementary Table 1). Multiplying 4 published points by each predictor's retention gives the points a predictor can justify when combined with position-based criteria:

| Predictor | Points |
|---|---|
| VEST4, MutPred, PROVEAN | About 2 |
| AlphaMissense, REVEL, and most others | About 1 |
| MPC, which is built from regional constraint | About 0 |

**ClinVar.** In 1,987,084 ClinVar missense variants:

- 30% of PP3-positive (AlphaMissense) P/LP variants lie at PM1-like sites.
- This enrichment is expected, partly because classifications use PP3 and PM1 themselves, so it is not used as evidence of overlap.
- If PP3 and PM1 were combined without ClinGen's cap, PP3 + PM1 alone would reach Likely Pathogenic for 14,343 VUS (AlphaMissense). The cap prevents this.
- **The overlap is concentrated in NDD genes.** In DDG2P monoallelic genes, 37.3% of AlphaMissense PP3-positive VUS lie at PM1-like sites, against 13.2% across all genes. Without the cap, PP3 + PM1 at Moderate alone would reach Likely Pathogenic for 7,101 VUS in these genes (13,245 with PM1 at Strong).
- **Stricter PM1 definition** (o/e < 0.2, Chao et al.'s Moderate level): the corresponding numbers are 9.7% and 2,148.
- **The cap and the conditional rule** both prevent PP3 + PM1 alone from reaching Likely Pathogenic in every setting (Supplementary Table 3).

## Discussion

**Main finding.** For de novo missense variants, most of what predictors contribute is information about the gene and the position of the variant in the protein. The evidence lost to position is not evidence that predictors are wrong: a predictor that scores a hotspot highly is correct, but it is restating what location already says. The finding concerns **overlap**.

**Not a statistical artefact.** It survives the comparison that would expose an artefact of conditioning (random partitions). It replicates across cohorts, phenotypes and predictor generations. It is stable across score thresholds.

**Implications for combining criteria.** The ACMG points framework treats PP3, PM1 and PP2 as additive with a single prior for every gene, so gene-level knowledge enters only through PP2 and PM1. Two consistent responses exist:

1. **Credit PP3 conditionally** when PM1 or PP2 is applied. The direct overlap test and the per-predictor retention shares give the size of the adjustment.
2. **Use gene-specific priors or calibrations** and condition predictor evidence on them. Chen et al. (2026) take this route with ClinVar.

**Additivity of PP3 and PM1.** The direct test bears on a disagreement in current guidance.

- CanVIG-UK permits PM1 and PP3 together because they are regarded as largely different evidence types.
- For de novo missense variants, our data do not support that premise. The joint evidence of AlphaMissense PP3 and a PM1-like flag was about half the product of their separate LRs, in three cohorts.
- **The overlap is consistent with the cap.** ClinGen limits PP3 + PM1 to Strong (Pejaver et al. 2022; restated as a no-double-counting rule by Stenton et al. 2024).
  - **Relative terms**, which do not depend on the scale of the de novo LR: overlap at the Strong tier removed 6–27% of the summed log-evidence, and the cap removes 33%. The cap is therefore large enough to absorb the overlap, and may be conservative for REVEL and VEST4.
  - **Uncapped co-use** overstated the joint evidence in every DD cell. In the classification analysis it produced Pathogenic calls in a class that was only 84% disease-associated.
  - **What we do not claim:** that the cap's exact size is optimal. The absolute joint points (2.5–3.9 against the cap's 4) are secondary, because the de novo LR is not on ClinGen's pathogenicity scale.
  - Expert-panel specifications that fold domain location into the computational criterion, as the ENIGMA BRCA1/BRCA2 VCEP does, are also consistent with the result.

**Relation to previous work.**

- **Prior benchmarks used de novos to rank predictors across genes:** AlphaMissense against DD de novo cases versus controls, and popEVE against DD de novos and unaffected siblings.
- **Gene identity dominates ClinVar benchmarks, and within-gene evaluation reorders predictors** (Harrizi et al. 2026). Our label-free results agree: models trained against a genome-wide null fell below AlphaMissense within gene, and models trained with a within-gene objective matched it without exceeding it significantly.
- **Calibration against gene-level depletion:** a thesis calibrated predictors to gene-level population depletion (Pathak 2026).
- **What is new here:**
  - predictor evidence expressed as LRs from the de novo excess, without variant-level labels, on the ACMG points scale;
  - its conditional decomposition by gene, position and constraint;
  - a direct test of the additivity the points framework assumes.

**Limitations.**

1. **Pre-specification and forking paths.**
   - The primary pre-specified DD criterion, a binary count of predictors losing most evidence, was not met (62% against 75%).
   - The fixed-window position control was adopted after that result.
   - We therefore present the continuous retention estimates and the direct overlap test as the findings, not the binary criterion.
   - The plan has thirteen amendments, each timestamped before the analysis it governs (Supplementary Note 1). The documents were not registered independently, and the amendments followed earlier results.
2. **Labels are not entirely absent.** No variant-level clinical labels are used, but strata use curated gene lists (DDG2P) and mechanism labels. Several of these genes were discovered in the same cohorts, so the DDG2P-stratum estimates are not independent of gene discovery. The all-genes analyses do not depend on them.
3. **Mutation-rate model.** Residual error is about 9% at CpG transitions. The δ margin and the perturbation analysis address it, but finer-scale error could bias small units.
4. **Penetrance and phenotype.** The de novo excess reflects variants of high penetrance for severe early-onset phenotypes. Retention measured here may not transfer to lower-penetrance or inherited variants.
5. **Transfer to ClinVar** (the conditional points) is an extrapolation beyond dominant de novo NDD variants.
6. **Ancestry.** The cohorts are predominantly of European ancestry.
7. **Precision.** The public cohort passed the pre-specified power check but gives wide intervals, and the model comparisons cannot exclude modest differences.

**Conclusion.** Missense predictors are most informative about which genes and regions matter, and much less about which variants matter within them. When PP3 is combined with PM1-like evidence, the joint evidence is substantially less than the sum, including within known disease genes, and the overlap is consistent with ClinGen's cap. Per-predictor retention shares and the direct overlap estimates indicate which predictors justify more or less than the cap.

---

## Tables and figures

- **Table 1.** Pre-specified tests and their outcomes, with denominators and intervals.
- **Table 2.** Direct overlap test: R and implied PP3 points for AlphaMissense, REVEL and VEST4 × PM1-like, PM1-strict and PP2-like flags, in DD, autism and the public cohort (shown in Results, Q4).
- **Figure 1.** Design: de novo excess, conditioning levels, retention.
- **Figure 2.** Evidence ladder: top-2% LR by level for AlphaMissense (DD, autism), popEVE (DDD) and MetaRNN.
- **Figure 3.** Per-predictor retention in three settings.
- **Figure 4.** Direct overlap test (`fig4_overlap.png`). R with 95% bootstrap intervals for AlphaMissense, REVEL and VEST4 PP3 combined with PM1-like, PM1-strict and PP2-like flags, in DD (primary), autism (replication) and the public cohort (descriptive). The dashed line is additivity (R = 1), which the points framework assumes.

**Supplementary:**

- Note 1: pre-registration timeline, amendments and deviation log.
- Note 2: robustness.
- Note 3: estimator simulation.
- Note 4: models trained on the de novo excess (Supplementary Figure 3: marginal against within-gene evidence).
- Table 1: per-predictor retention with intervals and implied points.
- Table 2: training-data ledger.
- Table 3: ClinVar sensitivity.
- Data 1: per-variant ClinVar PP3 + PM1 points.

---

## References (verified against source unless marked)
- Bergquist T, et al. *Genet Med* 2025: calibration of AlphaMissense, ESM1b and VARITY for PP3/BP4.
- Chao KR, et al. The landscape of regional missense mutational intolerance quantified from 125,748 exomes. bioRxiv 2024, doi:10.1101/2024.04.11.588920.
- Chen et al. Gene- and domain-aware calibration increases the clinical utility of variant effect predictors. Research Square preprint, May 2026.
- Cheng J, et al. AlphaMissense. *Science* 2023.
- Harrizi et al. Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. bioRxiv 2026, doi:10.64898/2026.08.21.746180.
- Havrilla JM, et al. Constrained coding regions. *Nat Genet* 2019.
- Jin SC, et al. Contribution of rare inherited and de novo variants in 2,871 congenital heart disease probands. *Nat Genet* 2017;49:1593–1601. doi:10.1038/ng.3970.
- Kaplanis J, et al. *Nature* 2020.
- Liu X, et al. dbNSFP v4. *Genome Med* 2020;12:103. doi:10.1186/s13073-020-00803-9.
- Orenbuch R, et al. popEVE. *Nat Genet* 2025.
- Pathak. PhD thesis, University of Edinburgh, 2026 (Edinburgh Research Archive).
- Pejaver V, et al. *Am J Hum Genet* 2022.
- Richards S, et al. *Genet Med* 2015; Tavtigian SV, et al. *Genet Med* 2018 (Bayesian framework); Tavtigian SV, et al. *Hum Mutat* 2020 (points).
- Stenton SL, et al. Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med* 2024 (states the PP3 + PM1 cap at Strong).
- Garrett A, et al. Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. *J Med Genet* 2021;58:297.
- Tejura M, et al. Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes. *Am J Hum Genet* 2024.
- Badonyi M, Marsh JA. acmgscaler: standardized gene-level variant effect score calibration within the ACMG/AMP framework. 2025 (PMC12496131).
- Grimm DG, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat* 2015.
- Seplyarskiy V, et al. Roulette. *Nat Genet* 2023.
- Turner TN, et al. denovo-db. *Nucleic Acids Res* 2017.
- Zhou X, et al. *Nat Genet* 2022. doi:10.1038/s41588-022-01148-2.
- ENIGMA BRCA1/BRCA2 VCEP specifications, *Genet Med Open* 2025 (to be confirmed against the published specification before submission).

---

## Supplementary Note 4. Models trained on the de novo excess

**Methods.**
**Training.** Gradient-boosted classifiers separate real de novo missense variants from sites drawn in proportion to mutation rate. Their inputs:

- 53 predictor scores;
- amino acids;
- regional constraint;
- AlphaFold structure (v2 onward).

Gene identity and mutation-rate features are excluded.

**Versions:**

| Model | Null sites drawn from | Training de novos | Status |
|---|---|---|---|
| v1, v2 | All genes | GeneDx | Pre-specified |
| v3 | The same transcript as each real de novo | GeneDx | Post hoc |
| v4 | The same transcript as each real de novo | All three DD centres | Post hoc |

v4 was evaluated only outside DD.

**Comparison with AlphaMissense:**

- Within-gene top-2% LR, paired bootstrap.
- Equivalence assessed with two one-sided tests, margin [0.80, 1.25].

**Summary.** Gradient-boosted models trained to separate real de novo missense variants from mutation-rate-matched sites had larger marginal LRs than any published predictor, but did not exceed AlphaMissense within gene. With a gene-matched null, they were statistically equivalent to AlphaMissense within gene in DD and autism. In the public cohort the point estimates were higher (v4: 1.30, 90% CI 1.03–1.80); this was the fifth, post hoc comparison and is exploratory.

**Genome-wide null (v1, v2).** Models trained against a genome-wide mutational null had the largest marginal LRs in the study (v2: 27.3), but less within-gene evidence than AlphaMissense:

| Comparison (v2 / AlphaMissense, within gene) | Ratio [95% CI] |
|---|---|
| DD cold genes | 0.69 [0.57, 0.84] |
| Autism | 0.81 [0.68, 0.95] |

The objective rewarded gene-level information.

**Gene-matched null (v3, v4).** Matching the null on gene raised within-gene evidence by 34–51% in every evaluation set. Within gene relative to AlphaMissense:

| Model | Cohort | Ratio [95% CI] |
|---|---|---|
| v3 | DD cold genes | 0.96 [0.82, 1.14] |
| v3 | Autism | 1.03 [0.89, 1.19] |
| v3 | Public cohort | 1.24 [0.91, 1.89] |
| v4 (all DD centres; evaluated outside DD only) | Autism | 1.09 [0.94, 1.29] |
| v4 | Public cohort | 1.30 [0.98, 1.91] |

- **Equivalence** (two one-sided tests, 90% CI within [0.80, 1.25]) was shown for v3 in DD, all genes (0.93 [0.89, 0.99]; the 0.96 above is for cold genes only, with a 95% interval) and in autism (1.03 [0.91, 1.17]), and for v4 in autism (1.09 [0.96, 1.25]).
- **In the public cohort,** v4's 90% interval lay above 1 (1.30 [1.03, 1.80]). As the fifth, post hoc model comparison, this is exploratory.
- The point estimates in the two independent cohorts favour the models, but no comparison was significant. These are the fourth and fifth model comparisons, two of them post hoc.
- **Conclusion:** a label-free model trained with a within-gene objective is at least comparable to AlphaMissense. Whether it is better requires a pre-registered test in a large independent NDD cohort.

---
