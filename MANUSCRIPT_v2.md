# Missense predictors add little evidence beyond gene and protein position for de novo variants: implications for combining PP3 with PM1 and PP2

*Draft 2, 6 October 2026. This draft is restructured in response to two external reviews; the point-by-point responses are in `REVIEW_RESPONSE.md`. Values marked ⟨pending⟩ come from analyses that are running and will be filled in from the committed outputs.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** ACMG/AMP classification adds computational evidence (PP3) to gene- and region-level evidence (PP2, PM1). We asked how much evidence missense predictors carry once the gene and the protein position of a variant are known, and what that implies for combining these criteria.

**Methods.** Likelihood ratios (LRs) were estimated from the excess of de novo missense variants over a per-site mutation-rate expectation, without variant-level clinical labels. Three cohorts were used:

- 31,058 developmental-disorder (DD) trios;
- 12,316 autism probands;
- a public multi-phenotype set of 4,256 trios.

We measured, for 55 predictors and for four models we trained on the same signal, the share of each predictor's log-evidence that remains after conditioning on gene, 50-residue protein windows and regional constraint ("retention"). We tested directly whether evidence from PP3 and PM1-like or PP2-like flags is additive. The analysis plan and its amendments were pre-specified and timestamped.

**Results.**

- **Gene identity and position account for most of the evidence.** Informative predictors retained a median of 0.56 of their log-evidence within gene, and 0.33 within 50-residue windows. AlphaMissense retained 0.30 [95% CI 0.28–0.32].
- **This is position, not small-unit conditioning.**
  - A random partition with the same number of units retained 0.51; the paired difference was 0.21 [0.19–0.23], positive for 24 of 24 predictors.
  - Curated InterPro domains retained 0.33.
- **Retention replicated** in the autism cohort (median 0.37), in 17 newer predictors (0.31) and in the public cohort (0.28).
- **Robustness:** retention was stable across score thresholds (and ⟨pending: under mutation-rate perturbation⟩), and the sibling negative control was null.
- **Supporting evidence after full adjustment:** none of the newer predictors (0/17) retained calibrated Supporting evidence after gene, position and constraint, and none in either replication cohort (0/36 each). In DD, 2/36 did under the originally pre-specified position control.
- **Direct overlap test:** in DD, the LR of AlphaMissense PP3 and a PM1-like flag together was ⟨pending⟩ of the product of their separate LRs. Implied PP3 evidence at PM1-like sites was ⟨pending⟩ points.
- **Models trained on the de novo excess** matched AlphaMissense within gene but did not exceed it significantly in any cohort.

**Conclusion.** For de novo missense variants, most predictor evidence is information that gene- and region-level criteria already carry. PP3 should be credited conditionally when PM1 or PP2 is applied. Per-predictor retention shares provide the size of the adjustment.

---

## Introduction

**The setting.** Missense variants are the largest class of variants of uncertain significance. Their classification under ACMG/AMP combines several kinds of evidence:

- computational predictions (PP3);
- location in a mutational hotspot or critical functional domain (PM1);
- membership of a gene with a low rate of benign missense variation (PP2).

**How the pieces are combined.** In the Bayesian points framework the evidence is additive in log-likelihood-ratio units, with a single prior probability for every gene. The ClinGen Sequence Variant Interpretation working group calibrated predictor thresholds to LRs using ClinVar (Pejaver et al. 2022; Bergquist et al. 2025) and, recognising that PP3 and PM1 can describe the same property, capped their combination at Strong.

**What is not known: how much the evidence overlaps.** A predictor that has learned which genes are intolerant, or which domains are hotspots, separates pathogenic from benign variants pooled across genes. Yet it may add little once the gene and the domain are known. ClinVar-based calibration cannot easily measure this:

- its labels are concentrated in a few genes;
- several predictors were trained on them;
- the evaluation is dominated by gene identity (Harrizi et al. 2026).

**Our approach.** Here we use de novo variants, whose excess over mutation-rate expectation in affected children defines a likelihood ratio without variant-level labels, and which can be estimated *within* a gene, window or domain. We ask four questions:

1. How much of each predictor's evidence survives conditioning on gene, protein position and regional constraint?
2. Is the loss due to position, or to the act of conditioning?
3. Does it replicate across cohorts, phenotypes and predictor generations?
4. What does it imply for combining PP3 with PM1 and PP2?

We also ask whether a model trained directly on this signal can do better.

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
- **Regional missense constraint:** gnomAD v2.1.1, lifted to GRCh38 (Chao et al. 2024).
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

### Models trained on the de novo excess
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
**Random partitions retained nearly all the within-gene evidence.** This control used all analysed sites; the main analyses used the common site set, which is why its values differ slightly from Q1 (0.31 against 0.33). A random partition of each gene's residues into the same number of groups as the 50-residue windows retained a median 0.51, against 0.56 for gene alone. Real windows retained 0.31.

- The paired difference was 0.21 [0.19, 0.23], positive for all 24 predictors.
- Shifting window boundaries by 25 residues changed nothing (−0.02).
- Curated InterPro domains retained 0.33, close to windows.

So the within-gene loss reflects where variants fall in the protein, and the information lost is of the kind PM1 encodes.

**Estimator behaviour** under a planted within-gene signal is in Supplementary Note 3 (⟨pending⟩). Retention across 25-, 50- and 100-residue units is in Supplementary Figure 2 (⟨pending⟩).

### Q3. Does the result replicate?
**Pre-specified tests:**

| Cohort / predictors | Informative | r(D2) < 0.5 | 95% CI | Supporting after full adjustment |
|---|---|---|---|---|
| DD, 36 predictors (primary; D1 position control) | 24 | 15 (62%): criterion not met | — | 2/36 (E1) |
| DD, D2 position control | 24 | 19 (79%) | 71–88% | 0/36 (E2) |
| Autism, 36 predictors (D2 pre-specified) | 20 | 16 (80%) | 60–95% | 0/36 |
| DD, 17 newer predictors (D2 pre-specified) | 14 | 13 (93%) | 93–100% | 0/17 |
| Public cohort, 36 predictors (descriptive; adequately powered) | 15 | 11 (73%) | — | 0/36 (E2 and E1) |

**Reading the binary criterion honestly.**

- The primary DD criterion failed under its pre-specified position control (D1). In genes without a constraint region, D1 does not control position.
- Under D2 the point estimates exceed 75% in every cohort, but the intervals for DD and autism include 75%.
- **The continuous measure is the more informative result, and it is stable.** Median retention after gene and position was 0.33 (DD), 0.37 (autism), 0.31 (newer predictors) and 0.28 (public cohort). The same four predictors (ClinPred, VEST4, MutationAssessor, MutPred) retained the most evidence in the DD, autism and public cohorts. Per-predictor retention correlated across cohorts (Spearman ρ 0.76).

**Negative control.** In 4,051 unaffected siblings, no predictor showed calibrated evidence at any level (0/37).

**Leakage.** Retention did not differ between predictors trained on clinical labels and label-free predictors (median r(D2) 0.29 against 0.37; Mann–Whitney p = 0.48). Leakage would in any case inflate, not reduce, within-gene evidence.

**Robustness** (Supplementary Note 2, ⟨pending⟩):

- perturbing CpG-transition rates by ±9%;
- excluding all ClinVar P/LP sites;
- stratifying by gene constraint, mechanism and constraint-region coverage.

### Q4. What does this imply for combining PP3 with PM1 and PP2?
**Direct overlap test** (DD, AlphaMissense, PM1-like flag): R = ⟨pending⟩ [⟨pending⟩].

- Applying PP3 at a PM1-like site corresponds to ⟨pending⟩ points of evidence, against 1–4 points from the published intervals.
- The autism cohort ⟨pending: replicated / did not replicate⟩.
- Results for REVEL, VEST4, the strict PM1 definition and the PP2-like flag are in Table 2 (⟨pending⟩).

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
- Sensitivity to PM1 strength, to the strict PM1 definition and to restriction to DDG2P monoallelic genes is in Supplementary Table 3 (⟨pending⟩).

### Q5. Can a model trained on the de novo excess do better?
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

- Equivalence within [0.80, 1.25] was ⟨pending⟩.
- The point estimates in the two independent cohorts favour the models, but no comparison was significant. These are the fourth and fifth model comparisons, two of them post hoc.
- **Conclusion for Q5:** a label-free model trained with a within-gene objective is at least comparable to AlphaMissense. Whether it is better requires a pre-registered test in a large independent NDD cohort.

---

## Discussion

**Main finding.** For de novo missense variants, most of what predictors contribute is information about the gene and the position of the variant in the protein. The evidence lost to position is not evidence that predictors are wrong: a predictor that scores a hotspot highly is correct, but it is restating what location already says. The finding concerns **overlap**.

**Not a statistical artefact.** It survives the comparison that would expose an artefact of conditioning (random partitions). It replicates across cohorts, phenotypes and predictor generations. It is stable across score thresholds.

**Implications for combining criteria.** The ACMG points framework treats PP3, PM1 and PP2 as additive with a single prior for every gene, so gene-level knowledge enters only through PP2 and PM1. Two consistent responses exist:

1. **Credit PP3 conditionally** when PM1 or PP2 is applied. The direct overlap test and the per-predictor retention shares give the size of the adjustment.
2. **Use gene-specific priors or calibrations** and condition predictor evidence on them. Chen et al. (2026) take this route with ClinVar.

ClinGen's cap of PP3 + PM1 at Strong is directionally consistent with our results. Our estimates suggest that for most predictors the combination justifies less than the cap allows, and for a few predictors (VEST4, MutPred) about as much. Expert-panel specifications that fold domain location into the computational criterion, as the ENIGMA BRCA1/BRCA2 VCEP does, are also consistent with the result.

**Relation to previous work.**

- **Prior benchmarks used de novos to rank predictors across genes:** AlphaMissense against DD de novo cases versus controls, and popEVE against DD de novos and unaffected siblings.
- **Gene identity dominates ClinVar benchmarks, and within-gene evaluation reorders predictors** (Harrizi et al. 2026). Our label-free results agree: models trained against a genome-wide null fell below AlphaMissense within gene, and models trained with a within-gene objective matched it without exceeding it significantly.
- **Calibration against gene-level depletion:** a thesis calibrated predictors to gene-level population depletion (Pathak 2026).
- **What is new here:**
  - predictor evidence expressed as LRs from the de novo excess, without variant-level labels, on the ACMG points scale;
  - its conditional decomposition by gene, position and constraint;
  - a direct test of the additivity the points framework assumes.

**Limitations.**

1. **Pre-specification.** The primary pre-specified DD test did not meet its criterion. The position control used afterwards was adopted after that result. Documents were timestamped but not registered independently.
2. **Labels are not entirely absent.** No variant-level clinical labels are used, but strata use curated gene lists (DDG2P) and mechanism labels. Several of these genes were discovered in the same cohorts, so the DDG2P-stratum estimates are not independent of gene discovery. The all-genes analyses do not depend on them.
3. **Mutation-rate model.** Residual error is about 9% at CpG transitions. The δ margin and the perturbation analysis address it, but finer-scale error could bias small units.
4. **Penetrance and phenotype.** The de novo excess reflects variants of high penetrance for severe early-onset phenotypes. Retention measured here may not transfer to lower-penetrance or inherited variants.
5. **Transfer to ClinVar** (the conditional points) is an extrapolation beyond dominant de novo NDD variants.
6. **Ancestry.** The cohorts are predominantly of European ancestry.
7. **Power.** The public cohort has limited power, and the model comparisons are underpowered to detect modest differences.

**Conclusion.** Missense predictors are most informative about which genes and regions matter, and much less about which variants matter within them. Classification frameworks should credit predictor evidence conditionally on the gene- and region-level criteria applied alongside it. Per-predictor retention shares provide a practical, label-free way to do so.

---

## Tables and figures

- **Table 1.** Pre-specified tests and their outcomes, with denominators and intervals.
- **Table 2.** Direct overlap test: R and implied PP3 points for AlphaMissense, REVEL and VEST4 × PM1-like, PM1-strict and PP2-like flags, in DD, autism and the public cohort. ⟨pending⟩
- **Figure 1.** Design: de novo excess, conditioning levels, retention.
- **Figure 2.** Evidence ladder: top-2% LR by level for AlphaMissense (DD, autism), popEVE (DDD) and MetaRNN.
- **Figure 3.** Per-predictor retention in three settings.
- **Figure 4.** Direct overlap test (replaces the ClinVar counterfactual as the main clinical figure). ⟨pending⟩
- **Figure 5.** Models: marginal against within-gene evidence.

**Supplementary:**

- Note 1: pre-registration timeline, amendments and deviation log.
- Note 2: robustness.
- Note 3: estimator simulation.
- Table 1: per-predictor retention with intervals and implied points.
- Table 2: training-data ledger.
- Table 3: ClinVar sensitivity.
- Data 1: per-variant ClinVar PP3 + PM1 points.

---

## References (all verified against source)
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
- Richards S, et al. *Genet Med* 2015; Tavtigian SV, et al. *Hum Mutat* 2020.
- Seplyarskiy V, et al. Roulette. *Nat Genet* 2023.
- Turner TN, et al. denovo-db. *Nucleic Acids Res* 2017.
- Zhou X, et al. *Nat Genet* 2022. doi:10.1038/s41588-022-01148-2.
- ENIGMA BRCA1/BRCA2 VCEP specifications, *Genet Med Open* 2025 (to be confirmed against the published specification before submission).
