# How much evidence do missense predictors add beyond gene and protein position? A pre-registered, label-free calibration using de novo mutations in 43,000 trios

*Draft 1, 5 October 2026. Target: Genetics in Medicine (structured abstract) or AJHG. Items marked **[PENDING]** are produced by analyses that are frozen in the pre-registration and currently running; every other number comes from committed outputs in `results_fullstudy/`.*

**Authors:** [names, affiliations]

---

## Abstract

**Purpose.** Computational predictions of missense variant effect (ACMG/AMP criterion PP3) are calibrated against ClinVar and then combined with gene- and region-level criteria such as PP2 and PM1. If a predictor's signal is largely "which gene" and "where in the protein", that combination counts the same evidence twice. ClinVar cannot measure this cleanly, because its labels are concentrated in a few genes and many predictors were trained on them. We measured, without clinical labels, how much evidence predictors add once gene, protein position and regional constraint are held fixed.

**Methods.** We estimated likelihood ratios (LRs) from the excess of de novo missense variants over a per-site mutation-rate expectation in 31,058 developmental-disorder (DD) trios. LRs were computed at five nested levels:

- A: none;
- B: gene;
- C: regional missense constraint;
- D: 50-residue protein windows;
- E: position plus constraint.

We applied this to 55 predictors: 36 from dbNSFP 4.1a, AlphaMissense, 17 dbNSFP 5.4a tools, and popEVE, the last evaluated on held-out DDD trios only. The plan and five amendments were timestamped before the corresponding results. We replicated the analysis in 12,316 autism probands, with 4,051 unaffected siblings as a negative control.

**Results.**

- **Gene identity:** informative predictors kept a median of 0.56 of their log-evidence after conditioning on gene.
- **Position:** 0.33 after conditioning on gene and 50-residue window. A random partition of each gene into the same number of non-contiguous groups left 0.51, and real windows removed more evidence than the random partition for all 24 predictors, so the loss reflects protein position rather than conditioning on small units.
- **Constraint:** 0.32 after adding regional constraint.
- **Specific predictors:** AlphaMissense's top-2% LR fell from 13.4 to 2.08, and popEVE's from 17.9 to 1.69.
- **Pre-registered tests.** The primary pre-registered test, in the DD cohort, did not meet its criterion under the originally specified position control: 15 of 24 informative predictors (62%) lost most of their evidence to gene and position, against a pre-specified 75%. Under fixed 50-residue windows the share was 19/24 (79%). With fixed windows made primary before the data were analysed, the criterion was met in the autism replication (16/20, 80%) and among the 17 newer predictors in the DD cohort (13/14, 93%).
- **Supporting evidence after full adjustment:**
  - DD cohort: 2 of 36 predictors kept calibrated Supporting-level evidence under the originally specified level (VEST4, MutPred; E1), and 0 of 36 under fixed windows (E2).
  - Replication: 0 of 36.
  - Newer predictors: 0 of 17.
- **Within known genes:** no predictor reached Moderate. De novo-derived thresholds were stricter than ClinVar-based ones for every tool examined.
- **ClinVar impact.** Among 1.99 million ClinVar missense variants, 30% of AlphaMissense-PP3-positive pathogenic variants lie in PM1-like regions. Applying PP3 and PM1 together at full strength would by themselves classify 14,343 VUS as Likely Pathogenic (AlphaMissense; 4,065 with REVEL). With each predictor's measured retention share, the combined PP3 + PM1 evidence falls below ClinGen's Strong cap for low-retention predictors (AlphaMissense, REVEL), but not for all VEST4 variants.

**Conclusion.** Most of what current missense predictors contribute to de novo variant interpretation is information about the gene and the protein region. Applying PP3 at full strength alongside PP2 or PM1 risks double counting. Predictor evidence should be calibrated conditionally, or capped when combined with gene- and region-level criteria.

---

## Introduction

Missense variants are the most common class of variant of uncertain significance. Their classification under the ACMG/AMP framework leans heavily on computational evidence. The ClinGen Sequence Variant Interpretation working group calibrated predictor score thresholds to likelihood ratios for each PP3/BP4 evidence strength (Pejaver et al. 2022), and later extended this to AlphaMissense, ESM1b and VARITY (Bergquist et al. 2025). These calibrations allow predictors to contribute up to Strong evidence.

The calibrations use ClinVar pathogenic and benign variants pooled across genes. This creates two problems for interpretation.

**First, gene and region effects.** Pathogenic variants are concentrated in particular genes and protein regions. A predictor that learns "this gene is intolerant" or "this domain is a hotspot" will separate pathogenic from benign variants across genes. Yet it may add little once the gene and region are already known, which is the situation in which a clinician applies PP3. ACMG/AMP already scores gene- and region-level evidence separately, through PP2 (missense-constrained genes) and PM1 (mutational hotspots and critical domains). ClinGen recognised the overlap with PM1 by capping PP3 + PM1 at Strong (Pejaver et al. 2022). The size of the overlap has not been measured.

**Second, circularity.** Many predictors are trained on ClinVar or HGMD labels, so ClinVar-based evaluation is partly in-sample.

**De novo mutations offer a label-free alternative.** In children with severe developmental disorders, de novo missense variants are enriched over the number expected from mutation alone. The size of that excess at a given score, relative to expectation, defines a likelihood ratio without any clinical classification. Because expectations can be computed per site, the excess can be measured *within* a gene, region or protein window. That directly asks how much evidence a predictor adds beyond those units.

We used this design to measure the conditional evidence of 55 missense predictors. The analysis was pre-registered, replicated in an independent autism cohort, and tested on a set of newer predictors that had not been examined when the main hypotheses were fixed.

---

## Methods

### Pre-registration and analysis governance

**Freezing procedure.** A feasibility gate, the full study plan and five amendments were each committed to a git repository and OpenTimestamped before the analyses they govern were run:

| Document | Commit | Content |
|---|---|---|
| Gate | `5c60cf7` | Feasibility gate |
| Plan | `6e2ee03` | Full study plan |
| A1 | `c719226` | Nested levels |
| A2 | `c5b94e5` | Autism replication |
| A3 | `4d58eb2` | dbNSFP 5.4a tools |
| A4 and A5 | `fffc57c` | ClinVar impact and position controls; null-contrast model v2 |

- Each amendment records what had and had not been seen when it was written.
- Deviations are logged with dates and reasons in `FULL_STUDY_PLAN.md`.

**Implementation check.** The analysis code reproduces every estimate of the frozen feasibility-gate implementation to a relative difference of ≤ 1e-11. This is a consistency check between implementations, not a validation of the method. The evidence for calibration comes from three checks:

- synonymous and mutation-class O/E (gate G1);
- consistency across centres (G3);
- the unaffected-sibling negative control.

### Cohorts

**Primary DD cohort.** De novo variants from 31,058 trios with severe developmental disorders (Kaplanis et al. 2020), from three centres: DDD, GeneDx and Radboudumc (RUMC).

- DDD was designated the held-out centre for anything fitted.
- Variants were lifted to GRCh38 and matched to a per-site universe built on canonical GENCODE v32 transcripts.

**Replication cohort.** Published de novo SNVs from SPARK, SSC, ASC and MSSNG (Zhou et al. 2022).

- 12,316 autism probands, with unaffected siblings (n = 4,051) as a negative control.
- Individuals sharing ≥ 2 identical de novo SNVs with a DD proband were removed as likely duplicates (n = 64).

**Site filters.** Sites within 3 bp of an exon boundary were excluded on both the observed and expected side. This is a pre-registered fix for a splice-edge artefact found at the gate. After filters, the DD analysis set contains 23,851 missense and 8,222 synonymous de novos.

### Expected counts and likelihood ratios

**Expected counts.**

- Each site's mutation rate comes from the Roulette model.
- A cohort-specific scale factor *c* = (observed synonymous de novos) / (summed synonymous mutation rate) converts rates to expected counts, E = *c*·ΣMR.

**Marginal LR.** For a score group *b*:

LR_b = (X_b / ΣX) / (E_b / ΣE), where X = O − E is the excess.

**Within-unit LR.** For a unit type *u* (gene, region, window), each observed de novo contributes its unit's opportunity profile across score groups. Writing *O_u* and *E_u* for the observed and expected counts of unit *u*:

LR_b = X_b / (Σ_u O_u · E_ub / E_u − E_b)

This compares the score distribution of de novos with that of other sites in the same unit.

**Uncertainty.**

- 95% intervals come from 1,000 Poisson bootstrap replicates over probands (seed 20261004).
- A pre-specified systematic margin δ = 0.091 (the largest calibration error of the mutation-rate model at the gate) is applied to give one-sided lower bounds.
- Evidence strengths use the LR cut-offs of Pejaver et al. 2022: Supporting 2.406, Moderate 5.79, Strong 33.53.

**Score groups.** Scores are binned into 1,000 tie-robust quantiles of mutational opportunity.

- The headline readout for every predictor is its **top 2%**: the smallest set of top bins holding at least 2% of opportunity.
- Local-window thresholds use 51-bin windows with a monotone rule.

### Nested conditioning levels

| Level | Held fixed | Units |
|---|---|---|
| A | nothing | marginal |
| B | gene | canonical transcript |
| C | regional constraint | gene × gnomAD regional missense constraint region (RMC; lifted to GRCh38) |
| D1 | position, matched count | gene × evenly spaced chunks, as many as the gene's RMC regions |
| D2 | position | gene × fixed 50-residue window |
| E1 / E2 | position + constraint | D1 or D2 unit × RMC region |

**Retention share.** For each predictor, r(level) = ln LR(level) / ln LR(A). This is the share of its log-evidence (equivalently, of its ACMG points) that survives conditioning.

**Pre-specified readouts (amendments A1–A3):**

1. **Generalisation.** The share of informative predictors (marginal top-2% LR ≥ 2.406) with r < 0.5 after gene + position. "Systematic" was defined as ≥ 75%.
2. **Constraint beyond position:** r(D) − r(E).
3. **Supporting after full adjustment:** the number of predictors whose level-E one-sided lower bound reaches 2.406.

The position control D1 was primary in the DD cohort (A1). D2 was made primary for the replication (A2) and the newer predictors (A3), for a reason recorded before either was run: D1 has a single chunk in the ~13,500 genes without an RMC region, so it does not control position there.

### Predictors

- **dbNSFP 4.1a (36 tools):** rankscores.
- **AlphaMissense:** the published scores.
- **dbNSFP 5.4a (17 tools):** MetaRNN, MutPred2, gMVP, MisFit_D/S, VARITY_R/ER and their LOO versions, ESM1b, PHACTboost, MutFormer, MutScore, GPN-MSA, GERP (92 mammals), phyloP and phastCons (470 mammals).
- **popEVE (5.4a):** it was developed and evaluated using DD cohorts, so it was assessed on DDD de novos only, with DDD's own scale factor.
- **Training-data ledger:** for every tool we recorded whether its training labels could include variants from these cohorts (`aim3_training_ledger.tsv`).

### Clinical impact on ClinVar (amendment A4.1)

**Variants.** ClinVar GRCh38 missense SNVs (≥ 1 review star; P/LP, B/LB or VUS) matched to analysed sites.

**PM1-like flag.** RMC missense o/e < 0.4 or CCR ≥ 95. A PM1-like site is assumed to receive PM1 at Moderate (+2).

**Points for AlphaMissense, REVEL and VEST4, under three rules:**

1. **Published:** points from the published ClinGen intervals, with PP3 + PM1 capped at 4.
2. **Conditional:** published points × the predictor's r_E2, i.e. the overlapping share removed, with no cap.
3. **De novo thresholds:** this study's within-unit E2 thresholds.

### Position-unit controls (amendment A4.2)

- **win50s:** 50-residue windows shifted by 25 residues.
- **rand50:** a random partition of each gene's residues into groups of the same number as D2, with no contiguity.
- **dom:** gene × InterPro domain, a direct PM1 proxy.

All three use AlphaMissense and the 24 informative DD predictors, all genes, and the unchanged estimators. The random partition assigns each residue to one of ⌈L/50⌉ groups of its gene with a seeded hash, so it keeps D2's unit count while removing position. If the drop from gene to window were an artefact of conditioning on many small units, r(rand50) would be close to r(D2). If it reflects position, r(rand50) would be close to r(B).

### Null-contrast models (Aim 4 and amendment A5)

**Training data.** Real GeneDx missense de novos against 10× as many sites sampled in proportion to mutation rate, restricted to 80% of genes; the other 20% are "cold" genes.

**Features.** Predictor scores, amino acids and regional constraint. Gene identity and mutation-rate features are deliberately excluded.

- **v1:** gradient-boosted trees, fixed hyperparameters.
- **v2:** adds the dbNSFP 5.4a tools (excluding popEVE) and AlphaFold structure features: pLDDT and Cα contact density, mapped by sequence identity.

**Models compared:** gradient boosting (primary), a five-seed MLP ensemble, and their average.

**LR conversion.** Classifier output is converted to an LR via the density ratio and the GeneDx pathogenic fraction.

**Evaluation.** On held-out DDD de novos, with the paired bootstrap, on all genes and on cold genes.

### Software and data availability

- Code: `gate.py`, `fullstudy.py`, `af_structure.py`.
- Pre-registration documents and OpenTimestamps proofs: in the repository.
- All inputs are public, except dbNSFP 5.4a (registration required) and the DD and autism de novo calls, which are as published.

---

## Results

### A label-free calibration that reproduces known structure

**Feasibility gate:**

- The pipeline passed its pre-registered checks:
  - transcripts 96.7%;
  - de novos mapped 93.4%;
  - REF mismatch 0.003%;
  - consistency across the three centres: Cochran's Q p ≥ 0.16 in the top bins.
- Mutation-rate calibration was within ±10% in every group. CpG transitions were over-predicted (O/E 0.919), which sets δ = 0.091.

**Excess and priors:**

- The de novo missense excess implies that 20.0% [18.1, 22.0] of missense de novos in DD probands are pathogenic.
- That share is 72.9% [71.8, 73.9] in DDG2P monoallelic genes, and 94.1% at AlphaMissense ≥ 0.99 in those genes.

**Marginal LRs** reproduce the expected ranking: AlphaMissense's top-2% LR is 13.4 (Aim 1 bin 0.99–1.00: 9.04 [8.26, 9.86]).

**Negative control.** In autism siblings, the missense O/E was 1.106 [1.045, 1.175]. This is the same ~10% offset seen at the gate and is absorbed by δ. Among AlphaMissense ≥ 0.99 sites, siblings showed no excess (O/E 0.97 [0.81, 1.14]), against 1.86 [1.73, 1.99] in autism probands.

### Gene identity accounts for about half of predictor evidence

Of 36 predictors, 24 were informative in the DD cohort (marginal top-2% LR ≥ 2.406).

- **Within gene:** they kept a median of 0.56 of their log-evidence (Figure 2, Figure 3).
- **Loss scales with how much a tool encodes gene- or region-level information:**
  - MPC, built from regional constraint: marginal 14.2, within gene 2.0;
  - DEOGEN2: 11.8 → 2.8;
  - AlphaMissense: 15.4 → 4.3;
  - REVEL: 9.0 → 3.6.
- **Conservation-only scores** added almost nothing within gene (LR 0.8–1.9).

### Position explains most of what remains; regional constraint adds little

**Effect of each level** (median r across informative DD tools):

| Level | Median r |
|---|---|
| B: gene | 0.56 |
| C: RMC region | 0.38 |
| D2: 50-residue window | 0.33 |
| E2: window + constraint | 0.32 |

**Constraint beyond position.** The median gain from adding constraint once position is fixed, r(D) − r(E), was 0.073 under D1 and 0.01–0.02 under D2.

**Arbitrary partitions do most of the work.** The pre-registered primary readout for AlphaMissense (bin 0.99–1.00) went from 3.37 within gene to 2.25 [2.16, 2.35] within gene × RMC region. Its lower bound of 2.10 falls short of Supporting. Evenly spaced pseudo-regions with the same unit count gave 2.50, so most of the step from gene to region does not depend on where the constraint boundaries lie.

**The drop is position, not unit size** (amendment A4.2). Conditioning on many small units could in principle shrink any predictor's within-unit signal, so we compared real 50-residue windows with a random partition of each gene's residues into the same number of groups. Median r over the 24 informative predictors:

| Unit | Median r |
|---|---|
| Gene | 0.56 |
| Random partition (rand50) | 0.51 |
| Real windows | 0.31 |
| Windows shifted by 25 residues | 0.29 |

- **Small units cost little.** Random partitions removed only a small part of the evidence. The paired difference r(rand50) − r(D2) was positive for **all 24 predictors** (median 0.21, IQR 0.18–0.22).
- **AlphaMissense:** within-gene top-2% LR 4.29, random partition 3.80, real windows 2.20.
- **Not a boundary artefact:** shifting window boundaries by half a window changed nothing (median difference −0.016).
- **Conclusion:** most of the loss from gene to window is due to where variants fall in the protein, not to the conditioning itself.

[PENDING: InterPro domain units, run last]

### Pre-registered generalisation tests

| Test | Data | Informative tools | r < 0.5 after gene + position | Pre-specified verdict | Tools at Supporting after full adjustment |
|---|---|---|---|---|---|
| A1 | DD, 36 tools | 24 | 62% (D1, primary); 79% (D2) | Not met | 2/36 (E1: VEST4, MutPred) |
| A2 | Autism, 36 tools | 20 | 80% (D2, primary); 75% (D1) | **Systematic** (adequately powered) | 0/36 (E2 and E1) |
| A3 | DD, 17 new tools | 14 | 93% (D2, primary); 57% (D1) | **Systematic** | 0/17 |

**Reading the DD result honestly.** The DD readout under its pre-registered position control did not meet the threshold.

- D1's limitation was understood before either later test: in genes without an RMC region it is not a position control.
- Under D2, every analysis reached 79–93%.

**Consistency across cohorts.**

- Retention shares agreed between cohorts (r(D2), Spearman ρ = 0.76 over the 20 tools informative in both).
- Raw rankings also agreed (ρ = 0.87).
- The tools that kept the most evidence were the same in both: MutationAssessor, VEST4, MutPred and ClinPred, with PROVEAN in DD.
- Median r was 0.56 / 0.33 / 0.32 (B / D2 / E2) in DD, 0.57 / 0.37 / 0.35 in autism, and 0.54 / 0.31 / 0.30 for the newer tools.

### The newest predictors behave the same way

**Strongest raw scorers lose the most.**

- MisFit_S (marginal top-2% LR 16.2), MutScore (15.1), PHACTboost (15.0) and gMVP (12.1) kept 20–24% of their log-evidence within 50-residue windows.
- MetaRNN was the only new tool to keep half (r(D2) = 0.505). It ends at E2 = 2.41, with a one-sided lower bound of 2.12, so it misses Supporting.

**popEVE** (DDD only) had the largest raw LR of any predictor (17.9). That fell to:

- 3.13 within gene;
- 1.76 within window;
- 1.69 after constraint (lower bound 1.58).

### Calibrated within-gene thresholds

**No predictor reached Moderate within gene.**

- Using local windows within gene in DDG2P genes, 15 of 37 tools (the 36 dbNSFP 4.1a tools and AlphaMissense) reached Supporting.
- VEST4 came closest to Moderate: its top-2% lower bound was 5.63, against the 5.79 cut-off.

**De novo thresholds were stricter than ClinVar-based ones for every tool examined:**

| Tool | De novo, marginal Supporting | De novo, within DDG2P genes | Published ClinGen Supporting |
|---|---|---|---|
| REVEL | ≥ 0.773 | ≥ 0.876 | 0.644 |
| BayesDel | ≥ 0.330 | ≥ 0.436 | 0.13 |
| CADD | 29.2 | not reached | 25.3 |

**DDD-only sensitivity for the newer predictors.** Restricting to DDD de novos, with DDD's own scale factor, gave the same result:

- 14 of 15 informative predictors kept < 50% of their evidence after gene and position (median r(D2) 0.33);
- none reached Supporting at E2;
- within-gene LRs were slightly higher than in the pooled cohort (e.g. MetaRNN 4.72 against 4.33).

### Mechanism dependence

- **Most tools depend on mechanism.** 17 of 37 tools gave significantly different within-gene evidence in loss-of-function genes than in dominant-negative or gain-of-function genes (Benjamini–Hochberg q < 0.05), using literature-curated mechanism labels.
  - 15 gave more evidence in LoF genes (e.g. PROVEAN, Δ log-LR 0.89 [0.51, 1.36]; REVEL, 0.45 [0.21, 0.72]).
- **AlphaMissense did not** (0.11 [−0.02, 0.24], q = 0.16). Its pre-registered LoF vs non-LoF test showed no difference.
- **Gain-of-function genes:** AlphaMissense's top-2% LR was flatter, 3.18 against 6.25 in dominant-negative genes.

### Training-label leakage does not explain the results

**The risk.** 28 of 55 tools were trained on clinical labels (ClinVar, HGMD, UniProt) that could include variants from these cohorts. Leakage would inflate within-unit evidence, i.e. bias *against* our conclusion.

**The check** (exploratory, not pre-specified). Retention after position did not differ between clinically trained and label-free tools: median r(D2) 0.29 against 0.37 (Mann–Whitney p = 0.48).

### Clinical impact on ClinVar missense classifications

**Variants.** 1,987,084 ClinVar missense variants (release 2026-10-04, ≥ 1 star) mapped to analysed sites:

- 45,997 P/LP;
- 1,825,639 VUS;
- 115,448 B/LB.

5.0% of them are at PM1-like sites (Figure 4).

**PP3 and PM1 overlap.**

- Among P/LP variants that receive PP3 from AlphaMissense under the published intervals, **30.3%** lie at PM1-like sites, six times the background rate (REVEL 24.3%, VEST4 24.4%).
- The same holds for 13.2% of PP3-positive VUS.

**Double-counting risk without the cap.** If PP3 at its published strength and PM1 at Moderate are both applied without ClinGen's cap, PP3 + PM1 alone reach Likely Pathogenic (≥ 6 points) for:

| Tool | VUS | P/LP | B/LB |
|---|---|---|---|
| AlphaMissense | 14,343 | 7,035 | 100 |
| REVEL | 4,065 | 3,471 | 21 |
| VEST4 | 3,004 | 2,441 | 17 |

ClinGen's cap removes this by construction (0 variants), as does the conditional rule (0).

**The cap is generous for de novo evidence.** At PM1-like sites, the conditional rule (published points × the predictor's r_E2) gives mean PP3 + PM1 of:

- AlphaMissense: 2.89 points for P/LP, against 3.77 under the capped published rule;
- REVEL: 2.69 against 3.50;
- VEST4: 3.26 against 3.52.

**How far the cap is from the conditional estimate depends on retention.**

- For predictors with r_E2 < 0.5, any PP3 + PM1 combination that reaches Strong under the cap falls below Strong under the conditional rule. This is arithmetic, not a separate finding: 4 × r_E2 + 2 < 4 when r_E2 < 0.5. It applies to all affected AlphaMissense variants (r_E2 = 0.28; 9,614 P/LP and 28,882 VUS) and REVEL variants (0.30; 7,701 and 16,453).
- VEST4 (r_E2 = 0.56) keeps enough independent evidence that 33% of its capped P/LP combinations remain at Strong (2,441 of 7,485), and 16% of its VUS.
- So the clinically useful output is the per-predictor retention share, which says how much of the cap each tool can justify.
- Under this study's own E2 thresholds, no predictor contributes PP3 at all once position and constraint are held fixed, leaving PM1 alone (2 points).

**Discrimination check (pre-specified).** Discrimination of P/LP from B/LB by PP3 points is unchanged by the conditional rule, which rescales points without reordering variants. AUC at PM1-like sites: AlphaMissense 0.911, REVEL 0.874, VEST4 0.875. So the adjustment changes how much evidence is credited, not which variants are ranked as more damaging.

### Can a model trained on the de novo signal do better?

[PENDING Aim 4 v1 and A5: calibration of predicted LRs on DDD; within-gene top-2% LR of v2-GBDT vs AlphaMissense on cold genes (primary), with the paired bootstrap; v2-MLP and ensemble; nested levels A–E for the model]

---

## Discussion

**Main finding.** Across 55 predictors, two phenotypes and two generations of tools, the same pattern holds:

- about half of a missense predictor's evidence for a de novo variant in a DD or autism proband is the identity of the gene;
- most of the remainder is the variant's position in the protein;
- once gene and position are fixed, regional constraint adds little.

After all three are held fixed, no predictor keeps calibrated Supporting-level evidence in the replication cohort or among the newest tools. This includes popEVE and AlphaMissense, whose raw LRs are the largest in the study.

**Why this happens.** It follows from how predictors are built and evaluated:

- Supervised tools learn from labels clustered in particular genes and domains.
- Unsupervised tools learn from evolutionary or population constraint, which is itself organised by gene and region.
- Both perform well in benchmarks that pool variants across genes, because between-gene differences dominate those benchmarks (cf. Harrizi et al. 2026).
- What a clinician needs is the within-gene, within-region increment. Our design measures that increment directly, and it is much smaller.

**Clinical implications.**

- **PP3 with PM1 or PP2.** When a variant already receives PM1 for location in a hotspot or critical domain, or PP2 for a missense-constrained gene, PP3 at full ClinVar-calibrated strength counts much of the same information again.
- **ClinGen's cap.** ClinGen's cap of PP3 + PM1 at Strong is directionally right. For de novo variants in DD genes our data suggest the overlap is larger than the cap allows: inside PM1-like regions, AlphaMissense's remaining evidence is below Supporting (within-region LR 2.04 for RMC o/e < 0.4).
- **What to do instead.** One option is to calibrate PP3 conditionally, within gene and region, as here. The other is to discount its points by the overlapping share when it is combined with gene- or region-level criteria.
- **Magnitude in ClinVar.** The overlap is not a corner case:
  - 30% of AlphaMissense-PP3-positive pathogenic variants sit in PM1-like regions;
  - without the cap, PP3 + PM1 alone would push 14,343 VUS to Likely Pathogenic.
  - These counts describe the risk if the cap is not applied, not observed misclassifications. How often laboratories apply the cap is not known from ClinVar.
  - Our retention estimates suggest the cap is generous for low-retention predictors: with r_E2 ≈ 0.3, the overlapping share leaves AlphaMissense or REVEL PP3 + PM1 at about 3 points, not 4. This rests on applying top-2% retention to all score intervals, which is an approximation.
  - A rule of the form "PP3 + PM1 ≤ Moderate + Supporting" (3 points) for predictors with low retention would be closer to the data than the current cap at Strong (4 points). Applied per tool through its retention share, it would credit high-retention predictors such as VEST4 more.

**Relation to previous work.**

- ClinGen calibrations (Pejaver et al. 2022; Bergquist et al. 2025) and gene- and domain-aware ClinVar calibrations (Chen et al. 2026, preprint) all estimate predictor evidence against ClinVar labels.
- A recent thesis calibrated predictors to gene-level population depletion and used DD de novos for evaluation (Pathak 2026).
- **De novo variants have been used to benchmark predictors before:**
  - AlphaMissense was evaluated on DD de novos against controls (Cheng et al. 2023);
  - popEVE built its test set from DD de novos and unaffected autism siblings (Orenbuch et al. 2025).
  - Those studies measure discrimination across genes.
- **Pooled benchmarks are dominated by gene identity.** It has been shown that gene identity dominates pooled ClinVar benchmarks (Harrizi et al. 2026).
- **What is new here** is narrower: predictor evidence expressed as calibrated likelihood ratios derived from the de novo excess without clinical labels, and measured conditionally on gene, protein position and constraint, on the scale ACMG/AMP uses.

**Estimand.** Our LR compares the score distribution of the de novo *excess* (observed minus expected, i.e. the disease-associated de novos) with that of the mutational background.

- **Not diluted.** Because the excess is used, the LR is not diluted by the ~80% of de novos that are background. A mixture-model correction would correct twice.
- **Different comparison class.** It is the evidence a score adds for a de novo missense variant in an affected proband against the background expected from mutation alone. That is not the ClinVar pathogenic-vs-benign LR:
  - the comparison class is "random new mutations", including damaging variants in genes unrelated to the phenotype;
  - ClinVar's is benign-labelled variants, which are selected to be tolerated and score low;
  - ClinVar-based LRs at high scores are therefore expected to be larger, and our absolute values conservative.
- **The robust quantity.** For this reason the retention shares, which compare a predictor with itself across levels, are the primary quantity. They replicate closely across cohorts. Comparisons with ClinGen cut-offs are secondary.

**Limitations.**

1. **Primary readout.** The primary pre-registered DD readout did not meet its threshold under D1. The conclusion rests on D2, adopted afterwards with a recorded reason, then confirmed prospectively in two independent tests.
2. **Mutation-rate model.** It has residual miscalibration of about 10%, mostly over-prediction at CpG transitions. δ is designed to absorb it, but finer-scale model error could still bias within-window estimates.
3. **Small units.** 50-residue windows contain few de novos individually. The within-unit estimator pools across units, but small-stratum estimates (some mechanism classes) are unstable.
4. **Ancestry.** The cohorts are predominantly of European ancestry.
5. **Training labels.** Labels for several tools are recorded from their papers' descriptions rather than re-derived.
6. **Top 2% only.** Retention is measured on each tool's top 2%. Applying it to all score intervals (A4.1) is an approximation.

**Conclusion.** Missense predictors are powerful at telling which genes and regions matter, and much weaker at telling which variants matter within them. Variant classification should calibrate predictor evidence conditionally on the gene- and region-level criteria it is combined with.

---

## Figure legends

**Figure 1. Study design** (`fig1_design.png`).

- **(a)** Likelihood ratios are estimated from the excess of observed de novo missense variants over the number expected from mutation alone. The expected count uses the Roulette per-site rate, scaled by each cohort's synonymous de novos, so no clinical labels are used.
  - The marginal LR (level A) compares sites across the genome.
  - Within-unit LRs (levels B–E) compare each de novo only with other sites in the same unit.
  - The retention share r is the fraction of a predictor's log-evidence (equivalently, of its ACMG points) that survives conditioning.
- **(b)** The nested units on a schematic protein:
  - B, the gene;
  - C, gnomAD regional missense constraint regions;
  - D2, fixed 50-residue windows;
  - E2, window × region.

**Figure 2. The evidence ladder** (`fig2_evidence_ladder.png`). The LR of the top 2% of scores at each conditioning level (log scale):

- AlphaMissense in the DD cohort and in the autism replication;
- popEVE in held-out DDD trios;
- MetaRNN, the predictor with the highest retention among the dbNSFP 5.4a tools.

Dashed lines mark the ClinGen LR cut-offs for PP3 Supporting (2.406) and Moderate (5.79). Every predictor starts above Supporting, and all but MetaRNN end below it. MetaRNN's E2 lower bound (2.12) is also below.

**Figure 3. Per-predictor retention in three settings** (`fig3_retention.png`).

- Each row is an informative predictor (marginal top-2% LR ≥ 2.406).
- Points show the retention share after conditioning on gene (B, grey), gene + 50-residue window (D2, colour) and window × constraint region (E2, open).
- The dashed line at 0.5 is the pre-specified criterion.
- Panels: the DD cohort with 36 dbNSFP 4.1a tools and AlphaMissense; the autism replication; the 17 dbNSFP 5.4a tools in the DD cohort.
- Two values above 1 in the autism panel (MutationAssessor, MutPred) are within-gene LRs exceeding the marginal LR.

**Figure 4. Clinical impact on ClinVar missense variants** (`fig4_clinvar_impact.png`).

- **(Left)** Mean PP3 + PM1 points for ClinVar P/LP and VUS missense variants at PM1-like sites, under three rules:
  - published ClinGen intervals with PP3 + PM1 capped at Strong;
  - conditional points (published × r_E2);
  - this study's within-unit E2 thresholds.
- **(Right)** The number of ClinVar VUS whose PP3 + PM1 points alone reach Likely Pathogenic (≥ 6) when both are applied at full strength without the cap, against under conditional points.

**Figure 5. Null-contrast models** [PENDING]. Calibration of predicted against empirical LR on held-out DDD trios, and the within-gene top-2% LR of each model against AlphaMissense on cold genes, with the paired bootstrap.

**Supplementary:** training-data ledger; per-tool thresholds; the mechanism contrast; the pre-registration timeline and deviation log.

---

## References (to be completed and checked against sources)

- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet* 2022.
- Bergquist T, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 2025.
- Kaplanis J, et al. Evidence for 28 genetic disorders discovered by combining healthcare and research data. *Nature* 2020.
- Zhou X, et al. Integrating de novo and inherited variants in 42,607 autism cases identifies mutations in new moderate-risk genes. *Nat Genet* 2022. doi:10.1038/s41588-022-01148-2.
- Cheng J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 2023.
- Orenbuch R, et al. Proteome-wide model for human disease genetics (popEVE). *Nat Genet* 2025.
- Seplyarskiy V, et al. Roulette mutation-rate model. [citation to be checked]
- Chao KR, et al. Regional missense constraint (gnomAD). [citation to be checked]
- Havrilla JM, et al. Constrained coding regions. *Nat Genet* 2019.
- Badonyi M, Marsh JA. *PLOS ONE* 2024. doi:10.1371/journal.pone.0307312.
- Liu X, et al. dbNSFP v4 / v5. [citations to be checked]
- Pathak 2026, PhD thesis, University of Edinburgh.
- Chen et al. 2026, Research Square preprint.
- Harrizi et al. 2026. [citation to be checked]
- Richards S, et al. ACMG/AMP standards. *Genet Med* 2015.
- Tavtigian SV, et al. Bayesian points framework. *Hum Mutat* 2020.
