# Full study: pre-registered analysis plan

## 0. Status

This plan was written after the feasibility gate (GO WITH CAUTION; see GATE_RESULT.md) and **before any full-study result was computed**. It will be committed and timestamped with OpenTimestamps before any analysis below is run. After that point:

- only bug fixes are allowed, each with a dated deviation-log entry at the bottom;
- definitions, strata, windows, thresholds and model settings are fixed.

Where the study lead had not specified a choice, a default was chosen and is marked **[default]**. Defaults may be changed only by a dated amendment made before the analysis it affects is run.

**What is already known.** The gate's AlphaMissense results (all genes and strata in GATE_RESULT.md) were seen before this plan was written. Aim 1's first conditioning, on gene, is therefore a known quantity: within-gene LR 3.45 at AM ≥ 0.99. Everything else below is new.

## 1. Shared data and universe

- **Universe, de novos, rates and liftover:** identical to the gate (GATE_PLAN.md §2–3). `fullstudy.py universe` rebuilds the gate's per-site universe with the same functions, and the build stops unless the re-aggregated counts equal the gate's saved aggregates exactly.
- **Exon-edge exclusion (new; motivated by the gate's exploratory finding).** All primary analyses exclude sites within 3 bp of an internal exon edge (`near`), on both sides:
  - synonymous sites in the scale factor *c*;
  - missense sites in every expected count;
  - de novos at those sites.
  - The gate's "include" behaviour is reported as a sensitivity analysis.
- **Other site filters, as in the gate:** Roulette FILTER ≠ low, rate model MR, pooled synonymous scaling.
- **Centres:** all three centres are pooled for every calibration of a fixed, published score (Aims 1–3, 5), except as in §6 for tools tuned on these cohorts. Aim 4 follows §6.

## 2. Shared estimands and inference

- **Expected counts.** E(S) = c · Σ_{s∈S} MR_s, with c = O_syn / Σ MR over synonymous universe sites.
- **Marginal LR** of a set S of missense sites: [X_S / ΣX] / [E_S / ΣE], with X = O − E. Here ΣX and ΣE run over the missense sites in the analysis (sites scored by the tool).
- **Within-unit LR:** LR_S = X_S / (Σ_u O_u · E_uS / E_u − E_S). A *unit* is the thing held fixed:
  - the gene (canonical transcript), as in the gate;
  - or a gene × constraint region (Aim 1).
- **Uncertainty.** Poisson bootstrap over probands, 1,000 replicates, seed 20261004 (as in the gate), with percentile intervals.
- **Systematic margin.** δ = 0.091 is **fixed from the gate run**; it is not re-estimated. The PP3 bound uses O − (1+δ)E and the BP4 bound uses O − (1−δ)E, as in GATE_PLAN.md §4.
- **ACMG points.** points = ln(LR) / ln(2.406). This reproduces Pejaver et al. 2022's cut-offs: 2.406, 5.790, 33.53 and 1124 → 1, 2, 4 and 8 points.
- **Local-window thresholds (Aims 1–3).**
  - Each score is cut into 1,000 fine bins at quantiles of the score weighted by MR (mutational opportunity) over the analysed missense universe.
  - A window is 51 consecutive fine bins (≈ 5% of opportunity), centred every 5 bins.
  - The **PP3 threshold** for a level is the smallest window-centre score such that every window centred at or above it has a one-sided 95% lower bound (with δ) ≥ the cut-off. BP4 mirrors this with upper bounds and reciprocal cut-offs.
  - Thresholds are reported as the score at the window centre.
- **Coarse bins** (the gate's 13 AlphaMissense bins, or score deciles plus the top 2%/1% for other tools) are reported alongside as tables.

## 3. Aim 1 (headline): is PP3 evidence independent of gene and regional constraint?

**Question.** How much of a predictor's evidence survives once the variant's (a) gene and (b) regional constraint are held fixed?

**Constraint sources**, chosen before seeing results:

1. **Primary [default]: gnomAD v2.1.1 regional missense constraint.** Source: `gnomAD_v2.1.1_transcripts_with_rmc.tsv` and `…_without_rmc.tsv`.
   - Each RMC region is a unit; for a transcript without RMC, the whole transcript is one region.
   - Region GRCh37 start and stop coordinates are lifted to GRCh38 with the UCSC chain. Regions whose endpoints fail to lift, change chromosome or invert are dropped.
   - A universe site is assigned to the RMC region of the **same gene symbol** that contains it.
   - Sites with no region are excluded from Aim 1's within-region analyses; their count is reported.
2. **Secondary: CCR** (Havrilla et al. 2019, `ccrs.autosomes.v2.20180420.bed.gz`, GRCh37).
   - Each CCR interval is lifted to GRCh38, and each site takes the `ccr_pct` of the interval containing it.
   - Units are gene × CCR class, with classes ccr_pct < 50, 50–80, 80–90, 90–95, 95–99 and ≥ 99.

**Estimands**, for AlphaMissense (primary) and every Aim 3 tool:

- (i) the marginal LR;
- (ii) the within-gene LR;
- (iii) the within-(gene × RMC region) LR;
- (iv) the within-(gene × CCR class) LR.

Each is computed in coarse bins and local windows, in two strata: all genes, and DDG2P monoallelic confirmed/probable genes.

**Primary readout.** The within-(gene × RMC region) LR in AM 0.99–1.00, and its PP3 level with δ.

**Pre-specified interpretation:**

- If the one-sided 95% lower bound of the readout is **< 2.406**, the variant-level evidence that survives conditioning on regional constraint does not reach Supporting: PP3 largely double-counts PM1-like regional constraint.
- If it is ≥ 2.406, PP3 retains independent evidence after regional constraint.
- The ratio ln(LR_within-region) / ln(LR_marginal) is reported as "fraction of log-evidence that is variant-level beyond region".

**Also reported:** the overlap of PP3 with a PM1-like flag.

- A site is flagged PM1-like if its RMC region has missense o/e < 0.4 **[default]**, or if its CCR ≥ 95.
- For each flag value, the AM LR is reported within the stratum of sites with that value.

## 4. Aim 2: mechanism-specific calibration

**Mechanism labels.** All are external and never inferred from these de novos. They are taken in this precedence:

1. **Badonyi & Marsh 2024 (PLOS ONE, doi:10.1371/journal.pone.0307312), S3 literature-curated labels** (`dn`, `gof`, `lof`, with PubMed evidence).
2. **DDG2P 2019 "mutation consequence".**
   - "loss of function" → LoF.
   - "activating" or "dominant negative" → non-LoF.
   - Other values are unlabelled.
3. **Badonyi & Marsh 2024 S5 predictions** (`verdict` column, t50 threshold). Genes with "-" are unlabelled.

**Primary [default]:**

- The analysis covers DDG2P monoallelic confirmed/probable genes, labelled by precedence 1 → 2.
- The contrast is **LoF vs non-LoF**, where non-LoF = dn ∪ gof.
- Secondary analyses: dn vs gof vs lof separately; and all genes labelled by precedence 1 → 2 → 3.

**Estimands:** marginal and within-gene LR curves per mechanism class, for AlphaMissense and every Aim 3 tool. They are reported in LR and in ACMG points, with local-window thresholds per class.

**Test.** For each tool, the difference in within-gene log-LR between LoF and non-LoF genes in the top 2% of the score. It is tested by a bootstrap percentile interval.

## 5. Aim 3: about 30 predictors

**Tools, v1.** These are dbNSFP v4.1a (Zenodo 4323592), plus AlphaMissense from the gate's file:

- **Missense predictors:** SIFT, SIFT4G, PolyPhen-2 HDIV, PolyPhen-2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, PROVEAN, VEST4, MetaSVM, MetaLR, M-CAP, REVEL, MutPred, MVP, MPC, PrimateAI, DEOGEN2, BayesDel addAF, BayesDel noAF, ClinPred, LIST-S2, CADD, DANN, fathmm-MKL, fathmm-XF, Eigen, Eigen-PC;
- **Conservation:** GERP++ RS, phyloP 100way, phyloP 30way, phastCons 100way, SiPhy;
- **Other:** GenoCanyon, integrated fitCons.

**Score definition.** Each tool's dbNSFP *rankscore* is used.

- It is single-valued, its most damaging value across transcripts is defined by dbNSFP, and it is oriented so that higher means more damaging.
- Thresholds are also reported in raw score units where the raw score is single-valued.
- Sites a tool does not score are excluded from that tool's analysis, and coverage (the share of opportunity scored) is reported.

**v2 (amendment, once dbNSFP 5.4a access is granted):** ESM1b, EVE, VARITY, gMVP, MetaRNN, PrimateAI-3D where licensed, and popEVE.

**Training-data ledger.** One row per tool: the citation, the training labels it used (ClinVar / HGMD / population data / none), and whether those labels could include these cohorts' variants.

- Ledger entries are taken from the dbNSFP readme and each tool's paper, and marked unverified until checked against the paper.
- A tool tuned on these de novo cohorts is evaluated only on DDD (§6).

**Comparison with ClinVar calibration.** Thresholds are compared side by side with published ClinVar-based thresholds:

- Pejaver et al. 2022 for REVEL, BayesDel, VEST4, MutPred2, CADD, PrimateAI, EVE, SIFT, PolyPhen-2, FATHMM and GERP++, where reported;
- Bergquist et al. 2025 for AlphaMissense, ESM1b and VARITY.

Numbers are read from the papers. The comparison is stated as different estimands (GATE_RESULT.md).

**Ranking.** Tools are ranked by the within-gene LR's local-window PP3 level in DDG2P genes, then by the score's ACMG points at the top 1%.

## 6. Held-out centre

DDD is the held-out test centre for anything fitted (Aim 4) and for any tool tuned on these cohorts (none in v1). No model selection, feature choice or hyper-parameter uses DDD de novos.

## 7. Aim 4: null-contrast learning (exploratory-confirmatory split)

**Task.** Distinguish real missense de novos from missense de novos simulated under the mutation-rate null.

**Training set:**

- **Real:** GeneDx missense de novos in the universe (primary filters).
- **Null:** 10 × as many universe missense sites, sampled with probability ∝ MR, from the same gene split.

**Features [default]:**

- AlphaMissense;
- every Aim 3 rankscore, with each missing value set to the tool's median and given a missingness indicator;
- one-hot reference and alternative amino acids;
- RMC region o/e and CCR percentile.

**No gene identity and no mutation-rate or sequence-context features**, so the model cannot learn the null's own mutability.

**Model [default]:** scikit-learn `HistGradientBoostingClassifier`, with fixed settings:

- max_iter 300, learning_rate 0.05, max_leaf_nodes 31, l2_regularization 1.0;
- early stopping on a 10% validation split of GeneDx;
- random_state 20261004.

**No tuning on RUMC or DDD.**

**From classifier to LR:**

- The density ratio is r(x) = [p/(1−p)] · (n_null / n_real).
- With π the GeneDx pathogenic fraction (1 − E/O over the training genes), LR(x) = (r(x) − (1−π)) / π, floored at 0.

**Gene split.** Genes (canonical transcripts) are split 80/20 into train and cold-test sets by a hash of the gene symbol. v1 is **not** paralog-aware; paralog-family grouping is a v2 amendment.

**Evaluation (confirmatory)** on DDD de novos:

1. Calibration: bin DDD sites by predicted LR, and compare predicted against empirical marginal LR. The empirical LR uses DDD de novos and DDD's own scale factor, by the gate's estimator.
2. The empirical within-gene LR of the model's top 1% / 2% / 5%, compared with AlphaMissense's.
3. Steps 1 and 2 repeated on cold-test genes only.

RUMC is reported as internal validation.

## 8. Aim 5: data-derived priors

- **Pathogenic fraction** of de novo missense, π = 1 − E/O, with bootstrap 95% CIs, in these classes:
  - all genes;
  - DDG2P monoallelic;
  - s_het ≥ 0.15 and < 0.15;
  - each Aim 2 mechanism class;
  - each AlphaMissense coarse bin.
- **Prior comparison.** π is reported as the prior for the scenario "a de novo missense variant in an NDD proband in gene class G", against Pejaver's prior of 0.0441.
- **Converting a posterior.** For a classification framework that already counts de novo status (PS2/PM6), the posterior odds implied by π and the prior already used are reported, so evidence is not double-counted.
- **DistCurve.** A comparison with DistCurve is included only after its citation and estimand are verified. Until then it is omitted.

## 9. Orthogonal check: MaveDB

- Count the genes that have both a MaveDB score set (human, missense) and at least 20 missense de novos in the universe.
- If there are 10 or more such genes **[default]**, compare the de novo within-gene LR in the top 2% of assay scores with the AlphaMissense equivalent. Otherwise, report only the count.

## 10. Small genes

- Only about 65 genes have ≥ 20 missense de novos; this is to be confirmed in Aim 5's output.
- Per-gene LRs are therefore never reported alone. Gene classes (s_het, DDG2P, mechanism) are the reporting unit.
- A hierarchical model is a v2 amendment.

## Deviation log

| Date | Change | Reason |
|---|---|---|
