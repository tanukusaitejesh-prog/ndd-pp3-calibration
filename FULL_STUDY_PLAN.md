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

## 3b. Amendment A1 (2026-10-05, before any per-tool Aim 1 result): nested conditioning on gene, position and constraint

**Status at the time of this amendment:**

- The AlphaMissense results for levels A–C, and the exploratory pseudo-region control (D1 for AlphaMissense only), had been seen.
- **No result for any other tool at any level had been computed.**
- For AlphaMissense, D and E are therefore exploratory. For the other 36 tools, everything in this section is confirmatory.

**Nested levels.** Every one is a within-unit LR, using the §2 estimator, on Aim 1's common site set:

| Level | Held fixed | Unit |
|---|---|---|
| A | nothing | marginal LR |
| B | gene | gene |
| C | gene + regional constraint | gene × RMC region (primary C); gene × CCR class (secondary) |
| D1 | gene + generic position | gene × matched-count position chunk: each gene is cut into as many evenly spaced protein-position chunks as it has RMC regions |
| D2 | gene + generic position | gene × fixed 50-residue window, by protein position: residues 1–50, 51–100, … |
| E1 | gene + position + constraint | gene × D1 chunk × RMC region |
| E2 | gene + position + constraint | gene × D2 window × RMC region |

**Readout per tool.** The marginal and within-unit LR over the tool's top 2% of mutational opportunity (fine bins 980–999), with bootstrap 95% CIs and the δ one-sided bound. The coarse tables and window thresholds are as in §2.

**Strata:** all genes; DDG2P monoallelic.

**Decomposition of log-evidence**, all genes, per tool, with r(X) = ln LR_X / ln LR_A:

- gene share: 1 − r(B);
- position beyond gene: r(B) − r(D1);
- constraint beyond position: r(D1) − r(E1).

The same is reported with D2/E2.

**Pre-specified readouts:**

1. **Generalisation.** Restrict to tools whose marginal top-2% LR is ≥ 2.406. Attenuation is declared "systematic" if at least 75% of those tools have r(D1) < 0.5, i.e. less than half the log-evidence remains once gene and generic position are held fixed.
2. **Constraint beyond position.** The median across tools of r(D1) − r(E1), with the across-tool interquartile range, reported for each tool and as a summary.
3. **Clinical readout.** For each tool, whether the one-sided lower bound of the E1 LR in the top 2% reaches 2.406 (Supporting). This is the count of tools that keep Supporting-level evidence once gene, position and constraint are all held fixed.

**Wording.** Results describe the "disease-associated de novo excess", not "pathogenic variants". RMC/CCR regions are "regional constraint", not PM1. PM1 is discussed as the clinical interpretation (§3 PM1-like strata are retained as defined).

## 3c. Amendment A2 (2026-10-05): independent replication in autism de novo data

**Status at the time of this amendment:**

- The amendment-A1 readouts on the DD cohort had been seen. Readout 1 was 62% under D1 (not met) and 79% under D2.
- The replication data had been downloaded and **only counted** (variant classes and phenotype labels). **No likelihood ratio, excess or annotation against the universe had been computed for them.**

**Data.** Zhou et al. 2022, *Nat Genet* (doi:10.1038/s41588-022-01148-2), Supplementary Data 1 and 2, the published de novo coding variants.

- Cohorts: SPARK, SSC, ASC and MSSNG.
- Supplementary Data 3, "other NDD studies", is **not** used, since it may include DDD, GeneDx or RUMC.
- Phenotype "Affected" (autism) is the replication set. "Unaffected" (siblings) is a **negative control**. "Uncertain" is excluded.
- Only single-nucleotide variants are used. They are lifted hg19 → hg38 with the gate's UCSC chain, and must match a universe site, including the REF base.
- The primary filters are those of §1: exon-edge sites excluded, Roulette FILTER ≠ low, MR > 0.

**Overlap with the DD cohort.** An autism individual is treated as a duplicate of a Kaplanis 2020 proband, and **removed**, if the two share **≥ 2 identical de novo SNVs** (same chromosome, position and allele; GRCh37 coordinates in both sources). Sharing two identical de novos by chance is negligible. The number removed is reported.

**Estimation.**

- The same estimators, δ = 0.091, the same bins and windows, and the same units (§2, §3, §3b).
- The scale factor *c* is fitted on the replication cohort's own synonymous de novos.
- The bootstrap is a Poisson bootstrap over replication individuals, 1,000 replicates, seed 20261004.

**Pre-specified power check.** It is declared before any result.

- The replication is **adequately powered** if both of these hold:
  - (i) the one-sided 95% lower bound of the total missense O/E is > 1;
  - (ii) the AlphaMissense within-gene LR over the top 2% has a lower 95% bound > 1.
- Otherwise the replication is reported as **underpowered**, and readouts (1)–(3) below are descriptive only.

**Readouts on the replication cohort**, all genes:

1. **Generalisation (primary).** The amendment-A1 rule, with the position control **D2** (fixed 50-residue windows) as primary.
   - **Justification, recorded before any replication result:** in the DD cohort, D1 has a single chunk for every gene without an RMC region (about 13,500 genes), so it does not control position there. D2 controls position uniformly.
   - "Systematic" means at least 75% of informative tools (replication marginal top-2% LR ≥ 2.406) have r(D2) < 0.5. D1 is reported as secondary.
2. **Constraint beyond position:** r(D2) − r(E2) (median, IQR), plus the D1/E1 version.
3. **The number of tools whose E2 one-sided lower bound (δ) reaches 2.406**, and the same for E1.
4. **AlphaMissense nested LRs A–E**, top 2% and gate bins.

**Negative control (unaffected siblings).** Total missense O/E with 95% CI; it is expected to include 1. No LR is interpreted if the excess's lower 95% bound is ≤ 0.

**Unchanged.** The DD-cohort A1 verdict (D1, 62%, not met) is unchanged and is reported as pre-registered. A2 is a separate, independent test.

## 3d. Amendment A3 (2026-10-05): dbNSFP 5.4a predictors (plan §5, v2)

**Status at the time of this amendment.** dbNSFP 5.4a (GRCh38 BGZF, md5 `d2ad53e1a974208fea6050117be8fee1`) was downloading. **No 5.4a score had been read or analysed.**

**New tools.** Each is the dbNSFP 5.4a rankscore, defined as in §5:

- MetaRNN, MutPred2, gMVP, MisFit_D, MisFit_S;
- VARITY_R, VARITY_ER, and their LOO versions;
- ESM1b (converted rankscore), PHACTboost, MutFormer, MutScore, popEVE (converted rankscore), GPN-MSA;
- GERP (92 mammals), phyloP 470-way, phastCons 470-way.

Tools already analysed from dbNSFP 4.1a are not re-analysed. AlphaMissense keeps the gate's own file.

**Held-out rule (§6).**

- popEVE was developed and evaluated using developmental-disorder de novo cohorts. It is evaluated on **DDD de novos only**, with DDD's own synonymous scale factor, in both Aim 3 and the A1 levels.
- For every other new tool, the full-cohort result is primary. A DDD-only run is a sensitivity analysis.

**Analyses.** Unchanged definitions throughout:

- Aim 3 table and ranking;
- Aim 2 LoF vs non-LoF contrast;
- the amendment-A1 nested levels A–E.

**Pre-specified confirmatory readout.** This is a new set of tools, not seen before this amendment.

- The A1 generalisation rule, with **D2 as the primary position control**, using the A2 justification.
- "Systematic" means at least 75% of the informative new tools (marginal top-2% LR ≥ 2.406) have r(D2) < 0.5. D1 is secondary.
- Also reported: the number of new tools whose E2 one-sided lower bound reaches 2.406.

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
| 2026-10-05 | `validate` command: the engine's reproduction check against the gate uses the gate's unrounded δ (0.091189). All analyses use the pre-registered δ = 0.091. | Engineering. The plan's 0.091 is the gate's δ rounded. No effect on any full-study estimate. |
| 2026-10-05 | Bug fixes before the analysis they affect was run: the Badonyi & Marsh S3 labels file is read as latin-1 (it mixes encodings); the dbNSFP join keeps raw scores for 13 tools, alongside the rankscores, so thresholds can be converted to Pejaver 2022 raw-score units; Aim 4 builds features in 400k-row slices and maps de novos during the pass (memory). | Engineering; no estimate affected. |
| 2026-10-05 | **Exploratory control added after the Aim 1 result:** `aim1x`, the within-unit LR when each gene is cut into the same number of evenly spaced pseudo-regions as its RMC regions, on Aim 1's common site set. | To tell real constraint information from the generic effect of conditioning on finer positional units. Not pre-registered and labelled exploratory. The pre-registered Aim 1 readout is unchanged. |
| 2026-10-05 | Engineering, no estimate affected:<br>• The dbNSFP join pre-selects columns with `zcat \| cut`. Its output on chr13 is identical to the straightforward reader's.<br>• A column store (`colstore`) holds one array per column per chromosome, so a streaming pass reads only the columns it needs. Its content is identical to the pickle path (checked on chr21).<br>• Tools are batched 8 per pass for gene-unit runs and 4 otherwise.<br>• Aim 3 and the per-tool Aim 2 strata run in one set of passes (`aim3`).<br>• Aim 1's estimands are run for every tool (`aim1tools`), as §3 specifies. | The first Aim 3 attempt re-read 17 GB of tables per batch from the USB drive (more than 25 min per batch of 4 tools). It was stopped before any of its output had been read. |
| 2026-10-05 | **Bug fix: quantile bins with tied scores.** dbNSFP rankscores of several tools carry large tied blocks; for example, SIFT = 0 is about 30% of opportunity at one rankscore. A float32/float64 grid mismatch put such blocks in the wrong fine bin and left the top bins empty, so the top-2% readout was undefined. The fix:<br>• each site goes to bin floor(1000 · F(<score)), the opportunity strictly below its score, so a tied block lands whole in one bin;<br>• the "top 2%" readout is the smallest set of top bins holding ≥ 2% of opportunity (bins 980–999 when there are no ties), and its actual share is reported (`top_share`);<br>• windows that contain no opportunity are dropped before the threshold rule.<br>Found when SIFT's first per-tool row came back empty; the runs were stopped and restarted. AlphaMissense fine-bin results were recomputed with the fix (its ties are negligible). The gate-bin readouts (`AM_gate`) are unaffected. | Implementation of the pre-registered quantile bins; the intent of §2 is unchanged. |
| 2026-10-05 | Aim 1 implementation detail: estimands (i)–(iv) are computed on one common site set, the missense sites with both an RMC region and a CCR value, so the four LRs are directly comparable. | Clarifies §3 ("sites with no region are excluded"); decided before Aim 1 was run. |
