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

## 3e. Amendment A4 (2026-10-05): clinical impact on ClinVar, and position-unit controls

**Status at the time of this amendment.**

- Seen: the A1 (DD), A2 (autism) and A3 (new tools, popEVE on DDD) readouts.
- The ClinVar file below is being **downloaded only**. No ClinVar variant has been matched, scored or counted.
- The AlphaFold archive (A5) is being downloaded. The parser `af_structure.py` has been run only on the first ~3 MB of the archive, as a format test.
- No InterPro domain or shifted-window unit has been built.

### A4.1 Clinical impact on ClinVar missense variants

**Data.**

- The ClinVar GRCh38 VCF current on 2026-10-05. The file name, release date and md5 are recorded in `results_fullstudy/provenance/`.
- Variants are SNVs whose `MC` field includes `missense_variant`, matched to a universe site (chromosome, position, REF, ALT) under the primary filters of §1.
- Review status must be ≥ 1 star: criteria provided, multiple submitters, expert panel, or practice guideline.
- Classes, from `CLNSIG`:
  - **P/LP:** Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic.
  - **B/LB:** Benign, Likely_benign, Benign/Likely_benign.
  - **VUS:** Uncertain_significance.
  - Conflicting and other values are excluded.

**Tools.** These are the tools with published ClinGen-calibrated thresholds whose exact score we hold, with intervals as in Bergquist et al. 2025 *Genet Med*, Table 1, which repeats Pejaver et al. 2022 for REVEL and VEST4. Each interval is given as points: +1 Supporting, +2 Moderate, +3, +4 Strong.

| Tool | +1 | +2 | +3 | +4 |
|---|---|---|---|---|
| AlphaMissense (primary; the gate's file) | [0.792, 0.905] | [0.906, 0.971] | [0.972, 0.989] | ≥ 0.990 |
| REVEL (dbNSFP raw) | [0.644, 0.772] | [0.773, 0.878] | [0.879, 0.931] | ≥ 0.932 |
| VEST4 (dbNSFP raw) | [0.764, 0.860] | [0.861, 0.908] | [0.909, 0.964] | ≥ 0.965 |

**PM1-like flag.** The §3 definition is used unchanged: an RMC region with missense o/e < 0.4, or CCR ≥ 95. A PM1-like site is assumed to receive PM1 at Moderate (+2).

**Points under three rules, for each tool:**

1. **Published.** PP3 points from the published intervals. Where the site is PM1-like, PP3 + PM1 is capped at 4 points, following ClinGen's recommendation (Pejaver et al. 2022).
2. **Conditional (primary).** Published PP3 points × r_E2(tool), where r_E2 is the tool's DD-cohort share of top-2% log-evidence retained at level E2. These are A1 for AlphaMissense and REVEL; VEST4 is taken from `a1_nested.tsv`, all genes. The result is added to PM1 (+2 if PM1-like) with no cap.
   - Rationale: points are proportional to ln LR. r_E2 is the share of a tool's log-evidence that is not explained by gene, position and constraint, so this keeps the ClinVar calibration and removes only the overlapping share.
   - **Stated limitation:** r is measured on the top 2% and is applied to every interval.
3. **De novo thresholds.** PP3 points from this study's own within-unit thresholds at level E2, all genes (`*_thresholds.tsv`). If E2 has no threshold for a strength, that strength is not available.

**Readouts.** Each is given for P/LP, VUS and B/LB separately:

1. The share of variants with published PP3 ≥ Supporting that are PM1-like.
2. For PM1-like variants, the mean PP3 + PM1 points under rules 1, 2 and 3, and the number of variants whose PP3 + PM1 points fall from ≥ 4 (rule 1) to < 4 (rule 2). Four points is Strong; on its own it reaches neither LP nor P.
3. **The number of VUS whose PP3 + PM1 total alone is ≥ 6 points (Likely Pathogenic) without the cap.** This is the double-counting risk when a lab applies both at full strength. It is reported under rule 1 without the cap, rule 1 with the cap, and rule 2.
4. **Discrimination check:** the AUC of PP3 points for P/LP against B/LB, under rules 1 and 2, within PM1-like sites and within other sites.

These are descriptive. No pass/fail threshold is set, because the clinical consequence is a magnitude rather than a hypothesis test.

### A4.2 Position-unit controls

These are new unit kinds, run through the unchanged estimators, for AlphaMissense and the informative A1 tools, on the DD cohort, all genes:

- **dom (curated domains, a direct PM1 proxy).** Gene × InterPro domain, from the dbNSFP 5.4a `Interpro_domain` field: the first listed domain at the site, with sites without a domain pooled as gene × "none". Level D_dom = within dom, and E_dom = dom × RMC.
- **win50s (shifted windows).** 50-residue windows offset by 25 residues. If D2 is a real position effect, r(win50s) ≈ r(win50).
- **rand50 (random partition control).** Within each gene, residues are assigned to ⌈L/50⌉ groups at random (seed 20261004), with no contiguity. This keeps D2's unit size but removes position. If the D2 drop is due to position, r(rand50) should be close to r(B), the gene level, not r(D2).

**Readout.** The median r over informative tools at B, D2, win50s, rand50 and D_dom, with the paired difference r(rand50) − r(D2).

## 3f. Amendment A5 (2026-10-05): Aim 4 v2, the null-contrast model with structure features

**Status at the time of this amendment.** Aim 4 v1 has **never produced a result**: two runs ended out of memory before the model was trained. No Aim 4 output exists. v1 (§7) is unchanged and is still run and reported.

**Motivation, written before any result.**

- v1 was frozen before the dbNSFP 5.4a tools existed in this study.
- None of v1's features encodes 3D protein structure.
- v2 adds new information. It does not add model capacity for its own sake.

**New features.** These are added to the v1 set.

1. **The 16 A3 tools**, excluding popEVE, which was developed using DD cohorts (§6). Each takes the median when missing, plus a missingness indicator.
2. **AlphaFold DB human proteome v6, single-fragment models only** (`af_structure.py`). For each residue:
   - pLDDT;
   - the 11-residue mean pLDDT;
   - the number of Cα atoms within 10 Å, with |i − j| ≥ 3 (burial).
3. **Mapping AlphaFold residues to analysed transcripts:**
   - Use the UniProt Swiss-Prot entry whose Ensembl cross-reference lists the transcript, ignoring versions; otherwise the entry with the same gene symbol.
   - Accept the mapping only if ≥ 95% of the transcript's analysed sites have a reference amino acid identical to the AlphaFold residue at the same position.
   - Otherwise, and for proteins split into multiple fragments, the features are missing, with an indicator.

**Models.** All use the same training data, LR conversion, gene split and evaluation as §7. Nothing is tuned on RUMC or DDD.

- **v2-GBDT (primary):** v1's exact `HistGradientBoostingClassifier` settings, with the extended features.
- **v2-MLP:** scikit-learn `MLPClassifier`, hidden layers (256, 128), ReLU, Adam, alpha 1e-4, early stopping on a 10% validation split, max_iter 200. Features are standardised on the training data. Five seeds (20261004 + k, k = 0…4), with probabilities averaged.
- **v2-ENS:** the mean of the v2-GBDT and v2-MLP log-odds.

**Evaluation.**

- All of §7's evaluations are run on DDD de novos (DDD scale factor), on all genes and on cold-test genes.
- In addition, each model's score is written as a column and run through the A1 nested levels A–E on DDD de novos only. This is the same machinery as popEVE in A3.

**Pre-specified comparisons.** These use a paired Poisson bootstrap over DDD probands (1,000 replicates, seed 20261004).

1. **Primary:** within-gene top-2% LR, v2-GBDT vs AlphaMissense, on **cold-test genes**. The model is reported as beating AlphaMissense if the 95% CI of the log ratio excludes 0.
2. The same comparison on all genes, and v2-GBDT vs v1.
3. r(D2) and E2 for each model, alongside the same values for AlphaMissense on DDD. This asks whether a model trained on the de novo contrast keeps more evidence beyond position.
4. v2-MLP and v2-ENS against v2-GBDT. These are secondary and descriptive.

## 3g. Amendment A6 (2026-10-05): uncertainty on retention, sibling negative control, models in the autism cohort

**Status at the time of this amendment.**

- **Seen:** the A1, A2 and A3 readouts, and a chr22-only code check of A4.1 (deviation log).
- **Not seen:** no Aim 4 (v1 or v2), A4.2 or genome-wide A4.1 output exists.
- **Origin:** this amendment answers an external review of the study pitch, which raised four points: no uncertainty on retention shares, a sibling control limited to the raw level, no evaluation of the models outside the DD cohort, and unchecked external facts.

### A6.1 Confidence intervals on retention shares

**Why the replicates pair across runs.** `run_scores` now saves the top-2% bootstrap replicates (`*_top2reps.npz`). The Poisson bootstrap weights are a fixed function of the proband set and the seed, so replicate *b* is the same resample of probands in every run on the same cohort.

**Reruns.** The runs that produced the A1 (`aim1tools`), A2 (`rep`) and A3 (`a1v54`) tables are re-run unchanged, to obtain replicates. Their point estimates must equal the committed tables; any difference is reported as a defect.

**Estimates.**

- For each tool and level, r_b = ln LRw_b(level) / ln LRm_b(A). Its 95% percentile interval is reported.
- The readout-1 share (informative tools with r(D2) < 0.5) also gets a bootstrap interval, computed per replicate with the informative set fixed at the point estimate.
- For A4.2, the paired difference r(rand50) − r(D2) gets an interval per tool and for the median over tools.

### A6.2 Sibling negative control at every level

**Analysis.** The A2 analysis (AlphaMissense + 36 tools, levels A–E, common site set, `rep` specification) on the unaffected siblings (Zhou et al. 2022), using the siblings' own synonymous scale factor.

**Expectation.** Sibling de novos carry essentially no disease-associated excess beyond the mutation-model offset, so no predictor should show calibrated evidence.

**Pre-specified pass criterion.** At most 2 of 37 predictors (5%) have a marginal top-2% LR whose one-sided lower bound (δ) reaches 2.406.

**Also reported:** the number of predictors reaching 2.406 at levels B and D2.

### A6.3 Aim 4 models in the autism cohort

**Analysis.** The fitted Aim 4 models are applied unchanged as columns of the column store: v1 and v2-GBDT/MLP/ENS (A5). Nothing is refitted. They are evaluated in autism probands and in siblings, with the A2 specification (levels A–E, common site set).

**Comparisons**, by paired Poisson bootstrap over autism individuals:

1. **Primary (A6):** within-gene top-2% LR, v2-GBDT vs AlphaMissense, all genes. The model is reported as better if the 95% interval of the log ratio excludes 0 in its favour.
2. The same comparison on the A5 cold-test genes.
3. r(D2) and E2 for each model, against AlphaMissense.
4. In siblings, every model's marginal top-2% lower bound (δ) should be below 2.406.

### A6.4 External facts and reporting

These involve no analysis. Before submission they are checked against primary sources and reported:

- ClinGen's recommendation on combining PP3 with PM1, quoted from Pejaver et al. 2022;
- the Zhou et al. 2022 trio counts, alongside our count of individuals with ≥ 1 usable de novo SNV and each filter step;
- informative-tool denominators next to every percentage.

## 3h. Amendment A7 (2026-10-06): Aim 4 v3, a gene-matched null-contrast model, designed after the v1/v2 results

**Status: written after the A5 results were seen.** This amendment is **not** independent of the data it follows. It was prompted by the A5 result:

- A5 primary: v2-GBDT / AlphaMissense within-gene top-2% LR on DDD cold genes = 0.69 [0.57, 0.84].
- The models had the largest marginal LRs in the study (v2-GBDT 27.3 against AlphaMissense 16.8 on DDD) and lower within-gene LRs (3.17 against 4.56).
- **Interpretation recorded before v3 is fitted.** v1/v2 were trained against an MR-weighted null drawn from all training genes, so their objective rewards between-gene information.

v3 changes only the null, so that the objective is within-gene. Its result is reported as a **post hoc, pre-registered-before-fitting** analysis. A5's confirmatory verdict stands as reported.

**Not yet seen:** A6 (models in autism, sibling control, CIs), the A5 nested levels for the models (`m4_nested`, running), and the InterPro control.

**Model v3.** Everything as v2-GBDT (§3f: features, hyperparameters, GeneDx training genes, seed, LR conversion), except the null:

- For each real GeneDx training de novo, 10 null sites are drawn **from the same transcript**, with probability ∝ MR.
- Transcripts have no other weighting; the null therefore matches the real de novos' gene distribution exactly.
- Gene-level features are then nearly constant within each matched set and carry little information. They are not removed.

**Evaluation.**

- As in A5: DDD and RUMC (all genes and cold genes), GeneDx cold genes, calibration, and top 1/2/5% marginal and within-gene, all with the paired bootstrap.
- The A1 nested levels on DDD de novos.
- In A6.3, the autism and sibling evaluation includes v3.

**Pre-specified comparisons:**

1. Within-gene top-2% LR, v3 / AlphaMissense, on DDD cold genes. "Better" if the 95% interval excludes 1 in v3's favour.
2. The same in the autism cohort, all genes (A6.3).
3. v3 / v2-GBDT within gene.
4. r(D2) for v3 against AlphaMissense on DDD.

## 3i. Amendment A8 (2026-10-06): Aim 4 v4, the gene-matched model trained on all DD centres and tested only in autism

**Status: written after every other result was seen.** This includes:

- the A5, A6 and A7 model results: v3 / AlphaMissense within gene was 0.96 [0.82, 1.14] on DDD cold genes and 1.03 [0.89, 1.19] in autism;
- the A4 and A6 results;
- the exploratory threshold sensitivity.

**This is the fourth model**, prompted by the study lead's request for one further attempt. It is **one shot**: it is fitted once, reported regardless of outcome, and no variant of it will be fitted. It is labelled post hoc. Its result cannot change any earlier verdict.

**Model v4.** Identical to v3 (§3h: features, gene-matched null of 10 MR-weighted same-transcript sites per real de novo, GBDT settings, seed, LR conversion), with two changes:

1. **Training de novos:** all three DD centres (GeneDx, RUMC, DDD), not GeneDx only.
2. **Training genes:** all genes. There is no cold split, because evaluation is in a different cohort. π uses all centres' synonymous scale factor.

**Evaluation**, in the autism cohort only (Zhou et al. 2022), which played no part in training, with the A2/A6.3 specification (common site set, levels A–E):

1. **Primary (A8):** within-gene top-2% LR, v4 / AlphaMissense, autism probands, all genes, paired Poisson bootstrap over individuals. "Better" if the 95% interval excludes 1 in v4's favour.
2. v4 / v3 within gene, in autism.
3. r(D2) and E2 for v4 against AlphaMissense, in autism.
4. Siblings: v4's marginal top-2% one-sided lower bound (δ) should be below 2.406.

**Multiplicity.** v4 is the fourth model compared with AlphaMissense in this study, and the result is reported in that light. A significant primary result would be described as an exploratory finding motivating an independent, pre-registered validation, not as a confirmed improvement.

## 3j. Amendment A9 (2026-10-06): a third, public, multi-phenotype cohort

**Status.** Written after all A1–A8 results were seen. The cohort has been assembled (variant lists only). **No analysis has been run on it:** no matching to the universe, no counts by class, no LR.

**Data** (`public3_dnv_snv_grch37.tsv`, md5 `f82386b55c643e74abae50caeb9c3fea`, built by `assemble_public3.py`):

| Source | Phenotype | Trios |
|---|---|---|
| Jin et al. 2017 *Nat Genet*, Table S9 (Springer static file, md5 `e54667ef…`) | Congenital heart disease | 2,645 |
| denovo-db v1.6.1 (non-SSC; md5 `97f8eb9c…`): Fromer 2014, Gulsuner 2013, McCarthy 2014 | Schizophrenia | 785 |
| denovo-db: Willsey 2017 | Tourette syndrome | 511 |
| denovo-db: Epi4K 2013 | Epileptic encephalopathy | 264 |
| denovo-db: Rauch 2012 | Intellectual disability | 51 |
| **Total** | | **4,256** |

**Exclusions** (decided from study metadata, before any analysis):

- Studies overlapping Kaplanis 2020 or Zhou 2022: DDD 2017, Sifrim 2016, de Ligt 2012, Lelieveld 2016, and every autism study.
- Homsy 2015 and Lifton, which are subsets of Jin 2017.
- Studies reporting < 0.15 synonymous de novos per proband, which indicates incomplete calls: Helbig 2016, McMichael 2015, Yu 2015, Halvardson 2016, Lemay 2015, Veeramah 2013.

**Processing.** As A2:

- GRCh37 → GRCh38 with the gate's chain; REF confirmed; primary filters.
- Duplicates are removed: any individual sharing ≥ 2 identical de novo SNVs with a Kaplanis proband **or a Zhou individual**.
- The scale factor comes from this cohort's own synonymous de novos.
- Poisson bootstrap over individuals, seed 20261004.

**Power check (pre-specified, as A2).** Adequately powered if the total missense O/E one-sided lower bound is > 1 **and** the AlphaMissense within-gene top-2% LR lower bound is > 1. Otherwise every readout is descriptive. Low power is expected, because most of these phenotypes have weak de novo missense excess.

**Readouts:**

1. **Main finding:** the A2 readouts (D2 primary) for AlphaMissense and the 36 dbNSFP 4.1a tools, with median r(B), r(D2) and r(E2).
2. **Models:** v3 and v4 against AlphaMissense, within gene, top 2%, paired bootstrap. The models were trained only on Kaplanis DD centres, so this cohort is independent of their training. Descriptive whatever the power, because this is the fifth model comparison.
3. **Descriptive:** the same within Jin 2017 probands flagged NDD = Yes.

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
| 2026-10-05 | Engineering, after amendment A3, no estimate affected:<br>• **dbNSFP 5.4a acquisition.** The file was downloaded in 24 byte-range pieces. One piece (`part.1`) began with 2,226,191,936 bytes from the file's start: during a single request the server ignored the byte range. This was found by the in-stream md5 and a BGZF read error, localised by comparison with server byte ranges, re-fetched (each 64 MB chunk confirmed as HTTP 206 of the right length), and re-verified.<br>• **A3 columns are extracted by streaming the pieces** (`scores54_stream`) rather than joining a 55 GB file. They are committed only if the md5 of the streamed bytes equals the published checksum.<br>• **`cmd_v54` run order** follows importance (study lead's instruction): the A3 confirmatory test, then popEVE (DDD only), then the Aim 3 table, then the DDD-only sensitivity. Every analysis is as specified. | Engineering and data integrity. |
| 2026-10-05 | Engineering: the laptop froze and had to be force-restarted (~19:50). chkdsk fixed only the X: bitmap (0 bad sectors). A full scan found the 18 A3 column files of chr22 empty; they were the last written before the freeze. Every other column-store file and the git repository (`git fsck`) were intact. chr22 was re-extracted from the trailing pieces of the already md5-verified download (`scores54_tail`; BGZF block CRCs checked; the stream began in chr19), with coverage identical to the first extraction (0.996). | Data integrity, no estimate affected. Analyses now run strictly one at a time. |
| 2026-10-05 | Aim 1 implementation detail: estimands (i)–(iv) are computed on one common site set, the missense sites with both an RMC region and a CCR value, so the four LRs are directly comparable. | Clarifies §3 ("sites with no region are excluded"); decided before Aim 1 was run. |
| 2026-10-05 | Engineering, after amendments A4/A5 were frozen (`fffc57c`), no specification changed:<br>• **Code checks.** The A4.1 ClinVar command was smoke-tested on chr22 only, to check that the code runs; the reported results come from the full run. The A5 feature builder was checked on 5,000 chr22 sites: the v1 features are an exact prefix of v2's. The domain-unit builder and the random-partition hash were checked on synthetic data. The first hash was visibly periodic along the protein; it was replaced by a splitmix64 mixer before any A4.2 run.<br>• **ClinVar release** `clinvar_20261004.vcf.gz` (md5 `c665dba2aa11f745fc38a6973af8025e`, verified against NCBI's checksum).<br>• **AlphaFold DB human proteome v6** was downloaded in 24 byte-range pieces and stream-parsed by `af_structure.py`; the archive was not kept.<br>• **Run order** (study lead's instruction: most important first, one job at a time): A3 sensitivity, then A4.1, then Aim 4 v1, then A4.2 (domain extraction, controls), then A5 (structure columns, v2 models, nested levels). | Engineering. |
| 2026-10-05 | **Count correction, no specification changed.** The A3 tool list (§3d) names 18 columns: 17 tools plus popEVE. The A3 commit message ("16 new tools incl. popEVE") and the A5 text ("the 16 A3 tools, excluding popEVE") miscount it. The analyses follow the stated lists and rules: A3 analysed all 17 non-popEVE tools in the full cohort and popEVE on DDD only, and A5 uses all 17 non-popEVE A3 tools as features. The study therefore covers 55 predictors (36 + AlphaMissense + 17 + popEVE). | Correction of a count in the text; the rules as written are unchanged. |
| 2026-10-05 | **Scheduling, to finish by 08:00 IST (study lead's deadline). No specification changed.**<br>• **Aim 4 v1 is not run as a separate job.** The `aim4v2` job refits the identical v1 model: same GeneDx real de novos, the same MR-weighted null drawn with the same random-number sequence, the same medians and the same v1 feature columns (an exact prefix of v2's), and the same hyperparameters. It evaluates v1 with every §7 evaluation: calibration and top 1/2/5% on DDD and RUMC, all genes and cold genes, plus GeneDx cold genes. v1 results are the `Aim4v1` rows of `aim4v2_*`.<br>• **A4.2 is run in two parts with identical specification:** the shifted-window and random-partition units first (`a4ctrl`), and the InterPro units (`a4ctrl_dom`) after the 1 h domain extraction, at the end of the night. | Scheduling. |
| 2026-10-05 | **Engineering: the first `aim4v2` run was killed for memory (out-of-memory killer, 7.6 GB) while writing predictions for chr7, after training had completed.** The training log is retained: 11,203 GeneDx de novos vs 115,100 null sites; v1 118 and v2 158 features; GBDT iterations v1 54 and v2 115; MLP epochs 15/15/17/14/17; π = 0.185. Prediction now reads memory-mapped columns in 100k-row slices. Training is deterministic (fixed seeds) and is re-run unchanged. The chain was reordered so the models run before A6, which needs the model columns for A6.3. The A6 run that had started (8 min into A6.2) was stopped and restarts in full after the models. | Engineering; no specification changed. |
| 2026-10-06 | Engineering, after A7 was frozen (`b6f6805`): `aim4v3` was checked end to end on a scratch copy of chr22 only. Its chr22 numbers rest on a few dozen de novos, are uninterpretable, and are not reported. Fix found by the check: `rep_report` now records a level with no top group as missing instead of failing. The night chain was rebuilt to run `aim4v3` between the A5 nested levels and A6, so A6.3 includes v3. The InterPro control (A4.2 part 2) runs last and may complete after 08:00 IST. | Engineering; no specification changed. |
| 2026-10-06 | Engineering, after A8 was frozen (`d96780e`): `aim4v4` was checked end to end on a scratch copy of chr22 (514 de novo sites). Its numbers are not interpretable and are not reported. The check showed that a cohort with no de novos in the top group writes no top-2% table, so the v4 report now records that case instead of failing. No specification changed. Separately, an exploratory threshold-sensitivity analysis (top 1/2/10%) was added before A8, rebuilt from committed coarse tables and labelled exploratory. | Engineering; no specification changed. |
