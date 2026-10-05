# Full study: results

Plan: `FULL_STUDY_PLAN.md`, frozen at commit `6e2ee03` (2026-10-04 23:04 UTC) before any full-study result; the manifest is OpenTimestamped in `prereg_fullstudy/`. The code is `fullstudy.py`, which imports the frozen `gate.py`. Tables are in `results_fullstudy/`.

**Analysis set.** These are the gate's sites and de novos, with sites within 3 bp of an exon edge excluded on both sides (pre-registered fix for the gate's splice-edge artefact). That leaves 23,851 missense and 8,222 synonymous de novos. δ = 0.091, the bootstrap is 1,000 Poisson replicates over probands, and the seed is 20261004.

**Correctness checks:**

- The per-site universe re-aggregates to the gate's saved aggregates exactly, on all 22 chromosomes.
- The new LR engine, run with the gate's filters, reproduces every O, E, LR, CI and δ-bound in the gate's `lr_bins.tsv` to a relative difference of ≤ 1e-11, in all four strata.
- The dbNSFP join's fast path gives output identical to a straightforward reference implementation (checked on chr13).

## Aim 1 (headline): does AlphaMissense's evidence survive conditioning on regional constraint?

**Pre-registered primary readout.** The within-(gene × gnomAD RMC region) LR at AlphaMissense 0.99–1.00.

All rows below use the common site set: missense sites with both an RMC region and a CCR value.

| AM 0.99–1.00, held fixed | All genes: LR [95% CI] | One-sided lower bound (δ) | Level | DDG2P monoallelic: LR [95% CI] | Level |
|---|---|---|---|---|---|
| nothing (marginal) | 9.04 [8.26, 9.86] | 8.09 | Moderate | 5.75 [5.53, 5.98] | Supporting |
| gene | 3.37 [3.21, 3.54] | 3.13 | Supporting | 3.74 [3.60, 3.88] | Supporting |
| **gene × RMC region (primary)** | **2.25 [2.16, 2.35]** | **2.10** | **not reached** | 2.26 [2.19, 2.34] | not reached |
| gene × CCR class | 2.30 [2.20, 2.39] | 2.14 | not reached | 2.34 [2.25, 2.43] | not reached |

**Pre-specified reading.** The lower bound (2.10) is below 2.406. With regional constraint held fixed, AlphaMissense's remaining variant-level evidence does not reach PP3_Supporting.

**PM1-like strata** (within-region LR at AM ≥ 0.99):

| Stratum | Within-region LR [95% CI] | Level |
|---|---|---|
| Inside constrained regions, RMC o/e < 0.4 | 2.04 [1.98, 2.12] | not reached |
| Inside constrained regions, CCR ≥ 95 | 1.84 [1.76, 1.91] | not reached |
| Outside constrained regions, by RMC | 2.91 [2.62, 3.24] | Supporting |
| Outside constrained regions, by CCR | 2.97 [2.76, 3.18] | Supporting |

**Share of log-evidence that is variant-level beyond region:** 0.37 for all genes and 0.47 for DDG2P genes, against 0.55 and 0.75 beyond gene.

**Local-window thresholds (AlphaMissense):**

- PP3_Supporting within gene: ≥ 0.9936.
- Within gene × region: no window reaches Supporting.
- Marginal: Supporting ≥ 0.978, Moderate ≥ 0.995.

### Exploratory control (added after the result; deviation log)

Holding any finer positional unit fixed might shrink the LR on its own. The control cuts each gene into the same number of **evenly spaced** pseudo-regions as its RMC regions.

| AM 0.99–1.00, within-unit | All genes | DDG2P |
|---|---|---|
| gene | 3.37 [3.21, 3.54] | 3.74 [3.60, 3.88] |
| evenly spaced pseudo-regions | 2.50 [2.40, 2.60] | 2.57 [2.47, 2.66] |
| RMC regions | 2.25 [2.16, 2.35] | 2.26 [2.19, 2.34] |

**Interpretation:**

- Most of the drop from gene to region happens with arbitrary positional partitions as well: 3.37 → 2.50, against 2.50 → 2.25 for the constraint-specific step.
- So AlphaMissense's within-gene evidence for de novo missense is largely positional: it encodes *where in the protein* a variant falls. That is the kind of information PM1 uses, which supports Pejaver et al. 2022's recommendation to cap PP3 + PM1 ("limit the sum of the evidence strength of PP3 and PM1 to strong").
- The data suggest the overlap is larger than a cap at Strong allows for de novo NDD variants: inside PM1-like regions, AlphaMissense adds less than Supporting.
- It would be an overstatement to say AlphaMissense "re-encodes constraint scores". The constraint-specific part is modest.

## Aim 2: mechanism-specific calibration (AlphaMissense; per-tool results in Aim 3)

Labels come from Badonyi & Marsh 2024 S3 (literature), then DDG2P. Within-gene LRs, DDG2P monoallelic genes:

| Mechanism | Genes | AM 0.95–0.98 | AM 0.98–0.99 | AM 0.99–1.00 | Top 2% (window) |
|---|---|---|---|---|---|
| LoF | 176 | 1.36 | 1.85 | 4.84 [4.48, 5.23] | 5.82 [5.27, 6.43] |
| non-LoF (dn ∪ gof) | 125 | 1.59 | 1.67 | 4.04 [3.76, 4.32] | 5.21 [4.76, 5.74] |
| dominant-negative | — | 0.71 | 0.92 | 4.72 [4.42, 5.04] | 6.25 [5.71, 6.92] |
| gain-of-function | — | 3.30 | 3.09 | 2.77 [2.24, 3.30] | 3.18 [2.38, 4.04] |

**Pre-registered test.** LoF vs non-LoF, within-gene log-LR in the top 2%: Δ = 0.11 [−0.02, 0.24]. **No detectable difference.** All genes (precedence 1 → 2 → 3): Δ = −0.10 [−0.23, 0.05].

**Secondary.** In gain-of-function genes, AlphaMissense's evidence is flatter: about 3 across 0.95–1.00, with a top-2% LR of 3.18 against 6.25 for dominant-negative genes. The pathogenic de novos there are less concentrated at the very top of the score.

**Caveat.** In small mechanism strata, some windows have zero or very few observed de novos, and the proband bootstrap cannot represent uncertainty at zero counts. The BP4 thresholds reported for GOF/DN genes (as low as 0.06) are artefacts of tiny counts, not evidence. This is a limit of the pre-registered estimator in small strata; a Poisson or hierarchical treatment is the fix for v2.

## Aim 5: data-derived priors

The pathogenic fraction of de novo missense is π = 1 − E/O, with 95% CIs. The Pejaver 2022 prior is 0.0441 for comparison.

| Class | π [95% CI] |
|---|---|
| all genes | 0.200 [0.181, 0.220] |
| DDG2P monoallelic | 0.729 [0.718, 0.739] |
| s_het ≥ 0.15 | 0.469 [0.449, 0.490] |
| s_het < 0.15 | 0.163 [0.143, 0.185] |
| DDG2P LoF | 0.631 [0.605, 0.656] |
| DDG2P non-LoF | 0.753 [0.736, 0.769] |
| all genes, AM 0.99–1.00 | 0.726 [0.715, 0.737] |
| DDG2P, AM 0.99–1.00 | 0.941 [0.938, 0.944] |

- With the exon-edge exclusion the all-genes π drops from the gate's 0.221 to 0.200, as the splice-edge artefact predicted: *c* is no longer biased low.
- These are priors for "a de novo missense variant in an NDD proband". They already include de novo status, so a classification that also applies PS2/PM6 must not use them as the baseline prior.
- The DistCurve comparison is omitted: the citation is not verified (plan §8).
- **Small genes:** only **48** genes have ≥ 20 missense de novos after filtering, not "about 65"; 188 have ≥ 10. Gene classes are the reporting unit (plan §10).

## Orthogonal check: MaveDB (plan §9)

- MaveDB maps 557 score sets to 394 human genes.
- Genes with ≥ 20 missense de novos that also have a MaveDB score set: **7** (BRAF, CHD7, KCNQ2, MAP2K1, NR2F1, SATB2, SPTAN1), with 16 at ≥ 10 de novos.
- That is below the pre-registered 10, so only the count is reported.

## Aim 3: 36 dbNSFP 4.1a predictors + AlphaMissense

**Readout.** The LR over each tool's top 2% of mutational opportunity, from `aim3_top2pct.tsv`.

Some tools have large tied blocks at their top score, so their "top" group is bigger than 2%. The share is shown in the table, and those tools are not directly comparable with the rest:

| Tool | Top-group share |
|---|---|
| SIFT | 17% |
| SIFT4G | 14% |
| PolyPhen-2 HDIV | 24% |
| LRT | 31% |
| MutationTaster | 35% |
| phastCons | 57% |

**Ranking.** Pre-registered order: DDG2P within-gene window PP3 level, then ACMG points of the within-gene LR in the top 2%.

| Tool | All genes: marginal LR | All genes: within-gene LR | DDG2P: within-gene LR [one-sided lower bound, δ] | DDG2P within-gene window level |
|---|---|---|---|---|
| VEST4 | 6.55 | 4.59 | 6.15 [5.63] | Supporting (≥ rankscore 0.926) |
| MutPred (v1) | 5.74 | 4.36 | 4.98 [4.51] | Supporting |
| AlphaMissense | 15.42 | 4.29 | 4.87 [4.64] | Supporting |
| MutationAssessor | 4.98 | 3.70 | 4.60 [4.15] | Supporting |
| ClinPred | 5.84 | 3.73 | 4.60 [4.22] | Supporting |
| REVEL | 9.03 | 3.63 | 4.31 [4.03] | Supporting |
| PROVEAN | 3.49 | 3.20 | 4.24 [3.75] | Supporting |
| BayesDel noAF / addAF | 8.13 / 8.94 | 3.28 / 3.36 | 4.00 / 3.98 | Supporting |
| Eigen, Eigen-PC, PolyPhen-2 HVAR, LIST-S2, SIFT4G, SIFT | — | 2.3–3.5 | 2.5–3.7 | Supporting |
| MetaLR, MVP, DEOGEN2, MetaSVM, M-CAP | 6.5–11.8 | 2.5–3.0 | 2.7–3.2 | not reached (in windows) |
| CADD, PrimateAI, fathmm-XF, FATHMM, **MPC** | 3.7–14.2 | 1.9–2.2 | 2.2–2.6 | not reached |
| conservation-only and other: phyloP, GERP++, phastCons, SiPhy, DANN, fathmm-MKL, GenoCanyon, fitCons, LRT, MutationTaster | 0.9–4.7 | 0.8–1.9 | 0.8–2.0 | not reached |

**Main points:**

- **No predictor reaches PP3_Moderate once the gene is known.** Using local windows within gene in DDG2P genes, 15 of 37 tools reach Supporting, and none reaches Moderate. VEST4 comes closest: its top-2% lower bound is 5.63, against the 5.79 cut-off.
- **A large marginal LR is mostly gene information, and the size of the drop depends on how the tool was built:**
  - **MPC**, which is built from regional constraint: 14.2 → 2.0.
  - **DEOGEN2**, which includes gene-level features: 11.8 → 2.8.
  - **AlphaMissense:** 15.4 → 4.3.
  - **REVEL:** 9.0 → 3.6.
  - The pure constraint score losing the most is what the framework predicts, so it serves as an internal check.
- **Conservation-only scores add almost nothing within gene:** LRs of about 1–2.

**Thresholds against ClinVar-based calibration.** Pejaver et al. 2022, Table 2 (read from the paper) against our marginal, all-genes thresholds, in raw-score units (`aim3_thresholds_long.tsv`):

| Tool | Pejaver PP3 Supporting / Moderate / Strong | De novo marginal PP3 | De novo DDG2P within-gene PP3 |
|---|---|---|---|
| REVEL | ≥ 0.644 / ≥ 0.773 / ≥ 0.932 | Supporting ≥ 0.773 | Supporting ≥ 0.876 |
| BayesDel addAF | ≥ 0.13 / ≥ 0.27 / ≥ 0.50 | Supporting ≥ 0.330 | Supporting ≥ 0.436 |
| VEST4 | ≥ 0.764 / ≥ 0.861 / ≥ 0.965 | Supporting ≥ 0.932 | Supporting ≥ 0.923 |
| CADD | ≥ 25.3 / ≥ 28.1 / — | Supporting ≥ 29.2 | not reached |
| PrimateAI | ≥ 0.790 / ≥ 0.867 / — | Supporting ≥ 0.858 | not reached |
| MPC | ≥ 1.360 / ≥ 1.828 / — | Supporting ≥ 1.945, Moderate ≥ 2.288 | not reached |
| PolyPhen-2 (HVAR) | ≥ 0.978 / ≥ 0.999 / — | Supporting ≥ 0.999 | Supporting ≥ 0.999 |
| phyloP (100-way) | ≥ 7.367 / ≥ 9.741 / — | Supporting ≥ 9.471 | not reached |

For every tool, the de novo calibration is stricter than the ClinVar-based one, and markedly so within gene. The same caveat applies as for AlphaMissense: the two calibrations measure different estimands (GATE_RESULT.md).

**Comparability caveats:**

- Our MutPred is MutPred v1 (dbNSFP); Pejaver calibrated MutPred2.
- Our PolyPhen-2 row is HVAR; check which version Pejaver used before claiming a difference.
- The phyloP alignment set must be matched to Pejaver's before any claim is made.
- Thresholds are window centres, at the resolution of 1,000 opportunity quantiles.

**Aim 2 per tool.** LoF vs non-LoF, DDG2P genes, within-gene log-LR in the top 2% (`aim3_contrast.tsv`):

- **Nominal:** 20 of 37 tools have a 95% CI excluding 0. Of these, 17 give more evidence in LoF genes and 3 in non-LoF genes.
- **After Benjamini–Hochberg correction across 37 tools:** 17 have q < 0.05. **15 give more within-gene evidence in LoF genes**; 2 (phyloP, VEST4) give more in non-LoF genes. The p-values are approximated from the bootstrap CI by a normal approximation (`aim3_contrast_BH.tsv`).
- Examples, Δ log-LR [95% CI]:
  - PROVEAN 0.89 [0.51, 1.36]
  - MutationAssessor 0.59
  - MetaLR 0.58
  - MetaSVM 0.57
  - BayesDel noAF 0.56
  - REVEL 0.45 [0.21, 0.72]
- **AlphaMissense shows no mechanism dependence:** 0.11 [−0.02, 0.24], q = 0.16.
- The reverse direction: VEST4 −0.36, phyloP −0.57.

**Training-data ledger** (`aim3_training_ledger.tsv`). The training-label column is **unverified** until it is checked against each tool's paper. Several top-ranked tools (ClinPred, REVEL, VEST4, MutPred) were trained on ClinVar/HGMD-type labels that could include variants from these cohorts.

**Not yet included** (plan §5 v2): ESM1b, EVE, VARITY, gMVP, MetaRNN (dbNSFP 5.4a, access requested 2026-10-05); popEVE (to be evaluated on DDD only); PrimateAI-3D (licensed).

## Amendment A2: independent replication in autism de novos (Zhou et al. 2022)

Pre-registered in plan §3c (commit `c5b94e5`, OpenTimestamped) before any replication LR was computed. Outputs are in `results_fullstudy/replication_affected/` and `replication_unaffected/`.

**Data.**

- 21,603 SNVs from autism probands (SPARK, SSC, ASC, MSSNG).
- **64 individuals removed** because they share ≥ 2 identical de novos with a Kaplanis 2020 proband.
- After liftover and matching to the universe: 12,805 missense and 4,960 synonymous de novos in 12,316 individuals. Of these, 11,809 missense de novos pass the §1 filters.

**Power check (pre-specified): adequately powered.**

- The total missense O/E is 1.126 [1.089, 1.166], with a one-sided lower bound of 1.095 (> 1 required).
- The AlphaMissense within-gene top-2% LR has a lower bound of 2.55 (> 1 required).

**Readout 1 (primary, D2): SYSTEMATIC.**

- 16 of 20 informative tools (80%) keep < 50% of their log-evidence after gene and position. The pre-specified threshold is 75%.
- D1 (secondary) gives 75%.
- For comparison, the DD cohort gave 79% under D2. Its pre-registered D1 verdict of 62% (not met) stands as reported.

**Readout 2: constraint beyond position is small.** The median r(D2) − r(E2) is 0.017 [IQR 0.003, 0.026], and r(D1) − r(E1) is 0.065. Once position is controlled, regional constraint adds little.

**Readout 3: no tool keeps Supporting.** 0 of 36 tools have an E2 one-sided lower bound ≥ 2.406, and 0 under E1. The highest E2 point estimates are:

| Tool | E2 LR | Lower bound (δ) |
|---|---|---|
| MutPred | 2.17 | 0.89 |
| ClinPred | 2.12 | 1.23 |
| MutationAssessor | 2.10 | 0.92 |
| AlphaMissense | 1.90 | 1.50 |

**AlphaMissense, nested levels (top 2%), autism replication vs DD:**

| Level | A raw | B gene | C RMC | D1 | D2 | E1 | E2 |
|---|---|---|---|---|---|---|---|
| Autism | 8.67 | 3.07 | 2.18 | 2.33 | 1.93 | 2.03 | 1.90 |
| DD | 13.43 | 4.10 | 2.57 | 2.91 | 2.19 | 2.45 | 2.08 |

**Consistency across cohorts.** These are descriptive, not pre-specified.

- The tool ranking by raw LR agrees: Spearman ρ = 0.87 over 37 tools.
- The share of evidence that survives position agrees: r(D2) has ρ = 0.76 over the 20 tools informative in both cohorts.
- The tools that keep ≥ 50% of their evidence under D2 are almost the same in both cohorts:
  - autism: MutationAssessor, VEST4, MutPred, ClinPred;
  - DD: the same four, plus PROVEAN.

**Negative control (unaffected siblings).**

- 4,246 missense de novos in 4,051 siblings. The total missense O/E is 1.106 [1.045, 1.175], so the CI excludes 1.
- This is the same roughly 10% offset seen in the gate's negative control: the Roulette expectation under-predicts missense relative to synonymous by a constant factor. The δ margin is designed to absorb this.
- What matters: among AlphaMissense 0.99–1.00 sites, the siblings' O/E is **0.97 [0.81, 1.14]**, compared with **1.86 [1.73, 1.99]** in autism probands. The high-score excess is specific to affected individuals.

**Interpretation.**

- In an independent cohort with a different phenotype, most predictor evidence is again explained by gene and position.
- Constraint beyond position is again small.
- No tool keeps calibrated Supporting evidence once gene, position and constraint are all accounted for.

## Related work, verified from source (2026-10-05)

**Pathak, *Harnessing the human mutational spectrum to advance missense variant effect prediction*, PhD thesis, University of Edinburgh, February 2026** (ERA).

- It calibrates VEP scores to gene-level **population-depletion** observed/expected, using gnomAD variation with Roulette expectations.
- The 31,058 DD de novos (Kaplanis 2020) and 9,060 control trios are used to **evaluate** case/control enrichment. It also estimates the de novo missense excess in diagnosed cases, about 56% causal, for interpreting recall.
- **Its ACMG likelihood ratios are derived from ClinVar** P/LP (autosomal dominant genes) against B/LB, following Tavtigian.
- There is no PM1 analysis, and no conditioning of predictor evidence on gene, position or constraint. The word PM1 does not occur in the thesis.

**Chen et al. 2026, *Gene- and domain-aware calibration increases the clinical utility of variant effect predictors*, Research Square preprint (8 May 2026; Fowler and Pejaver among the authors).**

- Gene-specific calibration for 132 genes, plus domain-pooled calibration for 2,637 genes, for REVEL, MutPred2 and AlphaMissense.
- Built on **ClinVar** controls and gnomAD, with a DistCurve-derived prior.
- No analysis of PP3/PM1 overlap.

**How this study differs from both:**

- Its likelihood ratios come from the **de novo excess over a mutation-rate null, with no ClinVar labels**.
- It **measures how much predictor evidence remains after conditioning** on gene, protein position and regional constraint.

Both earlier works calibrate *with* ClinVar, and neither quantifies conditional or incremental evidence.

## Aim 4: null-contrast learning

*(pending)*
