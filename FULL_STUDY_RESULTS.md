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

## Aim 3: about 30 predictors

*(pending: dbNSFP join)*

## Aim 4: null-contrast learning

*(pending)*
