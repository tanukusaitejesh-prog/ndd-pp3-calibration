# Supplementary Information

*RNA-family profiles provide blind-validated computational evidence for spliceosomal snRNA variants. Every number below is from the dated deviation log (`snrna/PREREG.md` §14) or from committed result files.*

## S1. Timeline of freezes and results
See `paper/TIMELINE.md` (Supplementary Table S1). The v4 code and output hashes (amendment A5) are anchored in Bitcoin block 971362 (2026-10-14 05:58 UTC). The A6 rule is anchored in block 971363 (06:04 UTC). Both precede the blind results, which were committed at 06:44 UTC.

## S2. Primary model (v2) and its endpoints
- **Aim 2 (pre-registered):** Spearman ρ with RNU4-2 SGE function, strict hold-out of the U4 family: 0.144. The bar was 0.346 (Evo 2 SNV baseline 0.246 + 0.10).
- **G3 (held-out genes vs CADD, RNU4ATAC unit):** ΔAUROC +0.083 against a bar of 0.15.
- **v3 (exploratory; v2 retrained with gnomAD + UK Biobank counts, 6,324 vs 4,100 observed SNVs):**
  - Aim 2 ρ 0.142;
  - RNU4 cluster held out 0.252 (v2 0.179);
  - RNU4ATAC selection f 0.687 (v2 0.646), framework S 0.859 (v2 0.852).
- **v2 on the sealed genes:**
  - RNU2-2 framework S 0.784, function component F 0.851, selection f 0.659; CADD 0.501.
  - Weighted Δ vs CADD +0.267 [0.134, 0.399].
- **Pre-registered Aim 4 calibration (v2; prior 0.10) — met at Supporting:**
  - RNU2-2: F ≥ 0.402;
  - RNU4-2: S ≥ 4.79;
  - RNU4ATAC: PP3 S ≥ 0.241 and BP4 S ≤ −0.046;
  - dominant pool: f ≥ 4.16.
- **The same calibration at prior 0.0441:**
  - RNU2-2 F ≥ 0.547;
  - RNU4ATAC S ≥ 0.241;
  - dominant pool f ≥ 4.21;
  - RNU4-2: no Supporting threshold.

## S3. Robustness of the blind RNU2-2 result (post hoc, specified before running)
**Circularity.** No ClinVar submitter applied PP3/BP4 to any blind positive. 21/22 list criteria (PS2, PS4, PM1, PM2, PM3, PP1), and one (UDN) lists none. Every positive has PM1 and 17/22 have PM2; both overlap the L3 definition and the local-constraint feature. AUROC by subset:

| Subset | AUROC |
|---|---|
| Positives with listed criteria and no PP3/BP4 (21) | 0.896 |
| De novo, PS2 (9) | 0.935 |
| Without PM2 (4) | 0.947 |
| v4 with leave-own-position local constraint, frozen coefficients | 0.893 (frozen 0.892) |

**Position-only baselines.** v4 exceeds each by +0.14 to +0.16, with all intervals excluding 0.

| Baseline | AUROC |
|---|---|
| Guidance region membership | 0.740 |
| Distance to the region | 0.728 |
| 5′ proximity | 0.752 |
| Position unobserved in gnomAD | 0.743 |

Inside n.1–67 (13 vs 14 variants), v4 reached 0.780 [0.582, 0.940] and 5′ proximity 0.860, a tie.

**Stricter negatives (AC ≥ 50).** v4 0.968, phyloP447 0.912, GPN-MSA 0.893, Rfam alternative-base rarity 0.912, CADD 0.669. AC ≥ 100 gave similar results.

**Ablation (post hoc).** Every drop-one v4 model reached 0.875–0.908. The best single features were Rfam alternative-base rarity (0.860) and reference-column frequency (0.847); CADD reached 0.501. v4 − alternative-base rarity: +0.031, interval including 0.

**Multiplicity.** The sealed-set rarity-null tests were added after unsealing and span 12 scores. They are unadjusted; v4's p = 0.014 corresponds to about 0.17 after Bonferroni correction. The A8 excess test for v4 was pre-registered.

**Independence of the two blind tests.** The A8 test uses 16 new positives against 125 L3, of which 125 are shared with the sealed test's 127. It replicates on new positives in the same gene against nearly the same negatives.

**Label-rule sensitivity (A8; exploratory).** These runs used 122 L3 for all rules; extras.py removed the curator's STRICT IDs from the negatives, so its "frozen rule" row differs from the primary A8 run (125 L3; main Table 1) by up to 0.007 in excess (+0.126 vs +0.133) and 0.002 in Δ vs CADD (+0.300 vs +0.298). The audit log first described this as "≤ 0.004"; the figures here are the exact ones.

| Label rule | New positives | v4 − CADD | v4 excess over null | Rfam reference-column frequency excess |
|---|---|---|---|---|
| Frozen rule | 16 | +0.300 [0.166, 0.429] | +0.126 (p 0.009) | +0.148 (p 0.018) |
| Curator STRICT | 20 | +0.267 [0.148, 0.387] | +0.079 (p 0.069) | +0.096 (p 0.051) |
| Union | 27 | +0.237 [0.124, 0.349] | +0.068 (p 0.068) | +0.084 (p 0.049) |

## S4. Functional readouts
**RNU4-2 SGE, 435 SNVs** (Spearman with −function score; 1,000 bootstraps, seed 20261104; `sge_spearman.py` → `results/sge_spearman.tsv`; point estimates identical to the earlier runs):

| Score | ρ [95% CI] |
|---|---|
| Rfam reference-column frequency | 0.459 [0.377, 0.539] |
| Rfam alternative-base rarity | 0.389 [0.311, 0.464] |
| ViennaRNA intermolecular ΔΔG | 0.237 [0.142, 0.320] |
| GPN-MSA | 0.143 [0.053, 0.236] |
| CADD | 0.088 [−0.004, 0.178] |
| ReMM | 0.062 [−0.025, 0.152] |
| phyloP447 | −0.024 [−0.117, 0.067] |

v4 refitted without the RNU4 cluster reached 0.382, computed on CPU without the Evo 2 feature.

**RNU4ATAC splicing assay** (Benoit-Pilven et al. 2020, S3 "Results cell.assay"; 23 single SNVs; Spearman with −efficiency):

| Score | ρ | p |
|---|---|---|
| Homology h (n = 11) | 0.745 | 0.008 |
| GPN-MSA | 0.581 | 0.004 |
| CADD | 0.472 | 0.023 |
| ReMM | 0.461 | |
| phyloP447 | 0.382 | |
| Rfam reference-column frequency | 0.344 | 0.11 |
| Rfam alternative-base rarity | 0.168 | 0.44 |
| ViennaRNA | 0.147 | |

17/23 are training labels, so v4 is not evaluated on this assay.

**RNU4-2 SGE classification** (post hoc, specified before running: `sge_classes.py` committed in e2085a0 before the result in 87fc8d8; `results/sge_classes.tsv`). Class "depleted" = function score < −0.302 (the study's threshold); "strong" = < −0.9. AUROC [95% CI], 2,000 class-stratified bootstraps; paired differences use the same resamples.

| Score | Depleted, all (106 vs 329) | Depleted, critical region (31 vs 8) | Depleted, outside (75 vs 321) | Strong, all (43 vs 392) | Strong, outside (25 vs 371) |
|---|---|---|---|---|---|
| Rfam reference-column frequency | 0.802 [0.753, 0.848] | 0.681 [0.454, 0.879] | 0.774 [0.711, 0.834] | 0.865 [0.822, 0.904] | 0.872 [0.817, 0.922] |
| Rfam alternative-base rarity | 0.733 [0.681, 0.785] | 0.702 [0.468, 0.899] | 0.722 [0.662, 0.779] | 0.771 [0.717, 0.823] | 0.775 [0.708, 0.837] |
| ViennaRNA ΔΔG (intermolecular) | 0.629 [0.568, 0.690] | 0.484 [0.310, 0.669] | 0.624 [0.556, 0.690] | 0.637 [0.540, 0.733] | 0.665 [0.548, 0.782] |
| GPN-MSA | 0.604 [0.547, 0.662] | 0.478 [0.232, 0.730] | 0.572 [0.505, 0.641] | 0.617 [0.543, 0.690] | 0.607 [0.507, 0.701] |
| CADD v1.7 | 0.564 [0.503, 0.624] | 0.683 [0.464, 0.875] | 0.510 [0.435, 0.583] | 0.633 [0.555, 0.714] | 0.583 [0.474, 0.697] |
| ReMM v0.4 | 0.513 [0.451, 0.577] | 0.736 [0.526, 0.917] | 0.495 [0.413, 0.577] | 0.636 [0.566, 0.705] | 0.662 [0.547, 0.778] |
| phyloP447 | 0.489 [0.435, 0.543] | 0.500 (constant) | 0.473 [0.403, 0.542] | 0.502 [0.421, 0.584] | 0.479 [0.345, 0.610] |

Paired differences, Rfam reference-column frequency minus each genome-wide score: depleted/all +0.197 to +0.313 and depleted/outside +0.202 to +0.302 (every interval excludes 0); strong/all +0.228 to +0.363 (every interval excludes 0); inside the critical region −0.054 to +0.204 (every interval includes 0). For alternative-base rarity, depleted/all +0.129 to +0.244 (every interval excludes 0). phyloP447 is constant inside the critical region (see S5b), so its AUROC there is exactly 0.5. These analyses use the same 435 variants as the Spearman correlations above.

## S5. Cross-family test (A7; pre-registered; not met)
Scores were frozen at commit 55033c3 before labels were fetched.

| Unit | Rfam alternative-base rarity | Reference-column frequency | CADD | phyloP447 | GPN-MSA | ReMM | v4 |
|---|---|---|---|---|---|---|---|
| RMRP (28 P vs 40 N) | 0.488 [0.350, 0.636] | 0.829 | 0.863 | 0.844 | 0.795 | 0.790 | 0.813 |
| SNORD118 + RNU7-1 (11 vs 134) | 0.879 [0.793, 0.948] | 0.883 | 0.816 | 0.860 | 0.875 | 0.840 | 0.923 |

TERC had 15 P vs 3 N and was not evaluable. Weighted Δ (primary − CADD) −0.251 [−0.368, −0.133]; v4 − CADD −0.005 [−0.060, 0.050].

## S5b. Alignment depth at snRNA loci (exploratory, specified before running)
`alignment_depth.py` was committed (b768536) before the result (0b2117f); `results/alignment_depth.tsv`. Per named gene: the 95th percentile of phyloP447 and phyloP100 over positions (ceiling), the fraction of positions sharing the modal phyloP447 value (flat), and the per-position Spearman between phyloP447 and Rfam reference-column frequency (agreement). CADD AUROC is over all labels per gene (deviation log).

| Gene | phyloP447 ceiling | Flat fraction | Agreement with Rfam | phyloP100 ceiling | CADD AUROC |
|---|---|---|---|---|---|
| RNU4-2 | 0.90 | 0.59 | −0.01 | 3.49 | 0.52 |
| RNU2-2 | 4.46 | 0.07 | 0.39 | 2.62 | 0.55 |
| RNU12 | 7.71 | 0.22 | 0.36 | 9.68 | 0.45 |
| RNU6ATAC | 7.75 | 0.16 | 0.40 | 9.32 | 0.71 |
| RNU4ATAC | 7.67 | 0.19 | 0.55 | 8.65 | 0.77 |
| RNU4-1 (paralog of RNU4-2) | 8.03 | 0.34 | 0.32 | 3.61 | — |

Across the five genes with CADD estimates, agreement tracked CADD AUROC (Spearman +0.90); the ceiling did not (+0.20). At RNU4-2, phyloP447 is flat at 0.687 across n.60–82, one of the most conserved stretches of U4 snRNA. RNU12 does not fit the low-ceiling pattern. The UCSC alignment itself could not be inspected (API unreachable), so the paralog explanation is a candidate, not a finding. Five genes; descriptive only.

## S6. Rarity-only null across benchmarks (A9–A15)
**Definitions.** Population-observed = gnomAD v4.1 genome PASS AC ≥ 9. Pool = SNVs in the positives' ±100-bp windows that are neither positives nor population-observed. Inflation = null − 0.5. Rarity-matched AUROC = positives vs the pool.

**Caveats.**
- The A9 and A10 pools were not cleaned of unsampled ClinVar P/LP variants or of coding and splice SNVs inside the windows, and they are not matched for consequence. The null is therefore an upper bound on what rarity alone contributes, and the rarity-matched AUROC is a lower bound.
- For A11 and A13 the pool comes from the positives' windows, while the published negatives lie elsewhere, so their inflation combines rarity with local context.

**A9: ClinVar non-coding (1,787 P/LP vs 8,666 population-observed).**

| Score | Raw | Null | Inflation | Rarity-matched |
|---|---|---|---|---|
| CADD | 0.828 | 0.632 | +0.132 | 0.704 |
| GPN-MSA | 0.805 | 0.676 | +0.176 | 0.659 |
| phyloP447 | 0.776 | 0.647 | +0.147 | 0.655 |
| ReMM | 0.730 | 0.618 | +0.118 | 0.633 |

Per-class inflation was +0.10 to +0.19.

**A10: ClinVar missense (1,484 P/LP vs 1,592 population-observed; 1,078 ClinVar B/LB).** AlphaMissense and CADD are indistinguishable among rare variants.

| Score | Raw | Null | Rarity-matched | vs ClinVar B/LB |
|---|---|---|---|---|
| AlphaMissense | 0.920 | 0.762 | 0.694 | 0.916 |
| CADD | 0.813 | 0.654 | 0.690 | 0.857 |
| GPN-MSA | 0.870 | 0.802 | 0.555 | 0.827 |
| phyloP447 | 0.721 | 0.621 | 0.622 | 0.725 |

**A11: NCBoost (655 P vs 6,494 published negatives, MAF > 5%).** GPN-MSA falls from first (raw) to last (rarity-matched).

| Score | Raw [95% CI] | Null | Rarity-matched |
|---|---|---|---|
| CADD | 0.825 [0.805, 0.845] | 0.708 | 0.660 |
| ReMM | 0.817 [0.799, 0.837] | 0.774 | 0.616 |
| GPN-MSA | 0.837 [0.820, 0.853] | 0.812 | 0.573 |
| phyloP447 | 0.772 [0.749, 0.794] | 0.695 | 0.610 |

**A13: TraitGym Mendelian (338 vs 3,042; 97% of controls population-observed).** GPN-MSA falls from first to third. Within-match-group AUROC equals raw AUROC.

| Score | Raw | Null | Rarity-matched |
|---|---|---|---|
| CADD | 0.893 | 0.744 | 0.717 |
| ReMM | 0.897 | 0.817 | 0.660 |
| GPN-MSA | 0.905 | 0.827 | 0.663 |
| phyloP447 | 0.848 | 0.706 | 0.697 |

**A14 (curated benign).** AlphaMissense separated P/LP from curated rare ClinVar B/LB at 0.967. Curated benign sets are therefore much easier than the rare-variant population.

**A15 simulation** (grid of pathogenicity signal β, rarity tracking γ and damaging fraction π; `results/a15_simulation.tsv`):
- Raw AUROC overstates the target T (pathogenic vs rare benign), by up to +0.42 when β = 0.
- The rarity-matched AUROC M is a lower bound, and (M − π/2)/(1 − π) recovers T.
- raw − null is 0 without signal but over-corrects with signal, so it is used only as a test statistic.

**Intervals.** The window-bootstrap intervals logged on 15 October collapsed repeated negatives and did not resample published negatives (bug found on 4 November). They are too narrow, and rankings in the main text are therefore described without significance claims. Corrected intervals come from `a15_null_validity.py ci <analysis> --cluster`.

## S7. Inheritance-mode rule (A6): not informative
The frozen A6 rule (local-constraint score, threshold set on non-blind genes) gave the following on literature-annotated blind-gene P/LP variants:
- RNU2-2 (4 dominant vs 31 recessive): AUROC 0.839 [0.685, 0.968];
- both blind genes (6 vs 31): 0.892 [0.763, 0.989].

The local-constraint score counts the variant's own allele, however. In this set, 24 of 31 recessive and 0 of 6 dominant variants are observed in gnomAD. The variant's own allele count alone separates the classes (AUROC 0.887, reproduced on 4 Nov 2026), and with the own position removed the RNU2-2 AUROC is about 0.70. The rule mainly detects whether a pathogenic variant is itself seen in carriers, which is expected for recessive alleles. We do not present it as a predictor.

## S8. Calibration details
- **Method:** local LR, nearest 20%, ≥ 10 per class, 1,000 bootstraps; thresholds from the 5th/95th percentile; prior 0.10 primary, with 0.01, 0.0441, 0.2 and 0.5 as sensitivity.
- **v4 (post hoc):** out-of-sample scores, using frozen v4 for the blind genes and v4 refitted without each gene's paralog cluster elsewhere. Pooled thresholds therefore mix scores from different cross-fitted models.
- **Negatives:** L3 (population-observed), not curated benign variants. A likelihood ratio against population-observed variants can overstate evidence for a rare variant.
- **Rarity check:** `v4_rare_passrate.py`, pre-specified on 4 Nov 2026, computes for each threshold the pass rate among P/LP and among never-labelled rare SNVs. Their ratio is a lower bound on the LR against rare benign variants; the rule is "robust" if ≥ 2.08.
- **VUS:** 279 ClinVar VUS/conflicting SNVs in the named genes. With pooled v4 thresholds, 61 reach PP3 and 35 BP4; with own-gene thresholds, 37 reach PP3. CADD pooled: 73 PP3, 36 BP4.

## S9. Gene-level constraint (exploratory; replication of Chen et al. 2024 and Messaoud et al. 2026)
A within-locus 20-nt window-depletion score, with the mutation scale fitted within the locus, ranked the five NDD snRNA genes first among 55 loci with gnomAD + UK Biobank counts (AUROC 1.0, permutation p = 5 × 10⁻⁵; named loci only p = 2.5 × 10⁻⁴). Whole-locus observed/expected ratios ranked them near the bottom (AUROC 0.08–0.12), because the expressed genes are mutation hotspots. The score is post hoc, and its minimum over windows grows with locus length and observation count.

## S10. Healthy-carrier splicing readout (pre-specified; null)
1000 Genomes 30× genotypes and Geuvadis LCL junction counts (462 samples). Carriers of rare SNVs in minor-spliceosome snRNAs (n = 72) did not differ in minor-intron junction usage from non-carriers (n = 195): one-sided p = 0.29. Predicted-damaging carriers (n = 39) gave p = 0.33, and major-snRNA carriers (n = 195) p = 0.65. About 58% of samples carried a "rare" variant at these loci, which suggests unreliable short-read calls at multicopy snRNA loci.

## S11. Prospective test (A12; ongoing)
ClinVar scans (11 October and 1 November 2026) and a literature tally (Nakano 2026, Sevilla-Porras 2026, Park 2026; 29 SNVs mapped) left every pre-registered unit below its minimum. Remaining candidates: RNU2-2 n.155C>G and n.104T>G, and RNU5B-1 n.37G>C, pending the sealed-set check. Rius et al. 2026 (Sup. Table 7) lists 9 RNU4-2 P/LP SNVs, but all 9 were training labels before the freeze, so none qualifies (4 Nov 2026).

| Unit | Role | New positives | Minimum |
|---|---|---|---|
| U1 RNU6 adRP | primary | 0 | 8 |
| U2 U5 (RNU5A-1, RNU5B-1) | primary | 0 | 5 |
| U3 RNU2-2 (beyond sealed and A8) | primary | 0 | 8 |
| U4 other untrained loci | secondary | 4 (RNU4-1; flagged: the source calls interpretation limited) | 5 |
| U5 trained genes (new variants) | secondary | 4 (RNU4-2 ×3, RNU4ATAC ×1) | 5 |

Scans continue monthly. Each unit is evaluated once, when it reaches its minimum.
