# Rarity inflates non-coding variant benchmarks; RNA-family profiles predict pathogenic snRNA variants beyond it: a pre-registered evaluation with two blind tests

*Draft, 14 October 2026. Every number is traceable to `PREREG.md` section 14 (dated deviation log) and `results/`.*

## Abstract

**Background.** Variants in spliceosomal small nuclear RNA (snRNA) genes cause neurodevelopmental disorders (RNU4-2, RNU2-2, RNU5B-1, RNU4ATAC, RNU12), but clinical guidance published in 2026 notes that no in silico tool has been designed or calibrated for these genes. Genome-wide scores perform at chance on RNU2-2.

**Methods.** We pre-registered and timestamped an evaluation of variant-effect scores across 55 snRNA loci (17,385 possible SNVs), with patient labels for two genes sealed in an encrypted file before any model was fitted. A selection-supervised graph model (primary, v2) and a regularised stacked model (v4, frozen before unsealing) were evaluated blind. Because "benign" controls drawn from common gnomAD variants let any rarity-tracking score separate random rare variants from controls, we report each score's excess over a rarity-only null. Post hoc and pre-registered follow-up analyses tested circularity, position-only baselines, measured function (RNU4-2 saturation genome editing), a second blind label set curated from the 2025–2026 literature, generalisation to other disease RNA genes, the size of the rarity effect across ClinVar's non-coding variants, and ClinGen PP3/BP4 calibration.

**Results.** The pre-registered primary model failed its per-variant endpoints. The frozen v4 model separated RNU2-2 patient variants from population variants in the sealed test (AUROC 0.892 vs CADD 0.501; excess over its rarity null +0.108, p = 0.014) and replicated on 16 patient variants reported only in the literature (0.910; excess +0.133, p = 0.012; Δ vs CADD +0.298 [0.153, 0.426]). The signal beyond rarity was carried by parameter-free RNA-family (Rfam) profile features (excess +0.17, p = 0.001), not by CADD, GPN-MSA, ReMM or phyloP. On measured RNU4-2 function, Rfam features correlated at ρ = 0.39–0.46 versus 0.09 (CADD) and −0.02 (phyloP). A pre-registered test on four non-snRNA disease RNA genes was not met: genome-wide scores perform well there (CADD 0.86 on RMRP), so their failure is specific to the multicopy, hypermutable snRNA loci. Calibrated v4 thresholds reach PP3 Supporting within RNU2-2, RNU4-2 and RNU4ATAC, where CADD reaches none, and 37 ClinVar VUS reach PP3 under thresholds from their own gene. Across 1,787 ClinVar non-coding pathogenic SNVs genome-wide, the rarity-only null was 0.62–0.68 for CADD, GPN-MSA, phyloP447 and ReMM, so 40–58% of their above-chance AUROC against common controls reflects rarity alone (pre-registered; met).

**Conclusions.** Benchmarks that use population-common variants as negatives overstate non-coding predictors by a score-specific rarity effect, which should be reported. Measured beyond it, RNA-family profiles provide calibrated in silico evidence for spliceosomal snRNA variants where genome-wide scores do not.

## Introduction

Monoallelic variants in RNU4-2 cause ReNU syndrome, one of the most common single-gene causes of neurodevelopmental disorder identified by genome sequencing (Chen et al. 2024; Greene et al. 2024). RNU2-2, RNU5B-1 and RNU5A-1 followed (Greene 2025; Jackson 2025; Nava 2025), with recessive forms of RNU4-2 and RNU2-2 (Rius 2026; Greene 2026), and longer-known recessive disorders of the minor spliceosome (RNU4ATAC, RNU12, RNU6ATAC). Laboratories now report snRNA variants routinely, yet the 2026 classification guidance (D'Souza et al. 2026) recommends PP3/BP4 only once a tool has been calibrated for snRNAs, and notes that CADD gives benign-range scores to recurrent pathogenic RNU2-2 variants.

Two features of these genes make generic scores unreliable. First, they are multicopy: each snRNA has paralogs and hundreds of pseudogene copies, which complicate genome alignment and conservation estimates. Second, they are hypermutable, with de novo rates far above background (Seplyarskiy 2023; D'Souza 2026), so population depletion is weaker and harder to model than in protein-coding genes. RNA-family covariance models (Rfam) describe each snRNA family's sequence and structure across many species and paralogs without relying on whole-genome alignment.

We designed the study to be falsifiable. The analysis plan, operational definitions and decision gates were frozen and timestamped before analysis; patient labels for two genes were sealed before any model was fitted; every later change is listed with what had been seen when it was made. We report the failures alongside the successes.

## Methods (summary; full detail in Supplementary Methods and `PREREG.md`)

**Loci and variants.** GENCODE v49 spliceosomal snRNA genes passing coverage and paralog QC (55 loci); all possible SNVs (17,385). Comparators: CADD v1.7, phyloP447, ReMM v0.4, GPN-MSA, Evo 2 7B, ViennaRNA ΔΔG.

**Labels.** Positives: ClinVar P/LP (≥ 1 star), RNUdb curation and a rule-based literature tally. Negatives (L3): gnomAD v4.1 genome PASS variants with AC ≥ 9, minus positives.

**Sealing.** RNU2-2 and RNU5B-1 labels were fetched directly into an AES-256 file before any model was fitted, and decrypted once (logged).

**Models.** v2 (primary, pre-registered): a graph-transformer selection model trained on gnomAD counts under a Roulette mutation model, with an inheritance-aware framework score. v4 (amendment A5, designed after the v2 endpoints were seen, frozen and hashed before unsealing): L2-regularised logistic regression on within-gene percentile ranks of 13 features (CADD, phyloP447, ReMM, GPN-MSA, Evo 2, ViennaRNA ΔΔG, Rfam alternative-base rarity and reference-column frequency, pairing, intermolecular pairing, local constraint, homology transfer of RNU4-2 SGE to RNU4ATAC), trained on non-blind genes only.

**Rarity-only null.** For a gene with *n* positives, we drew *n* random SNVs that are neither positives nor L3 (rare or absent in gnomAD) 1,000 times and computed each score's AUROC against the same L3. The excess (observed − null mean) and its one-sided p-value are the effect sizes we report.

**Follow-up analyses.** Each was written as code and logged before it was run: circularity (ACMG criteria read from every ClinVar submission; leave-position-out constraint), position-only baselines, stricter negatives (AC ≥ 20, 50, 100), RNU4-2 SGE (435 SNVs), a second blind test from the literature (A8), a zero-shot test on non-snRNA RNA genes (A7), calibration (ClinGen local likelihood ratio; Pejaver et al. 2022) and genome-wide rarity inflation (A9).

## Results

### 1. The pre-registered primary model failed per variant

The v2 selection model did not meet its pre-registered RNU4-2 SGE bar (ρ 0.144; bar 0.346), nor the held-out-gene gate G3 (ΔAUROC vs CADD +0.083; bar 0.15). Adding UK Biobank counts (v3) did not help (ρ 0.142). On the blind genes, v2's framework score beat CADD (RNU2-2 0.784 vs 0.501), but its selection component was weak (0.659). We report v2 as an informative negative: population counts at hypermutable snRNA loci do not resolve individual variants.

### 2. The frozen v4 model separates blind RNU2-2 patient variants from population variants

On the sealed RNU2-2 labels (22 P/LP vs 127 L3), v4 reached AUROC 0.892 [0.817, 0.951]; CADD 0.501, ReMM 0.501, phyloP447 0.724, GPN-MSA 0.724 (Fig. 1). Internally, leave-one-paralog-cluster-out ΔAUROC vs CADD was +0.186 [0.108, 0.267].

### 3. Most of the apparent gap is rarity; the RNA-family profile carries the rest

L3 variants are common, so any score that tracks rarity or constraint separates random rare variants from them. On RNU2-2 the rarity-only null was 0.78 for v4, 0.80 for GPN-MSA, 0.69 for phyloP and 0.53 for CADD (Table 1). Measured as excess over that null:
- The Rfam features retained a large signal: alternative-base rarity +0.175 (p = 0.001) and reference-column frequency +0.172 (p = 0.002).
- v4 retained +0.108 (p = 0.014).
- CADD (−0.033), GPN-MSA (−0.072), ReMM (−0.018) and phyloP447 (+0.038, p = 0.27) retained none.

v4 did not improve on the parameter-free Rfam feature alone (+0.031, interval crossing 0). **The headline AUROC gap therefore overstates the advantage, and the useful signal is the RNA-family profile.**

**Table 1.** Blind RNU2-2: AUROC, rarity-only null, and excess.

| Score | Sealed (22 P) AUROC | Null | Excess (p) | Literature-only (16 P) AUROC | Null | Excess (p) |
|---|---|---|---|---|---|---|
| v4 (frozen) | 0.892 | 0.783 | +0.108 (0.014) | 0.910 | 0.777 | +0.133 (0.012) |
| Rfam alt rarity | 0.860 | 0.685 | +0.175 (0.001) | 0.826 | 0.681 | +0.145 (0.010) |
| Rfam ref column freq | 0.847 | 0.675 | +0.172 (0.002) | 0.819 | 0.669 | +0.150 (0.010) |
| phyloP447 | 0.724 | 0.686 | +0.038 (0.27) | 0.787 | 0.683 | +0.104 (0.054) |
| GPN-MSA | 0.724 | 0.796 | −0.072 (0.93) | 0.859 | 0.793 | +0.066 (0.13) |
| CADD | 0.501 | 0.534 | −0.033 (0.74) | 0.613 | 0.531 | +0.081 (0.09) |

### 4. The result is not explained by circularity or position

- **Positive labels.** No ClinVar submitter applied PP3 or BP4 to any blind RNU2-2 positive. 21 of 22 list their criteria, all drawn from PS2, PS4, PM1, PM2, PM3 and PP1. Restricting to de novo (PS2) positives gave 0.935.
- **Negative labels.** v4's local-constraint feature counts a variant's own position. Recomputing it without the position, under the frozen coefficients, gave 0.893.
- **Position-only baselines.** Table 2 region membership (0.740), distance to it (0.728), 5′ proximity (0.752) and position unobserved in gnomAD (0.743) were all below v4 (Δ +0.14 to +0.16, every interval excluding 0). Within the n.1–67 region (13 vs 14 variants), v4 still ranked variants (0.780 [0.582, 0.940]).
- **Stricter negatives.** At AC ≥ 50, v4 reached 0.968 and CADD 0.669, and phyloP447 rose to 0.912. Genome-wide conservation separates patient variants from very common variants, but not from moderately common ones.

### 5. A second blind test replicates the result

We pre-registered a second test (A8) before curating the 2025–2026 RNU2-2 literature (Greene 2025/2026, Jackson 2025/2026, Leitão 2026, and others). A frozen wording rule accepted 35 RNU2-2 SNVs; 16 were absent from the sealed set (14 recessive). On these, frozen v4 reached 0.910 [0.825, 0.979], with excess +0.133 (p = 0.012) and Δ vs CADD +0.298 [+0.153, +0.426]: both pre-registered conditions were met. Pooled with the sealed set (38 positives), v4 reached 0.901 (excess +0.122, p = 0.001), the Rfam features +0.17 (p = 0.001), and CADD and GPN-MSA about 0. The first test was mostly dominant variants and this one mostly recessive, so the result holds across both inheritance classes.

### 6. Measured function agrees

On RNU4-2 saturation genome editing (435 SNVs), which involves no population label, Rfam reference-column frequency correlated with functional damage at ρ = 0.459 [0.376, 0.533] and alternative-base rarity at 0.389 [0.305, 0.464]. GPN-MSA reached 0.143, CADD 0.088, ReMM 0.062 and phyloP447 −0.024. A v4 refitted without the RNU4 paralog cluster reached ρ 0.382, above the pre-registered bar of 0.346 (post hoc).

### 7. The advantage is specific to the snRNA loci (pre-registered negative)

We tested the Rfam alternative-base rarity zero-shot on four other disease RNA genes (A7; scores frozen before labels). The test was not met:
- **RMRP (28 P vs 40 N):** Rfam alt rarity 0.488; CADD 0.863; phyloP447 0.844.
- **SNORD118 + RNU7-1 (11 vs 134):** Rfam alt rarity 0.879; CADD 0.816.
- **TERC:** too few common variants to evaluate.

The frozen v4 tied CADD (Δ −0.005). Genome-wide scores work on these single-copy RNA genes, so their failure at the snRNA loci reflects those loci (multicopy, hypermutable) rather than RNA genes in general.

### 8. Genome-wide rarity inflation of non-coding benchmarks (pre-registered, A9)

Is the RNU2-2 rarity effect peculiar to snRNAs? We pre-registered a genome-wide test before fetching any score (A9).
- **Design.** 1,787 ClinVar non-coding P/LP SNVs (intronic 1,500 of 3,438, sampled; 5′UTR 148; ncRNA exon 85; 3′UTR 51) against 8,666 gnomAD PASS SNVs with AC ≥ 9 within ±100 bp. The window-matched rarity null drew, for each positive, a random non-positive, non-common SNV from its own window.
- **Result.** Pooled, the null was 0.632 for CADD, 0.676 for GPN-MSA, 0.647 for phyloP447 and 0.618 for ReMM, against raw AUROCs of 0.828, 0.805, 0.776 and 0.730 (Table 2). Rarity alone therefore accounts for 40%, 58%, 53% and 51% of each tool's above-chance AUROC.
- **By class.** The inflation was +0.10 to +0.19 in every class, and largest in ncRNA exons.
- **Rankings.** The two reorderings after correction (phyloP447 vs GPN-MSA in intronic variants, ReMM vs phyloP447 in 5′UTRs) are within overlapping intervals and are ties. The gaps did change: GPN-MSA, the most rarity-sensitive score, lost most.

**Table 2.** Rarity-only inflation of non-coding predictors (pooled, 1,787 P/LP vs 8,666 common SNVs; phyloP447 via the UCSC API).

| Tool | Raw AUROC | Rarity null | Inflation | Corrected (raw − null) |
|---|---|---|---|---|
| CADD v1.7 | 0.828 | 0.632 | +0.132 | +0.196 |
| GPN-MSA | 0.805 | 0.676 | +0.176 | +0.129 |
| ReMM v0.4 | 0.730 | 0.618 | +0.118 | +0.112 |
| phyloP447 | 0.776 | 0.647 | +0.147 | +0.128 |

The null is not pure noise: rare variants are themselves enriched for deleterious ones, so the corrected value is a conservative lower bound on pathogenicity-specific signal. What the test shows is that AUROC against common controls mixes two quantities, selection on rarity and pathogenicity, in score-specific proportions.

### 9. Inheritance mode from local constraint

The A6 rule (frozen local-constraint score and threshold) was evaluable for the first time on literature-annotated blind-gene variants:
- **RNU2-2 (4 dominant vs 31 recessive):** AUROC 0.839 [0.685, 0.968], threshold accuracy 86%.
- **RNU2-2 and RNU5B-1 together:** 0.892.

Dominant-acting variants sit in heterozygote-depleted windows and recessive-acting ones do not. The dominant class is small, so this is a secondary finding.

### 10. Calibrated evidence for clinical use

Out-of-sample v4 (frozen for the blind genes, cluster-cross-fitted elsewhere; 99 P vs 1,117 B pooled) reached AUROC 0.788 versus 0.746 for CADD. With the ClinGen local-LR method at prior 0.10:
- **Pooled thresholds:** PP3 Supporting at v4 ≥ 0.362, BP4 Supporting at ≤ −1.35. No Moderate threshold was reached.
- **Within-gene PP3 Supporting:** v4 reached it in RNU2-2 (≥ 0.216), RNU4-2 (≥ 0.189) and RNU4ATAC (≥ 0.88). CADD reached it in none of these.

Of 279 ClinVar VUS or conflicting SNVs in named snRNA genes, 37 reach PP3 Supporting under their own gene's threshold. Under pooled thresholds, 61 reach PP3 and 35 BP4.

### 11. Gene-level constraint (exploratory; replication)

A within-locus window-depletion score ranked the five NDD snRNA genes above all other loci (exploratory, post hoc). Whole-locus observed/expected ratios ranked them near the bottom, because the disease genes are mutation hotspots. Both observations replicate and quantify prior work (Chen 2024; Messaoud 2026), and we do not claim them as discoveries.

## Discussion

Three results carry the paper:
- Two blind tests on independent label sets agree that a frozen snRNA model separates RNU2-2 patient variants from population variants, with the signal beyond rarity coming from RNA-family profiles.
- Measured function agrees, with no population labels involved.
- The advantage is specific to the snRNA loci; a pre-registered test on other disease RNA genes showed genome-wide scores working there.

The likely mechanism is that whole-genome alignments and genome-wide training sets handle multicopy, hypermutable loci poorly. Rfam covariance models are built from curated RNA families and avoid both problems.

The rarity-only null matters beyond snRNAs. When negatives are common population variants, a score that tracks constraint or conservation earns a large AUROC from rarity alone. On RNU2-2 the null was 0.80 for GPN-MSA and 0.53 for CADD, so raw AUROC comparisons between such scores mislead. Genome-wide, the same null is 0.62–0.68 for CADD, GPN-MSA, phyloP447 and ReMM across ClinVar's non-coding pathogenic variants, so this is a general property of the benchmark design rather than of snRNAs. We suggest that benchmarks using population controls report each score's rarity null.

**Limitations.**
- The blind evidence comes from one gene in two label sets (38 positives), and RNU5B-1 could not be evaluated.
- v4 was designed after the primary model's results were seen, although it was frozen before unsealing and replicated on new labels.
- Calibration reaches only Supporting strength.
- The literature curation was AI-assisted under a fixed rule, and the curated table is released for audit.
- There is no new functional data, and the study has a single author.

## Data and code availability

- Code, the pre-registration, OpenTimestamps proofs and the deviation log are in the project repository (`snrna/`).
- The score release (all 17,385 SNVs, v2 and v4 scores, PP3/BP4 labels) is on Zenodo **[DOI]**.
- Curated literature labels: `labels/A8_literature.tsv`, with rejected rows and reasons in `labels/A8_literature_rejected.tsv`.

## Figures
1. Blind RNU2-2 ROC for v4, Rfam alt rarity, phyloP447 and CADD, sealed and literature sets (`fig1_blind_roc`).
2. Excess over the rarity-only null per score, sealed and literature sets; position-only baselines (`fig2_robustness`).
3. RNU4-2 SGE correlation by score, with the pre-registered bar; v2/v3/v4 (`fig7_aim2`).
4. A7 cross-family results (`fig5_a7`) and A9 inflation by tool and class (`fig8_a9`).
5. Calibration: pooled local LR with evidence bars (`fig4_calibration`).
Supplementary: ablation (`fig3_ablation`), gene-level constraint (`fig6_gene_discovery`).
