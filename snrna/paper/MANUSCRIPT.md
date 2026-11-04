# RNA-family profiles provide blind-validated, calibrated computational evidence for spliceosomal snRNA variants: a pre-registered study with a rarity-adjusted benchmark

**Running title:** Calibrated in silico evidence for snRNA variants

**Author:** [Author name]¹ — ¹[Affiliation]. Correspondence: [email]

*Manuscript draft, version of 5 November 2026. Every number is traceable to the dated deviation log (`snrna/PREREG.md` §14), the result files in `snrna/results/`, and the verification report `paper/VERIFICATION_2026-11-05.md`.*

---

## Abstract

**Purpose.** Variants in spliceosomal small nuclear RNA (snRNA) genes are a frequent cause of neurodevelopmental disorders, but no computational tool has been calibrated for these genes. Without one, laboratories cannot apply the ACMG/AMP computational criteria (PP3/BP4).

**Methods.** We pre-registered and cryptographically timestamped an evaluation of variant-effect scores at 55 snRNA loci (17,385 possible single-nucleotide variants). Patient labels for RNU2-2 and RNU5B-1 were sealed before any model was fitted. A model was frozen and anchored in the Bitcoin blockchain before the labels were opened. We scored performance as the excess over a score-specific rarity-only null: random unlabelled rare variants scored against the same population-observed negatives.

**Results.** The pre-registered primary model, trained on population allele counts, failed its per-variant endpoints. On the sealed RNU2-2 labels (22 pathogenic vs 127 population-observed variants), the frozen model reached an AUROC of 0.892 versus 0.501 for CADD. It replicated on 16 patient variants reported only in the literature (0.910; Δ vs CADD +0.298 [0.153, 0.426]). Beyond rarity, the signal came from RNA-family (Rfam) alignment profiles (excess +0.17, p = 0.001) and not from CADD, ReMM, GPN-MSA or phyloP. The same profiles tracked measured RNU4-2 function (Spearman ρ 0.46 vs 0.09 for CADD). The advantage did not extend to other disease RNA genes (pre-registered negative). Calibrated thresholds gave PP3 Supporting evidence within RNU2-2, RNU4-2 and RNU4ATAC. Across 1,787 ClinVar non-coding pathogenic variants, rarity alone gave genome-wide scores AUROCs of 0.62–0.68.

**Conclusion.** RNA-family profiles provide blind-validated computational evidence for spliceosomal snRNA variants where genome-wide scores fail. Benchmarks with population-derived negatives should report each score's rarity-only null.

**Keywords:** snRNA; RNU2-2; RNU4-2; variant interpretation; PP3/BP4 calibration; benchmark bias; pre-registration

---

## Introduction

De novo variants in RNU4-2, which encodes the U4 small nuclear RNA (snRNA) of the major spliceosome, cause ReNU syndrome, one of the most frequent monogenic neurodevelopmental disorders (NDDs) identified by genome sequencing.¹,² Variants in other snRNA genes followed: RNU2-2 (dominant and recessive), RNU5B-1, recessive forms of RNU4-2,³⁻⁷ and the longer-known recessive disorders of the minor spliceosome (RNU4ATAC, RNU12, RNU6ATAC).⁸⁻¹¹ Variants in RNU4-2 and RNU6 genes also cause autosomal dominant retinitis pigmentosa.¹² Diagnostic laboratories now encounter snRNA variants routinely.

Interpretation is limited by the lack of computational evidence. Proposed snRNA classification guidance recommends computational evidence (PP3/BP4) only after a tool has been calibrated for snRNA genes. It notes that genome-wide scores such as CADD give benign-range values to recurrent pathogenic RNU2-2 variants.¹³ Saturation genome editing of RNU4-2 now provides functional measurements for one gene,¹⁴ but no comparable data exist for the others. snRNA genes are unusual targets for genome-wide predictors. They occur in multiple paralogous copies with hundreds of pseudogenes, and they are germline mutation hotspots.¹⁵,¹³ RNA-family covariance models (Rfam) instead describe each snRNA family's sequence and secondary structure across species and paralogs, independent of whole-genome alignment.¹⁶,¹⁷

A second obstacle is evaluation. Benchmarks often use population variants as negatives. Because rarity correlates with conservation and constraint, such negatives can flatter a predictor. This has been described for missense predictors and benign sets of common variants,¹⁸⁻²¹ and ClinGen calibration therefore restricts benign variants to allele frequency below 1%.²² In small genes, negatives are usually defined by population observation; here, gnomAD allele count ≥ 9. We show that most such negatives are rare (median allele frequency 1.8 × 10⁻⁴), yet any score that tracks rarity still separates a random unlabelled variant from them. We therefore report each score's excess over a score-specific rarity-only null.

We designed the study to be falsifiable. The plan, definitions and decision rules were frozen and timestamped before analysis. Patient labels for two genes were sealed in an encrypted file before any model was fitted, and every later change is logged with what had been seen when it was made. We report the failed endpoints alongside the successful ones.

## Materials and Methods

### Study design and pre-registration
The pre-registration (`PREREG.md`), analysis code and each amendment were hashed and timestamped with OpenTimestamps. The proofs are now anchored in Bitcoin blocks (Supplementary Table S1). Changes after data were seen are logged in a dated deviation log that records what had been viewed at the time (`PREREG.md` §14). Analyses are labelled pre-registered (with amendment number), post hoc, or exploratory throughout.

### Loci, variants and scores
We included the GENCODE v49 spliceosomal snRNA genes and paralogs that passed coverage and paralog quality control: 55 loci, 20 of them named genes, with all 17,385 possible SNVs. Comparators were:
- CADD v1.7;²³
- phyloP from the 447-mammal alignment;²⁴
- ReMM v0.4;²⁵
- GPN-MSA;²⁶
- Evo 2 (7B) log-likelihood ratios;²⁷
- ViennaRNA 2 folding and duplex free-energy changes for intermolecular snRNA pairs.²⁸

RNA-family features came from aligning each locus to its Rfam 15.1 covariance model with Infernal 1.1.5 (cmalign).¹⁶,¹⁷ We used the frequency of the reference base in its alignment column ("reference-column frequency") and one minus the frequency of the alternative base ("alternative-base rarity"). Both are parameter-free: no label is used to compute them.

### Labels
**Positives** were pathogenic or likely pathogenic SNVs from three sources: ClinVar (≥ 1 review star),²⁹ RNUdb curation, and a literature tally under a frozen wording rule (a source must call the variant pathogenic, likely pathogenic or disease-causing; any source calling it a VUS excludes it).

**Negatives (L3)** were gnomAD v4.1 genome PASS SNVs with allele count ≥ 9, excluding positives. L3 holds 1,233 SNVs. They are population-observed but mostly rare: median allele frequency 1.8 × 10⁻⁴, with 5.4% at ≥ 1%.

**Sealing.** RNU2-2 and RNU5B-1 labels were fetched directly into an AES-256-encrypted file before any model was fitted. They were decrypted once, and each event was logged.

### Models
**v2 (primary, pre-registered).** A graph-transformer selection model trained on gnomAD allele counts under the Roulette mutation-rate model,¹⁵ with paralog clusters as the hold-out unit. It used an inheritance-aware framework score: selection for dominant genes, and lab-derived function (RNU4-2 SGE transferred to RNU4ATAC by structural alignment) for recessive genes. A version with UK Biobank allele counts added (v3) was exploratory.

**v4 (amendment A5).** Designed after the v2 endpoints were seen. It is an L2-regularised logistic regression (C = 0.1, balanced class weights, no tuning) on within-gene percentile ranks of 13 features:
- the six comparator scores;
- the two Rfam features;
- reference pairing and intermolecular-pairing loss;
- absolute folding ΔΔG;
- local constraint: the negative z-score of observed versus expected variation in a 21-nt window, with gnomAD and UK Biobank counts;
- homology-transferred RNU4-2 SGE function (RNU4ATAC only).

It was trained on non-blind genes only; the blind genes never entered training. Its code and output hashes were committed and anchored in Bitcoin block 971362 (14 October 2026, 05:58 UTC), before the blind labels were decrypted. Its blind evaluation script checks those hashes before scoring.

### Evaluation and the rarity-only null
The primary metric was the area under the ROC curve (AUROC) for positives against L3. We report 95% class-stratified bootstrap intervals (2,000 replicates).

Because L3 is defined by population observation, any score that tracks rarity separates random unlabelled variants from L3. For a set of *n* positives in a gene, the **rarity-only null** draws *n* SNVs at random (without replacement) from that gene's SNVs that are neither positives nor L3. We score each draw against the same L3 and repeat 1,000 times. We report the null mean and the **excess** (observed AUROC − null mean), with a one-sided permutation p-value, (k + 1)/(B + 1).

Because the null pool contains unlabelled damaging variants, the excess is a conservative test statistic. In a ground-truth simulation (A15; Supplement), raw AUROC overstated the true discrimination between pathogenic and rare benign variants. The rarity-matched AUROC (positives against random rare variants) was a lower bound, which a correction for the damaging fraction recovered exactly.

### Blind tests
**Test 1: sealed labels.** RNU2-2 had 22 P/LP vs 127 L3 SNVs. RNU5B-1 had 1 P/LP and was not evaluable.

**Test 2: literature (amendment A8).** The design was fixed before curation. RNU2-2 and RNU5B-1 variants reported as pathogenic in the 2025–2026 literature were curated under the frozen wording rule, by an AI-assisted search instructed not to look up any score. Positives absent from the sealed set and from any training file formed the test set, and the labels were committed before scoring. The pre-registered conditions were:
- v4 excess over its rarity null with p < 0.05;
- v4 − CADD AUROC interval excluding 0.

### Robustness analyses (post hoc, specified before running)
- **Circularity:** the ACMG/AMP criteria applied by every ClinVar submitter to each blind positive, and v4 recomputed with leave-position-out local constraint.
- **Position-only baselines:** membership of guidance-defined disease regions,¹³ distance to them, 5′ proximity, and whether the position is observed in gnomAD.
- **Stricter negatives:** L3 restricted to AC ≥ 20, 50 or 100.

### Functional readout
- **RNU4-2 SGE:** Spearman correlation with function scores for 435 RNU4-2 SNVs.¹⁴ No population label is involved.
- **RNU4ATAC splicing assay (exploratory):** the published splicing efficiencies of 23 RNU4ATAC SNVs.³⁰

### Cross-family test (amendment A7)
To test whether the Rfam advantage generalises, alternative-base rarity was scored zero-shot on RMRP, TERC, SNORD118 and RNU7-1, with scores frozen before labels were fetched. The primary feature was chosen from the post hoc RNU2-2 ablation, before any A7 data existed.

### Benchmark analyses (amendments A9–A15)
- **A9:** 1,787 ClinVar non-coding P/LP SNVs (intronic 1,500 of 3,438 sampled, 5′UTR 148, ncRNA exon 85, 3′UTR 51, upstream 3) against 8,666 gnomAD PASS SNVs with AC ≥ 9 within ±100 bp. The window-matched null drew, for each positive, a random non-positive SNV with AC < 9 from the same window.
- **A10:** the same design for missense variants with AlphaMissense.³¹
- **A11:** the NCBoost published benchmark (common dbSNP negatives, MAF > 5%).³²
- **A13:** the TraitGym Mendelian benchmark (published matched controls).³³
- **A15:** a ground-truth simulation of the null, with bootstrap intervals over windows.

### Calibration
We used the ClinGen local likelihood-ratio method,²² as pre-registered:
- nearest 20% of labelled variants, at least 10 per class;
- 1,000 bootstraps;
- thresholds from the 5th (PP3) or 95th (BP4) bootstrap percentile of the local LR;
- primary prior 0.10;³⁴
- sensitivity priors 0.01, 0.0441, 0.2 and 0.5.

The pre-registered calibration (Aim 4) used the v2 framework score. Calibration of v4 was post hoc and used out-of-sample scores: frozen v4 for the blind genes, cluster-cross-fitted v4 elsewhere. Negatives were L3. ClinVar VUS and conflicting SNVs in the named genes (279 SNVs) were scored against the thresholds.

### Statistics and software
Analyses used Python 3 (NumPy, pandas, SciPy, scikit-learn). All code, inputs and logs are in the repository.

## Results

### The pre-registered primary model failed per variant
v2 did not meet its pre-registered endpoint against RNU4-2 SGE function:
- Spearman ρ 0.144 against a bar of 0.346, which was Evo 2's ρ of 0.246 plus 0.10;
- adding UK Biobank counts (v3) did not help (ρ 0.142).

v2 also failed the held-out-gene gate: ΔAUROC vs CADD +0.083 against a bar of 0.15. Population allele counts at these hypermutable loci therefore did not resolve individual variants. On the sealed RNU2-2 labels, v2's framework score (0.784) exceeded CADD (0.501). This came from its lab-transfer function component (0.851), not its selection component (0.659).

### The frozen v4 model separates blind RNU2-2 patient variants from population variants
In leave-one-paralog-cluster-out validation on non-blind genes, v4 exceeded CADD by +0.186 [0.108, 0.267]. On the sealed RNU2-2 labels (22 P/LP vs 127 L3; Fig. 1A), v4 reached an AUROC of 0.892 [0.817, 0.951]. The comparators were:
- CADD 0.501 [0.393, 0.612];
- ReMM 0.501;
- phyloP447 0.724;
- GPN-MSA 0.724.

### Most of the apparent gap is rarity; the RNA-family profile carries the rest
Because L3 is defined by population observation, scores differ in how far rarity alone carries them. On RNU2-2, the rarity-only null was 0.783 for v4, 0.796 for GPN-MSA, 0.686 for phyloP447 and 0.534 for CADD (Table 1; Fig. 2). Measured as excess over the null:
- **Rfam features** kept a large signal: alternative-base rarity +0.175 (p = 0.001) and reference-column frequency +0.172 (p = 0.002).
- **v4** kept +0.108 (p = 0.014).
- **Genome-wide scores** kept none: phyloP447 +0.038 (p = 0.27), CADD −0.033, ReMM −0.018 and GPN-MSA −0.072.

v4 did not improve on alternative-base rarity alone (+0.031; interval including 0). The headline gap between v4 and CADD therefore overstates the advantage: the useful signal is the RNA-family profile. We note that the two Rfam features were singled out by a post hoc ablation on the sealed labels. Their claim to replicate rests on the second blind test, the functional data and the cross-family test below.

**Table 1.** Blind RNU2-2: AUROC, rarity-only null (1,000 draws) and excess, for the sealed labels and the literature-only (A8) labels.

| Score | Sealed (22 P vs 127 L3) AUROC | Null | Excess (p) | A8 literature-only (16 P vs 125 L3) AUROC | Null | Excess (p) |
|---|---|---|---|---|---|---|
| v4 (frozen) | 0.892 | 0.783 | +0.108 (0.014) | 0.910 | 0.777 | +0.133 (0.012) |
| Rfam alternative-base rarity | 0.860 | 0.685 | +0.175 (0.001) | 0.826 | 0.681 | +0.145 (0.010) |
| Rfam reference-column frequency | 0.847 | 0.675 | +0.172 (0.002) | 0.819 | 0.669 | +0.150 (0.010) |
| phyloP447 | 0.724 | 0.686 | +0.038 (0.27) | 0.787 | 0.683 | +0.104 (0.054) |
| GPN-MSA | 0.724 | 0.796 | −0.072 | 0.859 | 0.793 | +0.066 (0.13) |
| CADD v1.7 | 0.501 | 0.534 | −0.033 | 0.613 | 0.532 | +0.081 (0.09) |
| ReMM v0.4 | 0.501 | 0.519 | −0.018 | — | — | — |

### The result is not explained by circularity or by position alone
- **Labels.** No ClinVar submitter applied PP3 or BP4 to any blind RNU2-2 positive. 21 of 22 positives list criteria, all from PS2, PS4, PM1, PM2, PM3 and PP1. AUROC was 0.896 for positives with listed criteria and no PP3/BP4, and 0.935 for de novo (PS2) positives.
- **Features.** v4's local-constraint feature counts each variant's own position. Recomputing it with that position left out, under the frozen coefficients, gave 0.893.
- **Position.** Four position-only baselines reached 0.728–0.752: disease-region membership, distance to it, 5′ proximity, and position unobserved in gnomAD. v4 exceeded each by +0.14 to +0.16, with every interval excluding 0. Within the guidance-defined n.1–67 region (13 vs 14 variants), v4 reached 0.780 [0.582, 0.940]; 5′ proximity reached 0.860, a tie within uncertainty.
- **Stricter negatives.** At AC ≥ 50, v4 reached 0.968, phyloP447 0.912, GPN-MSA 0.893 and CADD 0.669. Genome-wide conservation separates patient variants from frequent variants, but not from rarer population variants.

### A second blind test replicates the result
Amendment A8 was frozen before curation. The frozen rule accepted 35 RNU2-2 SNVs from the 2025–2026 literature, of which 16 were in neither the sealed set nor any training file (14 recessive, 2 dominant).
- **Pre-registered result:** frozen v4 reached 0.910 [0.825, 0.979] against 125 L3, with excess +0.133 (p = 0.012) and Δ vs CADD +0.298 [0.153, 0.426]. Both conditions were met (Table 1).
- **Pooled with the sealed set (38 positives):** v4 reached 0.901 (excess +0.122, p = 0.001), and the Rfam features had excess +0.17 (p = 0.001). CADD (+0.016) and GPN-MSA (−0.011) had no excess.
- **Broader label rules (post hoc):** the curator's STRICT set (20 new positives) and the union of both rules (27). v4 still exceeded CADD (+0.267 and +0.237, intervals excluding 0). Its excess over the rarity null fell to a trend (+0.079, p = 0.069; +0.068, p = 0.068).

The sealed set included de novo dominant variants (9 with PS2) and biallelic variants, whereas the second set was mostly recessive. The result therefore holds across both inheritance classes.

### Measured function agrees in RNU4-2, but not in RNU4ATAC
**RNU4-2 SGE** (435 SNVs; no population label involved; Fig. 3):
- reference-column frequency correlated with functional damage at ρ = 0.459 [0.376, 0.533];
- alternative-base rarity at 0.389 [0.305, 0.464];
- GPN-MSA 0.143, CADD 0.088, ReMM 0.062 and phyloP447 −0.024.

A v4 refitted without the RNU4 paralog cluster reached ρ = 0.382, above the 0.346 bar (post hoc; computed without the Evo 2 feature).

**RNU4ATAC assay (exploratory).** In the published splicing assay of 23 RNU4ATAC SNVs,³⁰ the pattern reversed:
- homology-transferred SGE function had ρ = 0.745 (n = 11, p = 0.008);
- GPN-MSA 0.581 (p = 0.004) and CADD 0.472 (p = 0.023);
- the Rfam features 0.344 (p = 0.11) and 0.168 (p = 0.44).

### The advantage is specific to the spliceosomal snRNA genes (pre-registered negative)
In the cross-family test (A7), Rfam alternative-base rarity did not beat genome-wide scores. The pre-registered test was not met: weighted ΔAUROC vs CADD −0.251 [−0.368, −0.133].
- **RMRP (28 P vs 40 N):** Rfam alternative-base rarity 0.488 [0.350, 0.636]; reference-column frequency 0.829; CADD 0.863; phyloP447 0.844.
- **SNORD118 and RNU7-1 (11 vs 134):** Rfam alternative-base rarity 0.879; CADD 0.816.
- **TERC:** too few negatives to evaluate.
- **Frozen v4:** tied CADD (Δ −0.005 [−0.060, 0.050]).

Genome-wide scores perform well on these single-copy RNA genes, and in RNU4ATAC (CADD 0.77) and RNU6ATAC (0.71). They fail in RNU2-2, RNU4-2 (0.52) and RNU12 (0.45). Across eight genes, the fraction of possible SNVs observed in gnomAD did not track CADD's performance (Spearman +0.21), so hypermutability does not explain the failure. Its cause remains open.

### Rarity inflation is general and large (pre-registered, A9–A13)
**ClinVar non-coding (A9).** A random unlabelled rare variant from the same window, scored against population-observed negatives, reached:
- CADD 0.632;
- GPN-MSA 0.676;
- phyloP447 0.647;
- ReMM 0.618.

The raw AUROCs of the ClinVar positives were 0.828, 0.805, 0.776 and 0.730 (Table 2). The inflation (null − 0.5) was +0.10 to +0.19 in every variant class. Measured against random rare variants instead, the rarity-matched AUROCs were:
- CADD 0.704 [0.695, 0.711];
- GPN-MSA 0.659 [0.646, 0.670];
- phyloP447 0.655 [0.643, 0.666];
- ReMM 0.633 [0.622, 0.639].

CADD exceeded every other score (paired differences +0.045 to +0.071; intervals excluding 0).

**Missense (A10).** AlphaMissense reached 0.920 against population-observed negatives (0.916 against ClinVar benign), but 0.694 [0.685, 0.704] against random rare missense variants. That is indistinguishable from CADD (0.691; difference +0.003 [−0.011, +0.016]). Against curated rare benign variants, AlphaMissense reached 0.967 (A14), so curated benign sets, even rare ones, are much easier than the rare-variant population.

**Published benchmarks.**
- **NCBoost (A11)**³²: every tool was inflated by +0.20 to +0.31. GPN-MSA ranked first on raw AUROC (0.837) and last after rarity matching (0.573 [0.551, 0.593]); CADD − GPN-MSA +0.087 [0.067, 0.102].
- **TraitGym Mendelian (A13)**³³: 97% of controls are population-observed. GPN-MSA ranked first on raw AUROC (0.905) and third after rarity matching (0.663; order CADD 0.717 > phyloP447 0.697 > GPN-MSA > ReMM 0.660). Intervals for this reordering were not computed, so it is descriptive.

The ground-truth simulation (A15) showed three things:
- raw AUROC overstates discrimination whenever a score tracks rarity, by up to +0.42 with no true signal;
- the rarity-matched AUROC is a lower bound;
- the difference raw − null is a valid test of "any signal beyond rarity" but over-corrects as a magnitude.

**Table 2.** Rarity-only null on ClinVar non-coding variants (A9; 1,787 P/LP vs 8,666 gnomAD AC ≥ 9 SNVs in ±100-bp windows). Rarity-matched AUROC intervals are from 200 window bootstraps (A15).

| Score | Raw AUROC | Rarity-only null | Inflation (null − 0.5) | Rarity-matched AUROC [95% CI] |
|---|---|---|---|---|
| CADD v1.7 | 0.828 | 0.632 | +0.132 | 0.704 [0.695, 0.711] |
| GPN-MSA | 0.805 | 0.676 | +0.176 | 0.659 [0.646, 0.670] |
| phyloP447 | 0.776 | 0.647 | +0.147 | 0.655 [0.643, 0.666] |
| ReMM v0.4 | 0.730 | 0.618 | +0.118 | 0.633 [0.622, 0.639] |

### Calibrated evidence for clinical use
**Pre-registered (Aim 4, v2 framework score; met).** At prior 0.10, PP3 Supporting was reached in:
- RNU2-2 (blind; function component ≥ 0.402);
- RNU4-2 (≥ 4.79);
- RNU4ATAC (PP3 ≥ 0.241; BP4 Supporting ≤ −0.046);
- the dominant pool (≥ 4.16).

No Moderate threshold was reached. CADD reached no PP3 threshold in RNU2-2 or RNU4-2.

**Post hoc (v4, out-of-sample; Table 3; Fig. 4).** Pooled over genes (99 P vs 1,117 L3), v4 reached an AUROC of 0.788 vs 0.746 for CADD. Thresholds were PP3 Supporting at v4 ≥ 0.362 and BP4 Supporting at ≤ −1.35. Within genes, v4 reached PP3 Supporting in RNU2-2 (≥ 0.216), RNU4-2 (≥ 0.189) and RNU4ATAC (≥ 0.88); CADD reached none in RNU2-2 or RNU4-2.

Of 279 ClinVar VUS or conflicting SNVs, 37 reach PP3 Supporting under their own gene's v4 threshold. Under pooled thresholds, 61 reach PP3 and 35 BP4. Because the negatives are population-observed variants, these likelihood ratios share the rarity contrast described above. A check of each threshold against random rare variants (positive-unlabelled lower bound) was pre-specified (`v4_rare_passrate.py`) and is reported in the Supplement.

**Table 3.** v4 PP3/BP4 Supporting thresholds at prior 0.10 (post hoc, out-of-sample scores; negatives L3).

| Unit | P / L3 | AUROC v4 (CADD) | PP3 Supporting | BP4 Supporting |
|---|---|---|---|---|
| Pooled (all genes) | 99 / 1,117 | 0.788 (0.746) | v4 ≥ 0.362 | v4 ≤ −1.35 |
| RNU2-2 | — | — | v4 ≥ 0.216 | — |
| RNU4-2 | — | — | v4 ≥ 0.189 | — |
| RNU4ATAC | — | — | v4 ≥ 0.88 | — |

### Inheritance mode from local constraint (secondary)
The frozen A6 rule (local-constraint score and threshold, set on non-blind genes) was applied to the literature-annotated blind-gene variants:
- **RNU2-2:** AUROC 0.839 [0.685, 0.968] for dominant (n = 4) vs recessive (n = 31) action, with threshold accuracy 86%;
- **RNU2-2 and RNU5B-1 together (6 vs 31):** 0.892 [0.763, 0.989].

Dominant-acting variants sat in heterozygote-depleted windows and recessive-acting variants did not. The dominant class is small.

### Gene-level constraint (exploratory; replication)
A within-locus window-depletion score ranked the five NDD snRNA genes above all 50 other loci (AUROC 1.0, permutation p = 5 × 10⁻⁵). Whole-locus observed/expected ratios ranked them near the bottom, because these genes are mutation hotspots. Both observations replicate prior work,¹,³⁵ and we do not claim them as discoveries.

## Discussion
Two blind tests on independent label sets agree that a frozen model separates RNU2-2 patient variants from population variants where CADD and ReMM are at chance. Measured beyond rarity, the signal comes from RNA-family alignment profiles. The same profiles track measured RNU4-2 function, where genome-wide scores do not. These results support RNA-family profiles as a source of computational evidence for the spliceosomal snRNA genes in which genome-wide predictors fail. The calibration supplies Supporting-strength thresholds that a laboratory can apply today.

The advantage has clear limits. It did not extend to RMRP, SNORD118 or RNU7-1, where genome-wide scores perform well. In RNU4ATAC, CADD, GPN-MSA and homology-transferred SGE tracked a functional assay better than the Rfam features. Genome-wide scores therefore fail in a subset of snRNA genes, RNU2-2, RNU4-2 and RNU12, and not in RNA genes generally. Hypermutability did not explain this pattern. One untested possibility is that paralogous copies confound the whole-genome alignments that underlie conservation scores at some loci. Rfam profiles are built from curated family alignments and do not depend on whole-genome alignment.

The pre-registered primary model failed. Population allele counts, even with UK Biobank added, did not resolve individual variants at these loci. The model that succeeded was designed after that failure was seen. It was frozen and anchored in a Bitcoin block before the blind labels were opened, and it then replicated on labels curated afterwards. We report both, because the failure carries information. Selection signals at hypermutable snRNA loci are too sparse for per-variant prediction, whereas family-level sequence information is not.

The rarity-only null matters beyond snRNAs. Frequency bias in benchmarks with common benign variants is recognised,¹⁸⁻²¹ and ClinGen calibration excludes common benign variants.²² We show that the bias persists when negatives are rare but population-observed, as is typical for small genes. Its size is score-specific, which makes raw comparisons between scores misleading. On RNU2-2, the null was 0.80 for GPN-MSA and 0.53 for CADD. Genome-wide, a random rare variant reached 0.62–0.68 against population-observed negatives. In published benchmarks, rarity matching reversed tool rankings: GPN-MSA fell from first to last on NCBoost, with an interval excluding a tie. The closest prior work stratified benchmarks by allele frequency.²¹ Our contribution is a per-score null that any benchmark with population-derived negatives can report, together with a simulation of how to interpret it.

Rarity matching also removed AlphaMissense's lead over CADD among rare ClinVar missense variants. This does not contradict trait-based rare-variant benchmarks in which AlphaMissense ranks first,³⁶ which measure a different endpoint. It does indicate that the advantage measured against population negatives overstates the advantage among the rare variants that laboratories classify.

**Limitations.**
- **One gene.** The blind evidence comes from one gene in two label sets (38 positives). RNU5B-1 could not be evaluated, and prospective units (A12) have not reached their pre-registered minimum.
- **Label rule.** In the second blind set, v4's excess over the rarity null was significant under the pre-registered label rule, but only a trend under broader rules. Its advantage over CADD held under all three.
- **Feature selection.** The Rfam features were singled out by a post hoc ablation, and their first confirmatory test, A7, was not met outside the spliceosomal snRNAs.
- **Calibration.** Negatives were population-observed variants rather than curated benign variants, which are too few in these genes. Thresholds reach Supporting strength only. The v4 calibration was post hoc.
- **Curation.** Literature curation was AI-assisted under a fixed rule, and the tables are released for audit.
- **Design.** No new functional data were generated, and the study has a single author.

## Conclusions
RNA-family alignment profiles, combined in a frozen model, provide blind-validated, ClinGen-calibrated computational evidence for variants in the spliceosomal snRNA genes where genome-wide scores fail, and the signal survives adjustment for rarity. Benchmarks that use population-derived negatives should report each score's rarity-only null, because its size is large and score-specific.

## Data availability
- The scores for all 17,385 SNVs (v2 and v4, with PP3/BP4 labels) will be deposited at Zenodo **[DOI to be added]**.
- The curated literature labels (`labels/A8_literature.tsv`, with rejected rows and reasons) and all benchmark selections (`a9/`–`a13/`) are in the repository.
- The sealed-label event log (`sealed/SEALED_MANIFEST.jsonl`) will be released with the scores.

## Code availability
All code, the pre-registration, amendments, OpenTimestamps proofs and the deviation log are available at https://github.com/tanukusaitejesh-prog/ndd-pp3-calibration (folder `snrna/`), and will be archived with the score release.

## Ethics declaration
The study used only public, de-identified data: ClinVar, gnomAD, UK Biobank allele-frequency summaries, published variant reports and public functional data. No ethics approval was required.

## Author contributions, competing interests, funding
[To be completed by the author.] AI assistance: analysis code, literature curation support and drafting were assisted by an AI system under the author's direction; all decisions and the final text are the author's responsibility.

## References
[Numbered reference list in `paper/REFERENCES.md`; every entry verified against Crossref/Europe PMC on 5 November 2026.]

## Figure legends
**Figure 1.** Blind RNU2-2 receiver-operating characteristic curves for frozen v4, Rfam alternative-base rarity, phyloP447 and CADD. (A) Sealed labels (22 P/LP vs 127 L3). (B) Literature-only labels (16 vs 125).

**Figure 2.** Excess AUROC over the rarity-only null for each score, on the sealed and literature-only RNU2-2 sets, with position-only baselines.

**Figure 3.** Spearman correlation with RNU4-2 SGE functional damage (435 SNVs), by score; v2, v3 and v4 relative to the pre-registered bar.

**Figure 4.** Pooled local likelihood ratio of out-of-sample v4 with bootstrap bands and PP3/BP4 evidence thresholds at prior 0.10.

**Figure 5.** Rarity-only inflation by score and benchmark: ClinVar non-coding (A9), ClinVar missense (A10), NCBoost (A11) and TraitGym Mendelian (A13).
