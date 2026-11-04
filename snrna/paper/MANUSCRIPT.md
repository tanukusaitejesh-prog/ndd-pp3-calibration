# RNA-family profiles provide blind-validated computational evidence for spliceosomal snRNA variants: a pre-registered study with a rarity-adjusted benchmark

**Running title:** Computational evidence for snRNA variants

**Author:** [Author name]¹ — ¹[Affiliation]. Correspondence: [email]

*Manuscript, version of 4 November 2026. Every number is traceable to the dated deviation log (`snrna/PREREG.md` §14), the result files in `snrna/results/`, and the verification report `paper/VERIFICATION_2026-11-04.md`. Citation numbers are generated from `paper/src/` by `paper/build_refs.py`.*

---

## Abstract

**Purpose.** Variants in spliceosomal small nuclear RNA (snRNA) genes are a frequent cause of neurodevelopmental disorders. Current guidance advises against computational evidence (ACMG/AMP PP3/BP4) for these genes until a tool has been calibrated on them, because genome-wide predictors perform poorly there.

**Methods.** We pre-registered and cryptographically timestamped an evaluation of variant-effect scores at 55 snRNA loci (17,385 possible single-nucleotide variants). Patient labels for RNU2-2 and RNU5B-1 were sealed before any model was fitted. A model was frozen and anchored in the Bitcoin blockchain before the labels were opened. Because negatives in small genes are population-observed variants, we scored each predictor against a score-specific rarity-only null.

**Results.** The pre-registered primary model, trained on population allele counts, failed its endpoints. A frozen model designed afterwards reached an AUROC of 0.892 on 22 sealed RNU2-2 pathogenic variants versus 127 population-observed variants (CADD 0.501). It reached 0.910 on 16 patient variants curated afterwards from the literature (Δ vs CADD +0.298 [0.153, 0.426]). Beyond rarity, the signal came from parameter-free RNA-family (Rfam) alignment profiles (excess +0.17, p = 0.001), not from CADD, ReMM, GPN-MSA or phyloP. In RNU4-2 saturation genome editing, with no population data involved, the Rfam reference-column frequency classified function-depleted variants with an AUROC of 0.80 [0.75, 0.85], against 0.49–0.60 for genome-wide scores (every paired interval excluded 0). The advantage did not extend to other RNA gene families (pre-registered negative). Provisional Supporting-strength thresholds were derived for RNU2-2, RNU4-2 and RNU4ATAC. Across 1,787 ClinVar non-coding pathogenic variants, random rare variants alone reached AUROCs of 0.62–0.68 against population-observed negatives.

**Conclusion.** RNA-family profiles provide blind-validated computational evidence for spliceosomal snRNA variants where genome-wide predictors fail. Benchmarks with population-derived negatives should report each score's rarity-only null.

**Keywords:** snRNA; RNU2-2; RNU4-2; variant interpretation; PP3/BP4; Rfam; benchmark bias; pre-registration

---

## Introduction

De novo variants in RNU4-2, which encodes the U4 small nuclear RNA (snRNA) of the major spliceosome, cause ReNU syndrome, one of the most frequent monogenic neurodevelopmental disorders (NDDs) found by genome sequencing.<sup>1,2</sup> Other snRNA disorders followed:
- dominant and recessive disorders of RNU2-2;<sup>3–6</sup>
- dominant disorders of RNU5B-1;<sup>4,7</sup>
- a recessive RNU4-2 syndrome.<sup>8,9</sup>

The minor-spliceosome snRNA disorders were known earlier: RNU4ATAC,<sup>10,11</sup> RNU12,<sup>12,13</sup> and RNU6ATAC.<sup>14,15</sup> Dominant variants in RNU4-2 and RNU6 genes also cause retinitis pigmentosa.<sup>16</sup> Diagnostic laboratories now encounter snRNA variants routinely.<sup>17</sup>

Interpretation lacks computational evidence. The proposed guidance for snRNA variant classification states that no tools yet capture features relevant to snRNA variant deleteriousness and that genome-wide tools perform poorly.<sup>18</sup> It therefore recommends using in silico evidence (PP3/BP4) only after a tool has been calibrated for snRNAs. Under thresholds calibrated on missense variants,<sup>19</sup> CADD scores the first two pathogenic RNU2-2 variants identified in the benign range.<sup>18</sup> More generally, genome-wide thresholds calibrated on missense variants give false benign evidence for pathogenic non-coding variants.<sup>20</sup> Saturation genome editing (SGE) of RNU4-2 now provides functional measurements for one gene,<sup>9</sup> but not for the others.

snRNA genes are unusual targets for genome-wide predictors. They occur in multiple paralogous copies among many pseudogenes,<sup>18,21</sup> and the spliceosomal snRNA genes are germline mutation hotspots.<sup>4,22</sup> RNA-family covariance models (Rfam) take a different approach. They describe each snRNA family's sequence and secondary structure across species and paralogs, independently of whole-genome alignment.<sup>23,24</sup>

Evaluation is a second obstacle. Benchmarks often use population variants as negatives, and rarity correlates with conservation and constraint. Such negatives can therefore flatter a predictor. This has been described for missense predictors evaluated with common benign variants<sup>25–28</sup> and for predictors evaluated across heterogeneous genomic contexts.<sup>29</sup> ClinGen calibration accordingly restricts benign variants to allele frequency below 1%.<sup>19</sup> Small genes rarely have curated benign variants, so their negatives are usually defined by population observation; here, a gnomAD allele count of at least 9. Most such negatives are rare (median allele frequency 1.8 × 10⁻⁴), yet any score that tracks rarity still separates a random unlabelled variant from them. We therefore report each score's excess over a score-specific rarity-only null, and we test the main claim on functional data that involve no population labels.

We designed the study to be falsifiable. The plan, definitions and decision rules were frozen and timestamped before analysis. Patient labels for two genes were sealed in an encrypted file before any model was fitted. Every later change is logged together with what had been seen when it was made. We report the failed endpoints alongside the successful ones.

## Materials and Methods

### Study design and pre-registration
The pre-registration (`PREREG.md`), analysis code and each amendment were hashed and timestamped with OpenTimestamps; the proofs are anchored in Bitcoin blocks (Supplementary Table S1). Changes made after data were seen are recorded in a dated deviation log (`PREREG.md` §14), with what had been viewed at the time. Throughout, analyses are labelled as pre-registered (with amendment number), post hoc, or exploratory. Post hoc analyses whose code was committed before they were run are labelled "specified before running"; git commit times are the authoritative record.

### Loci, variants and scores
We included the GENCODE v49 spliceosomal snRNA genes and paralogs that passed coverage and paralog quality control. This gave 55 loci, 20 of them named genes, with all 17,385 possible SNVs. Comparators were:
- CADD v1.7;<sup>30</sup>
- phyloP from the 447-way alignment;<sup>31,32</sup>
- ReMM v0.4;<sup>33</sup>
- GPN-MSA;<sup>34</sup>
- Evo 2 (7B) log-likelihood ratios;<sup>35</sup>
- ViennaRNA 2 folding and duplex free-energy changes for intermolecular snRNA pairs.<sup>36</sup>

RNA-family features came from aligning each locus to its Rfam 15.1 covariance model with Infernal 1.1.5 (cmalign).<sup>23,24</sup> We used two features:
- **reference-column frequency:** the frequency of the reference base in its alignment column;
- **alternative-base rarity:** one minus the frequency of the alternative base.

Neither feature uses a label, and neither has a fitted parameter.

### Labels
**Positives** were pathogenic or likely pathogenic (P/LP) SNVs from three sources:
- ClinVar, with at least one review star;<sup>37</sup>
- RNUdb curation;<sup>18</sup>
- a literature tally under a frozen wording rule: a source must call the variant pathogenic, likely pathogenic or disease-causing, and any source calling it a VUS excludes it.

**Negatives (L3)** were gnomAD v4.1 genome PASS SNVs<sup>38</sup> with allele count ≥ 9, excluding positives. L3 holds 1,233 SNVs. They are population-observed but mostly rare: median allele frequency 1.8 × 10⁻⁴, with 5.4% at ≥ 1%.

**Sealing.** RNU2-2 and RNU5B-1 labels were fetched directly into an AES-256-encrypted file before any model was fitted. They were decrypted once, and every access was logged.

### Models
**v2 (primary, pre-registered).** A graph-transformer selection model trained on gnomAD allele counts under the Roulette mutation-rate model,<sup>39</sup> with paralog clusters as the hold-out unit. Its framework score was inheritance-aware:
- dominant genes used selection;
- recessive genes used lab-derived function (RNU4-2 SGE, transferred to RNU4ATAC by structural alignment).

A version with UK Biobank allele counts added (v3) was exploratory.<sup>40</sup>

**v4 (amendment A5).** v4 was designed after the v2 endpoints were seen. It is an L2-regularised logistic regression (C = 0.1, balanced class weights, no tuning) on within-gene percentile ranks of 13 features:
- the six comparator scores;
- the two Rfam features;
- reference pairing and intermolecular-pairing loss;
- absolute folding ΔΔG;
- local constraint: the negative z-score of observed versus expected variation in a 21-nt window, with gnomAD and UK Biobank counts;
- homology-transferred RNU4-2 SGE function (RNU4ATAC only).

It was trained on non-blind genes only; the blind genes never entered training. Its code and output hashes were committed and anchored in Bitcoin block 971362 (14 October 2026, 05:58 UTC), before the blind labels were decrypted. Its blind evaluation script checks those hashes before scoring.

### Evaluation and the rarity-only null
The primary metric was the area under the ROC curve (AUROC) for positives against L3, with 95% class-stratified bootstrap intervals (2,000 replicates).

Because L3 is defined by population observation, any score that tracks rarity separates random unlabelled variants from L3. For *n* positives in a gene, the **rarity-only null** works as follows:
- draw *n* SNVs at random, without replacement, from that gene's SNVs that are neither positives nor L3;
- score the draw against the same L3;
- repeat 1,000 times.

We report the null mean and the **excess** (observed AUROC − null mean), with a one-sided permutation p-value, (k + 1)/(B + 1).

The null pool contains unlabelled damaging variants, so the excess is a conservative test statistic rather than an effect size. In a ground-truth simulation (A15; Supplement S6):
- the raw AUROC overstated the true discrimination between pathogenic and rare benign variants, so it is an upper bound;
- the rarity-matched AUROC (positives against random rare variants) was a lower bound;
- a correction for the damaging fraction π, (M − π/2)/(1 − π), recovered the true value.

### Blind tests
**Test 1: sealed labels.** RNU2-2 had 22 P/LP vs 127 L3 SNVs. RNU5B-1 had 1 P/LP and could not be evaluated.

**Test 2: literature (amendment A8).** The design was fixed before curation. RNU2-2 and RNU5B-1 variants reported as pathogenic in the 2025–2026 literature were curated under the frozen wording rule. Curation was AI-assisted, and the curator was instructed not to look up any score. Positives absent from the sealed set and from every training file formed the test set, and the labels were committed before scoring. The pre-registered conditions were:
- v4 excess over its rarity null with p < 0.05;
- a v4 − CADD AUROC interval excluding 0.

**Test 3: prospective (amendment A12).** Units of new P/LP variants reported after the freeze are evaluated once each, when they reach a pre-registered minimum (Supplement S11).

### Robustness analyses (post hoc, specified before running)
- **Circularity:** the ACMG/AMP criteria that each ClinVar submitter applied to each blind positive, and v4 recomputed with leave-own-position local constraint.
- **Position-only baselines:** membership of guidance-defined disease regions,<sup>18</sup> distance to them, 5′ proximity, and whether the position is observed in gnomAD.
- **Stricter negatives:** L3 restricted to AC ≥ 20, 50 or 100.

### Functional readouts
- **RNU4-2 SGE correlation:** Spearman correlation with functional damage (−function score) for the 435 in-transcript RNU4-2 SNVs with a score.<sup>9</sup>
- **RNU4-2 SGE classification (post hoc, specified before running):** AUROC for "depleted" variants (function score < −0.302, the study's threshold) and "strong" variants (< −0.9). We used 2,000 class-stratified bootstraps and paired differences between each Rfam feature and each genome-wide score. Analyses were repeated inside and outside the guidance-defined critical region, because position alone carries much of the signal in RNU4-2.
- **RNU4ATAC splicing assay (exploratory):** published splicing efficiencies of 23 RNU4ATAC SNVs.<sup>41</sup>

No population label or negative enters the SGE analyses.

### Cross-family test (amendment A7)
To test whether the Rfam advantage generalises, alternative-base rarity was scored zero-shot on RMRP, TERC, SNORD118 and RNU7-1, with scores frozen before labels were fetched. The primary feature had been chosen from the post hoc RNU2-2 ablation before any A7 data existed.

### Alignment depth (exploratory, specified before running)
For each named gene we summarised phyloP447 and phyloP100 across positions:
- the 95th-percentile value (a ceiling);
- the fraction of positions sharing the gene's most frequent value;
- the per-position Spearman correlation between phyloP447 and Rfam reference-column frequency.

### Benchmark analyses (amendments A9–A15)
- **A9:** 1,787 ClinVar non-coding P/LP SNVs against 8,666 gnomAD PASS SNVs with AC ≥ 9 within ±100 bp. The positives were intronic (1,500, sampled from 3,438), 5′UTR (148), ncRNA exon (85), 3′UTR (51) and upstream (3). The window-matched null drew, for each positive, a random non-positive SNV with AC < 9 from the same window.
- **A10:** the same design for missense variants, adding AlphaMissense.<sup>42</sup>
- **A11:** the NCBoost published benchmark (common dbSNP negatives).<sup>43</sup>
- **A13:** the TraitGym Mendelian benchmark (published matched controls).<sup>44</sup>
- **A15:** a ground-truth simulation of the null, with window-bootstrap intervals.

### Calibration
We used the ClinGen local likelihood-ratio method, as pre-registered:<sup>19</sup>
- nearest 20% of labelled variants, at least 10 per class;
- 1,000 bootstraps;
- thresholds from the 5th (PP3) or 95th (BP4) bootstrap percentile of the local LR;
- primary prior 0.10,<sup>45</sup> with sensitivity priors 0.01, 0.0441, 0.2 and 0.5.

The pre-registered calibration (Aim 4) used the v2 framework score. Calibration of v4 was post hoc and used out-of-sample scores: frozen v4 for the blind genes, and v4 refitted without the gene's paralog cluster elsewhere. Negatives were L3. ClinVar VUS and conflicting SNVs in the named genes (279 SNVs) were scored against the thresholds.

### Statistics and software
Analyses used Python 3 (NumPy, pandas, SciPy, scikit-learn). P-values are one-sided where the pre-registration specified a direction. The rarity-null tests on the sealed set were added after unsealing and span 12 scores; we report them unadjusted and give Bonferroni-adjusted values for the main claims. All code, inputs and logs are in the repository.

## Results

### The pre-registered primary model failed per variant
v2 did not meet its pre-registered endpoint against RNU4-2 SGE function. Its Spearman ρ was 0.144, against a bar of 0.346 (Evo 2's ρ of 0.246 plus 0.10). Adding UK Biobank counts (v3) did not help (ρ 0.142).

v2 also failed the held-out-gene gate, with ΔAUROC vs CADD of +0.083 against a bar of 0.15. Population allele counts at these hypermutable loci therefore did not resolve individual variants. On the sealed RNU2-2 labels, v2's framework score (0.784) exceeded CADD (0.501). This came from its lab-transfer function component (0.851), not its selection component (0.659).

### The frozen v4 model separates blind RNU2-2 patient variants from population variants
In leave-one-paralog-cluster-out validation on non-blind genes, v4 exceeded CADD by +0.186 [0.108, 0.267]. On the sealed RNU2-2 labels (22 P/LP vs 127 L3; Fig. 1A), v4 reached an AUROC of 0.892 [0.817, 0.951]. The comparators were:
- CADD 0.501 [0.393, 0.612];
- ReMM 0.501;
- phyloP447 0.724;
- GPN-MSA 0.724.

### Most of the apparent gap is rarity; the RNA-family profile carries the rest
Because L3 is defined by population observation, scores differ in how far rarity alone carries them. On RNU2-2, the rarity-only null was 0.783 for v4, 0.796 for GPN-MSA, 0.686 for phyloP447 and 0.534 for CADD (Table 1; Fig. 2). Measured as excess over the null:
- **Rfam features** kept a large signal: alternative-base rarity +0.175 (p = 0.001) and reference-column frequency +0.172 (p = 0.002). These remain significant after Bonferroni correction over the 12 scores (0.012 and 0.024).
- **v4** kept +0.108 (p = 0.014; about 0.17 after correction).
- **Genome-wide scores** kept none: phyloP447 +0.038 (p = 0.27), CADD −0.033, ReMM −0.018 and GPN-MSA −0.072.

v4 did not improve on alternative-base rarity alone (+0.031; interval including 0). The headline gap between v4 and CADD therefore overstates the advantage, and the useful signal is the RNA-family profile. The two Rfam features were singled out by a post hoc ablation on the sealed labels. Their claim to replicate rests on the second blind test and on the functional data below.

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

Sealed-set p-values are unadjusted (post hoc, 12 scores). The A8 test for v4 was pre-registered.

### The result is not explained by circularity or by position alone
- **Labels.** No ClinVar submitter applied PP3 or BP4 to any blind RNU2-2 positive. 21 of 22 positives list criteria, all from PS2, PS4, PM1, PM2, PM3 and PP1. v4 reached 0.896 for the positives with listed criteria and no PP3/BP4, and 0.935 for de novo (PS2) positives.
- **Features.** v4's local-constraint feature counts each variant's own position. Recomputing it with that position left out, under the frozen coefficients, gave 0.893.
- **Position.** Four position-only baselines reached 0.728–0.752: disease-region membership, distance to it, 5′ proximity, and position unobserved in gnomAD. v4 exceeded each by +0.14 to +0.16, with every interval excluding 0. Within the guidance-defined n.1–67 region (13 vs 14 variants), v4 reached 0.780 [0.582, 0.940] and 5′ proximity 0.860, a tie within uncertainty.
- **Stricter negatives.** At AC ≥ 50, v4 reached 0.968, phyloP447 0.912, GPN-MSA 0.893 and CADD 0.669. Genome-wide conservation separates patient variants from frequent variants, but not from rarer population variants.

### A second blind test replicates the result
Amendment A8 was frozen before curation. The frozen rule accepted 35 RNU2-2 SNVs from the 2025–2026 literature. Of these, 16 were in neither the sealed set nor any training file (14 recessive, 2 dominant).
- **Pre-registered result:** frozen v4 reached 0.910 [0.825, 0.979] against 125 L3, with excess +0.133 (p = 0.012) and Δ vs CADD +0.298 [0.153, 0.426]. Both conditions were met (Table 1).
- **Pooled with the sealed set (38 positives):** v4 reached 0.901 (excess +0.122, p = 0.001), and the Rfam features had excess +0.17 (p = 0.001). CADD (+0.016) and GPN-MSA (−0.011) had no excess.
- **Broader label rules (post hoc):** the curator's STRICT set (20 new positives) and the union of both rules (27). v4 still exceeded CADD (+0.267 and +0.237, intervals excluding 0), but its excess over the rarity null fell to a trend (+0.079, p = 0.069; +0.068, p = 0.068).

The two tests share their negatives: 125 of the sealed test's 127 L3 variants. The replication is therefore on new patient variants in the same gene, against nearly the same population variants. The sealed set included de novo dominant variants (9 with PS2), whereas the second set was mostly recessive.

### Measured function: the family profile classifies RNU4-2 variants without population labels
**SGE correlation** (435 RNU4-2 SNVs; Fig. 3A):
- reference-column frequency correlated with functional damage at ρ = 0.459 [0.377, 0.539];
- alternative-base rarity at 0.389 [0.311, 0.464];
- ViennaRNA intermolecular ΔΔG at 0.237;
- GPN-MSA 0.143, CADD 0.088, ReMM 0.062 and phyloP447 −0.024.

A v4 refitted without the RNU4 paralog cluster reached ρ = 0.382, above the 0.346 bar (post hoc; computed without the Evo 2 feature).

**SGE classification** (Table 2; Fig. 3B). Among the 435 SNVs, 106 were function-depleted.
- Reference-column frequency classified them with an AUROC of 0.802 [0.753, 0.848]. Genome-wide scores reached 0.489–0.604.
- Each paired difference was +0.197 to +0.313, with every interval excluding 0.
- For strongly depleted variants (43), the AUROC was 0.865, against at most 0.636 for genome-wide scores.
- Outside the guidance-defined critical region, the AUROC was 0.774 against at most 0.572, with every paired interval excluding 0.
- Inside the critical region (31 vs 8 variants), scores did not differ.

These analyses involve no population data, so the rarity contrast that affects every patient-versus-population comparison cannot arise. They use the same SGE variants as the correlation above, expressed as a classification.

**Table 2.** RNU4-2 SGE classification of function-depleted variants (function score < −0.302; post hoc, specified before running). AUROC [95% CI]; paired differences against reference-column frequency.

| Score | All (106 vs 329) | Outside critical region (75 vs 321) | Strong, < −0.9 (43 vs 392) | Rfam reference − score, all |
|---|---|---|---|---|
| Rfam reference-column frequency | 0.802 [0.753, 0.848] | 0.774 [0.711, 0.834] | 0.865 [0.822, 0.904] | — |
| Rfam alternative-base rarity | 0.733 [0.681, 0.785] | 0.722 [0.662, 0.779] | 0.771 [0.717, 0.823] | — |
| ViennaRNA ΔΔG (intermolecular) | 0.629 [0.568, 0.690] | 0.624 [0.556, 0.690] | 0.637 [0.540, 0.733] | — |
| GPN-MSA | 0.604 [0.547, 0.662] | 0.572 [0.505, 0.641] | 0.617 [0.543, 0.690] | +0.197 [0.122, 0.268] |
| CADD v1.7 | 0.564 [0.503, 0.624] | 0.510 [0.435, 0.583] | 0.633 [0.555, 0.714] | +0.237 [0.166, 0.309] |
| ReMM v0.4 | 0.513 [0.451, 0.577] | 0.495 [0.413, 0.577] | 0.636 [0.566, 0.705] | +0.288 [0.224, 0.352] |
| phyloP447 | 0.489 [0.435, 0.543] | 0.473 [0.403, 0.542] | 0.502 [0.421, 0.584] | +0.313 [0.237, 0.386] |

**RNU4ATAC assay (exploratory).** The pattern reversed in the published splicing assay of 23 RNU4ATAC SNVs:<sup>41</sup>
- homology-transferred SGE function had ρ = 0.745 (n = 11, p = 0.008);
- GPN-MSA had 0.581 (p = 0.004) and CADD 0.472 (p = 0.023);
- the Rfam features had 0.344 (p = 0.11) and 0.168 (p = 0.44).

### The advantage is specific to the spliceosomal snRNA genes (pre-registered negative)
In the cross-family test (A7), Rfam alternative-base rarity did not beat genome-wide scores. The pre-registered test was not met, with a weighted ΔAUROC vs CADD of −0.251 [−0.368, −0.133].
- **RMRP (28 P vs 40 N):** Rfam alternative-base rarity 0.488 [0.350, 0.636]; reference-column frequency 0.829; CADD 0.863; phyloP447 0.844.
- **SNORD118 and RNU7-1 (11 vs 134):** Rfam alternative-base rarity 0.879; CADD 0.816.
- **TERC:** too few negatives to evaluate.
- **Frozen v4:** tied CADD (Δ −0.005 [−0.060, 0.050]).

Genome-wide scores perform well on these single-copy RNA genes. Over all labels per gene, they also work in RNU4ATAC (CADD 0.77) and RNU6ATAC (0.71), but fail in RNU2-2 (0.55), RNU4-2 (0.52) and RNU12 (0.45). Across eight genes, the fraction of possible SNVs observed in gnomAD did not track CADD's performance (Spearman +0.21), so hypermutability alone does not explain the failure.

The alignment-depth check found one candidate mechanism, for RNU4-2:
- phyloP447 reaches at most 0.90 at RNU4-2, and 59% of its positions share a single value;
- its paralog RNU4-1 reaches 8.03;
- phyloP447 is uncorrelated with the Rfam profile across RNU4-2 positions (ρ −0.01);
- in RNU4ATAC it agrees with the profile (ρ 0.55).

The 447-way alignment therefore carries almost no informative conservation at RNU4-2, consistent with the SGE result (phyloP447 AUROC 0.489). Across the five genes with CADD estimates, the agreement between phyloP447 and the Rfam profile tracked CADD's performance (Spearman +0.90; descriptive, five genes). RNU12 has a high phyloP447 ceiling (7.71) yet CADD fails there, so alignment depth is at most a partial explanation.

### Rarity inflation is general and large (pre-registered, A9–A13)
**ClinVar non-coding (A9).** A random unlabelled rare variant from the same window, scored against population-observed negatives, reached:
- CADD 0.632;
- GPN-MSA 0.676;
- phyloP447 0.647;
- ReMM 0.618.

The raw AUROCs of the ClinVar positives were 0.828, 0.805, 0.776 and 0.730 (Table 3). The inflation (null − 0.5) was +0.10 to +0.19 in every variant class. Against random rare variants instead, the rarity-matched AUROCs (a lower bound) were:
- CADD 0.704;
- GPN-MSA 0.659;
- phyloP447 0.655;
- ReMM 0.633.

CADD ranked first both ways. Intervals for these rankings await the corrected bootstrap (Supplement S6), so we describe them without significance claims.

**Missense (A10).** AlphaMissense reached 0.920 against population-observed negatives and 0.916 against ClinVar benign variants. Against random rare missense variants it reached 0.694, close to CADD (0.690). Against curated rare benign variants it reached 0.967 (A14). Curated benign sets, even rare ones, are therefore much easier to separate from pathogenic variants than the rare-variant population is.

**Published benchmarks (descriptive).**
- **NCBoost (A11):**<sup>43</sup> the null inflated every tool by +0.20 to +0.31. GPN-MSA ranked first on raw AUROC (0.837) and last after rarity matching (0.573). CADD ranked first after rarity matching (0.660).
- **TraitGym Mendelian (A13):**<sup>44</sup> 97% of the controls are population-observed. GPN-MSA ranked first on raw AUROC (0.905) and third after rarity matching (0.663). The rarity-matched order was CADD 0.717, phyloP447 0.697, GPN-MSA, then ReMM 0.660.

In both benchmarks the null draws from the positives' windows, whereas the published negatives lie elsewhere. The inflation therefore combines rarity with local sequence context.

**Ground-truth simulation (A15).** The simulation showed three things:
- raw AUROC overstates discrimination whenever a score tracks rarity, by up to +0.42 with no true signal;
- the rarity-matched AUROC is a lower bound;
- the difference raw − null is a valid test of "any signal beyond rarity" but over-corrects as a magnitude.

**Table 3.** Rarity-only null on ClinVar non-coding variants (A9; 1,787 P/LP vs 8,666 gnomAD AC ≥ 9 SNVs in ±100-bp windows).

| Score | Raw AUROC (upper bound) | Rarity-only null | Inflation (null − 0.5) | Rarity-matched AUROC (lower bound) |
|---|---|---|---|---|
| CADD v1.7 | 0.828 | 0.632 | +0.132 | 0.704 |
| GPN-MSA | 0.805 | 0.676 | +0.176 | 0.659 |
| phyloP447 | 0.776 | 0.647 | +0.147 | 0.655 |
| ReMM v0.4 | 0.730 | 0.618 | +0.118 | 0.633 |

### Provisional calibrated evidence
**Pre-registered (Aim 4, v2 framework score; met).** At prior 0.10, PP3 Supporting was reached in:
- RNU2-2 (blind; function component ≥ 0.402);
- RNU4-2 (≥ 4.79);
- RNU4ATAC (PP3 ≥ 0.241; BP4 Supporting ≤ −0.046);
- the dominant pool (≥ 4.16).

No Moderate threshold was reached. CADD reached no PP3 threshold in RNU2-2 or RNU4-2.

**Post hoc (v4, out-of-sample; Table 4; Fig. 4).** Pooled over genes (99 P vs 1,117 L3), v4 reached an AUROC of 0.788, against 0.746 for CADD. Its pooled thresholds were PP3 Supporting at v4 ≥ 0.362 and BP4 Supporting at v4 ≤ −1.35. Within genes, v4 reached PP3 Supporting in RNU2-2 (≥ 0.216), RNU4-2 (≥ 0.189) and RNU4ATAC (≥ 0.88). CADD reached none in RNU2-2 or RNU4-2.

Of 279 ClinVar VUS or conflicting SNVs, 37 reach PP3 Supporting under their own gene's v4 threshold. Under pooled thresholds, 61 reach PP3 and 35 reach BP4.

These thresholds are provisional. The negatives are population-observed variants, so the likelihood ratios share the rarity contrast described above and may overstate evidence against a rare benign variant. A pre-specified check of each threshold against random rare variants (`v4_rare_passrate.py`; a lower bound on the likelihood ratio) requires the sealed scores and is pending.

**Table 4.** v4 PP3/BP4 Supporting thresholds at prior 0.10 (post hoc, out-of-sample scores; negatives L3; provisional).

| Unit | P / L3 | AUROC v4 (CADD) | PP3 Supporting | BP4 Supporting |
|---|---|---|---|---|
| Pooled (all genes) | 99 / 1,117 | 0.788 (0.746) | v4 ≥ 0.362 | v4 ≤ −1.35 |
| RNU2-2 | — | — | v4 ≥ 0.216 | — |
| RNU4-2 | — | — | v4 ≥ 0.189 | — |
| RNU4ATAC | — | — | v4 ≥ 0.88 | — |

### Prospective test (A12)
ClinVar re-scans and a literature tally of four 2026 reports<sup>8,17,46,47</sup> have not yet brought any pre-registered unit to its minimum (Supplement S11). Units in genes never seen in training have 0–4 new positives, and the unit in the trained genes has 4 of the 5 required. A12 therefore remains descriptive, and each unit will be evaluated once when it qualifies.

## Discussion
Two blind tests on independent patient-variant sets agree that a frozen model separates RNU2-2 patient variants from population variants where CADD and ReMM are at chance. Measured beyond rarity, the signal comes from RNA-family alignment profiles. The same profiles classify measured RNU4-2 function without any population data, where genome-wide scores reach AUROCs of only 0.49–0.60. RNA-family profiles are therefore a source of computational evidence for the spliceosomal snRNA genes in which genome-wide predictors fail, the gap identified by the snRNA classification guidance.<sup>18</sup> The calibration supplies provisional Supporting-strength thresholds.

The advantage has clear limits. It did not extend to RMRP, SNORD118 or RNU7-1, where genome-wide scores perform well. In RNU4ATAC, CADD, GPN-MSA and homology-transferred SGE tracked a functional assay better than the Rfam features. Genome-wide scores fail in a subset of snRNA genes (RNU2-2, RNU4-2 and RNU12), not in RNA genes generally. Hypermutability did not explain this pattern. For RNU4-2, the whole-genome alignment underlying phyloP carries almost no informative conservation, while the paralog RNU4-1 does. This is consistent with paralogous copies disrupting orthology in whole-genome alignments. Rfam profiles are built from curated family alignments and avoid this problem. The explanation does not fit RNU12, and it was checked descriptively across few genes.

The pre-registered primary model failed. Population allele counts, even with UK Biobank added, did not resolve individual variants at these loci. The model that succeeded was designed after that failure was seen. It was frozen and anchored in a Bitcoin block before the blind labels were opened, and it then replicated on labels curated afterwards. We report both results, because the failure is informative: selection signals at hypermutable snRNA loci are too sparse for per-variant prediction, whereas family-level sequence information is not.

**The rarity-only null.** The null matters beyond snRNAs. Frequency bias in benchmarks with common benign variants is well recognised,<sup>25–28</sup> and ClinGen calibration excludes common benign variants.<sup>19</sup> Lu et al. recently showed that heterogeneity in genomic context inflates predictor performance.<sup>29</sup> We show that a related bias persists when negatives are rare but population-observed, as is typical for small genes. Its size is score-specific, so raw comparisons between scores can mislead. On RNU2-2, the null was 0.80 for GPN-MSA and 0.53 for CADD. Genome-wide, a random rare variant reached 0.62–0.68 against population-observed negatives. On two published benchmarks, rarity matching changed the order of tools; this is descriptive until corrected intervals are available. Our contribution is a per-score null that any benchmark with population-derived negatives can report, together with a simulation showing how to interpret it: raw AUROC as an upper bound and rarity-matched AUROC as a lower bound.

The same concern applies to calibration. Recent non-coding calibrations have used population variants as controls,<sup>48</sup> and our v4 thresholds do too. Such likelihood ratios carry the rarity contrast and should be checked against rare variants, as we have pre-specified. Predictor-informed classifications are also feeding back into ClinVar,<sup>49</sup> which favours classification-independent benchmarks such as the SGE classification reported here. Rarity matching also removed AlphaMissense's lead over CADD among rare ClinVar missense variants. This does not contradict trait-based benchmarks in which AlphaMissense ranks first,<sup>50</sup> which measure a different endpoint. It does indicate that an advantage measured against population negatives can overstate the advantage among the rare variants that laboratories classify.

**Limitations.**
- **One gene for the blind tests.** The blind evidence comes from one gene in two label sets (38 positives) that share their negatives. RNU5B-1 could not be evaluated, and no prospective unit (A12) has reached its pre-registered minimum.
- **Multiplicity and label rules.** The sealed-set rarity tests were post hoc across 12 scores; after correction the Rfam features remain significant and v4 does not. In the second blind set, v4's excess over the rarity null was significant under the pre-registered label rule but only a trend under broader rules. Its advantage over CADD held under all three rules.
- **Feature selection.** The Rfam features were singled out by a post hoc ablation. Their first confirmatory test (A7) was not met outside the spliceosomal snRNAs. Their replication rests on the second blind test and on the SGE data, which are post hoc analyses of public data.
- **Function in one gene.** The population-free evidence comes from RNU4-2 alone. In RNU4ATAC, a small assay favoured other scores.
- **Calibration.** Negatives were population-observed variants rather than curated benign variants, which are too few in these genes. Thresholds reach Supporting strength only. The v4 calibration was post hoc, and its rare-variant check is pending.
- **Benchmarks.** Corrected intervals for the A9–A13 rankings are pending. AlphaGenome and other recent sequence models were not evaluated.
- **Curation.** Literature curation was AI-assisted under a fixed rule; the tables are released for audit.
- **Design.** No new functional data were generated, and the study has a single author.

## Conclusions
RNA-family alignment profiles provide blind-validated computational evidence for variants in the spliceosomal snRNA genes where genome-wide scores fail. The signal survives adjustment for rarity and classifies measured RNU4-2 function without population labels. Combined in a frozen model, the profiles yield provisional ClinGen-style Supporting thresholds. Benchmarks that use population-derived negatives should report each score's rarity-only null, because its size is large and score-specific.

## Data availability
- The scores for all 17,385 SNVs (v2 and v4, with PP3/BP4 labels) will be deposited at Zenodo **[DOI to be added]**.
- The curated literature labels (`labels/A8_literature.tsv` and `labels/A12_*.tsv`, with rejected rows and reasons) and all benchmark selections (`a9/`–`a13/`) are in the repository.
- The sealed-label event log (`sealed/SEALED_MANIFEST.jsonl`) will be released with the scores.

## Code availability
All code, the pre-registration, amendments, OpenTimestamps proofs and the deviation log are available at https://github.com/tanukusaitejesh-prog/ndd-pp3-calibration (folder `snrna/`), and will be archived with the score release.

## Ethics declaration
The study used only public, de-identified data: ClinVar, gnomAD, UK Biobank allele-frequency summaries, published variant reports and public functional data. No ethics approval was required.

## Author contributions, competing interests, funding
[To be completed by the author.] AI assistance: analysis code, literature curation support and drafting were assisted by an AI system under the author's direction; all decisions and the final text are the author's responsibility.

## References
1. Chen Y, Dawes R, Kim HC, Ljungdahl A, Stenton SL, Walker S, et al. De novo variants in the RNU4-2 snRNA cause a frequent neurodevelopmental syndrome. Nature. 2024;632(8026):832-840. doi:10.1038/s41586-024-07773-7. PMID: 38991538.
2. Greene D, Thys C, Berry IR, Jarvis J, Ortibus E, Mumford AD, et al. Mutations in the U4 snRNA gene RNU4-2 cause one of the most prevalent monogenic neurodevelopmental disorders. Nat Med. 2024;30(8):2165-2169. doi:10.1038/s41591-024-03085-5. PMID: 38821540.
3. Greene D, De Wispelaere K, Lees J, Codina-Solà M, Jensson BO, Hales E, et al. Mutations in the small nuclear RNA gene RNU2-2 cause a severe neurodevelopmental disorder with prominent epilepsy. Nat Genet. 2025;57(6):1367-1373. doi:10.1038/s41588-025-02159-5. PMID: 40210679.
4. Jackson A, Thaker N, Blakes A, Rice G, Griffiths-Jones S, Balasubramanian M, et al. Analysis of R-loop forming regions identifies RNU2-2 and RNU5B-1 as neurodevelopmental disorder genes. Nat Genet. 2025;57(6):1362-1366. doi:10.1038/s41588-025-02209-y. PMID: 40442284. Publisher Correction: Nat Genet. 2025;57(7):1792.
5. Greene D, Mendez R, Lees J, Barbosa M, Bruselles A, Chiriatti L, et al. Biallelic variants in RNU2-2 cause the most prevalent known recessive neurodevelopmental disorder. Nat Genet. 2026;58(4):774-781. doi:10.1038/s41588-026-02539-5. PMID: 41912932. Publisher Correction: Nat Genet. 2026;58(5):1187.
6. Leitão E, Santini A, Cogne B, Essid M, Athanasiadou M, LaFlamme CW, et al. Systematic analysis of snRNA genes reveals frequent RNU2-2 variants in dominant and recessive developmental and epileptic encephalopathies. Nat Genet. 2026;58(4):782-797. doi:10.1038/s41588-026-02547-5. PMID: 41912934.
7. Nava C, Cogne B, Santini A, Leitão E, Lecoquierre F, Chen Y, et al. Dominant variants in major spliceosome U4 and U5 small nuclear RNA genes cause neurodevelopmental disorders through splicing disruption. Nat Genet. 2025;57(6):1374-1388. doi:10.1038/s41588-025-02184-4. PMID: 40379786.
8. Rius R, Blakes AJM, Chen Y, De Jonghe J, Lecoquierre F, Dawes R, et al. Biallelic variants in the noncoding RNA gene RNU4-2 cause a recessive neurodevelopmental syndrome with distinct white matter changes. Nat Genet. 2026;58(4):761-773. doi:10.1038/s41588-026-02554-6. PMID: 41951959.
9. De Jonghe J, Kim HC, Adedeji A, Leitão E, Dawes R, Kajba CM, et al. Saturation editing of RNU4-2 reveals distinct dominant and recessive disorders. Nature. 2026;654(8118):429-436. doi:10.1038/s41586-026-10334-9. PMID: 41951737.
10. He H, Liyanarachchi S, Akagi K, Nagy R, Li J, Dietrich RC, et al. Mutations in U4atac snRNA, a component of the minor spliceosome, in the developmental disorder MOPD I. Science. 2011;332(6026):238-240. doi:10.1126/science.1200587. PMID: 21474760.
11. Edery P, Marcaillou C, Sahbatou M, Labalme A, Chastang J, Touraine R, et al. Association of TALS developmental disorder with defect in minor splicing component U4atac snRNA. Science. 2011;332(6026):240-243. doi:10.1126/science.1202205. PMID: 21474761.
12. Elsaid MF, Chalhoub N, Ben-Omran T, Kumar P, Kamel H, Ibrahim K, et al. Mutation in noncoding RNA RNU12 causes early onset cerebellar ataxia. Ann Neurol. 2017;81(1):68-78. doi:10.1002/ana.24826. PMID: 27863452.
13. Xing C, Kanchwala M, Rios JJ, Hyatt T, Wang RC, Tran A, et al. Biallelic variants in RNU12 cause CDAGS syndrome. Hum Mutat. 2021;42(8):1042-1052. doi:10.1002/humu.24239. PMID: 34085356.
14. Johnson MB, Russ-Silsby J, Blair PA, Govier M, Bonfield G, Domingo-Vila C, et al. Bi-allelic variants in the non-protein-coding minor spliceosome components RNU6ATAC and RNU4ATAC cause syndromic monogenic autoimmune diabetes. Am J Hum Genet. 2026;113(4):877-887. doi:10.1016/j.ajhg.2026.02.017. PMID: 41864208. Erratum: Am J Hum Genet. 2026;113(8):1776-1777.
15. Mendez R, Arriaga TM, Ma J, Bonner DE, Emami S, Levy RJ, et al. Bi-allelic RNU6ATAC variants cause a minor spliceopathy characterized by transcriptome-wide minor intron retention and multisystem manifestations. HGG Adv. 2026;7(2):100588. doi:10.1016/j.xhgg.2026.100588. PMID: 41808409.
16. Quinodoz M, Rodenburg K, Cvackova Z, Kaminska K, de Bruijn SE, Iglesias-Romero AB, et al. De novo and inherited dominant variants in U4 and U6 snRNA genes cause retinitis pigmentosa. Nat Genet. 2026;58(1):169-179. doi:10.1038/s41588-025-02451-4. PMID: 41513982.
17. Nakano Y, Suzuki H, Kuroda Y, Yoshihashi H, Okamoto N, Kondo A, et al. Whole-genome discovery of pathogenic snRNA variants and efficient extended-exome screening. iScience. 2026;29(8):116814. doi:10.1016/j.isci.2026.116814. PMID: 42502410.
18. D'Souza EN, Blakes AJ, Paluch R, Banka S, Chopra M, Coffey AJ, et al. Guidance for clinical variant classification in genes for spliceosomal small nuclear RNAs. medRxiv [Preprint]. 2026 Aug 4. doi:10.64898/2026.08.03.26359558.
19. Pejaver V, Byrne AB, Feng BJ, Pagel KA, Mooney SD, Karchin R, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. Am J Hum Genet. 2022;109(12):2163-2177. doi:10.1016/j.ajhg.2022.10.013. PMID: 36413997.
20. Miranda JP, Carvajal L, Aldana I, Alarcon P, Saavedra C. Beyond a single threshold: calibration of CADD v1.7 and AlphaGenome for ACMG/AMP computational evidence across variant classes, model versions, and reference data. bioRxiv [Preprint]. 2026 Nov 2. doi:10.64898/2026.10.31.763309.
21. Messaoud O, DiTroia S, Tarawneh R, Marten D, O'Heir E, O'Leary M, et al. Analysis of spliceosome-related coding and noncoding genes and pseudogenes reveals novel candidates. medRxiv [Preprint]. 2026 Aug 10. doi:10.64898/2026.08.06.26358951.
22. Xue J, Jiang Z, Zhang J, Li C. Cross-species analysis identifies APOBEC3-associated spliceosomal snRNA hypermutability. bioRxiv [Preprint]. 2026 Oct 13. doi:10.64898/2026.10.12.758683.
23. Ontiveros-Palacios N, Cooke E, Nawrocki EP, Triebel S, Marz M, Rivas E, et al. Rfam 15: RNA families database in 2025. Nucleic Acids Res. 2025;53(D1):D258-D267. doi:10.1093/nar/gkae1023. PMID: 39526405.
24. Nawrocki EP, Eddy SR. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics. 2013;29(22):2933-2935. doi:10.1093/bioinformatics/btt509. PMID: 24008419.
25. Grimm DG, Azencott CA, Aicheler F, Gieraths U, MacArthur DG, Samocha KE, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. Hum Mutat. 2015;36(5):513-523. doi:10.1002/humu.22768. PMID: 25684150.
26. Niroula A, Vihinen M. How good are pathogenicity predictors in detecting benign variants? PLoS Comput Biol. 2019;15(2):e1006481. doi:10.1371/journal.pcbi.1006481. PMID: 30742610.
27. Gunning AC, Fryer V, Fasham J, Crosby AH, Ellard S, Baple EL, et al. Assessing performance of pathogenicity predictors using clinically relevant variant datasets. J Med Genet. 2021;58(8):547-555. doi:10.1136/jmedgenet-2020-107003. PMID: 32843488.
28. Rastogi R, Chung R, Li S, Li C, Lee K, Woo J, et al. Critical assessment of missense variant effect predictors on disease-relevant variant data. Hum Genet. 2025;144(2-3):281-293. doi:10.1007/s00439-025-02732-2. PMID: 40113603.
29. Lu B, Liu X, Lin PY, Brandes N. Genomic heterogeneity inflates the performance of variant pathogenicity predictions. Nat Commun. 2026 [published online 2026 Oct 30]. doi:10.1038/s41467-026-78623-5.
30. Schubach M, Maass T, Nazaretyan L, Röner S, Kircher M. CADD v1.7: using protein language models, regulatory CNNs and other nucleotide-level scores to improve genome-wide variant predictions. Nucleic Acids Res. 2024;52(D1):D1143-D1154. doi:10.1093/nar/gkad989. PMID: 38183205.
31. Kuderna LFK, Ulirsch JC, Rashid S, Ameen M, Sundaram L, Hickey G, et al. Identification of constrained sequence elements across 239 primate genomes. Nature. 2024;625(7996):735-742. doi:10.1038/s41586-023-06798-8. PMID: 38030727.
32. Pollard KS, Hubisz MJ, Rosenbloom KR, Siepel A. Detection of nonneutral substitution rates on mammalian phylogenies. Genome Res. 2010;20(1):110-121. doi:10.1101/gr.097857.109. PMID: 19858363.
33. Smedley D, Schubach M, Jacobsen JOB, Köhler S, Zemojtel T, Spielmann M, et al. A whole-genome analysis framework for effective identification of pathogenic regulatory variants in Mendelian disease. Am J Hum Genet. 2016;99(3):595-606. doi:10.1016/j.ajhg.2016.07.005. PMID: 27569544.
34. Benegas G, Albors C, Aw AJ, Ye C, Song YS. A DNA language model based on multispecies alignment predicts the effects of genome-wide variants. Nat Biotechnol. 2025;43(12):1960-1965. doi:10.1038/s41587-024-02511-w. PMID: 39747647.
35. Brixi G, Durrant MG, Ku J, Naghipourfar M, Poli M, Sun G, et al. Genome modelling and design across all domains of life with Evo 2. Nature. 2026;652(8112):1349-1361. doi:10.1038/s41586-026-10176-5. PMID: 41781614.
36. Lorenz R, Bernhart SH, Höner zu Siederdissen C, Tafer H, Flamm C, Stadler PF, et al. ViennaRNA Package 2.0. Algorithms Mol Biol. 2011;6:26. doi:10.1186/1748-7188-6-26. PMID: 22115189.
37. Landrum MJ, Lee JM, Benson M, Brown GR, Chao C, Chitipiralla S, et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Res. 2018;46(D1):D1062-D1067. doi:10.1093/nar/gkx1153. PMID: 29165669.
38. Chen S, Francioli LC, Goodrich JK, Collins RL, Kanai M, Wang Q, et al. A genomic mutational constraint map using variation in 76,156 human genomes. Nature. 2024;625(7993):92-100. doi:10.1038/s41586-023-06045-0. PMID: 38057664. [gnomAD v4.1, 807,162 individuals, released 19 April 2024.]
39. Seplyarskiy V, Koch EM, Lee DJ, Lichtman JS, Luan HH, Sunyaev SR. A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription. Nat Genet. 2023;55(12):2235-2242. doi:10.1038/s41588-023-01562-0. PMID: 38036792.
40. UK Biobank Whole-Genome Sequencing Consortium. Whole-genome sequencing of 490,640 UK Biobank participants. Nature. 2025;645(8081):692-701. doi:10.1038/s41586-025-09272-9. PMID: 40770095.
41. Benoit-Pilven C, Besson A, Putoux A, Benetollo C, Saccaro C, Guguin J, et al. Clinical interpretation of variants identified in RNU4ATAC, a non-coding spliceosomal gene. PLoS One. 2020;15(7):e0235655. doi:10.1371/journal.pone.0235655. PMID: 32628740.
42. Cheng J, Novati G, Pan J, Bycroft C, Žemgulytė A, Applebaum T, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. 2023;381(6664):eadg7492. doi:10.1126/science.adg7492. PMID: 37733863.
43. Caron B, Luo Y, Rausell A. NCBoost classifies pathogenic non-coding variants in Mendelian diseases through supervised learning on purifying selection signals in humans. Genome Biol. 2019;20(1):32. doi:10.1186/s13059-019-1634-2. PMID: 30744685.
44. Benegas G, Eraslan G, Song YS. Benchmarking DNA sequence models for causal regulatory variant prediction in human genetics. bioRxiv [Preprint]. 2025 Feb 12. doi:10.1101/2025.02.11.637758.
45. Tavtigian SV, Greenblatt MS, Harrison SM, Nussbaum RL, Prabhu SA, Boucher KM, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. Genet Med. 2018;20(9):1054-1060. doi:10.1038/gim.2017.210. PMID: 29300386.
46. Sevilla-Porras M, Nieto-Molina E, Morte B, Perez-Jurado L. Further delineation of snRNA neurodevelopmental disorders: novel genetic contributors and refined phenotypic characterization. Research Square [Preprint]. 2026 Sep 29. doi:10.21203/rs.3.rs-9609995/v1.
47. Park J, Eddafir Z, Kegele J, Mei D, Smal N, Vossler-Wolf C, et al. Long-read genome sequencing resolves complex structural variants in exome-negative developmental and epileptic encephalopathy. medRxiv [Preprint]. 2026 Oct 28. doi:10.64898/2026.10.26.26365902.
48. Babazade A, Dolas Y, Duzenli T, Vural O, Furtana C, Ergun MA. Calibration of AlphaGenome Variant Impact and SpliceAI2 scores as PP3/BP4 evidence for missense, non-coding, and synonymous variants. bioRxiv [Preprint]. 2026 Oct 24. doi:10.64898/2026.10.21.760715.
49. Tao M, Marsh JA. Variant effect predictors are shaping the ClinVar classifications used to train and evaluate them. bioRxiv [Preprint]. 2026 Oct 27. doi:10.64898/2026.10.24.761512.
50. Tabet DR, Kuang D, Lancaster MC, Li R, Liu K, Weile J, et al. Benchmarking computational variant effect predictors by their ability to infer human traits. Genome Biol. 2024;25(1):172. doi:10.1186/s13059-024-03314-7. PMID: 38951922.

## Figure legends
**Figure 1.** Blind RNU2-2 receiver-operating characteristic curves for frozen v4, Rfam alternative-base rarity, phyloP447 and CADD. (A) Sealed labels (22 P/LP vs 127 L3). (B) Literature-only labels (16 vs 125).

**Figure 2.** Excess AUROC over the rarity-only null for each score, on the sealed and literature-only RNU2-2 sets, with position-only baselines.

**Figure 3.** RNU4-2 saturation genome editing. (A) Spearman correlation with functional damage (435 SNVs), by score, with v2, v3 and v4 relative to the pre-registered bar. (B) AUROC for function-depleted variants, overall and outside the critical region.

**Figure 4.** Pooled local likelihood ratio of out-of-sample v4, with bootstrap bands and PP3/BP4 evidence thresholds at prior 0.10.

**Figure 5.** Rarity-only inflation by score and benchmark: ClinVar non-coding (A9), ClinVar missense (A10), NCBoost (A11) and TraitGym Mendelian (A13).
