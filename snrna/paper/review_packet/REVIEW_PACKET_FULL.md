# Instructions for the external reviewer

You are reviewing a manuscript before submission. Act as two people at once:
- a senior referee for a Q1 human-genetics journal (Genetics in Medicine, Genome Medicine, American Journal of Human Genetics, Human Genetics, Human Mutation);
- a methods statistician who specialises in variant-effect predictor benchmarking and ClinGen PP3/BP4 calibration.

Be brutally honest. Do not take the authors' word for anything: separate what the materials demonstrate from what they merely assert. The author wants an independent rating. Their own rating has been removed from this packet so as not to anchor you.

## What you are given (below this prompt)
1. **MANUSCRIPT**: the full paper. Figures are not yet drawn; only the legends exist.
2. **SUPPLEMENT**: supplementary sections S1–S11.
3. **TIMELINE**: when each plan, model and result was frozen, with Bitcoin block times from OpenTimestamps proofs.
4. **VERIFICATION REPORT**: the author's own audit, listing what was checked, the flaws found and fixed, the scoop search, and open items. The rating section is removed.
5. **REFERENCES**: the numbered reference list.
6. **PRE-REGISTRATION** (full packet only): the frozen plan, the dated deviation log (§14), and the amendments (§15).

The code and data are in a GitHub repository named in the manuscript (folder `snrna/`). It may not be accessible to you; if it is, you may inspect it.

## Your tasks
1. **Claims.** State the paper's main claims in at most five bullets, as you understand them.
2. **Ratings.** Give two scores from 1 to 5, in steps of 0.5, each with a one-paragraph justification:
   - **Novelty**:
     - 1 = already done;
     - 2 = incremental;
     - 3 = solid new application of known ideas;
     - 4 = clearly new result or method that fills a recognised gap;
     - 5 = field-changing.
   - **Q1 potential** (likelihood of acceptance at a Q1 journal in this field, after normal revision):
     - 1 = reject everywhere Q1;
     - 2 = unlikely;
     - 3 = possible at a lower-Q1 journal;
     - 4 = likely at a good Q1 specialist journal;
     - 5 = competitive at the top journals.
3. **Flaws.** List flaws as FATAL, MAJOR or MINOR. Check in particular:
   - leakage or circularity between labels, features and the blind tests;
   - post hoc choices presented as confirmatory;
   - multiplicity;
   - whether the "rarity-only null" is logically sound and correctly interpreted;
   - whether the two blind tests are really independent;
   - whether the population-free RNU4-2 SGE analysis really answers the rarity objection;
   - whether the PP3/BP4 calibration is valid with population-observed negatives;
   - whether the pre-registration record supports the claimed order of events;
   - overclaiming in the title, abstract or conclusions.
4. **Consistency.** Look for numbers that disagree between the manuscript, the supplement and the deviation log, and for statements in one place that another place contradicts.
5. **Scoop check.** If you can search the web, look for prior or concurrent work, especially:
   - variant-effect scoring with Rfam or covariance models;
   - computational predictors or calibration for snRNA genes (RNU4-2, RNU2-2, RNU5B-1, RNU4ATAC, RNU12, RNU6ATAC);
   - benchmarks that correct for rarity or allele frequency when negatives are population variants.

   List what you searched and any paper that weakens a novelty claim. If you cannot search, say so.
6. **Improvements.** Name the specific changes or analyses that would raise each score by 0.5, separating things possible without new data from things needing new data or collaborators.
7. **Journals.** Rank the three best-fit journals. For each, predict the decision after first review (accept, minor revision, major revision or reject) and the two objections reviewers there are most likely to raise.
8. **Bottom line.** Finish with three lines:
   - "Novelty: x/5"
   - "Q1 potential: y/5"
   - "Single most important fix: …"

## Context for fair judgement
- The study has a single author with no wet-lab work. Every input is public data (ClinVar, gnomAD, UK Biobank summaries, published SGE and assay data, the literature).
- The pre-registered primary model (v2) failed its endpoints, and the paper says so. The model that succeeded (v4) was designed after that failure, then frozen and timestamped before the sealed labels were opened.
- Some items are pending because they run on a separate GPU server; they are listed in the verification report, section 6. Judge the paper as it stands, and say how those pending results could change your view.

---

# MATERIALS UNDER REVIEW
Six documents follow: manuscript, supplement, timeline, verification report, references and the full pre-registration.

---

# DOCUMENT 1: MANUSCRIPT

## RNA-family profiles provide blind-validated computational evidence for spliceosomal snRNA variants: a pre-registered study with a rarity-adjusted benchmark

**Running title:** Computational evidence for snRNA variants

**Author:** [Author name]¹ — ¹[Affiliation]. Correspondence: [email]

*Manuscript, version of 4 November 2026. Every number is traceable to the dated deviation log (`snrna/PREREG.md` §14), the result files in `snrna/results/`, and the verification report `paper/VERIFICATION_2026-11-04.md`. Citation numbers are generated from `paper/src/` by `paper/build_refs.py`.*

---

### Abstract

**Purpose.** Variants in spliceosomal small nuclear RNA (snRNA) genes are a frequent cause of neurodevelopmental disorders. Current guidance advises against computational evidence (ACMG/AMP PP3/BP4) for these genes until a tool has been calibrated on them, because genome-wide predictors perform poorly there.

**Methods.** We pre-registered and cryptographically timestamped an evaluation of variant-effect scores at 55 snRNA loci (17,385 possible single-nucleotide variants). Patient labels for RNU2-2 and RNU5B-1 were sealed before any model was fitted. A model was frozen and anchored in the Bitcoin blockchain before the labels were opened. Because negatives in small genes are population-observed variants, we scored each predictor against a score-specific rarity-only null.

**Results.** The pre-registered primary model, trained on population allele counts, failed its endpoints. A frozen model designed afterwards reached an AUROC of 0.892 on 22 sealed RNU2-2 pathogenic variants versus 127 population-observed variants (CADD 0.501). It reached 0.910 on 16 patient variants curated afterwards from the literature (Δ vs CADD +0.298 [0.153, 0.426]). Beyond rarity, the signal came from parameter-free RNA-family (Rfam) alignment profiles (excess +0.17, p = 0.001), not from CADD, ReMM, GPN-MSA or phyloP. In RNU4-2 saturation genome editing, with no population data involved, the Rfam reference-column frequency classified function-depleted variants with an AUROC of 0.80 [0.75, 0.85], against 0.49–0.60 for genome-wide scores (every paired interval excluded 0). The advantage did not extend to other RNA gene families (pre-registered negative). Provisional Supporting-strength thresholds were derived for RNU2-2, RNU4-2 and RNU4ATAC. Across 1,787 ClinVar non-coding pathogenic variants, random rare variants alone reached AUROCs of 0.62–0.68 against population-observed negatives.

**Conclusion.** RNA-family profiles provide blind-validated computational evidence for spliceosomal snRNA variants where genome-wide predictors fail. Benchmarks with population-derived negatives should report each score's rarity-only null.

**Keywords:** snRNA; RNU2-2; RNU4-2; variant interpretation; PP3/BP4; Rfam; benchmark bias; pre-registration

---

### Introduction

De novo variants in RNU4-2, which encodes the U4 small nuclear RNA (snRNA) of the major spliceosome, cause ReNU syndrome, one of the most frequent monogenic neurodevelopmental disorders (NDDs) found by genome sequencing.<sup>1,2</sup> Other snRNA disorders followed:
- dominant and recessive disorders of RNU2-2;<sup>3–6</sup>
- dominant disorders of RNU5B-1;<sup>4,7</sup>
- a recessive RNU4-2 syndrome.<sup>8,9</sup>

The minor-spliceosome snRNA disorders were known earlier: RNU4ATAC,<sup>10,11</sup> RNU12,<sup>12,13</sup> and RNU6ATAC.<sup>14,15</sup> Dominant variants in RNU4-2 and RNU6 genes also cause retinitis pigmentosa.<sup>16</sup> Diagnostic laboratories now encounter snRNA variants routinely.<sup>17</sup>

Interpretation lacks computational evidence. The proposed guidance for snRNA variant classification states that no tools yet capture features relevant to snRNA variant deleteriousness and that genome-wide tools perform poorly.<sup>18</sup> It therefore recommends using in silico evidence (PP3/BP4) only after a tool has been calibrated for snRNAs. Under thresholds calibrated on missense variants,<sup>19</sup> CADD scores the first two pathogenic RNU2-2 variants identified in the benign range.<sup>18</sup> More generally, genome-wide thresholds calibrated on missense variants give false benign evidence for pathogenic non-coding variants.<sup>20</sup> Saturation genome editing (SGE) of RNU4-2 now provides functional measurements for one gene,<sup>9</sup> but not for the others.

snRNA genes are unusual targets for genome-wide predictors. They occur in multiple paralogous copies among many pseudogenes,<sup>18,21</sup> and the spliceosomal snRNA genes are germline mutation hotspots.<sup>4,22</sup> RNA-family covariance models (Rfam) take a different approach. They describe each snRNA family's sequence and secondary structure across species and paralogs, independently of whole-genome alignment.<sup>23,24</sup>

Evaluation is a second obstacle. Benchmarks often use population variants as negatives, and rarity correlates with conservation and constraint. Such negatives can therefore flatter a predictor. This has been described for missense predictors evaluated with common benign variants<sup>25–28</sup> and for predictors evaluated across heterogeneous genomic contexts.<sup>29</sup> ClinGen calibration accordingly restricts benign variants to allele frequency below 1%.<sup>19</sup> Small genes rarely have curated benign variants, so their negatives are usually defined by population observation; here, a gnomAD allele count of at least 9. Most such negatives are rare (median allele frequency 1.8 × 10⁻⁴), yet any score that tracks rarity still separates a random unlabelled variant from them. We therefore report each score's excess over a score-specific rarity-only null, and we test the main claim on functional data that involve no population labels.

We designed the study to be falsifiable. The plan, definitions and decision rules were frozen and timestamped before analysis. Patient labels for two genes were sealed in an encrypted file before any model was fitted. Every later change is logged together with what had been seen when it was made. We report the failed endpoints alongside the successful ones.

### Materials and Methods

#### Study design and pre-registration
The pre-registration (`PREREG.md`), analysis code and each amendment were hashed and timestamped with OpenTimestamps; the proofs are anchored in Bitcoin blocks (Supplementary Table S1). Changes made after data were seen are recorded in a dated deviation log (`PREREG.md` §14), with what had been viewed at the time. Throughout, analyses are labelled as pre-registered (with amendment number), post hoc, or exploratory. Post hoc analyses whose code was committed before they were run are labelled "specified before running"; git commit times are the authoritative record.

#### Loci, variants and scores
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

#### Labels
**Positives** were pathogenic or likely pathogenic (P/LP) SNVs from three sources:
- ClinVar, with at least one review star;<sup>37</sup>
- RNUdb curation;<sup>18</sup>
- a literature tally under a frozen wording rule: a source must call the variant pathogenic, likely pathogenic or disease-causing, and any source calling it a VUS excludes it.

**Negatives (L3)** were gnomAD v4.1 genome PASS SNVs<sup>38</sup> with allele count ≥ 9, excluding positives. L3 holds 1,233 SNVs. They are population-observed but mostly rare: median allele frequency 1.8 × 10⁻⁴, with 5.4% at ≥ 1%.

**Sealing.** RNU2-2 and RNU5B-1 labels were fetched directly into an AES-256-encrypted file before any model was fitted. They were decrypted once, and every access was logged.

#### Models
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

#### Evaluation and the rarity-only null
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

#### Blind tests
**Test 1: sealed labels.** RNU2-2 had 22 P/LP vs 127 L3 SNVs. RNU5B-1 had 1 P/LP and could not be evaluated.

**Test 2: literature (amendment A8).** The design was fixed before curation. RNU2-2 and RNU5B-1 variants reported as pathogenic in the 2025–2026 literature were curated under the frozen wording rule. Curation was AI-assisted, and the curator was instructed not to look up any score. Positives absent from the sealed set and from every training file formed the test set, and the labels were committed before scoring. The pre-registered conditions were:
- v4 excess over its rarity null with p < 0.05;
- a v4 − CADD AUROC interval excluding 0.

**Test 3: prospective (amendment A12).** Units of new P/LP variants reported after the freeze are evaluated once each, when they reach a pre-registered minimum (Supplement S11).

#### Robustness analyses (post hoc, specified before running)
- **Circularity:** the ACMG/AMP criteria that each ClinVar submitter applied to each blind positive, and v4 recomputed with leave-own-position local constraint.
- **Position-only baselines:** membership of guidance-defined disease regions,<sup>18</sup> distance to them, 5′ proximity, and whether the position is observed in gnomAD.
- **Stricter negatives:** L3 restricted to AC ≥ 20, 50 or 100.

#### Functional readouts
- **RNU4-2 SGE correlation:** Spearman correlation with functional damage (−function score) for the 435 in-transcript RNU4-2 SNVs with a score.<sup>9</sup>
- **RNU4-2 SGE classification (post hoc, specified before running):** AUROC for "depleted" variants (function score < −0.302, the study's threshold) and "strong" variants (< −0.9). We used 2,000 class-stratified bootstraps and paired differences between each Rfam feature and each genome-wide score. Analyses were repeated inside and outside the guidance-defined critical region, because position alone carries much of the signal in RNU4-2.
- **RNU4ATAC splicing assay (exploratory):** published splicing efficiencies of 23 RNU4ATAC SNVs.<sup>41</sup>

No population label or negative enters the SGE analyses.

#### Cross-family test (amendment A7)
To test whether the Rfam advantage generalises, alternative-base rarity was scored zero-shot on RMRP, TERC, SNORD118 and RNU7-1, with scores frozen before labels were fetched. The primary feature had been chosen from the post hoc RNU2-2 ablation before any A7 data existed.

#### Alignment depth (exploratory, specified before running)
For each named gene we summarised phyloP447 and phyloP100 across positions:
- the 95th-percentile value (a ceiling);
- the fraction of positions sharing the gene's most frequent value;
- the per-position Spearman correlation between phyloP447 and Rfam reference-column frequency.

#### Benchmark analyses (amendments A9–A15)
- **A9:** 1,787 ClinVar non-coding P/LP SNVs against 8,666 gnomAD PASS SNVs with AC ≥ 9 within ±100 bp. The positives were intronic (1,500, sampled from 3,438), 5′UTR (148), ncRNA exon (85), 3′UTR (51) and upstream (3). The window-matched null drew, for each positive, a random non-positive SNV with AC < 9 from the same window.
- **A10:** the same design for missense variants, adding AlphaMissense.<sup>42</sup>
- **A11:** the NCBoost published benchmark (common dbSNP negatives).<sup>43</sup>
- **A13:** the TraitGym Mendelian benchmark (published matched controls).<sup>44</sup>
- **A15:** a ground-truth simulation of the null, with window-bootstrap intervals.

#### Calibration
We used the ClinGen local likelihood-ratio method, as pre-registered:<sup>19</sup>
- nearest 20% of labelled variants, at least 10 per class;
- 1,000 bootstraps;
- thresholds from the 5th (PP3) or 95th (BP4) bootstrap percentile of the local LR;
- primary prior 0.10,<sup>45</sup> with sensitivity priors 0.01, 0.0441, 0.2 and 0.5.

The pre-registered calibration (Aim 4) used the v2 framework score. Calibration of v4 was post hoc and used out-of-sample scores: frozen v4 for the blind genes, and v4 refitted without the gene's paralog cluster elsewhere. Negatives were L3. ClinVar VUS and conflicting SNVs in the named genes (279 SNVs) were scored against the thresholds.

#### Statistics and software
Analyses used Python 3 (NumPy, pandas, SciPy, scikit-learn). P-values are one-sided where the pre-registration specified a direction. The rarity-null tests on the sealed set were added after unsealing and span 12 scores; we report them unadjusted and give Bonferroni-adjusted values for the main claims. All code, inputs and logs are in the repository.

### Results

#### The pre-registered primary model failed per variant
v2 did not meet its pre-registered endpoint against RNU4-2 SGE function. Its Spearman ρ was 0.144, against a bar of 0.346 (Evo 2's ρ of 0.246 plus 0.10). Adding UK Biobank counts (v3) did not help (ρ 0.142).

v2 also failed the held-out-gene gate, with ΔAUROC vs CADD of +0.083 against a bar of 0.15. Population allele counts at these hypermutable loci therefore did not resolve individual variants. On the sealed RNU2-2 labels, v2's framework score (0.784) exceeded CADD (0.501). This came from its lab-transfer function component (0.851), not its selection component (0.659).

#### The frozen v4 model separates blind RNU2-2 patient variants from population variants
In leave-one-paralog-cluster-out validation on non-blind genes, v4 exceeded CADD by +0.186 [0.108, 0.267]. On the sealed RNU2-2 labels (22 P/LP vs 127 L3; Fig. 1A), v4 reached an AUROC of 0.892 [0.817, 0.951]. The comparators were:
- CADD 0.501 [0.393, 0.612];
- ReMM 0.501;
- phyloP447 0.724;
- GPN-MSA 0.724.

#### Most of the apparent gap is rarity; the RNA-family profile carries the rest
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

#### The result is not explained by circularity or by position alone
- **Labels.** No ClinVar submitter applied PP3 or BP4 to any blind RNU2-2 positive. 21 of 22 positives list criteria, all from PS2, PS4, PM1, PM2, PM3 and PP1. v4 reached 0.896 for the positives with listed criteria and no PP3/BP4, and 0.935 for de novo (PS2) positives.
- **Features.** v4's local-constraint feature counts each variant's own position. Recomputing it with that position left out, under the frozen coefficients, gave 0.893.
- **Position.** Four position-only baselines reached 0.728–0.752: disease-region membership, distance to it, 5′ proximity, and position unobserved in gnomAD. v4 exceeded each by +0.14 to +0.16, with every interval excluding 0. Within the guidance-defined n.1–67 region (13 vs 14 variants), v4 reached 0.780 [0.582, 0.940] and 5′ proximity 0.860, a tie within uncertainty.
- **Stricter negatives.** At AC ≥ 50, v4 reached 0.968, phyloP447 0.912, GPN-MSA 0.893 and CADD 0.669. Genome-wide conservation separates patient variants from frequent variants, but not from rarer population variants.

#### A second blind test replicates the result
Amendment A8 was frozen before curation. The frozen rule accepted 35 RNU2-2 SNVs from the 2025–2026 literature. Of these, 16 were in neither the sealed set nor any training file (14 recessive, 2 dominant).
- **Pre-registered result:** frozen v4 reached 0.910 [0.825, 0.979] against 125 L3, with excess +0.133 (p = 0.012) and Δ vs CADD +0.298 [0.153, 0.426]. Both conditions were met (Table 1).
- **Pooled with the sealed set (38 positives):** v4 reached 0.901 (excess +0.122, p = 0.001), and the Rfam features had excess +0.17 (p = 0.001). CADD (+0.016) and GPN-MSA (−0.011) had no excess.
- **Broader label rules (post hoc):** the curator's STRICT set (20 new positives) and the union of both rules (27). v4 still exceeded CADD (+0.267 and +0.237, intervals excluding 0), but its excess over the rarity null fell to a trend (+0.079, p = 0.069; +0.068, p = 0.068).

The two tests share their negatives: 125 of the sealed test's 127 L3 variants. The replication is therefore on new patient variants in the same gene, against nearly the same population variants. The sealed set included de novo dominant variants (9 with PS2), whereas the second set was mostly recessive.

#### Measured function: the family profile classifies RNU4-2 variants without population labels
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

#### The advantage is specific to the spliceosomal snRNA genes (pre-registered negative)
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

#### Rarity inflation is general and large (pre-registered, A9–A13)
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

#### Provisional calibrated evidence
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

#### Prospective test (A12)
ClinVar re-scans and a literature tally of four 2026 reports<sup>8,17,46,47</sup> have not yet brought any pre-registered unit to its minimum (Supplement S11). Units in genes never seen in training have 0–4 new positives, and the unit in the trained genes has 4 of the 5 required. A12 therefore remains descriptive, and each unit will be evaluated once when it qualifies.

### Discussion
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

### Conclusions
RNA-family alignment profiles provide blind-validated computational evidence for variants in the spliceosomal snRNA genes where genome-wide scores fail. The signal survives adjustment for rarity and classifies measured RNU4-2 function without population labels. Combined in a frozen model, the profiles yield provisional ClinGen-style Supporting thresholds. Benchmarks that use population-derived negatives should report each score's rarity-only null, because its size is large and score-specific.

### Data availability
- The scores for all 17,385 SNVs (v2 and v4, with PP3/BP4 labels) will be deposited at Zenodo **[DOI to be added]**.
- The curated literature labels (`labels/A8_literature.tsv` and `labels/A12_*.tsv`, with rejected rows and reasons) and all benchmark selections (`a9/`–`a13/`) are in the repository.
- The sealed-label event log (`sealed/SEALED_MANIFEST.jsonl`) will be released with the scores.

### Code availability
All code, the pre-registration, amendments, OpenTimestamps proofs and the deviation log are available at https://github.com/tanukusaitejesh-prog/ndd-pp3-calibration (folder `snrna/`), and will be archived with the score release.

### Ethics declaration
The study used only public, de-identified data: ClinVar, gnomAD, UK Biobank allele-frequency summaries, published variant reports and public functional data. No ethics approval was required.

### Author contributions, competing interests, funding
[To be completed by the author.] AI assistance: analysis code, literature curation support and drafting were assisted by an AI system under the author's direction; all decisions and the final text are the author's responsibility.

### References
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

### Figure legends
**Figure 1.** Blind RNU2-2 receiver-operating characteristic curves for frozen v4, Rfam alternative-base rarity, phyloP447 and CADD. (A) Sealed labels (22 P/LP vs 127 L3). (B) Literature-only labels (16 vs 125).

**Figure 2.** Excess AUROC over the rarity-only null for each score, on the sealed and literature-only RNU2-2 sets, with position-only baselines.

**Figure 3.** RNU4-2 saturation genome editing. (A) Spearman correlation with functional damage (435 SNVs), by score, with v2, v3 and v4 relative to the pre-registered bar. (B) AUROC for function-depleted variants, overall and outside the critical region.

**Figure 4.** Pooled local likelihood ratio of out-of-sample v4, with bootstrap bands and PP3/BP4 evidence thresholds at prior 0.10.

**Figure 5.** Rarity-only inflation by score and benchmark: ClinVar non-coding (A9), ClinVar missense (A10), NCBoost (A11) and TraitGym Mendelian (A13).

---

# DOCUMENT 2: SUPPLEMENT

## Supplementary Information

*RNA-family profiles provide blind-validated computational evidence for spliceosomal snRNA variants. Every number below is from the dated deviation log (`snrna/PREREG.md` §14) or from committed result files.*

### S1. Timeline of freezes and results
See `paper/TIMELINE.md` (Supplementary Table S1). The v4 code and output hashes (amendment A5) are anchored in Bitcoin block 971362 (2026-10-14 05:58 UTC). The A6 rule is anchored in block 971363 (06:04 UTC). Both precede the blind results, which were committed at 06:44 UTC.

### S2. Primary model (v2) and its endpoints
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

### S3. Robustness of the blind RNU2-2 result (post hoc, specified before running)
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

### S4. Functional readouts
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

### S5. Cross-family test (A7; pre-registered; not met)
Scores were frozen at commit 55033c3 before labels were fetched.

| Unit | Rfam alternative-base rarity | Reference-column frequency | CADD | phyloP447 | GPN-MSA | ReMM | v4 |
|---|---|---|---|---|---|---|---|
| RMRP (28 P vs 40 N) | 0.488 [0.350, 0.636] | 0.829 | 0.863 | 0.844 | 0.795 | 0.790 | 0.813 |
| SNORD118 + RNU7-1 (11 vs 134) | 0.879 [0.793, 0.948] | 0.883 | 0.816 | 0.860 | 0.875 | 0.840 | 0.923 |

TERC had 15 P vs 3 N and was not evaluable. Weighted Δ (primary − CADD) −0.251 [−0.368, −0.133]; v4 − CADD −0.005 [−0.060, 0.050].

### S5b. Alignment depth at snRNA loci (exploratory, specified before running)
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

### S6. Rarity-only null across benchmarks (A9–A15)
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

### S7. Inheritance-mode rule (A6): not informative
The frozen A6 rule (local-constraint score, threshold set on non-blind genes) gave the following on literature-annotated blind-gene P/LP variants:
- RNU2-2 (4 dominant vs 31 recessive): AUROC 0.839 [0.685, 0.968];
- both blind genes (6 vs 31): 0.892 [0.763, 0.989].

The local-constraint score counts the variant's own allele, however. In this set, 24 of 31 recessive and 0 of 6 dominant variants are observed in gnomAD. The variant's own allele count alone separates the classes (AUROC 0.887, reproduced on 4 Nov 2026), and with the own position removed the RNU2-2 AUROC is about 0.70. The rule mainly detects whether a pathogenic variant is itself seen in carriers, which is expected for recessive alleles. We do not present it as a predictor.

### S8. Calibration details
- **Method:** local LR, nearest 20%, ≥ 10 per class, 1,000 bootstraps; thresholds from the 5th/95th percentile; prior 0.10 primary, with 0.01, 0.0441, 0.2 and 0.5 as sensitivity.
- **v4 (post hoc):** out-of-sample scores, using frozen v4 for the blind genes and v4 refitted without each gene's paralog cluster elsewhere. Pooled thresholds therefore mix scores from different cross-fitted models.
- **Negatives:** L3 (population-observed), not curated benign variants. A likelihood ratio against population-observed variants can overstate evidence for a rare variant.
- **Rarity check:** `v4_rare_passrate.py`, pre-specified on 4 Nov 2026, computes for each threshold the pass rate among P/LP and among never-labelled rare SNVs. Their ratio is a lower bound on the LR against rare benign variants; the rule is "robust" if ≥ 2.08.
- **VUS:** 279 ClinVar VUS/conflicting SNVs in the named genes. With pooled v4 thresholds, 61 reach PP3 and 35 BP4; with own-gene thresholds, 37 reach PP3. CADD pooled: 73 PP3, 36 BP4.

### S9. Gene-level constraint (exploratory; replication of Chen et al. 2024 and Messaoud et al. 2026)
A within-locus 20-nt window-depletion score, with the mutation scale fitted within the locus, ranked the five NDD snRNA genes first among 55 loci with gnomAD + UK Biobank counts (AUROC 1.0, permutation p = 5 × 10⁻⁵; named loci only p = 2.5 × 10⁻⁴). Whole-locus observed/expected ratios ranked them near the bottom (AUROC 0.08–0.12), because the expressed genes are mutation hotspots. The score is post hoc, and its minimum over windows grows with locus length and observation count.

### S10. Healthy-carrier splicing readout (pre-specified; null)
1000 Genomes 30× genotypes and Geuvadis LCL junction counts (462 samples). Carriers of rare SNVs in minor-spliceosome snRNAs (n = 72) did not differ in minor-intron junction usage from non-carriers (n = 195): one-sided p = 0.29. Predicted-damaging carriers (n = 39) gave p = 0.33, and major-snRNA carriers (n = 195) p = 0.65. About 58% of samples carried a "rare" variant at these loci, which suggests unreliable short-read calls at multicopy snRNA loci.

### S11. Prospective test (A12; ongoing)
ClinVar scans (11 October and 1 November 2026) and a literature tally (Nakano 2026, Sevilla-Porras 2026, Park 2026; 29 SNVs mapped) left every pre-registered unit below its minimum. Remaining candidates: RNU2-2 n.155C>G and n.104T>G, and RNU5B-1 n.37G>C, pending the sealed-set check. Rius et al. 2026 (Sup. Table 7) lists 9 RNU4-2 P/LP SNVs, but all 9 were training labels before the freeze, so none qualifies (4 Nov 2026).

| Unit | Role | New positives | Minimum |
|---|---|---|---|
| U1 RNU6 adRP | primary | 0 | 8 |
| U2 U5 (RNU5A-1, RNU5B-1) | primary | 0 | 5 |
| U3 RNU2-2 (beyond sealed and A8) | primary | 0 | 8 |
| U4 other untrained loci | secondary | 4 (RNU4-1; flagged: the source calls interpretation limited) | 5 |
| U5 trained genes (new variants) | secondary | 4 (RNU4-2 ×3, RNU4ATAC ×1) | 5 |

Scans continue monthly. Each unit is evaluated once, when it reaches its minimum.

---

# DOCUMENT 3: TIMELINE

## Verified timeline of freezes and results (Supplementary Table S1)

Verified on 4 November 2026. Bitcoin times are the timestamps of the earliest block attesting each manifest's OpenTimestamps proof. Proofs were upgraded with `ots upgrade` on 4–4 November 2026, and block times were read from blockstream.info. Git times are commit times in UTC. A Bitcoin anchor proves the manifest existed by that time. A git time alone is weaker evidence: commit times are set by the committer, and independent confirmation rests on GitHub's push record.

| Event | Evidence | Time (UTC) |
|---|---|---|
| Pre-registration (`PREREG.md`) frozen | `prereg_manifest.txt.ots`, Bitcoin block 971237 | 2026-10-13 11:59 |
| Analysis code frozen | `code_manifest.txt.ots`, block 971241 | 2026-10-13 12:20 |
| Amendment A1 (G2 hotspot rule, simulation only) | block 971253 | 2026-10-13 14:07 |
| G2 mutation model frozen | `g2_freeze_manifest.txt.ots`, block 971255 | 2026-10-13 14:11 |
| Amendment A2 (v2 model specification) and A2b | block 971283 | 2026-10-13 18:33 |
| Amendment A4 (analysis plan for open endpoints) | block 971354 | 2026-10-14 05:07 |
| v4 code committed | git f4bd48a | 2026-10-14 05:49 |
| **A5: v4 code and DGX output hashes frozen** | `amendment_A5_manifest.txt.ots`, **block 971362** | **2026-10-14 05:58** |
| **A6: inheritance rule frozen** | `amendment_A6_manifest.txt.ots`, **block 971363** | **2026-10-14 06:04** |
| Blind (sealed) results committed | git 4081c92 | 2026-10-14 06:44 |
| A7 frozen (code) → scores frozen → evaluated | git 93682d0 → 55033c3 → 1d4086c | 09:47 → 10:10 → 10:20 |
| A8 frozen → labels committed → evaluated | git 4fdf25b → e1bb6e4 → 68710b0 | 10:23 → 14:13 → 14:16 |
| A9 frozen (positives) → evaluated | git 5b2145e → 4de7c80 | 10:32 → 14:22 |
| A10 frozen → missense evaluated | git a411288 → 13608ab | 15:16 → 17:15 |
| A11 frozen → evaluated | git e61a8f6 → 0de2637 | 10-14 15:18 → 10-15 06:13 |
| A12/A13 frozen → A13 evaluated | git cc399a8 → b9cec28 | 10-15 04:29 → 17:50 |
| A14/A15 frozen → A15 simulation run | git 35fdd9d → 4096c06 | 10-15 06:28 → 08:10 |
| A16 frozen (not yet run) | git b2818c8 | 10-15 18:07 |

### Notes
- **Sealing and decryption** of the RNU2-2 and RNU5B-1 labels are logged on the A100 machine in `sealed/SEALED_MANIFEST.jsonl`. That file must be added to the release before submission. It holds event times and hashes, not labels or the passphrase.
- **The A5 manifest records the DGX output hashes** (`v4_scores.tsv.gz`, `v4_model.json`) as 62 of 64 hex characters, transcribed from a screenshot. The full hashes are in `results_v4/v4_manifest.txt` on the DGX, and `v4_blind_eval.py` checks them before scoring. The 62-character suffix is still a unique fingerprint (248 bits).
- **The code hashes in the A5 and A6 manifests still match** the committed `v4.py`, `v4_blind_eval.py`, `inheritance_blind_eval.py` and `final_eval.py` (checked 4 November 2026).
- **Amendments A7–A16 have no OpenTimestamps proof in the repository**, so their order rests on git and GitHub push history.

---

# DOCUMENT 4: VERIFICATION REPORT (rating removed)

## Verification report, 4 November 2026

Scope: every claim in `paper/MANUSCRIPT.md` and `paper/SUPPLEMENT.md`, the integrity of the pre-registration record, the analysis code, the references, and novelty against work posted up to 4 November 2026. This report lists what was checked, how, what was wrong, and what remains open. Git commit times (UTC) are the authoritative record of order.

### 1. Integrity of the pre-registration record

| Check | Method | Result |
|---|---|---|
| Pre-registration, code, A1, G2 freeze, A2/A2b, A4, A5, A6 timestamps | `ots upgrade` + `ots verify`; Bitcoin block times (Supplementary Table S1, `paper/TIMELINE.md`) | All attested in Bitcoin blocks 971237–971363, 13–14 Oct 2026 |
| v4 frozen before blind labels were decrypted | A5 manifest in block 971362 (14 Oct 05:58 UTC) vs the decryption event log | Order holds |
| Frozen code unchanged since freeze | SHA-256 of `v4.py`, `v4_blind_eval.py` and related files vs the A5/A6 manifests | Match |
| Later amendments (A7–A16) specified before their results | git commit order of each specification vs each result | Order holds for every amendment |
| Post hoc checks of 4 Nov "specified before running" | `sge_classes.py` e2085a0 before 87fc8d8; `alignment_depth.py` b768536 before 0b2117f; `v4_rare_passrate.py` c183f40 (not yet run) | Order holds |
| Commit trailers | `git log` since 4 Nov | No Co-Authored-By lines |

Caveats recorded in the paper: the sealed-label event log (`sealed/SEALED_MANIFEST.jsonl`) is on the DGX and will be released with the scores; the A5 manifest stores output hashes truncated to a 62-character suffix.

### 2. Numbers

- **Recomputed here and matching:** RNU4-2 SGE Spearman correlations (all point estimates identical; intervals now from one saved run, `results/sge_spearman.tsv`); the A8 pattern (Rfam excess +0.16/+0.18, p = 0.001; CADD excess ≈ 0); the A10 and A11 report values; the alignment-depth Spearman (+0.90, recomputed by hand from ranks).
- **Checked against committed result files and the deviation log:** all 268 decimal numbers in the manuscript main text were traced to `paper/SUPPLEMENT.md`, `PREREG.md` §14, `results/*`, or the audited draft (automated check: 0 untraced).
- **Not recomputable here (DGX outputs):** the sealed-set AUROCs and nulls, v4 calibration thresholds, and the raw A9/A13 values. These come from DGX logs recorded in the deviation log at the time they were produced.

### 3. Code audit

An independent audit found no leakage (test suite 49/49 passing on 4 Nov, including the new cluster-bootstrap test) into the blind RNU2-2 result:
- RNU2-2 is alone in its paralog cluster;
- no training ID is in the blind genes;
- within-gene ranks use no labels;
- local-constraint self-counting changes the feature AUROC only from 0.745 to 0.748;
- the phyloP447 = GPN-MSA = 0.724 tie on the sealed set is a coincidence: the same Mann-Whitney U arises from different rank vectors.

### 4. Flaws found and fixed

| # | Flaw | Fix |
|---|---|---|
| 1 | Title and abstract said "calibrated" although the v4 calibration is post hoc with population-observed negatives | "Calibrated" dropped from the title; thresholds called provisional; the pre-specified rare-variant check is named as pending |
| 2 | L3 negatives were called "common" | Corrected: population-observed, median AF 1.8 × 10⁻⁴, 5.4% ≥ 1% |
| 3 | Claim that "40–58% of above-chance AUROC reflects rarity" | Withdrawn: the null pool is not cleaned or consequence-matched, so the share is an upper bound (A15 grid overstates it by +0.32 on average) |
| 4 | A6 inheritance rule presented as a predictor | Circular (the variant's own allele count alone gives AUROC 0.887); moved to Supplement S7 as non-informative |
| 5 | A15 window bootstrap collapsed repeated negatives and did not resample published negatives | Code fixed (`--cluster`) with a test; intervals removed from the main text; all A9–A13 rankings now descriptive; corrected intervals pending on the DGX |
| 6 | Discussion claimed the NCBoost reordering had "an interval excluding a tie" | Removed (it relied on the buggy bootstrap) |
| 7 | The A8 label-rule sensitivity table used 122 L3, not 125; it was described as differing by "≤ 0.004" | Explained in Supplement S3, with the exact differences: up to 0.007 in excess, 0.002 in Δ vs CADD |
| 8 | Sealed-set rarity tests were post hoc across 12 scores, with no multiplicity statement | Bonferroni values given: the Rfam features stay significant (0.012, 0.024); v4 does not (about 0.17); the A8 v4 test was pre-registered |
| 9 | The two blind tests share negatives (125 of 127) without saying so | Stated in Results and Limitations |
| 10 | Overclaims in the earlier draft | Fixed: "mostly dominant"; "CADD reached none" (true only in RNU2-2 and RNU4-2); "genome-wide scores at chance" (GPN-MSA reaches 0.60 on SGE); the tie inside n.1–67 (5′ proximity 0.860 vs v4 0.780) is now reported; the v4 SGE ρ of 0.382 is disclosed as computed without Evo 2 |
| 11 | Reference errors | Fixed: Greene 2026 is recessive RNU2-2; Johnson 2026 is RNU6ATAC/RNU4ATAC autoimmune diabetes; phyloP447 is cited as Kuderna 2024, not Zoonomia; Rfam 15 (Ontiveros-Palacios 2025); Evo 2 is Nature 2026; Jackson 2025 is RNU2-2 and RNU5B-1; Nava 2025 is RNU4-2, RNU5B-1 and the RNU5A-1 candidate; Lu et al. 2026 has 4 authors and is now in Nat Commun |
| 12 | Introduction cited the Roulette model (Pol III mutagenesis) for snRNA hotspots and claimed "hundreds of pseudogenes" | Cited instead: Jackson 2025 and Xue 2026 (Pol II snRNA hypermutability); "many pseudogenes", citing the guidance and Messaoud 2026 |
| 13 | Deviation-log rows and docstrings written on 4 Nov were dated "5 Nov" | Corrected text-only; erratum logged (PREREG §14); git times authoritative |
| 14 | SGE Spearman intervals differed slightly between the manuscript and the Supplement (separate bootstrap runs) | One seeded run saved to `results/sge_spearman.tsv`, cited by both |
| 15 | Rarity-matched A9 intervals in Table 3 came from the buggy bootstrap | Interval column removed; columns labelled upper bound (raw) and lower bound (rarity-matched) |

### 5. Novelty and scoop audit (to 4 Nov 2026)

**Search.**
- full bioRxiv listing 1 Sep–4 Nov 2026 (15,185 records), filtered for snRNA genes and variant-effect or calibration terms;
- medRxiv;
- Europe PMC, including full-text body searches;
- web search;
- full texts of the closest papers.

**Findings.**
- **No scoop.** No study proposes, blind-tests or calibrates an RNA-family (Rfam or covariance-model) profile as variant evidence for snRNA genes or other nuclear non-coding RNA genes. A Europe PMC body search for Rfam + snRNA + pathogenic returns 0 hits.
- **The gap is explicit.** The snRNA classification guidance (D'Souza et al., medRxiv, Aug 2026) states that "no tools have yet been developed that encompass features relevant to snRNA variant deleteriousness, and current genome-wide tools perform poorly". It calls for tools calibrated on snRNA genes and advises against PP3/BP4 until then.
- **Closest related work** (cited and positioned in the manuscript):
  - Miranda 2026: missense-calibrated CADD and AlphaGenome thresholds fail for non-coding variants;
  - Babazade 2026: AlphaGenome and SpliceAI2 calibrated in six non-coding classes, including population controls;
  - Tao & Marsh 2026: predictor evidence feeds back into ClinVar;
  - Lu et al. 2026: genomic heterogeneity inflates predictor performance;
  - Rastogi 2025: frequency-aware missense benchmark;
  - Xue 2026: germline hypermutability of snRNA genes;
  - Messaoud 2026: spliceosome gene and pseudogene prioritisation;
  - nAPOGEE: mitochondrial tRNA/rRNA, no Rfam, no nuclear genes.
- **The rarity-only null** quantifies and corrects a known caveat; it is not a new phenomenon. The manuscript says so.

### 6. Open items (DGX, morning)

1. `python v4_rare_passrate.py --blind sealed/blind_labels.tsv`. If any v4 threshold has LR_pool < 2.08, the paper reports that threshold as "not shown robust to rarity".
2. `python a15_null_validity.py ci a9 --cluster`, then the same for `a10`, `a11` and `a13:traitgym_mendelian`. These give corrected intervals for the descriptive rankings.
3. `sync_dgx_outputs.sh`. Include `results_v4/`, the A9 window files, `sealed/SEALED_MANIFEST.jsonl`, and the v4 thresholds at prior 0.0441.
4. `python a12_new_variants.py build --blind sealed/blind_labels.tsv`, to resolve the three pending RNU2-2/RNU5B-1 candidates.
5. Zenodo deposit of the scores; figures 1–5 from the synced outputs.

None of these can strengthen a claim beyond what the manuscript states. Items 1–2 can only weaken the calibration or ranking statements, which are already worded as provisional or descriptive.

### 7. Rating
[Removed from the review packet so that the external rating is independent.]


---

# DOCUMENT 5: REFERENCES

## References (numbered in order of first citation in MANUSCRIPT.md)

Generated by `paper/build_refs.py` from `paper/src/MANUSCRIPT.keys.md`. Entries are copied from `REFERENCES_VERIFIED.md` (Crossref, Europe PMC and bioRxiv APIs, 4 Nov 2026); entries marked † were added on 4 Nov 2026 after a bioRxiv API check of title, authors, date and DOI.

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
20. Miranda JP, Carvajal L, Aldana I, Alarcon P, Saavedra C. Beyond a single threshold: calibration of CADD v1.7 and AlphaGenome for ACMG/AMP computational evidence across variant classes, model versions, and reference data. bioRxiv [Preprint]. 2026 Nov 2. doi:10.64898/2026.10.31.763309. †
21. Messaoud O, DiTroia S, Tarawneh R, Marten D, O'Heir E, O'Leary M, et al. Analysis of spliceosome-related coding and noncoding genes and pseudogenes reveals novel candidates. medRxiv [Preprint]. 2026 Aug 10. doi:10.64898/2026.08.06.26358951.
22. Xue J, Jiang Z, Zhang J, Li C. Cross-species analysis identifies APOBEC3-associated spliceosomal snRNA hypermutability. bioRxiv [Preprint]. 2026 Oct 13. doi:10.64898/2026.10.12.758683. †
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
48. Babazade A, Dolas Y, Duzenli T, Vural O, Furtana C, Ergun MA. Calibration of AlphaGenome Variant Impact and SpliceAI2 scores as PP3/BP4 evidence for missense, non-coding, and synonymous variants. bioRxiv [Preprint]. 2026 Oct 24. doi:10.64898/2026.10.21.760715. †
49. Tao M, Marsh JA. Variant effect predictors are shaping the ClinVar classifications used to train and evaluate them. bioRxiv [Preprint]. 2026 Oct 27. doi:10.64898/2026.10.24.761512. †
50. Tabet DR, Kuang D, Lancaster MC, Li R, Liu K, Weile J, et al. Benchmarking computational variant effect predictors by their ability to infer human traits. Genome Biol. 2024;25(1):172. doi:10.1186/s13059-024-03314-7. PMID: 38951922.

---

# DOCUMENT 6: PRE-REGISTRATION (PREREG.md: frozen plan, deviation log §14, amendments §15)

## snRNA-VEP v2: pre-registration (week 0)

**Status.** Frozen on 13 October 2026, before any endpoint below was computed. The plan behind it is `PLAN.md`; this document turns that plan into exact, executable definitions. Where `PLAN.md` was ambiguous or internally inconsistent, the resolution is stated here and listed in section 13, so it can be checked before the analysis it affects is run.

**Timestamp.** The SHA-256 of this file and of the analysis code at freeze time are listed in `prereg_manifest.txt`, which is stamped with OpenTimestamps (`prereg_manifest.txt.ots`, Bitcoin-anchored) and committed to the repository. An OSF deposit of the same files is an author action still to do; if it happens after any analysis, the OSF record is a later copy and must not be described as predating the analysis.

**After freezing:**

- Only bug fixes, and the clarifications that section 13 explicitly allows, may change anything below.
- Each change gets a dated entry in the deviation log (section 14), saying whether any result had been seen when it was made.
- A threshold, split, label rule or test changed after its data were seen is reported as post hoc.

### 0. What had been seen before freezing

Honesty about prior exposure matters more than its absence. Before this document was written:

- **SGE table.** Supplementary Table 1 of the SGE paper was opened. Counts and categories were tabulated: 539 variants, 151 below −0.302, 72 below −0.9, 521 in-transcript, 18 ReNU-syndrome variants (12 SNVs), 353 population-observed (346 SNVs). The positions of the ReNU variants were listed; they define the critical region and are published. The paper itself (read in full by the plan's author) reports the CADD, ΔΔG and SGE AUCs. **No score was correlated with any label by this project.**
- **Crick repository** (commit `b404bf5`). Region boundaries and the published AUC code were read.
- **Access tests.** gnomAD v4.1 at RNU4-2 (3 records printed, from an access test), CADD v1.7 at one position, Roulette at two sites, gnomAD coverage at a few RNU4-2 bases.
- **Gene sets.** The GENCODE v49 gene list (111 spliceosomal genes reproduced exactly) and gene lengths.
- **RNUdb.** Only the gene list (RNU2-2, RNU4-2) and empty RNU4-2 region tracks. No variant endpoint was called for any gene.
- **Never accessed.** No RNU2-2 or RNU5B-1 patient-variant data has been downloaded, opened or tallied by this project. The two blind hotspot coordinates are known, because they are in `PLAN.md`, quoted from Table 2 of the guidance.

### 1. Data, pinned

| Input | Version / location | Use |
|---|---|---|
| RNU4-2 SGE | Nature 654:429–436 (2026), doi 10.1038/s41586-026-10334-9, Supplementary Table 1 (`41586_2026_10334_MOESM2_ESM.xlsx`, sheet `ST1_20251212_RNU4_2_merged_gf`); MD5 recorded by `snrna.py fetch` | L1 labels; paper's CADD and ΔΔG columns |
| SGE code | github.com/FrancisCrickInstitute/RNU4-2_Saturation_Genome_Editing, commit `b404bf5e98e11e0ed03aa186344b6c191591f05c` | Region map, published AUC definitions |
| gnomAD v4.1 genomes | `gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr*.vcf.bgz`, read by remote tabix | Allele counts |
| gnomAD genome coverage | gnomAD GraphQL API, `coverage(dataset: gnomad_r4).genome`, per-base median | Gene QC |
| Roulette | `N_rate_v5.2_TFBS_correction_all.vcf.bgz` (genetics.bwh.harvard.edu/downloads/Vova/Roulette), remote tabix with `.csi` | Mutation rates (`MR`) |
| CADD | v1.7 GRCh38 `whole_genome_SNVs.tsv.gz` (krishna.gs.washington.edu/download/CADD/v1.7/GRCh38), remote tabix | Baseline |
| GENCODE | v49 comprehensive annotation GTF | Gene set and coordinates |
| Reference | GRCh38 bases from the UCSC hg38 sequence API, cross-checked against the REF column of gnomAD, Roulette and CADD | REF confirmation |
| GIAB | Genome stratifications v3.3, `GRCh38_alllowmapandsegdupregions.bed.gz` | Paralog QC |
| ViennaRNA | 2.7.0, as in the SGE paper | Baseline |
| Evo 2 | `evo2_7b` checkpoint, `evo2` package version recorded at run time | G1 |
| Rfam | Release current at first download, recorded; family CMs for U1, U2, U4, U4atac, U5, U6, U6atac, U11, U12 | Consensus structure |

Restricted data (UK Biobank, All of Us, Genomics England, EGA RNA-seq) are not used. The aggregate UK Biobank and All of Us allele counts published in the CC BY 4.0 SGE table are used only to reproduce the paper's own AUCs (section 6.4). They do not enter any label set or model.

### 2. Gene universe and QC

1. **Spliceosomal genes.** GENCODE v49 genes with `gene_type "snRNA"`, on chr1–22, X or Y, whose `gene_name` matches `^(RNU1-\d+|RNVU1-\w+|RNU2-\d+|RNU4-\d+|RNU4ATAC|RNU5[A-F]-\d+|RNU6-\d+|RNU6ATAC|RNU11|RNU12|U1|U2|U4|U6|U4atac|U6atac|U11|U12)$` and does not end in `P`. This gives 111 genes, matching the guidance's count.
2. **Coverage.** gnomAD genome per-base median depth over every transcript base. Keep a gene if the median of those per-base medians is ≥ 20 and their minimum is ≥ 10. The guidance reports 87 genes. A different count is reported, not tuned.
3. **Paralog QC.** Exclude a gene if either holds:
   - more than 50% of its bases overlap the GIAB low-mappability-and-segmental-duplication union;
   - it is mostly low-VAF: summed over its PASS variants, more than 50% of heterozygous genotypes fall in allele-balance bins below 0.25 (gnomAD `ab_hist_alt_bin_freq`).

   RNU4-2/RNU4-1 is flagged as a high-risk pair but not excluded on that ground alone.
4. **Analysis set.** The genes passing 2 and 3. Every exclusion is listed with its reason.

### 3. Variants and coordinates

- **Key.** Every variant is keyed by `chrom-pos-ref-alt` in GRCh38. Transcript numbers are derived from GENCODE v49 exon coordinates and strand, and used only for display and for region definitions.
- **Universe.** All possible SNVs at every transcript base of every analysis-set gene: 3 per base.
- **Indels.** Used only in the SGE comparisons, where the SGE table supplies them.
- **REF confirmation.** Every published variant is converted to GRCh38 and its REF confirmed against the reference. A mismatch is listed and the variant dropped. A mismatch rate above 1% in any source stops the pipeline (QC failure, not a result).
- **Known-fact checks.** These must reproduce before any endpoint is computed (QC, not hypotheses):
  - 539 SGE variants;
  - 151 with score < −0.302;
  - 521 in-transcript;
  - 12 ReNU SNVs and 346 population SNVs;
  - the SGE library span chr12:120291753–120291903.

### 4. Label sets

| Set | Exact definition | Use |
|---|---|---|
| L1 SGE | `function_score_mean` (haploid) from SGE ST1; **depleted** = score < −0.302; **strong** = score < −0.9; `function_score_dip` (diploid) for sensitivity only. In-transcript = `nucleotide` 1–145 (521 variants: 435 SNVs, 86 indels) | Evaluation only |
| L2 Patient P/LP | Distinct variants (not individuals) reported as pathogenic or likely pathogenic by the reporting study, or by ClinVar with ≥ 1 review star, as of the tally date. Tagged by gene, Table 2 region and inheritance (dominant, recessive, both, unknown). RNU2-2 and RNU5B-1 entries go only into the sealed file (section 5) | Evaluation, calibration |
| L3 Benign-like | gnomAD v4.1 genomes, FILTER PASS, AC ≥ 9 (≈ AF 5.5e-5 × AN ≈ 152,000, scaling the SGE paper's rule). Sensitivity at AC ≥ 5, ≥ 20, ≥ 50. A variant in both L2 and L3 is removed from L3 and listed | Calibration |
| L4 Unlabeled | Every other possible SNV. Never called benign | Selection training |

Presence in gnomAD alone never makes a variant benign.

#### RNU4-2 regions

Transcript positions follow NR_003137.3 n.1–145, using the Crick code's `assign_region`:

| Region | Positions |
|---|---|
| Stem II | 3–16 |
| k-turn | 27–35 and 41–46 |
| Stem I | 56–62 |
| T-loop | 63–70 |
| Stem III | 75–79 |
| 3′ stem loop | 85–117 |
| Sm site | 118–126 |
| Terminal stem loop | 127–144 |
| Other | all remaining positions |

The **critical region** C is n.62–70 plus n.75–78, the span of the ReNU SNVs.

#### Blind hotspots

| Gene | Hotspot (as quoted in `PLAN.md` from guidance Table 2) |
|---|---|
| RNU2-2 | n.1–67 |
| RNU5B-1 | n.37–44 |

Section 13 allows one pre-specified correction if Table 2 itself gives different coordinates.

### 5. Splits and the sealed labels

- **Sealed genes.** RNU2-2 and RNU5B-1. Their L2 variants are compiled into `sealed/blind_labels.tsv` and encrypted with `snrna.py seal`, which uses AES-256 with a passphrase held by someone other than the analyst where possible. The SHA-256 of the plaintext and of the ciphertext is recorded in a dated amendment **before any selection model is trained**.
  - The pipeline refuses to read the file without `--unseal`, which logs the date, user and reason.
  - The file is unsealed once, for the final evaluation (section 9).
- **Other splits:**
  - leave-one-gene-out over the analysis set;
  - leave-family-out, U4 + U4atac: RNU4-1, RNU4-2, RNU4ATAC and the unnamed U4 loci held out together;
  - leave-region-out within RNU4-2, using the regions in section 4.
- **Training exclusions.** The two sealed genes are excluded from every training run until the final run. Hyperparameters are chosen only on gene-held-out folds of the unlabeled (L4 plus counts) data.

### 6. Gate G1 (week 1): is a foundation model already enough?

#### 6.1 Evo 2 score

- **Model and window.** `evo2_7b`. For each SGE variant, take the GRCh38 window of 8,192 nt centred on the variant's first reference base, chr12:c−4096 to c+4095. The alternative window is the same span with the variant applied, so indels change its length.
- **Score.** LL is the sum of per-token log-likelihoods. Δ = LL(alt) − LL(ref), computed on the forward strand and on the reverse complement of both windows, then averaged.
- **Damage score.** d = −Δ.
- **Sensitivity.** The length-normalised mean log-likelihood.

#### 6.2 Damage scores of the comparators (higher = more damaging)

| Predictor | Damage score d |
|---|---|
| CADD | v1.7 PHRED. Primary: the SGE table's `CADD_score` column. QC check: it must agree with the CADD v1.7 file to ±0.01 for SNVs |
| ViennaRNA U4/U6 ΔΔG | The SGE table's `ddG_bind` column, primary. Our 2.7.0 recomputation with the Crick parameters is a QC check |
| Allele-count floor | −log10(gnomAD v4.1 AC + 1), with AC = 0 if absent or non-PASS |
| Position-only floor | Minus the leave-one-position-out prediction of `function_score_mean` from an OLS model on two features: a one-hot of the nine RNU4-2 regions, and log(1 + distance in nt to the nearest critical-region position). All variants at a position are held out together |
| EVEE | Looked up for the ClinVar-overlapping subset only; reported, never used in a gate |

#### 6.3 Metrics

On the 521 in-transcript variants:

- **ρ.** Spearman correlation between −d and `function_score_mean`. Positive means concordant.
- **AUROC.** d for depleted versus not depleted.
- **Uncertainty.** 95% intervals from 2,000 bootstrap replicates, seed 20261013:
  - for AUROC, resampling within each class;
  - for ρ, resampling variants;
  - a position-cluster bootstrap (resampling positions with all their variants) is reported as a sensitivity analysis.
- **Also reported:**
  - SNV-only versions (435 variants);
  - per-region ρ;
  - diploid-score versions;
  - the replicate noise ceiling (pairwise Pearson r of `function_score_r1..r3`).

#### 6.4 ReNU-versus-population AUROC

- **Primary (as `PLAN.md` Appendix A specifies).** The 12 ReNU SNVs against in-transcript SGE SNVs that are not ReNU and are present in gnomAD v4.1 genomes (FILTER PASS, AC ≥ 1). The AC floor is undefined here, because controls are selected on AC, and is not reported.
- **Reproduction check.** The 12 ReNU SNVs against the paper's 346 population SNVs (`category == "UKBB/AllofUs"`, `Type == "SNV"`), using the paper's code definitions. The published values are CADD 0.65, ΔΔG 0.72 and SGE 0.93; failing to reproduce them to ±0.01 is a pipeline error.

#### 6.5 Decision

| Outcome | Rule | Consequence |
|---|---|---|
| **G1 FIRES** | Evo 2 ρ ≥ 0.60 **and** AUROC (depleted vs not) ≥ 0.85, point estimates, primary definitions | Drop the new-model framing; the paper becomes calibration plus failure analysis |
| **G1 PASSES** | Otherwise | Continue |
| Position floor flag | If the position-only floor's ρ is within 0.05 of the best non-floor predictor's ρ, the deep-model framing inside RNU4-2 ends (`PLAN.md`), whatever G1 says | Reported |
| **Pending** | If Evo 2 cannot be run, G1 stays pending | Everything else in section 6 is still reported |

### 7. Gate G2 (by week 4): is there a selection signal?

#### 7.1 Observed

For each possible SNV v in an analysis-set gene, O_v = 1 if gnomAD v4.1 genomes has v with FILTER PASS, AC ≥ 1 and AF < 0.001; otherwise O_v = 0.

Excluded from both O and E:

- sites whose gnomAD record for that allele is non-PASS (sensitivity: count them as O_v = 0);
- sites with AF ≥ 0.001;
- sites whose Roulette FILTER is `low`, or which have no Roulette record.

#### 7.2 Expected

p_v = 1 − exp(−κ_g · m_v · c_{k(v)}), where:

- m_v is Roulette `MR`;
- κ_g is a per-gene scale fitted by maximum likelihood (Bernoulli) on all included sites of gene g; it absorbs the Pol III rate, coverage and sample size;
- c_k is an element-class correction with two classes, paired and unpaired positions in the Rfam consensus structure (from cmalign of the gene to its family CM).

Two candidate mutation models:

- **M0:** c ≡ 1.
- **M1:** c_k fitted on unpaired positions outside every Table 2 region, in training genes only.

**Model choice.** Leave-one-gene-out Bernoulli log-loss over analysis-set genes, excluding the U4 family (RNU4-1, RNU4-2, RNU4ATAC and the unnamed U4 loci) and the two sealed genes. M1 is chosen only if it lowers the mean held-out log-loss by ≥ 0.002 nats per site; otherwise M0.

The chosen model is frozen, by a dated log entry with the model and its fitted c, before any U4-family or sealed-gene O/E is computed.

The site-rate covariate (number of other alleles seen at the same position) is **not** used for G2, because a position under strong selection also lacks other alleles, so the covariate would absorb the signal G2 tests. It is used only in the selection model of section 8, and as an ablation there.

**Mutation-model QC** (reported, not gating):

- observed against expected by gene and by trinucleotide context;
- a log-linear trend of O/E against position outside Table 2 regions;
- a trend p < 0.001 is flagged in the report.

#### 7.3 Windows and statistics

- **O/E.** For a set of sites S, O/E_S = Σ_S O_v / Σ_S p_v.
- **Windows.** Every run of 5 consecutive transcript positions, step 1.
- **Window depletion.** D_w = P(X ≤ O_w), where X follows the Poisson-binomial distribution with the window's p_v. Smaller D_w means more depleted.
- **Top decile.** The 10% of windows with the smallest D_w.
- **Median.** m_g is the median window O/E of gene g.

#### 7.4 Tests

**(a) RNU4-2 critical region C.**

- Primary: one-sided exact test, P(X ≤ O_C) with X ~ Poisson-binomial(min(1, m_g · p_v)) over the sites of C.
- **Recovered** if O/E_C < m_g and P < 0.01.
- Also reported, not gating: the circular-shift permutation p, which shifts C's position mask by every non-zero offset modulo the transcript length.

**(b) Each blind hotspot H.** *Superseded by amendment A1 (adopted 13 October 2026; see sections 14 and 15). The frozen rule below is still reported, as a secondary result.* **Recovered** if both hold:

1. at least one window lying wholly inside H is in the gene's top decile;
2. the circular-shift permutation p < 0.01, where:
   - the statistic S_H is the mean, over windows wholly inside H, of the window's depletion percentile (rank of D_w divided by the number of windows);
   - p = (1 + #{s ≠ 0 : S_{H+s} ≤ S_H}) / (1 + (L − 1)), for L transcript positions.

RNU2-2 (191 nt) and RNU5B-1 (117 nt) both allow p < 0.01.

#### 7.5 Decision

| Outcome | Rule | Consequence |
|---|---|---|
| **G2 PASS** | (a) recovered and both hotspots recovered | Continue with selection supervision; Aim 1 met |
| **G2 PARTIAL** | (a) recovered and exactly one hotspot recovered | Continue; Aim 1 reported as partly met; the novelty rating is not raised |
| **G2 KILL** | (a) not recovered, or neither hotspot recovered | Drop selection supervision; the paper becomes calibration-only |

**Sensitivity analyses** (reported, not gating):

- counting non-PASS records as O_v = 0;
- AF cut-off 0.01;
- M0 if M1 was chosen, and the reverse;
- a negative-binomial allele-count version, with element O/E from Σ AC.

### 8. Selection model and Aims 2–3 (weeks 5–8)

#### 8.1 Training

**Count model.** A_v ~ NegBin(λ_v, φ), with λ_v = μ_v T e^{−f(v)} and μ_v = m_v κ_g c_k, over L4 and observed sites of training genes. f(v) ≥ 0 is predicted by the encoder of `PLAN.md` Methods 3.

**Inputs:**

- sequence ±100 nt;
- Rfam profile and covariation;
- the base-pairing graph, including U4/U6 and U4atac/U6atac intermolecular pairs;
- element class;
- the site-rate covariate, computed excluding v itself.

**What training never sees.** SGE, L2 and L3 labels are never used in training, model selection or early stopping. Early stopping and hyperparameters use gene-held-out NB log-likelihood only.

**Uncertainty.** Five seeds; the score is the ensemble mean of f, and the interval is the range across seeds.

#### 8.2 Aim 2: SGE concordance

**Strict model.** Trained with the whole U4 family removed. The loose version removes RNU4-2 only and is reported as secondary.

**Met if either holds:**

- ρ (section 6.3) ≥ the best of CADD, ΔΔG and Evo 2 by at least 0.10;
- AUROC for ReNU against population (section 6.4 primary) ≥ 0.85.

**If Evo 2 is pending,** the comparison uses CADD and ΔΔG and says so.

#### 8.3 Aim 3 and gate G3 (week 6)

**Held-out sets.** Each gene is scored by a model trained without it. The week-6 set is RNU4ATAC (leave-one-gene-out and leave-family-out), plus RNU6ATAC, RNU12 and the named RNU6 paralogs (pooled as one unit, since their sequences are identical), each only if it has ≥ 10 distinct L2 SNVs. The sealed genes are not in the week-6 set (section 13, item 3).

**Metric.** Per-gene AUROC of L2 SNVs against L3 SNVs. The summary is the mean over genes weighted by the number of L2 SNVs, for the model and for CADD. ΔAUROC is model minus CADD.

**Uncertainty.** 2,000 bootstrap replicates resampling distinct variants within gene × class.

| Outcome | Rule | Consequence |
|---|---|---|
| **G3 PASS** | Weighted ΔAUROC ≥ 0.15 and the 95% interval excludes 0 | Full framework paper |
| **G3 FAIL** | Otherwise | Pivot to calibration-only or stop |
| **G3 NOT EVALUABLE** | Fewer than 10 L2 SNVs across all week-6 genes | G3 is decided at the final evaluation instead, and this is said |

**Leave-family-out.** Trained without the U4 family. RNU4ATAC ClinVar P/LP (≥ 1 star) against ClinVar B/LB, AUROC with an interval. Reported, not gating.

### 9. Final evaluation (single unsealing)

After the model, mutation model and calibration code are frozen by a dated log entry:

1. unseal `sealed/blind_labels.tsv` once;
2. score RNU2-2 and RNU5B-1 with models trained without them;
3. report each gene's AUROC (L2 against L3) and ΔAUROC against CADD, and add both genes to the Aim 3 summary as a second, blind estimate.

Nothing is refitted afterwards. Any later change is post hoc.

### 10. Aim 4: calibration (weeks 9–11)

- **Scores.** CADD, ΔΔG (U4/U6 and U4atac/U6atac only), Evo 2 if available, and the new model (cross-fitted).
- **Labels.** L2 positives, L3 negatives, as distinct SNVs. A positive is flagged if its source classification used CADD or conservation; results are shown with and without flagged positives.
- **Local LR.** For each evaluation point s on the score grid (every observed score), the local window is the nearest 20% of all labelled variants by score rank, with a minimum of 10 positives and 10 negatives. Local LR(s) = (n_P/N_P) / (n_B/N_B), with a 0.5 pseudocount per class count.
- **Bootstrap.** 1,000 replicates resampling distinct variants within class (seed 20261013).
- **Evidence thresholds.** At prior P, O_VSt(P) = (9(1−P)/P)^{4/3}. This gives 350 at P = 0.10 and ≈ 1,124 at P = 0.0441, matching Tavtigian 2018 and Pejaver 2022 to rounding. Strength k ∈ {Supporting, Moderate, Strong} needs LR ≥ O_VSt^{1/8, 1/4, 1/2}.
  - **Primary prior:** P = 0.10, the Tavtigian 2018 default.
  - **Sensitivity:** P ∈ {0.01, 0.0441, 0.2, 0.5}.
- **PP3 threshold at strength k.** The smallest grid score s* such that the 5th bootstrap percentile of local LR is ≥ the threshold for k at every grid point ≥ s*. **BP4** mirrors this: the 95th percentile must be ≤ the reciprocal at every grid point ≤ s*.
- **Aim 4 met** if the new model has a Supporting PP3 or BP4 threshold in at least one gene, or pooled across genes (pooled is stated as such). Otherwise the score is reported as "not yet usable".
- **Splits and drift.**
  - Dominant and recessive regions are calibrated separately where each has ≥ 10 positives.
  - Ancestry drift: L3 is redefined per gnomAD genetic ancestry group (AC ≥ 9 scaled to that group's AN), and thresholds and LRs are reported per group. Drift is reported, not corrected.

### 11. Secondary endpoints

**Mechanism probe.** Logistic regression on frozen encoder embeddings, three classes (dominant hotspot, recessive region, neutral), fitted on RNU4-2 and RNU4ATAC Table 2 regions only. Evaluated as:

- AUROC for dominant against recessive class on held-out positions;
- AUROC for SGE strong (< −0.9) against moderate (−0.9 to −0.302) depletion.

A research aid only; never calibrated.

**Other secondary endpoints:**

- agreement between the selection map (window O/E) and the SGE depletion map in RNU4-2: Spearman between window O/E and the window's mean SGE score;
- negative controls, each of which must give |ρ| < 0.1 and an AUROC interval containing 0.5:
  - shuffled training counts within gene;
  - a random base-pairing graph;
  - a model trained on pseudogenes only;
- pseudo-prospective: scores frozen in a dated Zenodo release, then evaluated on L2 variants first published after that date.

### 12. Analysis hygiene

- Seeds are fixed (20261013 unless stated).
- Every number in the paper is regenerated by `snrna.py` from pinned inputs.
- Each command writes a `run_info.json` with input checksums, versions and the git commit.
- Results files carry the command that made them.

### 13. Resolutions of `PLAN.md` ambiguities (check before the affected analysis runs)

1. **G2 decision when only some targets are recovered.** `PLAN.md` says G2 fails if "the known RNU4-2 depletion and two blind hotspots are not recovered". Here: PASS needs all three; PARTIAL (continue) needs RNU4-2 plus one hotspot; KILL is anything less (section 7.5).
2. **The critical-region p-value.** `PLAN.md` asks for "p < 0.01" without naming the test. Here the primary is an exact Poisson-binomial test against the gene-median-scaled expectation. A circular-shift permutation with 145 positions can reach p < 0.01 only if C is the single most depleted placement, so it is reported but not gating.
3. **Sealed genes at G3.** `PLAN.md` lists RNU2-2 and RNU5B-1 among G3's held-out genes, but also says their labels open only once, at the final evaluation. Here G3 uses non-sealed genes only; the sealed genes are evaluated once, at the end (section 9).
4. **Site-rate covariate in G2.** It is excluded from G2 for the circularity reason given in section 7.2.
5. **"ClinGen default prior".** Taken as P = 0.10 (Tavtigian 2018); Pejaver's 0.0441 is a sensitivity value.
6. **Population controls in week 1.** The gnomAD-based control set is primary, as `PLAN.md` Appendix A says. The paper's 346-SNV set is used only to reproduce published numbers.
7. **Blind hotspot coordinates.** If guidance Table 2 (read directly) gives coordinates different from `PLAN.md`'s, Table 2 wins. This is recorded as a dated clarification before G2 is run. It is the only coordinate change allowed.

### 14. Deviation log

| Date | Change | Data seen when made? | Effect |
|---|---|---|---|
| 2026-10-13 | **Bug fix (code to match frozen text).** The G1 position-floor flag was coded as floor ρ ≥ best ρ − 0.05; section 6.5 says "within 0.05". Corrected to \|floor ρ − best ρ\| ≤ 0.05; when the floor exceeds every predictor by more than 0.05, the report now adds a note. The frozen rule is unchanged | Yes: found on reading the first G1 output (floor ρ 0.570, best baseline 0.269) | The flag reads "not within 0.05" instead of firing. The case "floor far above every baseline" is not covered by the frozen rule and is left for the authors to judge |
| 2026-10-13 | **Bug fix.** Roulette has autosomes only; the reader now returns no rates for chrX/chrY instead of failing. Section 7.1 already excludes sites without a Roulette record | No endpoint data seen (week-0 site table) | One chrX U6 locus (309 SNVs) has no rates and drops out of G2 |
| 2026-10-13 | **Clarification under section 13 item 7: no change.** Guidance Table 2 was read directly (Europe PMC full text). It confirms RNU2-2 n.1–67 (recessive, 5′ region) and RNU5B-1 n.37–44 (dominant, 5′ stem loop I), and the RNU4-2 critical region n.62–70 (inclusive of insertions at n.61_62) and n.75–78. All Table 2 spans are in `table2_regions.tsv`, which is the source for "Table 2 regions" in sections 7.2, 10 and 11 | No G2 data seen | None. Note for G2: the RNU2-2 blind region is a **recessive** region, so depletion among gnomAD heterozygotes may be modest. This is the 2-fold scenario in `g2_power_simulation.txt`, where the frozen test has 6–12% power |
| 2026-10-13 | **Bug fix before any Evo 2 output existed.** `evo2_score.py` unpacked the Evo 2 forward output one level too shallow and read a version attribute the package does not export, so it would have crashed. It now calls Evo 2's own `score_sequences` (sum of per-token log-likelihoods, no BOS token), which is the section 6.1 definition. `g1` now applies the G1 rule only when all 521 variants are scored | No | None on the definition. A partial (smoke-test) scores file can no longer trigger a G1 decision |
| 2026-10-13 | **Amendment A1 adopted by the authors.** The G2 blind-hotspot rule (section 7.4b) becomes: recovered if (i) the exact test of section 7.4a, applied to the hotspot's sites, gives O/E below the gene median with p < 0.01, (ii) a window wholly inside the hotspot is in the top decile, and (iii) the circular-shift p < 0.05. The frozen rule is still reported. The critical-region test (section 7.4a) is unchanged: the optional circular-shift addition was not adopted, because simulation showed no robustness gain for it | **No.** No G2 statistic had been computed on real data; the basis is simulation only (`results/g2_power_simulation.txt`). Timestamped in `amendment_A1_manifest.txt.ots` | The decision table in section 7.5 now uses the amended rule |
| 2026-10-13 | **Clarification of M1 (section 7.2): the frozen wording is not estimable.** A two-class correction (paired, unpaired) cannot be fitted "on unpaired positions" only. Read as follows. Unpaired positions are the reference class (c_unpaired ≡ 1). c_paired is fitted by maximum likelihood, jointly with the per-gene κ, on training-gene positions outside every Table 2 region. Training genes are the analysis-set genes excluding the U4 family and the sealed genes; none of them has a Table 2 region. In leave-one-gene-out, c_paired is fitted without the held-out gene, and the held-out gene's κ is fitted on its own sites (as in M0). Family assignment for cmalign is by gene name: U1 RF00003, U2 RF00004, U4 RF00015, U4atac RF00618, U5 RF00020, U6 RF00026, U6atac RF00619, U11 RF00548, U12 RF00007. Positions in insert states, or without a consensus partner, count as unpaired | **No.** No G2 statistic computed; recorded before structure annotation and model fitting | If M1 is chosen, expected counts allow for the average stem-versus-loop difference, so G2 depletion means depletion beyond the structural class. This is the more conservative test |
| 2026-10-13 | **G2 mutation model frozen (section 7.2): M0.** Leave-one-gene-out log-loss on 42 training genes (11573 sites): M0 0.33278, M1 0.33224 nats per site. The improvement, +0.00054, is below the 0.002 bar, so M0. Fitted c_paired = 1.221 (unused). QC: per-gene O/E 0.92–1.18; per-context χ² p = 0.09 across 95 classes; position trend p = 0.67 (no flag). Structure from cmalign (Infernal 1.1.5, Rfam 15.1) | **No.** No U4-family or sealed-gene O/E had been computed. This entry and `results/g2_mutation_model.json` are committed before `g2-test` runs | G2 uses M0: p_v = 1 − exp(−κ_g m_v) |
| 2026-10-13 | **Amendment A2 adopted: specification of the selection model and its evaluation (sections 5, 8, 11).** Full text in section 15. It fixes the architecture, inputs, likelihood, hyperparameter grid, seeds, splits (paralog clusters as the holdout unit), ablations, controls, the SNV-only reading of Aim 2, the exploratory upgrades U1 and U3, and sealing as stage 0 of the A100 run | **No model had been fitted to real counts.** Known at the time: G1 and G2 results; the non-blind L2 tally (counts by gene and inheritance); fits to synthetic planted counts (tests and a CPU dry run). The dry run's evaluator printed the comparator's values on real non-blind labels: CADD AUROC 0.755 for RNU4ATAC L2 (26 SNVs) against L3 (68), and 0.880 for ClinVar P/LP against B/LB (26 vs 9). No model number on real counts existed. The sealed genes' labels had not been fetched or seen. Code and inputs are hashed in `amendment_A2_manifest.txt` (OpenTimestamps) | Sections 8 and 11 are run as A2 specifies; `train.py` refuses real counts until the blind labels are sealed |
| 2026-10-13 | **Amendment A3: evaluation-only code, comparators and literature tally** (text in section 15). Comparators (phyloP, ReMM, GPN-MSA) reported only; Aim 4 calibration, the single-unsealing final evaluation and the release packager written as code before any result; the literature tally added to L2 by a fixed rule | **No model result on real counts existed.** Comparator ReNU AUROCs were seen in a code check (listed in A3) | None on the model or the frozen bars; the evaluation is fixed as code before the results |
| 2026-10-13 | **A2b: progress logging only** in evo2_features.py, model.py, train.py, evaluate.py and run_a100.sh (which now checks the newest A2 manifest). Proof of no computational change: before and after the patch, mock Evo 2 deltas and embeddings and a synthetic-count model fit were bit-identical | No model fitted to real counts | None; hashes in `amendment_A2b_manifest.txt` (OpenTimestamps) |
| 2026-10-14 | **A4 step 1: automatic end-of-run results set aside unopened.** The v2 A100 run finished; `results/` was moved to `results_autorun_unopened/` without being opened and hashed (`autorun_unopened.sha256` on the A100). First three hashes: g1_by_region.tsv 740566c9172e2d88ffe799a6a0a707ff0aba99a6f48f6006c0bb5825bde81914; g1_evo2_sensitivity.txt 970cb9b07ad272afc81b38eb6c0b663f07ef3c1c3c1eb5c5077800d3853e8c69; g1_metrics.tsv c4a0a9f5b9e486b480eb7c9c9bc2a53849a4b64a7b1c915fc09755c248a5a056. The official v2 result is produced later by `evaluate.py` after the literature tally (A4) | No: model outputs not viewed | None on the analysis; documents the order |
| 2026-10-14 | **Bug fix (code to match the frozen L2 text, section 4 / A3) before the official evaluation.** `curate_literature.py` matched 'pathogenic' anywhere in the stated class, so rows whose per-variant class was unavailable (Rius 2026: 'class not accessible per variant (cohort: 1 pathogenic, 11 likely pathogenic, 18 VUS)') were accepted, while explicit disease-causing wording ('mutation causing MOPD1', 'cause of adRP') was rejected. Now 'causing' and 'cause of' count as disease-causing, and 'not accessible', 'unverified', 'VUS', 'uncertain' or 'benign' reject. Tally `lit/tally_2026-10-14.tsv` (AI-assisted search, LIT_TALLY_PROMPT.md; sealed genes absent) gives 40 SNVs P/LP (18 new: RNU4ATAC 11, RNU6ATAC 5, RNU12 2), 8 indels, 27 rejected (listed in labels/L2_curated_rejected.tsv). Wording that only says 'identified in individuals with NDD' or 'novel' stays rejected | No: no model output or endpoint viewed (autorun results unopened) | Labels for G3 and secondary analyses fixed before evaluation |
| 2026-10-14 | **Exploratory, post hoc: gene-level disease-gene finder (`gene_discovery.py`)**, added after the v2 endpoints were seen. It uses population counts only, no patient variant labels; disease status is public gene-level knowledge. First version: whole-locus and window O/E under one pooled mutation scale. It ranked every disease gene near the bottom (NDD AUROC 0.08-0.13), because the expressed snRNA genes are mutation hotspots (whole-locus O/E 2-3). Second version: most depleted 20-nt window relative to the locus's own fitted scale (`rel_z`). Both are kept in the output. With rel_z, gnomAD only: NDD genes AUROC 0.95 (permutation p 1e-4; named loci only 0.86, p 0.013). With gnomAD + UK Biobank: the top five loci are exactly the five NDD genes (AUROC 1.0; named loci only p 2.5e-4). The windows coincide with the reported disease regions (RNU4-2 n.62-81 vs ReNU n.62-79). Next-ranked locus: RNU4-1 n.60-79, the RNU4-2 paralogous region. adRP RNU6 genes are not enriched | Yes: all v2 endpoints viewed; G2 per-gene O/E seen earlier. Metric changed after the first run | Exploratory; reported with both versions and this history; no confirmatory claim |
| 2026-10-14 | **Exploratory analysis fixed before running: in-vivo carrier readout (`geuvadis_carriers.py`).** Uses open data only: 1000 Genomes 30x genotypes and Geuvadis LCL junction counts (462 samples). The hypotheses and readout are written in the script docstring and committed before the test is run. H1: carriers of rare SNVs (panel AF < 0.005) in RNU4ATAC, RNU6ATAC, RNU12 or RNU11 have lower minor-intron junction usage (MI - U2, z within population) than non-carriers (one-sided Mann-Whitney). H1b: the same for predicted-damaging variants (CADD >= 15 or RNU4ATAC h in the top half). N1: carriers in major-spliceosome snRNAs only show no shift (two-sided). Reported whatever the outcome | Yes for v2 endpoints; no for this readout (not computed yet) | Exploratory |
| 2026-10-14 | **Amendment A5 (v4) code written before unsealing.** `v4.py`: L2-logistic stack (C 0.1, balanced) of within-gene percentile ranks of CADD, phyloP447, ReMM, GPN-MSA, Evo 2, ViennaRNA, Rfam profile, pairing, local constraint and RNU4ATAC homology h. Trained on non-blind genes only; RNU2-2 and RNU5B-1 are excluded. `v4_blind_eval.py` is called by `final_eval.py` in the same single unsealing, after a hash check against `results_v4/v4_manifest.txt`. A CPU test without Evo 2 gave a leave-one-cluster-out weighted delta vs CADD of +0.165 [+0.087, +0.249] (RNU4-2, RNU4ATAC). One non-functional cleanup (an RNU4-2 SGE lookup that matched nothing) was made after that look. The DGX run (with Evo 2) is frozen by committing its manifest before final_eval | Yes: v2 endpoints and the v4 internal CV seen; blind labels never | v4 is a pre-specified secondary on the blind genes; v2 stays primary |
| 2026-10-14 | **A5 frozen.** The DGX v4 run (with Evo 2 and UK Biobank counts) is hashed in `amendment_A5_manifest.txt` (OTS-stamped) before any unsealing. Internal LOCO weighted delta vs CADD: +0.186 [+0.108, +0.267] | Blind labels never viewed | v4 is evaluated once, in final_eval.py |
| 2026-10-14 | **Amendment A6: inheritance-mode prediction, frozen before unsealing** (`inheritance_blind_eval.py`, run by final_eval.py; hashes in `amendment_A6_manifest.txt`, OTS). The score is local constraint from the A5-frozen v4 file. Non-blind check: dominant vs recessive L2 P/LP AUROC 1.0 (14 vs 62; within RNU4-2 14 vs 11); conservation does not separate them. Blind test: AUROC and threshold accuracy on the RNU2-2/RNU5B-1 P/LP SNVs with an inheritance annotation. Caveat: in RNU4-2 the two classes occupy different regions | Non-blind inheritance labels seen; blind never | New secondary endpoint, evaluated once |
| 2026-10-14 | **Single unsealing done (final_eval.py; logged in sealed/SEALED_MANIFEST.jsonl).** Blind genes RNU2-2 (22 P/LP vs 127 L3) and RNU5B-1 (1 P/LP vs 56 L3; not informative). v2 framework S vs CADD, weighted over blind genes: +0.267 [+0.134, +0.399] (RNU2-2: S 0.784, F 0.851, f 0.659). Combined with the week-6 units: +0.154 [+0.066, +0.241]; this is the second, blind estimate, and the G3 decision stays as made at evaluate.py. **v4 (A5, hash verified):** RNU2-2 AUROC 0.892 [0.817, 0.951] vs CADD 0.501 [0.393, 0.612]; weighted delta +0.374 [+0.262, +0.481]. **A6:** not evaluable. The blind P/LP SNVs carry no dominant/recessive annotation except the single RNU5B-1 variant (gene default dominant; predicted correctly). The v4 and A6 scripts were run by hand on the already-decrypted labels while final_eval's calibration step was still running (same frozen code) | Blind labels now viewed (once) | Final; nothing refitted |
| 2026-10-14 | **Final Aim 4 calibration completed** (the calibrate.py step of final_eval.py had stopped; re-run with the decrypted labels, same code, no refit). At prior 0.1, the model reaches PP3 Supporting in RNU2-2 (blind; F >= 0.402, CADD none), RNU4-2 (S >= 4.79), RNU4ATAC (PP3 S >= 0.241 and BP4 S <= -0.046) and the dominant pool (f >= 4.16). No Moderate at prior 0.1. **Exploratory, post hoc:** v4 feature ablation on blind RNU2-2: every drop-one model 0.875-0.908, best single features alt rarity 0.860 and Rfam ref_col_freq 0.847, CADD 0.501. v4 refitted without the RNU4 cluster vs RNU4-2 SGE (435 SNVs): Spearman 0.382 (CPU, without Evo 2; the Aim 2 bar was 0.346); Rfam ref_col_freq alone 0.459. Neither replaces a pre-registered endpoint | Yes | Exploratory |
| 2026-10-14 | **v3 (exploratory) results**: the same v2 architecture retrained on gnomAD + UK Biobank counts (6,324 vs 4,100 observed SNVs). Aim 2 strict rho 0.142 [0.045, 0.234] (v2 0.144). Loose 0.251 (v2 0.261). RNU4 cluster held out 0.252 (v2 0.179). ReNU primary AUROC 0.478 (v2 0.669). G3 RNU4ATAC: f 0.687 (v2 0.646), S 0.859 (v2 0.852), delta vs CADD +0.090 (v2 +0.083). About 54% more observed variants did not improve per-variant resolution of the selection model; the gene-level and family-profile signals carry the results | Yes | Exploratory |
| 2026-10-14 | **Exploratory, post hoc: circularity check of the blind RNU2-2 result (`circularity_check.py`), label side done before any re-scoring.** All ClinVar submissions for the 22 blind RNU2-2 P/LP SNVs were fetched (VCV XML) and the ACMG criteria each submitter applied were read. **No submitter applied PP3 or BP4 to any of them.** 21/22 list criteria (all from PS2, PS4, PM1, PM2, PM3, PP1; mostly Leitão et al. 2026 via Essen); 1 (UDN, chr11-62841662-G-T) states none. The two in-silico keyword hits are mechanism prose, not evidence codes. Residual overlap to report: every positive has PM1 (critical region) and 17/22 have PM2 (absent/rare in gnomAD), which overlap the L3 definition and v4's local-constraint feature; 4 clean positives lack PM2 and 9 have PS2 (de novo). Planned re-scoring (v4, S, f, F, CADD, comparators): label sets all / no PP3-BP4 / criteria listed / without PM2 / PS2 only; plus v4 with leave-position-out local constraint under the frozen coefficients (each L3 variant otherwise counts towards its own window). Evidence table in `sealed/blind_evidence.tsv` (not committed) | Yes (blind labels already unsealed) | Exploratory; AUROCs to be added when run on the DGX |
| 2026-10-14 | **Further post hoc checks, code written before running:** `blind_robustness.py` (position-only baselines on blind RNU2-2: Table 2 region membership, distance to it, 5′ proximity, per-position gnomAD observation; within-region AUROC; negatives at AC ≥ 20/50/100 and at low-mutability sites) and `v4_calibrate.py` (PP3/BP4 for out-of-sample v4: frozen scores for the blind genes, cluster-cross-fitted v4 elsewhere; pooled and per gene; CADD calibrated the same way; ClinVar VUS/conflicting SNVs in named snRNA genes, `labels/clinvar_vus_snrna.tsv`, scored against the thresholds). `release.py` now adds the hash-checked v4 score and its pooled evidence label | Yes | Exploratory |
| 2026-10-14 | **Rarity-only null added to `blind_robustness.py` before it was run on real labels.** A code check with randomly drawn pseudo-positives (rare RNU2-2 SNVs, no patient labels) showed that L3 (gnomAD AC ≥ 9) is separable from any rare SNV by scores that track rarity: on random pseudo-positives GPN-MSA reached ~0.79, the Rfam alt rarity ~0.68 and phyloP447 ~0.69, CADD ~0.54. So AUROC against L3 has a score-specific null well above 0.5. The check reports, for every score, the real blind AUROC, its rarity-only null (1,000 draws of as many random rare SNVs as positives) and the excess; the excess, not the raw AUROC, is the honest effect size | No real label used for the null | Exploratory; framing of every L2-vs-L3 AUROC in the paper |
| 2026-10-14 | **Post hoc results on blind RNU2-2 (22 P/LP vs 127 L3; frozen v4 hash verified).** *Circularity:* v4 with leave-position-out constraint 0.893 (frozen 0.892; reproduction max diff 1e-4); criteria-listed positives without PP3/BP4 (21) 0.896; PS2 de novo subset (9) 0.935; PM2-free (4) 0.947. *Position-only baselines:* Table 2 region 0.740, distance 0.728, 5′ proximity 0.752, position unobserved 0.743; v4 minus each +0.14 to +0.16, all intervals exclude 0. Inside n.1–67 (13 vs 14): v4 0.780 [0.582, 0.940], 5′ proximity 0.860 (tie). *Rarity-only null (1,000 draws):* excess over null — Rfam alt rarity +0.175 (p 0.001), Rfam ref column freq +0.172 (p 0.002), v4 +0.108 (p 0.014; null mean 0.783), local constraint +0.137 (p 0.021), phyloP447 +0.038 (p 0.27), CADD −0.033, ReMM −0.018, GPN-MSA −0.072. *Stricter negatives:* AC ≥ 50: v4 0.968, phyloP447 0.912, GPN-MSA 0.893, Rfam alt rarity 0.912, CADD 0.669; AC ≥ 100 similar. *v4 calibration (out-of-sample):* pooled 99 P vs 1,117 B, AUROC v4 0.788 vs CADD 0.746; PP3 Supporting v4 ≥ 0.362, BP4 Supporting ≤ −1.35; per gene PP3 Supporting in RNU2-2 (≥ 0.216), RNU4-2 (≥ 0.189), RNU4ATAC (≥ 0.88), CADD none in RNU2-2/RNU4-2. *VUS (279 SNVs):* v4 pooled PP3 61 / BP4 35; own-gene thresholds PP3 37; CADD pooled PP3 73 / BP4 36 (CADD has no within-gene threshold in RNU2-2 or RNU4-2). **Reading:** the patient-specific signal beyond rarity sits in the Rfam family features; v4 adds nothing over Rfam alt rarity alone (+0.031, interval crosses 0); "v4 0.89 vs CADD 0.50" overstates the gap because the two scores have different rarity nulls (0.78 vs 0.53). Genome-wide conservation separates patient variants from very common variants (AC ≥ 50) but not from moderately common ones | Yes | Exploratory; the paper reports excess over the rarity null as the effect size and leads with the family profile, not the stack |
| 2026-10-14 | **A7 evaluated (`a7_ncrna.py eval --acmg`, scores frozen at 55033c3 before labels): NOT MET.** RMRP (28 P vs 40 N): primary Rfam alt rarity 0.488 [0.350, 0.636]; Rfam ref column freq 0.829; CADD 0.863; phyloP447 0.844; GPN-MSA 0.795; ReMM 0.790; v4 frozen 0.813. SNORD118 + RNU7-1 (11 vs 134): primary 0.879 [0.793, 0.948]; ref column freq 0.883; CADD 0.816; phyloP447 0.860; GPN-MSA 0.875; ReMM 0.840; v4 frozen 0.923. TERC: 15 P vs 3 N, below the 5-per-class minimum, so not evaluable; the decision used the 2 evaluable units (pre-registration said 3 units; this is disclosed, not changed). Weighted ΔAUROC primary − CADD −0.251 [−0.368, −0.133], − phyloP447 −0.250 [−0.362, −0.134]; v4 frozen − CADD −0.005 [−0.060, +0.050]. Sensitivities agree (consensus-only, AC ≥ 20, no-PP3/BP4 positives). **Reading:** the family-profile advantage does not generalise to these RNA genes; genome-wide scores perform well on single-copy RNA genes (RMRP, SNORD118), so their failure is specific to the spliceosomal snRNA loci. The primary feature (alt rarity, chosen from the RNU2-2 ablation) fails on RMRP while ref column frequency does not; this is reported, and no post hoc switch of the primary is claimed | Yes (A7 labels, once) | Pre-registered negative; the cross-family claim is dropped. Motivates an exploratory test of why genome-wide scores fail at snRNA loci (paralog copy number, alignment, hypermutability), with the A7 genes as single-copy contrasts |
| 2026-10-14 | **Amendment A8 recorded before curation** (text in section 15; `a8_literature.py`, mapping checked on a synthetic row only) | No literature variant of the blind genes curated | Second blind test of v4 |
| 2026-10-14 | **Amendment A9 recorded; positives selected and committed before any score was fetched** (text in section 15) | No A9 score viewed | Genome-wide test of benchmark inflation |
| 2026-10-14 | **A8 labels curated and mapped, committed before scoring.** Two AI-assisted curations (RNU2-2: 105 variants from Greene 2025/2026, Jackson 2025/2026, Leitão 2026, Krygier 2026, Nakano 2026, Henry 2026, Messaoud 2026, Sevilla-Porras 2025/2026 and others; RNU5B-1: 7 variants; no RNU5A-1 variant called P/LP anywhere). Curators were told not to look up any score. **Bug fix in `a8_literature.py map` before any scoring:** the source-coordinate check misparsed `chr11:62841806 C>T`-style coordinates and rejected 18 correctly mapped rows; it now parses chromosome and position with a regex (checked on all four formats). Frozen wording rule applied unchanged: 37 SNVs accepted (RNU2-2 35: 4 dominant, 31 recessive; RNU5B-1 2), 6 indels, 69 rejected, mostly because one source called the variant VUS (strict by design). Of the RNU2-2 SNVs, 16 are not in the sealed blind set (14 recessive, 2 dominant), above the pre-registered minimum of 8. For transparency, 28 of the 35 RNU2-2 acceptances are also tagged STRICT by the curator and 13 STRICT SNVs were rejected by the frozen rule; the frozen rule stays primary | Counts and IDs only; no score of any A8 variant | A8 primary is evaluable; A6 becomes evaluable |
| 2026-10-14 | **A8 evaluated (`a8_literature.py eval`, frozen v4 hash verified): MET.** New literature-only RNU2-2 positives (16: 14 recessive, 2 dominant) vs 125 L3: v4 0.910 [0.825, 0.979], rarity null 0.777, excess +0.133 (p 0.012); v4 − CADD +0.298 [+0.153, +0.426]. Rfam ref column freq 0.819 (excess +0.150, p 0.010); Rfam alt rarity 0.826 (+0.145, p 0.010); phyloP447 0.787 (+0.104, p 0.054); CADD 0.613 (+0.081, p 0.090); GPN-MSA 0.859 (+0.066, p 0.128). Sealed + new pooled (38 vs 125): v4 0.901, excess +0.122 (p 0.001), v4 − CADD +0.350 [+0.251, +0.448]; Rfam features excess +0.166/+0.168 (p 0.001); CADD +0.016; GPN-MSA −0.011; phyloP447 +0.067 (p 0.057). RNU5B-1: 1 new positive, not evaluable. **A6 (frozen score and threshold) on literature-annotated blind-gene P/LP SNVs:** RNU2-2 4 dominant vs 31 recessive, AUROC 0.839 [0.685, 0.968], threshold accuracy 0.857; both genes 6 vs 31, AUROC 0.892 [0.763, 0.989], accuracy 0.865 | Yes (A8 labels, once) | Second blind test replicates the first; A6 evaluated for the first time (small dominant class) |
| 2026-10-14 | **A9 fetch made robust (no analysis change).** The first DGX fetch hung on one remote read inside htslib, where no Python timeout reaches. `a9_inflation.py fetch` now runs each ±100 bp window as its own task, abandons a task after 180 s, terminates the pool after each round and retries unfinished windows (finished ones are cached) for up to 3 rounds; windows still failing are skipped and listed in `a9/fetch_info.json`. Positives, windows, scores and the evaluation are unchanged | No A9 score viewed | Operational only |
| 2026-10-14 | **A9 evaluated: H1 MET, H2 NOT MET.** 1,787 ClinVar non-coding P/LP SNVs vs 8,666 gnomAD PASS AC ≥ 9 SNVs within ±100 bp. Pooled raw AUROC / rarity-only null / inflation / corrected: CADD 0.828 / 0.632 / +0.132 / +0.196; GPN-MSA 0.805 / 0.676 / +0.176 / +0.129; ReMM 0.730 / 0.618 / +0.118 / +0.112 (phyloP447 not readable remotely, skipped as pre-specified). Inflation per class: intronic (1,500 P) +0.124 to +0.188; 5′UTR (148) +0.110 to +0.123; 3′UTR (51) +0.099 to +0.139; ncRNA exon (85) +0.148 to +0.194. Rarity alone accounts for 40% (CADD), 58% (GPN-MSA) and 51% (ReMM) of the pooled above-chance AUROC. Tool rankings are unchanged after correction. **Reading:** benchmarks that use population-common variants as negatives measure selection/constraint signal and pathogenicity together; the rarity-only null separates them. The null is not pure noise: rare variants are themselves enriched for deleterious ones, so the corrected AUROC is a conservative estimate of pathogenicity-specific signal | Yes (A9, once) | Genome-wide; the snRNA inflation is one instance of a general effect |
| 2026-10-14 | **Exploratory, post hoc analyses written before running (`extras.py`).** (1) *A8 label-rule sensitivity:* the frozen rule rejects any variant one source calls VUS; A8 is re-run with the curator's STRICT rows (41 RNU2-2 SNVs) and with the union of both rules (48), new positives only, same metrics. (2) *phyloP447 in A9:* phyloP from the UCSC REST API per window (the bigWig was unreadable remotely), then the frozen A9 evaluation on the augmented table. (3) *Mechanism (descriptive):* per gene with ≥ 5 positives and ≥ 5 L3, fraction of possible SNVs observed in gnomAD and mean Roulette rate against the AUROC of CADD, phyloP447, GPN-MSA and the Rfam features, with Spearman across genes (few genes; reported as description) | A8 and A9 results seen | Exploratory |
| 2026-10-14 | **phyloP447 added to A9 (UCSC API, coverage 1.00; frozen evaluation re-run).** Pooled raw / null / inflation / corrected: phyloP447 0.776 / 0.647 / +0.147 / +0.128 (CADD, GPN-MSA, ReMM unchanged to ±0.001). Per class phyloP inflation +0.126 to +0.160. H1 unchanged (met). With phyloP in the panel the evaluator reports H2 met: in intronic variants (1,500 P) phyloP447 moves above GPN-MSA after correction (0.133 vs 0.132) and in 5′UTR ReMM above phyloP447 (0.102 vs 0.099). Both reorderings are within overlapping intervals; per the pre-registered H2 wording they are reported as ties, not as a finding. The substantive observation is that GPN-MSA, the most rarity-sensitive score (inflation +0.175), loses most | Yes | Exploratory completion of the pre-specified phyloP comparator |
| 2026-10-14 | **A8 label-rule sensitivity and mechanism table (exploratory).** *A8:* new RNU2-2 positives vs 122 L3 under the frozen rule (16), curator STRICT (20) and union (27): v4 − CADD +0.300 [+0.166, +0.429], +0.267 [+0.148, +0.387], +0.237 [+0.124, +0.349]; v4 excess over the rarity null +0.126 (p 0.009), +0.079 (p 0.069), +0.068 (p 0.068); Rfam ref column freq +0.148 (p 0.018), +0.096 (p 0.051), +0.084 (p 0.049). The advantage over CADD is robust to the label rule; the excess over rarity is significant only under the frozen rule and a trend under broader rules, which add weaker labels. *Mechanism:* across 8 genes, the fraction of possible SNVs observed in gnomAD did not track genome-wide score failure (Spearman with CADD AUROC +0.21; phyloP447 +0.48). CADD fails in RNU2-2 (0.55), RNU4-2 (0.52) and RNU12 (0.45) but not RNU4ATAC (0.77), RNU6ATAC (0.71), RMRP (0.86), SNORD118 (0.78) or RNU7-1 (0.85). **Hypermutability does not explain the failure; the cause is untested.** The manuscript no longer attributes it to multicopy/hypermutable loci | Yes | Two claims softened |
| 2026-10-14 | **Amendment A10 recorded; missense positives and ClinVar benign sets selected and committed before any A10 score was fetched** (text in section 15) | No A10 score viewed | Extends A9 to missense and SpliceAI; tests two fixes |
| 2026-10-14 | **Amendment A11 recorded before any A11 score was fetched** (NCBoost benchmark tables frozen in `a11/`; text in section 15) | No A11 score viewed | Shows the effect on a published benchmark |
| 2026-10-14 | **A10 missense evaluated: H1 MET.** 1,484 ClinVar P/LP missense SNVs with an AlphaMissense score vs 1,592 gnomAD-common missense SNVs in ±100 bp windows (pool 369,504 rare/absent missense SNVs; ClinVar B/LB 1,078). Raw / null / inflation / corrected / rarity-matched / vs ClinVar B/LB: AlphaMissense 0.920 / 0.762 / +0.262 / +0.158 / 0.694 / 0.916; CADD 0.813 / 0.654 / +0.154 / +0.159 / 0.690 / 0.857; GPN-MSA 0.870 / 0.802 / +0.302 / +0.068 / 0.555 / 0.827; phyloP447 0.721 / 0.621 / +0.121 / +0.100 / 0.622 / 0.725. H3: the rarity-matched AUROC lies closer to corrected + 0.5 than the raw AUROC in 100% of rows. ClinVar B/LB negatives do not remove the inflation (as stated in advance: they are largely frequency-defined). The rarity pool contains damaging rare variants, so rarity-matched and corrected AUROCs are lower bounds on pathogenicity-specific signal. SpliceAI step failed on the environment (setuptools ≥ 70 lacks pkg_resources); fixed by pinning setuptools < 70 and re-run | Yes (A10 missense, once) | The rarity effect extends to missense and to AlphaMissense |
| 2026-10-14 | **Exploratory, post hoc: rarity-free readout.** Spearman with RNU4-2 SGE damage (−function score, 435 SNVs; no population label involved): Rfam ref column freq +0.459 [0.376, 0.533] (known since the v4 SGE rematch), Rfam alt rarity +0.389 [0.305, 0.464], GPN-MSA +0.143, CADD +0.088, ReMM +0.062, phyloP447 −0.024 (1,000 bootstraps). The family profile carries functional signal on measured function where genome-wide scores do not, consistent with the blind RNU2-2 excess over the rarity null | Yes | Exploratory; second readout within the snRNAs |
| 2026-10-14 | **Amendment A7 proposed and frozen before any A7 score or label was computed** (text in section 15). Seen beforehand: ClinVar P/LP SNV counts per candidate gene only (RMRP 28, TERC 15, SNORD118 6, RNU7-1 5; RPPH1, RN7SK, RN7SL1, SNORD116-1, SNORD3A, SCARNA2, VTRNA1-1, MIR140 0; MIR184 1; MIR96 0), used to choose the units. The primary score was chosen from the post hoc v4 ablation on blind RNU2-2 (best single feature); the A7 genes have never been scored by any model in this project | No A7 score or label viewed | Confirmatory test of the family-profile hypothesis on untouched genes |
| 2026-10-14 | **A7 scores frozen before any A7 label was fetched.** `a7_ncrna.py score` ran on the DGX (Infernal from bioconda, ViennaRNA 2.7.0): 2,763 SNVs in RMRP, TERC, SNORD118 and RNU7-1, with 97%, 94%, 97% and 98% of positions at Rfam consensus columns. The run manifest is copied to `a7_freeze_manifest.txt` (the a7_ncrna.py hash matches commit 90701ea) and pushed to GitHub before `a7_ncrna.py eval`; an OpenTimestamps stamp of `a7/a7_manifest.txt` is added on the DGX | No A7 label viewed | Freeze record for A7 |
| 2026-10-15 | **Amendment A12 recorded before any A12 label was built or scored** (text in section 15; `a12_new_variants.py`, checked on synthetic data only, `tests/test_a12_a13.py`). Prospective blind tests of the frozen v4 on patient variants not in any label file committed by the A5 freeze, in five pre-registered units (three primary: RNU6 adRP genes, U5 genes, RNU2-2 beyond the sealed and A8 sets). Seen beforehand: the label files already in the repository (gene counts: no training label in RNU5A-1, RNU6-2, RNU6-7, RNU11, RNU4-1, RNU5D/E/F-1, RNVU1-*); the A8 curator's note that no RNU5A-1 variant was called P/LP. No ClinVar scan for A12 and no literature tally had been run | No A12 label or score viewed | Prospective test; each unit evaluated once, when it reaches its minimum |
| 2026-10-15 | **Amendment A13 recorded before any A13 score was fetched** (text in section 15; `a13_benchmarks.py`, synthetic check in `tests/test_a12_a13.py`; `raritybench.py` summary). Seen beforehand: TraitGym Mendelian file columns (chrom, pos, ref, alt, OMIM, consequence, label, tss_dist, match_group) and class counts (338 vs 3,042); no score or allele frequency of any TraitGym variant | No A13 score viewed | Second published benchmark for the rarity-only null |
| 2026-10-15 | **A12 first ClinVar scan (fileDate 2026-10-11, the release current at the freeze) and A13 benchmark frozen.** Scan: 397 SNV records in 12 loci; P/LP ≥ 1 star not in any earlier label file: 0 outside RNU2-2, 3 in RNU2-2 (still to be checked against the sealed set with `build --blind`); no VUS reclassified yet. So the ClinVar arm starts empty, as expected one day after the freeze; it is re-scanned monthly. Literature note (abstracts only, no variant list curated): the recurrent RNU6 adRP variants are insertions (n.55_56insG, n.56_57insG; Quinodoz 2026), which A12 does not score (SNVs only), so U1 may stay below its minimum; this is recorded now and the unit definition is not changed. TraitGym frozen: `a13/traitgym_mendelian/benchmark.tsv` (338 vs 3,042; md5 10e527e23718ebd7530b625eea22ecc2, LF line endings); no score fetched | Counts only; no score of any A12 or A13 variant | Baseline for the prospective scans |
| 2026-10-15 | **Bug fix in `a13_benchmarks.py` before any A13 score (operational only).** The freeze md5 was computed on a Windows checkout with CRLF line endings, so a Linux checkout of the same committed table would fail the fetch-time check. The md5 now ignores line endings and the table is written with LF; re-frozen from the same source file: content unchanged, md5 10e527e23718ebd7530b625eea22ecc2. New `a13_benchmarks.py` sha256 7c052f0ea03a2900ea3b92d9a63029cce28fb89530d66491cb0996efc2484a8c (the A12/A13 manifest holds the earlier one) | No A13 score viewed | None on the analysis |
| 2026-10-15 | **A11 evaluated (`a11_ncboost.py eval`, DGX; results transcribed from the DGX terminal into `results/a11_*`, originals stay on the DGX): H1 MET.** NCBoost published benchmark, 655 pathogenic vs 6,494 published common negatives (GRCh38; pool 205,825 rare/absent SNVs). Raw / null / inflation / corrected / rarity-matched: CADD 0.825 [0.805, 0.845] / 0.708 / +0.208 / +0.117 / 0.660; ReMM 0.817 [0.799, 0.837] / 0.774 / +0.274 / +0.043 / 0.616; GPN-MSA 0.837 [0.820, 0.853] / 0.812 / +0.312 / +0.025 / 0.573; phyloP447 0.772 [0.749, 0.794] / 0.695 / +0.195 / +0.077 / 0.610. Inflation on a published benchmark is larger than in A9 (+0.20 to +0.31 vs +0.12 to +0.18). **H2 (descriptive):** the published order ReMM > CADD does not hold even on raw AUROC (0.817 vs 0.825, overlapping intervals: a tie); GPN-MSA is first on raw AUROC and last after correction and on rarity-matched AUROC; no interval was computed for corrected AUROC, so the reorderings are reported as descriptive. Rarity alone accounts for 64% (CADD), 86% (ReMM), 93% (GPN-MSA) and 72% (phyloP447) of the above-chance raw AUROC | Yes (A11, once) | The rarity effect is shown on a published benchmark; it is large enough to reorder tools |
| 2026-10-15 | **A10 status.** The DGX `results/a10_report.txt` contains the missense arm only (as logged on 14 October); the SpliceAI arm (H2) has not produced output yet (no SpliceAI section in the report). H2 stays open until `spliceai_parse` and `eval` are re-run after SpliceAI finishes | No SpliceAI score viewed | A10 H2 pending |
| 2026-10-15 | **Amendments A14 and A15 recorded before any A14 variant was selected or scored and before the A15 simulation or any interval was computed** (text in section 15; `a14_recalibration.py`, `a15_null_validity.py`; synthetic checks in `tests/test_a14_a15.py`). Seen beforehand: A9, A10 (missense) and A11 results as logged above, including AlphaMissense AUROC 0.916 against ClinVar B/LB. Also added `sync_dgx_outputs.sh` to put the DGX-only outputs into git (operational) | No A14 score, no A15 output | Clinical consequence of the rarity effect (A14); validity and precision of the correction (A15) |
| 2026-10-15 | **Bug fix in `a15_null_validity.py` before any A15 output was seen.** The summary read `m.T` (pandas transpose) instead of the column `T`, so the first run stopped before printing or saving anything. Fixed to `m["T"]`; no other change. New sha256 a59c48552422eb061cda9b71a358c344f5442f660ce3f1edcbcf82c072b39523 | No | None |
| 2026-10-15 | **A15 simulation run: E1 FAILS, E2 and E3 HOLD.** Known truth T = AUROC of pathogenic vs rare benign. (1) Raw AUROC against common controls overstates T whenever the score tracks rarity (by up to +0.42 when the score has no pathogenicity signal; +0.04 to +0.08 when the signal is strong). (2) E1 fails because when the true signal is strong and the rare pool is 20% damaging, the rarity-matched AUROC M understates T by more than raw overstates it (e.g. beta 2, gamma 0.5, pi 0.2: raw +0.040, M -0.084); M is a lower bound, diluted by pi. (3) E3 holds: deattenuated M = (M - pi/2)/(1 - pi) recovers T (mean absolute error 0.0001) when pi is known. (4) **Corrected AUROC (raw - null) + 0.5 is biased downward whenever there is real signal and the score also tracks rarity (mean -0.088, as low as -0.360): it over-corrects and is not an estimate of T.** E2 holds: with no signal, M and corrected + 0.5 are at 0.5 while raw reaches 0.92. **Reading:** the paper's effect size for pathogenicity-specific discrimination should be the rarity-matched AUROC, reported as a lower bound and deattenuated over a plausible pi range; corrected AUROC is kept only as a test statistic for "any signal beyond rarity" (it is 0 when there is none), not as a magnitude. This changes how A9-A13 are reported, not what they found (inflation estimates are unaffected) | Yes (A15 simulation) | Reporting of A9-A13 changes: rarity-matched (deattenuated) AUROC becomes the magnitude; corrected AUROC is no longer presented as one |
| 2026-10-15 | **A15 window-bootstrap intervals (200 replicates, 50 null draws; DGX; numbers transcribed from the DGX screen, files to follow via `sync_dgx_outputs.sh`).** Rarity-matched AUROC [95% CI]. **A9 (1,787 non-coding):** CADD 0.704 [0.695, 0.711], GPN-MSA 0.659 [0.646, 0.670], phyloP447 0.655 [0.643, 0.666], ReMM 0.633 [0.622, 0.639]; CADD exceeds every other tool (paired differences +0.045 to +0.071, all intervals exclude 0); GPN-MSA vs phyloP447 a tie. **A10 (1,484 missense):** AlphaMissense 0.694 [0.685, 0.704], CADD 0.691 [0.681, 0.702], phyloP447 0.622 [0.609, 0.635], GPN-MSA 0.555 [0.545, 0.567]; **AlphaMissense - CADD +0.003 [-0.011, +0.016]: the 0.107 raw-AUROC advantage of AlphaMissense over CADD (0.920 vs 0.813) disappears among rare variants.** **A11 (NCBoost, 655):** CADD 0.660 [0.638, 0.676], ReMM 0.616 [0.587, 0.630], phyloP447 0.610 [0.585, 0.631], GPN-MSA 0.573 [0.551, 0.593]; GPN-MSA, first on raw AUROC (0.837), is last after rarity matching (CADD - GPN-MSA +0.087 [+0.067, +0.102]), a reordering whose interval excludes 0; ReMM vs phyloP447 a tie (+0.006 [-0.018, +0.018]) although raw AUROC separates them by 0.045. Every corrected AUROC interval excludes 0, so every tool carries signal beyond rarity. **Reading:** the reorderings are real, not ties; rarity matching removes AlphaMissense's lead over CADD on ClinVar missense | Yes (A15 intervals) | Rankings in published and ClinVar benchmarks depend on rarity; reported with intervals |
| 2026-10-15 | **A13 evaluated (TraitGym Mendelian; DGX; transcribed from the screen, files to follow): H1 MET.** 338 causal non-coding variants vs 3,042 published matched controls; 97% of controls are gnomAD-common, so H1 applies; pool 80,376. Raw / null / inflation / rarity-matched / within-match-group: CADD 0.893 [0.871, 0.915] / 0.744 / +0.244 / 0.717 / 0.887; ReMM 0.897 / 0.817 / +0.317 / 0.660 / 0.893; GPN-MSA 0.905 [0.886, 0.922] / 0.827 / +0.327 / 0.663 / 0.906; phyloP447 0.848 / 0.706 / +0.206 / 0.697 / 0.852. **H2:** GPN-MSA is first on raw AUROC and third after rarity matching (order CADD > phyloP447 > GPN-MSA > ReMM); the same reversal as A11, now on a benchmark built by GPN-MSA's own authors. Within-match-group AUROC equals raw AUROC (matching on consequence and TSS distance does not remove the rarity contrast). Intervals pending (`a15 ci a13:traitgym_mendelian`) | Yes (A13, once) | Second published benchmark: inflation +0.21 to +0.33, ranking reverses |
| 2026-10-15 | **A14 evaluated (DGX; transcribed; files to follow). H1 formally SUPPORTED but not informative about rarity.** ClinVar missense, ALL: 5,921 P/LP vs 5,745 B/LB (AlphaMissense scored); RARE (gnomAD genome AC <= 5): 5,727 vs 2,842. AlphaMissense AUROC 0.967 in both sets. Prior 0.0441 thresholds, ALL -> RARE: AlphaMissense PP3 Supporting 0.627 -> 0.657, Moderate 0.841 -> 0.858, Strong 0.985 -> 0.991; CADD Supporting 25.1 -> 25.1, Moderate 27.0 -> 26.6, Strong none in both; GPN-MSA Supporting 8.23 -> 8.39, Moderate 9.77 -> 10.1. Published AlphaMissense thresholds: Supporting and Moderate met in both sets; **Strong (>= 0.972) not met in either set** (minimum 5th-percentile local LR 23.0 in ALL and 18.2 in RARE vs 33.6 needed; ~1% of B/LB score above 0.972 in both). H2: under RARE thresholds 8.4% of rare P/LP lose one AlphaMissense strength level, 7.7% for GPN-MSA, 0% for CADD (4.1% gain). Validity check: our ALL-set AlphaMissense thresholds (0.627 / 0.841 / 0.985) do not reproduce the published ones (0.792 / 0.906 / 0.972); the local-LR procedure and data differ from Bergquist et al. **Reading:** H1 is met as worded, but the Strong failure occurs in the ALL set too, so it reflects the precision of our procedure at high scores (few benign variants there), not a rarity effect; we do not claim that the published Strong threshold is miscalibrated. The rarity-specific shift with curated rare B/LB is in the predicted direction but small. Contrast to report: AlphaMissense separates P/LP from curated rare B/LB at AUROC 0.967 but from random rare missense variants at 0.694 (A10/A15; deattenuated 0.70-0.74), so curated benign sets, even rare ones, are much easier than the rare-variant population (consistent with the stated caveat that rare B/LB may rest on BP4) | Yes (A14, once) | Clinical-threshold claim not made; the curated-benign vs rare-population gap is the finding |
| 2026-10-15 | **Amendment A16 recorded after the A14 results and before any A16 window was sampled or scored** (text in section 15; `a16_population_calibration.py`; synthetic checks in `tests/test_a16_tables.py`). Motivated by A14: curated rare B/LB are far easier than random rare variants (AlphaMissense AUROC 0.967 vs 0.694), and A14's "not met" verdicts reflected precision, so A16 uses a three-way verdict (demonstrably met / demonstrably not met / undetermined) and an assumption-free break-even pi. Primary pi (0.02) fixed now, before data. Also added `paper_tables.py` (tables and forest plot generated from result files; no estimation) | A14 seen; no A16 data | A16 is a follow-up designed after A14 and is labelled as such |
| 2026-11-04 | **Exploratory, post hoc: RNU4ATAC functional readout (Benoit-Pilven et al. 2020, PLoS One, doi:10.1371/journal.pone.0235655, S3 'Results cell.assay').** 23 single SNVs with U12-reporter splicing efficiency in TALS fibroblasts, all mapped (numbering as in the curated tally). Spearman with damage (-efficiency): homology h +0.745 (n 11, p 0.008), GPN-MSA +0.581 (p 0.004), CADD +0.472 (p 0.023), ReMM +0.461, phyloP447 +0.382, Rfam ref_col_freq +0.344 (p 0.11), Rfam alt rarity +0.168 (p 0.44), ViennaRNA +0.147. 17 of the 23 are training labels, so v4 is not evaluated on it. This is consistent with A7 and the RNU4ATAC blind results: in RNU4ATAC the genome-wide scores work, and the Rfam advantage is confined to RNU2-2, RNU4-2 and RNU12. It is reported as a limit of the Rfam claim | Yes | Exploratory |
| 2026-11-04 | **A12 literature tally 1 (analyst-curated from full texts; frozen wording rule).** Sources: Nakano 2026 iScience (doi:10.1016/j.isci.2026.116814, Table 1), Sevilla-Porras 2026 Research Square (doi:10.21203/rs.3.rs-9609995/v1, Table 1), Park 2026 medRxiv (doi:10.64898/2026.10.26.26365902, Table 1; n.73C>T and n.119G>A are VUS and were not entered). 29 SNVs mapped, 0 rejected (`lit/a12_batch_2026-11-04.tsv`). After removing anything in earlier label files: U5_trained 4 (RNU4-2 3, RNU4ATAC 1), U4_other 4 (all RNU4-1 from one paper, whose Discussion calls the gene unestablished; flagged). Candidates still unseen pending the sealed-set check on the DGX: RNU2-2 n.155C>G and n.104T>G; RNU5B-1 n.37G>C. Every unit is below its minimum, so no evaluation. No score was looked at | No | A12 not yet evaluable |
| 2026-11-04 | **Exploratory check written before running: do the v4 PP3/BP4 thresholds hold against random rare variants? (`v4_rare_passrate.py`)** The v4 calibration used L3 negatives (gnomAD PASS AC ≥ 9). These are population-observed but mostly rare: median AF 1.8e-4, and 5.4% with AF ≥ 1%. Per unit (pooled, RNU2-2, RNU4-2, RNU4ATAC), with the frozen prior-0.10 Supporting thresholds and the same out-of-sample scores, the check computes the fraction of P/LP, L3 and never-labelled rare SNVs that pass. LR_pool = pass_P / pass_pool is a lower bound on the LR against rare benign variants, because the pool contains damaging variants. Rule fixed now: 'robust to rarity' if LR_pool ≥ 350^(1/8) = 2.08 for PP3 (≤ 1/2.08 for BP4), else 'not shown'. To be run on the DGX | No: no output, v4 scores not available here | Qualifies the clinical-threshold claim |
| 2026-11-04 | **Independent code audit before submission; two bug fixes and four reporting corrections** (post hoc; no number recomputed here unless stated). Tests 48/48 pass. No leakage into the blind RNU2-2 result was found: RNU2-2 is alone in cluster C03; no training ID is in the blind genes; within-gene ranks are label-free; local-constraint self-count changes the feature AUROC .745 → .748; the phyloP447 = GPN-MSA = 0.724 tie on the sealed set is a genuine coincidence (same Mann-Whitney U from different rank vectors; columns distinct, Spearman 0.39). **(1) Bug in the A15 window bootstrap:** `metrics` resampled positives with multiplicity but collapsed repeated negatives and pool variants with np.unique, and never resampled published negatives (A11/A13). The intervals logged on 15 Oct are therefore too narrow. Fix: `a15_null_validity.py ci <x> --cluster` gives each SNV one home window and resamples windows with multiplicity and published negatives independently (test added; point estimates identical). Ranking claims are descriptive until it is re-run. **(2) The 'frozen rule' row of the A8 label-rule sensitivity (extras.py) is not a reproduction of A8:** it removes the curator's STRICT IDs from the negatives for all rules (122 vs 125 L3), and the seed differs (effects ≤ 0.004). The A8 primary result (125 L3) stands. **Reporting corrections:** (a) L3 is population-observed, not 'common' (median AF 1.8e-4; 5.4% ≥ 1%); (b) '40-58% of above-chance AUROC reflects rarity alone' is withdrawn. The A9/A10 null pool was not cleaned of unsampled P/LP or exonic SNVs in the windows and is not consequence-matched, so the null is an upper bound on rarity's contribution (the A15 grid shows the share overstated by +0.32 on average). (c) A11/A13 nulls draw from the positives' windows while the published negatives lie elsewhere, so their inflation combines rarity with local sequence context. (d) A6 is largely circular: local constraint counts the variant's own allele, 24/31 recessive vs 0/6 dominant A8 variants are observed in gnomAD, and the variant's own AC alone separates them (AUROC 0.887); with the own position removed, RNU2-2 0.70 (gnomAD counts). A6 moves to the Supplement as non-informative. Also stated in the paper: A8 shares 125 of 127 negatives with the sealed test (new positives, same gene and negatives); the sealed-set rarity-null tests are post hoc across 12 scores (v4 p 0.014, about 0.17 after Bonferroni), whereas the A8 v4 excess test was pre-registered | Yes | Paper text corrected; A15 intervals to be re-run with --cluster on the DGX |
| 2026-11-04 | **Exploratory test written before running: population-free classification on RNU4-2 SGE (`sge_classes.py`).** Labels: SGE 'depleted' (function < −0.302; secondary 'strong' < −0.9) among the 435 in-transcript RNU4-2 SNVs. No population-derived label or negative is involved. Scores: Rfam reference-column frequency and alternative-base rarity (parameter-free), CADD, GPN-MSA, phyloP447, ReMM, ViennaRNA ddG. AUROC with 2,000 class-stratified bootstraps, paired differences Rfam − genome-wide, all / critical region / outside. Expected from the known Spearman values (Rfam 0.39–0.46 vs ≤ 0.14): Rfam exceeds every genome-wide score. Reported whatever the result | SGE Spearman correlations already seen; AUROCs not computed | Rarity-free evidence for the family-profile claim |
| 2026-11-04 | **SGE classification result (sge_classes.py, as specified): population-free evidence for the family profile.** RNU4-2 SGE 'depleted' (106 vs 329 in-transcript SNVs): Rfam ref column freq AUROC 0.802 [0.753, 0.848], alt rarity 0.733; CADD 0.564, GPN-MSA 0.604, phyloP447 0.489, ReMM 0.513, ViennaRNA 0.629; Rfam ref − each genome-wide score +0.197 to +0.313, all intervals exclude 0. 'Strong' (43 vs 392): 0.865 vs ≤ 0.636. Outside the critical region (75 vs 321): 0.774 vs ≤ 0.572, intervals exclude 0. Inside the critical region (31 vs 8): ties. **Observation leading to a mechanism test, written before running (`alignment_depth.py`):** phyloP447 is constant (0.687) across RNU4-2 n.60–82, one of the most conserved stretches of U4 snRNA; a column conserved across hundreds of aligned mammals would score far higher. Hypothesis, stated in the draft discussion before this check: whole-genome alignments capture few species at some snRNA loci (paralogous copies), so conservation-based scores are capped there. Descriptive test per named gene: phyloP447 and phyloP100 ceiling (95th percentile over positions), the fraction of positions sharing the modal value (flat runs), and the per-position Spearman between phyloP447 and Rfam ref column freq. Prediction: low ceilings and weak phyloP–Rfam agreement in RNU2-2, RNU4-2 and RNU12 (where CADD fails), higher in RNU4ATAC and RNU6ATAC (where it works). Few genes; descriptive only | Yes (SGE values; phyloP flat run noticed in that output) | Mechanism candidate for why genome-wide scores fail at some snRNA loci |
| 2026-11-04 | **Alignment-depth check result (alignment_depth.py, as specified; descriptive).** RNU4-2: phyloP447 95th percentile 0.90 (max 0.90); 59% of positions share one value (0.687); per-position Spearman with Rfam reference-column frequency −0.01. Its tandem paralog RNU4-1 reaches 8.03 (agreement 0.32). phyloP100 (vertebrate multiz) is similar for both: 3.49 and 3.61. So the 447-way alignment gives RNU4-2 almost no informative conservation, consistent with orthologous sequence being aligned to the paralog. Other genes: RNU2-2 ceiling 4.46 (agreement 0.39); RNU4ATAC 7.67 (0.55); RNU6ATAC 7.75 (0.40); RNU12 7.71 (0.36). RNU12 does not fit the low-ceiling pattern, although CADD fails there (8 positives). Across the 5 genes with CADD AUROCs, agreement of phyloP447 with Rfam tracked CADD performance (Spearman +0.90); the ceiling did not (+0.20). Descriptive only; not confirmed against the alignment itself because the UCSC API was unreachable | Yes | Candidate mechanism for RNU4-2, reported as such |
| 2026-11-04 | **Date erratum (text only).** The five rows above from the v4 rarity check onward, and the docstrings of `v4_rare_passrate.py`, `sge_classes.py`, `alignment_depth.py` and the A15 bug-fix note in `a15_null_validity.py`, were first written with the date '5 Nov 2026'. Every one of these commits is dated 4 Nov 2026 UTC (commits c183f40 to 0b2117f, 18:02–18:35 UTC). The dates were corrected to 4 Nov; no other character changed (git diff). Git commit times are the authoritative record; the order 'specified before run' holds for every row (e2085a0 before 87fc8d8; b768536 before 0b2117f) | n/a | None on any result |
| 2026-11-04 | **A12 U5_trained: Rius et al. 2026 checked; no new positive.** Rius 2026 Nat Genet Sup. Table 7 lists 9 RNU4-2 SNVs classified P/LP (n.119A>G, n.7G>A, n.8C>G, n.11A>C, n.32G>A, n.45G>C, n.46G>A, n.118A>G, n.122T>G). All 9 are already in `labels/L2_nonblind.tsv` (training labels before the A5 freeze), so none is eligible for A12. U5_trained stays at 4 of 5; U4_other 4 of 5; primary units 0. No score was viewed | No | A12 remains descriptive-only |
| 2026-11-04 | **Scoop and novelty audit (Sep–Nov 2026 preprints; bioRxiv full scan 1 Sep–4 Nov, medRxiv, Europe PMC, web).** No study proposes or calibrates an RNA-family (Rfam/covariance-model) profile as variant evidence for snRNA or other nuclear non-coding RNA genes (Europe PMC body search Rfam + snRNA + pathogenic: 0 hits). The snRNA classification guidance (D'Souza et al., medRxiv 2026) states that 'no tools have yet been developed that encompass features relevant to snRNA variant deleteriousness, and current genome-wide tools perform poorly' and calls for tools calibrated on snRNA genes. Related, complementary work to cite: Miranda et al. 2026 (CADD/AlphaGenome thresholds calibrated on missense fail for non-coding pathogenic variants), Babazade et al. 2026 (AlphaGenome and SpliceAI2 calibrated in six non-coding classes, including population variants as controls), Tao & Marsh 2026 (VEP evidence feeds back into ClinVar labels; calls for classification-independent benchmarks), Lu et al. 2026 (genomic heterogeneity inflates VEP performance), Xue et al. 2026 (APOBEC3-associated germline hypermutability of Pol II spliceosomal snRNA genes), Messaoud et al. 2026 (spliceosome gene and pseudogene prioritisation), nAPOGEE (mitochondrial tRNA/rRNA predictor; no Rfam, no nuclear genes). The rarity-only null remains a new quantification and correction of a known caveat (closest: Rastogi et al. 2025; Lu et al. 2026) | No | Novelty claims scoped accordingly in the manuscript |
| 2026-11-04 | **Reproduction and wording correction (no new analysis).** (1) `sge_spearman.py` re-ran the RNU4-2 SGE Spearman correlations with a fixed seed (20261104) so that the manuscript and supplement cite one saved file (`results/sge_spearman.tsv`); point estimates are identical to the earlier runs, and interval bounds differ by ≤ 0.007 from the two earlier unsaved runs. ViennaRNA intermolecular ΔΔG (ρ 0.237) is now listed. (2) The audit row above described the extras.py 'frozen rule' row as differing from the primary A8 run by '≤ 0.004'; the exact differences are 0.007 in excess (+0.126 vs +0.133) and 0.002 in Δ vs CADD (+0.300 vs +0.298). Supplement S3 states the exact values | Yes (same data as before) | None on any conclusion |

### 15. Proposed amendments (the authors decide before the affected analysis runs)

**A1, the G2 blind-hotspot test (section 7.4b). ADOPTED 13 October 2026 (hotspot rule only; see the deviation log).** Proposed 13 October 2026, before any G2 statistic was computed on real data, on the basis of simulation only (`g2_power.py`, `results/g2_power_simulation.txt`).

- **The problem.** The frozen rule (top-decile window plus circular-shift p < 0.01) has low power, because with 117 or 191 placements, p < 0.01 needs the hotspot to be the single best placement. With a correct mutation model, a hotspot depleted 50-fold is recovered in 57% of simulations (8 nt in 117), and a 67-nt region depleted 2-fold in 12% (191 nt).
- **The exact test alone is not the fix.** The exact Poisson-binomial test of section 7.4a has far more power (97% and 70%), but trusts the mutation model. Under a deliberately misspecified, spatially autocorrelated rate error, its false-positive rate rises to 2–4%.
- **Proposed hybrid.** A hotspot is recovered if (i) the exact test of section 7.4a gives p < 0.01, (ii) at least one window lying wholly inside the hotspot is in the top decile, and (iii) the circular-shift p < 0.05.
  - False-positive rate ≤ 1.5% under both mutation models.
  - Power: 95%/90% (8 nt, 50-fold, correct/misspecified model) and 52%/27% (67 nt, 2-fold).
- **For consistency.** The critical-region test (section 7.4a) could add circular-shift p < 0.05; this costs little power (100% → 100% at 20-fold depletion, 66% → 65% at 3-fold).

Adopting A1 means a dated deviation-log entry stating that no G2 data had been seen.

**A2, the selection model and its evaluation (sections 5, 8, 11). ADOPTED 13 October 2026, before any model was fitted to real counts** (the authors asked for the model to be built and run as specified here). The code is the specification; `amendment_A2_manifest.txt` hashes it. Anything not changed here stays as frozen.

*A2.1 Inputs (replaces the input list of section 8.1).*

- **Gene universe.** The 55 analysis-set loci of `model_inputs/genes.tsv.gz`, built by `model_features.py`.
- **Graph per gene.** Nodes are the target snRNA's nucleotides plus partner RNAs:
  - U4: U6;
  - U6: U4 and U2;
  - U2: U6 and the branch-site consensus UACUAAC;
  - U1: the 5′ splice-site consensus CAGGUAAGU;
  - U4atac: U6atac;
  - U6atac: U4atac and U12;
  - U12: U6atac and UCCUUAAC;
  - U11: AUAUCCUUU;
  - U5: none.

  Partner sequences are RNU6-1, RNU4-2, RNU2-1 and the minor-spliceosome genes. Edges are the backbone, intramolecular pairs (Rfam 15.1 consensus via cmalign) and intermolecular pairs (RNAcofold MFE dimer, ViennaRNA 2.7.0).
- **Node features.** Base; molecule; Rfam profile (cmemit, 20,000 samples) and entropy; insert state; stack length; distance to helix end; paired flags; relative position; Sm motif; loop class.
- **Variant features.** Alternative base; transition; CpG; consensus column and pair frequencies; intermolecular pair broken or wobble; ΔΔG of folding; ΔΔG of binding to each partner (maximum and sum); reference and alternative pairing class.
- **Gene features.** log κ_g (from the gene's own counts; gene-level only, so it cannot reorder variants within a gene); Rfam bits per column; insert fraction; length ratio; HGNC-named flag; family.
- **Evo 2 features** (`evo2_features.py`, evo2_7b, 8,192-nt window centred on the gene):
  - layer `blocks.28.mlp.l3` embeddings at transcript positions, forward and reverse-complement averaged, reduced by PCA to 64 dimensions over all target positions (label-free);
  - per-SNV Δ log-likelihood, the mean over strands, as in section 6.1. Reference and alternatives are scored in identically shaped batches.

  "Sequence ±100 nt" in section 8.1 is replaced by these features plus base identity.
- **Site-rate covariate.** It enters only the offset (β·s_v), never f.

*A2.2 Likelihood (section 8.1 made concrete).*

- **Model.** AC_v ~ NegBin(μ_v, r), with log μ_v = α_g + log m_v + β s_v − γ_g f(v), where m_v is the Roulette rate. α_g is profiled out per gene by the Poisson moment estimator; this replaces κ_g T, and the structural correction c_k is M0 (none), as frozen in section 14.
- **Counts.** gnomAD v4.1 genomes FILTER-PASS AC.
- **Sites.** Those with the G2 `include` flag and a Roulette rate. Genes with no observed site carry no within-gene information and are skipped.
- **Selection and gate.** f(v) ≥ 0 (softplus). The gate γ_g ∈ [0, 1] is a logistic function of the gene features.

*A2.3 Encoder.*

- L layers, each:
  - message passing over typed edges (GINE-style);
  - global self-attention with learned per-head biases for edge type and shortest-path distance (capped at 16);
  - a feed-forward block.
- LayerNorm and dropout throughout; 4 heads.
- The variant head is an MLP on [node state, variant features]. Its output layer is zero-initialised, so f starts constant (absorbed by α_g) and any variation in f is learned from the counts. Without this, a CPU dry run on synthetic counts showed that a barely trained network echoes its input features (shuffled-counts control ρ −0.23), which would make the negative controls uninterpretable.

*A2.4 Training and model selection* (no SGE, L2 or L3 label is used for f; section 8.1 holds).

- **Optimisation.** Full batch; AdamW at lr 3·10⁻³; gradient clip 5; up to 600 epochs.
- **Early stopping.** NB log-likelihood is evaluated every 5 epochs on validation clusters (15% of the fit's training clusters, drawn per seed), with patience 12 evaluations; the best state is restored.
- **Hyperparameters.** Grid d ∈ {64, 128} × layers ∈ {2, 3} × dropout ∈ {0.1, 0.3} × weight decay ∈ {10⁻⁴, 10⁻²}. The configuration is chosen by the mean held-out NB log-likelihood in 3-fold cluster cross-validation on the strict training set (fold shuffle seed 20261013, model seed 11).
- **Seeds.** 11, 22, 33, 44 and 55. The score is the ensemble mean of f; the range across seeds is reported. GPU scatter-adds are not bit-reproducible, so seeds differ slightly between reruns.

*A2.5 Splits.*

- **Holdout unit.** Paralog clusters are the unit: single linkage at ≥ 90% identity within a family, giving 46 clusters. This is stricter than leave-one-gene-out. For example, RNU5B-1 is held out with RNU5A-1 and RNU5E-1 (C06), and RNU4-2 with RNU4-1 (C04).
- **Fits:**
  - strict (no U4 family);
  - loose (no RNU4-2; secondary, as frozen);
  - one leave-cluster-out fit for each of the 14 clusters containing a named gene. These give every named gene a cross-fitted score, and the sealed genes' final scores.
- **Sealed genes.** They are in no training set.

*A2.6 Aim 2 (section 8.2): SNVs only.*

- **Population.** The model scores SNVs, so Aim 2 uses the 435 SGE SNVs. The bar compares against the baselines' SNV values from G1: Evo 2 0.246, ΔΔG 0.237, CADD 0.088. Aim 2 is met if the strict model's ρ ≥ 0.346 or its ReNU AUROC (section 6.4 primary) ≥ 0.85.
- **Also reported:**
  - the loose fit and the C04 leave-cluster-out fit;
  - per-region ρ;
  - the position-only floor, as a supervised reference outside the bar.

*A2.7 Ablations and controls (section 11).*

- **Ablations.** Twelve strict-set ablations: no partner molecules; no intramolecular pairs; neither; no covariation features; no Evo 2 embeddings; no Evo 2 Δ; no Evo 2; no gate; no site-rate covariate; Poisson likelihood; no attention; MLP only. Each reports Aim 2 metrics and the held-out NB log-likelihood gain on the U4-family genes over the mutation-only null.
- **Controls.** Shuffled counts within gene; a random graph (pair and intermolecular edges rewired within gene, seed 7); pseudogenes only (the strict set's unnamed loci). The frozen criterion applies: |ρ| < 0.1 and an AUROC interval containing 0.5.
- **Score resolution.** Ensemble f is rounded to 0.01 (log units; about 3% of the spread of f in a fitted model, and a 1% change in expected count) before any metric, so a fit that learned nothing reads as tied (AUROC 0.5; undefined ρ counts as 0) instead of being ranked on optimiser noise.
- **Expectation, stated now.** The random-graph control keeps conservation and Evo 2 inputs, so it may not collapse. Its result is reported as is, and the contribution of structure is read from the ablations.

*A2.8 G3 (section 8.3).*

- **Scoring.** Each unit is scored by the fit that held out its cluster.
- **Units:**
  - RNU4ATAC;
  - RNU6ATAC;
  - RNU12;
  - the named RNU6 paralogs, pooled.

  A unit is eligible with ≥ 10 distinct L2 SNVs.
- **Labels.** L2 is `labels/L2_nonblind.tsv`, built by `labels.py nonblind` on 13 October 2026 from ClinVar P/LP with ≥ 1 star, RNUdb curated P/LP and the SGE paper's ReNU SNVs, plus any later curated additions in a dated file. L3 is gnomAD PASS AC ≥ 9 minus L2.
- **Comparator.** CADD for non-RNU4-2 genes is the v1.7 lookup.
- **Rules.** Metric, bootstrap and decision rule are unchanged.

*A2.8b The G3 gating score: inheritance-aware framework score S (upgrades U1 and U3, proposed before any model existed).*

- **Why.** Selection in gnomAD only sees variants that harm heterozygotes. G2 showed this directly: dominant regions are depleted (RNU4-2 T-loop O/E 0.10), recessive ones are not (RNU4ATAC 0.83–1.09). A selection score alone cannot rank recessive disease variants, and every G3 unit except the pooled RNU6 paralogs is recessive.
- **Definition.** For each G3 unit, by its established inheritance mode (guidance Table 1):
  - dominant (RNU6 paralogs): S = f, the cross-fitted selection score;
  - recessive (RNU4ATAC, RNU6ATAC, RNU12): S = F, the function score. F is the homology transfer h where the position maps to RNU4-2, else the function head d;
  - other or mixed: the mean within-unit percentile rank of f and F.
- **Homology transfer h** (`homology_transfer.py`, label-free). RNU4-2 and RNU4ATAC are aligned with Infernal 1.1.5 to both the U4 (RF00015) and U4atac (RF00618) covariance models, match states only. A RNU4ATAC position is mapped only where both alignments agree. That is 38 positions, RNU4ATAC n.1 and n.26–64 onto RNU4-2 n.1 and n.19–61: the 5′ stem-loop with the k-turn bound by Snu13/Prp31 in both spliceosomes, and stem I. h = −(SGE score of the same alternative base) when the reference bases match, else −(mean SGE score at the RNU4-2 position). 114 of 390 RNU4ATAC SNVs get h.
- **Decision.** The authors chose S as the G3 gating score on 13 October 2026, before any model was fitted to real counts.
- **Reporting.** S gates G3 with the frozen bar (weighted ΔAUROC vs CADD ≥ 0.15, interval excluding 0). The selection-only f result is reported next to it in every table.
- **Disclosure.** This choice was made after CPU dry runs printed, on the real non-blind RNU4ATAC labels, CADD 0.755 (L2 vs L3) and the function head d at 0.768 and 0.827, with f+d at 0.746 and 0.807. Those values came from representations of models fitted to synthetic counts. h itself had never been evaluated on any label. From this amendment on, dry runs permute patient labels within each unit, so no further label-based result can be seen before the frozen run.

*A2.9 Exploratory upgrades (never gating).*

- **U3, function head.** Ridge regression of the RNU4-2 SGE score on [frozen node state, variant features], with λ chosen by leave-one-region-out R².
  - It is fitted only in fits that do not hold out RNU4-2, and never scores RNU4-2.
  - Reported on the G3 units alone, and as the mean percentile rank with f.
  - This uses SGE labels in a separate head, never for f, early stopping or hyperparameters.
- **U1, inheritance (descriptive part).** In RNU4-2 (strict model), f separates dominant from recessive P/LP SNVs. Hypothesis: AUROC > 0.5, with a bootstrap interval.
- **Preliminary calibration.** Section 10 rules, using γ_g·f pooled over the non-blind minor and U6 units. The final Aim 4 comes after the single unsealing.
- **Additional comparators.** GPN-MSA, ReMM, MobiDeep (if released) and an RNA language model are added later as comparators. The frozen bars keep CADD.

*A2.10 Sealing as stage 0 (section 5).*

- **Command.** `labels.py seal` fetches the RNU2-2 and RNU5B-1 labels (ClinVar P/LP and B/LB with ≥ 1 star; RNUdb curated P/LP) straight into `sealed/blind_labels.tsv.enc`. Encryption is AES-256-CBC with PBKDF2 at 200,000 iterations. The passphrase is held by the authors, not the analyst.
- **What it prints.** Only counts and SHA-256s, which are written to `sealed/SEALED_MANIFEST.jsonl` and OpenTimestamps-stamped on the machine that seals.
- **Deviation from section 5.** The hashes are timestamped there instead of in a PREREG edit before training, because the GPU machine cannot push. They are copied into this log when the run's results are committed.
- **Guard.** `train.py` refuses real counts until the manifest records a SEALED event.

**A3, evaluation-only code, comparators and the literature tally (sections 4, 9, 10, 11). Recorded 13 October 2026, after A2 was frozen and before any result of a model fitted to real counts existed.** Nothing in A3 changes the model, its training, the A2 gating rules or any frozen bar. `amendment_A3_manifest.txt` hashes the files and is OpenTimestamps-stamped.

- **Comparators** (`comparators.py`, reported only, never gating; the bars stay against CADD).
  - phyloP for 100 vertebrates and 447 mammals (UCSC tracks phyloP100way, phyloP447way).
  - ReMM v0.4 (hg38).
  - GPN-MSA: minus its precomputed log-likelihood ratio.

  Coverage is 100% (phyloP100 98.9%). `comparators_report.py` reports them on the Aim 2 and per-unit endpoints.
- **Aim 4** (`calibrate.py`, implementing section 10 as written).
  - Scored as follows:
    - per gene: the framework score S (A2.8b), with f where F is unavailable;
    - pools: positives grouped by their own inheritance, dominant (scored by f) and recessive (scored by F), each against the L3 variants of the genes contributing positives.
  - Dominant and recessive positives are also calibrated separately within a gene when each has ≥ 10.
  - Ancestry-specific L3: PASS and AC_g ≥ 9·AN_g/AN at that site (at least 1).
  - Flagged positives are those whose literature source states in-silico or conservation evidence was used. ClinVar and RNUdb positives are unflagged (unknown).
  - The clinician table holds the per-gene Supporting/Moderate/Strong PP3 and BP4 thresholds of S at prior 0.10.
- **Final evaluation** (`final_eval.py`, section 9). It unseals once, logged, and refuses dry runs. RNU2-2 (mixed: mean rank of f and F) and RNU5B-1 (dominant: f) are scored by the fits that held out their clusters, C03 and C06. It reports AUROC and Δ vs CADD for S, every component and every comparator; then Aim 3 with the blind genes as a second estimate; then the final calibration.
- **Release** (`release.py`, section 11). A scores-only Zenodo release right after the run starts the pseudo-prospective clock. The final release adds the thresholds and per-variant evidence labels.
- **Literature tally** (`curate_literature.py`). Three searches by gene group, never touching RNU2-2 or RNU5B-1. Variants count if the reporting study calls them pathogenic, likely pathogenic or disease-causing. HGVS n. positions are mapped to GRCh38, and the reference base must match the transcript; mismatches are listed and dropped. Output goes to `labels/L2_curated_additions.tsv`, which evaluate.py and calibrate.py already merge.
- **Seen before freezing A3.** A comparator-report code check with shuffled SGE scores still showed real ReNU-vs-population AUROCs for the comparators, because the ReNU category had not been shuffled: CADD 0.621 (known from G1), phyloP100 0.546, phyloP447 0.547, ReMM 0.513, GPN-MSA 0.673. That check now shuffles the category too. No other label-based number was computed for the comparators.

**A4, analysis plan for the open endpoints (recorded 14 October 2026, before any endpoint of a real-count model was viewed).** The end-of-run evaluation written by `run_a100.sh` is moved aside unopened (`results_autorun_unopened/`, hashes in `autorun_unopened.sha256`). The **official** v2 evaluation is a single run of `evaluate.py` after the literature tally (`labels/L2_curated_additions.tsv`, dated) is committed. Additional pre-specified analyses, all secondary:
1. *Versus CADD:* (a) dominant genes (RNU4-2 Aim 2 metrics; blind RNU5B-1) model vs CADD; (b) combined score = mean within-gene rank of S and CADD vs CADD alone, per G3 unit and blind gene; (c) PP3/BP4 evidence strength reached by S vs CADD per gene; (d) number of ClinVar VUS SNVs reaching PP3 or BP4 under S vs CADD.
2. *Inheritance prediction:* positions with high F and low f predicted recessive-acting; high f predicted dominant-acting; tested on RNU4-2 (and RNU2-2 at unsealing) dominant vs recessive P/LP SNVs, AUROC of (rank f − rank F).
3. *3D selection map:* per-position cross-fitted f mapped onto cryo-EM spliceosome structures; hypothesis: positions within 4 Å of protein or partner RNA have higher f than other positions (one-sided Mann-Whitney).
4. *Candidate genes:* named genes ranked by the mean of their top-decile f; reported as predictions, timestamped by the Zenodo release.
5. *v3 (exploratory follow-up):* same pipeline with UK Biobank counts added (`v3_counts.py`, `run_v3.sh`), run only after the official v2 evaluation; claims of improvement require variants unseen by both versions.

**A7, zero-shot test of RNA-family profiles on non-snRNA disease RNA genes (recorded 14 October 2026, before any A7 score or label was computed).** The code is the specification (`a7_ncrna.py`).

- **Hypothesis.** At structured non-coding RNA genes, the Rfam family profile ranks patient variants above population variants better than genome-wide scores (CADD; phyloP447 as genome-wide conservation).
- **Genes and units.** RMRP (RF00030, recessive); TERC (RF00024, dominant); SNORD118 (U8, RF00096) and RNU7-1 (RF00066), pooled as one unit because each has fewer than 10 P/LP SNVs. GRCh38 coordinates from Ensembl (in the script). No model, threshold or analysis in this project has used these genes.
- **Primary score.** `rfam_alt_rarity` = 1 − frequency of the alternative base in the variant's Rfam consensus column (cmemit profile, 20,000 samples, as `model_features.py`); positions in insert states score 1, as in v4. Parameter-free. Chosen because it was the best single feature in the post hoc v4 ablation on blind RNU2-2 (disclosed).
- **Secondary.** The A5-frozen v4 coefficients applied to the same feature definitions; Evo 2 and lab function take the frozen missing value; partner-binding features are 0.
- **Labels.** Positives: ClinVar P/LP SNVs with ≥ 1 review star. Negatives: gnomAD v4.1 genome PASS AC ≥ 9, minus positives.
- **Two stages.** `a7_ncrna.py score` computes every feature and score without reading any label and writes `a7/a7_manifest.txt`. The manifest and scores are committed, pushed and OpenTimestamps-stamped. Only then does `a7_ncrna.py eval` run; it refuses unless `a7/` is committed and unchanged.
- **Decision (pre-registered).** A7 is met if (i) the weighted (by positives) mean per-unit ΔAUROC of the primary score vs CADD and vs phyloP447 are both > 0, with 95% class-stratified bootstrap intervals (2,000 replicates, seed 20261015) excluding 0 and Holm-adjusted p < 0.05; and (ii) the primary AUROC is ≥ 0.70 in at least 2 of the 3 units.
- **Sensitivity (reported, not gating).** Consensus-column positions only; negatives at AC ≥ 20; positives without PP3/BP4 in their ClinVar submissions (ACMG codes read as in `circularity_check.py`); GPN-MSA and ReMM as further comparators; v4 frozen vs CADD.
- **Either outcome is reported.** If not met, the family-profile result stays specific to the spliceosomal snRNAs.

**A8, a second blind test of the frozen v4 on literature-only variants of the blind genes (recorded 14 October 2026, before any literature variant of RNU2-2 or RNU5B-1 was curated).** The code is the specification (`a8_literature.py`).

- **Why.** The first blind test rests on 22 RNU2-2 SNVs. RNU2-2 and RNU5B-1 were never in any v4 training set and the v4 scores are hash-frozen (A5), so patient variants reported only in the 2025–2026 literature are a second, independent blind set.
- **Curation.** A tally of every RNU2-2 and RNU5B-1 variant reported as pathogenic, likely pathogenic or disease-causing in peer-reviewed or preprint papers, supplements, OMIM or LOVD (same columns and wording rule as the A3 literature tally; `curate_literature.py` mapping with reference-base check). Curators never see any score.
- **Positives (primary).** Literature SNVs in RNU2-2 not in the sealed blind P/LP set and not in any training label file. **Negatives.** L3 (gnomAD PASS AC ≥ 9) minus every positive of either set.
- **Order.** Curate → `a8_literature.py map` → commit and push `labels/A8_literature.tsv` → `a8_literature.py eval` (refuses unless the labels are committed and unchanged).
- **Decision (pre-registered; needs ≥ 8 new RNU2-2 positives, else descriptive only).** Met if (i) the frozen v4 AUROC exceeds its rarity-only null (1,000 draws of as many random rare RNU2-2 SNVs) with one-sided p < 0.05, and (ii) v4 − CADD AUROC > 0 with the 95% class-stratified bootstrap interval excluding 0.
- **Secondary.** The same for RNU5B-1 (≥ 5 new positives); sealed + new pooled; Rfam alt rarity, Rfam ref column freq, phyloP447 and GPN-MSA reported with their rarity nulls; A6 run unchanged (frozen score and threshold) on every blind-gene P/LP SNV whose inheritance the literature states.
- **Seen beforehand.** The 22 sealed RNU2-2 positives and the ACMG codes of their ClinVar submissions; titles and abstracts of the 2025–2026 RNU2-2 papers (read during the novelty audit, including that Greene 2026 lists 32 variants it calls destabilising); no literature variant list, and no per-variant v4 score of any RNU2-2 SNV.

**A9, rarity-only inflation of non-coding predictor benchmarks (recorded 14 October 2026, before any score of an A9 variant was fetched).** The code is the specification (`a9_inflation.py`).

- **Why.** On RNU2-2, any score that tracks rarity separates random rare SNVs from gnomAD-common "benign" variants (null AUROC 0.53–0.80 depending on the score). If this holds genome-wide, published non-coding benchmarks overstate predictors by a score-specific amount.
- **Positives.** ClinVar GRCh38 VCF (fileDate 2026-10-11) SNVs classified Pathogenic, Likely pathogenic or Pathogenic/Likely pathogenic with ≥ 1 review star, whose molecular consequences are all non-coding (any coding or canonical splice-site term excludes). Classes by first match: 5′UTR, 3′UTR, ncRNA exon, intronic, upstream; at most 1,500 per class (seed 20261017). Selected and committed before any score was fetched: intronic 1,500 (of 3,438), 5′UTR 148, ncRNA exon 85, 3′UTR 51, upstream 3 (`a9/positives.tsv`, `a9/select_info.json`).
- **Negatives.** gnomAD v4.1 genome PASS SNVs with AC ≥ 9 within ±100 bp of a positive, not ClinVar P/LP.
- **Rarity-only null.** For each positive, one SNV from its own ±100 bp window that is neither a negative nor a positive; 200 draws; null AUROC against the same negatives.
- **Scores.** CADD v1.7, GPN-MSA, ReMM v0.4; phyloP447 if the bigWig is readable.
- **Metrics.** Per class with ≥ 50 positives and pooled: raw AUROC (1,000 bootstraps), null mean (2.5–97.5 percentiles), inflation = null − 0.5, corrected = raw − null.
- **H1 (primary).** Pooled inflation ≥ 0.05 for at least 2 of CADD, GPN-MSA, ReMM. **H2 (descriptive).** Whether the ranking of tools by raw AUROC differs from the ranking by corrected AUROC in a class with ≥ 200 positives; a reordering within overlapping intervals is reported as such, not as a finding.
- **Code check before freezing.** The evaluation was run on synthetic window data with random scores: every inflation was within ±0.03 of 0 and H1 was not met, as expected under the null.
- **Either outcome is reported.** If H1 is not met, the RNU2-2 inflation is reported as locus-specific (consistent with snRNA hypermutability).

**A10, inflation in missense and SpliceAI benchmarks, and two fixes (recorded 14 October 2026, after A9 and before any A10 score was fetched).** The code is the specification (`a10_extend.py`).

- **Missense.** ClinVar (fileDate 2026-10-11) P/LP SNVs with ≥ 1 star whose consequences include missense_variant and no splice-site, nonsense, stop or start term: 59,582 available, 1,500 kept (seed 20261019; `a10/missense_positives.tsv`, committed before any fetch). ±100 bp windows; the analysis is restricted to SNVs with an AlphaMissense score, so positives, gnomAD-common negatives (PASS AC ≥ 9) and the rarity pool are all missense. Tools: AlphaMissense, CADD, GPN-MSA, phyloP447.
- **SpliceAI.** The A9 positives and common negatives plus, per positive, up to 10 pool SNVs from its own window (seed 20261019), scored with SpliceAI 1.3 (GRCh38 annotation, distance 50; score = maximum delta; unscored = 0). The null is drawn from that sample.
- **Fixes, for every tool.** (a) Rarity-matched AUROC: positives against all pool SNVs in their windows (rare or absent variants, no frequency contrast). (b) AUROC against ClinVar Benign/Likely benign SNVs (≥ 1 star) in the same windows and class (non-coding: 9,845; missense: 1,126). Caveat stated in advance: ClinVar benign classifications often rest on frequency (BA1/BS1), so (b) may retain part of the rarity effect; (a) is the cleaner fix.
- **H1.** AlphaMissense pooled inflation ≥ 0.05. **H2.** SpliceAI inflation ≥ 0.05 on intronic positives. **H3 (descriptive).** The share of tool × unit rows in which the rarity-matched AUROC lies closer than the raw AUROC to the null-corrected value + 0.5.
- **Code check.** The evaluation was run on synthetic windows with random scores before freezing: inflation within ±0.02 and AUROCs near 0.5.

**A11, re-analysis of a published benchmark (recorded 14 October 2026, before any A11 score was fetched).** The code is the specification (`a11_ncboost.py`).

- **Benchmark.** NCBoost (Caron, Luo & Rausell, Genome Biology 2019, doi 10.1186/s13059-019-1634-2): 655 public curated pathogenic non-coding SNVs (Table S1) against 6,550 region-matched common dbSNP SNVs (MAF > 5%, no clinical assertion; Table S3). The paper ranked ReMM above Eigen, DeepSEA and CADD. Tables (CC BY 4.0) are copied to `a11/` as TSV (GRCh37) and lifted to GRCh38 with pyliftover; variants that fail liftover or whose GRCh38 reference disagrees are dropped and counted. NCBoost's own score is not openly downloadable and is not re-scored.
- **Metrics.** As A9, against the published negatives: raw AUROC, window-matched rarity-only null (one rare/absent SNV from each positive's ±100 bp window, 200 draws), inflation, corrected AUROC, rarity-matched AUROC. Tools: CADD v1.7, ReMM v0.4, GPN-MSA, phyloP447.
- **H1.** Inflation ≥ 0.05 for at least 2 of the 4 tools. **H2 (descriptive).** Whether the published order ReMM > CADD holds for raw, corrected and rarity-matched AUROC.
- **Code check.** Synthetic data with random scores: every inflation within ±0.002 and H1 not met.

**A12, prospective blind tests on new patient variants, including never-labelled genes (recorded 15 October 2026, before any A12 label was built or scored).** The code is the specification (`a12_new_variants.py`).

- **Why.** The blind evidence rests on one gene (RNU2-2). v4 and its scores for all 17,385 SNVs were hash-frozen on 14 October 2026 (A5), so any P/LP variant not in a label file committed by then is unseen by v4, and several loci never had any label (RNU5A-1, RNU6-2, RNU6-7, RNU11, RNU4-1, RNU5D/E/F-1, RNVU1-*). They are blind genes by construction.
- **Sources.** ClinVar GRCh38 VCF, re-scanned monthly (`scan`; each scan's fileDate and per-gene counts are committed in `a12/scan_log.tsv`); literature tallies under `lit/A12_TALLY_PROMPT.md`, mapped with the frozen wording rule of `curate_literature.py` (`map`).
- **Positives.** ClinVar P/LP with ≥ 1 review star, or literature P/LP / disease-causing, SNVs, not in L2_nonblind, L2_curated_additions, A8_literature or the sealed blind set. A literature positive that ClinVar classifies B/LB (≥ 1 star) is dropped and listed. RNU2-2 and RNU5B-1 positives are built only on the machine holding the sealed labels (`--blind`).
- **Units.** U1 RNU6 adRP (RNU6-1, -2, -7, -8; primary; ≥ 8 positives). U2 U5 (RNU5A-1, RNU5B-1; primary; ≥ 5). U3 RNU2-2 beyond the sealed and A8 sets (primary; ≥ 8). U4 every other locus without training labels (secondary; ≥ 5). U5 RNU4-2, RNU4ATAC, RNU12, RNU6ATAC, new variants only (secondary; these genes were seen in training). RNU6-1 and RNU6-8 contributed three training labels; their genes are therefore not fully unseen, which is reported.
- **Negatives and null.** L3 (gnomAD v4.1 PASS AC ≥ 9) in the unit's loci minus every positive ever used; rarity-only null from 1,000 draws of as many SNVs from the unit's loci that are neither L3 nor any positive.
- **Scores.** Frozen v4 (primary; hash-checked), Rfam reference column frequency (key secondary: the Rfam feature that held up in A7), Rfam alternative-base rarity, CADD, phyloP447, GPN-MSA.
- **Decision per primary unit.** MET if (i) the v4 excess over its rarity null has one-sided p < 0.05/3 (Bonferroni over the three primary units, fixed now so that units evaluated at different times share one alpha) and (ii) v4 − CADD AUROC has a 95% class-stratified bootstrap interval above 0. Below the minimum: descriptive only.
- **One look per unit.** A unit is evaluated once, after `labels/A12_labels.tsv` is committed; the look is logged in `a12/eval_log.tsv`. A repeat needs `--again` and is exploratory.
- **Reclassification arm (descriptive).** ClinVar VUS committed on 14 October 2026 (`labels/clinvar_vus_snrna.tsv`) that later become P/LP or B/LB: v4 and CADD AUROC once each class has ≥ 5.
- **Code check.** Synthetic tests: random scores give no excess; a planted v4 signal gives MET; prior, sealed and ClinVar-benign rows are excluded.

**A13, the A11 re-analysis applied to other published benchmarks, first TraitGym (recorded 15 October 2026, before any A13 score was fetched).** The code is the specification (`a13_benchmarks.py`).

- **Benchmark.** TraitGym Mendelian (Benegas et al. 2025; `songlab/TraitGym`, `mendelian_traits_matched_9/test.parquet`, GRCh38): 338 causal non-coding variants for Mendelian traits, each with 9 matched controls (3,042). Frozen into `a13/traitgym_mendelian/benchmark.tsv` (md5 recorded) and committed before fetching. Any further benchmark is added by the same `freeze` → commit → `fetch` → `eval` steps, each with its own dated deviation-log entry before `fetch`.
- **Applicability check first.** The share of negatives that are gnomAD-common (PASS AC ≥ 9) is reported. H1 is evaluated only if ≥ 50% are; otherwise H1 is "not applicable" and the null is descriptive.
- **Metrics.** As A11: raw AUROC, window-matched rarity-only null (200 draws), inflation, corrected, rarity-matched AUROC; plus the within-match-group AUROC (pairs only inside a match group).
- **H1.** Inflation ≥ 0.05 for at least 2 of CADD, ReMM, GPN-MSA, phyloP447. **H2 (descriptive).** Whether the tool order changes between raw, corrected and rarity-matched AUROC (overlapping intervals are ties).
- **Summary.** `raritybench.py` tabulates A9, A10, A11 and A13 (share of above-chance AUROC explained by rarity); exploratory.
- **Code check.** Synthetic windows: a random score shows no inflation; a score that tracks rarity shows inflation and ~0 corrected AUROC.

**A14, ClinGen PP3/BP4 recalibration with rare benign variants (recorded 15 October 2026, before any A14 variant was selected or scored).** The code is the specification (`a14_recalibration.py`).

- **Why.** ClinVar B/LB negatives did not remove the rarity-only inflation (A10). ClinGen calibrations use such B/LB sets, but PP3 is applied to variants that are already rare (PM2), so the evidence a score carries is its separation of pathogenic from rare benign variants.
- **Variants.** ClinVar missense SNVs (A10 consequence rule), P/LP or B/LB with ≥ 1 star; at most 30 per gene per class and 6,000 per class (seed 20261023); committed before fetching.
- **Sets.** ALL: P/LP vs all B/LB (the ClinGen design). RARE: P/LP vs B/LB, both absent or AC ≤ 5 in gnomAD v4.1 genomes (PASS).
- **Method.** `calibrate.py` local-LR procedure unchanged; prior 0.0441 primary (Pejaver 2022), 0.10 secondary. Tools: AlphaMissense (primary), CADD v1.7, GPN-MSA, REVEL if a local REVEL v1.3 file is supplied.
- **Validity check first.** Our ALL-set AlphaMissense thresholds against the published ones (0.792 / 0.906 / 0.972).
- **H1 (primary).** At prior 0.0441, AlphaMissense's published PP3 Strong threshold (≥ 0.972) is not met in the RARE set (the 5th-percentile local LR falls below the Strong LR somewhere at or above 0.972). **H2 (descriptive).** Per tool and strength, ALL vs RARE thresholds (stricter, unchanged, unreachable), and the share of rare P/LP whose PP3 strength drops under RARE thresholds.
- **Caveat stated in advance.** Rare B/LB were classified without frequency evidence and may rest partly on BP4, which makes them easier to separate; a shift towards stricter thresholds is therefore conservative.

**A15, validity of the rarity-only null and intervals for every corrected AUROC (recorded 15 October 2026, before the simulation or any interval was computed).** The code is the specification (`a15_null_validity.py`).

- **Estimand.** T = AUROC of pathogenic vs rare benign variants (clinical use: the variant is already rare). The rarity-matched AUROC M estimates T diluted by the damaging fraction π of the rare pool: T ≈ (M − π/2)/(1 − π).
- **Simulation.** Windows with 1 pathogenic positive, 60 rare SNVs (fraction π damaging) and 6 common benign SNVs; score = β·damaging + γ·rare + N(0, 1); β ∈ {0, 0.5, 1, 2}, γ ∈ {0, 0.5, 1, 2}, π ∈ {0.05, 0.2}; 1,000 windows, 100 replicates.
- **Expectations.** E1: whenever γ > 0, raw AUROC overstates T and |M − T| < |raw − T|. E2: when β = 0, M and corrected + 0.5 are near 0.5 (within 0.03) while raw is not. E3: deattenuated M is within 0.02 of T on average. Each is reported as holding or failing.
- **Intervals.** Window (cluster) bootstrap over positives, 200 replicates, 50 null draws each: raw AUROC, null mean, corrected and rarity-matched AUROC per tool, and paired differences between tools for A9, A10, A11 and A13. A reordering of tools is called a reordering only if its paired interval excludes 0.

**A16, PP3 evidence against the rare-variant population (recorded 15 October 2026, after the A14 results and before any A16 window was sampled or scored).** The code is the specification (`a16_population_calibration.py`).

- **Why.** Clinical PP3 is applied to a rare variant drawn from the population of rare variants, not from ClinVar's benign set, and A14 showed curated rare B/LB are much easier to separate than random rare variants.
- **Sets.** P: rare ClinVar P/LP missense from A14 (gnomAD genome AC ≤ 5). U: every missense SNV observed in gnomAD v4.1 genomes (PASS, AC 1–5) in 3,000 random 1-kb windows centred on uniformly sampled AlphaMissense positions (seed 20261025), minus ClinVar-labelled variants (≥ 1 star). Curated B/LB (A14 ALL and RARE) are reported side by side.
- **Positive-unlabelled correction.** U = (1 − π)B + πP; true LR = (1 − π)/(1/LR_obs − π). π grid {0, 0.01, 0.02, 0.05}; **primary π = 0.02**; π = 0 (U benign, most conservative) and 0.05 (stress test) alongside. π_max = min over scores of 1/LR_obs, with a bootstrap interval; verdicts near π_max are not interpreted.
- **Band LRs (primary readout).** For each published band (AlphaMissense 0.792 / 0.906 / 0.972; REVEL 0.644 / 0.773 / 0.932; CADD 25.3 / 28.1), P(band | P)/P(band | comparison), 0.5 pseudocounts, 1,000 bootstraps, π-corrected. Verdict at prior 0.0441: DEMONSTRABLY MET (5th percentile ≥ ClinGen LR), DEMONSTRABLY NOT MET (95th percentile < it), otherwise UNDETERMINED.
- **Break-even π.** π* = (need/LR_obs − 1)/(need − 1), from the point and the 95th-percentile LR; if π* > π_max the band cannot reach its strength for any π the data allow.
- **H1 (primary).** AlphaMissense's published Strong band (≥ 0.972) is DEMONSTRABLY NOT MET against U at π = 0.02. **H2 (descriptive).** All verdicts by tool, band, π and comparison set; the population false-PP3 burden (share of U in each band, π-free); local-LR thresholds against U.
- **Status.** Designed after A14 was seen; reported as a pre-registered follow-up, not as an original aim.

---

# END OF MATERIALS
Now complete tasks 1–8 from the instructions at the top.
