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
- This is a **revised** version. An earlier version was reviewed once; section 6b of the verification report lists what that review raised and what changed. Judge the current version on its merits, and say whether the changes resolve the earlier points.
- The study has a single author with no wet-lab work. Every input is public data (ClinVar, gnomAD, UK Biobank summaries, published SGE and assay data, the literature).
- The pre-registered primary model (v2) failed its endpoints, and the paper says so. The model that succeeded (v4) was designed after that failure, then frozen and timestamped before the sealed labels were opened.
- Some items are pending because they run on a separate GPU server; they are listed in the verification report, section 6. Judge the paper as it stands, and say how those pending results could change your view.

---

# MATERIALS UNDER REVIEW
Five documents follow: manuscript, supplement, timeline, verification report and references. The full pre-registration and deviation log are omitted for length; ask for them if you need a specific entry.

---

# DOCUMENT 1: MANUSCRIPT

## Whole-genome conservation fails at key spliceosomal snRNA genes, and RNA-family alignments recover the signal: a pre-registered evaluation with sealed patient labels

**Running title:** RNA-family conservation for snRNA variants

**Author:** [Author name]¹ — ¹[Affiliation]. Correspondence: [email]

*Manuscript, version of 5 November 2026. Numbers come from committed result files in `snrna/results/` or, where marked ‡, from outputs of the analysis server transcribed into the dated deviation log (`snrna/PREREG.md` §14) at the time they were produced. Citation numbers are generated by `paper/build_refs.py`.*

---

### Abstract

**Purpose.** Variants in spliceosomal small nuclear RNA (snRNA) genes are a frequent cause of neurodevelopmental disorders. Current guidance advises against computational evidence (ACMG/AMP PP3/BP4) for these genes until a tool has been calibrated on them, because genome-wide predictors perform poorly there. We asked where and why genome-wide scores fail, and whether RNA-family (Rfam) alignment profiles do better.

**Methods.** We pre-registered and timestamped an evaluation of variant-effect scores at 55 snRNA loci (17,385 possible single-nucleotide variants). Patient labels for RNU2-2 and RNU5B-1 were sealed before any model was fitted. Because negatives in small genes are population-observed variants, scores were compared against rarity-only nulls matched on gene, disease region and RNA structure.

**Results.** The pre-registered model, built on population allele counts, failed its endpoints. A model designed afterwards and frozen before unsealing reached an AUROC of 0.892 on 22 sealed RNU2-2 patient variants (CADD 0.501). Its best single inputs were two parameter-free Rfam features, identified post hoc. Across 35 public RNU2-2 patient variants, the disease region explained most of the Rfam signal beyond rarity. A residual excess survived matching on region and base pairing (+0.07, p = 0.009), which no genome-wide score showed. In RNU4-2 saturation genome editing, with no population data involved, the Rfam reference-column frequency classified function-depleted variants with an AUROC of 0.80 (95% CI 0.72–0.88, resampling positions), against 0.49–0.60 for genome-wide scores. The advantage concerned which positions matter, not which base substitution. At RNU4-2, the 447-way whole-genome alignment contains only primate sequences (186 species and no non-primate mammal, against 442 species at its paralog RNU4-1), which caps phyloP447 at 0.90. Across five disease genes, Rfam matched or exceeded CADD, with the gains where CADD fails (RNU2-2, RNU4-2, RNU12). The advantage did not extend to other RNA families (pre-registered negative).

**Conclusion.** At several spliceosomal snRNA genes, whole-genome alignments carry little usable conservation, and curated RNA-family alignments recover position-level signal. snRNA-specific computational evidence should be built on family alignments. Clinical calibration remains exploratory.

**Keywords:** snRNA; RNU2-2; RNU4-2; variant interpretation; conservation; Rfam; PP3/BP4; pre-registration

---

### Introduction

De novo variants in RNU4-2, which encodes the U4 small nuclear RNA (snRNA) of the major spliceosome, cause ReNU syndrome. It is one of the most frequent monogenic neurodevelopmental disorders (NDDs) found by genome sequencing.<sup>1,2</sup> Other snRNA disorders followed:
- dominant and recessive disorders of RNU2-2;<sup>3–6</sup>
- dominant disorders of RNU5B-1;<sup>4,7</sup>
- a recessive RNU4-2 syndrome.<sup>8,9</sup>

The minor-spliceosome snRNA disorders were known earlier: RNU4ATAC,<sup>10,11</sup> RNU12<sup>12,13</sup> and RNU6ATAC.<sup>14,15</sup> Dominant variants in RNU4-2 and RNU6 genes also cause retinitis pigmentosa.<sup>16</sup> Diagnostic genome sequencing now finds snRNA variants routinely.<sup>17</sup>

Interpretation lacks computational evidence. The proposed guidance for snRNA variant classification states that no tool yet captures features relevant to snRNA variant deleteriousness and that genome-wide tools perform poorly.<sup>18</sup> It recommends in silico evidence (PP3/BP4) only after a tool has been calibrated for snRNAs. Under thresholds calibrated on missense variants,<sup>19</sup> CADD scores the first two pathogenic RNU2-2 variants identified in the benign range.<sup>18</sup> More generally, genome-wide thresholds calibrated on missense variants give false benign evidence for pathogenic non-coding variants.<sup>20</sup> Saturation genome editing (SGE) of RNU4-2 shows that measured function outperforms in silico prediction in that gene,<sup>9</sup> but no comparable data exist for the other genes.

Why genome-wide scores fail here is unclear. snRNA genes occur in multiple paralogous copies among many pseudogenes,<sup>18,21</sup> and the spliceosomal snRNA genes are germline mutation hotspots.<sup>4,22</sup> Either could disturb the whole-genome alignments and population data on which genome-wide scores rely. RNA-family covariance models (Rfam) describe each snRNA family's sequence and secondary structure from curated alignments, independently of whole-genome alignment.<sup>23,24</sup> They are a natural alternative source of conservation.

Evaluation is a second obstacle. In small genes, negatives are usually population-observed variants. Rarity correlates with conservation, so such negatives flatter any conservation-like score.<sup>25–29</sup> Patient variants are also ascertained in disease regions, so any score that encodes position benefits. We therefore compared scores against rarity-only nulls matched on region and structure, and we tested the main claim on functional data that involve no population labels.

The study was designed to be falsifiable. The plan, definitions and decision rules were frozen and timestamped before analysis. Patient labels for two genes were sealed before any model was fitted. Every later change is logged together with what had been seen when it was made. We report the failed endpoints alongside the successful ones.

### Materials and Methods

#### Study design and record of order
The pre-registration (`PREREG.md`), the analysis code and amendments A1–A6 were hashed and timestamped with OpenTimestamps; the proofs are anchored in Bitcoin blocks 971237–971363 (13–14 October 2026; Supplementary Table S1). Later amendments (A7–A16) and all post hoc analyses rest on git commit times pushed to a public repository. Each post hoc analysis labelled "specified before running" had its code committed before its result. Changes made after data were seen are recorded in a dated deviation log (`PREREG.md` §14), with what had been viewed at the time.

The sealed labels were fetched into an AES-256-encrypted file before any model was fitted. They were decrypted once, and every access was logged in `sealed/SEALED_MANIFEST.jsonl`. That log is held on the analysis server and will be released with an OpenTimestamps proof. Until then, the sealing claim rests on the procedure and the anchored freeze of the model (block 971362).

#### Loci, variants and scores
We included the GENCODE v49 spliceosomal snRNA genes and paralogs that passed coverage and paralog quality control. This gave 55 loci, 20 of them named genes, with all 17,385 possible SNVs. Comparators were:
- CADD v1.7;<sup>30</sup>
- phyloP from the 447-way alignment;<sup>31,32</sup>
- ReMM v0.4;<sup>33</sup>
- GPN-MSA;<sup>34</sup>
- Evo 2 (7B) log-likelihood ratios;<sup>35</sup>
- ViennaRNA 2 folding and duplex free-energy changes for intermolecular snRNA pairs.<sup>36</sup>

Each locus was aligned to its Rfam 15.1 covariance model with Infernal 1.1.5 (cmalign).<sup>23,24</sup> Two features came from the alignment:
- **reference-column frequency:** the frequency of the reference base in its alignment column; this is constant across the three alternatives at a position;
- **alternative-base rarity:** one minus the frequency of the alternative base.

Neither feature uses a label or has a fitted parameter.

#### Labels
**Positives** were pathogenic or likely pathogenic (P/LP) SNVs from three sources:
- ClinVar, with at least one review star;<sup>37</sup>
- RNUdb curation;<sup>18</sup>
- a literature tally under a frozen wording rule: a source must call the variant pathogenic, likely pathogenic or disease-causing, and any source calling it a VUS excludes it.

Curation was AI-assisted and released for audit.

**Negatives (L3)** were gnomAD v4.1 genome PASS SNVs with allele count ≥ 9,<sup>38</sup> excluding positives. L3 holds 1,233 SNVs. They are population-observed but mostly rare: median allele frequency 1.8 × 10⁻⁴, with 5.4% at ≥ 1%. Population-observed is not benign for recessive genes, where L3 can contain carrier alleles.

#### Pre-registered gates and models
**G1** asked whether a foundation model alone (Evo 2) already predicts RNU4-2 function. **G2** asked whether population data show selection at known disease regions. **v2**, the primary model, was a graph-transformer selection model trained on gnomAD allele counts under the Roulette mutation-rate model,<sup>39</sup> with paralog clusters as the hold-out unit. Its framework score was inheritance-aware:
- selection (f) for dominant genes;
- lab-derived function (F) for recessive genes;
- the mean rank of f and F for mixed-mode genes such as RNU2-2.

A version with UK Biobank counts (v3) was exploratory.<sup>40</sup>

**v4 (amendment A5)** was designed after the v2 endpoints were seen. It is an L2-regularised logistic regression (C = 0.1, balanced class weights, no tuning) on within-gene percentile ranks of 13 features:
- the six comparators;
- the two Rfam features;
- reference pairing and intermolecular-pairing loss;
- absolute folding ΔΔG;
- local constraint (observed versus expected variation in a 21-nt window);
- homology-transferred RNU4-2 SGE function (RNU4ATAC only).

It was trained on non-blind genes only. Its code and output hashes were anchored in Bitcoin block 971362 (14 October 2026, 05:58 UTC), before the blind labels were decrypted.

#### Blind tests
- **Test 1, sealed labels:** RNU2-2, 22 P/LP vs 127 L3. RNU5B-1 had 1 P/LP and could not be evaluated.
- **Test 2, literature (amendment A8):** RNU2-2 variants reported as pathogenic in the 2025–2026 literature, curated after unsealing under the frozen rule, with the curator instructed not to look up any score. The rule accepted 35 SNVs; 19 were already in the sealed set and 16 were new. Pre-registered conditions: the v4 excess over its rarity null with p < 0.05, and a v4 − CADD AUROC interval excluding 0.
- **Test 3, prospective (A12):** new P/LP variants reported after the freeze, evaluated once per unit when the unit reaches its pre-registered minimum.

#### Rarity-only nulls matched on region and structure
Because L3 is defined by population observation, any score that tracks rarity separates random unlabelled variants from L3. For *n* positives, a null draws *n* SNVs that are neither positives nor L3, scores them against the same L3, and repeats the draw. Three nulls were used:
- **N0, gene-wide:** draws from the whole gene (the original definition);
- **N1, region × pairing:** each draw is matched to its positive's disease-region status (guidance Table 2<sup>18</sup>) and reference base-pairing status;
- **N2, local window:** each draw comes from within ±5 nt of its positive.

All three used 10,000 draws. We report the excess (AUROC − null mean) with a one-sided p = (k + 1)/(B + 1). The same nulls were applied to four position-only baselines: region membership, distance to the region, 5′ proximity, and position unobserved in gnomAD. A score that only encodes position should lose its excess under N1 or N2. Because these analyses need only public labels, they used the 35 public RNU2-2 positives. They were specified after an external review and before running.

#### Functional readouts (no population labels)
- **RNU4-2 SGE:**<sup>9</sup>
  - Spearman correlation with functional damage across the 435 in-transcript SNVs;
  - classification of "depleted" variants (function score < −0.302, the study's threshold), overall and outside the guidance-defined critical region;
  - 95% intervals from resampling positions (145 clusters), because reference-column frequency is constant within a position;
  - a within-position test: centring each score and the function score within position, to ask whether a score ranks the alternatives at the same position.
- **RNU4ATAC splicing assay (exploratory):** published splicing efficiencies of 23 SNVs.<sup>41</sup>

#### Gene-level analyses
- **Per gene** (post hoc, specified before running): Rfam, CADD, phyloP447 and GPN-MSA AUROC in every gene with ≥ 5 public P/LP SNVs (RNU2-2 literature set, RNU4-2, RNU4ATAC, RNU12, RNU6ATAC).
- **Alignment depth** (exploratory, specified before running): the phyloP447 ceiling (95th percentile), the fraction of positions sharing one value, and per-position agreement with the Rfam profile. We predicted that the Rfam advantage over CADD would be larger where phyloP447 agrees less with the profile.
- **Direct alignment check** (post hoc, specified before running): the UCSC Cactus 447-way alignment (bigMaf), read by HTTP range requests, gave the number of species aligned at each human position of seven genes. We predicted that fewer species would be aligned at RNU4-2 than at RNU4-1 and the other genes. A descriptive extension, written after the counts were seen, split aligned species into primates and non-primates using UCSC's species tree.
- **Cross-family test (A7, pre-registered):** zero-shot scoring of RMRP, TERC, SNORD118 and RNU7-1, with scores frozen before labels were fetched.

#### Calibration (exploratory for v4)
We used the ClinGen local likelihood-ratio method as pre-registered:<sup>19,42</sup>
- nearest 20% of labelled variants, at least 10 per class;
- 1,000 bootstraps;
- PP3 thresholds at the 5th percentile;
- prior 0.10, with sensitivity priors.

The pre-registered calibration (Aim 4) used the v2 framework score. The v4 calibration was post hoc and used out-of-sample scores. Negatives were L3. Because L3 can contain carrier alleles of recessive genes, we report PP3 only and no BP4. Our implementation did not reproduce published AlphaMissense thresholds in a validity check (A14): it gave 0.627/0.841/0.985 against the published 0.792/0.906/0.972. Thresholds from it are therefore approximate.

#### Population-negative benchmarks
The rarity-only null was also applied to ClinVar non-coding (A9) and missense (A10) variants, to two published benchmarks (A11 NCBoost,<sup>43</sup> A13 TraitGym<sup>44</sup>), and to a ground-truth simulation (A15). These are reported in Supplement S6.

#### Statistics and software
Analyses used Python 3 (NumPy, pandas, SciPy, scikit-learn). P-values are one-sided where a direction was specified. Sealed-set rarity tests were added after unsealing across 12 scores; we give Bonferroni-adjusted values. Selection of the Rfam features from the ablation is not covered by that correction.

### Results

#### Pre-registered gates and the primary model
- **G1 passed (continue).** Evo 2 alone did not predict RNU4-2 function (ρ 0.303; the bar for "foundation model sufficient" was 0.60). A supervised position-only floor reached ρ 0.570, fitted on neighbouring SGE scores with leave-one-position-out. The frozen rule did not cover a floor that far above every predictor; we report it here as the ceiling that knowledge of position alone attains.
- **G2 was partial.** Population depletion recovered the RNU4-2 critical region (O/E 0.13, p 9 × 10⁻⁹). RNU2-2 n.1–67 was recovered only under the amended rule A1 (O/E 0.62; circular-shift p 0.031). RNU5B-1 n.37–44 was not recovered (circular-shift p 0.051).
- **The primary model failed.** v2 reached ρ 0.144 against the Aim 2 bar of 0.346, and +0.083 against the held-out-gene bar of 0.15 (ΔAUROC vs CADD). Adding UK Biobank counts (v3) did not help (ρ 0.142).
- **The pre-specified RNU2-2 score.** On the sealed RNU2-2 labels, the pre-specified mixed-mode score S reached 0.784‡ (CADD 0.501). Its function component F reached 0.851 and its selection component f 0.659.

#### A frozen post-failure model on sealed RNU2-2 labels
v4 exceeded CADD by +0.186 [0.108, 0.267] in leave-one-cluster-out validation on non-blind genes. On the sealed RNU2-2 labels, v4 reached 0.892 [0.817, 0.951]‡, against:
- CADD 0.501;
- ReMM 0.501;
- phyloP447 0.724;
- GPN-MSA 0.724.

A post hoc ablation found that the best single inputs were the two Rfam features (alternative-base rarity 0.860, reference-column frequency 0.847), and v4 did not improve on alternative-base rarity alone (+0.031, interval including 0). The result is therefore an Rfam result, and the Rfam features were identified after unsealing.

The literature test (A8) met both pre-registered conditions:
- v4 reached 0.910 [0.825, 0.979] on 16 new positives against 125 L3;
- its excess over the rarity null was +0.133 (p = 0.012);
- Δ vs CADD was +0.298 [0.153, 0.426].

The two tests are not independent. Both are RNU2-2, they share 125 of 127 negatives and draw on overlapping literature, and 14 of the 16 new positives are recessive. Under broader label rules (20 and 27 positives), v4 still exceeded CADD (+0.267, +0.237), but its rarity-null excess fell to a trend (p ≈ 0.07).

#### Most of the signal beyond rarity is region; a residual survives matching
On the 35 public RNU2-2 positives (Table 1), the gene-wide null (N0) reproduced the earlier pattern: Rfam excess +0.18 (p = 0.0001), CADD +0.02, GPN-MSA −0.02. Region membership alone, however, also reached an N0 excess of +0.16, and so did 5′ proximity. So the N0 excess does not distinguish pathogenicity from where patient variants are ascertained.

Matching the null on region and base pairing (N1) removed the baselines' excess, by construction. The Rfam features kept +0.071 (p = 0.009) and +0.073 (p = 0.018). No genome-wide score kept any: CADD +0.013, phyloP447 +0.001, GPN-MSA −0.050. Under the strictest null (N2, ±5 nt), the Rfam excess was +0.041 (p = 0.069) and +0.044 (p = 0.088), a trend. N2 removes most positional information, which is what a conservation score supplies. Inside n.1–67 (25 vs 14 variants), Rfam separated poorly (0.59 and 0.66). Outside it (10 vs 111), Rfam and phyloP447 were similar (0.81 each).

**Table 1.** RNU2-2, 35 public patient variants vs 125 L3: AUROC and excess over three rarity-only nulls (post hoc, specified before running; `results/rnu2_position_control.tsv`).

| Score | AUROC | N0 gene-wide excess (p) | N1 region × pairing excess (p) | N2 ±5-nt excess (p) |
|---|---|---|---|---|
| Rfam reference-column frequency | 0.846 | +0.179 (0.0001) | +0.071 (0.009) | +0.041 (0.069) |
| Rfam alternative-base rarity | 0.844 | +0.161 (0.0001) | +0.073 (0.018) | +0.044 (0.088) |
| GPN-MSA | 0.771 | −0.024 (0.74) | −0.050 (0.91) | −0.052 (0.96) |
| phyloP447 | 0.745 | +0.061 (0.09) | +0.001 (0.50) | −0.010 (0.62) |
| CADD v1.7 | 0.555 | +0.019 (0.33) | +0.013 (0.37) | +0.038 (0.10) |
| ReMM v0.4 | 0.514 | −0.013 (0.60) | +0.024 (0.31) | +0.064 (0.06) |
| Baseline: in region n.1–67 | 0.801 | +0.161 (0.0001) | +0.000 (1.0) | +0.005 (0.66) |
| Baseline: 5′ proximity | 0.813 | +0.160 (0.0008) | −0.011 (0.68) | −0.000 (0.54) |

#### RNU4-2 function: the family profile identifies the positions that matter
In RNU4-2 SGE (435 SNVs at 145 positions; no population data):
- **Correlation with functional damage:** reference-column frequency ρ = 0.459 [0.377, 0.539] and alternative-base rarity 0.389, against GPN-MSA 0.143, CADD 0.088, ReMM 0.062 and phyloP447 −0.024.
- **Classification:** reference-column frequency classified the 106 function-depleted variants with an AUROC of 0.802 [0.721, 0.876] (position-cluster interval), against 0.489–0.604 for genome-wide scores. Every paired difference excluded 0, ranging from +0.197 [0.101, 0.284] against GPN-MSA to +0.313 [0.199, 0.430] against phyloP447.
- **Outside the critical region** (132 positions): 0.774 [0.664, 0.868], with every paired difference still excluding 0.
- **Inside the critical region** (31 vs 8): no score differed.

**Within position** (Table 2), alternative-base rarity barely ranked the alternatives (ρ 0.12 [0.01, 0.23]), similar to GPN-MSA (0.14) and above CADD (−0.07). The advantage is therefore position-level: the family profile identifies which positions are functionally constrained, which whole-genome scores fail to do here. It does not identify which substitution at a position is worse.

For reference, the supervised position floor of G1, fitted on neighbouring SGE scores, reached 0.859 on the same SNVs. The parameter-free profile recovers most of that ceiling without using any functional data.

**Table 2.** RNU4-2 SGE, function-depleted variants (< −0.302): AUROC with 95% position-cluster intervals, and within-position Spearman (`results/sge_position_control.tsv`, `results/sge_classes.tsv`).

| Score | AUROC, all (106 vs 329) | AUROC, outside critical region | Within-position ρ |
|---|---|---|---|
| Rfam reference-column frequency | 0.802 [0.721, 0.876] | 0.774 [0.664, 0.868] | constant within position |
| Rfam alternative-base rarity | 0.733 [0.664, 0.804] | 0.722 [0.639, 0.801] | 0.119 [0.009, 0.225] |
| GPN-MSA | 0.604 [0.546, 0.660] | 0.572 [0.504, 0.640] | 0.135 [0.012, 0.249] |
| CADD v1.7 | 0.564 [0.471, 0.654] | 0.510 [0.401, 0.611] | −0.073 [−0.189, 0.046] |
| ReMM v0.4 | 0.513 [0.414, 0.611] | 0.495 [0.367, 0.627] | −0.012 [−0.042, 0.019] |
| phyloP447 | 0.489 [0.408, 0.568] | 0.473 [0.370, 0.579] | constant within position |
| Supervised position floor (G1; uses SGE) | 0.859 | — | — |

A v4 refitted without the RNU4 paralog cluster reached ρ = 0.382 on SGE, computed without the Evo 2 feature. In the published RNU4ATAC splicing assay (23 SNVs),<sup>41</sup> the pattern reversed:
- homology-transferred SGE function ρ = 0.745 (n = 11);
- GPN-MSA 0.581 and CADD 0.472;
- the Rfam features 0.344 and 0.168.

#### Where and why genome-wide scores fail
**Per gene** (Table 3), reference-column frequency matched or exceeded CADD in all five disease genes with public labels. The gain was large where CADD fails, and small and uncertain where CADD works.

**Table 3.** Per-gene AUROC, public P/LP SNVs vs L3 (post hoc, specified before running; `results/per_gene_rfam.tsv`). Agreement = per-position Spearman of phyloP447 with the Rfam profile.

| Gene | P vs L3 | Rfam ref (excess) | CADD (excess) | phyloP447 (excess) | GPN-MSA (excess) | Rfam ref − CADD [95% CI] | Agreement |
|---|---|---|---|---|---|---|---|
| RNU4-2 | 26 vs 52 | 0.913 (+0.29) | 0.519 (+0.08) | 0.515 (−0.05) | 0.825 (+0.03) | +0.394 [0.250, 0.537] | −0.01 |
| RNU2-2 | 35 vs 125 | 0.846 (+0.18) | 0.555 (+0.02) | 0.745 (+0.06) | 0.771 (−0.02) | +0.290 [0.164, 0.409] | 0.39 |
| RNU12 | 8 vs 79 | 0.695 (+0.12) | 0.453 (+0.03) | 0.577 (+0.00) | 0.748 (+0.04) | +0.241 [0.005, 0.432] | 0.36 |
| RNU4ATAC | 37 vs 65 | 0.832 (+0.28) | 0.769 (+0.20) | 0.712 (+0.13) | 0.813 (+0.13) | +0.063 [−0.026, 0.154] | 0.55 |
| RNU6ATAC | 5 vs 18 | 0.744 (+0.12) | 0.706 (+0.15) | 0.644 (+0.06) | 0.672 (−0.05) | +0.039 [−0.350, 0.467] | 0.40 |

These are AUROCs against population-observed negatives, so they carry the rarity and region contrasts shown in Table 1; excess is over the gene-wide null (N0), which does not remove region. GPN-MSA's high AUROCs in RNU4-2 and RNU12 are almost entirely its rarity null (excess +0.03 and +0.04). The Rfam excess is significant in RNU2-2, RNU4-2 and RNU4ATAC (p ≤ 0.001) but not in RNU12 or RNU6ATAC, which have 8 and 5 positives.

**Alignment depth.** At RNU4-2, phyloP447 reaches at most 0.90, and 59% of positions share a single value. Its tandem paralog RNU4-1 reaches 8.03. The alignment itself explains this (Fig. 4A):
- **RNU4-2:** 186 species are aligned, every one a primate, and none of the 204 non-primate mammals in the alignment. The aligned sequences match human at 99.5% identity, so the gene is conserved, but only across short primate branches. That caps phyloP.
- **RNU4-1:** 442 species are aligned, 202 of them non-primate mammals.
- **Other disease genes:** 425–439 species (189–202 non-primates). RNU2-2 has 363 (148 non-primates).

No non-primate information therefore reaches RNU4-2 in the whole-genome alignment. We did not resolve whether non-primate orthologs align to RNU4-1 instead or the duplication is primate-specific.

Across the five genes, the Rfam advantage over CADD was larger where phyloP447 agreed less with the Rfam profile (Spearman −0.80, as predicted; exact one-sided p = 0.067; Fig. 4B). RNU12 is fully aligned (202 non-primates) yet CADD fails there, so alignment explains RNU4-2, in part RNU2-2, and not RNU12. Across eight genes, the fraction of possible SNVs observed in gnomAD did not track CADD's performance (Spearman +0.21).

**Other RNA families (A7, pre-registered; not met).** Rfam alternative-base rarity did not beat genome-wide scores on RMRP, SNORD118, RNU7-1 or TERC: weighted ΔAUROC vs CADD −0.251 [−0.368, −0.133]. On these single-copy RNA genes, genome-wide scores already work (RMRP: CADD 0.863, phyloP447 0.844).

#### Exploratory PP3 thresholds
**Pre-registered (Aim 4, v2).** At prior 0.10, PP3 Supporting was reached:
- in RNU2-2, by the F component (≥ 0.402)‡;
- in RNU4-2 (S ≥ 4.79)‡;
- in RNU4ATAC (S ≥ 0.241)‡.

No Moderate threshold was reached. The calibration of the pre-specified RNU2-2 score S, and the results at prior 0.0441, are on the analysis server and are not yet reported.

**Post hoc (v4, out-of-sample).** PP3 Supporting was reached:
- pooled over genes at v4 ≥ 0.362‡ (99 P vs 1,117 L3; AUROC 0.788 vs CADD 0.746);
- within RNU2-2 at ≥ 0.216‡, RNU4-2 at ≥ 0.189‡ and RNU4ATAC at ≥ 0.88‡.

CADD reached none in RNU2-2 or RNU4-2. We do not report BP4, because population-observed variants of recessive genes can be pathogenic carrier alleles.

These thresholds are exploratory for four reasons:
- they rest on population negatives;
- our implementation did not reproduce published thresholds (A14);
- a pre-specified check against random rare variants (`v4_rare_passrate.py`) is pending;
- likelihood ratios with leave-own-position local constraint have not been computed.

#### Prospective test (A12)
No pre-registered unit has reached its minimum. Units in genes never seen in training have 0–4 new positives; the trained-gene unit has 4 of 5. Rius et al. 2026<sup>8</sup> added no eligible variant, because all its P/LP SNVs were already training labels.

#### Population-negative benchmarks (Supplement S6)
Across 1,787 ClinVar non-coding P/LP SNVs, a random rare variant from the same window reached AUROCs of 0.62–0.68 against population-observed negatives (CADD 0.632, GPN-MSA 0.676). Rarity matching reordered tools on two published benchmarks. Intervals for these orderings, and the cleaning of the null pools, are incomplete, so this part is descriptive. It quantifies a known caveat<sup>25,28,29</sup> and is not a main claim of this paper.

### Discussion
At RNU4-2, RNU2-2 and RNU12, genome-wide predictors carry little information beyond rarity about which variants are pathogenic or functionally damaging. A curated RNA-family alignment recovers position-level conservation that they miss. The clearest evidence is functional and free of population labels. In RNU4-2 SGE, the Rfam profile classified function-depleted variants at 0.80, against at most 0.60, with intervals that respect position clustering. The reason at RNU4-2 is concrete. In the 447-way whole-genome alignment, only primates are aligned there, so phyloP is computed over short branches and capped at 0.90. Its paralog RNU4-1 is aligned across 442 mammals. A curated family alignment does not depend on that orthology assignment.

The patient-variant evidence is weaker than our earlier framing implied. A frozen model separated sealed RNU2-2 patient variants from population variants (0.892 vs CADD 0.501) and met its pre-registered replication on literature variants. But its strength came from Rfam features identified after unsealing. Both label sets come from one gene with shared negatives. Most of the signal beyond rarity is explained by the disease region where patient variants are ascertained. A residual Rfam excess survived matching on region and structure (+0.07), and no genome-wide score showed one.

These results define what an snRNA-specific tool should be built on. Genome-wide scores work on single-copy RNA genes (A7) and in RNU4ATAC and RNU6ATAC, where the Rfam profile adds little. The Rfam profile is position-level conservation. It does not rank the substitutions at a position, and in a small RNU4ATAC assay other scores did better. A practical tool would therefore combine a family-alignment conservation term with substitution-level information, such as structure-aware models or functional data, and use it in the genes where whole-genome conservation is uninformative. An alignment-agreement diagnostic like ours could flag those genes in advance, although five genes cannot establish such a rule (p = 0.067).

The pre-registered primary model failed. Population allele counts at these hypermutable loci did not resolve individual variants, even with UK Biobank added. We report the failure because it is informative, and because the model that succeeded was designed afterwards. Its freeze before unsealing is anchored in Bitcoin, but the claim it supports is narrower than a pre-registered confirmation.

The calibration is not yet fit for clinical use. Population-observed negatives inflate likelihood ratios, and they include carrier alleles in recessive genes. Our implementation of the ClinGen method did not reproduce published thresholds, and the rare-variant check is pending. The PP3 thresholds are exploratory. A calibration suitable for laboratories will need curated benign variants or functional readouts in each gene. This is consistent with recent evidence that predictor-informed classifications feed back into ClinVar,<sup>45</sup> and that genome-wide thresholds do not transfer across variant classes.<sup>20,46</sup>

**Limitations.**
- **Blind evidence from one gene.** The blind evidence comes from one gene (38 positives in two overlapping sets with shared negatives). RNU5B-1 could not be evaluated, and no prospective unit has reached its minimum.
- **Post hoc features.** The Rfam features were identified post hoc, and their first confirmatory test outside the spliceosomal snRNAs (A7) failed. The region-matched and functional analyses are post hoc analyses of public data.
- **Function in one gene, at position level.** The population-free evidence comes from RNU4-2 alone, from a haploid-cell loss-of-function readout. Inside the critical region no score separated variants.
- **Incomplete mechanism.** The alignment explains RNU4-2 (primates only) and partly RNU2-2, but not RNU12, which is fully aligned yet poorly scored by CADD. The gene-level association is five genes (p = 0.067).
- **Record of order.** Amendments after A6 rest on git times, not Bitcoin anchoring, and the seal log is not yet released.
- **Population data at multicopy loci.** gnomAD calls at multicopy snRNA loci may include mapping artefacts (Supplement S10).
- **Missing comparators.** AlphaGenome, covariation scores and RNA structure-effect tools were not evaluated, and we did not build an ortholog-only conservation baseline.
- **Design.** No new functional data were generated, and the study has a single author.

### Conclusions
Whole-genome conservation carries little usable signal beyond rarity at RNU4-2, RNU2-2 and RNU12, and curated RNA-family alignments recover position-level conservation there. This holds against measured RNU4-2 function without population labels, and in part against patient variants after matching on region and structure. snRNA variant-effect tools should build on family alignments rather than whole-genome alignments in these genes. Their clinical calibration needs curated benign or functional data.

### Data availability
- Scores for all 17,385 SNVs (v2 and v4) will be deposited at Zenodo **[DOI to be added]**.
- Curated literature labels (`labels/A8_literature.tsv`, `labels/A12_*.tsv`, with rejected rows and reasons) and all benchmark selections are in the repository.
- The sealed-label event log (`sealed/SEALED_MANIFEST.jsonl`) will be released with an OpenTimestamps proof.

### Code availability
All code, the pre-registration, amendments, OpenTimestamps proofs and the deviation log are at https://github.com/tanukusaitejesh-prog/ndd-pp3-calibration (folder `snrna/`). Every analysis in Tables 1–3 runs from committed public inputs: `rnu2_position_control.py`, `sge_position_control.py`, `sge_classes.py` and `per_gene_rfam.py`.

### Ethics declaration
The study used only public, de-identified data: ClinVar, gnomAD, UK Biobank allele-frequency summaries, published variant reports and public functional data. No ethics approval was required.

### Author contributions, competing interests, funding
[To be completed by the author.] AI assistance: analysis code, literature curation support and drafting were assisted by an AI system under the author's direction. All decisions and the final text are the author's responsibility.

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
42. Tavtigian SV, Greenblatt MS, Harrison SM, Nussbaum RL, Prabhu SA, Boucher KM, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. Genet Med. 2018;20(9):1054-1060. doi:10.1038/gim.2017.210. PMID: 29300386.
43. Caron B, Luo Y, Rausell A. NCBoost classifies pathogenic non-coding variants in Mendelian diseases through supervised learning on purifying selection signals in humans. Genome Biol. 2019;20(1):32. doi:10.1186/s13059-019-1634-2. PMID: 30744685.
44. Benegas G, Eraslan G, Song YS. Benchmarking DNA sequence models for causal regulatory variant prediction in human genetics. bioRxiv [Preprint]. 2025 Feb 12. doi:10.1101/2025.02.11.637758.
45. Tao M, Marsh JA. Variant effect predictors are shaping the ClinVar classifications used to train and evaluate them. bioRxiv [Preprint]. 2026 Oct 27. doi:10.64898/2026.10.24.761512.
46. Babazade A, Dolas Y, Duzenli T, Vural O, Furtana C, Ergun MA. Calibration of AlphaGenome Variant Impact and SpliceAI2 scores as PP3/BP4 evidence for missense, non-coding, and synonymous variants. bioRxiv [Preprint]. 2026 Oct 24. doi:10.64898/2026.10.21.760715.

### Figure legends
**Figure 1.** RNU2-2 patient variants against population variants. ROC curves for Rfam reference-column frequency, GPN-MSA, phyloP447 and CADD on the 35 public patient SNVs vs 125 L3 (AUROC in the legend). The sealed-set curves for frozen v4 require the analysis server and will be added.

**Figure 2.** RNU2-2: excess AUROC over three rarity-only nulls (gene-wide; region × pairing; ±5-nt window) for each score and two position baselines (Table 1). * one-sided p < 0.05.

**Figure 3.** RNU4-2 saturation genome editing. (A) AUROC for function-depleted variants with 95% position-cluster intervals, overall and outside the critical region. (B) Along RNU4-2: Rfam reference-column frequency, phyloP447 (note the flat, low track) and the minimum SGE function score per position.

**Figure 4.** Why genome-wide conservation fails at RNU4-2. (A) Species aligned over at least half of each gene in the UCSC 447-way alignment, split into primates and non-primate mammals: RNU4-2 is aligned only in primates. (B) Per-gene Rfam − CADD AUROC difference (95% CI) against the agreement of phyloP447 with the Rfam profile (Table 3).

‡ Transcribed from analysis-server outputs into the dated deviation log at the time they were produced; the files will be released with the scores.

---

# DOCUMENT 2: SUPPLEMENT

## Supplementary Information

*RNA-family profiles provide blind-validated computational evidence for spliceosomal snRNA variants. Every number below is from the dated deviation log (`snrna/PREREG.md` §14) or from committed result files.*

### S1. Timeline of freezes and results
The pre-registration, code and amendments A1–A6 are anchored in Bitcoin blocks (`paper/TIMELINE.md`). Amendments A7–A16 and every post hoc analysis rest on git commit times only. The seal/unseal log (`sealed/SEALED_MANIFEST.jsonl`) is on the analysis server and will be released with an OpenTimestamps proof.

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

### S3b. Region- and structure-matched rarity nulls on RNU2-2 (post hoc, specified before running; after an external review)
`rnu2_position_control.py` (committed ed12221 before the result dad9b58); `results/rnu2_position_control.tsv`. Positives: the 35 public RNU2-2 SNVs accepted by the frozen A8 rule (16 A8-new, 19 also in the sealed set). Negatives: 125 RNU2-2 L3. Pool: 413 SNVs that are neither. 10,000 draws per null. N0 gene-wide; N1 matched on Table 2 region status × reference pairing; N2 within ±5 nt. Inside n.1–67: 25 vs 14; outside: 10 vs 111.

| Score | AUROC | Rarity-matched AUROC | N0 null / excess (p) | N1 null / excess (p) | N2 null / excess (p) | AUROC inside | AUROC outside |
|---|---|---|---|---|---|---|---|
| Rfam ref column freq | 0.846 | 0.710 | 0.667 / +0.179 (0.0001) | 0.775 / +0.071 (0.0088) | 0.805 / +0.041 (0.0686) | 0.587 | 0.811 |
| Rfam alt rarity | 0.844 | 0.696 | 0.683 / +0.161 (0.0001) | 0.771 / +0.073 (0.0176) | 0.800 / +0.044 (0.0884) | 0.660 | 0.811 |
| CADD v1.7 | 0.555 | 0.514 | 0.536 / +0.019 (0.3304) | 0.543 / +0.013 (0.3689) | 0.517 / +0.038 (0.1019) | 0.483 | 0.578 |
| phyloP447 | 0.745 | 0.532 | 0.684 / +0.061 (0.0896) | 0.744 / +0.001 (0.5027) | 0.755 / -0.010 (0.6202) | 0.609 | 0.808 |
| GPN-MSA | 0.771 | 0.452 | 0.795 / -0.024 (0.7423) | 0.821 / -0.050 (0.9129) | 0.823 / -0.052 (0.9552) | 0.591 | 0.769 |
| ReMM v0.4 | 0.514 | 0.485 | 0.526 / -0.013 (0.6001) | 0.489 / +0.024 (0.3072) | 0.450 / +0.064 (0.0607) | 0.396 | 0.772 |
| BASELINE in region n.1-67 | 0.801 | 0.661 | 0.640 / +0.161 (0.0001) | 0.801 / +0.000 (1.0000) | 0.796 / +0.005 (0.6567) | 0.500 | 0.500 |
| BASELINE distance to region | 0.790 | 0.652 | 0.641 / +0.149 (0.0012) | 0.802 / -0.011 (0.6904) | 0.790 / +0.000 (0.4579) | 0.500 | 0.458 |
| BASELINE 5' proximity | 0.813 | 0.648 | 0.653 / +0.160 (0.0008) | 0.823 / -0.011 (0.6791) | 0.813 / -0.000 (0.5406) | 0.776 | 0.458 |
| BASELINE position unobserved | 0.731 | 0.517 | 0.717 / +0.014 (0.3856) | 0.817 / -0.086 (0.9939) | 0.833 / -0.102 (0.9989) | 0.691 | 0.497 |

Reading: under N0 the region baseline gets as much excess as Rfam, so the N0 excess cannot separate pathogenicity from ascertainment. Under N1 only the Rfam features keep a significant excess. Under N2 every excess is ≤ +0.064 and none is significant. "Position unobserved" counts the negative's own allele and is optimistic by construction.

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

**Position-cluster bootstrap and within-position test** (post hoc, specified before running after the external review: `sge_position_control.py`, committed before its result; `results/sge_position_control.tsv`). Reference-column frequency is constant within a position, so the SNV-level intervals above are too narrow. Resampling the 145 positions (132 outside the critical region), 2,000 draws:

| Subset | Test | Score | Value [95% CI] |
|---|---|---|---|
| all | cluster AUROC | Rfam ref column freq | 0.802 [0.721, 0.876] |
| all | cluster AUROC | CADD v1.7 | 0.564 [0.471, 0.654] |
| all | cluster diff (ref - score) | CADD v1.7 | +0.237 [+0.130, +0.350] |
| all | cluster AUROC | GPN-MSA | 0.604 [0.546, 0.660] |
| all | cluster diff (ref - score) | GPN-MSA | +0.197 [+0.101, +0.284] |
| all | cluster AUROC | phyloP447 | 0.489 [0.408, 0.568] |
| all | cluster diff (ref - score) | phyloP447 | +0.313 [+0.199, +0.430] |
| all | cluster AUROC | ReMM v0.4 | 0.513 [0.414, 0.611] |
| all | cluster diff (ref - score) | ReMM v0.4 | +0.288 [+0.188, +0.386] |
| all | cluster AUROC | Rfam alt rarity | 0.733 [0.664, 0.804] |
| outside critical region | cluster AUROC | Rfam ref column freq | 0.774 [0.664, 0.868] |
| outside critical region | cluster AUROC | CADD v1.7 | 0.510 [0.401, 0.611] |
| outside critical region | cluster diff (ref - score) | CADD v1.7 | +0.265 [+0.129, +0.405] |
| outside critical region | cluster AUROC | GPN-MSA | 0.572 [0.504, 0.640] |
| outside critical region | cluster diff (ref - score) | GPN-MSA | +0.202 [+0.076, +0.319] |
| outside critical region | cluster AUROC | phyloP447 | 0.473 [0.370, 0.579] |
| outside critical region | cluster diff (ref - score) | phyloP447 | +0.302 [+0.157, +0.447] |
| outside critical region | cluster AUROC | ReMM v0.4 | 0.495 [0.367, 0.627] |
| outside critical region | cluster diff (ref - score) | ReMM v0.4 | +0.279 [+0.150, +0.412] |
| outside critical region | cluster AUROC | Rfam alt rarity | 0.722 [0.639, 0.801] |

Within position (score and −function score centred within each position; position-cluster intervals):

| Score | Within-position Spearman |
|---|---|
| Rfam alt rarity | +0.119 [+0.009, +0.225] |
| CADD v1.7 | -0.073 [-0.189, +0.046] |
| GPN-MSA | +0.135 [+0.012, +0.249] |
| ReMM v0.4 | -0.012 [-0.042, +0.019] |
| ViennaRNA ddG_bind | +0.117 [-0.002, +0.234] |

Reading: the Rfam advantage survives position clustering but is position-level; within a position, alternative-base rarity ties GPN-MSA. The G1 supervised position floor (fitted on neighbouring SGE scores, leave-one-position-out) reached ρ 0.518 and AUROC 0.859 on these SNVs (`results/g1_report.txt`).

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

Across the five genes with CADD estimates, agreement tracked CADD AUROC (Spearman +0.90); the ceiling did not (+0.20). At RNU4-2, phyloP447 is flat at 0.687 across n.60–82, one of the most conserved stretches of U4 snRNA. RNU12 does not fit the low-ceiling pattern. Five genes; descriptive only.

**Direct check of the 447-way alignment** (`alignment_species.py`, specified before running; `alignment_species_clades.py`, descriptive extension written after the counts; `results/alignment_species*.tsv`). Species aligned over ≥ 50% of the gene:

| Gene | Aligned species | Primates | Non-primates | Median species per position (min) | Identity to human |
|---|---|---|---|---|---|
| RNU4-2 | 186 | 186 | 0 | 186 (177) | 0.995 |
| RNU4-1 | 442 | 240 | 202 | 442 (435) | 0.975 |
| RNU2-2 | 363 | 215 | 148 | 359 (193) | 0.959 |
| RNU12 | 431 | 229 | 202 | 430 (389) | 0.968 |
| RNU4ATAC | 439 | 240 | 199 | 438.5 (0) | 0.923 |
| RNU6ATAC | 432 | 235 | 197 | 432 (14) | 0.970 |
| RNU5B-1 | 425 | 236 | 189 | 418 (399) | 0.949 |

The tree holds 242 primates (excluding human) and 204 non-primate mammals. The paralog explanation of the flat RNU4-2 phyloP track is therefore confirmed in the sense that no non-primate sequence is aligned to RNU4-2.

### S5c. Per-gene comparison on public labels (post hoc, specified before running)
`per_gene_rfam.py`; `results/per_gene_rfam.tsv`. AUROC (excess over the gene-wide rarity null, 2,000 draws) against the gene's L3.

| Gene | P vs L3 | Rfam ref | Rfam alt | CADD | phyloP447 | GPN-MSA | Rfam ref − CADD [95% CI] |
|---|---|---|---|---|---|---|---|
| RNU2-2 | 35 vs 125 | 0.846 (+0.180, p 0.0005) | 0.844 (+0.162) | 0.555 (+0.021) | 0.745 (+0.063) | 0.771 (-0.022) | +0.290 [+0.164, +0.409] |
| RNU4-2 | 26 vs 52 | 0.913 (+0.288, p 0.0005) | 0.828 (+0.110) | 0.519 (+0.077) | 0.515 (-0.053) | 0.825 (+0.027) | +0.394 [+0.250, +0.537] |
| RNU4ATAC | 37 vs 65 | 0.832 (+0.279, p 0.0005) | 0.803 (+0.129) | 0.769 (+0.204) | 0.712 (+0.133) | 0.813 (+0.131) | +0.063 [-0.026, +0.154] |
| RNU12 | 8 vs 79 | 0.695 (+0.122, p 0.1304) | 0.750 (+0.078) | 0.453 (+0.030) | 0.577 (+0.002) | 0.748 (+0.037) | +0.241 [+0.005, +0.432] |
| RNU6ATAC | 5 vs 18 | 0.744 (+0.120, p 0.1919) | 0.689 (+0.027) | 0.706 (+0.154) | 0.644 (+0.059) | 0.672 (-0.054) | +0.039 [-0.350, +0.467] |

Across the five genes, the Rfam − CADD advantage correlated with phyloP447–Rfam agreement at Spearman −0.80 (prediction: negative; exact one-sided permutation p = 0.067).

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
- **BP4 not reported in the main text:** population-observed variants of recessive genes (RNU2-2, RNU4ATAC, RNU12, RNU6ATAC) can be pathogenic carrier alleles, so BP4 against L3 is unsafe. The logged BP4 values (pooled v4 ≤ −1.35; v2 RNU4ATAC S ≤ −0.046) are kept here for completeness only.
- **Implementation check (A14):** our local-LR procedure gave AlphaMissense thresholds 0.627/0.841/0.985 against the published 0.792/0.906/0.972 (Bergquist et al.), so thresholds from it are approximate.
- **Pending (analysis server):** the calibration of the pre-specified RNU2-2 score S, all results at prior 0.0441, the rare-variant pass-rate check, and likelihood ratios with leave-own-position local constraint.
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

### 6b. External review (5 Nov 2026) and what changed

An external AI review of the 4 Nov packet listed fatal, major and minor issues. Its factual claims were checked against the log and hold. Changes made without the analysis server (each new analysis committed before it was run):

| Review point | Response |
|---|---|
| Title and abstract oversell "pre-registered, blind-validated evidence" | Retitled around what survives: whole-genome conservation fails at key snRNA genes and family alignments recover position-level signal; v4 described as a post-failure model; Rfam features stated as identified after unsealing |
| Rarity null confounded by region and position ascertainment | New `rnu2_position_control.py` on 35 public RNU2-2 positives. The region baseline alone matches the Rfam N0 excess (+0.16 vs +0.18), which confirms the point. After region × pairing matching, Rfam keeps +0.07 (p 0.009) and genome-wide scores ≈ 0; under ±5-nt matching it is a trend (+0.04, p 0.07). Reported in Table 1 |
| SGE bootstrap treats SNVs as independent | New `sge_position_control.py`: position-cluster intervals (Rfam 0.80 [0.72, 0.88]; every difference still excludes 0, also outside the critical region). Within position, Rfam ties GPN-MSA (0.12 vs 0.14), so the claim is stated as position-level |
| Supervised position floor (G1) omitted | Reported: AUROC 0.859 on the same SNVs, as a ceiling that uses SGE data |
| G1, G2, A1 outcomes absent | Reported in Results |
| Pre-specified RNU2-2 score S not shown | S = 0.784 reported; its calibration flagged as pending on the server |
| A14 validity failure not in paper | Reported in Methods and Supplement S8 |
| BP4 unsafe with population negatives in recessive genes | BP4 removed from the main text; PP3 thresholds labelled exploratory |
| Rfam never shown for RNU12, RNU6ATAC | New `per_gene_rfam.py`: Rfam ≥ CADD in all 5 genes; gain where CADD fails (RNU2-2, RNU4-2, RNU12); ties in RNU4ATAC, RNU6ATAC; alignment-agreement prediction ρ −0.80, p 0.067 |
| Benchmark-inflation section is weak and secondary | Moved to Supplement S6; one descriptive paragraph in the main text |
| OpenTimestamps claimed for all amendments; seal log absent | Methods now state A1–A6 Bitcoin-anchored, A7–A16 git only, seal log pending release with a proof |
| Header claimed every number traces to committed files | Server-transcribed numbers now marked ‡ |
| Missing comparators (AlphaGenome, covariation, structure tools, ortholog-only baseline) | Listed as limitations; not run |
| Paralog explanation for RNU4-2 unchecked against the alignment | New `alignment_species.py` reads the UCSC 447-way alignment directly: RNU4-2 is aligned in 186 species, all primates (0 of 204 non-primates), against 442 species at RNU4-1. The mechanism is now a demonstrated alignment fact for RNU4-2; it does not explain RNU12 |
| Figures not drawn, Table 4 empty cells | Figures 1–4 drawn from public inputs (`figures_paper.py`); Figure 1 uses the 35 public RNU2-2 variants because the sealed-set v4 scores are on the server |

Still requiring the analysis server: v4 under the matched nulls on the sealed set; the rare-variant threshold check; the calibration of S and prior-0.0441 results; leave-own-position likelihood ratios; the seal log with an OpenTimestamps proof; corrected benchmark intervals.

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
42. Tavtigian SV, Greenblatt MS, Harrison SM, Nussbaum RL, Prabhu SA, Boucher KM, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. Genet Med. 2018;20(9):1054-1060. doi:10.1038/gim.2017.210. PMID: 29300386.
43. Caron B, Luo Y, Rausell A. NCBoost classifies pathogenic non-coding variants in Mendelian diseases through supervised learning on purifying selection signals in humans. Genome Biol. 2019;20(1):32. doi:10.1186/s13059-019-1634-2. PMID: 30744685.
44. Benegas G, Eraslan G, Song YS. Benchmarking DNA sequence models for causal regulatory variant prediction in human genetics. bioRxiv [Preprint]. 2025 Feb 12. doi:10.1101/2025.02.11.637758.
45. Tao M, Marsh JA. Variant effect predictors are shaping the ClinVar classifications used to train and evaluate them. bioRxiv [Preprint]. 2026 Oct 27. doi:10.64898/2026.10.24.761512. †
46. Babazade A, Dolas Y, Duzenli T, Vural O, Furtana C, Ergun MA. Calibration of AlphaGenome Variant Impact and SpliceAI2 scores as PP3/BP4 evidence for missense, non-coding, and synonymous variants. bioRxiv [Preprint]. 2026 Oct 24. doi:10.64898/2026.10.21.760715. †

---

# END OF MATERIALS
Now complete tasks 1–8 from the instructions at the top.
