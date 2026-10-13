# Methods (draft, written before the model results; numbers in brackets are filled from `results/`)

## Study design and pre-registration

The analysis plan, operational definitions and decision gates were frozen before any endpoint was computed (`PREREG.md`). The frozen text and code were hashed, and the hashes were timestamped with OpenTimestamps (Bitcoin-anchored) at four points:
- before the first analysis (`prereg_manifest.txt`, `code_manifest.txt`);
- before the G2 test (amendment A1, and the mutation-model freeze `g2_freeze_manifest.txt`);
- before any model was fitted to real allele counts (amendment A2, `amendment_A2_manifest.txt`);
- before any model result existed, for the evaluation-only code (amendment A3).

Every later change is listed in a dated deviation log, with what had been seen when it was made. Three decision gates were pre-specified:
- **G1:** is a genomic foundation model already sufficient?
- **G2:** do population data carry a selection signal?
- **G3:** does the model beat CADD on held-out genes?

The two blind genes, RNU2-2 and RNU5B-1, had their patient labels compiled straight into an AES-256-encrypted file, on the GPU machine, by the authors. That happened before any selection model was fitted. The analyst never held the passphrase, and the file was decrypted once, for the final evaluation (logged).

## Data

- **Genes.** GENCODE v49 spliceosomal snRNA genes on chromosomes 1–22, X and Y gave 111 genes. Of these, [87] passed gnomAD v4.1 genome coverage (median depth ≥ 20, minimum ≥ 10) and [55] passed paralog QC. Paralog QC removed genes with more than half their bases in GIAB low-mappability or segmental-duplication regions, or with mostly low allele-balance heterozygous calls.
- **Variants.** All 3 possible SNVs at every transcript base ([17,385] SNVs), keyed by GRCh38 chrom-pos-ref-alt.
- **Mutation rates** came from Roulette.
- **Allele counts** came from gnomAD v4.1 genomes (FILTER PASS), with per-ancestry counts for the drift analysis.
- **CADD** v1.7.
- **RNU4-2 function scores** came from saturation genome editing (SGE): 521 in-transcript variants, 435 of them SNVs.
- **Comparators:**
  - phyloP for 100 vertebrates and 447 mammals;
  - ReMM v0.4;
  - GPN-MSA;
  - Evo 2 7B: the summed log-likelihood difference over an 8,192-nt window, averaged over both strands.

## Label sets

- **L2, patient variants:** distinct variants reported as pathogenic or likely pathogenic, either by ClinVar with at least one review star or by RNUdb curation, or by the reporting study (literature tally, mapped to GRCh38 with reference-base checks).
- **L3, benign-like:** gnomAD PASS variants with AC ≥ 9, minus any L2 variant.

Presence in gnomAD alone never made a variant benign.

## Selection signal (gate G2)

For each site, observation was presence in gnomAD at AF < 0.001. Expectation was p_v = 1 − exp(−κ_g m_v), with a per-gene scale κ_g. A structural correction was rejected by leave-one-gene-out log-loss on training genes (frozen before testing). Two tests were used:
- **RNU4-2 critical region:** an exact Poisson-binomial test against the gene-median-scaled expectation.
- **Blind hotspots:** a hybrid rule (amendment A1, adopted on simulation evidence before any G2 data): exact-test p < 0.01, a top-decile window inside the hotspot, and circular-shift p < 0.05.

## The model

**Graph per gene.** Each locus is a graph: the target snRNA's nucleotides plus its spliceosomal partner RNAs.
- U4 with U6.
- U6 with U4 and U2.
- U2 with U6 and the branch-site consensus.
- U1 with the 5′ splice-site consensus.
- The minor-spliceosome equivalents.

Edges are of three kinds:
- the backbone;
- intramolecular base pairs from the Rfam 15.1 consensus (Infernal cmalign);
- intermolecular pairs from RNAcofold (ViennaRNA 2.7.0).

**Inputs, all label-free:**
- node features: base, molecule, Rfam profile and entropy, insert state, helix geometry and loop class;
- variant features: alternative base, transition/CpG, consensus column and pair frequencies, broken intermolecular pairs, ΔΔG of folding and of partner binding;
- Evo 2 embeddings of each position (layer 28, both strands, PCA to 64 dimensions) and each SNV's Evo 2 log-likelihood difference;
- gene features: the gene's mutation-rate scale, Rfam fit, HGNC-named flag and family.

**Encoder.** L layers of typed message passing (GINE), then global self-attention with learned per-head biases for edge type and shortest-path distance. A variant head maps [node state, variant features] to a selection score f(v) ≥ 0 (softplus). Its output layer is zero-initialised, so any variation in f is learned from the counts.

**Gate.** A logistic gate γ_g ∈ [0, 1] scales f per gene, allowing unexpressed copies to show no selection.

**Likelihood.** Allele counts follow a negative binomial with log μ_v = α_g + log m_v + β s_v − γ_g f(v):
- α_g is profiled out per gene, so only within-gene departures from the mutation model carry signal;
- s_v, a site-rate covariate, enters only the offset.

f never sees allele counts as input.

**Training.**
- Full-batch AdamW (lr 3·10⁻³) with gradient clipping.
- Early stopping on the negative-binomial log-likelihood of held-out validation clusters.
- Hyperparameters chosen by 3-fold cluster cross-validation of that log-likelihood, over d ∈ {64, 128}, layers ∈ {2, 3}, dropout ∈ {0.1, 0.3} and weight decay ∈ {10⁻⁴, 10⁻²}.
- Five seeds; the score is their mean.

No SGE or patient label was used for f.

**Holdout.** The unit is the paralog cluster (single linkage at ≥ 90% identity within a family), which is stricter than leave-one-gene-out. Every named gene is scored by a model trained without its whole cluster. The strict model, trained without the U4 family, scores RNU4-2 for Aim 2.

**Function score F.** For recessive disease genes, selection on heterozygotes is weak, so we also defined F on the SGE function scale.
- For RNU4ATAC positions mapped to RNU4-2 by agreement of two structure-aware alignments (U4 and U4atac covariance models), F is the measured RNU4-2 SGE score of the same change. The mapping covers 38 positions in the 5′ stem-loop/k-turn and stem I.
- Elsewhere, F is a ridge regression from the model's frozen representations to the RNU4-2 SGE score, fitted only in models that never saw the gene scored.

**Framework score S.** S = f for dominant genes, F for recessive genes, and the mean within-gene rank of both for mixed genes (pre-specified; amendment A2).

## Evaluation

**Aim 2.** Concordance with RNU4-2 SGE on 435 SNVs, using Spearman ρ and AUROC for depleted variants. The ReNU-versus-population AUROC compares 12 ReNU SNVs against gnomAD-present SNVs. The strict model must exceed the best of CADD, ΔΔG and Evo 2 by 0.10 in ρ, or reach ReNU AUROC ≥ 0.85.

**Ablations and controls.**
- Ablations: twelve on the strict training set.
- Negative controls: shuffled counts within gene, a randomly rewired pairing graph, and pseudogenes only. Each must give |ρ| < 0.1 and an AUROC interval containing 0.5.

**Gate G3.** Per-gene AUROC of L2 against L3 SNVs on held-out genes, weighted by the number of L2 SNVs, compared with CADD. It passes if Δ ≥ 0.15 and the 95% bootstrap interval (2,000 replicates within gene × class) excludes 0.

**Blind genes.** RNU2-2 and RNU5B-1 were scored by models trained without their clusters (RNU5B-1 without RNU5A-1 and RNU5E-1). They were evaluated once, after unsealing.

## Calibration (ClinGen PP3/BP4)

We followed the ClinGen local-likelihood-ratio approach (Pejaver et al. 2022).

**Local LR.**
- Evaluated at every observed score, over the nearest 20% of labelled variants by score rank, with at least 10 per class and 0.5 pseudocounts.
- Bootstrap: 1,000 class-stratified replicates.

**Evidence thresholds.**
- Very-Strong odds: O_VSt = (9(1−P)/P)^(4/3) at prior P = 0.10, with sensitivity at 0.01, 0.0441, 0.2 and 0.5.
- Supporting, Moderate and Strong need LR ≥ O^(1/8), O^(1/4) and O^(1/2).
- PP3 thresholds require the 5th bootstrap percentile to clear the bar at every higher score; BP4 mirrors this with the 95th percentile.

**Calibration units.**
- per gene;
- dominant and recessive pools;
- dominant and recessive positives within a gene.

**Robustness.**
- with and without positives whose classification used in-silico evidence;
- per gnomAD ancestry group, with L3 redefined as AC_g ≥ 9·AN_g/AN.

## Code and data availability

All code, inputs, timestamps and the deviation log are in the GitHub repository (`snrna/`). Scores for every SNV were released on Zenodo [DOI, date] before the final evaluation. Patient variants published after that date form a pseudo-prospective test set.
