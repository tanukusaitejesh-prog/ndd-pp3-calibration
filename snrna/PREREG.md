# snRNA-VEP v2: pre-registration (week 0)

**Status.** Frozen on 13 October 2026, before any endpoint below was computed. The plan behind it is `PLAN.md`; this document turns that plan into exact, executable definitions. Where `PLAN.md` was ambiguous or internally inconsistent, the resolution is stated here and listed in section 13, so it can be checked before the analysis it affects is run.

**Timestamp.** The SHA-256 of this file and of the analysis code at freeze time are listed in `prereg_manifest.txt`, which is stamped with OpenTimestamps (`prereg_manifest.txt.ots`, Bitcoin-anchored) and committed to the repository. An OSF deposit of the same files is an author action still to do; if it happens after any analysis, the OSF record is a later copy and must not be described as predating the analysis.

**After freezing:**

- Only bug fixes, and the clarifications that section 13 explicitly allows, may change anything below.
- Each change gets a dated entry in the deviation log (section 14), saying whether any result had been seen when it was made.
- A threshold, split, label rule or test changed after its data were seen is reported as post hoc.

## 0. What had been seen before freezing

Honesty about prior exposure matters more than its absence. Before this document was written:

- **SGE table.** Supplementary Table 1 of the SGE paper was opened. Counts and categories were tabulated: 539 variants, 151 below −0.302, 72 below −0.9, 521 in-transcript, 18 ReNU-syndrome variants (12 SNVs), 353 population-observed (346 SNVs). The positions of the ReNU variants were listed; they define the critical region and are published. The paper itself (read in full by the plan's author) reports the CADD, ΔΔG and SGE AUCs. **No score was correlated with any label by this project.**
- **Crick repository** (commit `b404bf5`). Region boundaries and the published AUC code were read.
- **Access tests.** gnomAD v4.1 at RNU4-2 (3 records printed, from an access test), CADD v1.7 at one position, Roulette at two sites, gnomAD coverage at a few RNU4-2 bases.
- **Gene sets.** The GENCODE v49 gene list (111 spliceosomal genes reproduced exactly) and gene lengths.
- **RNUdb.** Only the gene list (RNU2-2, RNU4-2) and empty RNU4-2 region tracks. No variant endpoint was called for any gene.
- **Never accessed.** No RNU2-2 or RNU5B-1 patient-variant data has been downloaded, opened or tallied by this project. The two blind hotspot coordinates are known, because they are in `PLAN.md`, quoted from Table 2 of the guidance.

## 1. Data, pinned

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

## 2. Gene universe and QC

1. **Spliceosomal genes.** GENCODE v49 genes with `gene_type "snRNA"`, on chr1–22, X or Y, whose `gene_name` matches `^(RNU1-\d+|RNVU1-\w+|RNU2-\d+|RNU4-\d+|RNU4ATAC|RNU5[A-F]-\d+|RNU6-\d+|RNU6ATAC|RNU11|RNU12|U1|U2|U4|U6|U4atac|U6atac|U11|U12)$` and does not end in `P`. This gives 111 genes, matching the guidance's count.
2. **Coverage.** gnomAD genome per-base median depth over every transcript base. Keep a gene if the median of those per-base medians is ≥ 20 and their minimum is ≥ 10. The guidance reports 87 genes. A different count is reported, not tuned.
3. **Paralog QC.** Exclude a gene if either holds:
   - more than 50% of its bases overlap the GIAB low-mappability-and-segmental-duplication union;
   - it is mostly low-VAF: summed over its PASS variants, more than 50% of heterozygous genotypes fall in allele-balance bins below 0.25 (gnomAD `ab_hist_alt_bin_freq`).

   RNU4-2/RNU4-1 is flagged as a high-risk pair but not excluded on that ground alone.
4. **Analysis set.** The genes passing 2 and 3. Every exclusion is listed with its reason.

## 3. Variants and coordinates

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

## 4. Label sets

| Set | Exact definition | Use |
|---|---|---|
| L1 SGE | `function_score_mean` (haploid) from SGE ST1; **depleted** = score < −0.302; **strong** = score < −0.9; `function_score_dip` (diploid) for sensitivity only. In-transcript = `nucleotide` 1–145 (521 variants: 435 SNVs, 86 indels) | Evaluation only |
| L2 Patient P/LP | Distinct variants (not individuals) reported as pathogenic or likely pathogenic by the reporting study, or by ClinVar with ≥ 1 review star, as of the tally date. Tagged by gene, Table 2 region and inheritance (dominant, recessive, both, unknown). RNU2-2 and RNU5B-1 entries go only into the sealed file (section 5) | Evaluation, calibration |
| L3 Benign-like | gnomAD v4.1 genomes, FILTER PASS, AC ≥ 9 (≈ AF 5.5e-5 × AN ≈ 152,000, scaling the SGE paper's rule). Sensitivity at AC ≥ 5, ≥ 20, ≥ 50. A variant in both L2 and L3 is removed from L3 and listed | Calibration |
| L4 Unlabeled | Every other possible SNV. Never called benign | Selection training |

Presence in gnomAD alone never makes a variant benign.

### RNU4-2 regions

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

### Blind hotspots

| Gene | Hotspot (as quoted in `PLAN.md` from guidance Table 2) |
|---|---|
| RNU2-2 | n.1–67 |
| RNU5B-1 | n.37–44 |

Section 13 allows one pre-specified correction if Table 2 itself gives different coordinates.

## 5. Splits and the sealed labels

- **Sealed genes.** RNU2-2 and RNU5B-1. Their L2 variants are compiled into `sealed/blind_labels.tsv` and encrypted with `snrna.py seal`, which uses AES-256 with a passphrase held by someone other than the analyst where possible. The SHA-256 of the plaintext and of the ciphertext is recorded in a dated amendment **before any selection model is trained**.
  - The pipeline refuses to read the file without `--unseal`, which logs the date, user and reason.
  - The file is unsealed once, for the final evaluation (section 9).
- **Other splits:**
  - leave-one-gene-out over the analysis set;
  - leave-family-out, U4 + U4atac: RNU4-1, RNU4-2, RNU4ATAC and the unnamed U4 loci held out together;
  - leave-region-out within RNU4-2, using the regions in section 4.
- **Training exclusions.** The two sealed genes are excluded from every training run until the final run. Hyperparameters are chosen only on gene-held-out folds of the unlabeled (L4 plus counts) data.

## 6. Gate G1 (week 1): is a foundation model already enough?

### 6.1 Evo 2 score

- **Model and window.** `evo2_7b`. For each SGE variant, take the GRCh38 window of 8,192 nt centred on the variant's first reference base, chr12:c−4096 to c+4095. The alternative window is the same span with the variant applied, so indels change its length.
- **Score.** LL is the sum of per-token log-likelihoods. Δ = LL(alt) − LL(ref), computed on the forward strand and on the reverse complement of both windows, then averaged.
- **Damage score.** d = −Δ.
- **Sensitivity.** The length-normalised mean log-likelihood.

### 6.2 Damage scores of the comparators (higher = more damaging)

| Predictor | Damage score d |
|---|---|
| CADD | v1.7 PHRED. Primary: the SGE table's `CADD_score` column. QC check: it must agree with the CADD v1.7 file to ±0.01 for SNVs |
| ViennaRNA U4/U6 ΔΔG | The SGE table's `ddG_bind` column, primary. Our 2.7.0 recomputation with the Crick parameters is a QC check |
| Allele-count floor | −log10(gnomAD v4.1 AC + 1), with AC = 0 if absent or non-PASS |
| Position-only floor | Minus the leave-one-position-out prediction of `function_score_mean` from an OLS model on two features: a one-hot of the nine RNU4-2 regions, and log(1 + distance in nt to the nearest critical-region position). All variants at a position are held out together |
| EVEE | Looked up for the ClinVar-overlapping subset only; reported, never used in a gate |

### 6.3 Metrics

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

### 6.4 ReNU-versus-population AUROC

- **Primary (as `PLAN.md` Appendix A specifies).** The 12 ReNU SNVs against in-transcript SGE SNVs that are not ReNU and are present in gnomAD v4.1 genomes (FILTER PASS, AC ≥ 1). The AC floor is undefined here, because controls are selected on AC, and is not reported.
- **Reproduction check.** The 12 ReNU SNVs against the paper's 346 population SNVs (`category == "UKBB/AllofUs"`, `Type == "SNV"`), using the paper's code definitions. The published values are CADD 0.65, ΔΔG 0.72 and SGE 0.93; failing to reproduce them to ±0.01 is a pipeline error.

### 6.5 Decision

| Outcome | Rule | Consequence |
|---|---|---|
| **G1 FIRES** | Evo 2 ρ ≥ 0.60 **and** AUROC (depleted vs not) ≥ 0.85, point estimates, primary definitions | Drop the new-model framing; the paper becomes calibration plus failure analysis |
| **G1 PASSES** | Otherwise | Continue |
| Position floor flag | If the position-only floor's ρ is within 0.05 of the best non-floor predictor's ρ, the deep-model framing inside RNU4-2 ends (`PLAN.md`), whatever G1 says | Reported |
| **Pending** | If Evo 2 cannot be run, G1 stays pending | Everything else in section 6 is still reported |

## 7. Gate G2 (by week 4): is there a selection signal?

### 7.1 Observed

For each possible SNV v in an analysis-set gene, O_v = 1 if gnomAD v4.1 genomes has v with FILTER PASS, AC ≥ 1 and AF < 0.001; otherwise O_v = 0.

Excluded from both O and E:

- sites whose gnomAD record for that allele is non-PASS (sensitivity: count them as O_v = 0);
- sites with AF ≥ 0.001;
- sites whose Roulette FILTER is `low`, or which have no Roulette record.

### 7.2 Expected

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

### 7.3 Windows and statistics

- **O/E.** For a set of sites S, O/E_S = Σ_S O_v / Σ_S p_v.
- **Windows.** Every run of 5 consecutive transcript positions, step 1.
- **Window depletion.** D_w = P(X ≤ O_w), where X follows the Poisson-binomial distribution with the window's p_v. Smaller D_w means more depleted.
- **Top decile.** The 10% of windows with the smallest D_w.
- **Median.** m_g is the median window O/E of gene g.

### 7.4 Tests

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

### 7.5 Decision

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

## 8. Selection model and Aims 2–3 (weeks 5–8)

### 8.1 Training

**Count model.** A_v ~ NegBin(λ_v, φ), with λ_v = μ_v T e^{−f(v)} and μ_v = m_v κ_g c_k, over L4 and observed sites of training genes. f(v) ≥ 0 is predicted by the encoder of `PLAN.md` Methods 3.

**Inputs:**

- sequence ±100 nt;
- Rfam profile and covariation;
- the base-pairing graph, including U4/U6 and U4atac/U6atac intermolecular pairs;
- element class;
- the site-rate covariate, computed excluding v itself.

**What training never sees.** SGE, L2 and L3 labels are never used in training, model selection or early stopping. Early stopping and hyperparameters use gene-held-out NB log-likelihood only.

**Uncertainty.** Five seeds; the score is the ensemble mean of f, and the interval is the range across seeds.

### 8.2 Aim 2: SGE concordance

**Strict model.** Trained with the whole U4 family removed. The loose version removes RNU4-2 only and is reported as secondary.

**Met if either holds:**

- ρ (section 6.3) ≥ the best of CADD, ΔΔG and Evo 2 by at least 0.10;
- AUROC for ReNU against population (section 6.4 primary) ≥ 0.85.

**If Evo 2 is pending,** the comparison uses CADD and ΔΔG and says so.

### 8.3 Aim 3 and gate G3 (week 6)

**Held-out sets.** Each gene is scored by a model trained without it. The week-6 set is RNU4ATAC (leave-one-gene-out and leave-family-out), plus RNU6ATAC, RNU12 and the named RNU6 paralogs (pooled as one unit, since their sequences are identical), each only if it has ≥ 10 distinct L2 SNVs. The sealed genes are not in the week-6 set (section 13, item 3).

**Metric.** Per-gene AUROC of L2 SNVs against L3 SNVs. The summary is the mean over genes weighted by the number of L2 SNVs, for the model and for CADD. ΔAUROC is model minus CADD.

**Uncertainty.** 2,000 bootstrap replicates resampling distinct variants within gene × class.

| Outcome | Rule | Consequence |
|---|---|---|
| **G3 PASS** | Weighted ΔAUROC ≥ 0.15 and the 95% interval excludes 0 | Full framework paper |
| **G3 FAIL** | Otherwise | Pivot to calibration-only or stop |
| **G3 NOT EVALUABLE** | Fewer than 10 L2 SNVs across all week-6 genes | G3 is decided at the final evaluation instead, and this is said |

**Leave-family-out.** Trained without the U4 family. RNU4ATAC ClinVar P/LP (≥ 1 star) against ClinVar B/LB, AUROC with an interval. Reported, not gating.

## 9. Final evaluation (single unsealing)

After the model, mutation model and calibration code are frozen by a dated log entry:

1. unseal `sealed/blind_labels.tsv` once;
2. score RNU2-2 and RNU5B-1 with models trained without them;
3. report each gene's AUROC (L2 against L3) and ΔAUROC against CADD, and add both genes to the Aim 3 summary as a second, blind estimate.

Nothing is refitted afterwards. Any later change is post hoc.

## 10. Aim 4: calibration (weeks 9–11)

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

## 11. Secondary endpoints

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

## 12. Analysis hygiene

- Seeds are fixed (20261013 unless stated).
- Every number in the paper is regenerated by `snrna.py` from pinned inputs.
- Each command writes a `run_info.json` with input checksums, versions and the git commit.
- Results files carry the command that made them.

## 13. Resolutions of `PLAN.md` ambiguities (check before the affected analysis runs)

1. **G2 decision when only some targets are recovered.** `PLAN.md` says G2 fails if "the known RNU4-2 depletion and two blind hotspots are not recovered". Here: PASS needs all three; PARTIAL (continue) needs RNU4-2 plus one hotspot; KILL is anything less (section 7.5).
2. **The critical-region p-value.** `PLAN.md` asks for "p < 0.01" without naming the test. Here the primary is an exact Poisson-binomial test against the gene-median-scaled expectation. A circular-shift permutation with 145 positions can reach p < 0.01 only if C is the single most depleted placement, so it is reported but not gating.
3. **Sealed genes at G3.** `PLAN.md` lists RNU2-2 and RNU5B-1 among G3's held-out genes, but also says their labels open only once, at the final evaluation. Here G3 uses non-sealed genes only; the sealed genes are evaluated once, at the end (section 9).
4. **Site-rate covariate in G2.** It is excluded from G2 for the circularity reason given in section 7.2.
5. **"ClinGen default prior".** Taken as P = 0.10 (Tavtigian 2018); Pejaver's 0.0441 is a sensitivity value.
6. **Population controls in week 1.** The gnomAD-based control set is primary, as `PLAN.md` Appendix A says. The paper's 346-SNV set is used only to reproduce published numbers.
7. **Blind hotspot coordinates.** If guidance Table 2 (read directly) gives coordinates different from `PLAN.md`'s, Table 2 wins. This is recorded as a dated clarification before G2 is run. It is the only coordinate change allowed.

## 14. Deviation log

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

## 15. Proposed amendments (the authors decide before the affected analysis runs)

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
