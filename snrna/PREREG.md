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
