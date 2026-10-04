# Feasibility gate: pre-registered plan

**Project.** De novo excess in developmental-disorder (DD) trios as a ClinVar-free truth set for calibrating ACMG/AMP PP3/BP4 evidence for de novo missense variants.

**Status of this document.** Written, together with `gate.py`, before any calibration result was computed on real data. Commit it (with the scripts) to a timestamped place before running `python gate.py run` on the real data; see step 5 of README.md. After that point:

- only bug fixes are allowed;
- each fix gets a dated entry in the deviation log at the bottom;
- no threshold, bin, stratum, rate model or criterion below may change.

## 1. Question

Can the published de novo mutations (DNMs) from 31,058 DD trios support likelihood ratios for AlphaMissense that are precise enough to place it on the ClinGen PP3/BP4 scale? The DNMs come from DDD, GeneDx and Radboudumc (Kaplanis et al. 2020). Expected counts come from the Roulette mutation-rate model, and no ClinVar labels are used.

## 2. Data (all open and citable)

| Input | Version | Licence / source |
|---|---|---|
| De novo list | `DDD_RUMC_GDX_denovos_cadd_shet_wweights_minrep_indelfix_2025_09_18.txt.gz`, DeNovoWEST commit `49068922` (github.com/HurlesGroupSanger/DeNovoWEST, formerly queenjobo/DeNovoWEST) | Kaplanis et al. 2020 (Supplementary Table 1, re-normalised 2025-09-18 version); repository MIT |
| Gene table, s_het | `extended_denovoWEST_results.tab` (same commit) | repository MIT; s_het values from Cassa et al. 2017 (Supplementary Table 1) |
| DD gene list | `DDG2P_22_9_2019.csv` (same commit) | Gene2Phenotype |
| AlphaMissense | `AlphaMissense_hg38.tsv.gz`, Zenodo 10.5281/zenodo.10813168 (md5 `9fd167735f16a1b87da6eb3e4c25fcb5`) | CC BY 4.0 |
| Transcript models | GENCODE v32 (the release AlphaMissense hg38 was mapped with) | open |
| Mutation rates | Roulette per-site rates, files `N_rate_v5.2_TFBS_correction_all.vcf.bgz`, GRCh38 autosomes | Seplyarskiy et al. 2023; public download, no licence stated |
| Liftover | UCSC hg19ToHg38 chain, applied with pyliftover | open |

## 3. Universe and counts

The **universe** is every possible single-nucleotide variant (SNV) that meets all of the following:

- it lies in the coding sequence of an AlphaMissense canonical transcript on autosomes 1–22;
- it is either synonymous, or missense with an AlphaMissense score;
- the transcript passes QC: present in GENCODE v32, complete CDS, no internal stop, and at least 99% concordance between our codon enumeration and AlphaMissense in both directions;
- it does not fall in codon 1;
- it is not shared by more than one QC-passing transcript;
- its Roulette FILTER is not `low`.

**Observed (O).** O is the count of published SNV DNMs, lifted to GRCh38, that land on a universe site with the REF allele confirmed. X-chromosome DNMs are out of scope because Roulette has no X.

**Expected (E).** For a set of sites S:

- E(S) = c · Σ_{s∈S} μ_s
- μ is the raw Roulette rate (MR).
- c = O_syn / Σ μ over synonymous universe sites, using all centres for pooled analyses and each centre's own synonymous DNMs for per-centre analyses.

## 4. Estimands

They are computed over AlphaMissense bins with edges 0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .98, .99 and 1.

- **Marginal LR**, LR_b = [X_b / ΣX] / [E_b / ΣE] with X_b = O_b − E_b. Pathogenic means excess de novo missense; benign means background de novo missense.
- **Within-gene LR**, LR_b = X_b / (Σ_g O_g·E_gb/E_g − E_b). This is the benign distribution reweighted to the gene composition of the pathogenic set, so it measures variant-level information beyond gene identity.
- **Strata**:
  - all genes;
  - DDG2P monoallelic confirmed/probable genes (the clinical use case);
  - genes with s_het ≥ 0.15;
  - genes with s_het < 0.15;
  - each centre (GeneDx, DDD, RUMC).
- **Cumulative LR** for AlphaMissense ≥ 0.5, 0.8, 0.9, 0.95, 0.98 and 0.99.
- **Uncertainty**: Poisson bootstrap over probands (weights ~ Poisson(1)), 1,000 replicates, seed 20261004, with percentile intervals.
- **Evidence mapping**:
  - LR cut-offs from Pejaver et al. 2022: 2.406, 5.790, 33.53 and 1124, with prior 0.0441. The BP4 cut-offs are their reciprocals.
  - Threshold rule: the smallest score τ such that every bin at or above τ has a one-sided 95% lower bound ≥ the cut-off. BP4 mirrors this.
  - Systematic margin. Evidence bounds also allow for miscalibration of the expected counts:
    - the PP3 bound is computed with E_b inflated by δ, and the BP4 bound with E_b deflated by δ;
    - exactly: the PP3 bound is the 5th bootstrap percentile of LR_b with X_b replaced by O_b − (1+δ)·E_b, and the BP4 bound is the 95th percentile with X_b replaced by O_b − (1−δ)·E_b. ΣX and E_b/ΣE (marginal), and the within-gene denominator Σ_g O_g·E_gb/E_g − E_b, keep their unshifted values;
    - δ = max(0.05, the largest |O/E − 1| in the gated synonymous null groups).

## 5. Gate criteria (decided now)

| | Criterion | Rule |
|---|---|---|
| **G0** pipeline QC | AlphaMissense transcripts passing all checks ≥ 95%; VEP missense/synonymous DNMs mapped into the universe ≥ 85%; REF mismatch after liftover ≤ 1%. Definitions: the mapping rate's denominator is every autosomal SNV DNM whose source VEP consequence is missense or synonymous, including those that fail liftover, and its numerator counts those landing on a universe site as defined in §3 except for the Roulette FILTER/rate condition (a data exclusion, not a mapping failure; the rate after that exclusion is reported, informational). The REF-mismatch rate is over lifted DNMs whose GRCh38 position has a Roulette record, i.e. whose REF can be checked | FAIL → fix pipeline, no interpretation |
| **G1** null calibration (**kill**) | Synonymous O/E in each of three groupings: ~10 equal-mass Roulette mutability groups; ~10 equal-mass groups of the codon's mean AlphaMissense score; CpG transition / other transition / transversion | **KILL** if any group's 95% CI lies wholly outside 0.90–1.10, which means significant miscalibration of more than 10%. WARN if any point estimate is outside 0.90–1.10, the omnibus χ² p < 0.001, or (for ordered groupings) the log-linear trend p < 0.001. |
| **G2** precision | In the four bins with AlphaMissense ≥ 0.90, marginal LR hi95/lo95 ≤ 3, and the top-bin lower bound > 1 | FAIL → underpowered |
| **G3** centre consistency | Cochran's Q p ≥ 0.01 across GeneDx/DDD/RUMC in the top four bins; a p-value that cannot be computed (non-finite LR or no bootstrap variance) counts as not passing | WARN only |

Equal-mass groups for G1 are formed by cutting the cumulative synonymous mutational mass over a fine grid (log10 MR in steps of 0.05; codon mean AlphaMissense in steps of 0.01) into 10 parts. If one fine bin carries enough mass to skip a group, the groups are renumbered consecutively, so fewer than 10 groups can result.

**Pre-registered fallback for G1.**

1. If the mutation-class check fails under the single pooled scale factor, switch to a separate scale factor per mutation class: CpG transitions, other transitions and transversions, each fitted on synonymous de novos of that class.
2. Re-run the mutability and protein-context checks under the class-specific factors.
3. If neither check fails, the gate continues with class-specific scaling and G1 is recorded as WARN. Otherwise the verdict is KILL.

No other remedy is allowed after results are seen. A class miscalibration also shows up in the mutability groups, because CpG transitions dominate the most mutable groups. The fallback therefore applies whenever mutation class fails, whatever happens to the mutability groups.

**Negative control (WARN only).** Missense de novos with AlphaMissense < 0.10 in genes with s_het < 0.15 should show no excess. If the 95% CI of O/E lies wholly outside 0.90–1.10, G1 is downgraded to WARN, because the synonymous scaling may not transfer to missense.

**Overall verdict**:

- **GO**: G0, G1 and G2 PASS, and G3 PASS.
- **GO WITH CAUTION**: G1 or G3 is WARN.
- **UNDERPOWERED**: G2 FAIL.
- **KILL**: G1 FAIL.
- **FIX PIPELINE**: G0 FAIL.
- Precedence when several apply: FIX PIPELINE > KILL > UNDERPOWERED > GO WITH CAUTION > GO.

Informational only, not gating: synonymous O/E by s_het class, by DDG2P membership and at exon edges. The sensitivity analyses are also informational: rate models AR, MG and MC; excluding sites within 3 bp of splice junctions; including FILTER=low sites; the other scaling, i.e. class-specific when the primary is pooled and vice versa. Roulette writes AR only where its Pol III / TF-binding correction applies, so a missing AR is taken as AR = MR. MG and MC are used where present, and the fraction of universe sites with each rate is reported.

The within-gene LR is the key scientific readout, but it is **not** a gate criterion: a small within-gene LR is a finding, not a failure.

## 6. Commitments for the next phase (made before seeing any result)

- The gate fits no model, so nothing can be tuned to it.
- For the ML model in the full study, **DDD is the held-out test centre**. No model selection, feature choice or hyper-parameter will use DDD de novos. Development uses GeneDx, and RUMC is used for internal validation.
- The full study adds sliding local windows for thresholds, about 30 predictors, and a training-data ledger for each predictor. Predictors tuned on these cohorts, such as popEVE, are evaluated only on the held-out centre.

## Pre-freeze amendments (2026-10-04)

These were made before the pre-registration commit, after an independent review of the code against this plan and a check of every citation. No real-data computation had been run; the Roulette files were downloading but the pipeline had read nothing. The self-test passed before and after the changes.

1. **§2: s_het source corrected** from Weghorn et al. 2019 to Cassa et al. 2017. The DeNovoWEST `s_het` values equal Cassa's Supplementary Table 1 for all 15,251 shared genes, and Kaplanis 2020 cites Cassa. Other data-table details were added: repository rename, AlphaMissense md5, and Roulette naming and licence status.
2. **§5 G0: definitions made explicit, and the code changed to match.** The mapping-rate denominator now includes DNMs that fail liftover, which the code had left out. The REF-mismatch rate is over DNMs whose REF can be checked; the code had included unchecked non-coding positions.
3. **§5 G3:** a Cochran's Q p-value that cannot be computed now counts as WARN. The code had counted it as PASS.
4. **§5 G1 grouping:** equal-mass groups are renumbered consecutively. This fixes a crash when one fine bin is heavy, since real Roulette MR is quantised, and changes nothing otherwise.
5. **§5 sensitivity:** a missing Roulette AR is taken as AR = MR, and the coverage of each rate model is reported.
6. **§3:** "shared by more than one transcript" now reads "more than one QC-passing transcript", which is what the code does. **§4:** the systematic-margin formula is written out exactly as coded. **§5:** verdict precedence is stated.
7. **Engineering, no effect on any estimate:**
   - `sites` frees intermediate tables (chr1 peak RAM ~7.5 → ~4.2 GB);
   - `run_info.json` records tool versions and hashes of the Roulette extracts and region files;
   - `environment.lock.yml` (exact package versions) and `.gitattributes` (LF line endings) are added;
   - downloads use https, the renamed DeNovoWEST repository and an AlphaMissense md5 check;
   - Roulette extraction checks for truncated input.

## Deviation log

| Date | Change | Reason |
|---|---|---|
| 2026-10-05 | `get_roulette.sh`: the truncation check "last extracted position ≥ start of the last coding region" is replaced by a check that the downloaded VCF ends with the BGZF end-of-file block. The extraction command is unchanged. | Bug fix, no effect on any estimate. The old check gave a false alarm on chr5, whose last coding region (OR4F3, 5:181,367,287–181,368,222) has no Roulette records at all: `bcftools query -r` returns 0, and the flanking sequence and the previous region are fully present. Roulette omits some duplicated regions. The chr5 VCF was intact: its size equals the server's Content-Length and its EOF block is present. Found during processing; no calibration result had been computed. |
| 2026-10-05 | Engineering note, data acquisition only; `gate.py` and `get_roulette.sh` unchanged. The Roulette VCFs were stored on an external USB drive. A USB fault (Windows disk event 51 at 00:22 IST) corrupted the six files being written at that moment (chr7–12), although their sizes matched the server's Content-Length; `chkdsk` then repaired the file system (0 bad sectors). Every Roulette VCF and index was then checked block by block (header, decompression, CRC32 and length of every BGZF block, unbroken block chain, EOF block) with a helper script kept outside this repository. Damaged byte ranges were re-fetched from the server with HTTP Range requests and re-checked: chr7 163, chr8 143, chr9 139, chr10 295, chr11 147, chr12 318 ranges; 0 for the other 16 chromosomes. Extraction ran only on files that had passed this check. chr8 had been extracted before the fault was found. Re-extracting it from the verified file gave identical extract content and identical site tables. The genome-wide G0 check (`--qc-only`, no calibration result) was then run and passed: transcripts 96.69%, mapping 93.40%, REF mismatch 0.003%. | Disk fault during download, found by bcftools read errors on chr10/chr11. No effect on any estimate. Recorded before the calibration run. |
| 2026-10-05 | **Exploratory, recorded after the calibration run; no estimate changed.** The informational exon-edge check showed synonymous O/E 0.17 [0.13, 0.23] within 3 bp of a splice junction. Post-hoc explanation: the source de novo list contains no `splice_region_variant` records. VEP ranks splice_region above synonymous but below missense, so synonymous de novos at exon edges were dropped from the source list while edge missense de novos stayed. Consequence: the synonymous scale factor c is about 2.4% low and the primary LRs are conservative. The pre-registered "exclude sites near splice junctions" sensitivity gives a top bin of 10.61 marginal and 3.55 within-gene. For the full study: pre-register excluding exon-edge (±3 bp) sites from both the synonymous scaling and the missense universe. | Post-hoc explanation of an informational check, labelled exploratory; not a change to the gate. |
