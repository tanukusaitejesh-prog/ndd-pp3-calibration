# Feasibility gate result

**Verdict: GO WITH CAUTION.** G0, G2 and G3 pass. G1 is a WARN, not a kill.

## Pre-registration

- **Plan and code commit:** `5c60cf7b3576ff7661656642677c2a320acebaa6`, committed 2026-10-04 17:41:32 UTC.
- **Independent timestamp:** the manifest (`prereg/prereg_manifest.txt`) holds the SHA-256 of every frozen file. It is stamped with OpenTimestamps (`prereg_manifest.txt.ots`), and the proof is anchored in Bitcoin block **969884**, mined 2026-10-04 17:50:50 UTC; its merkle root matches the proof. A second anchor is in block 969890.
- **Calibration run:** `python gate.py run`, results commit `4d58ae0` at 2026-10-04 21:44 UTC, about 4 hours after the anchor.
- **Code identity:** `gate.py` is byte-identical to the pre-registered file (SHA-256 `53001a2c…`).
- **Changes since freezing:** recorded in GATE_PLAN.md's deviation log:
  - a bug fix to the Roulette truncation check;
  - an engineering note on download integrity, made before the run;
  - an exploratory note on the splice-edge deficit, made after the run.

## What this likelihood ratio measures

The LR here does not measure the same thing as ClinVar-based calibration, so read it before comparing the two.

- **This study.** "Pathogenic" is the de novo missense *excess* over the Roulette expectation. "Benign" is the *background* de novo missense expected from mutation alone. That background includes damaging variants in genes unrelated to NDD. The LR is the evidence a score adds **for a de novo missense variant found in an NDD proband**, which is the clinical question.
- **ClinVar-based calibration** (Pejaver 2022; Bergquist 2025). "Benign" is curated benign variants, mostly common polymorphisms. The LR answers "pathogenic versus a known-benign variant".
- **So a low LR here is not "AlphaMissense is bad".** A score can separate damaging from polymorphic variants well and still add little for a de novo variant in an affected child. In that setting many background variants are also damaging, and the gene already carries much of the information.

The stratum closest to how PP3 is applied in practice is **within-gene, in known NDD genes (DDG2P monoallelic)**, because the clinician already knows the gene. This report leads with that stratum.

## Headline

**With the gene held fixed, AlphaMissense reaches only PP3_Supporting, and only at ≥ 0.99.** Marginally, it reaches PP3_Moderate at ≥ 0.99, mostly through gene identity.

| Stratum | LR at AM 0.99–1.00 [95% CI] | PP3 | BP4 |
|---|---|---|---|
| **DDG2P monoallelic, within-gene** | **3.86 [3.72, 3.99]** | Supporting ≥ 0.99 | Supporting ≤ 0.60 |
| DDG2P monoallelic, marginal | 5.88 [5.66, 6.09] | Supporting ≥ 0.99 | Supporting ≤ 0.60, Moderate ≤ 0.30 |
| All genes, within-gene | 3.45 [3.29, 3.61] | Supporting ≥ 0.99 | not reached |
| All genes, marginal | 9.68 [8.86, 10.57] | Supporting ≥ 0.98, Moderate ≥ 0.99 | not reached |

- **The rule:** these are Pejaver et al. 2022's LR cut-offs (prior 0.0441). A level is reached at the smallest score where every bin at or above it has a one-sided 95% bound past the cut-off.
- **The δ margin:** bounds include the pre-registered systematic margin δ = 0.091 on expected counts.
- **Resolution:** thresholds are at bin resolution, with edges 0.90/0.95/0.98/0.99; the full study uses sliding local windows.
- **BP4:** BP4 in DDG2P genes depends on the accuracy of the expected counts in low-score bins, so it is reported with δ applied.

## Comparison with ClinVar-based calibration

The published values below are from Bergquist et al., *Genet Med* 27:101402 (2025), doi:10.1016/j.gim.2025.101402, Table 1, read from the paper itself.

| Evidence level | ClinVar-based (Bergquist 2025) | De novo NDD, all genes, marginal | De novo NDD, DDG2P, within-gene |
|---|---|---|---|
| PP3 Supporting | 0.792–0.905 | ≥ 0.98 | ≥ 0.99 |
| PP3 Moderate | 0.906–0.971 | ≥ 0.99 | not reached |
| PP3 +3 points | 0.972–0.989 | not assessed (the gate uses the 4 Pejaver levels) | not assessed |
| PP3 Strong | ≥ 0.990 | not reached | not reached |
| BP4 Supporting | 0.100–0.169 | not reached | ≤ 0.60 |
| BP4 Moderate | 0.071–0.099 | not reached | not reached (marginal: ≤ 0.30) |

**Reading it:**

- **At AlphaMissense ≥ 0.99**, ClinVar calibration assigns Strong. The de novo data supports Moderate marginally and Supporting within gene.
- **Between 0.79 and 0.98**, ClinVar calibration gives Supporting or Moderate. Here the within-gene LR is about 0.8–1.4, so essentially no variant-level evidence.
- **What the comparison does and doesn't show:** the two calibrations answer different questions (see above). What it does support is that applying ClinVar-derived AlphaMissense thresholds to *de novo missense variants in NDD probands* would overstate the variant-level evidence. It does not show that the ClinVar calibration is wrong for its own question.
- **This is from a feasibility gate,** at bin resolution, with one predictor. The confirmatory comparison is part of the full study.

## Gate criteria

| | Result | Detail |
|---|---|---|
| G0 pipeline QC | PASS | Transcripts pass 96.69% (17,788/18,397; threshold ≥ 95%). De novos mapped 93.40% (≥ 85%). REF mismatch 0.003% (≤ 1%). |
| G1 null calibration | WARN | See below. |
| G2 precision | PASS | Marginal hi/lo in the top 4 bins: 1.50, 1.36, 1.37, 1.19 (≤ 3). The top-bin lower bound is 8.86. |
| G3 centre consistency | PASS | Cochran's Q p in the top 4 bins: 0.92, 0.77, 0.26, 0.16. |

**G1 in detail.** Every group's point estimate lies inside 0.90–1.10, and no 95% CI lies wholly outside it, so there is no KILL. The WARN comes from p-values below 0.001:

| Grouping | Omnibus χ² p | Trend p | Result |
|---|---|---|---|
| Mutability groups | 1.4e-4 | 1.2e-5 (last/first ratio 0.86) | WARN |
| Mutation class | 1.0e-8 | — | WARN (no trend test, since the groups are not ordered) |
| AlphaMissense codon context | 0.13 | 0.038 | PASS |

- **Source:** the same effect drives both: CpG transitions are over-predicted (O/E 0.919 [0.892, 0.944]), and the three most mutable groups sit at 0.91–0.93.
- **The fallback was not triggered.** Mutation class did not FAIL, so the pre-registered per-class fallback stayed off and pooled scaling remains primary.
- **Effect on the result:** per-class scaling (a pre-registered sensitivity analysis) moves the top-bin LRs by ≤ 1% (9.58 vs 9.68 marginal; 3.41 vs 3.45 within-gene).
- **δ = 0.091** equals the largest |O/E − 1| in the gated groups, as the rule requires.

## Exploratory finding (after the run): splice-edge synonymous deficit

- **What:** synonymous O/E is 0.17 [0.13, 0.23] within 3 bp of a splice junction. This is an informational check, not a gate.
- **Explanation:** the published de novo list contains **no `splice_region_variant` records**. VEP ranks splice_region above synonymous but below missense. So synonymous de novos at exon edges were dropped when the source list was built, while edge missense de novos stayed.
- **Status:** this explanation is post hoc, so it is labelled exploratory. It is logged in GATE_PLAN.md.
- **Consequence:** the synonymous scale factor *c* is about 2.4% low. Expected counts are therefore slightly low, which makes the primary LRs **conservative**.
- **Sensitivity check:** the pre-registered sensitivity analysis that excludes splice-adjacent sites gives a top bin of 10.61 marginal and 3.55 within-gene.
- **Full study:** pre-register excluding exon-edge (±3 bp) sites from both the synonymous scaling and the missense universe.
- **Negative control:** O/E 1.096 [1.057, 1.138], which passes. It is consistent with this under-scaling plus a small genuine low-score excess.

## Other results

- **Missense excess:** O/E 1.284 [1.256, 1.318], an excess of 5,454 [4,982, 5,968]. About 22% of de novo missense are pathogenic.
- **Robustness:** across rate models AR, MG and MC, and when FILTER=low sites are included, the top-bin LRs stay at 8.0–10.0 marginal and 3.2–3.5 within-gene.

## Data-integrity incident (deviation log, before the run)

- A USB fault corrupted the Roulette VCFs for chr7–12 while they were downloading.
- Every VCF was then checked block by block (BGZF CRC32), and only the damaged ranges were re-fetched: 139–318 per affected chromosome, 0 elsewhere.
- chr8 had been extracted before the fault was found. Re-extracting it from the verified file gave identical output.
- No pipeline code changed.

## Files

- **Report:** `results/gate_report.txt`
- **Tables:** `lr_bins.tsv`, `null_checks.tsv`, `thresholds.tsv`, `cumulative_lr.tsv`, `qc.tsv`, `qc_by_chrom.tsv`, `run_info.json`
- **Figures:** `lr_curves.png`, `null_checks.png`
- **Repository:** `/root/ndd_pp3_gate` in WSL (distro `ndd-gate`)
- **Full bundle:** `prereg/ndd-pp3-gate-with-results.bundle`

## Next (GATE_PLAN.md §6, plus items from this gate)

The full study, with these design points:

- **DDD held out** as the test centre; GeneDx for development, RUMC for internal validation;
- sliding local-window thresholds;
- about 30 predictors;
- a side-by-side comparison with ClinVar calibration;
- pre-registered exclusion of exon-edge sites;
- results led by the DDG2P within-gene stratum.
