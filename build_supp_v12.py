"""Builds a clean SUPPLEMENT_v12.md (no console output; tables generated from committed result files). Run: python build_supp_v12.py"""
import re
import numpy as np, pandas as pd

R = "results_fullstudy_backup/"


def ci(x):
    a, b = [float(v) for v in re.findall(r"-?\d+\.\d+", str(x))[:2]]
    return f"[{a:.2f}, {b:.2f}]"


# Supplementary Table 5 rows: de novo dependence, gene-clustered
g = pd.read_csv(R + "a12_gclust2.tsv", sep="\t")
g["cohort"] = g["cohort"].replace({"public3": "public"})
t5 = ["| Cohort | Predictor | PP3 tier | Flag | R | f | Within gene, all | Within gene, DDG2P |", "|---|---|---|---|---|---|---|---|"]
for r in g.itertuples():
    t5.append(f"| {r.cohort} | {r.tool} | {r.tier} | {r.flag} | {r.R:.2f} {ci(r.R_ci)} | {r.f:.2f} {ci(r.f_ci)} | "
              f"{r.R_within_all:.2f} {ci(r.R_within_all_ci)} | {r.R_within_ddg2p:.2f} {ci(r.R_within_ddg2p_ci)} |")
# assays
m = pd.read_csv(R + "mave_all_sites_summary.tsv", sep="\t").dropna(subset=["R_sup"]).sort_values(["ddg2p", "gene"])
ta = ["| Assay | Gene set | Variants | PM1-like share | LR of PM1-like | R, Supporting [95% CI] | R, Strong |", "|---|---|---|---|---|---|---|"]
for r in m.itertuples():
    ta.append(f"| {r.ds} | {'DDG2P' if r.ddg2p else 'other'} | {int(r.n)} | {r.F_frac:.2f} | {r.LR_F:.2f} | {r.R_sup:.2f} [{r.lo_sup:.2f}, {r.hi_sup:.2f}] | {r.R_str:.2f} |")
# benchmark (bench2) top rows
b = pd.read_csv(R + "bench2.tsv", sep="\t")
b = b[(b.stratum == "all_genes") & (b.site_set == "a1") & (b.q == 0.02) & (b.predictor != "NULL") & (b.LR_marginal >= 2.406)].copy()
b["rank_m"] = b.LR_marginal.rank(ascending=False).astype(int)
b = b.sort_values("LR_window", ascending=False)
tb = ["| Predictor | Pooled LR (rank) | Within-window LR [95% CI] | Rank [95% interval] | Retention |", "|---|---|---|---|---|"]
for r in b.itertuples():
    tb.append(f"| {r.predictor} | {r.LR_marginal:.2f} ({r.rank_m}) | {r.LR_window:.2f} [{r.LR_window_lo:.2f}, {r.LR_window_hi:.2f}] | "
              f"{int(r.rank_window)} [{int(round(r.rank_window_lo))}, {int(round(r.rank_window_hi))}] | {r.r_window:.2f} |")
T5, TA, TB = "\n".join(t5), "\n".join(ta), "\n".join(tb)

S = f"""# Supplementary material — Dependence between computational and location evidence in missense variant classification, and its consequences in expert and laboratory practice

*Draft 12, 7 October 2026. Accompanies `MANUSCRIPT_v12.md`. Full outputs for every analysis are in `results_fullstudy/` in the study repository; file names are given in each note.*

---

## Supplementary Note 1. Analysis log and pre-specified tests

The analysis plan (`FULL_STUDY_PLAN.md`) and its seventeen amendments (A1–A17) were committed to the study repository, each before the analysis it governs. Each amendment records what had been seen when it was written. Changes made after each amendment are in the deviation log at the end of the plan.

None of these documents was deposited in an independent registry before its analysis.

**Pre-specified readout.** The share of informative predictors (marginal top-2% LR ≥ 2.406) whose evidence within 50-residue windows is less than half their marginal evidence (retention r < 0.5). The criterion was a share of ≥ 75%.

**Supplementary Table S1a. Pre-specified retention tests** (exact binomial intervals; these ignore correlation between predictors).

| Cohort and predictors | Position control | Informative | r < 0.5 | 95% CI |
|---|---|---|---|---|
| DD, dbNSFP 4.1a (primary) | D1 | 24 | 15 (62%): not met | 41–81% |
| DD, dbNSFP 4.1a | D2 | 24 | 19 (79%) | 58–93% |
| Autism, dbNSFP 4.1a (pre-specified D2) | D2 | 20 | 16 (80%) | 56–94% |
| DD, 17 newer predictors (pre-specified D2) | D2 | 14 | 13 (93%) | 66–100% |
| Public, dbNSFP 4.1a (descriptive) | D2 | 15 | 11 (73%) | 45–92% |

**Why the primary test failed.**
- **The original position control (D1)** split genes into as many evenly spaced chunks as they have regional-constraint regions. In genes without a constraint split, that is a single unit, so D1 does not control position.
- **The fixed-window control (D2)** was adopted after this result and pre-specified for all later tests.

Across all 39 predictors, with genes and predictors both resampled, the share with r < 0.5 was 0.87 [0.69, 0.97].

**What followed.**
- **Pre-specified before its results (amendment A10):** the PP3 × PM1 dependence test.
- **Not pre-specified, and descriptive:** the PM1 × PP2 and triple tests, the benchmark, the oracle, and the expert-panel and laboratory analyses.

---

## Supplementary Note 2. Estimator, controls and robustness

**Estimator.**
- **Marginal LR:** (X_b/ΣX)/(E_b/ΣE).
- **Within-unit LR:** X_b / (Σ_u O_u E_ub/E_u − E_b).

**Simulation.** One hundred datasets per scenario at the real sample size. Planted within-gene signals were recovered with relative bias ≤ 1%, and 95% interval coverage was 0.86–0.99. Because coverage dipped below 0.95, every R and f interval in the paper resamples genes.

**Controls** (24 dbNSFP 4.1a predictors in DD, median retention):
- **Random partitions.** A random partition of each gene into as many groups as it has windows kept 0.51, against 0.56 for gene alone and 0.31 for real windows. The paired difference was 0.21 [0.19, 0.23], positive for all 24 predictors.
- **Window size.** Real windows against equal-sized random partitions:
  - 25 residues: 0.27 vs 0.48;
  - 50 residues: 0.33 vs 0.52;
  - 100 residues: 0.39 vs 0.53.
- **InterPro domains** kept 0.33.
- **Unaffected siblings:** no predictor showed evidence (0 of 37).
- **Perturbations.** Retention was unchanged by ±9% CpG mutation-rate perturbation and by removing all 64,573 ClinVar P/LP sites (0.33–0.34).
- **Thresholds.** It was similar at the 1%, 2% and 10% thresholds (0.33, 0.32, 0.29).

**Null predictor** (a hash of the variant key): LR 0.84–1.11 at every threshold and site set.

---

## Supplementary Note 3. Predictor evidence incremental to location

**Retention.** Across 39 informative predictors in DD, the median top-2% LR was 6.06 marginally, 3.00 within gene and 1.79 within 50-residue windows. The median retention was 0.32 (`bench2.tsv`).

**The oracle** (`oracle_summary.txt`) is an analytic perfect predictor under the same estimators. It assumes pathogenic sites carry de novo rate k × mutation rate, with each window's excess taken from a fitted rate model.

**Supplementary Table S3a. Perfect-predictor evidence** (top 2% of opportunity, DD).

| k | Pathogenic share of opportunity | Marginal LR | Within gene | Within window | Retention within window |
|---|---|---|---|---|---|
| 10 | 3.5% | 28.4 | 5.7 | 2.8 | 0.30 |
| 30 | 1.1% | 50.0 | 15.6 | 8.2 | 0.54 |
| 100 | 0.3% | 50.0 | 30.4 | 20.1 | 0.77 |
| 300 | 0.1% | 50.0 | 41.2 | 33.5 | 0.90 |

k is not identified by the data. Retention therefore measures evidence incremental to location; it is not a measure of predictor skill.

**Benchmark.** Predictors were ranked by within-window LR (`bench2.tsv`).
- **Across all genes, the pooled and within-window rankings were not positively associated:** Spearman −0.17 [−0.26, 0.06] at 2%, and between −0.02 and −0.22 across thresholds and site sets.
- **Within DDG2P genes they agreed** (0.68 [0.53, 0.76]).
- **Validation.** Within-window LR correlated with within-gene ClinVar AUROC (ρ 0.55 [0.22, 0.77]) and with assay performance (0.63 [0.39, 0.79]). Marginal LR also correlated with both (0.37), and the differences were not significant (`benchval.tsv`).

**Supplementary Table S3b. Gene- and position-adjusted benchmark** (DD, top 2%; gene-clustered intervals).

{TB}

---

## Supplementary Note 4. Dependence on the de novo scale

**Supplementary Table 5. All dependence estimates** (`a12_gclust2.tsv`; gene-clustered 95% intervals). f is the fraction of summed log-evidence kept.

{T5}

**Other analyses** (AlphaMissense PP3 at Supporting × PM1-like, DD unless stated):

| Analysis | R [95% CI] | Source |
|---|---|---|
| All 40 informative predictors | median 0.51 (range 0.18–0.80); all upper bounds < 1 | `a12_allr_summary.txt` |
| Excluding the 50 genes with the largest excess | 0.48 [0.42, 0.53] | `a12_lto_summary.txt` |
| AlphaFold pLDDT ≥ 70 / < 70 | 0.53 [0.50, 0.57] / 0.39 [0.28, 0.51] | `a12_plddt_summary.txt` |
| Buried / exposed sites | 0.55 [0.52, 0.59] / 0.36 [0.29, 0.43] | `a12_burial_summary.txt` |
| PM1-like × PP2-like, DD / autism | 0.31 [0.28, 0.33] / 0.36 [0.30, 0.45] | `pm1pp2_summary.txt` |
| PP3 × PM1-like × PP2-like (R3), DD / autism | 0.13 [0.11, 0.14] / 0.18 [0.13, 0.24] | `triple_summary.txt` |

**Mechanism.** AlphaMissense flags 51% of mutation-expected sites at PM1-like positions, against 19% overall. Among excess de novos the rise is only from 79% to 94%.

---

## Supplementary Note 5. ClinVar and functional assays

**ClinVar** (P/LP against B/LB, DD de novo sites removed): AlphaMissense PP3 × PM1-like gave R 0.32 [0.28, 0.37] at Supporting and 0.23 [0.18, 0.28] at Strong in DDG2P genes. It gave 0.20 and 0.15 across all genes. ClinVar classifications use PP3 and PM1 themselves, so these estimates are partly circular.

**Functional assays** (`make_mave_all.py`; `mave_all_summary.txt`).
- **Data.** Human ProteinGym v1.3 single-substitution assays, with UniProt identifiers mapped to gene symbols through the UniProt REST API.
- **Scope.** 82 assays in 70 genes were mapped. Sixteen were informative (PM1-like share 5–95%, ≥ 10 damaging and ≥ 10 non-damaging variants).
- **Pooled R** (random effects):
  - outside DDG2P genes (13 assays): 0.94 [0.88, 0.99] at Supporting and 0.98 [0.92, 1.03] at Strong;
  - all 16 assays: 0.89 [0.83, 0.95].
- **The PM1-like flag** barely predicted damage within these assays: median LR 1.10 outside DDG2P genes (IQR 1.04–1.42).

**Supplementary Table S5. Informative assays.**

{TA}

---

## Supplementary Note 6. Expert-panel practice

**Data.** ClinGen Evidence Repository, export of 2026-10-07 (`make_stacking_audit.py`, `make_stacking_cv.py`).

**Co-application.** Of 6,177 missense classifications:
- PP3 was applied in 2,990, PM1 in 1,092 and PP2 in 592;
- PP3 + PM1 in 787 (Supporting 605, Moderate 178, Strong 4);
- PM1 + PP2 in 174: RASopathy 62, FBN1 61, Monogenic Diabetes 22, Brain Malformations 11, PTEN 11, Pulmonary Hypertension 5, ACADVL 2;
- all three in 138.

**Truth.** Points after removing PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5. P/LP ≥ 6, B/LB ≤ −1: 2,942 classifications in 130 genes.
- **Dependence.** R(panel PP3, panel PM1) was 0.21 [0.12, 0.31].
- **PM1 × PP2 and PP3 × PP2** were not estimable (1 and 2 B/LB in the joint cells).

**Supplementary Table S6a. PP3 credit beside panel PM1** (points [95% CI]; gene-clustered).

| Truth definition | Any strength | Supporting | Moderate |
|---|---|---|---|
| Primary | 0.60 [0.14, 1.05] (304; 17 B/LB) | 0.81 [0.24, 1.68] (235; 11) | 0.02 [−1.69, 0.89] (67; 6) |
| P/LP ≥ 4 | 0.61 [0.15, 1.07] | 0.82 [0.28, 1.77] | 0.05 [−1.13, 0.80] |
| P/LP ≥ 8 | 0.61 [0.13, 1.04] | 0.89 [0.30, 1.84] | −0.23 [−2.61, 0.59] |
| B/LB ≤ −3 | 0.65 [0.16, 1.30] | 0.83 [0.30, 2.06] | 0.14 [−2.08, 2.34] |
| PS2, PM6 also removed | 0.58 [0.08, 1.04] | 0.75 [0.17, 1.56] | 0.15 [−2.10, 0.97] |
| PS1, PP5, BP6, PS2, PM6 also removed | 0.57 [0.08, 1.04] | 0.73 [0.21, 1.60] | 0.17 [−2.06, 1.03] |
| Frequency criteria also removed | 0.49 [−0.06, 1.31] (6 B/LB) | — | — |

**Moderate against Supporting** (primary truth): −0.79 [−2.81, 0.27], bootstrap p = 0.14. The standard error was 1.55, so the minimum detectable difference at 80% power is 4.35 points.

**Out of sample.**
- **Gene splits.** Across 200 random 50/50 gene splits, the credit was 0.61 on the estimating half and 0.59 on the held-out half (medians). The held-out credit fell below the points applied in 98% of splits.
- **Temporal split.** Estimated on classifications approved before 2025, the credit was 0.70 [−0.14, 1.26]. On classifications approved in 2025–2026 it was 0.40 [−0.25, 1.07], against 1.40 applied. Among the 2,507 later classifications, limiting PP3 to Supporting changed 15 and the cap 16.

---

## Supplementary Note 7. PP3 defined by score rather than by panel choice

**Data and truth.** Panel classifications were mapped to analysed sites (`make_erepo_eligibility.py`). Truth excludes PS2 and PM6 in addition: 2,381 classifications, 103 genes.

**Selection into PP3.** Within panel-PM1 classifications, likelihood ratios for pathogenicity were:
- 2.16 where PP3 was applied (n 254);
- 1.05 where the variant was eligible by score but PP3was not applied (n 70);
- 0.04 where it was ineligible (n 21).

Panels apply PP3 preferentially to more convincing variants.

**Supplementary Table S7. Points PP3 earns beside panel PM1, by definition.**

| Definition | Points [95% CI] | n | B/LB |
|---|---|---|---|
| Panel-applied PP3 | 0.88 [0.35, 1.49] | 254 | 13 |
| REVEL ≥ 0.644 or AlphaMissense ≥ 0.792 | 0.71 [0.35, 1.20] | 322 | 19 |
| REVEL ≥ 0.644 | 0.68 [0.34, 1.19] | 315 | 19 |
| AlphaMissense ≥ 0.792 | 1.66 [0.98, 2.58] | 264 | 7 |
| REVEL band: Supporting / Moderate / Strong | −0.52 [−1.46, 1.17] / 0.58 [0.08, 1.46] / 1.84 [0.81, —] | 45 / 137 / 133 | 7 / 9 / 3 |
| AlphaMissense band: Supporting / Moderate / Strong | 1.63 [0.09, —] / 1.44 [0.34, —] / 1.75 [0.76, 3.74] | 37 / 63 / 164 | 1 / 2 / 4 |

"—" marks an unbounded interval (≤ 3 B/LB). The definitions disagree. AlphaMissense-defined PP3 shows no deficit, whereas REVEL-defined and panel-applied PP3 do.

---

## Supplementary Note 8. Laboratory practice

**Data.** ClinVar 2026-10 germline submissions for missense SNVs, excluding expert-panel and practice-guideline submissions, whose comments list any ACMG/AMP code: 87,873 submissions from 472 submitters (`make_lab_stacking.py`).

**Parser validation** (`make_parser_validation.py`). The free-text summaries of 13,091 Evidence Repository classifications were parsed and compared with their structured codes.

| Code | Precision | Recall | Strength agreement |
|---|---|---|---|
| All codes | 0.92 | 0.98 | — |
| PP3 | 0.95 | 0.99 | 0.99 |
| PM1 | 0.98 | 0.99 | 0.99 |
| PP2 | 0.98 | 0.99 | 1.00 |

Most remaining false positives are codes that the text lists as not met.

**Co-application.**
- PP3 + PM1 in 4,471 submissions (Supporting 4,207, Moderate 184, Strong 80);
- PM1 + PP2 in 1,481;
- all three in 1,035.

**Truth.** 39,015 submissions (P/LP 9,455, B/LB 29,560; 5,976 genes, 287 submitters).

**Supplementary Table S8. Laboratory truth.**

| Estimate | Value | Gene-clustered | Submitter-clustered | PS1, PP5, BP6, PS2, PM6 also removed |
|---|---|---|---|---|
| R(PP3, PM1) | 0.18 | [0.14, 0.24] | [0.08, 0.43] | 0.18 [0.14, 0.24] |
| PP3 Supporting credit beside PM1 | 0.84 | [0.52, 1.24] | [0.06, 1.61] | 0.82 [0.52, 1.20] |
| R(PM1, PP2) | 0.27 | [0.19, 0.46] | [0.14, 0.54] | 0.30 [0.21, 0.48] |
| PP2 credit beside PM1 | −0.22 | [−0.67, 0.42] | [−1.23, 0.79] | −0.41 [−0.89, 0.24] |
| R(PP3, PP2) | 0.42 | [0.35, 0.51] | — | — |

**Thin cells.** PP3 credit at Moderate (3 B/LB) and Strong (0 B/LB) is not estimable.

**Re-scoring** (changed classifications, of 87,873):

| Rule | Changed | P→LP | LP→VUS | Other |
|---|---|---|---|---|
| ClinGen cap | 40 | 30 | 10 | 0 |
| PP3 + PM1 + PP2 ≤ 4 | 63 | 39 | 22 | 2 |
| PP2 not added when PM1 applies | 352 | 124 | 223 | 5 |
| PP3 limited to Supporting when PM1 applies | 90 | 39 | 51 | 0 |

**An earlier version** restricted the data to submissions that mention PP3 or PM1. That selects on the exposure and biased the marginal LRs, so it is not used.

---

## Supplementary Note 9. Consequences in developmental disorders and ClinVar

**Developmental disorders.** In DD, uncapped co-use of PS2 + PM2 + PM1-like + AlphaMissense PP3 placed 2,438 de novo variants at Pathogenic.
- The de novo excess implies that 84% of that class is disease-associated.
- Taken literally, that is below even the 90% floor for Likely Pathogenic.
- The estimate is class-level, and refers to association with any severe developmental disorder rather than a per-variant posterior.

**ClinVar.** In DDG2P monoallelic genes, 37.3% of AlphaMissense PP3-positive VUS lie at PM1-like sites, against 13.2% across all genes. Without the cap, PP3 + PM1 at Moderate alone would reach Likely Pathogenic for 7,101 of them.

---

## Supplementary Note 10. Attempts to build better predictors or rules from de novo data (negative)

- **De novo-trained classifiers (A7, A8, A14, A15).** These were gradient-boosted models separating real de novos from mutation-matched sites. On ClinVar within DDG2P genes their within-gene AUROC was 0.923 (v3), 0.905 (v4) and 0.879 (v5, label-free inputs only), against 0.967 for AlphaMissense.
- **Hierarchical gene → window → AlphaMissense model (A16).** Specification: `prereg_fullstudy/A16_HIERARCHICAL_PREREG.md`. On ClinVar variants reclassified from VUS between 2023 and 2026, it calibrated worse than the cap: log-loss difference +0.147 [+0.106, +0.190]. Its gene-specific PP3 evidence did not track ClinVar gene-specific calibration in 142 DDG2P genes (Spearman 0.01).
- **Within-window model (A17, task 13).** A Poisson model learning only within-window variation added nothing to AlphaMissense on held-out de novos: +0.0008 [−0.0006, +0.0023] nats per de novo. It was worse on ClinVar (AUROC 0.894 against 0.967) and in 7 of 8 assays.

---

## Supplementary Tables (files)

| Table | Content | File |
|---|---|---|
| 1 | Per-predictor retention with intervals | `supp_table1_predictor_resource.tsv` |
| 2 | Training-data ledger for predictors | `aim3_training_ledger.tsv` |
| 3 | ClinVar sensitivity of the rule comparison | `a10_clinvar_sensitivity.tsv` |
| 4 | Conditional PP3 points on the de novo scale (scale-dependent; see Supplementary Note 4) | `supp_table4_conditional_pp3.tsv` |
| 5 | All dependence estimates (Supplementary Note 4) | `a12_gclust2.tsv` |
| 6 | Hierarchical-model calibration table | `a16_calibration_table.tsv` |
| 7 | Benchmark, all predictors, thresholds and site sets | `bench2.tsv` |
"""
S = S.replace("PP3was", "PP3 was")
open("SUPPLEMENT_v12.md", "w", encoding="utf-8").write(S)
print(len(S.split()), "words")
