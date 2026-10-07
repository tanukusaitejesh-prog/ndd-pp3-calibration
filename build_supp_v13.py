"""Builds SUPPLEMENT_v13.md: seven notes, unique table numbers (S1-S12), numbers identical to the manuscript.
Run: python build_supp_v13.py"""
import re
import pandas as pd

R = "results_fullstudy_backup/"


def ci(x):
    a, b = [float(v) for v in re.findall(r"-?\d+\.\d+", str(x))[:2]]
    return f"[{a:.2f}, {b:.2f}]"


g = pd.read_csv(R + "a12_gclust2.tsv", sep="\t")
g["cohort"] = g["cohort"].replace({"public3": "public"})
t_dep = ["| Cohort | Predictor | PP3 tier | Flag | R | f | Within gene, all | Within gene, DDG2P |", "|---|---|---|---|---|---|---|---|"]
for r in g.itertuples():
    t_dep.append(f"| {r.cohort} | {r.tool} | {r.tier} | {r.flag} | {r.R:.2f} {ci(r.R_ci)} | {r.f:.2f} {ci(r.f_ci)} | "
                 f"{r.R_within_all:.2f} {ci(r.R_within_all_ci)} | {r.R_within_ddg2p:.2f} {ci(r.R_within_ddg2p_ci)} |")
m = pd.read_csv(R + "mave_all_sites_summary.tsv", sep="\t").dropna(subset=["R_sup"]).sort_values(["ddg2p", "gene"])
t_as = ["| Assay | Gene set | Variants | PM1-like share | LR of PM1-like | R, Supporting [95% CI] | R, Strong |", "|---|---|---|---|---|---|---|"]
for r in m.itertuples():
    t_as.append(f"| {r.ds} | {'DDG2P' if r.ddg2p else 'other'} | {int(r.n)} | {r.F_frac:.2f} | {r.LR_F:.2f} | "
                f"{r.R_sup:.2f} [{r.lo_sup:.2f}, {r.hi_sup:.2f}] | {r.R_str:.2f} |")
TDEP, TAS = "\n".join(t_dep), "\n".join(t_as)

S = f"""# Supplementary material: Dependence between computational and location evidence in missense variant classification: a label-free measurement and an audit of practice

*Draft 13, 7 October 2026. Accompanies `MANUSCRIPT_v13.md`. Complete outputs are in `results_fullstudy/` in the study repository; each note names its files.*

---

## Supplementary Note 1. Analysis log and pre-specified tests

**The log.** The analysis plan (`FULL_STUDY_PLAN.md`) and seventeen amendments were committed to the study repository, each before the analysis it governs, with a deviation log. They were not deposited in an independent registry.

**Other analyses in the log, not reported in this paper:**
- a predictor benchmark;
- de novo-trained classifiers and a hierarchical rule (both negative).

They are in the repository.

**Primary pre-specified readout.** The share of informative predictors (marginal top-2% likelihood ratio ≥ 2.406) whose evidence within 50-residue windows is less than half their marginal evidence, with ≥ 75% as the criterion.

**Table S1. Pre-specified retention tests** (exact binomial intervals; these ignore correlation between predictors).

| Cohort and predictors | Position control | Informative | Retention < 0.5 | 95% CI |
|---|---|---|---|---|
| DD, dbNSFP 4.1a (primary) | D1 | 24 | 15 (62%): not met | 41–81% |
| DD, dbNSFP 4.1a | D2 | 24 | 19 (79%) | 58–93% |
| Autism, dbNSFP 4.1a (pre-specified D2) | D2 | 20 | 16 (80%) | 56–94% |
| DD, 17 newer predictors (pre-specified D2) | D2 | 14 | 13 (93%) | 66–100% |
| Public, dbNSFP 4.1a (descriptive) | D2 | 15 | 11 (73%) | 45–92% |

**Why the primary test failed.** D1 divides each gene into as many chunks as it has regional-constraint regions, which in genes without a constraint split is a single unit. The fixed-window control (D2) was adopted after this result.

---

## Supplementary Note 2. Estimator, controls and evidence incremental to location

**Estimators.**
- Marginal likelihood ratio: (X_b/ΣX)/(E_b/ΣE).
- Within-unit likelihood ratio: X_b/(Σ_u O_u E_ub/E_u − E_b).

**Simulation.** One hundred datasets per scenario at the real sample size. Planted within-gene signals were recovered with relative bias ≤ 1%, and 95% interval coverage was 0.86–0.99. Because coverage fell below 0.95, every R and f interval resamples genes.

**Controls** (24 dbNSFP 4.1a predictors in DD, median retention):
- **Random partitions.** A random partition into as many groups as windows kept 0.51, against 0.56 for gene alone and 0.31 for real windows. The paired difference was 0.21 [0.19, 0.23].
- **Window size.** Windows of 25, 50 and 100 residues kept 0.27, 0.33 and 0.39, against 0.48, 0.52 and 0.53 for equal-sized random partitions.
- **InterPro domains** kept 0.33.
- **Siblings.** No predictor showed evidence in siblings (0 of 37).
- **Perturbations.** Retention was unchanged by ±9% CpG mutation-rate perturbation and by removing all ClinVar P/LP sites.
- **Null predictor** (a hash of the variant key): likelihood ratio 0.84–1.11.

**Evidence incremental to location.** Across 39 informative predictors in DD, the median top-2% likelihood ratio was 6.06 marginally, 3.00 within gene and 1.79 within 50-residue windows. The median retention was 0.32.

**Table S2. Evidence available to a perfect predictor** (analytic; pathogenic sites carry de novo rate k × mutation rate; top 2% of opportunity; DD).

| k | Pathogenic share of opportunity | Marginal LR | Within gene | Within window | Retention within window |
|---|---|---|---|---|---|
| 10 | 3.5% | 28.4 | 5.7 | 2.8 | 0.30 |
| 30 | 1.1% | 50.0 | 15.6 | 8.2 | 0.54 |
| 100 | 0.3% | 50.0 | 30.4 | 20.1 | 0.77 |
| 300 | 0.1% | 50.0 | 41.2 | 33.5 | 0.90 |

k is not identified by the data. Retention therefore measures evidence incremental to location, not predictor skill.

---

## Supplementary Note 3. Dependence on the de novo scale

**Table S3. All dependence estimates** (`a12_gclust2.tsv`; gene-clustered 95% intervals; f is the fraction of summed log-evidence kept).

{TDEP}

**Table S4. Further de novo analyses** (AlphaMissense PP3 at Supporting × PM1-like, DD unless stated).

| Analysis | Result | File |
|---|---|---|
| All 40 informative predictors (all analysed sites) | median R 0.51 (range 0.18–0.80); all upper bounds < 1 | `a12_allr_summary.txt` |
| Excluding the 50 genes with the largest excess | 0.48 [0.42, 0.53] | `a12_lto_summary.txt` |
| AlphaFold pLDDT ≥ 70 / < 70 | 0.53 [0.50, 0.57] / 0.39 [0.28, 0.51] | `a12_plddt_summary.txt` |
| Buried / exposed sites | 0.55 [0.52, 0.59] / 0.36 [0.29, 0.43] | `a12_burial_summary.txt` |
| PP3 × PM1-like × PP2-like (R3), DD / autism | 0.13 [0.11, 0.14] / 0.18 [0.13, 0.24] | `triple_summary.txt` |

**Predictor sets.** The 40-predictor set uses all analysed sites and includes popEVE, phyloP470way and fathmm-MKL. The 39-predictor set in Supplementary Note 2 requires constraint annotation, which admits SIFT and PolyPhen2_HVAR instead. Thirty-six predictors are common to both.

**Mechanism.** AlphaMissense flags 51% of mutation-expected sites at PM1-like positions, against 19% overall. Among excess de novos the rise is only from 79% to 94%.

---

## Supplementary Note 4. ClinVar and functional assays

**ClinVar** (P/LP against B/LB, DD de novo sites removed). AlphaMissense PP3 × PM1-like gave R 0.32 [0.28, 0.37] at Supporting and 0.23 [0.18, 0.28] at Strong in DDG2P genes, and 0.20 and 0.15 across all genes. These estimates are partly circular.

**Functional assays** (`make_mave_all.py`).
- **Mapping.** UniProt identifiers were mapped to gene symbols through the UniProt REST API. Of 82 mapped assays in 70 genes, 16 were informative (PM1-like share 5–95%; ≥ 10 damaging and ≥ 10 non-damaging variants).
- **Pooled R** (random effects): 0.94 [0.88, 0.99] at Supporting and 0.98 [0.92, 1.03] at Strong outside DDG2P genes (13 assays), and 0.89 [0.83, 0.95] across all 16.
- **The PM1-like flag barely predicted damage within these assays:** median LR 1.10 outside DDG2P genes (IQR 1.04–1.42).

**Table S5. Informative assays.**

{TAS}

---

## Supplementary Note 5. Expert-panel practice

**Data.** ClinGen Evidence Repository, export of 7 October 2026 (`make_stacking_audit.py`, `make_stacking_cv.py`, `make_erepo_eligibility.py`).

**Co-application.** Of 6,177 missense classifications:
- PP3 was applied in 2,990, PM1 in 1,092 and PP2 in 592;
- PP3 + PM1 in 787;
- PM1 + PP2 in 174: RASopathy 62, FBN1 61, Monogenic Diabetes 22, Brain Malformations 11, PTEN 11, Pulmonary Hypertension 5, ACADVL 2;
- all three in 138.

**Truth.** Points without PP3, BP4, PM1, PP2, BP1, BP3, PM4 or PM5 gave 2,942 classifications in 130 genes. R(panel PP3, panel PM1) was 0.21 [0.12, 0.31]. PM1 × PP2 and PP3 × PP2 were not estimable (1 and 2 B/LB in the joint cells).

**Table S6. PP3 credit beside panel PM1** (points [95% CI]; gene-clustered).

| Truth definition | Any strength | Supporting | Moderate |
|---|---|---|---|
| Primary | 0.60 [0.14, 1.05] (304; 17 B/LB) | 0.81 [0.24, 1.68] (235; 11) | 0.02 [−1.69, 0.89] (67; 6) |
| P/LP ≥ 4 | 0.61 [0.15, 1.07] | 0.82 [0.28, 1.77] | 0.05 [−1.13, 0.80] |
| P/LP ≥ 8 | 0.61 [0.13, 1.04] | 0.89 [0.30, 1.84] | −0.23 [−2.61, 0.59] |
| B/LB ≤ −3 | 0.65 [0.16, 1.30] | 0.83 [0.30, 2.06] | 0.14 [−2.08, 2.34] |
| PS2, PM6 also removed | 0.58 [0.08, 1.04] | 0.75 [0.17, 1.56] | 0.15 [−2.10, 0.97] |
| PS1, PP5, BP6, PS2, PM6 also removed | 0.57 [0.08, 1.04] | 0.73 [0.21, 1.60] | 0.17 [−2.06, 1.03] |
| Frequency criteria also removed | 0.49 [−0.06, 1.31] (6 B/LB) | not estimable | not estimable |

**Moderate against Supporting.** The difference was −0.79 [−2.81, 0.27], bootstrap p = 0.14. With a standard error of 1.55, the minimum detectable difference at 80% power is 4.35 points.

**Out of sample.**
- **Gene splits.** Across 200 random gene splits, the credit was 0.61 on the estimating half and 0.59 on the held-out half. The held-out credit fell below the points applied in 98% of splits.
- **Temporal split.** Estimated before 2025 the credit was 0.70 [−0.14, 1.26]; in 2025–2026 it was 0.40 [−0.25, 1.07], against 1.40 applied.

**Table S7. PP3 defined by score rather than by panel choice** (classifications mapped to analysed sites; truth also without PS2 and PM6; 2,381 classifications, 103 genes).

| Definition | Points [95% CI] | n | B/LB |
|---|---|---|---|
| Panel-applied PP3 | 0.88 [0.35, 1.49] | 254 | 13 |
| REVEL ≥ 0.644 or AlphaMissense ≥ 0.792 | 0.71 [0.35, 1.20] | 322 | 19 |
| REVEL ≥ 0.644 | 0.68 [0.34, 1.19] | 315 | 19 |
| AlphaMissense ≥ 0.792 | 1.66 [0.98, 2.58] | 264 | 7 |

**Selection into PP3.** Within panel-PM1 classifications, the likelihood ratio for pathogenicity was:
- 2.16 where PP3 was applied (n 254);
- 1.05 where the variant was score-eligible but PP3 was not applied (n 70);
- 0.04 where it was ineligible (n 21).

---

## Supplementary Note 6. Laboratory practice

**Data.** ClinVar October 2026 germline submissions for missense SNVs, excluding expert-panel and practice-guideline submissions, whose comments list ACMG/AMP codes. This gave 87,873 submissions from 472 submitters (`make_lab_stacking.py`).

**Table S8. Parser validation** (`make_parser_validation.py`; free-text summaries of 13,091 Evidence Repository classifications against their structured codes).

| Code | Precision | Recall | Strength agreement |
|---|---|---|---|
| All codes | 0.92 | 0.98 | — |
| PP3 | 0.95 | 0.99 | 0.99 |
| PM1 | 0.98 | 0.99 | 0.99 |
| PP2 | 0.98 | 0.99 | 1.00 |

**Co-application.**
- PP3 + PM1 in 4,471 submissions (Supporting 4,207, Moderate 184, Strong 80);
- PM1 + PP2 in 1,481;
- all three in 1,035.

**Truth.** 39,015 submissions: P/LP 9,455 and B/LB 29,560, in 5,976 genes from 287 submitters.

**Table S9. Laboratory truth: estimates and clustering.**

| Estimate | Value, gene-clustered 95% CI | Submitter-clustered 95% CI | PS1, PP5, BP6, PS2, PM6 also removed |
|---|---|---|---|
| R(PP3, PM1) | 0.18 [0.14, 0.25] | [0.08, 0.43] | 0.18 [0.14, 0.24] |
| PP3 Supporting credit beside PM1 | 0.84 [0.52, 1.26] | [0.06, 1.61] | 0.82 [0.52, 1.20] |
| R(PM1, PP2) | 0.27 [0.19, 0.44] | [0.14, 0.54] | 0.30 [0.21, 0.48] |
| PP2 credit beside PM1 | −0.22 [−0.67, 0.40] | [−1.23, 0.79] | −0.41 [−0.89, 0.24] |
| R(PP3, PP2) | 0.42 [0.35, 0.51] | — | — |

The last column uses its own bootstrap, so its intervals differ slightly from the first.

**Table S10. Replication in independent halves of the laboratory data** (gene-clustered 95% CI; 500 replicates).

| Half | n (PM1 + PP2 joint; B/LB) | PP2 credit beside PM1 | R(PM1, PP2) | PP3 Supporting credit | R(PP3, PM1) |
|---|---|---|---|---|---|
| Evaluated up to 2025 | 30,257 (473; 38) | 0.07 [−0.56, 1.18] | 0.32 [0.20, 0.70] | 0.94 [0.58, 1.43] | 0.20 [0.15, 0.30] |
| Evaluated after 2025 | 8,590 (79; 18) | −1.23 [−1.93, −0.57] | 0.17 [0.11, 0.28] | 0.26 [−0.17, 0.78] | 0.12 [0.08, 0.18] |
| Submitter half A | 27,449 (323; 22) | −0.77 [−1.26, −0.10] | 0.19 [0.13, 0.33] | 0.24 [−0.09, 0.70] | 0.09 [0.07, 0.13] |
| Submitter half B | 11,566 (242; 35) | 0.14 [−0.52, 1.12] | 0.31 [0.19, 0.63] | 1.21 [0.72, 1.86] | 0.33 [0.23, 0.56] |

**Random submitter halves.** In 100 random halves of the submitters, the held-out PP2 credit had median −0.22, with an upper 95% bound below 1 in 62% of halves.

**Equivalence of Supporting PP3 credit to 1 point** (margin ±0.5; 90% CI):
- gene-clustered [0.58, 1.17]: equivalent;
- submitter-clustered [0.10, 1.47]: not shown.

**Thin cells.** Credit for PP3 at Moderate (3 B/LB) and Strong (0 B/LB) is not estimable in laboratory data.

**Table S11. Laboratory re-scoring** (changed classifications, of 87,873).

| Rule | Changed | P→LP | LP→VUS | Other |
|---|---|---|---|---|
| ClinGen limit | 40 | 30 | 10 | 0 |
| PP3 + PM1 + PP2 ≤ 4 | 63 | 39 | 22 | 2 |
| PP2 not added when PM1 applies | 352 | 124 | 223 | 5 |
| PP3 limited to Supporting when PM1 applies | 90 | 39 | 51 | 0 |

**An earlier version** restricted the data to submissions that mention PP3 or PM1. That selects on the exposure, so it is not used.

---

## Supplementary Note 7. Consequences in developmental disorders and ClinVar

**Developmental disorders.** In DD, uncapped co-use of PS2 + PM2 + PM1-like + AlphaMissense PP3 placed 2,438 de novo variants at Pathogenic, in a class the de novo excess implies is 84% disease-associated. This is a class-level estimate of association with any severe developmental disorder, not a per-variant posterior.

**ClinVar.** In DDG2P monoallelic genes, 37.3% of AlphaMissense PP3-positive VUS lie at PM1-like sites, against 13.2% across all genes. Without the ClinGen limit, PP3 + PM1 at Moderate alone would reach Likely Pathogenic for 7,101 of them.

---

## Supplementary Data files

| File | Content |
|---|---|
| Supplementary Data 1 (`a12_gclust2.tsv`) | All de novo dependence estimates (Table S3) |
| Supplementary Data 2 (`mave_all_sites_summary.tsv`) | All mapped functional assays |
| Supplementary Data 3 (`supp_table1_predictor_resource.tsv`) | Per-predictor retention with intervals |
| Supplementary Data 4 (`lab_stacking_summary.txt`, `stacking_cv_summary.txt`) | Expert-panel and laboratory outputs |
"""
open("SUPPLEMENT_v13.md", "w", encoding="utf-8").write(S)
print(len(S.split()), "words")
