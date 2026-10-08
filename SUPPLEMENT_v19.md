# Supplementary material: Dependence between computational and location evidence in missense variant classification: a label-free measurement, its cost, and a pre-registered test of practice

*Draft 19, 8 October 2026. Accompanies `MANUSCRIPT_v19.md`. Complete outputs are in `results_fullstudy/` in the study repository; each note names its files.*

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

| Cohort | Predictor | PP3 tier | Flag | R | f | Within gene, all | Within gene, DDG2P |
|---|---|---|---|---|---|---|---|
| DD | AlphaMissense | Supporting | PM1-like | 0.44 [0.41, 0.47] | 0.77 [0.76, 0.78] | 0.75 [0.71, 0.78] | 0.69 [0.65, 0.73] |
| DD | AlphaMissense | Strong | PM1-like | 0.33 [0.31, 0.35] | 0.73 [0.72, 0.75] | 0.70 [0.66, 0.74] | 0.65 [0.60, 0.69] |
| DD | AlphaMissense | Supporting | PM1-strict | 0.37 [0.35, 0.40] | 0.77 [0.75, 0.78] | 0.67 [0.63, 0.71] | 0.61 [0.57, 0.65] |
| DD | AlphaMissense | Strong | PM1-strict | 0.25 [0.23, 0.28] | 0.72 [0.70, 0.74] | 0.58 [0.54, 0.63] | 0.53 [0.48, 0.58] |
| DD | AlphaMissense | Supporting | PP2-like | 0.61 [0.57, 0.65] | 0.85 [0.83, 0.86] | 0.96 [0.93, 1.00] | 0.95 [0.92, 0.97] |
| DD | AlphaMissense | Strong | PP2-like | 0.52 [0.48, 0.56] | 0.83 [0.81, 0.84] | 0.97 [0.93, 1.02] | 0.96 [0.92, 0.99] |
| DD | REVEL | Supporting | PM1-like | 0.57 [0.53, 0.61] | 0.84 [0.82, 0.86] | 0.78 [0.73, 0.83] | 0.72 [0.68, 0.77] |
| DD | REVEL | Strong | PM1-like | 0.53 [0.48, 0.58] | 0.84 [0.82, 0.86] | 0.77 [0.71, 0.84] | 0.68 [0.63, 0.75] |
| DD | REVEL | Supporting | PM1-strict | 0.48 [0.43, 0.54] | 0.83 [0.80, 0.85] | 0.70 [0.64, 0.76] | 0.64 [0.58, 0.71] |
| DD | REVEL | Strong | PM1-strict | 0.45 [0.37, 0.54] | 0.83 [0.79, 0.87] | 0.69 [0.60, 0.79] | 0.62 [0.52, 0.74] |
| DD | REVEL | Supporting | PP2-like | 0.73 [0.67, 0.79] | 0.90 [0.88, 0.93] | 0.98 [0.93, 1.04] | 0.97 [0.92, 1.01] |
| DD | REVEL | Strong | PP2-like | 0.75 [0.67, 0.84] | 0.92 [0.90, 0.95] | 1.03 [0.96, 1.12] | 0.99 [0.93, 1.05] |
| DD | VEST4 | Supporting | PM1-like | 0.67 [0.63, 0.70] | 0.87 [0.86, 0.89] | 0.88 [0.83, 0.93] | 0.80 [0.75, 0.84] |
| DD | VEST4 | Strong | PM1-like | 0.72 [0.66, 0.78] | 0.91 [0.88, 0.93] | 0.90 [0.82, 0.97] | 0.77 [0.71, 0.84] |
| DD | VEST4 | Supporting | PM1-strict | 0.60 [0.55, 0.64] | 0.87 [0.84, 0.88] | 0.81 [0.75, 0.87] | 0.73 [0.67, 0.80] |
| DD | VEST4 | Strong | PM1-strict | 0.59 [0.50, 0.67] | 0.87 [0.84, 0.91] | 0.75 [0.66, 0.86] | 0.63 [0.53, 0.74] |
| DD | VEST4 | Supporting | PP2-like | 0.77 [0.72, 0.82] | 0.91 [0.89, 0.93] | 1.01 [0.97, 1.07] | 0.98 [0.94, 1.01] |
| DD | VEST4 | Strong | PP2-like | 0.93 [0.84, 1.01] | 0.98 [0.95, 1.00] | 1.12 [1.04, 1.20] | 1.01 [0.94, 1.07] |
| autism | AlphaMissense | Supporting | PM1-like | 0.56 [0.47, 0.66] | 0.79 [0.75, 0.84] | 0.96 [0.84, 1.10] | 0.74 [0.67, 0.82] |
| autism | AlphaMissense | Strong | PM1-like | 0.44 [0.36, 0.54] | 0.75 [0.71, 0.79] | 0.92 [0.80, 1.07] | 0.73 [0.64, 0.82] |
| autism | AlphaMissense | Supporting | PM1-strict | 0.50 [0.42, 0.61] | 0.79 [0.74, 0.84] | 0.88 [0.75, 1.03] | 0.68 [0.60, 0.77] |
| autism | AlphaMissense | Strong | PM1-strict | 0.39 [0.30, 0.50] | 0.75 [0.69, 0.80] | 0.83 [0.67, 1.02] | 0.68 [0.55, 0.83] |
| autism | AlphaMissense | Supporting | PP2-like | 0.65 [0.55, 0.76] | 0.84 [0.79, 0.89] | 1.05 [0.93, 1.18] | 0.90 [0.84, 0.97] |
| autism | AlphaMissense | Strong | PP2-like | 0.57 [0.47, 0.68] | 0.82 [0.78, 0.87] | 1.07 [0.93, 1.22] | 0.90 [0.81, 0.98] |
| autism | REVEL | Supporting | PM1-like | 0.82 [0.69, 0.99] | 0.92 [0.85, 1.00] | 1.08 [0.92, 1.30] | 0.74 [0.65, 0.84] |
| autism | REVEL | Strong | PM1-like | 0.74 [0.57, 0.97] | 0.90 [0.82, 0.99] | 1.05 [0.83, 1.45] | 0.66 [0.56, 0.80] |
| autism | REVEL | Supporting | PM1-strict | 0.73 [0.56, 0.96] | 0.89 [0.82, 0.98] | 0.95 [0.74, 1.22] | 0.66 [0.53, 0.80] |
| autism | REVEL | Strong | PM1-strict | 0.68 [0.42, 1.02] | 0.89 [0.77, 1.01] | 0.96 [0.63, 1.40] | 0.59 [0.38, 0.84] |
| autism | REVEL | Supporting | PP2-like | 0.93 [0.78, 1.13] | 0.97 [0.90, 1.06] | 1.19 [1.02, 1.42] | 0.90 [0.80, 0.99] |
| autism | REVEL | Strong | PP2-like | 0.94 [0.73, 1.21] | 0.98 [0.90, 1.07] | 1.19 [0.95, 1.53] | 0.90 [0.79, 1.02] |
| autism | VEST4 | Supporting | PM1-like | 0.77 [0.64, 0.91] | 0.89 [0.83, 0.96] | 1.06 [0.92, 1.25] | 0.86 [0.77, 0.97] |
| autism | VEST4 | Strong | PM1-like | 0.84 [0.67, 1.05] | 0.94 [0.86, 1.02] | 1.13 [0.92, 1.43] | 0.83 [0.71, 0.98] |
| autism | VEST4 | Supporting | PM1-strict | 0.79 [0.64, 0.97] | 0.92 [0.85, 0.99] | 1.13 [0.92, 1.37] | 0.89 [0.73, 1.06] |
| autism | VEST4 | Strong | PM1-strict | 0.76 [0.51, 1.08] | 0.92 [0.80, 1.03] | 1.04 [0.72, 1.43] | 0.78 [0.53, 1.05] |
| autism | VEST4 | Supporting | PP2-like | 0.80 [0.69, 0.95] | 0.91 [0.85, 0.97] | 1.09 [0.96, 1.25] | 0.97 [0.89, 1.06] |
| autism | VEST4 | Strong | PP2-like | 1.02 [0.82, 1.29] | 1.01 [0.93, 1.11] | 1.24 [1.03, 1.53] | 0.95 [0.83, 1.08] |
| public | AlphaMissense | Supporting | PM1-like | 0.67 [0.48, 0.96] | 0.83 [0.73, 0.98] | 0.96 [0.74, 1.28] | 0.75 [0.62, 0.89] |
| public | AlphaMissense | Strong | PM1-like | 0.48 [0.34, 0.70] | 0.75 [0.66, 0.85] | 0.83 [0.60, 1.13] | 0.64 [0.48, 0.79] |
| public | AlphaMissense | Supporting | PM1-strict | 0.55 [0.38, 0.83] | 0.81 [0.69, 0.92] | 0.84 [0.61, 1.20] | 0.73 [0.61, 0.95] |
| public | AlphaMissense | Strong | PM1-strict | 0.34 [0.20, 0.54] | 0.70 [0.57, 0.81] | 0.63 [0.40, 0.99] | 0.57 [0.40, 0.87] |
| public | AlphaMissense | Supporting | PP2-like | 0.78 [0.57, 1.07] | 0.89 [0.78, 1.04] | 1.00 [0.80, 1.29] | 0.93 [0.83, 1.05] |
| public | AlphaMissense | Strong | PP2-like | 0.69 [0.49, 0.97] | 0.87 [0.78, 0.99] | 1.01 [0.80, 1.30] | 0.93 [0.82, 1.06] |
| public | REVEL | Supporting | PM1-like | 0.78 [0.54, 1.15] | 0.90 [0.77, 1.07] | 0.89 [0.66, 1.25] | 0.78 [0.64, 0.94] |
| public | REVEL | Strong | PM1-like | 0.71 [0.45, 1.09] | 0.89 [0.74, 1.04] | 0.75 [0.50, 1.13] | 0.77 [0.57, 1.00] |
| public | REVEL | Supporting | PM1-strict | 0.54 [0.26, 0.88] | 0.80 [0.56, 0.95] | 0.66 [0.35, 0.99] | 0.65 [0.40, 0.90] |
| public | REVEL | Strong | PM1-strict | 0.48 [0.15, 0.94] | 0.81 [0.51, 0.98] | 0.55 [0.21, 1.06] | 0.62 [0.28, 1.04] |
| public | REVEL | Supporting | PP2-like | 1.02 [0.73, 1.43] | 1.01 [0.88, 1.19] | 1.18 [0.93, 1.59] | 0.96 [0.86, 1.07] |
| public | REVEL | Strong | PP2-like | 0.81 [0.48, 1.19] | 0.93 [0.76, 1.07] | 0.85 [0.57, 1.16] | 0.85 [0.64, 1.01] |
| public | VEST4 | Supporting | PM1-like | 0.75 [0.54, 1.04] | 0.87 [0.76, 1.02] | 1.02 [0.77, 1.36] | 0.91 [0.75, 1.11] |
| public | VEST4 | Strong | PM1-like | 0.84 [0.56, 1.26] | 0.94 [0.80, 1.11] | 1.09 [0.76, 1.58] | 0.85 [0.66, 1.10] |
| public | VEST4 | Supporting | PM1-strict | 0.63 [0.43, 0.95] | 0.85 [0.72, 0.98] | 0.96 [0.67, 1.45] | 0.86 [0.66, 1.26] |
| public | VEST4 | Strong | PM1-strict | 0.80 [0.45, 1.37] | 0.93 [0.77, 1.11] | 1.16 [0.68, 2.00] | 0.96 [0.52, 1.73] |
| public | VEST4 | Supporting | PP2-like | 0.74 [0.54, 0.97] | 0.87 [0.73, 0.99] | 0.95 [0.72, 1.22] | 0.92 [0.79, 1.05] |
| public | VEST4 | Strong | PP2-like | 0.97 [0.67, 1.41] | 0.99 [0.86, 1.15] | 1.12 [0.82, 1.55] | 0.96 [0.80, 1.15] |

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

| Assay | Gene set | Variants | PM1-like share | LR of PM1-like | R, Supporting [95% CI] | R, Strong |
|---|---|---|---|---|---|---|
| EPHB2_HUMAN_Tsuboyama_2023_1F0M | other | 376 | 0.21 | 1.27 | 0.40 [0.24, 0.62] | 0.33 |
| HXK4_HUMAN_Gersing_2022_activity | other | 2413 | 0.13 | 1.64 | 0.84 [0.59, 1.50] | 0.76 |
| HXK4_HUMAN_Gersing_2023_abundance | other | 2388 | 0.13 | 0.73 | 0.79 [0.59, 1.06] | 0.95 |
| HMDH_HUMAN_Jiang_2019 | other | 4994 | 0.46 | 3.14 | 1.40 [1.16, 1.66] | 1.56 |
| ILF3_HUMAN_Tsuboyama_2023_2L33 | other | 399 | 0.57 | 1.10 | 0.92 [0.82, 1.03] | 0.86 |
| MK01_HUMAN_Brenan_2016 | other | 2070 | 0.81 | 1.04 | 0.97 [0.93, 1.02] | 0.99 |
| MET_HUMAN_Estevam_2023 | other | 1600 | 0.58 | 1.42 | 0.83 [0.76, 0.92] | 0.85 |
| PPARG_HUMAN_Majithia_2016 | other | 2943 | 0.12 | 2.24 | 0.65 [0.52, 0.81] | 0.69 |
| RBP1_HUMAN_Tsuboyama_2023_2KWH | other | 278 | 0.23 | 1.25 | 0.60 [0.38, 0.89] | 0.71 |
| SRC_HUMAN_Ahler_2019 | other | 1038 | 0.89 | 1.01 | 1.00 [0.96, 1.04] | 1.01 |
| SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM | other | 1107 | 0.89 | 1.04 | 1.00 [0.96, 1.06] | 1.02 |
| SRC_HUMAN_Nguyen_2022 | other | 1036 | 0.89 | 1.00 | 1.02 [0.98, 1.07] | 1.02 |
| UBR5_HUMAN_Tsuboyama_2023_1I2T | other | 325 | 0.70 | 1.05 | 0.95 [0.88, 1.03] | 0.91 |
| RASH_HUMAN_Bandaru_2017 | DDG2P | 939 | 0.45 | 1.53 | 0.74 [0.66, 0.82] | 0.69 |
| PPM1D_HUMAN_Miller_2022 | DDG2P | 2386 | 0.73 | 1.45 | 0.72 [0.65, 0.79] | 0.71 |
| SHOC2_HUMAN_Kwon_2022 | DDG2P | 3294 | 0.19 | 1.28 | 0.86 [0.71, 1.02] | 0.82 |

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

## Supplementary Note 8. Curated PM1, within-submitter and frequency-handled estimates, a held-out test, and trends

**Curated PM1 on functional assays** (`make_mave_curated.py`; `mave_curated_summary.txt`). Curated PM1 residues are the protein positions at which any expert panel applied PM1 (any strength) to a missense variant in the Evidence Repository: 1,439 classifications, 841 residues, 67 genes. Six of these genes had mapped ProteinGym assays. GCK, TP53, HRAS and KRAS passed the reference check (≥ 90% of curated residues matching the transcript's reference amino acid); PTEN (0.86) and RAF1 (0.88) did not. Seven assays in three genes had at least 10 variants at curated residues. BRCA1 has no curated PM1 residues. With the PM1-like proxy, TP53 could not be tested because the proxy flags none of its assayed residues; BRCA1 had 3% flagged and PTEN nearly all.

**Table S12. AlphaMissense PP3 × curated PM1 on functional truth** (position-clustered 95% CI; 500 replicates).

| Assay | Gene | Variants | At curated PM1 | LR(curated PM1) | LR(PM1-like) | R, Supporting | R, Strong |
|---|---|---|---|---|---|---|---|
| Gersing 2022, activity | GCK | 2,413 | 72 (3%) | 4.31 | 1.64 | 0.40 [0.27, 1.06] | 0.51 |
| Gersing 2023, abundance | GCK | 2,388 | 72 (3%) | 1.00 | 0.73 | 0.38 [0.30, 0.46] | 0.36 |
| Bandaru 2017 | HRAS | 939 | 23 (2%) | 0.85 | 1.53 | 0.60 [0.53, 0.67] | 0.46 |
| Giacomelli 2018, null, etoposide | TP53 | 2,020 | 199 (10%) | 3.87 | — | 0.43 [0.32, 0.64] | 0.39 |
| Giacomelli 2018, null, nutlin | TP53 | 2,020 | 199 (10%) | 3.16 | — | 0.52 [0.35, 0.92] | 0.51 |
| Giacomelli 2018, WT, nutlin | TP53 | 2,020 | 199 (10%) | 4.71 | — | 0.40 [0.29, 0.59] | 0.33 |
| Kotler 2018 | TP53 | 933 | 175 (19%) | 2.32 | — | 0.64 [0.47, 0.98] | 0.65 |
| **Pooled (random effects)** | | | | median 3.16 | | **0.49 [0.40, 0.60]** | **0.43 [0.36, 0.51]** |

**PP3 credit at curated PM1 residues.** Within each assay, AlphaMissense Supporting PP3 had a median likelihood ratio for damage of 1.78 at curated PM1 residues, against 3.41 across the assay. The median share of log-evidence kept was 0.40 (GCK activity 0.31, GCK abundance −0.01, TP53 0.40–0.61, HRAS 0.00). At curated residues 77% of variants were PP3-positive (34% elsewhere) and 51% damaging (26% elsewhere). Where curated PM1 itself carried no information about damage (GCK abundance, HRAS), R < 1 arises because PP3 adds nothing at those residues. For classification this is the same loss of credit as duplication.

"—": no assayed residue is PM1-like. Assays within a gene share residues and are not independent; the three-gene result should be read with that in mind. Panels apply PM1 where pathogenic variants have been classified, so curated residues are enriched for known pathogenic positions; the curated flag uses no assay or predictor data.

**Within-submitter credit** (`make_review20_lab.py`; `review20_lab_summary.txt`). Among 2,812 PM1-applied laboratory submissions with truth, 35 submitters had both P/LP and B/LB submissions (2,043 submissions). Credit was pooled across them with the Mantel–Haenszel risk-ratio estimator (truth as exposure, criterion applied as outcome, submitter as stratum), with a submitter-clustered bootstrap (1,000 replicates).

**Table S13. Credit once PM1 applies: pooled against within submitter** (laboratory truth; points).

| Criterion | Pooled | Within submitter [95% CI] | 90% CI | Applied |
|---|---|---|---|---|
| PP2 | −0.22 | 0.37 [−0.36, 0.95] | −0.31, 0.89 | 1.00 |
| PP3 Supporting | 0.84 | 0.93 [−0.03, 1.62] | 0.03, 1.53 | 1.00 |
| PP3, any strength | 0.85 | 1.00 [0.05, 1.68] | 0.11, 1.60 | 1.09 |

**PP2 credit with frequency evidence handled** (`make_review21_lab.py`; `review21_lab_summary.txt`). Laboratory benign truth among PM1 submissions rests mostly on frequency: 238 B/LB, of which only 43 lack BA1, BS1 and BS2.

**Table S14. PP2 credit beside PM1 by handling of frequency evidence** (points; applied 1; pooled gene-clustered and within-submitter submitter-clustered 95% CI, 1,000 replicates).

| Subset | Submissions | B/LB | PP2 & B/LB | Pooled | Within submitter |
|---|---|---|---|---|---|
| All PM1 submissions | 2,812 | 238 | 57 | −0.22 [−0.69, 0.37] | 0.37 [−0.41, 0.96] |
| No BA1/BS1/BS2 applied | 2,612 | 43 | 15 | −0.65 [−1.39, 21.8] | 0.90 [−0.49, 2.18] |
| PM2 applied | 2,322 | 39 | 21 | −1.08 [−1.48, −0.63] | −0.42 [−0.90, 0.63] |
| B/LB requires a non-frequency benign criterion | 2,646 | 72 | 24 | −0.60 [−1.24, 0.79] | 0.16 [−0.71, 1.14] |

The first row repeats Table S13 with a fresh bootstrap. Intervals widen once frequency is handled, and none of the within-submitter intervals excludes 1.

**Held-out test of "PP2 not added when PM1 applies"** (`make_review21_lab.py`). ClinVar submission_summary of January 2023 was parsed with the same filters and parser as the 2026 laboratory set: 8,510 germline, non-expert-panel missense SNV submissions listing ACMG codes, 1,748 of them Likely Pathogenic. Variants in the 2026 PM1-applied truth set used to estimate PP2 credit (2,445 variants) were excluded.
- **Rule-affected:** Likely Pathogenic, with PM1 and PP2, and below 6 points once PP2 is removed (41).
- **Controls A:** Likely Pathogenic at 6 points without PP2 beside PM1 (450).
- **Controls B:** the subset of A with PM1 (125).
- **Primary outcome (O2):** among other submitters' October 2026 classifications of the variant, VUS/LB/B outnumber P/LP.
- **Secondary outcome (O1):** the same SCV accession is no longer P/LP. This cannot be held out, because it requires the submission to persist.
- **Statistics:** exact conditional odds ratio with exact central 95% CI and the matching central exact p, plus Mantel–Haenszel odds ratios (Robins–Breslow–Greenland CI) stratified by gene and by submitter.

**Table S15. Later reclassification of 2023 Likely Pathogenic calls.**

| Outcome, controls | Rule-affected | Controls | Exact OR [95% CI], p | MH by gene | MH by submitter |
|---|---|---|---|---|---|
| **O2 (primary), controls A** | 9/34 (26%) | 48/411 (12%) | 2.71 [1.05, 6.46], 0.039 | 7.43 [0.98, 56.5] (8 genes) | 4.69 [1.39, 15.9] (13 submitters) |
| O2, controls B (with PM1) | 9/34 (26%) | 17/113 (15%) | 2.02 [0.71, 5.51], 0.21 | 2.33 [0.26, 21.4] | 4.17 [0.88, 19.7] |
| O1, controls A | 2/41 (4.9%) | 2/442 (0.5%) | 11.2 [0.79, 158], 0.076 | not estimable | 10.8 [0.81, 144] |
| O1, controls B | 2/41 (4.9%) | 1/122 (0.8%) | 6.12 [0.31, 369], 0.31 | not estimable | not estimable |

The primary comparison is significant and the others point the same way. The comparison with PM1-matched controls is underpowered and does not separate PP2 from location evidence in general.

**Reach and trend** (`review21_lab_summary.txt`). Of the 1,481 laboratory classifications combining PM1 and PP2, 352 (23.8%) change class if PP2 is not added. Of their 624 Likely Pathogenic calls, 223 (35.7%) become VUS, and of their 536 Pathogenic calls, 124 become Likely Pathogenic. The changes span 160 genes and 66 laboratories. The leading genes are NOTCH3 (24), FBN1 (21), COL4A3 (12), CASR (10), POLG (10), CACNA1A (9), KCNH2 (9), SCN1A (9), COL3A1 (8), PTEN (6), COL4A4 (6) and MYH7 (6).

**Table S16. Co-application among laboratory PM1 submissions by year of last evaluation.**

| Year | PM1 submissions | PP2 beside PM1 | PP3 beside PM1 | PP3 ≥ Moderate beside PM1 |
|---|---|---|---|---|
| 2018 | 121 | 37% | 79% | 0.0% |
| 2019 | 182 | 15% | 65% | 0.5% |
| 2020 | 274 | 16% | 78% | 0.4% |
| 2021 | 315 | 17% | 81% | 0.0% |
| 2022 | 713 | 20% | 69% | 0.8% |
| 2023 | 1,013 | 18% | 66% | 0.9% |
| 2024 | 1,469 | 19% | 63% | 4.3% |
| 2025 | 1,712 | 24% | 62% | 8.5% |
| 2026 (to October) | 796 | 31% | 60% | 3.8% |

Submission-level Spearman trend: PP2 beside PM1, ρ = 0.075 (p = 1.2 × 10⁻⁹); PP3 ≥ Moderate beside PM1, ρ = 0.127 (p = 2.8 × 10⁻²⁵). The year is the date of last evaluation, so laboratories that re-evaluate often are over-represented in recent years.

---

## Supplementary Table S0. Hypotheses and when they were specified relative to data access

| Hypothesis | Specified | Data first tabulated for it | Status |
|---|---|---|---|
| PP3 × PM1 dependence on the de novo scale (R < 1) | analysis log, before results (Note 1) | after specification | pre-specified; confirmed |
| Predictor evidence incremental to gene and position (primary criterion) | analysis log | after specification | pre-specified; not met |
| Expert-panel and laboratory credit, within-submitter, frequency-handled, held-out and trend analyses (drafts 6–16) | not pre-specified | before specification | exploratory |
| A16 hierarchical label-free model | A16 pre-registration (commit 34bac58) | after | pre-registered; refuted |
| A18 temporal reclassification | A18 pre-registration (commit 8a57850) | after | pre-registered; null |
| H1 rules for PP3 beside PM1 (primary) | A19 (commit ceb1a11), 8 Oct 2026 | after | pre-registered; not met (cap retained) |
| H2 PP2 credit beside PM1 < 1 (primary) | A19 | after (the pooled −0.22 and within-submitter 0.37 estimates of drafts 13–16 were known; the mixed model, its decision rule and random effects were new) | pre-registered; confirmed |
| H2b PM1-rationale interaction | A19 | after | pre-registered; not testable (17 constraint-based) |
| H3 calibrate rather than cap | A19 | after | pre-registered; dependence persists |
| H4 replication across eras and negative controls | A19 | after | pre-registered; R < 1 in all eras; negative controls non-specific |
| H5a rules for PP3 beside PM1 on de novo truth, primary | A20 (commit 2e9692a) | after (cell estimates known from drafts 13–17; no rule evaluated) | pre-registered; not met (selected rule overfitted) |
| H5 secondary rule contrasts; H5b PP2; H5c triple | A20 | after | secondary / descriptive: one level lower beats cap (AlphaMissense only); no PP2 beside PM1 and cap + no PP2 best calibrated |
| A21 PP3 × curated PM1 on MaveDB, random-residue control | A21 (commit 8cdd9c2) | after (score-set list known; no estimate) | pre-registered; confirmed (R 0.59 [0.44, 0.80]); control passed (R 1.00) |
| Cost headline, validity analyses (Notes 8–9) | not pre-specified | — | exploratory / descriptive |

## Supplementary Note 9. Pre-registered tests (A19)

**H1 (`make_rules_cv.py`; `a19_h1_rules_2026-10_summary.txt`, `..._acmg_summary.txt`).** Held-out log-loss among 2,812
PM1-applied laboratory submissions with truth (gene-grouped 5-fold CV × 20): uncapped 0.2763, cap 0.2769, Supporting limit
0.2770, discount 0.2732, learned credit 0.2732, no credit 0.3147. Differences vs cap (submitter-clustered bootstrap):
learned credit −0.0037 [−0.0101, +0.0105]; discount −0.0037 [−0.0092, +0.0061]; uncapped −0.0006 [−0.0015, −0.0002];
no credit +0.0378 [+0.0239, +0.0499]. Deviation 1 (ACMG convention, absent code carries no evidence): all differences
within ±0.001 of the cap except no credit (+0.0088 [−0.0076, +0.0178]); learned credit beside PM1 0.83 (Supporting),
0.57 (Moderate), 1.93 (Strong) points; discount factor 0.78.

**H2 (`make_pp2_confirm.py`; `a19_h2_pp2_2026-10_summary.txt`).** Variational-Bayes binomial mixed model, PM1-applied
submissions with truth (n 2,812; 238 B/LB): PP2 +0.35 points [−0.12, 0.81]; PP3 +2.34 [2.01, 2.67] (present vs absent).
Random-effect SDs: submitter 2.15, gene 1.34 (log-odds). ACMG-convention credits (gene-clustered): PP2 beside PM1 −0.23
[−0.67, 0.37] (n 565; 57 B/LB), without PM1 +1.12 [0.90, 1.31]; PP3 beside PM1 +0.85 [0.54, 1.23], without +2.82
[2.72, 2.92]. PM1 rationale classes: other 2,325, hotspot/domain 470, constraint 17; PP2 beside hotspot/domain PM1 −1.24
[−1.95, −0.55] (14 B/LB; descriptive). H2c: PP2 − PP3 beside PM1 −1.08 [−1.69, −0.36] (submitter-clustered).

**H3 (`make_genecal_pp3.py`; `a19_h3_summary.txt`).** Published Supporting threshold 0.792; per-gene thresholds for
341 genes (median 0.132, IQR 0.099–0.192: the pre-registered cumulative-LR rule is permissive). R(PP3, PM1): published
0.31 [0.18, 0.88], gene-calibrated 0.41 [0.30, 0.61]; ratio 1.33 [0.65, 1.96]; in calibrated genes 0.15 vs 0.51.

**H4 (`make_vintage_replication.py`; `a19_h4_summary.txt`).** 88,202 unique submissions (first release carrying codes:
2020-01 1,979; 2022-01 3,053; 2024-01 27,248; 2026-10 55,922); 38,911 with truth. R(PP3, PM1) by era of last
evaluation: ≤ 2019 0.27 [0.14, 0.86]; 2020–2022 0.10 [0.07, 0.18]; 2023–2024 0.24 [0.16, 0.40]; 2025–2026 0.17 [0.13,
0.26]. PP2 credit beside PM1: −0.47, +0.82, +0.45, −0.58 (only 2025–2026 had ≥ 20 B/LB). Negative controls (truth also
excluding PP1, PP4; n 38,007): credit beside / without PM1 — PP1 +1.34 / +5.40; PP4 +0.76 / +5.00; PP2 −0.22 / +1.06;
PP3 +0.84 / +2.80. R(·, PM1): PP1 0.04 [0.02, 0.22]; PP4 0.03 [0.02, 0.12]; PP2 0.29 [0.20, 0.47]; PP3 0.18 [0.14, 0.25].

## Supplementary Note 10. Cost of stacking and validity of the laboratory analysis (exploratory)

**Cost (`make_cost_headline.py`; `cost_headline_summary.txt`).** DD classes under PS2 + PM2 + PM1-like + AlphaMissense
PP3: uncapped Pathogenic 2,438 (84% [84, 85] disease-associated vs ≥ 99% nominal; implied 5.4 points vs ≥ 10); uncapped
Likely Pathogenic 5,686 (33%); capped and conditional Likely Pathogenic 8,124 (48%). Points over-counted (−ln R / ln
2.406) for every predictor-flag pair are in the summary file.

**Validity (`make_validity_v14.py`; `validity_v14_summary.txt`).** PM1 strength when applied: Supporting 556, Moderate
6,077, Strong 146. PP3 + PM1 sums: 2: 282; 3: 3,823; 4: 155; 5: 138; 6: 69; 8: 4. Reported class equals the points-based
class in 88.2% of submissions; where the cap would change the class (38), laboratories reported the uncapped class in 32
and the capped class in 6 (Mayo Clinic Laboratories 4 of 20). Templated submitters (comment-opening redundancy ≥ 0.3,
≥ 20 PP3 & PM1 submissions): 18, covering 1,681 submissions; without them R(PP3, PM1) 0.11 [0.09, 0.15], PP3 credit
beside PM1 +0.28, PP2 credit −0.80.

**Parser validation by annotation.** `parser_annotation_sheet_v14.tsv`: 400 comments from 248 submitters (≤ 5 each,
PP3/PM1/PP2 mentions prioritised). Two annotators, blind to `parser_codes`, record each code and strength; per-code
precision, recall and Cohen's kappa are then computed. Not yet done.

## Supplementary Note 11. Previous positions, and Supplementary Figure S1

**Table S17. Previous positions on overlap between PP3, PM1 and PP2, and what this study adds** (formerly Table 7).

| Source | Position or finding | Basis | This study |
|---|---|---|---|
| ClinGen SVI (Pejaver et al. 2022; Stenton et al. 2024) | PP3 + PM1 limited to Strong | Expert judgement of shared attributes | Dependence measured; 98–99% of co-application lies below the level at which the limit acts |
| ACGS 2024 (Durkie et al. 2024) | Avoid counting constraint in both PM1 and PP2 | Expert judgement | Estimate of what PP2 adds once PM1 applies |
| Sherloc (Nykamp et al. 2017) | Rule dependencies to prevent one observation counting twice | Expert judgement | Dependence quantified for one such pairing set |
| ClinGen cardiomyopathy panel (Kelly et al. 2018) | PP2 not used beside PM1 in the MYH7 head domain | Expert judgement, one gene | Tested across 5,976 genes in laboratory practice |
| CanVIG-UK (Garrett et al. 2021) | PP3 and PM1 may be combined | Expert judgement | Supporting PP3 beside PM1 not over-credited |
| Amendola et al. 2016; Harrison et al. 2017 | Laboratories apply ACMG/AMP criteria inconsistently | Inter-laboratory comparison | Co-application audited in 46 panels and 472 laboratories |
| Grimm et al. 2015; Harrizi et al. 2026 | ClinVar evaluations of predictors are circular and gene-dominated | ClinVar | Likelihood ratios without variant-level labels |
| Tejura et al. 2024; Chen et al. 2026 | Predictor calibration varies between genes and domains | ClinVar and assays | Dependence between criteria rather than between genes |


**Supplementary Figure S1.** PP3 credit beside PM1, earned against applied, by data source and analysis (`fig3_credit.png`; formerly Figure 3).

## Supplementary Note 12. Rules on de novo truth (A20)

`make_denovo_rules.py`; trained in DD, tested in autism and the public cohort; cells with ≥ 50 de novos in the joint cell and a strictly positive R interval; error in points; MSE weighted by 1/SE²; parametric bootstrap 10,000 (cells treated as independent).

```
AlphaMissense: cells train (DD) 6, test (autism, public) 5; learned kappa constant 0.40; per band nominal 1: -16.99, nominal 2: -0.18, nominal 4: 0.41
  rule                    train MSE   test MSE   test MSE - cap [95% CI]
  R0 uncapped                 1.670      0.901   +0.801 [+0.509, +1.088]
  R1 ClinGen cap              0.234      0.100   +0.000 [+0.000, +0.000]
  R2 PP3 <= Supporting        0.222      0.166   +0.067 [-0.085, +0.219]
  R3 one level lower          0.143      0.038   -0.062 [-0.109, -0.016]
  R6 no PP3 credit            0.776      0.614   +0.514 [+0.211, +0.820]
  R4 learned constant         0.099      0.040   -0.060 [-0.141, +0.025]
  R5 learned per band         0.017      0.521   +0.422 [+0.049, +0.899]  <- selected on DD
  H5a (primary): R5 learned per band vs cap on test cohorts -> does not beat the cap (upper bound +0.899)
    test cell autism  PM1-like   PP3 Supporting only  n  104: error cap -0.23 pts, R5 learned per band +3.15 pts
    test cell autism  PM1-like   PP3 Moderate only    n  131: error cap -0.93 pts, R5 learned per band +0.02 pts
    test cell autism  PM1-like   PP3 >=Strong(3-4)    n  519: error cap -0.02 pts, R5 learned per band +0.15 pts
    test cell autism  PM1-strict PP3 >=Strong(3-4)    n  157: error cap -0.16 pts, R5 learned per band +0.00 pts
    test cell public3 PM1-like   PP3 >=Strong(3-4)    n  107: error cap +0.06 pts, R5 learned per band +0.22 pts
REVEL: cells train (DD) 6, test (autism, public) 4; learned kappa constant 0.59; per band nominal 1: 0.31, nominal 2: 0.52, nominal 4: 0.66
  rule                    train MSE   test MSE   test MSE - cap [95% CI]
  R0 uncapped                 0.614      0.106   -0.020 [-0.253, +0.225]
  R1 ClinGen cap              0.391      0.126   +0.000 [+0.000, +0.000]
  R2 PP3 <= Supporting        0.565      0.400   +0.274 [+0.124, +0.424]
  R3 one level lower          0.103      0.141   +0.015 [-0.113, +0.147]
  R6 no PP3 credit            1.263      0.990   +0.864 [+0.550, +1.176]
  R4 learned constant         0.054      0.046   -0.079 [-0.175, +0.019]
  R5 learned per band         0.014      0.035   -0.090 [-0.220, +0.045]  <- selected on DD
  H5a (secondary): R5 learned per band vs cap on test cohorts -> does not beat the cap (upper bound +0.045)
    test cell autism  PM1-like   PP3 Supporting only  n  148: error cap -0.15 pts, R5 learned per band +0.23 pts
    test cell autism  PM1-like   PP3 Moderate only    n  132: error cap -0.28 pts, R5 learned per band +0.18 pts
    test cell autism  PM1-like   PP3 >=Strong(3-4)    n  173: error cap +0.45 pts, R5 learned per band +0.20 pts
    test cell autism  PM1-strict PP3 >=Strong(3-4)    n   54: error cap +0.36 pts, R5 learned per band +0.11 pts
H5b PP2 beside PM1-like (error = ln R + (1 - kappa_PP2) ln LR(PP2), points):
  test autism: add PP2 (kappa 1)            error -1.16 points (SE 0.12)
  test autism: no PP2 beside PM1 (kappa 0)  error +0.62 points (SE 0.12)
  test autism: learned in DD (kappa 0.33)   error +0.02 points (SE 0.12)
H5c triple PP3 x PM1 x PP2 (descriptive; error in points, test = autism):
  autism Supporting x PM1-like  : cap + add PP2              error -1.97 points
  autism Supporting x PM1-like  : cap + no PP2               error -0.19 points
  autism Supporting x PM1-like  : Supporting limit + no PP2  error -0.19 points
  autism Supporting x PM1-like  : one level lower + no PP2   error +1.07 points
  autism Supporting x PM1-like  : uncapped                   error -1.97 points
  autism Strong     x PM1-like  : cap + add PP2              error -1.31 points
  autism Strong     x PM1-like  : cap + no PP2               error +0.47 points
  autism Strong     x PM1-like  : Supporting limit + no PP2  error +0.93 points
  autism Strong     x PM1-like  : one level lower + no PP2   error +0.47 points
  autism Strong     x PM1-like  : uncapped                   error -2.23 points
  autism Supporting x PM1-strict: cap + add PP2              error -2.09 points
  autism Supporting x PM1-strict: cap + no PP2               error -0.31 points
  autism Supporting x PM1-strict: Supporting limit + no PP2  error -0.31 points
  autism Supporting x PM1-strict: one level lower + no PP2   error +0.95 points
  autism Supporting x PM1-strict: uncapped                   error -2.09 points
  autism Strong     x PM1-strict: cap + add PP2              error -1.51 points
  autism Strong     x PM1-strict: cap + no PP2               error +0.27 points
  autism Strong     x PM1-strict: Supporting limit + no PP2  error +0.73 points
  autism Strong     x PM1-strict: one level lower + no PP2   error +0.27 points
  autism Strong     x PM1-strict: uncapped                   error -2.43 points
```

## Supplementary Note 13. Functional assays with curated PM1 in MaveDB (A21)

`make_mave_mavedb.py`; score-set status, estimates and negative control in `a21_mavedb_sets.tsv`, `a21_mavedb_estimates.tsv`, `a21_mavedb_negative_control.tsv`.

```
A21: score sets examined 81; status: too few nonsense/synonymous 54, mapping failed 15, informative 12
AlphaMissense Supporting x curated PM1: genes 5, score sets 12; random effects R 0.59 [0.44, 0.80]; genes with R < 1: 4/5
  leave-one-gene-out R range 0.53-0.71
  PRIMARY decision: dependence confirmed if upper bound < 1 -> CONFIRMED
    F9       R 104398.54 [3616.86, 3013406.68]
    GCK      R 0.49 [0.43, 0.54]
    LDLR     R 0.52 [0.44, 0.61]
    PTEN     R 0.57 [0.51, 0.62]
    TP53     R 0.61 [0.44, 0.84]
AlphaMissense Strong x curated PM1: genes 5, score sets 12; random effects R 0.48 [0.30, 0.76]; genes with R < 1: 4/5
PP3 share of log-evidence kept at curated residues: median 0.40 (score sets 12)
Negative control (random residue sets, same size): median R per gene F9 1.08, GCK 1.00, LDLR 0.99, PTEN 0.98, TP53 1.02; across genes median 1.00
Sensitivity (not pre-registered):
  score sets with inverted truth (LR(PP3) < 1, e.g. dominant-negative assays where nonsense behaves like wild type): urn:mavedb:00001235-a-1
  score sets with degenerate R (> 100; a joint cell without non-damaging variants): urn:mavedb:00001200-d-1, urn:mavedb:00001200-b-1, urn:mavedb:00001200-e-1
  excluding inverted-truth sets: genes 5, sets 11; random effects R 0.55 [0.40, 0.75]; F9 104398.54, GCK 0.49, LDLR 0.52, PTEN 0.57, TP53 0.42
  excluding inverted and degenerate sets: genes 5, sets 8; random effects R 0.52 [0.47, 0.57]; F9 2.88, GCK 0.49, LDLR 0.52, PTEN 0.57, TP53 0.42
```

## Supplementary Data files

| File | Content |
|---|---|
| Supplementary Data 1 (`a12_gclust2.tsv`) | All de novo dependence estimates (Table S3) |
| Supplementary Data 2 (`mave_all_sites_summary.tsv`) | All mapped functional assays |
| Supplementary Data 3 (`supp_table1_predictor_resource.tsv`) | Per-predictor retention with intervals |
| Supplementary Data 4 (`lab_stacking_summary.txt`, `stacking_cv_summary.txt`) | Expert-panel and laboratory outputs |
| Supplementary Data 5 (`mave_curated_sites.tsv`, `review20_lab_summary.txt`, `review21_lab_summary.txt`) | Curated-PM1 assays, within-submitter, frequency-handled, held-out and trend analyses |
| Supplementary Data 6 (`a19_*_summary.txt`, `a19_h4_eras.tsv`) | Pre-registered tests H1–H4 and negative controls |
| Supplementary Data 8 (`a20_denovo_rules_summary.txt`) | Rules on de novo truth |
| Supplementary Data 9 (`a21_mavedb_*.tsv`, `a21_mavedb_summary.txt`) | MaveDB functional assays with curated PM1 |
| Supplementary Data 7 (`cost_headline_summary.txt`, `validity_v14_summary.txt`, `parser_annotation_sheet_v14.tsv`) | Cost of stacking, validity checks, parser annotation sheet |
