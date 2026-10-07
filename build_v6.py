"""Builds MANUSCRIPT_v6.md / SUPPLEMENT_v6.md from draft 5: benchmark-led reframing (review 9). Run: python build_v6.py"""
import re
import pandas as pd

s = open("MANUSCRIPT_v5.md", encoding="utf-8").read()
b = pd.read_csv("supp_table7_benchmark.tsv", sep="\t", index_col=0)
for o, n in (("4", "5"), ("3", "4"), ("2", "3")):          # main-text tables shift down by one
    s = re.sub(r"(?<!Supplementary )Table " + o + r"\b", "Table " + n, s)


def rep(a, new):
    global s
    assert s.count(a) == 1, a[:70]
    s = s.replace(a, new)


rep("# Computational (PP3) and hotspot/domain (PM1) evidence overlap for de novo missense variants in neurodevelopmental "
    "disorders: a label-free test of additivity",
    "# De novo benchmarks overstate missense predictor skill: a gene- and position-adjusted evaluation and the double "
    "counting of computational and hotspot evidence (PP3 + PM1)")
rep("*Draft 5, 7 October 2026. Supersedes draft 4 (`MANUSCRIPT_v4.md`). Supplementary material: `SUPPLEMENT_v5.md`.",
    "*Draft 6, 7 October 2026. Supersedes draft 5 (`MANUSCRIPT_v5.md`). Supplementary material: `SUPPLEMENT_v6.md`.")

ABSTRACT = """## Abstract

**Purpose.** Missense predictors are benchmarked on de novo and clinical data that pool variants across genes. We asked how much of their evidence concerns variant effect within a gene, and what this implies for combining computational (PP3) with hotspot/domain (PM1) evidence under ACMG/AMP.

**Methods.** Likelihood ratios (LRs) were estimated without variant-level labels. They came from the excess of de novo missense variants over mutation-rate expectation in 31,058 developmental-disorder (DD) trios, 12,316 autism probands and 4,256 public trios. LRs were computed marginally, within genes and within 50-residue windows. Additivity of PP3 and PM1-like evidence was measured as R, the joint LR over the product of the separate LRs, and checked on ClinVar and expert-panel truth.

**Results.** Across 39 informative predictors, a median 0.32 of log-evidence survived conditioning on gene and window, and pooling inflated LRs 3.6-fold. Marginal and within-window rankings were unrelated (Spearman −0.18). The four strongest marginal predictors fell to ranks 17–39. The within-window ranking replicated in autism (ρ 0.85). PP3 and PM1-like evidence were therefore not additive. For AlphaMissense, R was 0.44 [0.41–0.47] in DD, 0.56 [0.47–0.66] in autism and 0.67 [0.48–0.96] in the public cohort, and 0.69–0.75 within DDG2P genes. On expert-panel truth stripped of computational and position criteria, Strong PP3 earned 1.83 points [0.84–3.54] once PM1 applied, against 4 under uncapped co-use. A pre-registered label-free correction did not beat ClinGen's cap.

**Conclusion.** De novo case–control benchmarks overstate predictors' variant-level skill. Gene- and position-adjusted evaluation re-ranks predictors and explains why PP3 double-counts PM1. The evidence is consistent with ClinGen's cap and excludes uncapped co-use.

"""
i = s.index("## Abstract")
s = s[:i] + ABSTRACT + s[s.index("---", i):]

rep("""De novo variants avoid these problems. Their excess over mutation-rate expectation in affected children defines a likelihood ratio without variant-level labels, and that ratio can be computed within a gene, a window or a domain.

We used the de novo excess to answer three questions:

1. How much predictor evidence survives conditioning on gene and protein position?
2. Is PP3 additive with PM1-like evidence, and does the answer hold on independent truth?
3. Does a label-free alternative calibrate better than the cap?""",
    """De novo variants avoid these problems, but published de novo evaluations do not exploit this. AlphaMissense, PrimateAI-3D, popEVE and MVP were all assessed by separating case from control de novos pooled across genes (Qi et al. 2021; Cheng et al. 2023; Gao et al. 2023; Orenbuch et al. 2025). Such comparisons reward a predictor for knowing which genes and regions are intolerant, which is information ACMG/AMP already credits through PP2 and PM1. The de novo excess over mutation-rate expectation defines a likelihood ratio without variant-level labels, and that ratio can be computed within a gene or a window.

We used it to answer four questions:

1. How much predictor evidence survives conditioning on gene and protein position?
2. Does a gene- and position-adjusted benchmark rank predictors differently?
3. Is PP3 additive with PM1-like evidence, and does the answer hold on independent truth?
4. Does a label-free alternative calibrate better than ClinGen's cap?""")

rep("### Direct overlap test\n**Flags at each site:**",
    """### Gene- and position-adjusted benchmark
Each informative predictor was ranked by its within-window LR for its top 2% of scores. This is the evidence it adds once gene and 50-residue window are fixed. The ranking was compared with the marginal ranking, and within-window LRs were compared across cohorts. Predictors from dbNSFP 5.4a were evaluated on their own site set, and intervals are from the individual-level bootstrap. The full resource is Supplementary Table 7 (`make_benchmark.py`).

### Direct overlap test
**Flags at each site:**""")

top = b.head(10)
big = b.sort_values("DD_marginal_LR", ascending=False).head(4)
rows = pd.concat([top, big[~big.index.isin(top.index)]])
L = ["| Predictor | Marginal LR (rank) | Within gene | Within window [95% CI] (rank) | Share retained |", "|---|---|---|---|---|"]
for k, r in rows.iterrows():
    L.append(f"| {k} | {r.DD_marginal_LR:.1f} ({int(r.rank_marginal)}) | {r.DD_within_gene_LR:.2f} | "
             f"{r.DD_within_window_LR:.2f} [{r.D2_lo:.2f}, {r.D2_hi:.2f}] ({int(r.rank_within_window)}) | {r.r_D2:.2f} |")
TAB = "\n".join(L)
rep("### PP3 and PM1-like evidence are not additive\n", f"""### A gene- and position-adjusted benchmark re-ranks predictors
**Table 2** ranks the 39 informative predictors (25 from dbNSFP 4.1a plus AlphaMissense, and 14 newer) by within-window LR. The ten best are shown, with the four strongest marginal predictors appended (all 39, with replications, in Supplementary Table 7).

**Table 2. Predictors ranked by evidence within 50-residue windows** (DD; top 2% of scores; ranks among 39).

{TAB}

**Inflation.**
- Pooling across genes inflated the top-2% LR a median 3.6-fold (range 1.5–10.6).
- The median LR fell from 6.08 (marginal) to 2.99 within gene and 1.78 within window.

**The rankings were unrelated** (Spearman −0.18, p = 0.28).
- The four strongest marginal predictors, MisFit_S, MutScore, PHACTboost and MPC (marginal LR 13.6–16.2), ranked 17th, 21st, 26th and 39th within windows.
- VEST4, MutPred, MetaRNN and ClinPred, ranked 20th–32nd marginally, led within windows.
- popEVE, evaluated in DDD only, fell from 17.9 to 1.76.
- AlphaMissense was the only predictor in the top five on both scales.
- Only VEST4's within-window LR had a lower 95% bound above the ClinGen Supporting LR (2.406).

**Replication.** Within-window LRs correlated across predictors between DD and autism (Spearman 0.85, 25 predictors) and, more weakly, with the smaller public cohort (0.49).

The predictors that look strongest in pooled de novo evaluations are largely those that encode gene and regional constraint. That is the information that also underlies PM1.

### PP3 and PM1-like evidence are not additive
""")

rep("""**Main finding.** For de novo missense variants in NDD, most of what missense predictors contribute is information about which gene and which part of the protein a variant falls in. A predictor that scores a hotspot highly is not wrong; it restates what location already says.""",
    """**Main finding.** De novo case–control benchmarks reward missense predictors for knowing which genes and regions matter. Once gene and position are held fixed:

- about two-thirds of the evidence disappears;
- the ranking of predictors changes completely;
- the strongest predictors in pooled evaluations, several of them built partly from constraint, add among the least within a window.

A predictor that scores a hotspot highly is not wrong; it restates what location already says. Pooled de novo benchmarks therefore overstate variant-level skill, and the gene- and position-adjusted ranking is the relevant one wherever PP3 is applied alongside gene- or region-level criteria. The same redundancy is the mechanism of PP3 + PM1 double counting.""")
rep("- **What is new here:**",
    """- **Predictor benchmarks.** Published de novo evaluations of AlphaMissense, PrimateAI-3D, popEVE and MVP compare cases with controls pooled across genes. Gene-level heterogeneity and circularity of predictor performance are documented on ClinVar (Grimm et al. 2015; Tejura et al. 2024; Harrizi et al. 2026). To our knowledge, no previous benchmark evaluates predictors within gene and window on label-free de novo data.
- **What is new here:**""")
rep("  - a label-free measurement of PP3 + PM1 additivity on the points scale;",
    "  - a label-free, gene- and position-adjusted benchmark of 39 predictors, released as a resource;\n"
    "  - a label-free measurement of PP3 + PM1 additivity on the points scale;")
rep("7. **Ancestry.** The cohorts are predominantly of European ancestry.",
    "7. **Ancestry.** The cohorts are predominantly of European ancestry.\n"
    "8. **Benchmark scope.** The benchmark scores each predictor's top 2% of scores against de novo excess in severe developmental "
    "disorders. Predictors that are informative mainly at lower thresholds, or for other phenotypes, may rank differently.")
rep("- **Table 1.** Pre-specified retention tests.",
    "- **Table 1.** Pre-specified retention tests.\n- **Table 2.** Predictors ranked by within-window evidence.")
rep("  - Table 6: hierarchical model calibration table (`results_fullstudy/a16_calibration_table.tsv`).",
    "  - Table 6: hierarchical model calibration table (`results_fullstudy/a16_calibration_table.tsv`);\n"
    "  - Table 7: gene- and position-adjusted predictor benchmark, all 39 predictors (`supp_table7_benchmark.tsv`).")
rep("- Garrett A, et al.",
    "- Gao H, et al. The landscape of tolerated genetic variation in humans and primates (PrimateAI-3D). *Science* 2023.\n- Garrett A, et al.")
rep("- Richards S, et al. *Genet Med* 2015.",
    "- Qi H, et al. MVP predicts the pathogenicity of missense variants by deep learning. *Nat Commun* 2021;12:510.\n"
    "- Richards S, et al. *Genet Med* 2015.")
s = s.replace("SUPPLEMENT_v5.md", "SUPPLEMENT_v6.md")
open("MANUSCRIPT_v6.md", "w", encoding="utf-8").write(s)
sup = open("SUPPLEMENT_v5.md", encoding="utf-8").read()
sup = sup.replace("*Draft 5, 7 October 2026. Accompanies `MANUSCRIPT_v5.md`.", "*Draft 6, 7 October 2026. Accompanies `MANUSCRIPT_v6.md`.")
open("SUPPLEMENT_v6.md", "w", encoding="utf-8").write(sup)
print("ok")
