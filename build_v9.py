"""Builds MANUSCRIPT_v9.md from draft 8.1, applying review 13. Run: python build_v9.py"""
s = open("MANUSCRIPT_v8.md", encoding="utf-8").read()


def rep(a, b):
    global s
    assert s.count(a) == 1, a[:80]
    s = s.replace(a, b)


rep("*Draft 8, 7 October 2026. Supplementary material: `SUPPLEMENT_v8.md`.*",
    "*Draft 9, 7 October 2026. Supplementary material: `SUPPLEMENT_v9.md`.*")
# ---------------------------------------------------------------- abstract (fewer numbers, nominal scale, honest status)
i = s.index("## Abstract")
j = s.index("---", i)
s = s[:i] + """## Abstract

**Purpose.** ACMG/AMP classification sums computational (PP3), hotspot/domain (PM1) and gene-level (PP2) evidence as if independent. ClinGen caps PP3 + PM1 at Strong. We measured, without variant-level labels, whether these criteria are independent and how expert panels combine them.

**Methods.** Likelihood ratios came from the excess of de novo missense variants over mutation-rate expectation in 31,058 developmental-disorder trios, 12,316 autism probands and 4,256 public trios. Dependence between criteria was measured on this scale and on 6,177 expert-panel missense classifications re-scored without computational or location criteria. The primary pre-specified retention criterion was not met; the expert-practice analyses were not pre-specified.

**Results.** Most predictor evidence concerned gene and position: a median 0.32 of log-evidence survived conditioning on 50-residue windows. No pair of criteria was independent (AlphaMissense PP3 × PM1-like R 0.44 [0.41–0.47]; PM1 × PP2 R 0.31 [0.28–0.33]). With Strong PP3, the three criteria nominally sum to 7 points (5 when capped) against 4.05 points of de novo evidence. At Supporting the sum, 4, approximately matched the evidence (3.65). Expert panels applied PP3 beside PM1 almost always below Strong (783 of 787), where the cap does not act. There, Moderate PP3 added no detectable evidence (0.02 points [−1.69–0.89]), whereas Supporting PP3 was approximately calibrated (0.81 [0.24–1.68]). This held in held-out genes and in later classifications.

**Conclusion.** Location-type criteria overlap, and the overlap that matters in expert practice lies above Supporting PP3, below the reach of the current cap. Limiting PP3 to Supporting when PM1 applies is consistent with these data.

""" + s[j:]
# ---------------------------------------------------------------- methods: status wording
rep("- **The benchmark, the oracle, the PM1 × PP2 test and the predictor-validation analyses** were requested in review and are descriptive.",
    "- **Not pre-specified:** the benchmark, the oracle, the PM1 × PP2 and triple tests, the predictor-validation analyses and the expert-practice audit were added after the main analyses. They are descriptive.")
rep("- **k is unknown**, so results are reported for k = 10–1,000. As a reference, k ≈ 30 corresponds to about 1% of all missense opportunity being DD-pathogenic.",
    "- **k is unknown**, and the de novo data do not identify it. Results are reported for k = 10–300. The oracle is used only to show that retention is not a measure of skill, not to place predictors on a skill scale.")
rep("""**Alternative rules.** Every classification was re-scored on the points scale under four rules:
1. the ClinGen cap;
2. a cap on all three location-type criteria;
3. PP2 not added when PM1 applies;
4. PP3 credited at its measured value beside PM1.""",
"""**Alternative rules.** Every classification was re-scored on the points scale under five rules:
1. the ClinGen cap;
2. a cap on all three location-type criteria;
3. PP2 not added when PM1 applies;
4. PP3 limited to Supporting when PM1 applies;
5. PP3 credited at its measured value beside PM1.

**Robustness and out-of-sample tests.**
- **Truth thresholds:** P/LP at ≥ 4, 6 or 8 points; B/LB at ≤ −1 or −3.
- **Frequency criteria** (BA1, BS1, BS2, PM2) also removed from the truth points.
- **Gene-split cross-validation:** 200 random halves of the genes.
- **Temporal split:** credit estimated on classifications approved before 2025 and tested on those approved in 2025–2026.""")
# ---------------------------------------------------------------- results: figure call-outs, oracle, benchmark framing
rep("**Retention.** Across the 39 informative predictors in DD,",
    "**Retention** (design in Figure 1). Across the 39 informative predictors in DD,")
rep("At k ≈ 30 a perfect predictor retains 0.54. Real predictors (median 0.32; AlphaMissense 0.30 [0.27, 0.33]; best, VEST4, 0.59) therefore keep roughly 60% of the attainable within-window evidence. Retention is not a measure of skill on its own scale: it measures evidence incremental to location, which is the quantity that matters when PP3 is combined with PM1 or PP2.",
    "Across plausible k, a perfect predictor retains 0.30–0.90, a range that contains the observed median (0.32; AlphaMissense 0.30 [0.27, 0.33]; VEST4 0.59). Because k is not identified, retention cannot be read as a fraction of attainable skill. It measures evidence incremental to location, which is the quantity that matters when PP3 is combined with PM1 or PP2.")
rep("### A gene- and position-adjusted benchmark\n**Figure 2** and **Table 3** compare pooled and within-window rankings.",
    "### A gene- and position-adjusted benchmark (resource)\n**Figure 2** and **Table 3** compare pooled and within-window rankings. We present the benchmark as a mechanism and a resource; the re-ranking is not a headline finding.")
rep("**Pooled and within-window rankings differed.**", "**Where the rankings differ, and where they do not.**")
rep("- Only VEST4's within-window LR had a lower bound above the ClinGen Supporting LR (2.406).",
    "- Only VEST4's within-window LR had a lower bound above the ClinGen Supporting LR (2.406). Below the top two, rank intervals span 5–19 places.")
rep("""**Within-window evidence tracks independent variant-level performance.** Across the 39 predictors, within-window LR correlated with within-gene ClinVar AUROC (ρ 0.55 [0.22, 0.77]; 0.56 in DDG2P genes) and with assay performance (ρ 0.63 [0.39, 0.79]). Marginal LR correlated less strongly with the same measures (ρ 0.37 for each). With 39 correlated predictors the differences were not significant: +0.18 [−0.32, 0.63] for ClinVar and +0.26 [−0.11, 0.65] for assays.""",
    """**Within-window evidence tracks independent variant-level performance, but not demonstrably better than pooled evidence.**
- Across the 39 predictors, within-window LR correlated with within-gene ClinVar AUROC (ρ 0.55 [0.22, 0.77]; 0.56 in DDG2P genes) and with assay performance (ρ 0.63 [0.39, 0.79]).
- Marginal LR also correlated with them (ρ 0.37 for each).
- The differences were not significant: +0.18 [−0.32, 0.63] for ClinVar and +0.26 [−0.11, 0.65] for assays.

Where PP3 is applied (within disease genes), pooled and within-window rankings largely agree.""")
rep("### PP3, PM1 and PP2 are not independent\n**PP3 and PM1.**",
    "### PP3, PM1 and PP2 are not independent\n**PP3 and PM1** (Figure 3).")
rep("- Across all 40 predictors informative in this analysis, R × PM1-like had an upper bound below 1 in DD (median 0.51). Predictors using population constraint or allele frequency were no more dependent than others (p = 0.87).",
    "- Across the 40 predictors informative on all analysed sites, R × PM1-like had an upper bound below 1 in DD (median 0.51). This set includes popEVE, phyloP470way and fathmm-MKL; the benchmark's site set (constraint annotated) instead admits SIFT and PolyPhen2_HVAR, and 36 predictors are in both. Predictors using population constraint or allele frequency were no more dependent than others (p = 0.87).")
rep("- PM1 and PP2 were the most dependent pair: R 0.31 [0.28, 0.33] in DD and 0.36 [0.30, 0.45] in autism for PM1-like, and 0.27 and 0.34 for PM1-strict. Constrained regions lie mostly in constrained genes, so applying both counts the same constraint twice.",
    "- PM1 and PP2 were the most dependent pair: R 0.31 [0.28, 0.33] in DD and 0.36 [0.30, 0.45] in autism for PM1-like, and 0.27 and 0.34 for PM1-strict.\n- Part of this is by construction. Both of our flags are built from missense constraint, which ACGS assigns to PP2 rather than PM1 (Durkie et al. 2024). Panel-defined PM1 rests also on hotspots and functional domains, so its dependence on PP2 may be smaller.")
# ---------------------------------------------------------------- stacking section: honest numbers, sensitivity, temporal
rep("  - **Moderate PP3** (67 classifications) earned 0.02 points [−1.69, 0.89] against 2.00 applied: no detectable independent evidence. The ClinGen cap allows Moderate PP3 with Moderate PM1 (4 points), so it never intervenes here.",
    """  - **Moderate PP3** (67 classifications; 6 B/LB in the joint cell) earned 0.02 points [−1.69, 0.89] against 2.00 applied. The interval excludes the points applied, though it is wide. The ClinGen cap allows Moderate PP3 with Moderate PM1 (4 points), so it never intervenes here.
- **Robustness.** The credit was stable across truth thresholds: 0.60–0.65 overall, 0.81–0.89 for Supporting PP3 and −0.23 to 0.14 for Moderate. With frequency criteria also removed it was 0.49 [−0.06, 1.31], on only 6 B/LB variants in the joint cell (Supplementary Note 11).
- **Truth that conditions on other evidence.** Truth requires the remaining evidence to reach ≥ 6 or ≤ −1 points, and panels choose PP3 strength knowing that evidence. Both can bias the credit; neither the threshold analysis nor the out-of-sample tests suggest a large effect.""")
rep("- **Out of sample.** Across 200 random splits of the genes, the credit estimated on one half was 0.61 (median) and on the held-out half 0.59. The held-out credit fell below the points applied in 98% of splits.",
    """- **Out of sample.**
  - **Held-out genes.** Across 200 random splits of the genes, the credit estimated on one half was 0.61 (median) and on the held-out half 0.59. The held-out credit fell below the points applied in 98% of splits.
  - **Later classifications.** Estimated on classifications approved before 2025, the credit was 0.70 [−0.14, 1.26]. On classifications approved in 2025–2026 it was 0.40 [−0.25, 1.07], against 1.40 points applied.""")
rep("""**The triple stack on the de novo scale.**
- **Supporting.** When AlphaMissense PP3, PM1-like and PP2-like all apply, their joint de novo LR in DD was 24.6 (3.65 points [3.50, 3.79]), against a product of 196.7 (6.02 points): R3 0.13 [0.11, 0.14].
- **Strong.** At Strong PP3 the joint evidence was 4.05 points against 6.75 summed (R3 0.09).
- **Autism.** The pattern replicated: R3 0.18 [0.13, 0.24] at Supporting and 0.14 at Strong.
- **PM1-strict.** With PM1-strict, R3 was 0.10 and 0.07.""",
"""**The triple stack on the de novo scale.** When AlphaMissense PP3, PM1-like and PP2-like all apply:

- **Dependence.** Their joint de novo evidence was far below the product of their separate de novo LRs (R3 0.13 [0.11, 0.14] at Supporting, 0.09 at Strong; autism 0.18 and 0.14).
- **Against the nominal ACMG weights** (PP3 1 or 4, PM1 2, PP2 1):

| PP3 tier | Nominal sum | Under the cap | De novo evidence, DD | De novo evidence, autism |
|---|---|---|---|---|
| Supporting | 4 | 4 | 3.65 [3.50, 3.79] | 2.96 [2.67, 3.22] |
| Strong | 7 | 5 | 4.05 [3.89, 4.20] | 3.27 [2.96, 3.53] |

At Supporting, the nominal weights approximately match the joint evidence. With Strong PP3 they exceed it by about 3 points uncapped and 1 point capped. The de novo scale is not the clinical scale, so these are comparisons of direction and size, not calibrated points.""")
rep("Applied to held-out genes, with the credit estimated on the other half, the measured-credit rule changed a median of 89 of 3,103 classifications per split (95% range 13–140), about 182 when scaled to the full set. Limiting PP3 to Supporting needs no estimated parameter.",
    """**Out of sample.**
- **Held-out genes.** With the credit estimated on the other half of the genes, the measured-credit rule changed a median of 89 of 3,103 classifications per split (2.9%; 95% range 13–140), equivalent to about 177 of 6,177.
- **Later classifications.** Among the 2,507 missense classifications approved in 2025–2026, limiting PP3 to Supporting would change 15 and the cap 16.

Limiting PP3 to Supporting needs no estimated parameter.""")
# ---------------------------------------------------------------- clinical consequences: remove stale A12.5 paragraph
rep("**Expert-panel classifications.** Among 13,293, 1,018 applied both PP3 and PM1. The cap would change 24 of them and the conditional rule 238.\n\n",
    "**Expert-panel classifications.** See *Stacking in expert practice*. The cap changes 24 classifications whether all 13,293 classifications or only the 6,177 missense ones are considered. Of the 72 classifications in which PP3 + PM1 exceed 4 points, 70 are missense, and the 2 that are not do not change class.\n\n")
# ---------------------------------------------------------------- discussion
rep("""- **About a third of the log-evidence is incremental.** A perfect predictor would keep more, but not all, so about two-thirds of what pooled evaluations credit to predictors restates location.
- **This matters where location is credited separately**, as ACMG/AMP does through PP2 and PM1.
- **It is not a statement about skill on an absolute scale.** The oracle shows that the attainable within-window evidence depends on how enriched pathogenic sites are. At plausible enrichment, current predictors capture roughly 60% of it.""",
"""- **About a third of the log-evidence is incremental to gene and window.** That is the quantity that matters where location is credited separately, as ACMG/AMP does through PP2 and PM1.
- **It is not a statement about skill.** A perfect predictor would also lose evidence within windows, by an amount that depends on an enrichment parameter the data do not identify.""")
rep("""- **The ranking changes.** Once gene and window are fixed, the ranking is unrelated to the pooled one, and the predictors leading pooled comparisons fall furthest. MPC, built from regional constraint for exactly this cohort type, falls to last.
- **Within known disease genes the rankings agree.** The divergence is about between-gene information, not a different notion of variant effect.
- **The adjusted ranking is supported by independent truth.** It tracks within-gene ClinVar and assay performance, as Harrizi et al. (2026) found for within-gene ClinVar evaluation arbitrated by assays. We could not show it does so better than the pooled ranking.
- **Our recommendation is modest:** when PP3 is combined with gene- or region-level criteria, within-window evidence is the relevant benchmark.""",
"""- **Across all genes the rankings differ**, consistent with Harrizi et al. (2026) on ClinVar. The predictors leading pooled de novo comparisons tend to encode constraint; MPC, built from regional constraint, falls to last.
- **Within known disease genes the rankings agree**, and within-window evidence tracks within-gene clinical and functional performance no better than pooled evidence does, as far as these data can tell.
- **We therefore offer the benchmark as a mechanism and a resource**, not as evidence that a different predictor should be used for PP3.""")
rep("- **Expert practice.** Panels stack these criteria, mostly at strengths the existing cap does not reach. On the panels' own classifications, PP3 beside PM1 was worth about half the points applied. When PP3, PM1 and PP2 all apply, the summed points overstate the de novo evidence by about 2.4 points, enough to move a variant across a class boundary.",
    "- **Expert practice.** Panels stack these criteria, mostly at strengths the existing cap does not reach.\n  - On the panels' own classifications, Supporting PP3 beside PM1 was approximately calibrated, while Moderate PP3 added no detectable evidence.\n  - On the de novo scale the nominal weights match the joint evidence when PP3 is Supporting, and exceed it when PP3 is stronger.\n  - The two scales agree that the overlap that matters lies above Supporting PP3.")
rep("""  - When PM1 applies, count PP3 at most as Supporting. Moderate PP3 beside PM1 added no detectable evidence, and the existing cap does not prevent this combination.
  - Do not sum PM1 with PP2 unadjusted. PM1 × PP2 was the most dependent pair on the de novo scale, and seven expert panels currently stack them.
  - Treat PP3 + PM1 + PP2 as one location-type evidence budget. Panels that wish to keep co-application should calibrate the combination jointly against independent truth, as done here.""",
"""  - **Limit PP3 to Supporting when PM1 applies.** This is consistent with the data and needs no estimated parameter. The supporting evidence rests on 67 Moderate-PP3 classifications, so we suggest panels confirm it on their own data.
  - **Check PM1 + PP2 for shared evidence**, as ACGS advises. Our flags quantify a strong dependence, partly by construction, and expert truth was too sparse in benign variants to test panel-defined PM1 + PP2.
  - **Calibrate co-application jointly.** Panels that keep co-applying location-type criteria should calibrate the combination against independent truth, as done here.""")
rep("UK guidance already cautions against this qualitatively (ACGS 2024); our data quantify it, and show that seven ClinGen panels still co-apply the two.",
    "UK guidance already advises against counting the same constraint evidence in both (Durkie et al. 2024). Our data quantify the dependence for constraint-based flags, and show that seven ClinGen panels co-apply the two.")
rep("**Conclusion.** Most of what missense predictors contribute restates where a variant lies. The criteria that credit location separately are not independent: PM1 and PP2 most of all, and PP3 and PM1 even within disease genes. Stacked together, the three overstate the de novo evidence by about 2.4 points. Expert panels stack them mainly at strengths the existing cap does not reach, and crediting PP3 at its measured value beside PM1 would change up to 187 expert classifications. Location-type evidence should be budgeted jointly rather than summed.",
    "**Conclusion.** Most of what missense predictors contribute concerns where a variant lies, and the criteria that credit location separately are not independent. In expert practice, PP3 is stacked beside PM1 mainly at strengths the existing cap does not reach. There, Supporting PP3 appears calibrated and Moderate PP3 does not. Limiting PP3 to Supporting when PM1 applies is a simple rule consistent with these data.")
# ---------------------------------------------------------------- references
rep("- ACGS. Best practice guidelines for variant classification in rare disease 2024. Association for Clinical Genomic Science, 2024.",
    "- Durkie M, Cassidy E-J, Berry I, et al. ACGS best practice guidelines for variant classification in rare disease 2024. Association for Clinical Genomic Science, 2024.")
s = s.replace("(ACGS 2024)", "(Durkie et al. 2024)")
s = s.replace("SUPPLEMENT_v8.md", "SUPPLEMENT_v9.md")
open("MANUSCRIPT_v9.md", "w", encoding="utf-8").write(s)
sup = open("SUPPLEMENT_v8.md", encoding="utf-8").read().replace(
    "*Draft 8, 7 October 2026. Accompanies `MANUSCRIPT_v8.md`.", "*Draft 9, 7 October 2026. Accompanies `MANUSCRIPT_v9.md`.")
cv = open("results_fullstudy_backup/stacking_cv_summary.txt", encoding="utf-8").read().strip()
i = sup.index("**Strength-specific credit and gene-split cross-validation**")
j = sup.index("**All variant types**")
sup = sup[:i] + "**Strength-specific credit, sensitivity analyses, gene-split cross-validation and temporal split** (`make_stacking_cv.py`):\n\n```\n" + cv + "\n```\n\n" + sup[j:]
open("SUPPLEMENT_v9.md", "w", encoding="utf-8").write(sup)
print("ok")
