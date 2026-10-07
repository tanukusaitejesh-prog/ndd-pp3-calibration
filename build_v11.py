"""Builds MANUSCRIPT_v11.md from draft 10.1 (review 17). Run: python build_v11.py"""
s = open("MANUSCRIPT_v10.md", encoding="utf-8").read()


def rep(a, b):
    global s
    assert s.count(a) == 1, a[:80]
    s = s.replace(a, b)


rep("# Computational and location evidence are not independent in missense variant classification, and current practice over-credits their combination",
    "# Location-type evidence is double counted in missense variant classification: PP2 adds nothing beside PM1, and PP3 above Supporting adds little")
rep("*Draft 10.1, 7 October 2026. Supplementary material: `SUPPLEMENT_v10.md`.*",
    "*Draft 11, 7 October 2026. Supplementary material: `SUPPLEMENT_v11.md`.*")
i = s.index("## Abstract"); j = s.index("---", i)
s = s[:i] + """## Abstract

**Purpose.** ACMG/AMP classification sums computational evidence (PP3) with location evidence: hotspot or domain (PM1) and gene-level (PP2). It treats them as independent, and ClinGen caps PP3 + PM1 only at Strong. We measured how much each criterion adds once another applies.

**Methods.** Likelihood ratios came from the excess of de novo missense variants over mutation-rate expectation in 47,630 trios and probands. We estimated the points each criterion earns beside PM1 on 6,177 expert-panel and 87,873 laboratory missense classifications, re-scored without computational or location criteria.

**Results.**
- **No pair of criteria was independent** on the de novo scale.
- **PP2 adds nothing beside PM1 in laboratory practice:** it earned −0.20 points [−0.66–0.42] against 1 applied. Laboratories stacked the two in 1,511 submissions, and not adding PP2 would change 357 classifications.
- **Supporting PP3 beside PM1 was close to its nominal point** on panel and laboratory truth (0.81, 0.83).
- **Above Supporting, PP3 beside PM1 earned well below its nominal points.** On the de novo scale, Moderate AlphaMissense earned −0.12 [−0.38–0.10] against 2 and Strong 1.10 [1.04–1.16] against 4. In expert panels, Moderate PP3 earned 0.02 [−1.69–0.89].
- **The cap acts only at Strong,** which accounts for 1–2% of PP3 + PM1 co-application.

**Conclusion.** Location-type criteria are double counted. Adding PP2 beside PM1 is not supported, and PP3 above Supporting adds little beside PM1, where the existing cap rarely acts.

""" + s[j:]
rep("We used the de novo excess to measure the dependence between PP3, PM1 and PP2, checked it on four independent kinds of truth, and audited how expert panels and clinical laboratories combine these criteria.",
    "We used the de novo excess to measure the dependence between PP3, PM1 and PP2, and checked it on four independent kinds of truth. We then measured how many points PP2 and PP3 earn once PM1 applies, in expert-panel and laboratory practice, and compared these with the points curators give them.")
rep("**Credit.** The credit PP3 earns once PM1 applies is ln[LR(PP3 ∧ PM1)/LR(PM1)]/ln 2.406, on the ClinGen points scale. It was computed by PP3 strength, with gene-clustered intervals, and Moderate and Supporting PP3 were compared directly (difference test and its minimum detectable difference).",
    """**Credit.**
- **Definition.** The points a criterion C earns once PM1 applies are ln[LR(C ∧ PM1)/LR(PM1)]/ln 2.406, on the ClinGen points scale. We report this credit rather than R because it is what classification uses. R can be well below 1 while the credit equals the nominal points, when C alone is worth more than its nominal points on that truth.
- **Where it was computed:**
  - PP2 and PP3 on panel and laboratory truth, with PP3 by applied strength and, for panels, by predictor score band;
  - PP3 on the de novo scale, by AlphaMissense and REVEL score band (Supporting, Moderate, Strong), with gene-clustered intervals.
- **Moderate against Supporting.** Panel-applied Moderate and Supporting PP3 were compared directly, with a difference test and its minimum detectable difference.""")
rep("| ClinVar P/LP vs B/LB | AlphaMissense; PM1-like | DDG2P genes | 0.32 [0.28, 0.37] |",
    "| ClinVar P/LP vs B/LB* | AlphaMissense; PM1-like | DDG2P genes | 0.32 [0.28, 0.37] |")
rep("| Functional assays | AlphaMissense; PM1-like | 13 assays, 10 genes outside DDG2P | 0.94 [0.88, 0.99] |",
    """| Functional assays | AlphaMissense; PM1-like | 13 assays, 10 genes outside DDG2P | 0.94 [0.88, 0.99] |
| Laboratory truth: PM1 × PP2 | Laboratory-applied | 39,075 submissions | 0.28 [0.19, 0.44] |
| Laboratory truth: PP3 × PP2 | Laboratory-applied | 39,075 submissions | 0.42 [0.35, 0.51] |

\\*ClinVar classifications partly use PP3 and PM1 themselves, so this row is partly circular. The bias runs towards R = 1, because curators who apply both enrich P/LP variants in the joint cell.""")
# ---------------------------------------------------------------- replace practice results
i = s.index("### How panels and laboratories stack the criteria")
j = s.index("---\n\n## Discussion")
s = s[:i] + """### How panels and laboratories stack the criteria
**Table 3. Co-application of PP3, PM1 and PP2.**

| | Expert panels | Laboratories |
|---|---|---|
| Missense classifications | 6,177 | 87,873 |
| PP3 beside PM1 | 787 | 4,477 |
| — PP3 Supporting / Moderate / Strong | 605 / 178 / 4 | 4,212 / 185 / 80 |
| Share of PP3 + PM1 at Strong (where the cap acts) | 0.5% | 1.8% |
| PP2 beside PM1 | 174 (7 panels) | 1,511 |
| All three | 138 | 1,059 |

### PP2 adds nothing beside PM1
- **Laboratory truth.** PM1 and PP2 were strongly dependent: R 0.28 [0.19, 0.44]. Once PM1 applied, PP2 earned −0.20 points [−0.66, 0.42] against the 1 point applied (571 submissions; 57 B/LB).
- **De novo scale.** The constraint-based flags were similarly dependent (R 0.31 [0.28, 0.33]); part of this is by construction.
- **Expert-panel truth** had only one B/LB variant in the joint cell, so it could not test the pair.

UK guidance already advises against counting the same constraint evidence in both criteria (Durkie et al. 2024). These data show that in laboratory practice PP2 beside PM1 carries no detectable independent evidence.

### PP3 beside PM1, by strength
**Table 4. Points PP3 earns once PM1 applies** (95% intervals, gene-clustered; nominal points: Supporting 1, Moderate 2, Strong 4).

| Source of truth and definition of PP3 | Supporting | Moderate | Strong |
|---|---|---|---|
| **De novo excess, DD** — AlphaMissense band | −0.92 [−1.37, −0.59] | −0.12 [−0.38, 0.10] | 1.10 [1.04, 1.16] |
| **De novo excess, DD** — REVEL band | 0.47 [0.27, 0.64] | 1.04 [0.91, 1.16] | 1.82 [1.61, 2.01] |
| **De novo excess, autism** — AlphaMissense band | −0.04 [−0.72, 0.45] | −0.13 [−0.94, 0.35] | 0.90 [0.75, 1.04] |
| **De novo excess, autism** — REVEL band | 0.39 [−0.09, 0.71] | 0.86 [0.58, 1.08] | 1.39 [0.86, 1.76] |
| **Expert panels** — strength applied by the panel | 0.81 [0.24, 1.68] (n 235) | 0.02 [−1.69, 0.89] (n 67) | n 4 |
| **Expert panels** — REVEL band | −0.52 [−1.46, 1.17] | 0.58 [0.08, 1.46] | 1.84 [0.81, —] |
| **Expert panels** — AlphaMissense band | 1.63 [0.09, —] | 1.44 [0.34, —] | 1.75 [0.76, 3.74] |
| **Laboratories** — strength applied | 0.83 [0.53, 1.24] (n 1,681) | 0.75 [—] (3 B/LB) | 0 B/LB |

"—" marks an unbounded interval (≤ 3 B/LB in the joint cell). De novo bands are relative to all other PM1-like sites, and the de novo scale is not the clinical scale.

**Two independent lines of evidence.** The de novo excess and the expert panels are independent of each other; the panel rows share classifications.
- **On the de novo scale, every estimate above Supporting is below its nominal points, with intervals excluding them,** for both predictors and in both cohorts. Moderate AlphaMissense beside PM1 earns nothing (DD −0.12 [−0.38, 0.10]). Strong earns about 1 point (1.10 [1.04, 1.16]) where 4 are given.
- **In expert panels, Moderate PP3 applied beside PM1 earned at most 0.89 points** against 2. Where the score-band estimates are bounded, they also fall below their nominal points: REVEL Moderate below 2 and AlphaMissense Strong below 4.
- **Not every definition agrees.** AlphaMissense bands on panel truth were imprecise and their point estimates were at or above nominal at Supporting and Moderate. Defining PP3 by AlphaMissense eligibility alone, it earned 1.66 points [0.98, 2.58] beside panel PM1. REVEL eligibility gave 0.68 [0.34, 1.19], and either predictor 0.71 [0.35, 1.20] (Supplementary Note 14).
- **Moderate against Supporting.** In panel-applied data, Moderate could not be distinguished from Supporting (difference −0.79 [−2.81, 0.27]; p = 0.14; minimum detectable difference about 4 points).

**Supporting PP3.**
- **On clinical truth it is calibrated:** 0.81 in panels and 0.83 in laboratories, both intervals including 1.
- **On the de novo scale it earned less than 1** (REVEL 0.47; AlphaMissense below 0). Because the scales differ, we do not conclude that Supporting PP3 is over-credited in practice.
- **Laboratories overall** did not over-credit PP3 beside PM1 (0.85 [0.53, 1.20] against 1.09 applied), because 94% of their PP3 + PM1 co-applications are at Supporting.

**Why R < 1 does not imply over-credit at Supporting.** On panel truth PP3 alone has LR 8.2, more than its nominal Supporting LR of 2.4. Beside PM1 that falls to about the nominal value. The dependence removes the excess, not the nominal point.

**Expert-panel robustness of the overall credit** (0.60 [0.14, 1.05] against 1.24 applied; 17 B/LB in the joint cell; Figure 2):
- held-out genes: 0.59 on held-out halves, below the points applied in 98% of 200 splits;
- classifications approved in 2025–2026: 0.40 [−0.25, 1.07] against 1.40;
- de novo criteria removed from truth: 0.58 [0.08, 1.04];
- truth thresholds: 0.60–0.65;
- frequency criteria removed: 0.49 [−0.06, 1.31] (Supplementary Notes 11 and 14).

### Consequences for classification
**Table 5. Classifications that change class under alternative rules.**

| Rule | Expert panels (6,177) | Laboratories (87,873) |
|---|---|---|
| ClinGen cap (PP3 + PM1 ≤ 4) | 24 (12 P→LP, 12 LP→VUS) | 40 (30 P→LP, 10 LP→VUS) |
| PP3 + PM1 + PP2 ≤ 4 | 35 (17, 18) | 64 (40, 22) |
| **PP2 not added when PM1 applies** | 39 (28, 11) | **357 (125, 227)** |
| **PP3 limited to Supporting when PM1 applies** | 43 (22, 21) | 91 (40, 51) |

- **Later panel classifications.** Among the 2,507 approved in 2025–2026, limiting PP3 to Supporting changed 15 and the cap 16.
- **Developmental disorders.** In DD, uncapped co-use with de novo evidence placed 2,438 variants at Pathogenic, in a class the de novo excess implies is only 84% disease-associated.

""" + s[j:]
# ---------------------------------------------------------------- discussion
i = s.index("**What we found.**")
j = s.index("**Relation to previous work.**")
s = s[:i] + """**What we found.** Location-type evidence is double counted in missense classification.
- **The dependence.** On a label-free de novo scale, computational evidence (PP3) is not independent of hotspot or domain evidence (PM1) or gene-level evidence (PP2), and PM1 is not independent of PP2.
- **Why.** Predictors have learned where pathogenic variation falls, which is what PM1 and PP2 describe.
- **Where it holds.** On expert-panel and laboratory truth the dependence holds. On functional assays it is weak, because within assayed regions PM1 carries little information about molecular damage.

**PP2 beside PM1.** This is the clearest practical finding.
- **The guidance.** UK guidance advises against counting constraint evidence in both criteria.
- **Practice.** Laboratories stack them in 1,511 submissions and seven ClinGen panels in 174 classifications.
- **The evidence.** On laboratory truth, PP2 beside PM1 earned no detectable evidence.
- **Recommendation.** Do not add PP2 when PM1 applies, unless PP2 rests on evidence independent of PM1. This would change 357 laboratory and 39 expert classifications.

**PP3 above Supporting beside PM1.**
- **The evidence.** On the de novo scale, Moderate and Strong PP3 beside PM1 earned well under their nominal points, with intervals excluding them, for both predictors and in two cohorts. Panel-applied Moderate PP3 points the same way.
- **The cap misses it.** ClinGen's cap acts only when Strong PP3 meets Moderate PM1, which accounts for 0.5% of panel and 1.8% of laboratory PP3 + PM1 co-application.
- **Recommendation.** Counting PP3 at most as Supporting when PM1 applies is consistent with the data and needs no estimated parameter. It would change 43 expert and 91 laboratory classifications.
- **Its limits.** In panel-applied data, Moderate could not be distinguished from Supporting, and AlphaMissense-defined PP3 on panel truth did not show the deficit. We therefore present it as a rule for panels to check against their own classifications, not as an established standard.

**Supporting PP3 beside PM1** was calibrated on clinical truth. It earned less on the de novo scale, which measures association with severe developmental disorders rather than pathogenicity in general. Practice at Supporting, which is 94–96% of PP3 + PM1 co-application, therefore looks sound.

""" + s[j:]
rep("5. **Assays** tested dependence in ten genes outside NDD; key clinical genes such as BRCA1 and TP53 could not be tested.",
    "5. **Assays** tested dependence in ten genes outside NDD; key clinical genes such as BRCA1 and TP53 could not be tested.\n6. **Thin cells.** Several estimates rest on few benign variants: 17 for the overall panel PP3 credit and 6 for panel Moderate PP3. The de novo estimates above Supporting carry the claim.")
rep("6. **Coverage.** Transcripts match", "7. **Coverage.** Transcripts match")
i = s.index("**Conclusion.** PP3, PM1 and PP2 are not independent")
j = s.index("---", i)
s = s[:i] + """**Conclusion.** PP3, PM1 and PP2 are not independent, and the dependence matters in practice in two places:
- **PP2 beside PM1** adds no detectable evidence, yet laboratories and some expert panels add it routinely.
- **PP3 above Supporting beside PM1** adds little, below the reach of the current cap.

Two simple rules follow, both consistent with these data: do not add PP2 when PM1 applies, and count PP3 at most as Supporting when PM1 applies.

""" + s[j:]
rep("  - Table 3: co-application in practice;\n  - Table 3b: PP3 credit above Supporting from three independent analyses;\n  - Table 4: classification changes.",
    "  - Table 3: co-application in practice;\n  - Table 4: PP3 credit beside PM1 by strength, de novo and clinical;\n  - Table 5: classification changes in panels and laboratories.")
s = s.replace("SUPPLEMENT_v10.md", "SUPPLEMENT_v11.md")
open("MANUSCRIPT_v11.md", "w", encoding="utf-8").write(s)
print("ok")
