"""Builds MANUSCRIPT_v8.md / SUPPLEMENT_v8.md from draft 7: evidence stacking in expert practice (review 11)."""
s = open("MANUSCRIPT_v7.md", encoding="utf-8").read()


def rep(a, b):
    global s
    assert s.count(a) == 1, a[:70]
    s = s.replace(a, b)


rep("# How much do missense predictors add beyond gene and position? A label-free audit of computational and location evidence for de novo variants in neurodevelopmental disorders",
    "# Computational, regional and gene-level evidence are stacked but not independent: a label-free audit of missense criteria and of expert classification practice")
rep("*Draft 7, 7 October 2026. Supplementary material: `SUPPLEMENT_v7.md`.*",
    "*Draft 8, 7 October 2026. Supplementary material: `SUPPLEMENT_v8.md`.*")
i = s.index("## Abstract")
j = s.index("---", i)
s = s[:i] + """## Abstract

**Purpose.** ACMG/AMP classification sums computational (PP3), hotspot/domain (PM1) and gene-level (PP2) evidence as if independent. ClinGen caps PP3 + PM1 at Strong. We measured, without variant-level labels, whether these criteria are independent, how much predictors add beyond location, and how expert panels stack them.

**Methods.** Likelihood ratios came from the excess of de novo missense variants over mutation-rate expectation in 31,058 developmental-disorder trios, 12,316 autism probands and 4,256 public trios, marginally and within genes and 50-residue windows. We compared them with ClinVar, assays, and 6,177 expert-panel missense classifications re-scored without computational or location criteria.

**Results.** Predictors' evidence was mostly location: a median 0.32 of log-evidence survived conditioning on gene and window (a perfect predictor: 0.30–0.90). No pair of criteria was independent. Together, AlphaMissense PP3, PM1 and PP2 carried 3.65 points of de novo evidence where the rules sum 6.02 (R3 0.13 [0.11–0.14]; autism 0.18). Expert panels applied PP3 beside PM1 in 787 classifications, 783 of them below Strong, where the cap does not apply. There PP3 earned 0.60 points [0.11–1.10] against 1.24 applied (R 0.21 [0.12–0.31]). Seven panels also stacked PM1 with PP2. Crediting PP3 at its measured value would change 43–187 expert classifications, against 24 under the cap.

**Conclusion.** Location-type criteria are stacked in expert practice at strengths the existing cap does not reach. PP3 should be conditioned on PM1 at every strength, and PM1 + PP2 should not be summed unadjusted.

""" + s[j:]
rep("""### Clinical consequences""", """### Stacking in expert practice
**How panels combine the criteria.** Among 6,177 missense classifications in the ClinGen Evidence Repository (46 panels):

| Criteria applied | Classifications |
|---|---|
| PP3 | 2,990 |
| PM1 | 1,092 |
| PP2 | 592 |
| PP3 + PM1 | 787 |
| PM1 + PP2 | 174 |
| All three | 138 |

- **The cap rarely applies.** Where PP3 accompanied PM1, it was applied at Supporting in 605 classifications, at Moderate in 178 and at Strong in 4. The ClinGen cap binds only at Strong, so it constrains almost none of this stacking.
- **PM1 + PP2 is concentrated in seven panels:** RASopathy 62, FBN1 61, Monogenic Diabetes 22, Brain Malformations 11, PTEN 11, Pulmonary Hypertension 5 and ACADVL 2.

**Dependence of the panels' own criteria on independent truth.**
- **The truth set.** Classifications were re-scored after removing PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5: 2,942 missense classifications in 130 genes reached P/LP (≥ 6 points) or B/LB (≤ −1). Unlike the Table 4 analysis, this uses the panels' own PP3 rather than AlphaMissense, so no site mapping is needed.
- **PP3 × PM1.** Panel-applied PP3 and PM1 were strongly dependent: R 0.21 [0.12, 0.31] (gene-clustered).
- **What PP3 earns beside PM1.** Once a panel applied PM1, its PP3 earned 0.60 points [0.11, 1.10] of independent evidence, against a mean 1.24 points applied.
- **PM1 × PP2 and PP3 × PP2** could not be estimated on this truth: there were 1 and 2 B/LB variants in the joint cells.

**The triple stack on the de novo scale.**
- **Supporting.** When AlphaMissense PP3, PM1-like and PP2-like all apply, their joint de novo LR in DD was 24.6 (3.65 points [3.50, 3.79]), against a product of 196.7 (6.02 points): R3 0.13 [0.11, 0.14].
- **Strong.** At Strong PP3 the joint evidence was 4.05 points against 6.75 summed (R3 0.09).
- **Autism.** The pattern replicated: R3 0.18 [0.13, 0.24] at Supporting and 0.14 at Strong.
- **PM1-strict.** With PM1-strict, R3 was 0.10 and 0.07.

**Effect on expert classifications.** Re-scoring the 6,177 missense classifications under alternative rules changed:

| Rule | Classifications changed | P → LP | LP → VUS |
|---|---|---|---|
| ClinGen cap (PP3 + PM1 ≤ 4) | 24 | 12 | 12 |
| Cap on all three location-type criteria (PP3 + PM1 + PP2 ≤ 4) | 35 | 17 | 18 |
| PP2 not added when PM1 applies | 39 | 28 | 11 |
| PP3 credited at its measured value beside PM1 (1.10 / 0.60 points) | 43 / 187 | 22 / 84 | 21 / 102 |

The PP3 credit was estimated on the same panels, so the last row illustrates the size of the effect rather than predicting it.

### Clinical consequences""")
rep("""| Reclassified ClinVar (87 P, 887 B) | PM1-like | Not estimable (3 B/LB in the joint cell) | Not estimable |""",
    """| Reclassified ClinVar (87 P, 887 B) | PM1-like | Not estimable (3 B/LB in the joint cell) | Not estimable |
| Expert panel, panel-applied PP3 × PM1 (2,942 missense, 130 genes) | Panel-applied PM1 | 0.21 [0.12, 0.31] (any PP3 strength) | — |""")
rep("""- **What it implies for rules.**""", """- **Expert practice.** Panels stack these criteria, mostly at strengths the existing cap does not reach. On the panels' own classifications, PP3 beside PM1 was worth about half the points applied. When PP3, PM1 and PP2 all apply, the summed points overstate the de novo evidence by about 2.4 points, enough to move a variant across a class boundary.
- **What it implies for rules.**""")
rep("""We therefore support the cap, or a conditional rule, over uncapped co-use, and do not claim the cap's size is optimal.""",
    """We therefore support the cap, or a conditional rule, over uncapped co-use, and do not claim the cap's size is optimal.
- **Recommendations.**
  - Condition PP3 on PM1 at every strength, not only at Strong. The measured credit beside PM1 is about half a Supporting point.
  - Do not sum PM1 with PP2 unadjusted. PM1 × PP2 was the most dependent pair on the de novo scale, and seven expert panels currently stack them.
  - Treat PP3 + PM1 + PP2 as one location-type evidence budget. Panels that wish to keep co-application should calibrate the combination jointly against independent truth, as done here.""")
rep("""**Conclusion.** Most of what missense predictors contribute to pooled evaluations restates the gene and position of a variant. Measured label-free, the evidence they add within windows ranks predictors differently from pooled benchmarks and tracks within-gene clinical and functional performance. Computational, regional and gene-level criteria are not independent: PM1 and PP2 most of all, and PP3 and PM1 even within disease genes. Uncapped summation is not supported.""",
    """**Conclusion.** Most of what missense predictors contribute restates where a variant lies. The criteria that credit location separately are not independent: PM1 and PP2 most of all, and PP3 and PM1 even within disease genes. Stacked together, the three overstate the de novo evidence by about 2.4 points. Expert panels stack them mainly at strengths the existing cap does not reach, and crediting PP3 at its measured value beside PM1 would change up to 187 expert classifications. Location-type evidence should be budgeted jointly rather than summed.""")
rep("""4. **Truth sets.** Expert-panel truth covers 108 genes and is ascertained through panels.""",
    """4. **Truth sets.** Expert-panel truth covers 108–130 genes, is ascertained through panels, and was used both to estimate PP3's credit beside PM1 and to show the effect of applying it. The classification-change counts are therefore illustrative.""")
rep("10. PM1 × PP2 dependence.", "10. PM1 × PP2 dependence and the triple stack.\n11. Stacking in expert practice.")
s = s.replace("SUPPLEMENT_v7.md", "SUPPLEMENT_v8.md")
open("MANUSCRIPT_v8.md", "w", encoding="utf-8").write(s)
sup = open("SUPPLEMENT_v7.md", encoding="utf-8").read().replace(
    "*Draft 7, 7 October 2026. Accompanies `MANUSCRIPT_v7.md`.", "*Draft 8, 7 October 2026. Accompanies `MANUSCRIPT_v8.md`.")
tri = "\n".join(l for l in open("results_fullstudy_backup/triple_summary_log.txt", encoding="utf-8").read().splitlines() if "R3" in l)
tri = tri.split("] ", 1)[-1] if tri.startswith("[") else tri
st_m = open("results_fullstudy_backup/stacking_audit_missense_summary.txt", encoding="utf-8").read().strip()
st_a = open("results_fullstudy_backup/stacking_audit_all_summary.txt", encoding="utf-8").read().strip()
sup = sup.replace("## Supplementary Note 10. PM1 × PP2 dependence", "## Supplementary Note 10. PM1 × PP2 dependence and the triple stack")
sup = sup.rstrip("\n") + f"""

**Triple stack** (`fullstudy.py triple`; AlphaMissense PP3 × PM1 × PP2-like; transcript-clustered bootstrap, 1,000 replicates):

```
{tri}
```

---

## Supplementary Note 11. Stacking in expert practice

**Code and data.** `make_stacking_audit.py` runs on the ClinGen Evidence Repository export of 2026-10-07 (13,293 classifications; retracted entries removed).

**Truth.** Points from the applied criteria, excluding PP3, BP4, PM1, PP2, BP1, BP3, PM4 and PM5. P/LP at ≥ 6, B/LB at ≤ −1. Gene-clustered bootstrap, 1,000 replicates.

**Rules.** Class thresholds follow the Tavtigian points scale: P ≥ 10, LP 6–9, VUS 0–5, LB −1 to −6, B ≤ −7.

**Missense classifications:**

```
{st_m}
```

**All variant types** (for reference):

```
{st_a}
```
"""
open("SUPPLEMENT_v8.md", "w", encoding="utf-8").write(sup + "\n")
print("ok")
