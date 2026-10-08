"""Draft 14, roadmap section 3 (exploratory; not pre-registered): validity gaps in the laboratory analysis.
1. PM1 strength and PP3+PM1 sums; how often the ClinGen cap binds.
2. Reported vs re-scored class: where the cap binds and changes the points-based class, which class did the laboratory
   report? Reported = capped class suggests the cap (or a platform enforcing it) is already applied.
3. Templated submitters: template ratio = 1 - distinct 60-character comment openings / SCVs (PP3 & PM1 SCVs); submitters
   with ratio >= 0.3 and >= 20 such SCVs are flagged; key estimates re-run without them.
4. Parser validation sheet: 400 comments stratified by submitter (<= 5 each), oversampling PP3/PM1/PP2, for two
   independent human annotators (kappa computed by make_parser_kappa.py once both columns are filled)."""
import math

import numpy as np
import pandas as pd

LN = math.log(2.406)
rng = np.random.default_rng(20261008)
s = pd.read_csv("work/cvt/lab_codes_2026-10.tsv.gz", sep="\t")
out = []
pm = s.PM1 > 0
both = pm & (s.PP3 > 0)
out.append("1. PM1 strength when applied (points: n): " + ", ".join(f"{k}: {v}" for k, v in s.loc[pm, "PM1"].value_counts().sort_index().items()))
out.append("   PP3 + PM1 sum when both applied: " + ", ".join(f"{k}: {v}" for k, v in (s.PP3 + s.PM1)[both].value_counts().sort_index().items())
           + f"; cap binds (sum > 4) in {int(((s.PP3 + s.PM1) > 4)[both].sum())} of {int(both.sum())}")


def cls_pts(t):
    return np.where(t >= 10, "P", np.where(t >= 6, "LP", np.where(t >= 0, "VUS", np.where(t >= -6, "LB", "B"))))


def cls_rep(x):
    x = str(x).lower()
    if "conflicting" in x:
        return "other"
    if x.startswith("pathogenic/likely") or x == "pathogenic":
        return "P" if x == "pathogenic" else "P/LP"
    if x.startswith("likely pathogenic"):
        return "LP"
    if "uncertain" in x:
        return "VUS"
    if "likely benign" in x:
        return "LB"
    if x.startswith("benign"):
        return "B"
    return "other"


s["cls_uncapped"] = cls_pts(s.total)
capped = s.total - s.PP3 + np.minimum(s.PP3, np.maximum(0, 4 - s.PM1)) * (s.PM1 > 0) + s.PP3 * (s.PM1 <= 0)
s["cls_capped"] = cls_pts(capped)
s["cls_reported"] = s.ClinicalSignificance.map(cls_rep)
agree = (s.cls_reported == s.cls_uncapped).mean()
out.append(f"2. Reported class equals points-based class (uncapped) in {agree:.1%} of SCVs listing codes")
d = s[both & (s.cls_uncapped != s.cls_capped)]
tab = d.groupby(["cls_uncapped", "cls_capped", "cls_reported"]).size().reset_index(name="n")
out.append(f"   SCVs where the cap changes the points-based class: {len(d)}; reported class:")
for r in tab.itertuples():
    out.append(f"     uncapped {r.cls_uncapped:3s} -> capped {r.cls_capped:3s}: reported {r.cls_reported:5s} n {r.n}")
rep_capped = (d.cls_reported == d.cls_capped).sum()
out.append(f"   reported = capped class in {rep_capped} of {len(d)}; reported = uncapped class in {(d.cls_reported == d.cls_uncapped).sum()}")
for sub, g in d.groupby("Submitter"):
    if (g.cls_reported == g.cls_capped).sum() >= 3:
        out.append(f"     submitter reporting the capped class >= 3 times: {sub} ({(g.cls_reported == g.cls_capped).sum()} of {len(g)})")

# 3. templated submitters
b = s[both]
g = b.groupby("Submitter").agg(n=("PP3", "size"), uniq=("Description", lambda x: x.str[:60].nunique()))
g["ratio"] = 1 - g.uniq / g.n
flag = g[(g.ratio >= 0.3) & (g.n >= 20)]
out.append(f"3. Templated submitters flagged (ratio >= 0.3, >= 20 PP3&PM1 SCVs): {len(flag)} covering {int(flag.n.sum())} of {len(b)} PP3&PM1 SCVs: "
           + "; ".join(f"{i} ({r.n}, {r.ratio:.2f})" for i, r in flag.iterrows()))
T = s[(s.indep >= 6) | (s.indep <= -1)].copy()
T["y"] = (T.indep >= 6).astype(int)
for c in ["PP3", "PM1", "PP2"]:
    T[c + "_on"] = (T[c] > 0).astype(int)


def lr(d, m):
    p, q = d[d.y == 1], d[d.y == 0]
    return ((m(p).sum() + 0.5) / (len(p) + 1)) / ((m(q).sum() + 0.5) / (len(q) + 1))


def R(d):
    return lr(d, lambda x: (x.PP3_on == 1) & (x.PM1_on == 1)) / (lr(d, lambda x: x.PP3_on == 1) * lr(d, lambda x: x.PM1_on == 1))


def credit(d, code):
    e = d[d.PM1_on == 1]
    a = e[e[code + "_on"] == 1]
    return math.log(((a.y.sum() + .5) / ((1 - a.y).sum() + .5)) / ((e.y.sum() + .5) / ((1 - e.y).sum() + .5))) / LN


def ci(d, f, n=1000):
    gi = d.groupby("gene").indices
    k = list(gi)
    return np.nanquantile([f(d.iloc[np.concatenate([gi[k[i]] for i in rng.choice(len(k), len(k))])]) for _ in range(n)], [.025, .975])


for lab, d in [("all submitters", T), ("templated submitters removed", T[~T.Submitter.isin(flag.index)])]:
    r, (rl, rh) = R(d), ci(d, R)
    c3, (al, ah) = credit(d, "PP3"), ci(d, lambda x: credit(x, "PP3"))
    c2, (bl, bh) = credit(d, "PP2"), ci(d, lambda x: credit(x, "PP2"))
    out.append(f"   {lab:29s}: n {len(d)}; R(PP3,PM1) {r:.2f} [{rl:.2f}, {rh:.2f}]; PP3 credit beside PM1 {c3:+.2f} [{al:+.2f}, {ah:+.2f}]; "
               f"PP2 credit beside PM1 {c2:+.2f} [{bl:+.2f}, {bh:+.2f}]")

# 4. parser annotation sheet
cand = s.assign(w=1 + 2 * ((s.PP3 > 0) | (s.PM1 > 0) | (s.PP2 > 0)))
rows = []
for sub, g2 in cand.groupby("Submitter"):
    k = min(5, len(g2))
    g2 = g2.assign(r=rng.random(len(g2))).sort_values(["w", "r"], ascending=[False, True])
    rows.append(g2.head(k))
sheet = pd.concat(rows)
sheet = sheet.assign(r=rng.random(len(sheet))).sort_values(["w", "r"], ascending=[False, True]).head(400)
sheet = sheet[["SCV", "Submitter", "gene", "Description", "codes"]].rename(columns={"codes": "parser_codes"})
for a in ["annotator1", "annotator2"]:
    sheet[a + "_codes"] = ""
sheet.to_csv("results_fullstudy/parser_annotation_sheet_v14.tsv", sep="\t", index=False)
out.append(f"4. Parser annotation sheet: {len(sheet)} comments from {sheet.Submitter.nunique()} submitters -> results_fullstudy/parser_annotation_sheet_v14.tsv "
           f"(annotators fill annotator1_codes / annotator2_codes as JSON code->points, blind to parser_codes)")
open("results_fullstudy/validity_v14_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
