"""Review 20, laboratory analyses (descriptive; not pre-specified).
(1) Within-submitter credit once PM1 applies: among PM1-applied submissions, the likelihood ratio of C (PP2; Supporting PP3)
    for P/LP vs B/LB truth is pooled across submitters with the Mantel-Haenszel risk-ratio estimator (truth = exposure,
    C applied = outcome, submitter = stratum), so only within-submitter contrasts contribute. Credit = ln(LR_MH)/ln 2.406.
    Submitter-clustered bootstrap.
(2) Held-out test of the rule "PP2 not added when PM1 applies": ClinVar 2023-01 laboratory Likely Pathogenic missense
    submissions whose class drops to VUS when PP2 is removed beside PM1 (rule-affected) are compared with Likely Pathogenic
    submissions on the same total points without PP2 beside PM1 (controls). Outcomes, read from ClinVar 2026-10:
    O1, the same SCV is no longer P/LP (submitter's own reclassification); O2, the variant's other submitters (2026)
    classify it VUS/LB/B by majority. Crude and submitter-stratified Mantel-Haenszel odds ratios.
Run in the WSL repo: python make_review20_lab.py"""
import gzip, math, re
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

CVT = "work/cvt/"
BASE = {"PVS": 8, "PS": 4, "PM": 2, "PP": 1, "BA": -8, "BS": -4, "BP": -1}
MOD = {"verystrong": 8, "strong": 4, "moderate": 2, "supporting": 1, "standalone": 8}
EXCL = {"PP3", "BP4", "PM1", "PP2", "BP1", "BP3", "PM4", "PM5"}
CODE = re.compile(r"\b(PVS|PS|PM|PP|BA|BS|BP)([1-9])(?:[_ ]?\(?(very[ _]?strong|strong|moderate|supporting|stand[ _-]?alone)\)?)?\b", re.I)
NEG_A = re.compile(r"^(?:/[A-Z]{2,3}\d)?\W*(?:criteria\s+|criterion\s+)?(?:was\s+|is\s+|were\s+|are\s+)?(?:not\b|n/a|unmet|absent\b)", re.I)
NEG_B = re.compile(r"(?:\b(?:not|no|without|neither|nor)\b|does\s+not\s+meet|do\s+not\s+meet|not\s+meet|insufficient\s+for|below\s+the\s+threshold\s+for)[^.;]{0,40}$", re.I)
MIS = re.compile(r"\(p\.(?:[A-Z][a-z]{2})\d+(?!Ter|=|del|ins|dup|fs|ext)[A-Z][a-z]{2}\)")
OUT = "results_fullstudy/review20_lab_summary.txt"


def parse(txt):
    ev = {}
    for m in CODE.finditer(txt or ""):
        if NEG_A.match(txt[m.end():m.end() + 25]) or NEG_B.search(txt[max(0, m.start() - 60):m.start()]):
            continue
        fam, num, mod = m.group(1).upper(), m.group(2), m.group(3)
        v = BASE[fam]
        if mod:
            v = MOD[re.sub(r"[ _-]", "", mod.lower())] * (1 if v > 0 else -1)
        ev[fam + num] = v
    return ev


def load(path, cols):
    with gzip.open(path, "rt") as fh:
        n = next(i for i, line in enumerate(fh) if line.startswith("#VariationID\t"))
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i == n:
                avail = line.rstrip("\n").split("\t")
                break
    use = [c for c in cols if c in avail]
    d = pd.read_csv(path, sep="\t", skiprows=n, header=0, dtype=str, quoting=3, compression="gzip", on_bad_lines="skip",
                    usecols=use).rename(columns={"#VariationID": "VariationID"})
    return d


COLS = ["#VariationID", "ClinicalSignificance", "Description", "OriginCounts", "ReviewStatus", "Submitter", "SCV"]
vs = pd.read_csv(CVT + "variant_summary_2026-10.txt.gz", sep="\t", dtype=str, compression="gzip",
                 usecols=["VariationID", "Name", "GeneSymbol", "Type", "Assembly"])
vs = vs[(vs.Assembly == "GRCh38") & (vs.Type == "single nucleotide variant")].drop_duplicates("VariationID")
vs = vs[vs.Name.str.contains(MIS, na=False)]


def lab_set(path):
    s = load(path, COLS)
    s = s[s.OriginCounts.fillna("").str.contains("germline", case=False)]
    s = s[~s.ReviewStatus.fillna("").str.contains("expert panel|practice guideline", case=False)]
    s = s.merge(vs[["VariationID", "GeneSymbol"]], on="VariationID")
    return s


import os
out = []
# ---------------- (1) within-submitter credit, 2026-10
if os.environ.get("SKIP1"):
    out = open(OUT).read().splitlines()[:4]
s = lab_set(CVT + "submission_summary_2026-10.txt.gz")
s = s[s.Description.fillna("").str.contains(r"\b(?:PVS|PS|PM|PP|BA|BS|BP)[1-9]\b", regex=True)].copy()
ev = s.Description.map(parse)
for c in ("PP3", "PM1", "PP2"):
    s[c] = ev.map(lambda e, c=c: e.get(c, 0))
s["indep"] = ev.map(lambda e: sum(v for c, v in e.items() if c not in EXCL))
s["gene"] = s.GeneSymbol
T = s[(s.indep >= 6) | (s.indep <= -1)].copy()
T["truth"] = (T.indep >= 6).astype(int)
P = T[T.PM1 > 0].copy()
out += [] if os.environ.get("SKIP1") else [f"(1) Within-submitter credit once PM1 applies (laboratory truth; {len(T)} submissions with truth, {len(P)} with PM1)"]


def lr_mh(z, cond):
    num = den = 0.0
    k = 0
    for _, g in z.groupby("Submitter"):
        y = g.truth.to_numpy(); a_ = cond(g).to_numpy()
        n1, n0 = (y == 1).sum(), (y == 0).sum()
        if n1 == 0 or n0 == 0:
            continue
        N = n1 + n0
        a = (a_ & (y == 1)).sum(); c = (a_ & (y == 0)).sum()
        num += a * n0 / N; den += c * n1 / N
        k += 1
    return (num / den if den > 0 else np.nan), k


def lr_crude(z, cond):
    y = z.truth.to_numpy(); a_ = cond(z).to_numpy()
    return a_[y == 1].mean() / max(a_[y == 0].mean(), 1e-9)


rng = np.random.default_rng(20261020)
su = P.Submitter.unique(); si = P.groupby("Submitter").indices
for lab, cond, applied in ([] if os.environ.get("SKIP1") else  (("PP2", lambda g: g.PP2 > 0, 1.0), ("PP3 Supporting", lambda g: g.PP3 == 1, 1.0),
                           ("PP3 any strength", lambda g: g.PP3 > 0, None))):
    lc = lr_crude(P, cond); lm, k = lr_mh(P, cond)
    bs = []
    for _ in range(1000):
        zz = P.iloc[np.concatenate([si[c] for c in rng.choice(su, len(su))])]
        v, _ = lr_mh(zz, cond)
        bs.append(math.log(v) / math.log(2.406) if v and v > 0 and np.isfinite(v) else np.nan)
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    lo90, hi90 = np.nanpercentile(bs, [5, 95])
    sel = cond(P)
    inf = P.groupby("Submitter").filter(lambda g: g.truth.nunique() == 2)
    out.append(f"  {lab:16s}: pooled (between + within) {math.log(lc) / math.log(2.406):.2f}; within-submitter (MH) "
               f"{math.log(lm) / math.log(2.406):.2f} [{lo:.2f}, {hi:.2f}] (90% {lo90:.2f}, {hi90:.2f}); informative submitters {k}, "
               f"their PM1 submissions {len(inf)} (C applied {int(cond(inf).sum())}; B/LB with C {int((cond(inf) & (inf.truth == 0)).sum())}); "
               f"applied {P.loc[sel, 'PP3' if 'PP3' in lab else 'PP2'].mean():.2f}")
open(OUT, "w").write("\n".join(out) + "\n")
print("\n".join(out))

# ---------------- (2) held-out reclassification test, 2023-01 -> 2026-10
o = lab_set(CVT + "submission_summary_2023-01.txt.gz")
o = o[o.Description.fillna("").str.contains(r"\b(?:PVS|PS|PM|PP|BA|BS|BP)[1-9]\b", regex=True)].copy()
eo = o.Description.map(parse)
for c in ("PM1", "PP2"):
    o[c] = eo.map(lambda e, c=c: e.get(c, 0))
o["pts"] = eo.map(lambda e: sum(e.values())).astype(float)
o["LPcall"] = o.ClinicalSignificance.fillna("").str.strip().str.lower().eq("likely pathogenic")
stack = (o.PM1 > 0) & (o.PP2 > 0)
o["affected"] = o.LPcall & stack & (o.pts >= 6) & (o.pts < 10) & (o.pts - o.PP2 < 6)
o["control"] = o.LPcall & ~stack & (o.pts >= 6) & (o.pts < 7)
out.append(f"(2) Held-out test, ClinVar 2023-01 laboratory missense submissions with ACMG codes: {len(o)}; Likely Pathogenic "
           f"{int(o.LPcall.sum())}; rule-affected (LP only through PP2 beside PM1) {int(o.affected.sum())}; controls (LP at the "
           f"same points, 6, without PP2 beside PM1) {int(o.control.sum())}")

n_all = lab_set(CVT + "submission_summary_2026-10.txt.gz")
n_all["acc"] = n_all.SCV.str.split(".").str[0]
o["acc"] = o.SCV.str.split(".").str[0]
plp = lambda x: x.fillna("").str.contains("pathogenic", case=False) & ~x.fillna("").str.contains("conflict|uncertain|benign", case=False)
nonplp = lambda x: x.fillna("").str.contains("uncertain|benign", case=False)
same = n_all.drop_duplicates("acc").set_index("acc").ClinicalSignificance
o["same_2026"] = o.acc.map(same)
o["O1"] = np.where(o.same_2026.isna(), np.nan, (~plp(o.same_2026.astype(str))).astype(float))
# other submitters' 2026 calls
g = n_all[["VariationID", "Submitter", "ClinicalSignificance"]].copy()
g["plp"] = plp(g.ClinicalSignificance).astype(int); g["non"] = nonplp(g.ClinicalSignificance).astype(int)
key = o[["VariationID", "Submitter"]].drop_duplicates()
m = key.merge(g, on="VariationID", suffixes=("", "_o"))
m = m[m.Submitter != m.Submitter_o]
agg = m.groupby(["VariationID", "Submitter"]).agg(n_o=("plp", "size"), plp_o=("plp", "sum"), non_o=("non", "sum")).reset_index()
o = o.merge(agg, on=["VariationID", "Submitter"], how="left")
o["O2"] = np.where(o.n_o.fillna(0) == 0, np.nan, (o.non_o > o.plp_o).astype(float))


def mh_or(z, outc):
    num = den = 0.0
    for _, g_ in z.groupby("Submitter"):
        e = g_.affected.to_numpy(); y = g_[outc].to_numpy()
        a = ((e) & (y == 1)).sum(); b = ((e) & (y == 0)).sum(); c = ((~e) & (y == 1)).sum(); d = ((~e) & (y == 0)).sum()
        N = a + b + c + d
        if N == 0 or (a + b) == 0 or (c + d) == 0:
            continue
        num += a * d / N; den += b * c / N
    return num / den if den > 0 else np.nan


Z = o[o.affected | o.control].copy()
for outc, lab in (("O1", "same SCV no longer P/LP in 2026 (submitter's own reclassification)"),
                  ("O2", "other submitters in 2026 classify VUS/LB/B by majority")):
    z = Z.dropna(subset=[outc])
    ra, rc = z.loc[z.affected, outc].mean(), z.loc[z.control, outc].mean()
    na, nc = int(z.affected.sum()), int(z.control.sum())
    xa, xc = int(z.loc[z.affected, outc].sum()), int(z.loc[z.control, outc].sum())
    orc = (xa / max(na - xa, 1e-9)) / max(xc / max(nc - xc, 1e-9), 1e-9)
    orm = mh_or(z, outc)
    tab = [[xa, na - xa], [xc, nc - xc]]
    ex = odds_ratio(tab); ci = ex.confidence_interval(0.95)
    out.append(f"  {outc} ({lab}): affected {xa}/{na} ({ra:.1%}) vs controls {xc}/{nc} ({rc:.1%}); conditional ML OR {ex.statistic:.2f} "
               f"[{ci.low:.2f}, {ci.high:.2f}] (exact); Fisher two-sided p {fisher_exact(tab)[1]:.4f}; submitter-stratified MH OR {orm:.2f}; "
               f"submitters with affected {z.loc[z.affected, 'Submitter'].nunique()}")
open(OUT, "w").write("\n".join(out) + "\n")
print("\n".join(out[-3:]))
