"""A19 H1 (pre-registered, prereg_fullstudy/A19_V14_PREREG.md): which rule for PP3 beside PM1 is best calibrated out of
sample on laboratory truth? Gene-grouped 5-fold CV x 20 repeats; prediction logit p = a + r(PP3, PM1) ln 2.406 with
a = logit truth rate among PM1-applied, PP3-absent training SCVs; held-out log-loss; submitter-clustered bootstrap of
per-SCV loss differences (averaged over repeats). Rules: R0 uncapped, R1 ClinGen cap (PP3+PM1 <= 4), R2 PP3 <= 1,
R3 discount k (training earned/applied), R4 learned credit per strength (training), R5 no PP3 credit.
Usage: python make_rules_cv.py [vintage]"""
import math
import sys

import numpy as np
import pandas as pd

V = next((x for x in sys.argv[1:] if not x.startswith("--")), "2026-10")
# Deviation (logged in prereg_fullstudy/A19_DEVIATIONS.md): '--acmg' uses the ACMG convention that an absent code carries
# no evidence: baseline a = logit truth rate among ALL PM1-applied SCVs, and learned credit = ln LR(PP3 | PM1) as in v13.
ACMG = "--acmg" in sys.argv
TAG = V + ("_acmg" if ACMG else "")
LN = math.log(2.406)
rng = np.random.default_rng(20261008)
s = pd.read_csv(f"work/cvt/lab_codes_{V}.tsv.gz", sep="\t", usecols=["SCV", "Submitter", "gene", "PP3", "PM1", "indep", "total"])
s = s[(s.PM1 > 0) & ((s.indep >= 6) | (s.indep <= -1))].copy()
s["y"] = (s.indep >= 6).astype(int)
s["PP3"] = s.PP3.clip(lower=0)
genes = s.gene.unique()


def logit(p):
    p = min(max(p, 1e-6), 1 - 1e-6)
    return math.log(p / (1 - p))


def credit(tr, strength_mask):
    """points earned by PP3 (cells PP3 in mask vs PP3 = 0) among PM1-applied training SCVs, 0.5 smoothing"""
    a = tr[strength_mask(tr.PP3)]
    b = tr if ACMG else tr[tr.PP3 == 0]
    o1 = (a.y.sum() + 0.5) / ((1 - a.y).sum() + 0.5)
    o0 = (b.y.sum() + 0.5) / ((1 - b.y).sum() + 0.5)
    return math.log(o1 / o0) / LN


def rules(tr):
    c = {k: credit(tr, lambda x, k=k: x == k) for k in [1, 2, 4, 8]}
    applied = tr.loc[tr.PP3 > 0, "PP3"].mean()
    k = credit(tr, lambda x: x > 0) / applied
    return {"R0 uncapped": lambda d: d.PP3,
            "R1 ClinGen cap": lambda d: np.minimum(d.PP3, np.maximum(0, 4 - d.PM1)),
            "R2 PP3 <= Supporting": lambda d: np.minimum(d.PP3, 1),
            "R3 discount": lambda d: d.PP3 * k,
            "R4 learned credit": lambda d: d.PP3.map(lambda x: 0.0 if x == 0 else c.get(int(x), c[1])),
            "R5 no PP3 credit": lambda d: d.PP3 * 0.0}


NAMES = list(rules(s).keys())
loss = {n: np.zeros(len(s)) for n in NAMES}
REP, K = 20, 5
idx = {g: i for i, g in enumerate(genes)}
for rep in range(REP):
    fold = rng.permutation(len(genes)) % K
    gf = s.gene.map(lambda g: fold[idx[g]]).to_numpy()
    for f in range(K):
        tr, te = s[gf != f], s[gf == f]
        a = logit(tr.y.mean() if ACMG else tr.loc[tr.PP3 == 0, "y"].mean())
        for n, r in rules(tr).items():
            z = a + np.asarray(r(te), float) * LN
            p = 1 / (1 + np.exp(-z))
            p = np.clip(p, 1e-6, 1 - 1e-6)
            loss[n][gf == f] += -(te.y * np.log(p) + (1 - te.y) * np.log(1 - p)) / REP
L = pd.DataFrame(loss)
L["Submitter"] = s.Submitter.to_numpy()
subs = L.Submitter.unique()
by = L.groupby("Submitter")
sums, cnt = by[NAMES].sum(), by.size()
out = [f"A19 H1, vintage {V}{' (ACMG-convention baseline: SENSITIVITY, not pre-registered)' if ACMG else ' (pre-registered specification)'}: PM1-applied SCVs with laboratory truth {len(s)} (P/LP {s.y.sum()}, B/LB {(1 - s.y).sum()}; "
       f"PP3 applied {int((s.PP3 > 0).sum())}; genes {len(genes)}; submitters {len(subs)})",
       "Held-out mean log-loss (gene-grouped 5-fold CV x 20):"]
for n in NAMES:
    out.append(f"  {n:22s} {L[n].mean():.4f}")
out.append("Difference vs R1 ClinGen cap (submitter-clustered bootstrap, 2,000):")
res = []
for n in NAMES:
    if n == "R1 ClinGen cap":
        continue
    diff = (sums[n] - sums["R1 ClinGen cap"])
    bs = []
    for _ in range(2000):
        pick = rng.choice(len(subs), len(subs))
        bs.append(diff.iloc[pick].sum() / cnt.iloc[pick].sum())
    lo, hi = np.quantile(bs, [0.025, 0.975])
    d = diff.sum() / cnt.sum()
    res.append(dict(rule=n, delta_logloss_vs_cap=d, lo=lo, hi=hi))
    out.append(f"  {n:22s} {d:+.4f} [{lo:+.4f}, {hi:+.4f}]")
r4 = [x for x in res if x["rule"] == "R4 learned credit"][0]
out.append(f"H1a decision (pre-specified): R4 better than the cap if upper bound < 0 -> "
           f"{'CONFIRMED' if r4['hi'] < 0 else 'NOT confirmed (cap retained)'}")
full = rules(s)
c_all = {k: credit(s, lambda x, k=k: x == k) for k in [1, 2, 4]}
out.append("Full-data learned credit (points) for PP3 beside PM1: " + ", ".join(f"nominal {k}: {v:.2f}" for k, v in c_all.items()))
k_all = credit(s, lambda x: x > 0) / s.loc[s.PP3 > 0, "PP3"].mean()
out.append(f"Full-data discount factor k = {k_all:.2f}")
# classifications changed per rule among all PM1 & PP3 SCVs (all SCVs, not only truth)
allx = pd.read_csv(f"work/cvt/lab_codes_{V}.tsv.gz", sep="\t", usecols=["PP3", "PM1", "total"])
allx = allx[(allx.PM1 > 0) & (allx.PP3 > 0)]


def cls(t):
    return np.where(t >= 10, "P", np.where(t >= 6, "LP", np.where(t >= 0, "VUS", "B")))


base = cls(allx.total)
for n, r in full.items():
    t = allx.total - allx.PP3 + np.asarray(r(allx), float)
    ch = cls(t) != base
    out.append(f"  classifications changed vs reported points under {n:22s}: {int(ch.sum())} of {len(allx)} "
               f"(P->LP {int(((base == 'P') & (cls(t) == 'LP')).sum())}, LP->VUS {int(((base == 'LP') & (cls(t) == 'VUS')).sum())})")
pd.DataFrame(res).to_csv(f"results_fullstudy/a19_h1_rules_{TAG}.tsv", sep="\t", index=False)
open(f"results_fullstudy/a19_h1_rules_{TAG}_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
