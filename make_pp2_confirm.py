"""A19 H2 (pre-registered): PP2 credit beside PM1 on laboratory truth.
H2  : Bayesian binomial mixed GLM (variational Bayes) among PM1-applied SCVs with truth: y ~ PP2 + PP3, random
      intercepts for submitter and gene; PP2 credit = coefficient / ln 2.406; confirmed if 95% CrI upper bound < 1.
H2b : PM1 rationale from the comment text (window of 250 characters around the PM1 mention; pre-registered keyword
      rules): constraint-based / hotspot-domain / other. PP2 credit by class (ACMG convention, gene-clustered bootstrap)
      and PP2 x constraint interaction in the mixed model.
H2c : within-submitter PP2 - PP3 credit difference beside PM1 (submitter-clustered bootstrap); descriptive.
Baselines (descriptive): PP2 and PP3 credit with and without PM1 (ACMG convention: ln LR(code | PM1 status) / ln 2.406).
Usage: python make_pp2_confirm.py [vintage]"""
import math
import re
import sys

import numpy as np
import pandas as pd
from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM

V = sys.argv[1] if len(sys.argv) > 1 else "2026-10"
LN = math.log(2.406)
rng = np.random.default_rng(20261008)
s = pd.read_csv(f"work/cvt/lab_codes_{V}.tsv.gz", sep="\t",
                usecols=["SCV", "Submitter", "gene", "PP3", "PM1", "PP2", "indep", "Description"])
s = s[(s.indep >= 6) | (s.indep <= -1)].copy()
s["y"] = (s.indep >= 6).astype(int)
for c in ["PP3", "PM1", "PP2"]:
    s[c + "_on"] = (s[c] > 0).astype(int)

CONSTRAINT = re.compile(r"constraint|missense[- ]depleted|depleted of|regional|o/e|z-?score|tolerat|intolerant|lack of benign|"
                        r"absence of benign|few benign|no benign", re.I)
HOTSPOT = re.compile(r"hot ?spot|domain|critical|residue|region|motif|binding|active site", re.I)


def pm1_class(txt):
    t = txt or ""
    m = re.search(r"\bPM1\b", t)
    if not m:
        return "other"
    w = t[max(0, m.start() - 250):m.end() + 250]
    if CONSTRAINT.search(w):
        return "constraint"
    if HOTSPOT.search(w):
        return "hotspot/domain"
    return "other"


s["pm1_class"] = np.where(s.PM1_on == 1, s.Description.map(pm1_class), "")


def credit(df, code, cond):
    """ACMG-convention credit of `code` among rows `cond`: ln[odds(y | code, cond) / odds(y | cond)] / ln 2.406"""
    d = df[cond(df)]
    a = d[d[code + "_on"] == 1]
    o1 = (a.y.sum() + 0.5) / ((1 - a.y).sum() + 0.5)
    o0 = (d.y.sum() + 0.5) / ((1 - d.y).sum() + 0.5)
    return math.log(o1 / o0) / LN, len(a), int((1 - a.y).sum())


def boot(df, f, unit, n=1000):
    units = df[unit].unique()
    g = {u: i for u, i in df.groupby(unit).groups.items()}
    vals = []
    for _ in range(n):
        pick = rng.choice(units, len(units))
        b = df.loc[np.concatenate([g[u] for u in pick])]
        try:
            vals.append(f(b))
        except Exception:
            pass
    return np.nanquantile(vals, [0.025, 0.975])


out = [f"A19 H2, vintage {V}: laboratory truth SCVs {len(s)} (P/LP {s.y.sum()}, B/LB {(1 - s.y).sum()}); PM1-applied {int(s.PM1_on.sum())}"]
pm1 = s[s.PM1_on == 1].copy()
# ---- H2 mixed model
X = pm1[["PP2_on", "PP3_on"]].astype(float).assign(const=1.0)[["const", "PP2_on", "PP3_on"]]
Z = pd.get_dummies(pm1[["Submitter", "gene"]].astype(str), prefix=["sub", "gene"]).astype(float)
ident = np.array([0 if c.startswith("sub_") else 1 for c in Z.columns])
m = BinomialBayesMixedGLM(pm1.y.to_numpy(), X.to_numpy(), Z.to_numpy(), ident, fep_names=list(X.columns),
                          vcp_names=["submitter", "gene"]).fit_vb()
mu, sd = m.fe_mean, m.fe_sd
for i, nme in enumerate(X.columns[1:], start=1):
    lo, hi = (mu[i] - 1.96 * sd[i]) / LN, (mu[i] + 1.96 * sd[i]) / LN
    out.append(f"H2 mixed model (PM1-applied, n {len(pm1)}, B/LB {(1 - pm1.y).sum()}): {nme} credit {mu[i] / LN:+.2f} points "
               f"[95% CrI {lo:+.2f}, {hi:+.2f}] (pre-registered: present vs absent)")
    if nme == "PP2_on":
        out.append(f"H2 decision: PP2 credit below 1 point if CrI upper bound < 1 -> {'CONFIRMED' if hi < 1 else 'NOT confirmed'}")
out.append(f"  random-effect SDs (log-odds): submitter {np.exp(m.vcp_mean[0]):.2f}, gene {np.exp(m.vcp_mean[1]):.2f}")
# ---- ACMG-convention credits, baselines
out.append("ACMG-convention credit (points [gene-clustered 95% CI]; n with code; B/LB with code):")
for code, lab, cond in [("PP2", "PP2 beside PM1", lambda d: d.PM1_on == 1), ("PP2", "PP2 without PM1", lambda d: d.PM1_on == 0),
                        ("PP3", "PP3 beside PM1", lambda d: d.PM1_on == 1), ("PP3", "PP3 without PM1", lambda d: d.PM1_on == 0)]:
    c, n, nb = credit(s, code, cond)
    lo, hi = boot(s, lambda b, code=code, cond=cond: credit(b, code, cond)[0], "gene")
    out.append(f"  {lab:18s}: {c:+.2f} [{lo:+.2f}, {hi:+.2f}]; n {n}; B/LB {nb}")
# ---- H2b PM1 rationale
out.append("H2b PM1 rationale classes (PM1-applied SCVs with truth): " + ", ".join(f"{k} {v}" for k, v in pm1.pm1_class.value_counts().items()))
for k in ["constraint", "hotspot/domain", "other"]:
    sub = pm1[pm1.pm1_class == k]
    if len(sub) < 30:
        continue
    c, n, nb = credit(sub, "PP2", lambda d: d.PM1_on == 1)
    lo, hi = boot(sub, lambda b: credit(b, "PP2", lambda d: d.PM1_on == 1)[0], "gene")
    flag = "" if nb >= 20 else " (descriptive: < 20 B/LB)"
    out.append(f"  PP2 credit beside {k:15s} PM1: {c:+.2f} [{lo:+.2f}, {hi:+.2f}]; n PP2 {n}; B/LB {nb}{flag}")
pm1["constraint"] = (pm1.pm1_class == "constraint").astype(float)
X2 = pm1[["PP2_on", "PP3_on", "constraint"]].astype(float).assign(PP2xC=lambda d: d.PP2_on * d.constraint, const=1.0)
X2 = X2[["const", "PP2_on", "PP3_on", "constraint", "PP2xC"]]
m2 = BinomialBayesMixedGLM(pm1.y.to_numpy(), X2.to_numpy(), Z.to_numpy(), ident, fep_names=list(X2.columns),
                           vcp_names=["submitter", "gene"]).fit_vb()
i = list(X2.columns).index("PP2xC")
out.append(f"  mixed-model interaction PP2 x constraint-based PM1: {m2.fe_mean[i] / LN:+.2f} points "
           f"[{(m2.fe_mean[i] - 1.96 * m2.fe_sd[i]) / LN:+.2f}, {(m2.fe_mean[i] + 1.96 * m2.fe_sd[i]) / LN:+.2f}]")
# ---- H2c within-submitter contrast
f = lambda b: credit(b, "PP2", lambda d: d.PM1_on == 1)[0] - credit(b, "PP3", lambda d: d.PM1_on == 1)[0]
d0 = f(s)
lo, hi = boot(s[s.Submitter.isin(pm1.Submitter.unique())], f, "Submitter")
out.append(f"H2c PP2 - PP3 credit beside PM1: {d0:+.2f} [submitter-clustered {lo:+.2f}, {hi:+.2f}] (descriptive)")
open(f"results_fullstudy/a19_h2_pp2_{V}_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
