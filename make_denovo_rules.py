"""A20 (pre-registered, prereg_fullstudy/A20_DENOVO_RULES_PREREG.md): rules for PP3 and PP2 beside PM1 on the de novo
scale, trained in DD and tested in autism and the public cohort. Cell-level (saved de novo estimates); error
e = ln R + (1 - kappa) ln LR(PP3) in points; weighted MSE (1/SE^2); parametric bootstrap 10,000."""
import math
import re

import numpy as np
import pandas as pd

LN = math.log(2.406)
rng = np.random.default_rng(20261008)
o = pd.read_csv("results_fullstudy/a10_overlap.tsv", sep="\t")
NOM = {"PP3 Supporting only": 1, "PP3 Moderate only": 2, "PP3 >=Strong(3-4)": 4}
c = o[o.tool.isin(["AlphaMissense", "REVEL"]) & o.flag.isin(["PM1-like", "PM1-strict"]) & o.pp3_level.isin(NOM)].copy()
c = c[(c.n_dnv_PP3_and_F >= 50) & (c.R_lo95 > 0)].copy()
c["lnR"] = np.log(c.R)
c["se"] = (np.log(c.R_hi95) - np.log(c.R_lo95)) / 3.92
c["a"] = np.log(c.LR_PP3)
c["s"] = c.pp3_level.map(NOM)
FIXED = {"R0 uncapped": lambda s: 1.0, "R1 ClinGen cap": lambda s: min(1.0, 2.0 / s), "R2 PP3 <= Supporting": lambda s: 1.0 / s,
         "R3 one level lower": lambda s: {1: 0.0, 2: 0.5, 4: 0.5}[s], "R6 no PP3 credit": lambda s: 0.0}


def fit_k(d, lnR):
    y = lnR + d.a.to_numpy()
    w = 1 / d.se.to_numpy() ** 2
    return float((w * d.a * y).sum() / (w * d.a ** 2).sum())


def mse(d, lnR, kap):
    e = (lnR + (1 - kap) * d.a.to_numpy()) / LN
    w = 1 / d.se.to_numpy() ** 2
    return float((w * e ** 2).sum() / w.sum())


def evaluate(tool, lnR_all):
    d = c[c.tool == tool]
    idx = d.index
    lr = pd.Series(lnR_all, index=c.index).loc[idx].to_numpy()
    tr, te = (d.cohort == "DD").to_numpy(), (d.cohort != "DD").to_numpy()
    k4 = fit_k(d[tr], lr[tr])
    k5 = {s: fit_k(d[tr & (d.s == s).to_numpy()], lr[tr & (d.s == s).to_numpy()]) for s in sorted(d.s.unique())}
    kap = {n: d.s.map(f).to_numpy(float) for n, f in FIXED.items()}
    kap["R4 learned constant"] = np.full(len(d), k4)
    kap["R5 learned per band"] = d.s.map(k5).to_numpy(float)
    train = {n: mse(d[tr], lr[tr], k[tr]) for n, k in kap.items()}
    test = {n: mse(d[te], lr[te], k[te]) for n, k in kap.items()}
    return train, test, k4, k5, d


out = []
for tool in ["AlphaMissense", "REVEL"]:
    train, test, k4, k5, d = evaluate(tool, c.lnR.to_numpy())
    out.append(f"{tool}: cells train (DD) {int((d.cohort == 'DD').sum())}, test (autism, public) {int((d.cohort != 'DD').sum())}; "
               f"learned kappa constant {k4:.2f}; per band " + ", ".join(f"nominal {s}: {v:.2f}" for s, v in k5.items()))
    cand = ["R2 PP3 <= Supporting", "R3 one level lower", "R4 learned constant", "R5 learned per band"]
    best = min(cand, key=lambda n: train[n])
    boot = {n: [] for n in train}
    for _ in range(10000):
        draw = rng.normal(c.lnR.to_numpy(), c.se.to_numpy())
        trb, teb, *_ = evaluate(tool, draw)
        for n in train:
            boot[n].append(teb[n] - teb["R1 ClinGen cap"])
    out.append(f"  {'rule':22s} {'train MSE':>10s} {'test MSE':>10s}   test MSE - cap [95% CI]")
    for n in train:
        b = np.array(boot[n])
        lo, hi = np.quantile(b, [0.025, 0.975])
        out.append(f"  {n:22s} {train[n]:10.3f} {test[n]:10.3f}   {test[n] - test['R1 ClinGen cap']:+.3f} [{lo:+.3f}, {hi:+.3f}]"
                   + ("  <- selected on DD" if n == best else ""))
    b = np.array(boot[best])
    hi = np.quantile(b, 0.975)
    out.append(f"  H5a ({'primary' if tool == 'AlphaMissense' else 'secondary'}): {best} vs cap on test cohorts -> "
               f"{'BEATS THE CAP' if hi < 0 else 'does not beat the cap'} (upper bound {hi:+.3f})")
    te = d[d.cohort != "DD"]
    for r in te.itertuples():
        e_cap = (r.lnR + (1 - min(1, 2 / r.s)) * r.a) / LN
        e_b = (r.lnR + (1 - (k4 if best == 'R4 learned constant' else k5.get(r.s, k4) if best == 'R5 learned per band' else FIXED[best](r.s))) * r.a) / LN
        out.append(f"    test cell {r.cohort:7s} {r.flag:10s} {r.pp3_level:20s} n {r.n_dnv_PP3_and_F:4d}: error cap {e_cap:+.2f} pts, {best} {e_b:+.2f} pts")

# ---- H5b PP2 beside PM1
pp = {}
for line in open("results_fullstudy/pm1pp2_summary.txt"):
    m = re.match(r"(\w+)\s+PM1-like\s+x PP2-like: R ([\d.]+) \[([\d.]+), ([\d.]+)\]", line)
    if m:
        pp[m.group(1)] = (math.log(float(m.group(2))), (math.log(float(m.group(4))) - math.log(float(m.group(3)))) / 3.92)
lr_pp2 = o[(o.tool == "AlphaMissense") & (o.flag == "PP2-like") & (o.pp3_level == "PP3>=Supporting")].set_index("cohort").LR_F
out.append("H5b PP2 beside PM1-like (error = ln R + (1 - kappa_PP2) ln LR(PP2), points):")
kd = 1 + pp["DD"][0] / math.log(lr_pp2["DD"])
for nm, k in [("add PP2 (kappa 1)", 1.0), ("no PP2 beside PM1 (kappa 0)", 0.0), (f"learned in DD (kappa {kd:.2f})", kd)]:
    e = (pp["autism"][0] + (1 - k) * math.log(lr_pp2["autism"])) / LN
    out.append(f"  test autism: {nm:28s} error {e:+.2f} points (SE {pp['autism'][1] / LN:.2f})")

# ---- H5c triple
out.append("H5c triple PP3 x PM1 x PP2 (descriptive; error in points, test = autism):")
lrF = o[(o.tool == "AlphaMissense") & (o.pp3_level == "PP3>=Supporting")].set_index(["cohort", "flag"]).LR_F
lrP = o[(o.tool == "AlphaMissense") & (o.flag == "PM1-like")].set_index(["cohort", "pp3_level"]).LR_PP3
for line in open("results_fullstudy/triple_summary.txt"):
    m = re.match(r"(\w+)\s+AlphaMissense (\w+)\s+x (PM1-\w+)\s+x PP2-like: R3 ([\d.]+) .*?joint LR ([\d.]+)", line)
    if not m or m.group(1) != "autism":
        continue
    coh, band, flag, joint = m.group(1), m.group(2), m.group(3), float(m.group(5))
    s = 1 if band == "Supporting" else 4
    lp = lrP[(coh, "PP3>=Supporting" if band == "Supporting" else "PP3 >=Strong(3-4)")]
    l1, l2 = math.log(lrF[(coh, flag)]), math.log(lrF[(coh, "PP2-like")])
    for nm, k3, k2 in [("cap + add PP2", min(1, 2 / s), 1.0), ("cap + no PP2", min(1, 2 / s), 0.0),
                       ("Supporting limit + no PP2", 1 / s, 0.0), ("one level lower + no PP2", {1: 0, 4: 0.5}[s], 0.0),
                       ("uncapped", 1.0, 1.0)]:
        e = (math.log(joint) - (l1 + k3 * math.log(lp) + k2 * l2)) / LN
        out.append(f"  autism {band:10s} x {flag:10s}: {nm:26s} error {e:+.2f} points")
open("results_fullstudy/a20_denovo_rules_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
