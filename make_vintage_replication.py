"""A19 H4 (pre-registered): replication across ClinVar vintages and evaluation eras, plus negative-control co-criteria.
SCVs pooled over vintages 2020-01, 2022-01, 2024-01, 2026-10; each SCV accession (version stripped) counted once at its
earliest vintage carrying codes. Eras by DateLastEvaluated. Per era: R(PP3, PM1) = LR(joint)/(LR(PP3) LR(PM1)) on truth,
and ACMG-convention PP2 credit beside PM1; gene-clustered bootstrap (1,000).
Negative control: PP1 and PP4 credit beside vs without PM1, with truth also excluding PP1 and PP4 (so the controls are not
part of their own truth)."""
import math

import numpy as np
import pandas as pd

LN = math.log(2.406)
rng = np.random.default_rng(20261008)
cols = ["SCV", "Submitter", "gene", "PP3", "PM1", "PP2", "PP1", "PP4", "indep", "DateLastEvaluated", "vintage"]
S = pd.concat([pd.read_csv(f"work/cvt/lab_codes_{v}.tsv.gz", sep="\t", usecols=cols) for v in ["2020-01", "2022-01", "2024-01", "2026-10"]])
S["acc"] = S.SCV.str.split(".").str[0]
S = S.sort_values("vintage").drop_duplicates("acc", keep="first").reset_index(drop=True)
S["year"] = pd.to_datetime(S.DateLastEvaluated, errors="coerce").dt.year
S["era"] = pd.cut(S.year, [0, 2019, 2022, 2024, 2100], labels=["<=2019", "2020-2022", "2023-2024", "2025-2026"])
for c in ["PP3", "PM1", "PP2", "PP1", "PP4"]:
    S[c + "_on"] = (S[c] > 0).astype(int)


def truth(df, extra_excl=()):
    pts = df.indep - sum(df[c].clip(lower=0) for c in extra_excl) if extra_excl else df.indep
    d = df[(pts >= 6) | (pts <= -1)].copy()
    d["y"] = (pts[d.index] >= 6).astype(int)
    return d


def lr(d, mask):
    p, b = d[d.y == 1], d[d.y == 0]
    return ((mask(p).sum() + 0.5) / (len(p) + 1)) / ((mask(b).sum() + 0.5) / (len(b) + 1))


def R(d):
    return lr(d, lambda x: (x.PP3_on == 1) & (x.PM1_on == 1)) / (lr(d, lambda x: x.PP3_on == 1) * lr(d, lambda x: x.PM1_on == 1))


def credit(d, code, pm1):
    e = d[d.PM1_on == pm1]
    a = e[e[code + "_on"] == 1]
    o1 = (a.y.sum() + 0.5) / ((1 - a.y).sum() + 0.5)
    o0 = (e.y.sum() + 0.5) / ((1 - e.y).sum() + 0.5)
    return math.log(o1 / o0) / LN


def ci(d, f, n=1000):
    g = d.groupby("gene").indices
    keys = list(g)
    v = []
    for _ in range(n):
        pick = rng.choice(len(keys), len(keys))
        b = d.iloc[np.concatenate([g[keys[i]] for i in pick])]
        try:
            v.append(f(b))
        except Exception:
            pass
    return np.nanquantile(v, [0.025, 0.975])


T = truth(S)
out = [f"A19 H4: pooled unique SCVs {len(S)} (by first vintage: " + ", ".join(f"{k} {v}" for k, v in S.vintage.value_counts().sort_index().items())
       + f"); with truth {len(T)}"]
rows = []
for era, d in T.groupby("era", observed=True):
    nb_joint = int(((d.PP3_on == 1) & (d.PM1_on == 1) & (d.y == 0)).sum())
    nb_pp2 = int(((d.PP2_on == 1) & (d.PM1_on == 1) & (d.y == 0)).sum())
    r, (rlo, rhi) = R(d), ci(d, R)
    c, (clo, chi) = credit(d, "PP2", 1), ci(d, lambda b: credit(b, "PP2", 1))
    rows.append(dict(era=era, n=len(d), R=r, R_lo=rlo, R_hi=rhi, B_joint=nb_joint, pp2=c, pp2_lo=clo, pp2_hi=chi, B_pp2=nb_pp2))
    out.append(f"  era {era:9s} n {len(d):6d}: R(PP3,PM1) {r:.2f} [{rlo:.2f}, {rhi:.2f}] (B/LB in joint {nb_joint}); "
               f"PP2 credit beside PM1 {c:+.2f} [{clo:+.2f}, {chi:+.2f}] (B/LB {nb_pp2})")
E = pd.DataFrame(rows)
okR = E[E.B_joint >= 20]
okP = E[E.B_pp2 >= 20]
out.append(f"Replication rule (pre-specified): R upper bound < 1 in every era with >= 20 B/LB -> "
           f"{'REPLICATED' if (okR.R_hi < 1).all() and len(okR) else 'not replicated'} ({len(okR)} eligible eras); "
           f"PP2 credit < 1 in every eligible era -> {'REPLICATED' if (okP.pp2 < 1).all() and len(okP) else 'not replicated'} ({len(okP)} eligible eras)")
# negative controls
T2 = truth(S, extra_excl=("PP1", "PP4"))
out.append(f"Negative control (truth also excluding PP1, PP4; n {len(T2)}): credit beside PM1 vs without PM1 (gene-clustered 95% CI)")
for code in ["PP1", "PP4", "PP2", "PP3"]:
    a, (alo, ahi) = credit(T2, code, 1), ci(T2, lambda b, c=code: credit(b, c, 1))
    w, (wlo, whi) = credit(T2, code, 0), ci(T2, lambda b, c=code: credit(b, c, 0))
    n1 = int(((T2[code + "_on"] == 1) & (T2.PM1_on == 1)).sum())
    out.append(f"  {code}: beside PM1 {a:+.2f} [{alo:+.2f}, {ahi:+.2f}] (n {n1}); without PM1 {w:+.2f} [{wlo:+.2f}, {whi:+.2f}]")
E.to_csv("results_fullstudy/a19_h4_eras.tsv", sep="\t", index=False)
open("results_fullstudy/a19_h4_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))

# R(code, PM1) for each co-criterion on the same truth (PP1/PP4 excluded from truth): specificity of the dependence
out.append("Specificity: R(code, PM1) with the negative-control truth (gene-clustered 95% CI)")
for code in ["PP1", "PP4", "PP2", "PP3"]:
    f = lambda b, c=code: lr(b, lambda x: (x[c + "_on"] == 1) & (x.PM1_on == 1)) / (lr(b, lambda x: x[c + "_on"] == 1) * lr(b, lambda x: x.PM1_on == 1))
    r, (lo, hi) = f(T2), ci(T2, f)
    out.append(f"  R({code}, PM1) {r:.2f} [{lo:.2f}, {hi:.2f}]")
    print(out[-1], flush=True)
for code in ["PP1", "PP4"]:
    f = lambda b, c=code: (lr(b, lambda x: (x.PP3_on == 1) & (x.PM1_on == 1)) / (lr(b, lambda x: x.PP3_on == 1) * lr(b, lambda x: x.PM1_on == 1))) / \
        (lr(b, lambda x: (x[c + "_on"] == 1) & (x.PM1_on == 1)) / (lr(b, lambda x: x[c + "_on"] == 1) * lr(b, lambda x: x.PM1_on == 1)))
    r, (lo, hi) = f(T2), ci(T2, f)
    out.append(f"  ratio R(PP3,PM1)/R({code},PM1) {r:.2f} [{lo:.2f}, {hi:.2f}] (< 1 = PP3 dependence beyond the control)")
    print(out[-1], flush=True)
open("results_fullstudy/a19_h4_summary.txt", "w").write("\n".join(out) + "\n")
