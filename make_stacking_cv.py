"""Out-of-sample check of the expert-practice stacking result (review 12, descriptive).
(1) Strength-specific PP3 credit beside panel PM1 (PP3 Supporting vs PP3 Moderate), gene-clustered bootstrap.
(2) Gene-split cross-validation (200 random 50/50 gene splits): estimate PP3 credit beside PM1 on half A; on half B
    check (a) whether the held-out credit is below the points applied and (b) how many half-B classifications change
    when PP3 is credited at the half-A estimate.
Missense classifications from the ClinGen Evidence Repository export. Run: python make_stacking_cv.py <erepo.tsv>"""
import math, re, sys
import numpy as np, pandas as pd

d = pd.read_csv(sys.argv[1], sep="\t", dtype=str).fillna("")
d = d[d["Retracted"].str.lower() != "true"]
MIS = re.compile(r"\(p\.(?:[A-Z][a-z]{2})\d+(?!Ter|=|del|ins|dup|fs|ext)[A-Z][a-z]{2}\)")
d = d[d["#Variation"].str.contains(MIS)].reset_index(drop=True)
BASE = {"PVS": 8, "PS": 4, "PM": 2, "PP": 1, "BA": -8, "BS": -4, "BP": -1}
MOD = {"verystrong": 8, "strong": 4, "moderate": 2, "supporting": 1, "standalone": 8}
EXCL = {"PP3", "BP4", "PM1", "PP2", "BP1", "BP3", "PM4", "PM5"}


def codes(txt):
    ev = {}
    for c in re.split(r"[,;]\s*", txt):
        m = re.match(r"(PVS|PS|PM|PP|BA|BS|BP)(\d+)(?:_(\w+))?", c.strip())
        if not m:
            continue
        v = BASE[m.group(1)]
        if m.group(3) and m.group(3).lower().replace("_", "") in MOD:
            v = MOD[m.group(3).lower().replace("_", "")] * (1 if v > 0 else -1)
        ev[m.group(1) + m.group(2)] = v
    return ev


def cls(p):
    return np.select([p >= 10, p >= 6, p >= 0, p >= -6], ["P", "LP", "VUS", "LB"], "B")


ev = d["Applied Evidence Codes (Met)"].map(codes)
d["pts"] = ev.map(lambda e: sum(e.values())).astype(float)
for k in ("PP3", "PM1"):
    d[k] = ev.map(lambda e, k=k: e.get(k, 0)).astype(float)
d["indep"] = ev.map(lambda e: sum(v for c, v in e.items() if c not in EXCL))
d["gene"] = d["HGNC Gene Symbol"]
T = d[(d.indep >= 6) | (d.indep <= -1)].copy()
T["truth"] = (T.indep >= 6).astype(int)


def lr(m, y):
    return m[y == 1].mean() / max(m[y == 0].mean(), 1e-9)


def credit(z, sel=None):
    """Points of independent evidence earned by PP3 (optionally only PP3 at strength `sel`) once PM1 applies."""
    y, B = z.truth.to_numpy(), (z.PM1 > 0).to_numpy()
    A = (z.PP3 > 0).to_numpy() if sel is None else (z.PP3 == sel).to_numpy()
    lb = lr(B, y)
    return math.log(max(lr(A & B, y) / max(lb, 1e-9), 1e-9)) / math.log(2.406)


rng = np.random.default_rng(20261014)
genes = T.gene.unique()
gi = T.groupby("gene").indices
out = [f"Missense classifications {len(d)}; independent truth {len(T)} (P/LP {T.truth.sum()}, B/LB {(1 - T.truth).sum()}, genes {len(genes)})"]
for lab, sel in (("PP3 any strength", None), ("PP3 Supporting (1 point)", 1.0), ("PP3 Moderate (2 points)", 2.0)):
    c0 = credit(T, sel)
    bs = [credit(T.iloc[np.concatenate([gi[g] for g in rng.choice(genes, len(genes))])], sel) for _ in range(1000)]
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    m = (T.PM1 > 0) & ((T.PP3 > 0) if sel is None else (T.PP3 == sel))
    out.append(f"{lab:26s} beside panel PM1: earns {c0:.2f} points [{lo:.2f}, {hi:.2f}]; applied {T.loc[m, 'PP3'].mean():.2f}; "
               f"n {int(m.sum())} (B/LB {int((m & (T.truth == 0)).sum())})")
# gene-split cross-validation
res = []
for _ in range(200):
    A = set(rng.choice(genes, len(genes) // 2, replace=False))
    tA, tB = T[T.gene.isin(A)], T[~T.gene.isin(A)]
    cA, cB = credit(tA), credit(tB)
    dB = d[~d.gene.isin(A)]
    base = cls(dB.pts.to_numpy())
    corr = dB.pts - np.where((dB.PM1 > 0) & (dB.PP3 > 0), np.maximum(dB.PP3 - max(cA, 0), 0), 0)
    ch = cls(corr.to_numpy()) != base
    mB = (tB.PM1 > 0) & (tB.PP3 > 0)
    res.append((cA, cB, tB.loc[mB, "PP3"].mean(), int(ch.sum()), len(dB)))
R = pd.DataFrame(res, columns=["credit_A", "credit_B", "applied_B", "changed_B", "n_B"])
out.append(f"Gene-split CV (200 splits): training-half credit median {R.credit_A.median():.2f} (IQR {R.credit_A.quantile(.25):.2f}-"
           f"{R.credit_A.quantile(.75):.2f}); held-out credit median {R.credit_B.median():.2f} (IQR {R.credit_B.quantile(.25):.2f}-"
           f"{R.credit_B.quantile(.75):.2f}); held-out credit below points applied in {(R.credit_B < R.applied_B).mean():.0%} of splits")
out.append(f"Held-out classifications changed when PP3 is credited at the training-half estimate: median {R.changed_B.median():.0f} "
           f"of {R.n_B.median():.0f} (2.5-97.5%: {R.changed_B.quantile(.025):.0f}-{R.changed_B.quantile(.975):.0f}); scaled to all "
           f"classifications: median {(R.changed_B / R.n_B).median() * len(d):.0f}")
open("results_fullstudy_backup/stacking_cv_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
