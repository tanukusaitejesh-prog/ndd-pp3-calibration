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

# ---- review 13: sensitivity of the panel-applied credit, and a temporal out-of-sample test
d["date"] = pd.to_datetime(d["Approval Date"], errors="coerce")
FREQ = {"BA1", "BS1", "BS2", "PM2"}


def truthset(hi=6, lo=-1, drop_freq=False):
    ex = EXCL | FREQ if drop_freq else EXCL
    ind = ev.map(lambda e: sum(v for c, v in e.items() if c not in ex))
    t = d[(ind >= hi) | (ind <= lo)].copy()
    t["truth"] = (ind[t.index] >= hi).astype(int)
    return t


def boot_credit(t, sel=None, nb=1000):
    g = t.gene.unique()
    ix = t.groupby("gene").indices
    bs = [credit(t.iloc[np.concatenate([ix[x] for x in rng.choice(g, len(g))])], sel) for _ in range(nb)]
    return credit(t, sel), *np.nanpercentile(bs, [2.5, 97.5])


out.append("Sensitivity of PP3 credit beside panel PM1 (points [95% CI]; n joint; B/LB in joint):")
for lab, kw in (("truth >=6 / <=-1 (primary)", {}), ("truth >=4 / <=-1", {"hi": 4}), ("truth >=8 / <=-1", {"hi": 8}),
                ("truth >=6 / <=-3", {"lo": -3}), ("frequency codes also removed", {"drop_freq": True})):
    t = truthset(**kw)
    for sl, sel in (("any", None), ("Supporting", 1.0), ("Moderate", 2.0)):
        c0, lo, hi = boot_credit(t, sel)
        m = (t.PM1 > 0) & ((t.PP3 > 0) if sel is None else (t.PP3 == sel))
        out.append(f"  {lab:30s} PP3 {sl:10s}: {c0:.2f} [{lo:.2f}, {hi:.2f}]; n {int(m.sum())}; B/LB {int((m & (t.truth == 0)).sum())}")
cut = pd.Timestamp("2025-01-01")
tr, te = T[d.loc[T.index, "date"] < cut], T[d.loc[T.index, "date"] >= cut]
ctr = boot_credit(tr)
cte = boot_credit(te)
mte = (te.PM1 > 0) & (te.PP3 > 0)
out.append(f"Temporal split (approved before vs from 2025-01-01): training credit {ctr[0]:.2f} [{ctr[1]:.2f}, {ctr[2]:.2f}] "
           f"(n truth {len(tr)}); held-out credit {cte[0]:.2f} [{cte[1]:.2f}, {cte[2]:.2f}] (n truth {len(te)}; joint {int(mte.sum())}, "
           f"B/LB {int((mte & (te.truth == 0)).sum())}); held-out applied {te.loc[mte, 'PP3'].mean():.2f}")
later = d[d.date >= cut]
base = cls(later.pts.to_numpy())
for lab, newpts in (("PP3 at most Supporting when PM1 applies", later.pts - np.where((later.PM1 > 0) & (later.PP3 > 1), later.PP3 - 1, 0)),
                    ("PP3 at training-period credit when PM1 applies", later.pts - np.where((later.PM1 > 0) & (later.PP3 > 0), np.maximum(later.PP3 - max(ctr[0], 0), 0), 0)),
                    ("ClinGen cap", later.pts - (later.PP3 + later.PM1) + np.minimum(later.PP3 + later.PM1, 4))):
    ch = cls(newpts.to_numpy()) != base
    out.append(f"  2025-2026 classifications ({len(later)}): rule [{lab}] changes {int(ch.sum())}")
open("results_fullstudy_backup/stacking_cv_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out[-25:]))

# ---- review 15: Moderate - Supporting difference test with power; de novo criteria (PS2, PM6) excluded from truth
DENOVO = {"PS2", "PM6"}


def truthset2(drop):
    ex = EXCL | drop
    ind = ev.map(lambda e: sum(v for c, v in e.items() if c not in ex))
    t = d[(ind >= 6) | (ind <= -1)].copy()
    t["truth"] = (ind[t.index] >= 6).astype(int)
    return t


out.append("Moderate - Supporting difference in PP3 credit beside panel PM1 (gene-clustered bootstrap, 2,000):")
for lab, drop in (("primary truth", set()), ("PS2/PM6 also excluded", DENOVO), ("PS2/PM6 and BA1/BS1/BS2/PM2 excluded", DENOVO | FREQ)):
    t = truthset2(drop)
    g = t.gene.unique()
    ix = t.groupby("gene").indices
    vals = []
    for _ in range(2000):
        z = t.iloc[np.concatenate([ix[x] for x in rng.choice(g, len(g))])]
        vals.append((credit(z, 1.0), credit(z, 2.0)))
    V = np.array(vals)
    dlt = V[:, 1] - V[:, 0]
    d0 = credit(t, 2.0) - credit(t, 1.0)
    lo, hi = np.nanpercentile(dlt, [2.5, 97.5])
    se = np.nanstd(dlt)
    p_ = 2 * min((dlt >= 0).mean(), (dlt <= 0).mean())
    nS = int(((t.PM1 > 0) & (t.PP3 == 1)).sum()); nM = int(((t.PM1 > 0) & (t.PP3 == 2)).sum())
    bS = int(((t.PM1 > 0) & (t.PP3 == 1) & (t.truth == 0)).sum()); bM = int(((t.PM1 > 0) & (t.PP3 == 2) & (t.truth == 0)).sum())
    out.append(f"  {lab:38s}: Supporting {credit(t, 1.0):.2f}, Moderate {credit(t, 2.0):.2f}; difference {d0:+.2f} [{lo:+.2f}, {hi:+.2f}], "
               f"bootstrap p {p_:.3f}; n Supporting {nS} (B/LB {bS}), Moderate {nM} (B/LB {bM}); SE {se:.2f}; "
               f"minimum detectable difference (80% power, two-sided 5%) {2.8 * se:.2f} points; "
               f"overall credit {credit(t):.2f}; truth n {len(t)}")
open("results_fullstudy_backup/stacking_cv_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out[-4:]))
out.append("PP3 credit beside PM1 with PS2/PM6 excluded from truth (gene-clustered 95% CI):")
t = truthset2(DENOVO)
for lab, sel in (("any", None), ("Supporting", 1.0), ("Moderate", 2.0)):
    c0, lo, hi = boot_credit(t, sel)
    m = (t.PM1 > 0) & ((t.PP3 > 0) if sel is None else (t.PP3 == sel))
    out.append(f"  PP3 {lab:10s}: {c0:.2f} [{lo:.2f}, {hi:.2f}]; applied {t.loc[m, 'PP3'].mean():.2f}; n {int(m.sum())}")
open("results_fullstudy_backup/stacking_cv_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out[-4:]))
