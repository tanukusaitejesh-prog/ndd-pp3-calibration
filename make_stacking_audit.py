"""Stacking audit on ClinGen Evidence Repository classifications (review 11, descriptive).
(1) How often panels apply PP3, PM1 and PP2 together; (2) classification changes if location-type evidence is capped;
(3) dependence of panel-applied PM1, PP2 and PP3 on independent expert truth (points excluding PP3, BP4, PM1, PP2, BP1,
BP3, PM4, PM5), gene-clustered bootstrap. Run: python make_stacking_audit.py work_erepo/erepo_2026-10-07.tsv"""
import math, re, sys
import numpy as np, pandas as pd

d = pd.read_csv(sys.argv[1], sep="\t", dtype=str).fillna("")
d = d[d["Retracted"].str.lower() != "true"]
MIS = re.compile(r"\(p\.(?:[A-Z][a-z]{2})\d+(?!Ter|=|del|ins|dup|fs|ext)[A-Z][a-z]{2}\)")
if len(sys.argv) > 2 and sys.argv[2] == "missense":
    d = d[d["#Variation"].str.contains(MIS)]
TAG = sys.argv[2] if len(sys.argv) > 2 else "all"
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
    return "P" if p >= 10 else "LP" if p >= 6 else "VUS" if p >= 0 else "LB" if p >= -6 else "B"


ev = d["Applied Evidence Codes (Met)"].map(codes)
d["pts"] = ev.map(lambda e: sum(e.values()))
for k in ("PP3", "PM1", "PP2"):
    d[k] = ev.map(lambda e, k=k: e.get(k, 0))
d["indep"] = ev.map(lambda e: sum(v for c, v in e.items() if c not in EXCL))
out = [f"Classifications {len(d)} (retracted removed); panels {d['Expert Panel'].nunique()}"]
f = {k: d[k] > 0 for k in ("PP3", "PM1", "PP2")}
out.append(f"applied: PP3 {f['PP3'].sum()}, PM1 {f['PM1'].sum()}, PP2 {f['PP2'].sum()}; PP3+PM1 {(f['PP3'] & f['PM1']).sum()}, "
           f"PM1+PP2 {(f['PM1'] & f['PP2']).sum()}, PP3+PP2 {(f['PP3'] & f['PP2']).sum()}, all three {(f['PP3'] & f['PM1'] & f['PP2']).sum()}")
stack = d[f["PM1"] & f["PP2"]]
out.append("PM1+PP2 by panel: " + "; ".join(f"{k} {v}" for k, v in stack["Expert Panel"].value_counts().items()))
out.append(f"PM1+PP2 classifications with points from PP3+PM1+PP2 > 4: {((stack.PP3 + stack.PM1 + stack.PP2) > 4).sum()} of {len(stack)}")
loc = d.PP3 + d.PM1 + d.PP2
rules = {
    "as applied": d.pts,
    "ClinGen cap only (PP3+PM1 <= 4)": d.pts - (d.PP3 + d.PM1) + np.minimum(d.PP3 + d.PM1, 4),
    "PM1+PP2 <= PM1 (PP2 dropped when PM1 applies)": d.pts - np.where((d.PM1 > 0) & (d.PP2 > 0), d.PP2, 0),
    "location cap: PP3+PM1+PP2 <= 4": d.pts - loc + np.minimum(loc, 4),
    "PP3 credited at measured 0.6 points when PM1 applies": d.pts - np.where((d.PM1 > 0) & (d.PP3 > 0), d.PP3 - 0.6, 0),
    "PP3 credited at 1.10 (upper 95% bound) when PM1 applies": d.pts - np.where((d.PM1 > 0) & (d.PP3 > 0), np.maximum(d.PP3 - 1.10, 0), 0),
    "PP3 credited at 0.11 (lower 95% bound) when PM1 applies": d.pts - np.where((d.PM1 > 0) & (d.PP3 > 0), d.PP3 - 0.11, 0),
}
base = d.pts.map(cls)
for name, p in rules.items():
    c = p.map(cls)
    ch = (c != base)
    tr = pd.crosstab(base[ch], c[ch]) if ch.any() else None
    moves = "; ".join(f"{a}->{b} {int(tr.loc[a, b])}" for a in tr.index for b in tr.columns if tr.loc[a, b] > 0) if tr is not None else "none"
    out.append(f"rule [{name}]: classifications changed {int(ch.sum())} ({moves})")
t = d[(d.indep >= 6) | (d.indep <= -1)].copy()
t["truth"] = (t.indep >= 6).astype(int)
out.append(f"independent expert truth: {len(t)} classifications (P/LP {t.truth.sum()}, B/LB {(1 - t.truth).sum()}, genes {t['HGNC Gene Symbol'].nunique()})")
rng = np.random.default_rng(20261012)
genes = t["HGNC Gene Symbol"].to_numpy()
ug = np.unique(genes)
gi = pd.Series(np.arange(len(t))).groupby(genes).indices


def lr(m, y):
    return m[y == 1].mean() / max(m[y == 0].mean(), 1e-9)


def R2(z, a, b):
    y, A, B = z.truth.to_numpy(), (z[a] > 0).to_numpy(), (z[b] > 0).to_numpy()
    la, lb, lj = lr(A, y), lr(B, y), lr(A & B, y)
    return lj / (la * lb) if la * lb > 0 else np.nan, lj, la, lb


def R3(z):
    y = z.truth.to_numpy()
    A, B, C = ((z[k] > 0).to_numpy() for k in ("PP3", "PM1", "PP2"))
    l = [lr(A, y), lr(B, y), lr(C, y)]
    lj = lr(A & B & C, y)
    return lj / np.prod(l) if np.prod(l) > 0 else np.nan, lj, l
for a, b in (("PM1", "PP2"), ("PP3", "PM1"), ("PP3", "PP2")):
    r0, lj, la, lb = R2(t, a, b)
    bs = [R2(t.iloc[np.concatenate([gi[g] for g in rng.choice(ug, len(ug))])], a, b)[0] for _ in range(1000)]
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    nj = int(((t[a] > 0) & (t[b] > 0)).sum())
    out.append(f"expert truth R({a},{b}) {r0:.2f} [{lo:.2f}, {hi:.2f}]; LR {a} {la:.2f}, {b} {lb:.2f}, joint {lj:.2f}; n joint {nj} "
               f"(B/LB in joint {int(((t[a] > 0) & (t[b] > 0) & (t.truth == 0)).sum())})")
def pts(z):
    y, A, B = z.truth.to_numpy(), (z.PP3 > 0).to_numpy(), (z.PM1 > 0).to_numpy()
    return math.log(max(lr(A & B, y) / max(lr(B, y), 1e-9), 1e-9)) / math.log(2.406)
p0 = pts(t)
bs = [pts(t.iloc[np.concatenate([gi[g] for g in rng.choice(ug, len(ug))])]) for _ in range(1000)]
lo, hi = np.nanpercentile(bs, [2.5, 97.5])
out.append(f"points earned by panel-applied PP3 (any strength) once panel PM1 applies: {p0:.2f} [{lo:.2f}, {hi:.2f}]; "
           f"mean PP3 points applied in that cell {t.loc[(t.PP3 > 0) & (t.PM1 > 0), 'PP3'].mean():.2f}")
r0, lj, l = R3(t)
bs = [R3(t.iloc[np.concatenate([gi[g] for g in rng.choice(ug, len(ug))])])[0] for _ in range(1000)]
lo, hi = np.nanpercentile(bs, [2.5, 97.5])
out.append(f"expert truth R3(PP3,PM1,PP2) {r0:.3f} [{lo:.3f}, {hi:.3f}]; joint LR {lj:.2f} vs product {np.prod(l):.2f}; "
           f"summed points {sum(math.log(x) for x in l) / math.log(2.406):.2f} vs joint points {math.log(max(lj, 1e-9)) / math.log(2.406):.2f}")
open(f"results_fullstudy_backup/stacking_audit_{TAG}_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
