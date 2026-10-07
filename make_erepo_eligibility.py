"""Selection into PP3 on expert-panel truth (review 15). Among Evidence Repository missense classifications mapped to
analysed sites, PP3 eligibility is computed from scores (REVEL >= 0.644 or AlphaMissense >= 0.792, the ClinGen Supporting
thresholds). Within classifications where the panel applied PM1, we compare PP3-applied, eligible-not-applied and
ineligible groups on independent truth, and estimate a selection-free credit: points earned by score-defined PP3 once
panel PM1 applies. Truth excludes PP3, BP4, PM1, PP2, BP1, BP3, PM4, PM5 and the de novo criteria PS2, PM6.
Run in the WSL repo: python make_erepo_eligibility.py <erepo.tsv>"""
import math, re, sys
import numpy as np, pandas as pd

sys.path.insert(0, ".")
import gate

d = pd.read_csv(sys.argv[1], sep="\t", dtype=str).fillna("")
d = d[d["Retracted"].str.lower() != "true"]
BASE = {"PVS": 8, "PS": 4, "PM": 2, "PP": 1, "BA": -8, "BS": -4, "BP": -1}
MOD = {"verystrong": 8, "strong": 4, "moderate": 2, "supporting": 1, "standalone": 8}
EXCL = {"PP3", "BP4", "PM1", "PP2", "BP1", "BP3", "PM4", "PM5", "PS2", "PM6"}


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


rows = []
for r in d.itertuples(index=False):
    best = {}
    for x in re.finditer(r"NC_0000(\d\d)\.(\d+):g\.(\d+)([ACGT])>([ACGT])", r[3]):
        ch, ver = int(x.group(1)), int(x.group(2))
        if ch not in best or ver > best[ch][0]:
            best[ch] = (ver, int(x.group(3)), x.group(5))
    if len(best) != 1:
        continue
    ch, (ver, pos, alt) = next(iter(best.items()))
    ev = codes(r[9])
    ind = sum(v for c, v in ev.items() if c not in EXCL)
    rows.append({"chrom": str(ch), "pos": pos, "alt": alt, "gene": r[4], "ind": ind, "PP3": ev.get("PP3", 0), "PM1": ev.get("PM1", 0)})
x = pd.DataFrame(rows)
x = x[x.chrom.isin(gate.AUTOSOMES)]
x["key"] = x.pos.astype(np.int64) * 4 + x.alt.map(gate.B2I)
parts = []
for c, g in x.groupby("chrom"):
    od = f"work/fullstudy/col/chr{c}"
    k = pd.Index(np.load(f"{od}/key.npy"))
    ix = k.get_indexer(g.key.to_numpy())
    g = g[ix >= 0].copy(); ix = ix[ix >= 0]
    g["am"] = np.load(f"{od}/am.npy", mmap_mode="r")[ix]
    g["revel"] = np.load(f"{od}/raw_REVEL_score.npy", mmap_mode="r")[ix]
    parts.append(g)
x = pd.concat(parts, ignore_index=True)
t = x[(x.ind >= 6) | (x.ind <= -1)].copy()
t["truth"] = (t.ind >= 6).astype(int)
t["elig"] = (t.revel >= 0.644) | (t.am >= 0.792)
out = [f"Mapped missense classifications {len(x)}; with truth (PS2/PM6 also excluded) {len(t)} (P/LP {t.truth.sum()}, B/LB {(1 - t.truth).sum()}, genes {t.gene.nunique()})"]
pm = t[t.PM1 > 0]
y = pm.truth.to_numpy()


def lr(m, yy):
    return m[yy == 1].mean() / max(m[yy == 0].mean(), 1e-9)


grp = {"PP3 applied": (pm.PP3 > 0).to_numpy(),
       "eligible, PP3 not applied": (pm.elig & (pm.PP3 <= 0)).to_numpy(),
       "not eligible, PP3 not applied": (~pm.elig & (pm.PP3 <= 0)).to_numpy()}
out.append(f"Within panel-PM1 classifications (n {len(pm)}; P/LP {int(y.sum())}, B/LB {int((1 - y).sum())}):")
for k_, m_ in grp.items():
    out.append(f"  {k_:32s}: n {int(m_.sum())} (P/LP {int(m_[y == 1].sum())}, B/LB {int(m_[y == 0].sum())}); LR {lr(m_, y):.2f}")
rng = np.random.default_rng(20261017)
genes = t.gene.unique(); gi = t.groupby("gene").indices


def credit(z, flag):
    yy, B = z.truth.to_numpy(), (z.PM1 > 0).to_numpy()
    A = flag(z)
    return math.log(max(lr(A & B, yy) / max(lr(B, yy), 1e-9), 1e-9)) / math.log(2.406)


for lab, flag in (("panel-applied PP3", lambda z: (z.PP3 > 0).to_numpy()),
                  ("score-defined PP3 (REVEL >= 0.644 or AM >= 0.792)", lambda z: z.elig.to_numpy()),
                  ("score-defined PP3 (REVEL >= 0.644)", lambda z: (z.revel >= 0.644).to_numpy()),
                  ("score-defined PP3 (AM >= 0.792)", lambda z: (z.am >= 0.792).to_numpy())):
    c0 = credit(t, flag)
    bs = [credit(t.iloc[np.concatenate([gi[g] for g in rng.choice(genes, len(genes))])], flag) for _ in range(1000)]
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    m_ = flag(t) & (t.PM1 > 0).to_numpy()
    out.append(f"credit beside panel PM1, {lab:52s}: {c0:.2f} points [{lo:.2f}, {hi:.2f}]; n {int(m_.sum())}; B/LB {int((m_ & (t.truth == 0).to_numpy()).sum())}")
open("results_fullstudy/erepo_eligibility_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))

# ---- review 16: score-defined PP3 strength bands beside panel PM1 (no curator choice)
BANDS = {"AlphaMissense": ("am", [(0.792, 0.906, "Supporting", 1), (0.906, 0.972, "Moderate", 2), (0.972, 9, "Strong", 4)]),
         "REVEL": ("revel", [(0.644, 0.773, "Supporting", 1), (0.773, 0.932, "Moderate", 2), (0.932, 9, "Strong", 4)])}
out.append("Score-defined PP3 bands beside panel PM1 (credit = points earned by the band, relative to all other PM1 classifications):")
for tool, (col, bands) in BANDS.items():
    res = {}
    for lo_, hi_, nm, pts in bands:
        flag = lambda z, lo_=lo_, hi_=hi_: ((z[col] >= lo_) & (z[col] < hi_)).to_numpy()
        c0 = credit(t, flag)
        bs = np.array([credit(t.iloc[np.concatenate([gi[g] for g in rng.choice(genes, len(genes))])], flag) for _ in range(1000)])
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        res[nm] = bs
        m_ = flag(t) & (t.PM1 > 0).to_numpy()
        out.append(f"  {tool:13s} {nm:10s} band (nominal {pts}): earns {c0:.2f} points [{lo:.2f}, {hi:.2f}]; n {int(m_.sum())}; "
                   f"B/LB {int((m_ & (t.truth == 0).to_numpy()).sum())}")
    dl = res["Moderate"] - res["Supporting"]
    out.append(f"  {tool:13s} Moderate - Supporting band: {np.nanmedian(dl):+.2f} [{np.nanpercentile(dl, 2.5):+.2f}, {np.nanpercentile(dl, 97.5):+.2f}]")
open("results_fullstudy/erepo_eligibility_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out[-9:]))
