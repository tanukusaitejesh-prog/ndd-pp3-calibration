"""Task 6: PP3 x PM1 overlap and rule calibration on ClinGen expert-panel truth that excludes computational and
position evidence. Truth class = Tavtigian points from all Met codes EXCEPT PP3, BP4, PM1, PP2, BP1, BP3, PM4, PM5
(PM5 and PM1 are position-based), from the eRepo export. AlphaMissense and the PM1-like flag come from the column store.
Gene-clustered bootstrap, 1,000, seed 20261008. Run in the WSL repo: python make_erepo_truth.py <erepo.tsv>"""
import math, re, sys
import numpy as np, pandas as pd

sys.path.insert(0, ".")
import gate

src = sys.argv[1]
d = pd.read_csv(src, sep="\t", dtype=str).fillna("")
BASE = {"PVS": 8, "PS": 4, "PM": 2, "PP": 1, "BA": -8, "BS": -4, "BP": -1}
MOD = {"verystrong": 8, "strong": 4, "moderate": 2, "supporting": 1}
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


rows = []
pat = re.compile(r"NC_0000(\d\d)\.\d+:g\.(\d+)([ACGT])>([ACGT])")
for r in d.itertuples(index=False):
    hg = r[3]
    m = [x for x in pat.findall(hg)]
    if not m:
        continue
    # GRCh38 RefSeq accession versions: chr1 .11, others .10-.12 ; take the accession whose version is the highest per chrom
    best = {}
    for x in re.finditer(r"NC_0000(\d\d)\.(\d+):g\.(\d+)([ACGT])>([ACGT])", hg):
        ch, ver = int(x.group(1)), int(x.group(2))
        if ch not in best or ver > best[ch][0]:
            best[ch] = (ver, int(x.group(3)), x.group(5))
    if len(best) != 1:
        continue
    ch, (ver, pos, alt) = next(iter(best.items()))
    ev = codes(r[9])
    pts = sum(v for k, v in ev.items() if k not in EXCL)
    truth = 1 if pts >= 6 else (0 if pts <= -1 else -1)          # independent LP/P vs LB/B; else undetermined
    rows.append({"chrom": str(ch), "pos": pos, "alt": alt, "gene": r[4], "truth": truth, "pts_indep": pts,
                 "PP3_met": "PP3" in ev, "PM1_met": "PM1" in ev, "assertion": r[8]})
x = pd.DataFrame(rows)
x = x[x["chrom"].isin(gate.AUTOSOMES)]
x["key"] = x["pos"].astype(np.int64) * 4 + x["alt"].map(gate.B2I)
parts = []
for c, g in x.groupby("chrom"):
    od = f"work/fullstudy/col/chr{c}"
    k = pd.Index(np.load(f"{od}/key.npy"))
    ix = k.get_indexer(g["key"].to_numpy())
    g = g[ix >= 0].copy(); ix = ix[ix >= 0]
    for col in ("am", "rmc_oe", "ccr_pct"):
        g[col] = np.load(f"{od}/{col}.npy", mmap_mode="r")[ix]
    parts.append(g)
x = pd.concat(parts, ignore_index=True)
x["F"] = (x["rmc_oe"] < 0.4) | (x["ccr_pct"] >= 95)
x.to_csv("results_fullstudy/erepo_truth_sites.tsv", sep="\t", index=False)
t = x[(x["truth"] >= 0) & x["am"].notna()]
out = [f"eRepo missense mapped to analysed sites: {len(x)}; with independent truth (pts excl. PP3/BP4/PM1/PP2/BP1/BP3/PM4/PM5 "
       f">= 6 or <= -1): {len(t)} (P/LP {int(t.truth.sum())}, B/LB {int((1 - t.truth).sum())}, genes {t.gene.nunique()})"]
rng = np.random.default_rng(20261008)


def lrs(z, cut, flag):
    def lr(mk):
        return mk[z.truth == 1].mean() / max(mk[z.truth == 0].mean(), 1e-9)
    pp = z.am >= cut
    f = z[flag]
    return lr(pp), lr(f), lr(pp & f), lr(pp & ~f)


for flag, lab in (("F", "PM1-like (o/e, CCR)"), ("PM1_met", "PM1 applied by panel")):
    for cut, tier in ((0.792, "Supporting"), (0.972, "Strong")):
        lp, lf, lj, lnf = lrs(t, cut, flag)
        R = lj / (lp * lf)
        genes = t.gene.unique(); gi = t.groupby("gene").indices
        bs = []
        for _ in range(1000):
            idx = np.concatenate([gi[g] for g in rng.choice(genes, len(genes))])
            b = lrs(t.iloc[idx], cut, flag)
            bs.append(b[2] / (b[0] * b[1]) if b[0] * b[1] > 0 else np.nan)
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        pts_given = math.log(max(lj / lf, 1e-9)) / math.log(2.406)
        out.append(f"{lab:22s} AM>={tier:10s} LR(PP3) {lp:.2f} LR(F) {lf:.2f} joint {lj:.2f} R {R:.2f} [{lo:.2f}, {hi:.2f}]; "
                   f"PP3 given F {lj / lf:.2f} ({pts_given:.2f} points) vs not F {lnf:.2f}; n(PP3&F) {int(((t.am >= cut) & t[flag]).sum())}")
open("results_fullstudy/erepo_truth_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
