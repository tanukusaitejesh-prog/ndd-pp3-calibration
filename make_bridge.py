"""Estimand bridge: AlphaMissense PP3 LR and R(PP3 x PM1-like) on the ClinVar scale (P/LP vs B/LB), DD de novo sites removed,
all genes and DDG2P monoallelic; gene-clustered bootstrap (1,000, seed 20261007). Run inside the WSL repo."""
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.getcwd())
import fullstudy as fs, gate

a = argparse.Namespace(data="data", work="work", fs="work/fullstudy", results="results_fullstudy")
m = pd.read_csv("results_fullstudy/a4_clinvar_matched.tsv.gz", sep="\t", dtype={"chrom": str})
m = m[m["cls"].isin(["P/LP", "B/LB"])].copy()
meta = pd.read_csv("work/dnv_grch38.tsv", sep="\t", dtype={"chrom": str}, usecols=["chrom", "pos38", "alt38"])
meta = meta[meta["alt38"].isin(list(gate.B2I))]
dk = set(zip(meta["chrom"].str.replace("chr", ""), meta["pos38"].astype(np.int64) * 4 + meta["alt38"].map(gate.B2I)))
m = m[[(c.replace("chr", ""), int(k)) not in dk for c, k in zip(m["chrom"], m["key"])]]
tx = {}
for c in gate.AUTOSOMES:
    od = f"work/fullstudy/col/chr{c}"
    tx[c] = pd.Series(np.load(f"{od}/tx.npy"), index=np.load(f"{od}/key.npy"))
m["tx2"] = [tx[c.replace("chr", "")].get(k, -1) if c.replace("chr", "") in tx else -1 for c, k in zip(m["chrom"], m["key"])]
ti = fs.load_tinfo(a)
m = m[m["tx2"] >= 0]
m["gene"] = ti["gene"].reindex(m["tx2"].to_numpy()).to_numpy()
m["dd"] = ti["ddg2p"].reindex(m["tx2"].to_numpy()).fillna(False).to_numpy(bool)
m["F"] = m["pm1like"].astype(str) == "True"
m["y"] = (m["cls"] == "P/LP").astype(int)
cuts = [0.792, 0.906, 0.972]
rng = np.random.default_rng(20261007)
out = []


def stats(x):
    P, B = x[x.y == 1], x[x.y == 0]
    r = {}
    for nm, cut in (("Supp", cuts[0]), ("Strong", cuts[2])):
        pp = x.am >= cut
        def lr(mask):
            return (mask[x.y == 1].mean()) / max(mask[x.y == 0].mean(), 1e-9)
        lrP, lrF, lrJ = lr(pp), lr(x.F), lr(pp & x.F)
        r[nm] = (lrP, lrF, lrJ, lrJ / (lrP * lrF))
    return r


for name, sub in (("all_genes", m), ("ddg2p_monoallelic", m[m.dd])):
    sub = sub.dropna(subset=["am"])
    st = stats(sub)
    genes = sub["gene"].unique()
    g_ix = {g: ix for g, ix in sub.groupby("gene").indices.items()}
    boots = {k: [] for k in st}
    for _ in range(1000):
        pick = rng.choice(genes, len(genes))
        idx = np.concatenate([g_ix[g] for g in pick])
        b = stats(sub.iloc[idx])
        for k in st:
            boots[k].append(b[k][3])
    for k, (lrP, lrF, lrJ, R) in st.items():
        lo, hi = np.nanpercentile(boots[k], [2.5, 97.5])
        out.append(f"ClinVar {name:18s} AM PP3>={k:6s} LR(PP3) {lrP:.1f} LR(PM1-like) {lrF:.2f} joint {lrJ:.1f} "
                   f"R {R:.2f} [{lo:.2f}, {hi:.2f}]  (P/LP {int(sub.y.sum())}, B/LB {int((1 - sub.y).sum())}, genes {len(genes)})")
open("results_fullstudy/bridge_clinvar_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
