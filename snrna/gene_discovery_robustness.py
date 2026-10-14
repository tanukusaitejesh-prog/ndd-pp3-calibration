# Robustness of gene_discovery.py rel_z: window size, per-gene site-permutation and circular-shift nulls, AN in window.
# Needs model_inputs_v3 (python v3_counts.py cohorts/ukb/*.csv). Output: results/gene_discovery_robustness.txt
import sys, numpy as np, pandas as pd
sys.path.insert(0, ".")
from gene_discovery import fit_k, DISEASE_GENES
rng = np.random.default_rng(1)
NDD = [g for g, (_, c) in DISEASE_GENES.items() if c == "NDD"]


def minz(n, obs, p, w):
    pos = np.arange(n.min(), n.max() + 1)
    i0 = n - n.min()
    O = np.bincount(i0, obs, len(pos)); E = np.bincount(i0, p, len(pos)); V = np.bincount(i0, p * (1 - p), len(pos))
    c = lambda x: np.concatenate([[0], np.cumsum(x)])
    cO, cE, cV = c(O), c(E), c(V)
    w = min(w, len(pos))
    k = np.arange(len(pos) - w + 1)
    o, e, v = cO[k + w] - cO[k], cE[k + w] - cE[k], cV[k + w] - cV[k]
    z = np.where(e >= 2, (o - e) / np.sqrt(np.maximum(v, 1e-9)), np.inf)
    i = int(np.argmin(z))
    return z[i], pos[0] + k[i]


for inp in ("model_inputs", "model_inputs_v3"):
    v = pd.read_csv(f"{inp}/variants.tsv.gz", sep="\t")
    v = v[(v.include == 1) & v.log_mr.notna()].copy(); v["m"] = np.exp(v.log_mr); v["obs"] = (v.AC_pass > 0).astype(float)
    named = pd.read_csv(f"{inp}/genes.tsv.gz", sep="\t").set_index("label").named
    rows = []
    for lab, g in v.groupby("label"):
        g = g.sort_values(["n", "alt"]); n = g.n.to_numpy(); obs = g.obs.to_numpy(); m = g.m.to_numpy()
        p = 1 - np.exp(-fit_k(obs, m) * m)
        r = {"label": lab[:22], "gene": g.gene.iloc[0], "named": bool(named[lab])}
        for w in (15, 20, 30):
            r[f"z{w}"] = minz(n, obs, p, w)[0]
        z0, st = minz(n, obs, p, 20)
        sh = np.array([minz(n, np.roll(obs, rng.integers(1, len(obs))), p, 20)[0] for _ in range(2000)])
        pm = np.array([minz(n, rng.permutation(obs), p, 20)[0] for _ in range(2000)])
        r["p_shift"] = (np.sum(sh <= z0) + 1) / 2001; r["p_perm"] = (np.sum(pm <= z0) + 1) / 2001
        win = (n >= st) & (n < st + 20); seen = obs > 0
        r["win"] = f"{st}-{st + 19}"
        r["AN_win"] = np.nanmedian(g.AN.to_numpy()[win & seen]) if (win & seen).any() else np.nan
        r["AN_rest"] = np.nanmedian(g.AN.to_numpy()[~win & seen])
        rows.append(r)
    t = pd.DataFrame(rows); t["dis"] = t.gene.map(lambda x: DISEASE_GENES.get(x, ("", ""))[1]); t = t.sort_values("z20")
    print(f"\n== {inp}"); print(t.drop(columns="gene").head(10).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    for col in ("z15", "z20", "z30", "p_shift", "p_perm"):
        P = t[t.gene.isin(NDD)][col].to_numpy(); N = t[~t.gene.isin(DISEASE_GENES) & t.named][col].to_numpy()
        print(f"  NDD vs other named, {col}: AUROC {np.mean((P[:, None] < N[None, :]) + 0.5 * (P[:, None] == N[None, :])):.3f}")
