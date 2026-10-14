#!/usr/bin/env python3
"""v4 (amendment A5): stacked, regularised per-variant classifier, frozen before the single unsealing.

Designed after the v2 endpoints were seen (they showed: position-level information and lab-data transfer carry the
signal, the graph network overfitted). Trained and tuned on non-blind genes only; RNU2-2 and RNU5B-1 (sealed) are
excluded from training entirely and scored once, blind, by v4_blind_eval.py inside final_eval.py.

Features (each converted to a within-gene percentile rank, missing -> 0.5, so gene-level shifts cancel and the model
transfers to genes it has not seen):
  cadd_phred, phylop447, remm, gpn_msa, evo2 (-(alt - ref) log-likelihood, mean of strands), ddG_bind_max, |ddG_fold|,
  1 - alt_col_freq (Rfam profile), ref_col_freq, ref_paired, alt_inter_broken, local constraint (-z of the 21-nt
  window centred on the site, k fitted within the locus; gnomAD + UK Biobank counts when model_inputs_v3 exists),
  lab function (homology transfer h, RNU4ATAC only; a gene's own SGE is not used: RNU4-2 labels may rest on it)
Labels: L2 P/LP SNVs (positives) vs L3 SNVs not in L2 (negatives), named non-sealed genes.
Model (fixed in advance, no tuning): L2-penalised logistic regression, C = 0.1, balanced class weights.
Internal validation: leave-one-paralog-cluster-out (genes.tsv cluster) for every gene with >= 5 L2 SNVs; AUROC of v4
vs CADD with 2000 bootstrap CI; weighted delta AUROC over units with >= 10 L2 SNVs (the G3 rule).

    python v4.py [--out results_v4]      # -> v4_cv.tsv, v4_report.txt, v4_model.json, v4_scores.tsv.gz, v4_manifest.txt
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evaluate as E  # noqa: E402
import snrna as S  # noqa: E402
from gene_discovery import fit_k  # noqa: E402

C_REG = 0.1
SEED = 20261014
FEATURES = ["cadd_phred", "phylop447", "remm", "gpn_msa", "evo2", "ddG_bind_max", "abs_ddG_fold", "alt_rarity",
            "ref_col_freq", "ref_paired", "alt_inter_broken", "local_constraint", "lab_function"]


def local_constraint(v, w=21):
    """-z of observed vs expected SNVs in the w-nt window centred on each site (k refitted within the locus)."""
    out = pd.Series(np.nan, index=v.index)
    ok = (v.include == 1) & v.log_mr.notna()
    for lab, g in v[ok].groupby("label"):
        m = np.exp(g.log_mr.to_numpy()); obs = (g.AC_pass > 0).to_numpy(float)
        p = 1 - np.exp(-fit_k(obs, m) * m)
        n = g.n.to_numpy(); n0 = n.min(); L = n.max() - n0 + 1
        O = np.bincount(n - n0, obs, L); Ex = np.bincount(n - n0, p, L); V = np.bincount(n - n0, p * (1 - p), L)
        h = w // 2
        cs = lambda x: np.concatenate([[0], np.cumsum(x)])
        cO, cE, cV = cs(O), cs(Ex), cs(V)
        lo = np.clip(np.arange(L) - h, 0, L); hi = np.clip(np.arange(L) + h + 1, 0, L)
        z = ((cO[hi] - cO[lo]) - (cE[hi] - cE[lo])) / np.sqrt(np.maximum(cV[hi] - cV[lo], 1e-9))
        allg = v.label == lab
        out[allg] = -z[(v.n[allg] - n0).clip(0, L - 1).to_numpy()]
    return out


def build_features(inputs, inputs_counts, evo2_dir):
    v = pd.read_csv(os.path.join(inputs, "variants.tsv.gz"), sep="\t")
    vc = pd.read_csv(os.path.join(inputs_counts, "variants.tsv.gz"), sep="\t") if inputs_counts != inputs else v
    v["local_constraint"] = local_constraint(v.assign(AC_pass=v.ID.map(dict(zip(vc.ID, vc.AC_pass)))))
    comp = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    v = v.merge(comp[["ID", "phylop447", "remm", "gpn_msa"]], on="ID", how="left")
    dp = os.path.join(evo2_dir, "deltas.tsv")
    if os.path.exists(dp):
        d = pd.read_csv(dp, sep="\t").drop_duplicates("ID")
        delta = ((d.ll_alt_fwd - d.ll_ref_fwd) + (d.ll_alt_rc - d.ll_ref_rc)) / 2
        v["evo2"] = v.ID.map(dict(zip(d.ID, -delta)))
    else:
        v["evo2"] = np.nan
    v["abs_ddG_fold"] = v.ddG_fold.abs()
    v["alt_rarity"] = 1 - v.alt_col_freq
    v["ref_paired"] = 1 - v.ref_intra_unpaired
    v["lab_function"] = v.ID.map(E.homology_scores(inputs))  # RNU4ATAC only (RNU4-2 SGE transferred); no own-gene SGE
    for f in FEATURES:  # within-gene percentile rank, missing -> 0.5
        r = v.groupby("label")[f].rank(pct=True)
        v[f + "_r"] = r.fillna(0.5)
    return v


def auc_ci(y, x, B, rng):
    a = S.auroc(y, x)
    i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
    bs = [S.auroc(np.r_[np.ones(len(i1), bool), np.zeros(len(i0), bool)],
                  x[np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]]) for _ in range(B)]
    return a, *np.nanpercentile(bs, [2.5, 97.5])


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    from sklearn.linear_model import LogisticRegression
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--inputs-counts", default=None, help="counts for local constraint (default model_inputs_v3 if present)")
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--labels", default=os.path.join(HERE, "labels"))
    ap.add_argument("--out", default=os.path.join(HERE, "results_v4"))
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()
    ic = a.inputs_counts or (os.path.join(HERE, "model_inputs_v3") if os.path.exists(os.path.join(HERE, "model_inputs_v3", "variants.tsv.gz")) else a.inputs)
    rng = np.random.default_rng(SEED)
    v = build_features(a.inputs, ic, a.evo2)
    genes = pd.read_csv(os.path.join(a.inputs, "genes.tsv.gz"), sep="\t")
    cl = genes.set_index("label").cluster
    named = set(genes.label[genes.named & ~genes.gene.isin(S.SEALED_GENES)])
    L = E.load_l2(a.labels)
    plp = set(L.ID[L.label == "P/LP"])
    v["y"] = v.ID.isin(plp)
    train = v[v.label.isin(named) & (v.y | (v.L3_benign_like.astype(bool) & ~v.y))].copy()
    X = [f + "_r" for f in FEATURES]
    mk = lambda: LogisticRegression(C=C_REG, class_weight="balanced", max_iter=5000)
    rows, sets = [], []
    for lab, g in train.groupby("label"):
        if g.y.sum() < 5 or (~g.y).sum() < 5:
            continue
        tr = train[train.label.map(cl) != cl[lab]]
        m = mk().fit(tr[X], tr.y)
        s4 = m.decision_function(g[X]); y = g.y.to_numpy(); cadd = g.cadd_phred.to_numpy(float)
        a4, lo4, hi4 = auc_ci(y, s4, a.boot, rng); ac, loc, hic = auc_ci(y, cadd, a.boot, rng)
        rows.append({"gene": g.gene.iloc[0], "cluster_held_out": cl[lab], "n_L2": int(y.sum()), "n_L3": int((~y).sum()),
                     "auroc_v4": a4, "v4_lo": lo4, "v4_hi": hi4, "auroc_cadd": ac, "cadd_lo": loc, "cadd_hi": hic,
                     "delta": a4 - ac})
        if y.sum() >= 10:
            sets.append((g.gene.iloc[0], y, s4, cadd))
    cv = pd.DataFrame(rows)
    rep = ["v4 (A5) internal validation: leave-one-paralog-cluster-out on non-blind genes", ""]
    rep += [cv.to_string(index=False, float_format=lambda x: f"{x:.3f}"), ""]
    if sets:
        w = np.array([y.sum() for _, y, _, _ in sets], float)
        pt = float((w * np.array([S.auroc(y, s) - S.auroc(y, c) for _, y, s, c in sets])).sum() / w.sum())
        bs = []
        for _ in range(a.boot):
            vals = []
            for _, y, s, c in sets:
                i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
                idx = np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]
                vals.append(S.auroc(y[idx], s[idx]) - S.auroc(y[idx], c[idx]))
            bs.append((w * np.array(vals)).sum() / w.sum())
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        rep.append(f"Weighted delta AUROC v4 - CADD over units with >= 10 L2 ({', '.join(u for u, *_ in sets)}): "
                   f"{pt:+.3f} [{lo:+.3f}, {hi:+.3f}]")
    final = mk().fit(train[X], train.y)
    coef = dict(zip(FEATURES, final.coef_[0].round(4).tolist()))
    rep += ["", "Final model (all non-blind training genes), coefficients on within-gene percentile ranks:",
            *[f"  {k:18s} {c:+.3f}" for k, c in sorted(coef.items(), key=lambda t: -abs(t[1]))],
            f"  intercept {final.intercept_[0]:+.3f}", f"Training: {int(train.y.sum())} L2 vs {int((~train.y).sum())} L3 SNVs; "
            f"counts for local constraint from {os.path.basename(ic)}; evo2 deltas {'yes' if v.evo2.notna().any() else 'NO'}"]
    os.makedirs(a.out, exist_ok=True)
    v["v4"] = final.decision_function(v[X])
    v[["ID", "label", "gene", "v4"] + FEATURES].to_csv(os.path.join(a.out, "v4_scores.tsv.gz"), sep="\t", index=False,
                                                      float_format="%.5g", compression="gzip")
    cv.to_csv(os.path.join(a.out, "v4_cv.tsv"), sep="\t", index=False, float_format="%.4f")
    json.dump({"features": FEATURES, "C": C_REG, "coef": coef, "intercept": float(final.intercept_[0]),
               "counts": ic, "seed": SEED}, open(os.path.join(a.out, "v4_model.json"), "w"), indent=1)
    open(os.path.join(a.out, "v4_report.txt"), "w").write("\n".join(rep) + "\n")
    files = [os.path.join(HERE, f) for f in ("v4.py", "v4_blind_eval.py")] + \
            [os.path.join(a.out, f) for f in ("v4_scores.tsv.gz", "v4_model.json")]
    with open(os.path.join(a.out, "v4_manifest.txt"), "w") as f:
        for p in files:
            if os.path.exists(p):
                f.write(f"{sha(p)}  {os.path.relpath(p, HERE)}\n")
    print("\n".join(rep))
    print("\nFREEZE: paste results_v4/v4_manifest.txt into the chat (amendment A5) BEFORE running final_eval.py")


if __name__ == "__main__":
    main()
