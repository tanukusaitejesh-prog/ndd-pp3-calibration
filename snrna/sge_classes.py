#!/usr/bin/env python3
"""Exploratory, post hoc (written 5 Nov 2026 BEFORE running; PREREG s14): population-free classification test on
RNU4-2 saturation genome editing (De Jonghe et al. 2026). No population-derived label or negative is involved, so the
rarity contrast that affects every L2-vs-L3 AUROC cannot arise here.

Labels: the 435 in-transcript RNU4-2 SNVs with an SGE function score. Primary class: "depleted" (function score
< -0.302, the study's threshold; 151 of 539 variants overall). Secondary: "strong" (< -0.9).
Scores (higher = more damaging): Rfam reference-column frequency, Rfam alternative-base rarity (parameter-free, no
labels), CADD v1.7, GPN-MSA, phyloP447, ReMM v0.4, ViennaRNA intermolecular ddG (ddG_bind_max).
Metrics: AUROC with 2,000 class-stratified bootstrap intervals; paired bootstrap difference of each Rfam feature minus
each genome-wide score (same resamples). Also within the guidance-defined critical region and outside it, because
position alone separates much of the signal in RNU4-2.
Expectation stated in advance from the known Spearman correlations (Rfam 0.39-0.46 vs genome-wide <= 0.14): Rfam
features exceed every genome-wide score; reported whatever the result.

    python sge_classes.py      # -> results/sge_classes.tsv, results/sge_classes.txt
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snrna as S  # noqa: E402

SCORES = {"Rfam ref column freq": "ref_col_freq", "Rfam alt rarity": "alt_rarity", "CADD v1.7": "cadd_phred",
          "GPN-MSA": "gpn_msa", "phyloP447": "phylop447", "ReMM v0.4": "remm", "ViennaRNA ddG_bind": "ddG_bind_max"}
B = 2000
SEED = 20261105


def main():
    v = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    comp = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    v = v.merge(comp[["ID", "phylop447", "remm", "gpn_msa"]], on="ID", how="left")
    v["alt_rarity"] = 1 - v.alt_col_freq
    s = pd.read_csv(os.path.join(HERE, "work", "sge.tsv"), sep="\t")
    s = s[(s.Type == "SNV")][["ID", "function_score_mean", "critical"]].dropna(subset=["function_score_mean"])
    t = v[v.gene == "RNU4-2"].merge(s, on="ID")
    rows, lines = [], []
    rng = np.random.default_rng(SEED)
    for cls, thr in (("depleted (< -0.302)", -0.302), ("strong (< -0.9)", -0.9)):
        for subset, m in (("all", np.ones(len(t), bool)), ("critical region", t.critical.astype(bool).to_numpy()),
                          ("outside critical region", ~t.critical.astype(bool).to_numpy())):
            x = t[m]
            y = (x.function_score_mean < thr).to_numpy()
            if y.sum() < 5 or (~y).sum() < 5:
                continue
            i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
            idx = [np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))] for _ in range(B)]
            boots = {}
            for name, col in SCORES.items():
                sc = x[col].to_numpy(float)
                sc = np.where(np.isnan(sc), np.nanmedian(sc), sc)
                a = S.auroc(y, sc)
                bs = np.array([S.auroc(y[j], sc[j]) for j in idx])
                boots[name] = bs
                lo, hi = np.percentile(bs, [2.5, 97.5])
                rows.append({"class": cls, "subset": subset, "n_pos": int(y.sum()), "n_neg": int((~y).sum()),
                             "score": name, "auroc": a, "lo": lo, "hi": hi})
            for rf in ("Rfam ref column freq", "Rfam alt rarity"):
                for gw in ("CADD v1.7", "GPN-MSA", "phyloP447", "ReMM v0.4"):
                    d = boots[rf] - boots[gw]
                    a_rf = [r["auroc"] for r in rows if r["class"] == cls and r["subset"] == subset and r["score"] == rf][0]
                    a_gw = [r["auroc"] for r in rows if r["class"] == cls and r["subset"] == subset and r["score"] == gw][0]
                    lo, hi = np.percentile(d, [2.5, 97.5])
                    rows.append({"class": cls, "subset": subset, "n_pos": int(y.sum()), "n_neg": int((~y).sum()),
                                 "score": f"{rf} - {gw}", "auroc": a_rf - a_gw, "lo": lo, "hi": hi})
    out = pd.DataFrame(rows)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    out.to_csv(os.path.join(HERE, "results", "sge_classes.tsv"), sep="\t", index=False, float_format="%.4f")
    for (c, sub), g in out.groupby(["class", "subset"], sort=False):
        lines.append(f"{c} | {sub} ({int(g.n_pos.iloc[0])} vs {int(g.n_neg.iloc[0])})")
        for r in g.itertuples():
            lines.append(f"  {r.score:38s} {r.auroc:+.3f} [{r.lo:+.3f}, {r.hi:+.3f}]" if " - " in r.score
                         else f"  {r.score:38s} {r.auroc:.3f} [{r.lo:.3f}, {r.hi:.3f}]")
    open(os.path.join(HERE, "results", "sge_classes.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
