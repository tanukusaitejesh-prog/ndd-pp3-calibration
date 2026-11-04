#!/usr/bin/env python3
"""Reproduction (4 Nov 2026; PREREG s14): Spearman correlation of each score with RNU4-2 SGE functional damage
(-function_score_mean) over the in-transcript RNU4-2 SNVs with an SGE score, with 1,000 bootstrap intervals and a fixed
seed, so that the manuscript and the supplement cite one saved file. Same inputs and merging as sge_classes.py.

    python sge_spearman.py      # -> results/sge_spearman.tsv
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
SCORES = {"Rfam ref column freq": "ref_col_freq", "Rfam alt rarity": "alt_rarity", "CADD v1.7": "cadd_phred",
          "GPN-MSA": "gpn_msa", "phyloP447": "phylop447", "ReMM v0.4": "remm", "ViennaRNA ddG_bind": "ddG_bind_max"}
B = 1000
SEED = 20261104


def main():
    v = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    comp = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    v = v.merge(comp[["ID", "phylop447", "remm", "gpn_msa"]], on="ID", how="left")
    v["alt_rarity"] = 1 - v.alt_col_freq
    s = pd.read_csv(os.path.join(HERE, "work", "sge.tsv"), sep="\t")
    s = s[s.Type == "SNV"][["ID", "function_score_mean"]].dropna()
    t = v[v.gene == "RNU4-2"].merge(s, on="ID")
    dmg = -t.function_score_mean.to_numpy(float)
    rng = np.random.default_rng(SEED)
    idx = [rng.integers(0, len(t), len(t)) for _ in range(B)]
    rows = []
    for name, col in SCORES.items():
        x = t[col].to_numpy(float)
        rho = spearmanr(x, dmg, nan_policy="omit").correlation
        bs = np.array([spearmanr(x[j], dmg[j], nan_policy="omit").correlation for j in idx])
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        rows.append({"score": name, "n": int(np.isfinite(x).sum()), "rho": rho, "lo": lo, "hi": hi})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(HERE, "results", "sge_spearman.tsv"), sep="\t", index=False, float_format="%.3f")
    print(f"{len(t)} RNU4-2 SNVs with an SGE score")
    print(out.to_string(index=False, float_format=lambda z: f"{z:.3f}"))


if __name__ == "__main__":
    main()
