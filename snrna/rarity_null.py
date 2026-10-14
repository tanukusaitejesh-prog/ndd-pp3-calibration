#!/usr/bin/env python3
"""Rarity-only null for variant-effect benchmarks that use population-common variants as negatives.

When negatives are common population variants (e.g. gnomAD AC >= 9), any score that tracks selection or conservation
separates a RANDOM rare variant from them, so a score's AUROC has a score-specific null above 0.5. Report the AUROC,
this null, and the excess (AUROC - null); the excess is the part not explained by rarity alone.

    from rarity_null import rarity_null
    res = rarity_null(score_pos, score_neg, score_pool, draws=1000, seed=0)
    # res = {"auroc": .., "null_mean": .., "null_lo": .., "null_hi": .., "excess": .., "p": ..}

score_pos   scores of the positives (e.g. ClinVar P/LP)
score_neg   scores of the common negatives
score_pool  scores of candidate pseudo-positives: variants from the same loci (or windows) that are neither positives
            nor negatives, i.e. rare or absent in the population. Matching the pool to the positives' loci matters:
            draw from the same genes, or pass `groups` to draw within each positive's own window/gene.
groups      optional: list with one array of pool scores per positive (window-matched null, as in A9)

Command line (TSV with columns score and class in {pos, neg, pool}):
    python rarity_null.py scores.tsv --score cadd_phred
"""
import argparse

import numpy as np


def _auroc(y, s):
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s))
    sorted_s = s[order]
    i = 0
    while i < len(s):  # average ranks for ties
        j = i
        while j + 1 < len(s) and sorted_s[j + 1] == sorted_s[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    n1, n0 = y.sum(), (~y).sum()
    if n1 == 0 or n0 == 0:
        return float("nan")
    return float((ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def rarity_null(score_pos, score_neg, score_pool=None, groups=None, draws=1000, seed=0):
    rng = np.random.default_rng(seed)
    pos, neg = np.asarray(score_pos, float), np.asarray(score_neg, float)
    y = np.r_[np.ones(len(pos), bool), np.zeros(len(neg), bool)]
    au = _auroc(y, np.r_[pos, neg])
    if groups is not None:
        draw = lambda: np.array([rng.choice(np.asarray(g, float)) for g in groups])
    else:
        pool = np.asarray(score_pool, float)
        if len(pool) < len(pos):
            raise ValueError("pool must hold at least as many variants as there are positives")
        draw = lambda: rng.choice(pool, len(pos), replace=False)
    nulls = np.array([_auroc(y, np.r_[draw(), neg]) for _ in range(draws)])
    lo, hi = np.percentile(nulls, [2.5, 97.5])
    return {"auroc": au, "null_mean": float(nulls.mean()), "null_lo": float(lo), "null_hi": float(hi),
            "excess": au - float(nulls.mean()), "p": float((np.sum(nulls >= au) + 1) / (draws + 1))}


def main():
    import pandas as pd
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv")
    ap.add_argument("--score", required=True)
    ap.add_argument("--draws", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    d = pd.read_csv(a.tsv, sep="\t").dropna(subset=[a.score])
    r = rarity_null(d[a.score][d["class"] == "pos"], d[a.score][d["class"] == "neg"], d[a.score][d["class"] == "pool"],
                    draws=a.draws, seed=a.seed)
    print(f"{a.score}: AUROC {r['auroc']:.3f}  rarity-only null {r['null_mean']:.3f} [{r['null_lo']:.3f}, {r['null_hi']:.3f}]"
          f"  excess {r['excess']:+.3f}  p {r['p']:.4f}")


if __name__ == "__main__":
    main()
