#!/usr/bin/env python3
"""Post hoc, written before running (5 Nov 2026, in response to an external review; PREREG s14).

The review noted that Rfam reference-column frequency is constant across the three alternatives at a position, so the
variant-level bootstrap in sge_classes.py treats ~435 SNVs as independent when ~145 positions carry the signal.
Two checks on RNU4-2 SGE (no population labels):
 (1) Position-cluster bootstrap: resample positions with replacement (all SNVs at a drawn position come along),
     2,000 draws; AUROC for 'depleted' (< -0.302) and paired differences Rfam ref col freq - each genome-wide score,
     overall and outside the critical region.
 (2) Within-position test (conditions on position entirely): among the alternatives at the same position, does the
     score rank SGE damage? Pooled within-position Spearman of each variant-varying score (Rfam alt rarity, CADD,
     GPN-MSA, ReMM, ViennaRNA ddG) after centring score and -function score within position; position-cluster
     bootstrap interval. Reference-column frequency and phyloP are constant within a position, so they are undefined
     here by construction.
Expectation stated in advance: (1) intervals widen; whether Rfam - genome-wide still excludes 0 is the question.
(2) unknown; reported whatever the result.

    python sge_position_control.py   # -> results/sge_position_control.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snrna as S  # noqa: E402

B = 2000
SEED = 20261105
GW = {"CADD v1.7": "cadd_phred", "GPN-MSA": "gpn_msa", "phyloP447": "phylop447", "ReMM v0.4": "remm"}
VARYING = {"Rfam alt rarity": "alt_rarity", "CADD v1.7": "cadd_phred", "GPN-MSA": "gpn_msa", "ReMM v0.4": "remm",
           "ViennaRNA ddG_bind": "ddG_bind_max"}


def load():
    v = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    comp = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    v = v.merge(comp[["ID", "phylop447", "remm", "gpn_msa"]], on="ID", how="left")
    v["alt_rarity"] = 1 - v.alt_col_freq
    s = pd.read_csv(os.path.join(HERE, "work", "sge.tsv"), sep="\t")
    s = s[s.Type == "SNV"][["ID", "function_score_mean", "critical"]].dropna(subset=["function_score_mean"])
    t = v[v.gene == "RNU4-2"].merge(s, on="ID")
    for c in set(GW.values()) | set(VARYING.values()) | {"ref_col_freq"}:
        t[c] = t[c].fillna(t[c].median())
    return t


def cluster_draws(pos, rng):
    up = np.unique(pos)
    groups = {p: np.flatnonzero(pos == p) for p in up}
    return [np.concatenate([groups[p] for p in rng.choice(up, len(up))]) for _ in range(B)]


def main():
    t = load()
    rng = np.random.default_rng(SEED)
    rows = []
    for subset, m in (("all", np.ones(len(t), bool)), ("outside critical region", ~t.critical.astype(bool).to_numpy())):
        x = t[m].reset_index(drop=True)
        y = (x.function_score_mean < -0.302).to_numpy()
        draws = [d for d in cluster_draws(x.n.to_numpy(), rng)]
        ok = [d for d in draws if 0 < y[d].sum() < len(d)]
        ref = x.ref_col_freq.to_numpy(float)
        bref = np.array([S.auroc(y[d], ref[d]) for d in ok])
        rows.append({"test": "cluster AUROC", "subset": subset, "score": "Rfam ref column freq", "n_pos_sites": x.n.nunique(),
                     "value": S.auroc(y, ref), "lo": np.percentile(bref, 2.5), "hi": np.percentile(bref, 97.5)})
        for name, col in {**GW, "Rfam alt rarity": "alt_rarity"}.items():
            sc = x[col].to_numpy(float)
            bs = np.array([S.auroc(y[d], sc[d]) for d in ok])
            rows.append({"test": "cluster AUROC", "subset": subset, "score": name, "n_pos_sites": x.n.nunique(),
                         "value": S.auroc(y, sc), "lo": np.percentile(bs, 2.5), "hi": np.percentile(bs, 97.5)})
            if name in GW:
                d = bref - bs
                rows.append({"test": "cluster diff (ref - score)", "subset": subset, "score": name,
                             "n_pos_sites": x.n.nunique(), "value": S.auroc(y, ref) - S.auroc(y, sc),
                             "lo": np.percentile(d, 2.5), "hi": np.percentile(d, 97.5)})
    # within-position
    x = t.copy()
    x["dmg"] = -x.function_score_mean
    x = x[x.groupby("n").ID.transform("size") >= 2]
    dm = (x.dmg - x.groupby("n").dmg.transform("mean")).to_numpy()
    draws = cluster_draws(x.n.to_numpy(), rng)
    for name, col in VARYING.items():
        sc = (x[col] - x.groupby("n")[col].transform("mean")).to_numpy(float)
        rho = spearmanr(sc, dm).correlation
        bs = np.array([spearmanr(sc[d], dm[d]).correlation for d in draws])
        rows.append({"test": "within-position Spearman", "subset": "all", "score": name, "n_pos_sites": x.n.nunique(),
                     "value": rho, "lo": np.nanpercentile(bs, 2.5), "hi": np.nanpercentile(bs, 97.5)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(HERE, "results", "sge_position_control.tsv"), sep="\t", index=False, float_format="%.3f")
    print(out.to_string(index=False, float_format=lambda z: f"{z:.3f}"))


if __name__ == "__main__":
    main()
