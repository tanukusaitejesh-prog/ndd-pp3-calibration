#!/usr/bin/env python3
"""Post hoc, written before running (5 Nov 2026, in response to an external review; PREREG s14).

Question: does the Rfam signal on RNU2-2 patient variants survive conditioning on region, structure and local
position? Uses only public, committed labels, so it runs without the sealed set: the 35 RNU2-2 SNVs accepted by the
frozen A8 wording rule (labels/A8_literature.tsv; 16 of them are the A8 test set, the rest overlap the sealed set).
Negatives: RNU2-2 L3 (gnomAD PASS AC >= 9) minus these positives. A few sealed-only positives (not public) may remain
among the negatives; this biases every score towards 0.5 equally.

Scores: Rfam reference-column frequency, Rfam alternative-base rarity, CADD, phyloP447, GPN-MSA, ReMM, and the four
position-only baselines of blind_robustness.py (Table 2 region membership, distance to region, 5' proximity, position
unobserved in gnomAD; the last counts the negative's own allele and is optimistic by construction).
Nulls (10,000 draws each; one-sided p = (k+1)/(B+1)):
  N0 gene-wide: for each positive, a random non-L3, non-positive RNU2-2 SNV (the original rarity-only null).
  N1 region x pairing: drawn from SNVs with the positive's own Table 2 region status and reference pairing status
     (paired vs unpaired in the reference structure).
  N2 local window: drawn from SNVs within +-5 nt of the positive.
Also: rarity-matched AUROC (positives vs the N0 pool) and AUROC inside and outside n.1-67.
Prediction stated in advance: excess shrinks from N0 to N2 for every score; the question is whether the Rfam excess
stays above 0 under N1 and N2 while the position baselines' excess goes to ~0 under N2. Reported whatever the result.

    python rnu2_position_control.py   # -> results/rnu2_position_control.tsv
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snrna as S  # noqa: E402
from rarity_null import rarity_null  # noqa: E402

GENE = "RNU2-2"
DRAWS = 10000
SEED = 20261105
SCORES = {"Rfam ref column freq": "rfam_ref", "Rfam alt rarity": "rfam_alt", "CADD v1.7": "cadd_phred",
          "phyloP447": "phylop447", "GPN-MSA": "gpn_msa", "ReMM v0.4": "remm",
          "BASELINE in region n.1-67": "in_region", "BASELINE distance to region": "minus_dist",
          "BASELINE 5' proximity": "five_prime", "BASELINE position unobserved": "pos_unobserved"}


def load():
    var = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    v = var[var.gene == GENE].copy()
    comp = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    v = v.merge(comp[["ID", "phylop447", "gpn_msa", "remm"]], on="ID", how="left")
    reg = pd.read_csv(os.path.join(HERE, "table2_regions.tsv"), sep="\t", comment="#")
    reg = reg[reg.gene == GENE]
    inreg = np.zeros(len(v), bool)
    dist = np.full(len(v), np.inf)
    for r in reg.itertuples():
        inreg |= v.n.between(r.start, r.end).to_numpy()
        dist = np.minimum(dist, np.where(v.n < r.start, r.start - v.n, np.where(v.n > r.end, v.n - r.end, 0)))
    v["in_region"] = inreg.astype(float)
    v["minus_dist"] = -dist
    v["five_prime"] = -v.n.astype(float)
    obs = v[v.AC_pass > 0].groupby("n").size()
    v["pos_unobserved"] = -v.n.map(obs).fillna(0)
    v["rfam_alt"] = 1 - v.alt_col_freq
    v["rfam_ref"] = v.ref_col_freq
    v["paired"] = (v.ref_intra_unpaired.fillna(1) == 0).astype(int)
    for c in SCORES.values():
        v[c] = v[c].astype(float).fillna(v[c].astype(float).median())
    lab = pd.read_csv(os.path.join(HERE, "labels", "A8_literature.tsv"), sep="\t")
    plp = set(lab.ID[(lab.gene == GENE) & (lab.label == "P/LP")])
    v["pos"] = v.ID.isin(plp)
    v["neg"] = (v.L3_benign_like.astype(bool)) & ~v.pos
    v["pool"] = ~v.pos & ~v.neg
    return v


def main():
    v = load()
    P, N, U = v[v.pos], v[v.neg], v[v.pool]
    print(f"{GENE}: {len(P)} positives, {len(N)} L3 negatives, {len(U)} pool SNVs; "
          f"{int(P.in_region.sum())} positives inside n.1-67")
    rows = []
    for name, col in SCORES.items():
        sp, sn = P[col].to_numpy(), N[col].to_numpy()
        g_n1 = [U[(U.in_region == r.in_region) & (U.paired == r.paired)][col].to_numpy() for r in P.itertuples()]
        g_n2 = [U[(U.n - r.n).abs() <= 5][col].to_numpy() for r in P.itertuples()]
        row = {"score": name, "auroc": S.auroc(np.r_[np.ones(len(sp)), np.zeros(len(sn))], np.r_[sp, sn]),
               "rarity_matched_auroc": S.auroc(np.r_[np.ones(len(sp)), np.zeros(len(U))], np.r_[sp, U[col].to_numpy()])}
        for tag, kw in (("N0", {"score_pool": U[col].to_numpy()}), ("N1", {"groups": g_n1}), ("N2", {"groups": g_n2})):
            r = rarity_null(sp, sn, draws=DRAWS, seed=SEED, **kw)
            row[f"{tag}_null"], row[f"{tag}_excess"], row[f"{tag}_p"] = r["null_mean"], r["excess"], r["p"]
        for tag, m in (("inside", v.in_region == 1), ("outside", v.in_region == 0)):
            x = v[m & (v.pos | v.neg)]
            row[f"auroc_{tag}"] = S.auroc(x.pos.to_numpy(), x[col].to_numpy()) if 0 < x.pos.sum() < len(x) else np.nan
            row[f"n_{tag}"] = f"{int(x.pos.sum())} vs {int((~x.pos).sum())}"
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(HERE, "results", "rnu2_position_control.tsv"), sep="\t", index=False, float_format="%.4f")
    with pd.option_context("display.width", 250):
        print(out.to_string(index=False, float_format=lambda z: f"{z:.3f}"))


if __name__ == "__main__":
    main()
