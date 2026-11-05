#!/usr/bin/env python3
"""Post hoc, written before running (5 Nov 2026, in response to an external review; PREREG s14).

The review asked for Rfam performance in every labelled gene (RNU12 and RNU6ATAC were never shown) and a gene-by-gene
test of the alignment-depth explanation. Per named gene with >= 5 public P/LP SNVs:
  positives: L2_nonblind + L2_curated_additions P/LP (RNU4-2, RNU4ATAC, RNU12, RNU6ATAC) and, for RNU2-2, the 35
  A8-rule literature SNVs (public; the sealed set is not used);
  negatives: the gene's L3 (gnomAD PASS AC >= 9) minus positives.
Scores: Rfam ref column freq, Rfam alt rarity, CADD, phyloP447, GPN-MSA. AUROC; Rfam ref - CADD with a 2,000-draw
class-stratified bootstrap interval; gene-wide rarity-only null excess (2,000 draws).
Alignment-depth prediction stated in advance: across genes, the Rfam-minus-CADD advantage is larger where phyloP447
agrees less with the Rfam profile (results/alignment_depth.tsv), i.e. a negative Spearman. Five genes: descriptive,
with the exact permutation p reported. Reported whatever the result.

    python per_gene_rfam.py   # -> results/per_gene_rfam.tsv
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snrna as S  # noqa: E402
from rarity_null import rarity_null  # noqa: E402

SCORES = {"Rfam ref column freq": "rfam_ref", "Rfam alt rarity": "rfam_alt", "CADD v1.7": "cadd_phred",
          "phyloP447": "phylop447", "GPN-MSA": "gpn_msa"}
B = 2000
SEED = 20261105


def main():
    v = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    comp = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    v = v.merge(comp[["ID", "phylop447", "gpn_msa"]], on="ID", how="left")
    v["rfam_alt"], v["rfam_ref"] = 1 - v.alt_col_freq, v.ref_col_freq
    lab = pd.concat([pd.read_csv(os.path.join(HERE, "labels", f), sep="\t")
                     for f in ("L2_nonblind.tsv", "L2_curated_additions.tsv")])
    lit = pd.read_csv(os.path.join(HERE, "labels", "A8_literature.tsv"), sep="\t")
    plp = set(lab.ID[lab.label == "P/LP"]) | set(lit.ID[(lit.gene == "RNU2-2") & (lit.label == "P/LP")])
    rng = np.random.default_rng(SEED)
    rows = []
    for g in ("RNU2-2", "RNU4-2", "RNU4ATAC", "RNU12", "RNU6ATAC"):
        x = v[v.gene == g].copy()
        for c in SCORES.values():
            x[c] = x[c].astype(float).fillna(x[c].astype(float).median())
        pos = x.ID.isin(plp).to_numpy()
        neg = x.L3_benign_like.astype(bool).to_numpy() & ~pos
        if pos.sum() < 5:
            continue
        m = pos | neg
        y = pos[m]
        i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
        idx = [np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))] for _ in range(B)]
        row = {"gene": g, "n_pos": int(pos.sum()), "n_neg": int(neg.sum())}
        for name, col in SCORES.items():
            s = x[col].to_numpy(float)
            row[name] = S.auroc(y, s[m])
            r = rarity_null(s[pos], s[neg], score_pool=s[~pos & ~neg], draws=B, seed=SEED)
            row[f"{name} excess"], row[f"{name} p"] = r["excess"], r["p"]
        sr, sc = x.rfam_ref.to_numpy(float)[m], x.cadd_phred.to_numpy(float)[m]
        d = np.array([S.auroc(y[j], sr[j]) - S.auroc(y[j], sc[j]) for j in idx])
        row["ref - CADD"], row["ref - CADD lo"], row["ref - CADD hi"] = (row["Rfam ref column freq"] - row["CADD v1.7"],
                                                                          *np.percentile(d, [2.5, 97.5]))
        rows.append(row)
    out = pd.DataFrame(rows)
    ad = pd.read_csv(os.path.join(HERE, "results", "alignment_depth.tsv"), sep="\t")[["gene", "rho_phyloP447_vs_rfam"]]
    out = out.merge(ad, on="gene", how="left")
    out.to_csv(os.path.join(HERE, "results", "per_gene_rfam.tsv"), sep="\t", index=False, float_format="%.4f")
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(out.to_string(index=False, float_format=lambda z: f"{z:.3f}"))
    a, b = out["ref - CADD"].to_numpy(), out.rho_phyloP447_vs_rfam.to_numpy()
    rho = spearmanr(a, b).correlation
    perms = [spearmanr(a, np.array(p)).correlation for p in itertools.permutations(b)]
    p_one = np.mean(np.array(perms) <= rho + 1e-12)
    print(f"\nAcross {len(out)} genes: Spearman(Rfam ref - CADD advantage, phyloP447-Rfam agreement) = {rho:+.2f}; "
          f"exact one-sided p (negative) = {p_one:.3f}")


if __name__ == "__main__":
    main()
