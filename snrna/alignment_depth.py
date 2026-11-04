#!/usr/bin/env python3
"""Exploratory mechanism check (written 5 Nov 2026 before running; PREREG s14): are genome-wide conservation scores
capped at some snRNA loci because whole-genome alignments capture few species there?

Per named snRNA gene (positions in the transcript), from comparators/comparators.tsv.gz:
  - phyloP447 and phyloP100 ceiling: 95th percentile over positions (a column conserved across many aligned species
    scores high; a low ceiling at a deeply conserved RNA means few species are aligned)
  - flat: fraction of positions whose phyloP447 equals the gene's most frequent value (rounded to 3 dp)
  - rho_rfam: per-position Spearman between phyloP447 and Rfam reference-column frequency (agreement of the genome
    alignment with the curated RNA-family alignment)
Prediction stated in advance: low ceilings, flat runs and weak agreement in RNU2-2, RNU4-2 and RNU12 (CADD fails
there), not in RNU4ATAC and RNU6ATAC (CADD works). Descriptive; few genes.

    python alignment_depth.py   # -> results/alignment_depth.tsv
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
CADD_AUROC = {"RNU2-2": 0.55, "RNU4-2": 0.52, "RNU12": 0.45, "RNU4ATAC": 0.77, "RNU6ATAC": 0.71}  # PREREG s14 (all labels)


def main():
    v = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    c = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    genes = pd.read_csv(os.path.join(HERE, "model_inputs", "genes.tsv.gz"), sep="\t")
    named = set(genes.gene[genes.named])
    v = v[v.gene.isin(named)].merge(c[["ID", "phylop447", "phylop100"]], on="ID", how="left")
    pos = v.groupby(["gene", "n"]).agg(p447=("phylop447", "first"), p100=("phylop100", "first"),
                                       rfam=("ref_col_freq", "first")).reset_index()
    rows = []
    for g, x in pos.groupby("gene"):
        p = x.p447.round(3)
        mode = p.mode().iloc[0]
        rows.append({"gene": g, "positions": len(x), "phyloP447_p95": np.nanpercentile(x.p447, 95),
                     "phyloP447_max": x.p447.max(), "phyloP100_p95": np.nanpercentile(x.p100, 95),
                     "flat_frac_447": float((p == mode).mean()), "modal_447": mode,
                     "rho_phyloP447_vs_rfam": spearmanr(x.p447, x.rfam, nan_policy="omit").correlation,
                     "rho_phyloP100_vs_rfam": spearmanr(x.p100, x.rfam, nan_policy="omit").correlation,
                     "cadd_auroc_all_labels": CADD_AUROC.get(g, np.nan)})
    out = pd.DataFrame(rows).sort_values("phyloP447_p95")
    out.to_csv(os.path.join(HERE, "results", "alignment_depth.tsv"), sep="\t", index=False, float_format="%.3f")
    print(out.to_string(index=False, float_format=lambda z: f"{z:.3f}"))
    k = out.dropna(subset=["cadd_auroc_all_labels"])
    for col in ("phyloP447_p95", "rho_phyloP447_vs_rfam", "flat_frac_447"):
        print(f"Spearman across {len(k)} genes, {col} vs CADD AUROC: {spearmanr(k[col], k.cadd_auroc_all_labels).correlation:+.2f}")


if __name__ == "__main__":
    main()
