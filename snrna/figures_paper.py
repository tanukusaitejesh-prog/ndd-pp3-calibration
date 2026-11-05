#!/usr/bin/env python3
"""Manuscript figures from committed public inputs and result files (no sealed data).

    python figures_paper.py   # -> paper/figures/fig{1..4}.{png,pdf}
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rnu2_position_control as R  # noqa: E402

OUT = os.path.join(HERE, "paper", "figures")
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "legend.frameon": False, "font.family": "DejaVu Sans"})


def roc(y, s):
    o = np.argsort(-s, kind="mergesort"); y = y[o]
    tp = np.r_[0, np.cumsum(y)] / y.sum(); fp = np.r_[0, np.cumsum(~y)] / (~y).sum()
    return fp, tp


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig1():
    v = R.load(); m = v.pos | v.neg; x = v[m]; y = x.pos.to_numpy()
    tab = pd.read_csv(os.path.join(HERE, "results", "rnu2_position_control.tsv"), sep="\t").set_index("score")
    fig, ax = plt.subplots(figsize=(3.4, 3.2))
    ax.plot([0, 1], [0, 1], color=GRID, lw=1, zorder=0)
    for name, col, c in (("Rfam reference-column frequency", "rfam_ref", BLUE), ("GPN-MSA", "gpn_msa", AQUA),
                         ("phyloP447", "phylop447", YELLOW), ("CADD v1.7", "cadd_phred", ORANGE)):
        fp, tp = roc(y, x[col].to_numpy(float))
        key = {"rfam_ref": "Rfam ref column freq", "gpn_msa": "GPN-MSA", "phylop447": "phyloP447", "cadd_phred": "CADD v1.7"}[col]
        ax.plot(fp, tp, color=c, lw=2, label=f"{name}  {tab.loc[key, 'auroc']:.2f}")
    ax.set_xlabel("False-positive rate"); ax.set_ylabel("True-positive rate")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.01); ax.set_aspect("equal")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=7, title="AUROC", title_fontsize=7)
    ax.set_title("RNU2-2: 35 public patient variants vs 125 population variants", fontsize=8, color=INK, loc="left")
    save(fig, "fig1_rnu2_roc")


def fig2():
    t = pd.read_csv(os.path.join(HERE, "results", "rnu2_position_control.tsv"), sep="\t")
    keep = ["Rfam ref column freq", "Rfam alt rarity", "GPN-MSA", "phyloP447", "CADD v1.7", "ReMM v0.4",
            "BASELINE in region n.1-67", "BASELINE 5' proximity"]
    lab = {"Rfam ref column freq": "Rfam reference-column freq.", "Rfam alt rarity": "Rfam alternative-base rarity",
           "BASELINE in region n.1-67": "Baseline: in region n.1–67", "BASELINE 5' proximity": "Baseline: 5′ proximity"}
    t = t.set_index("score").loc[keep]
    fig, ax = plt.subplots(figsize=(5.2, 3.3))
    ys = np.arange(len(keep))[::-1]
    for k, (tag, c, nm) in enumerate((("N0", BLUE, "Gene-wide null"), ("N1", ORANGE, "Region × pairing null"),
                                      ("N2", AQUA, "±5-nt window null"))):
        ax.scatter(t[f"{tag}_excess"], ys + (1 - k) * 0.22, s=28, color=c, edgecolor=SURF, linewidth=1.2, zorder=3, label=nm)
        for yy, (e, p) in zip(ys + (1 - k) * 0.22, zip(t[f"{tag}_excess"], t[f"{tag}_p"])):
            if p < 0.05:
                ax.text(e + 0.008, yy, "*", va="center", fontsize=8, color=INK)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_yticks(ys); ax.set_yticklabels([lab.get(s, s) for s in keep])
    ax.set_xlabel("AUROC excess over rarity-only null (* one-sided p < 0.05)")
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=3, fontsize=7)
    save(fig, "fig2_rnu2_matched_nulls")


def fig3():
    import sge_classes  # noqa: F401  (same inputs)
    v = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    comp = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    v = v.merge(comp[["ID", "phylop447"]], on="ID", how="left")
    s = pd.read_csv(os.path.join(HERE, "work", "sge.tsv"), sep="\t")
    s = s[s.Type == "SNV"][["ID", "function_score_mean"]].dropna()
    t = v[v.gene == "RNU4-2"].merge(s, on="ID")
    pos = t.groupby("n").agg(rfam=("ref_col_freq", "first"), p447=("phylop447", "first"),
                             fmin=("function_score_mean", "min")).reset_index()
    pc = pd.read_csv(os.path.join(HERE, "results", "sge_position_control.tsv"), sep="\t")
    pc = pc[pc.test == "cluster AUROC"]
    fig = plt.figure(figsize=(7.0, 4.2))
    gs = fig.add_gridspec(3, 2, width_ratios=[1, 1.6], hspace=0.35, wspace=0.35)
    ax = fig.add_subplot(gs[:, 0])
    order = ["Rfam ref column freq", "Rfam alt rarity", "GPN-MSA", "CADD v1.7", "ReMM v0.4", "phyloP447"]
    ys = np.arange(len(order))[::-1]
    for k, (sub, c, nm) in enumerate((("all", BLUE, "All variants"), ("outside critical region", ORANGE, "Outside critical region"))):
        d = pc[pc.subset == sub].set_index("score").loc[order]
        yy = ys + (0.5 - k) * 0.3
        ax.errorbar(d.value, yy, xerr=[d.value - d.lo, d.hi - d.value], fmt="o", ms=5, color=c, ecolor=c,
                    elinewidth=1.4, capsize=0, mec=SURF, mew=1, label=nm)
    ax.axvline(0.5, color=INK2, lw=0.8)
    ax.set_yticks(ys); ax.set_yticklabels(order); ax.set_xlim(0.3, 1.0)
    ax.set_xlabel("AUROC, function-depleted (95% CI, positions resampled)")
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=7)
    ax.text(-0.02, 1.09, "A", transform=ax.transAxes, fontweight="bold", fontsize=10, ha="right")
    for i, (col, lab, c) in enumerate((("rfam", "Rfam ref. column freq.", BLUE), ("p447", "phyloP447", YELLOW),
                                       ("fmin", "SGE function\n(min per position)", INK2))):
        a = fig.add_subplot(gs[i, 1])
        a.bar(pos.n, pos[col], width=0.8, color=c, linewidth=0)
        a.set_ylabel(lab, fontsize=7); a.set_xlim(0.5, pos.n.max() + 0.5)
        if i < 2:
            a.set_xticklabels([])
        else:
            a.set_xlabel("RNU4-2 position (n.)")
        if i == 0:
            a.set_title("B", loc="left", fontweight="bold")
    save(fig, "fig3_rnu4_2_sge")


def fig4():
    cl = pd.read_csv(os.path.join(HERE, "results", "alignment_species_clades.tsv"), sep="\t")
    pg = pd.read_csv(os.path.join(HERE, "results", "per_gene_rfam.tsv"), sep="\t")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw={"width_ratios": [1.2, 1], "wspace": 0.4})
    xs = np.arange(len(cl))
    a.bar(xs, cl.primates, width=0.6, color=BLUE, label="Primates")
    a.bar(xs, cl.non_primates, bottom=cl.primates + 3, width=0.6, color=ORANGE, label="Non-primate mammals")
    a.set_xticks(xs); a.set_xticklabels(cl.gene, rotation=35, ha="right")
    a.set_ylabel("Species aligned (≥ 50% of gene)\nin the 447-way alignment")
    a.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    a.text(-0.02, 1.12, "A", transform=a.transAxes, fontweight="bold", fontsize=10, ha="right")
    a.grid(axis="y", color=GRID, lw=0.6); a.set_axisbelow(True)
    b.scatter(pg.rho_phyloP447_vs_rfam, pg["ref - CADD"], s=40, color=BLUE, edgecolor=SURF, linewidth=1.2, zorder=3)
    b.errorbar(pg.rho_phyloP447_vs_rfam, pg["ref - CADD"], yerr=[pg["ref - CADD"] - pg["ref - CADD lo"],
               pg["ref - CADD hi"] - pg["ref - CADD"]], fmt="none", ecolor=BLUE, elinewidth=1, alpha=0.6)
    off = {"RNU4-2": (6, 4), "RNU2-2": (6, 6), "RNU12": (-36, -2), "RNU6ATAC": (6, -12), "RNU4ATAC": (-24, 9)}
    for r in pg.itertuples():
        b.annotate(r.gene, (r.rho_phyloP447_vs_rfam, getattr(r, "_19")), xytext=off.get(r.gene, (5, 3)),
                   textcoords="offset points", fontsize=7, color=INK)
    b.axhline(0, color=INK2, lw=0.8)
    b.set_xlabel("phyloP447 agreement with Rfam profile (ρ)"); b.set_ylabel("AUROC, Rfam ref − CADD (95% CI)")
    b.text(-0.02, 1.12, "B", transform=b.transAxes, fontweight="bold", fontsize=10, ha="right")
    b.grid(color=GRID, lw=0.6); b.set_axisbelow(True)
    save(fig, "fig4_alignment_and_genes")


def main():
    os.makedirs(OUT, exist_ok=True)
    fig1(); fig2(); fig3(); fig4()
    print("written:", sorted(os.listdir(OUT)))


if __name__ == "__main__":
    main()
