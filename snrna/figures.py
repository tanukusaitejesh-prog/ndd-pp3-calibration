#!/usr/bin/env python3
"""Paper figures (PDF + PNG in paper/figures/). Each figure is drawn only if its inputs exist; missing ones are listed.

    python figures.py [--blind sealed/blind_labels.tsv]

  fig1_blind_roc        ROC on blind RNU2-2: frozen v4, Rfam alt rarity (label-free), phyloP447, CADD
  fig2_robustness       blind RNU2-2 AUROC of v4 vs position-only baselines and comparators (blind_robustness.py)
  fig3_ablation         v4 drop-one / single-feature AUROC on blind RNU2-2 (v4_ablation.py)
  fig4_calibration      pooled local LR of out-of-sample v4 with the PP3/BP4 evidence bars (v4_calibrate.py)
  fig5_a7               A7 per-unit AUROC: Rfam alt rarity vs CADD vs phyloP447 on non-snRNA RNA genes
  fig6_gene_discovery   loci ranked by within-locus window depletion (rel_z), disease genes marked
  fig7_aim2             RNU4-2 SGE Spearman by model version against the pre-registered bar
  fig8_a9               genome-wide rarity-only inflation of CADD, GPN-MSA, ReMM by variant class (A9)
Colours: validated categorical slots 1-4 (light mode, print); every series also has its own line style or marker
and a direct label, so identity never rests on colour alone.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snrna as S  # noqa: E402

OUT = os.path.join(HERE, "paper", "figures")
RES = os.path.join(HERE, "results")
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
SERIES = {"v4": ("#2a78d6", "-", "v4 (frozen)"), "rfam": ("#1baf7a", "--", "Rfam alt rarity"),
          "phylop447": ("#4a3aa7", "-.", "phyloP447"), "cadd": ("#eb6834", ":", "CADD")}


def setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                         "legend.frameon": False, "pdf.fonttype": 42, "savefig.dpi": 300})
    return plt


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight")
    print(f"  wrote paper/figures/{name}.pdf/.png")


def roc_xy(y, s):
    o = np.argsort(-s, kind="mergesort")
    y = y[o]
    tp = np.r_[0, np.cumsum(y)] / max(y.sum(), 1)
    fp = np.r_[0, np.cumsum(~y)] / max((~y).sum(), 1)
    return fp, tp


def fig_blind_roc(plt, blind):
    sp = os.path.join(HERE, "results_v4", "v4_scores.tsv.gz")
    if not (blind and os.path.exists(blind) and os.path.exists(sp)):
        return "needs --blind and results_v4/v4_scores.tsv.gz"
    var = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    v = var[var.gene == "RNU2-2"].merge(pd.read_csv(sp, sep="\t")[["ID", "v4"]], on="ID")
    v = v.merge(pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")[["ID", "phylop447"]], on="ID", how="left")
    b = pd.read_csv(blind, sep="\t")
    plp = set(b.ID[(b.gene == "RNU2-2") & (b.label == "P/LP")])
    pos = v.ID.isin(plp)
    v = v[pos | (v.L3_benign_like.astype(bool) & ~pos)]
    y = v.ID.isin(plp).to_numpy()
    cols = {"v4": v.v4, "rfam": 1 - v.alt_col_freq, "phylop447": v.phylop447, "cadd": v.cadd_phred}
    fig, ax = plt.subplots(figsize=(3.4, 3.3))
    ax.plot([0, 1], [0, 1], color=MUTED, lw=0.8, ls=(0, (2, 2)))
    for k, s in cols.items():
        s = s.fillna(s.median()).to_numpy(float)
        fp, tp = roc_xy(y, s)
        col, ls, lab = SERIES[k]
        ax.plot(fp, tp, color=col, ls=ls, lw=2, label=f"{lab}  {S.auroc(y, s):.2f}")
    ax.set(xlabel="False-positive rate (L3 population SNVs)", ylabel="True-positive rate (P/LP SNVs)", xlim=(0, 1), ylim=(0, 1.01))
    ax.set_title(f"Blind RNU2-2: {int(y.sum())} P/LP vs {int((~y).sum())} L3", loc="left", color=INK, fontsize=9)
    ax.legend(loc="lower right", title="AUROC", fontsize=8, title_fontsize=8)
    save(fig, "fig1_blind_roc")


def dot_ci(plt, df, label_col, title, ref=None, highlight=()):
    df = df.sort_values("auroc")
    fig, ax = plt.subplots(figsize=(4.6, 0.28 * len(df) + 0.9))
    yy = np.arange(len(df))
    for i, r in enumerate(df.itertuples()):
        hi = any(h in getattr(r, label_col) for h in highlight)
        col = SERIES["v4"][0] if hi else INK2
        if pd.notna(getattr(r, "lo", np.nan)):
            ax.plot([r.lo, r.hi], [i, i], color=col, lw=1.5, solid_capstyle="round")
        ax.plot(r.auroc, i, "o", color=col, ms=5, mec="#fcfcfb", mew=1)
        ax.text(1.005, i, f"{r.auroc:.2f}", va="center", fontsize=7, color=INK2, transform=ax.get_yaxis_transform())
    ax.set_yticks(yy, df[label_col])
    ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (2, 2)))
    if ref is not None:
        ax.axvline(ref[0], color=SERIES["cadd"][0], lw=1, ls=":")
        ax.text(ref[0], len(df) - 0.4, ref[1], fontsize=7, color=INK2, ha="center")
    ax.set_xlabel("AUROC (95% CI)")
    ax.set_title(title, loc="left", color=INK, fontsize=9)
    ax.grid(axis="y", visible=False)
    return fig


def fig_robustness(plt):
    p = os.path.join(RES, "blind_robustness.tsv")
    if not os.path.exists(p):
        return "needs results/blind_robustness.tsv (blind_robustness.py)"
    t = pd.read_csv(p, sep="\t")
    t = t[(t.negatives == "L3 AC>=9 (frozen)") & (t.positions == "all positions")].dropna(subset=["auroc"])
    fig = dot_ci(plt, t, "score", "Blind RNU2-2: v4 vs position-only baselines", highlight=("v4 (frozen)",))
    save(fig, "fig2_robustness")


def fig_ablation(plt):
    p = os.path.join(RES, "v4_ablation.tsv")
    if not os.path.exists(p):
        return "needs results/v4_ablation.tsv (v4_ablation.py)"
    t = pd.read_csv(p, sep="\t").rename(columns={"auroc_RNU2-2": "auroc"})
    cadd = t.loc[t.model.str.startswith("CADD"), "auroc"]
    t = t[~t.model.str.startswith("CADD")]
    fig = dot_ci(plt, t, "model", "v4 ablation on blind RNU2-2 (post hoc)",
                 ref=(float(cadd.iloc[0]), "CADD") if len(cadd) else None, highlight=("v4 full", "col_freq", "alt_rarity"))
    save(fig, "fig3_ablation")


def fig_calibration(plt):
    p = os.path.join(RES, "v4_calibration_lr_curves.tsv.gz")
    if not os.path.exists(p):
        return "needs results/v4_calibration_lr_curves.tsv.gz (v4_calibrate.py)"
    c = pd.read_csv(p, sep="\t")
    c = c[(c.unit == "pooled (all genes)") & (c.score == "v4")].sort_values("s")
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    col = SERIES["v4"][0]
    ax.fill_between(c.s, c.lr_lo5.clip(lower=1e-3), c.lr_hi95, color=col, alpha=0.18, lw=0, label="5th-95th bootstrap")
    ax.plot(c.s, c.lr, color=col, lw=2, label="local LR")
    th = S.lr_thresholds(0.10)
    for k, ls in (("Supporting", "--"), ("Moderate", ":")):
        for val, side in ((th[k], "PP3"), (1 / th[k], "BP4")):
            ax.axhline(val, color=INK2, lw=0.8, ls=ls)
            ax.text(c.s.min(), val * 1.08, f"{side} {k}", fontsize=7, color=INK2)
    ax.set_yscale("log")
    ax.set(xlabel="v4 score (out-of-sample)", ylabel="Local likelihood ratio (P/LP : L3)")
    ax.set_title("ClinGen local-LR calibration of v4, pooled (prior 0.10)", loc="left", color=INK, fontsize=9)
    ax.legend(loc="upper left", fontsize=7)
    save(fig, "fig4_calibration")


def fig_a7(plt):
    p = os.path.join(RES, "a7_results.tsv")
    if not os.path.exists(p):
        return "needs results/a7_results.tsv (a7_ncrna.py eval)"
    t = pd.read_csv(p, sep="\t")
    t = t[t.set.str.startswith("PRIMARY")].dropna(subset=["auroc_rfam_alt_rarity"])
    keys = [("rfam_alt_rarity", "rfam"), ("phylop447", "phylop447"), ("cadd_phred", "cadd")]
    fig, ax = plt.subplots(figsize=(4.6, 0.55 * len(t) + 1.0))
    for i, r in enumerate(t.itertuples()):
        for j, (c, k) in enumerate(keys):
            col, _, lab = SERIES[k]
            y = i + (j - 1) * 0.22
            ax.plot([getattr(r, f"{c}_lo"), getattr(r, f"{c}_hi")], [y, y], color=col, lw=1.5)
            ax.plot(getattr(r, f"auroc_{c}"), y, "osD"[j], color=col, ms=5, mec="#fcfcfb", mew=1, label=lab if i == 0 else None)
    ax.set_yticks(range(len(t)), [f"{u} ({p} P/LP)" for u, p in zip(t.unit, t.n_P)])
    ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (2, 2)))
    ax.set_xlabel("AUROC (95% CI)")
    ax.set_title("A7: zero-shot on non-snRNA disease RNA genes", loc="left", color=INK, fontsize=9)
    ax.legend(loc="lower left", fontsize=7, bbox_to_anchor=(0, 1.02), ncol=3)
    ax.grid(axis="y", visible=False)
    save(fig, "fig5_a7")


def fig_gene_discovery(plt):
    p = os.path.join(RES, "gene_discovery_v3.tsv")
    p = p if os.path.exists(p) else os.path.join(RES, "gene_discovery.tsv")
    if not os.path.exists(p):
        return "needs results/gene_discovery*.tsv"
    t = pd.read_csv(p, sep="\t").sort_values("rel_z").reset_index(drop=True)
    dis = t.disease.fillna("").astype(str).str.contains("NDD")
    fig, ax = plt.subplots(figsize=(6.4, 2.6))
    x = np.arange(len(t))
    ax.vlines(x, 0, t.rel_z, color=np.where(dis, SERIES["v4"][0], GRID), lw=2)
    ax.plot(x[dis], t.rel_z[dis], "o", color=SERIES["v4"][0], ms=5, mec="#fcfcfb", mew=1)
    last = None
    for i in x[dis]:  # labels to the right; nudge down when two disease genes have near-equal depth
        yl = t.rel_z[i] if last is None or abs(t.rel_z[i] - last) > 0.45 else last + 0.45
        ax.text(i + 1.2, yl, t.gene[i], fontsize=7, color=INK, va="center")
        last = yl
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xticks([])
    ax.set(xlabel=f"{len(t)} snRNA loci, ranked by most-depleted 20-nt window", ylabel="Window depletion z (rel. to own locus)")
    ax.set_title("Within-locus constraint ranks the five NDD snRNA genes first (exploratory)", loc="left", color=INK, fontsize=9)
    ax.grid(axis="x", visible=False)
    save(fig, "fig6_gene_discovery")


def fig_aim2(plt):
    # values from PREREG s14 (v2 official evaluate.py; v3 exploratory; v4 post hoc, RNU4 cluster held out, CPU, no Evo 2)
    rows = [("v2 graph model\n(pre-registered)", 0.144), ("v3 (+UK Biobank)", 0.142), ("v4, RNU4 cluster\nheld out (post hoc)", 0.382),
            ("Rfam ref column\nfreq alone", 0.459)]
    fig, ax = plt.subplots(figsize=(4.4, 2.4))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, v) in zip(y, rows):
        ax.barh(yi, v, height=0.5, color=SERIES["v4"][0] if "Rfam" not in lab else SERIES["rfam"][0])
        ax.text(v + 0.01, yi, f"{v:.3f}", va="center", fontsize=7, color=INK2)
    ax.axvline(0.346, color=INK2, lw=1, ls="--")
    ax.text(0.352, y.max() + 0.42, "pre-registered bar 0.346", fontsize=7, color=INK2, ha="left", va="center")
    ax.set_yticks(y, [r[0] for r in rows], fontsize=7)
    ax.set(xlabel="Spearman ρ with RNU4-2 SGE (435 SNVs)", xlim=(0, 0.55), ylim=(-0.5, y.max() + 0.7))
    ax.set_title("Aim 2: concordance with measured function", loc="left", color=INK, fontsize=9)
    ax.grid(axis="y", visible=False)
    save(fig, "fig7_aim2")


def fig_a9(plt):
    p = os.path.join(RES, "a9_results.tsv")
    if not os.path.exists(p):
        return "needs results/a9_results.tsv (a9_inflation.py eval)"
    t = pd.read_csv(p, sep="	")
    units = [u for u in ("intronic", "5'UTR", "3'UTR", "ncRNA exon", "pooled") if u in set(t.unit)]
    tools = [("CADD v1.7", "cadd"), ("GPN-MSA", "v4"), ("ReMM v0.4", "phylop447")]
    fig, ax = plt.subplots(figsize=(5.2, 0.62 * len(units) + 1.0))
    for i, u in enumerate(units):
        for j, (lab, k) in enumerate(tools):
            r = t[(t.unit == u) & (t.tool == lab)]
            if not len(r):
                continue
            r = r.iloc[0]
            y = i + (j - 1) * 0.24
            col = SERIES[k][0]
            ax.plot([0.5, r.null_mean], [y, y], color=GRID, lw=5, solid_capstyle="butt")
            ax.plot([r.null_mean, r.raw_auroc], [y, y], color=col, lw=5, solid_capstyle="butt", label=lab if i == 0 else None)
            ax.plot(r.null_mean, y, "|", color=INK2, ms=9, mew=1.5)
    ax.set_yticks(range(len(units)), [f"{u} ({int(t[t.unit == u].n_P.iloc[0])} P/LP)" for u in units])
    ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (2, 2)))
    ax.set_xlim(0.45, 0.92)
    ax.set_xlabel("AUROC vs gnomAD-common SNVs (grey = rarity-only null)")
    ax.set_title("A9: rarity alone inflates non-coding predictors", loc="left", color=INK, fontsize=9, pad=22)
    ax.legend(loc="lower left", fontsize=7, bbox_to_anchor=(0, 1.0), ncol=3, borderaxespad=0.2)
    ax.grid(axis="y", visible=False)
    save(fig, "fig8_a9")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blind", default=os.path.join(HERE, "sealed", "blind_labels.tsv"))
    a = ap.parse_args()
    plt = setup()
    skipped = {}
    for name, fn in (("fig1_blind_roc", lambda: fig_blind_roc(plt, a.blind)), ("fig2_robustness", lambda: fig_robustness(plt)),
                     ("fig3_ablation", lambda: fig_ablation(plt)), ("fig4_calibration", lambda: fig_calibration(plt)),
                     ("fig5_a7", lambda: fig_a7(plt)), ("fig6_gene_discovery", lambda: fig_gene_discovery(plt)),
                     ("fig7_aim2", lambda: fig_aim2(plt)), ("fig8_a9", lambda: fig_a9(plt))):
        print(name, flush=True)
        why = fn()
        if why:
            skipped[name] = why
        plt.close("all")
    for k, v in skipped.items():
        print(f"  skipped {k}: {v}")


if __name__ == "__main__":
    main()
