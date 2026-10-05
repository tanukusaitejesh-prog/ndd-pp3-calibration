"""Manuscript figures 2 and 3 from committed result tables (run on Windows; reads the WSL results through \\\\wsl$)."""
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter  # noqa: E402

R = sys.argv[1] if len(sys.argv) > 1 else r"\\wsl$\ndd-gate\root\ndd_pp3_gate\results_fullstudy"
OUT = os.path.dirname(os.path.abspath(__file__))
S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7c5cc4"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
LEVELS = ["A", "B", "C", "D2", "E2"]
LEVEL_LAB = ["A\nraw", "B\ngene", "C\nconstraint\nregion", "D2\n50-residue\nwindow", "E2\nwindow ×\nconstraint"]


def nested(name, stratum=True):
    d = pd.read_csv(os.path.join(R, name), sep="\t")
    return d[d["stratum"] == "all_genes"] if (stratum and "stratum" in d) else d


def popeve():
    t = pd.read_csv(os.path.join(R, "a1v54_popeve_ddd_top2pct.tsv"), sep="\t")
    t = t[t["stratum"] == "all_genes"].set_index("unit")
    return {"A": t.loc["gene", "LRm_top2"], "B": t.loc["gene", "LRw_top2"], "C": t.loc["rmc", "LRw_top2"],
            "D2": t.loc["win50", "LRw_top2"], "E2": t.loc["win50_rmc", "LRw_top2"]}


dd = nested("a1_nested.tsv")
rep = nested(os.path.join("replication_affected", "rep_nested.tsv"), stratum=False)
v5 = nested("a1v54_nested.tsv", stratum=False)

# ---------------------------------------------------------------- Figure 2: evidence ladder
fig, ax = plt.subplots(figsize=(7.2, 4.2))
series = [("popEVE (DDD only)", popeve(), S[1]),
          ("AlphaMissense, DD", dd.set_index("score").loc["AlphaMissense"], S[0]),
          ("AlphaMissense, autism", rep.set_index("score").loc["AlphaMissense"], S[4]),
          ("MetaRNN, DD", v5.set_index("score").loc["MetaRNN"], S[2])]
x = np.arange(len(LEVELS))
for lab, r, col in series:
    y = [float(r[k]) for k in LEVELS]
    ax.plot(x, y, "-o", color=col, lw=2, ms=4.5, label=f"{lab}  ({y[0]:.1f} → {y[-1]:.2f})", zorder=3)
for v, n in [(2.406, "PP3 Supporting 2.406"), (5.79, "PP3 Moderate 5.79")]:
    ax.axhline(v, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(4.15, v * 1.05, n, color=MUTED, va="bottom", fontsize=7.5)
ax.axhline(1, color=GRID, lw=1, zorder=0)
ax.set_yscale("log")
ax.yaxis.set_major_locator(FixedLocator([1, 2, 5, 10, 20]))
ax.yaxis.set_minor_locator(FixedLocator([]))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
ax.set_ylim(0.9, 25)
ax.set_xlim(-0.2, 5.3)
ax.set_xticks(x)
ax.set_xticklabels(LEVEL_LAB, fontsize=8)
ax.set_ylabel("Likelihood ratio, top 2% of scores")
ax.legend(frameon=False, fontsize=7.5, loc="upper right", bbox_to_anchor=(0.86, 1.0))
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig2_evidence_ladder.png"), dpi=200)

# ---------------------------------------------------------------- Figure 3: per-tool retention, three settings
sets = [("DD (36 tools)", dd, S[0]), ("Autism replication", rep, S[4]), ("New tools, dbNSFP 5.4a", v5, S[2])]
fig, axs = plt.subplots(1, 3, figsize=(13, 5.2), sharey=False, gridspec_kw={"wspace": 0.55})
for ax, (title, d, col) in zip(axs, sets):
    inf = d[(d["A"] >= 2.406)].copy()
    inf = inf.sort_values("r_D2")
    yy = np.arange(len(inf))
    ax.hlines(yy, inf["r_E2"], inf["r_B"], color=GRID, lw=3, zorder=1)
    ax.scatter(inf["r_B"], yy, s=18, color=MUTED, label="gene (B)", zorder=2)
    ax.scatter(inf["r_D2"], yy, s=22, color=col, label="gene + position (D2)", zorder=3)
    ax.scatter(inf["r_E2"], yy, s=14, facecolors="white", edgecolors=INK, lw=0.8, label="+ constraint (E2)", zorder=4)
    ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (4, 3)))
    share = (inf.loc[inf["score"] != "AlphaMissense", "r_D2"] < 0.5).mean()
    ax.set_title(f"{title}\n{share:.0%} of informative tools keep < 50%\nafter gene + position", fontsize=9)
    ax.set_yticks(yy)
    ax.set_yticklabels(inf["score"], fontsize=7)
    ax.set_xlim(0, max(1.0, float(np.nanmax(inf["r_B"])) + 0.05))
    ax.set_xlabel("Share of log-evidence retained, r")
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=8, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
fig.subplots_adjust(left=0.08, right=0.98, top=0.84, bottom=0.17, wspace=0.55)
fig.savefig(os.path.join(OUT, "fig3_retention.png"), dpi=200)
print("written fig2_evidence_ladder.png, fig3_retention.png")

# ---------------------------------------------------------------- Figure 4: ClinVar clinical impact (A4.1)
imp_path = os.path.join(R, "a4_clinvar_impact.tsv")
if os.path.exists(imp_path):
    imp = pd.read_csv(imp_path, sep="\t")
    tools = list(dict.fromkeys(imp["tool"]))
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.2), gridspec_kw={"wspace": 0.35, "width_ratios": [1.45, 1]})
    ax = axs[0]
    rules = [("pm1like_mean_pub_capped", "published, PP3+PM1 capped at 4", S[3]),
             ("pm1like_mean_conditional", "conditional (published × r_E2)", S[0]),
             ("pm1like_mean_denovo_thr", "de novo E2 thresholds", MUTED)]
    groups = [(t, c) for t in tools for c in ("P/LP", "VUS")]
    xg = np.arange(len(groups))
    for j, (col, lab, colr) in enumerate(rules):
        vals = [float(imp[(imp.tool == t) & (imp.cls == c)][col].iloc[0]) for t, c in groups]
        ax.bar(xg + (j - 1) * 0.26, vals, 0.26, color=colr, label=lab)
    ax.axhline(4, color=MUTED, lw=0.8, ls=(0, (4, 3)))
    ax.text(len(groups) - 0.5, 4.05, "Strong (4)", ha="right", va="bottom", color=MUTED, fontsize=7.5)
    ax.set_xticks(xg)
    ax.set_xticklabels([f"{t}\n{c}" for t, c in groups], fontsize=7.5)
    ax.set_ylabel("Mean PP3 + PM1 points, PM1-like sites")
    ax.set_ylim(0, 4.6)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right", bbox_to_anchor=(1.0, 1.12), ncol=1)
    ax = axs[1]
    v = imp[imp.cls == "VUS"].set_index("tool").loc[tools]
    xt = np.arange(len(tools))
    ax.bar(xt - 0.2, v["ge6_uncapped"], 0.4, color=S[1], label="published, no cap")
    ax.bar(xt + 0.2, v["ge6_conditional"], 0.4, color=S[0], label="conditional")
    for i, n in enumerate(v["ge6_uncapped"]):
        ax.text(i - 0.2, n, f"{int(n):,}", ha="center", va="bottom", fontsize=7.5)
    for i, n in enumerate(v["ge6_conditional"]):
        ax.text(i + 0.2, n, f"{int(n):,}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(xt)
    ax.set_xticklabels(tools)
    ax.set_ylabel("ClinVar VUS reaching Likely Pathogenic\n(≥ 6 points) from PP3 + PM1 alone")
    ax.legend(frameon=False, fontsize=7.5)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.18)
    fig.savefig(os.path.join(OUT, "fig4_clinvar_impact.png"), dpi=200)
    print("written fig4_clinvar_impact.png")

# ---------------------------------------------------------------- Figure 5: null-contrast models (Aim 4, A5, A7)
top_path = os.path.join(R, "aim4v2_top.tsv")
if os.path.exists(top_path):
    tp = pd.read_csv(top_path, sep="\t")
    v3p = os.path.join(R, "aim4v3_top.tsv")
    if os.path.exists(v3p):
        t3 = pd.read_csv(v3p, sep="\t")
        tp = pd.concat([tp, t3[t3["score"] == "Aim4v3_GBDT"]], ignore_index=True)
    tp = tp[(tp["centre"] == "DDD") & (tp["top_pct"] == 2)]
    order = [s for s in ["AlphaMissense", "Aim4v1", "Aim4v2_GBDT", "Aim4v2_MLP", "Aim4v2_ENS", "Aim4v3_GBDT"]
             if s in set(tp["score"])]
    lab = {"AlphaMissense": "AlphaMissense", "Aim4v1": "v1", "Aim4v2_GBDT": "v2 GBDT", "Aim4v2_MLP": "v2 MLP",
           "Aim4v2_ENS": "v2 ensemble", "Aim4v3_GBDT": "v3 gene-matched"}
    colr = {"AlphaMissense": S[0], "Aim4v1": MUTED, "Aim4v2_GBDT": S[1], "Aim4v2_MLP": S[3], "Aim4v2_ENS": S[4],
            "Aim4v3_GBDT": S[2]}
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.3), gridspec_kw={"wspace": 0.3})
    for ax, (col, lo, hi, title) in zip(axs, [("LRm", "lo95", "hi95", "Marginal (across genes)"),
                                             ("LRw", "lo95_w", "hi95_w", "Within gene")]):
        for j, gs in enumerate(["all_genes", "cold_genes"]):
            d = tp[tp["genes"] == gs].set_index("score")
            xs = np.arange(len(order)) + (j - 0.5) * 0.36
            vals = [float(d.loc[s, col]) if s in d.index else np.nan for s in order]
            los = [float(d.loc[s, lo]) if s in d.index else np.nan for s in order]
            his = [float(d.loc[s, hi]) if s in d.index else np.nan for s in order]
            ax.bar(xs, vals, 0.34, color=[colr[s] for s in order], alpha=1.0 if j == 0 else 0.55,
                   label="all genes" if j == 0 else "cold (held-out) genes")
            ax.errorbar(xs, vals, yerr=[np.subtract(vals, los), np.subtract(his, vals)], fmt="none", ecolor=INK, lw=0.8,
                        capsize=2)
        ax.set_xticks(np.arange(len(order)))
        ax.set_xticklabels([lab[s] for s in order], rotation=25, ha="right", fontsize=8)
        ax.set_title(f"{title}: top-2% LR, held-out DDD trios", fontsize=9)
        ax.axhline(2.406, color=MUTED, lw=0.8, ls=(0, (4, 3)))
    axs[1].legend(frameon=False, fontsize=7.5)
    fig.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.22)
    fig.savefig(os.path.join(OUT, "fig5_models.png"), dpi=200)
    print("written fig5_models.png")
