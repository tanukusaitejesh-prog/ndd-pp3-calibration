import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = "/root/ndd_pp3_gate/results_fullstudy/"
S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # reference palette slots 1-4 (light)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
CUTS = [(2.406, "Supporting"), (5.79, "Moderate")]
from matplotlib.ticker import FuncFormatter, FixedLocator
def plainlog(ax, ticks):
    ax.yaxis.set_major_locator(FixedLocator(ticks)); ax.yaxis.set_minor_locator(FixedLocator([]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))

def refs(ax, xpos):
    for v, n in CUTS:
        ax.axhline(v, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.text(xpos, v * 1.04, f"PP3 {n}", color=MUTED, ha="left", va="bottom", fontsize=7.5)
    ax.axhline(1, color=GRID, lw=1, zorder=0)

fig, axs = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={"width_ratios": [1.25, 0.95, 1.15], "wspace": 0.55})
# A: Aim 1 local windows, all genes
w = pd.read_csv(R + "aim1_windows.tsv", sep="\t")
w = w[(w.score == "AlphaMissense") & (w.stratum == "all_genes")]
ax = axs[0]
series = [("gene", "marginal", "marginal (gene + variant)", S[3]), ("gene", "within_unit", "within gene", S[0]),
          ("rmc", "within_unit", "within gene × RMC region", S[1]), ("ccr", "within_unit", "within gene × CCR class", S[2])]
nudge = {"within gene × RMC region": 0.88, "within gene × CCR class": 1.12}
for unit, kind, lab, col in series:
    d = w[(w.unit == unit) & (w.kind == kind) & (w.centre_score >= 0.5)].sort_values("centre_score")
    ax.plot(d.centre_score, d.LR.clip(lower=0.05), color=col, lw=2, zorder=3)
    ax.text(d.centre_score.iloc[-1] + 0.006, d.LR.iloc[-1] * nudge.get(lab, 1), lab, color=INK, fontsize=7.5, va="center")
ax.set_yscale("log"); ax.set_xlim(0.5, 1.16); refs(ax, 0.52); plainlog(ax, [0.3, 0.5, 1, 2, 3, 5, 10]); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xlabel("AlphaMissense (local window centre)"); ax.set_ylabel("likelihood ratio (log scale)")
ax.set_title("A  Aim 1: evidence left after conditioning", loc="left", fontsize=10, color=INK)
# B: top bin, within-unit LR with CI, incl. exploratory pseudo-regions
c1 = pd.read_csv(R + "aim1_coarse.tsv", sep="\t"); cx = pd.read_csv(R + "aim1x_control_coarse.tsv", sep="\t")
rows = []
for st, lab_st in (("all_genes", "all genes"), ("ddg2p_monoallelic", "DDG2P")):
    for unit, src, lab in (("gene", c1, "gene"), ("pseudo", cx, "even pseudo-regions*"), ("rmc", c1, "RMC region"), ("ccr", c1, "CCR class")):
        r = src[(src.score == "AM_gate") & (src.score_lo == 0.99) & (src.kind == "within_unit") & (src.unit == unit) & (src.stratum == st)].iloc[0]
        rows.append((lab_st, lab, r.LR, r.lo95, r.hi95))
ax = axs[1]
labs = ["gene", "even pseudo-regions*", "RMC region", "CCR class"]
for j, (st, col) in enumerate((("all genes", S[0]), ("DDG2P", S[1]))):
    sub = [r for r in rows if r[0] == st]
    y = np.arange(len(labs)) + (j - 0.5) * 0.28
    for k, lab in enumerate(labs):
        r = [x for x in sub if x[1] == lab][0]
        ax.plot([r[3], r[4]], [y[k], y[k]], color=col, lw=2, zorder=2)
        ax.plot(r[2], y[k], "o", ms=8, color=col, mec="white", mew=2, zorder=3)
    ax.plot([], [], "o", color=col, ms=7, label=st)
ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs); ax.set_ylim(3.6, -0.75)
ax.axvline(2.406, color=MUTED, lw=0.8, ls=(0, (4, 3))); ax.text(2.43, -0.55, "PP3 Supporting", color=MUTED, fontsize=7.5)
ax.legend(frameon=False, loc="lower right", fontsize=8); ax.set_xlim(2.0, 4.1)
ax.set_xlabel("LR at AlphaMissense 0.99–1.00 (95% CI)")
ax.set_title("B  What drives the drop", loc="left", fontsize=10, color=INK)
ax.text(0.0, -0.2, "* exploratory control, added after the result", transform=ax.transAxes, fontsize=7, color=MUTED)
# C: Aim 2 mechanisms, within gene, DDG2P
w2 = pd.read_csv(R + "aim2_windows.tsv", sep="\t")
w2 = w2[(w2.score == "AlphaMissense") & (w2.kind == "within_unit")]
ax = axs[2]
for st, lab, col, nd in (("ddg2p_lof", "loss of function", S[0], 0.72), ("ddg2p_dn", "dominant negative", S[1], 1.3), ("ddg2p_gof", "gain of function", S[2], 1.22)):
    d = w2[(w2.stratum == st) & (w2.centre_score >= 0.5)].sort_values("centre_score")
    ax.plot(d.centre_score, d.LR.clip(lower=0.05), color=col, lw=2)
    ax.text(d.centre_score.iloc[-1] + 0.006, d.LR.iloc[-1] * nd, lab, color=INK, fontsize=7.5, va="center")
ax.set_yscale("log"); ax.set_xlim(0.5, 1.2); refs(ax, 0.52); plainlog(ax, [0.05, 0.1, 0.3, 1, 2, 3, 5]); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.text(0.0, -0.2, "small strata: windows with few de novos are noisy", transform=ax.transAxes, fontsize=7, color=MUTED)
ax.set_xlabel("AlphaMissense (local window centre)"); ax.set_ylabel("within-gene LR (log scale)")
ax.set_title("C  Aim 2: by mechanism (DDG2P genes)", loc="left", fontsize=10, color=INK)
fig.subplots_adjust(left=0.05, right=0.97, bottom=0.17, top=0.9, wspace=0.62); fig.savefig(R + "fig_aims_1_2.png", dpi=160, facecolor="#fcfcfb")
print("saved")
