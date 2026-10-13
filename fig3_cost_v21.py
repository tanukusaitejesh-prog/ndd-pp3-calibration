"""Figure 3 (draft 21): class-level calibration of PS2 + PM2 + PM1-like + AlphaMissense PP3 in DD, replotted from
results_fullstudy/cost_headline.tsv (written by make_cost_headline.py). The 'conditional' rule gives classes identical to
the capped rule and is omitted; the legend sits above the plot so it cannot overlap the data."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
C = pd.read_csv(os.path.join(HERE, "results_fullstudy", "cost_headline.tsv"), sep="\t")
C = C[C.rule != "conditional"].iloc[::-1].reset_index(drop=True)
INK, INK2, GRID, BLUE, ORANGE = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(7.2, 2.7))
lab = [f"{'No cap' if r.rule == 'uncapped' else 'ClinGen cap'}: {r.cls}\n(n = {r.n:,})" for r in C.itertuples()]
for i, r in enumerate(C.itertuples()):
    hi = 1.0 if r.cls == "Pathogenic" else 0.99
    ax.barh(i, hi - r.nominal_min_posterior, left=r.nominal_min_posterior, height=0.5, color="#c0d7f2", zorder=2)
    ax.errorbar(r.observed, i, xerr=[[r.observed - r.obs_lo], [r.obs_hi - r.observed]], fmt="o", color=ORANGE,
                ms=7, capsize=3, zorder=3, mec="white", mew=1.5)
    ax.text(r.obs_lo - 0.02, i, f"{r.observed:.0%}", va="center", ha="right", color=INK, fontsize=8.5)
ax.set_yticks(range(len(C)), lab)
ax.set_xlim(0, 1.0); ax.xaxis.grid(True, color=GRID, lw=0.6, zorder=0)
ax.set_xlabel("Disease-associated fraction of the class (de novo excess, DD)")
ax.scatter([], [], color=ORANGE, label="observed fraction (95% CI)")
ax.barh([], [], color="#c0d7f2", label="nominal posterior of the class")
ax.legend(frameon=False, fontsize=8.5, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
fig.tight_layout()
out = os.path.join(HERE, "figures", "fig3_cost_v21.png")
fig.savefig(out, dpi=300)
print("written", out)
