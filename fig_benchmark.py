"""Figure 2: predictor ranks by marginal (pooled) vs within-window de novo evidence (bench2, DD, top 2%).
Run from the project folder: python fig_benchmark.py -> fig2_benchmark.png"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

R = "results_fullstudy/" if os.path.exists("results_fullstudy/bench2.tsv") else "results_fullstudy_backup/"
b = pd.read_csv(R + "bench2.tsv", sep="\t")
s = b[(b.stratum == "all_genes") & (b.site_set == "a1") & (b.q == 0.02) & (b.predictor != "NULL")]
s = s[s.LR_marginal >= 2.406].copy()
s["rm"] = s.LR_marginal.rank(ascending=False).astype(int)
s["rw"] = s.rank_window.astype(int)
s["d"] = s.rm - s.rw                                   # positive = rises once gene and position are fixed
UP, DOWN, GREY = "#2a78d6", "#eb6834", "#c9c8c3"
INK, INK2 = "#0b0b0b", "#52514e"
n = len(s)
fig, ax = plt.subplots(figsize=(6.4, 7.6), dpi=200)
fig.patch.set_facecolor("#fcfcfb")
ax.set_facecolor("#fcfcfb")
lab_up = s.nsmallest(4, "rw").predictor.tolist()                 # four best within windows
lab_dn = s.nsmallest(4, "rm").predictor.tolist() + ["AlphaMissense"]   # four best pooled, plus AlphaMissense
for _, r in s.sort_values("d", key=abs).iterrows():
    col = UP if r.d >= 10 else (DOWN if r.d <= -10 else GREY)
    lw = 2 if col != GREY else 1
    ax.plot([0, 1], [r.rm, r.rw], color=col, lw=lw, solid_capstyle="round", zorder=3 if col != GREY else 2)
    ax.scatter([0, 1], [r.rm, r.rw], s=18, color=col, edgecolor="#fcfcfb", linewidth=1.5, zorder=4)
    if r.predictor in lab_up + lab_dn:
        ax.text(-0.04, r.rm, f"{r.predictor}  {r.LR_marginal:.1f}", ha="right", va="center", fontsize=7, color=INK)
        ax.text(1.04, r.rw, f"{r.LR_window:.2f}  {r.predictor}", ha="left", va="center", fontsize=7, color=INK)
ax.set_ylim(n + 1, 0)
ax.set_xlim(-0.75, 1.75)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Rank by pooled\n(marginal) LR", "Rank by LR within\n50-residue windows"], fontsize=8, color=INK)
ax.set_yticks([1, 10, 20, 30, n])
ax.tick_params(axis="y", labelsize=7, colors=INK2, length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.grid(axis="y", color="#ebeae6", lw=0.6, zorder=0)
ax.plot([], [], color=UP, lw=2, label="rises ≥ 10 ranks")
ax.plot([], [], color=DOWN, lw=2, label="falls ≥ 10 ranks")
ax.plot([], [], color=GREY, lw=1, label="other")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False, fontsize=7, labelcolor=INK2)
ax.set_title("39 missense predictors, DD de novo excess (top 2% of scores)\nlabels: LR at each end", fontsize=8.5, color=INK, loc="center")
fig.tight_layout()
fig.savefig("fig2_benchmark.png", facecolor=fig.get_facecolor())
print("saved", n)
