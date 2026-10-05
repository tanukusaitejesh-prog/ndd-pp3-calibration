"""Figure 1: study design schematic (nested conditioning on gene, constraint region and protein window)."""
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7c5cc4"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 9})
fig = plt.figure(figsize=(12.5, 5.2))

# ---- panel a: the estimator
ax = fig.add_axes([0.02, 0.06, 0.40, 0.88])
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
ax.text(0, 9.6, "a", fontsize=13, fontweight="bold")


def box(x, y, w, h, txt, col, fs=8.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.15", fc="white", ec=col, lw=1.4))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, color=INK)


box(0.3, 7.3, 4.2, 1.7, "Observed de novo missense\n31,058 DD trios\n(replication: 12,316 autism)", S[0])
box(5.4, 7.3, 4.3, 1.7, "Expected under mutation alone\nRoulette per-site rate × c\nc from synonymous de novos", MUTED)
box(2.2, 4.4, 5.6, 1.9, "Excess X = O − E in each score group\nLR = share of excess / share of expected\n(no clinical labels)", S[1])
box(0.3, 1.0, 4.4, 2.3, "Marginal (A)\nall sites compared\nacross the genome", S[3])
box(5.3, 1.0, 4.4, 2.3, "Within unit (B–E)\neach de novo compared only\nwith sites in the same\ngene / region / window", S[2])
for (x0, y0, x1, y1) in [(2.4, 7.3, 4.3, 6.35), (7.5, 7.3, 5.7, 6.35), (4.0, 4.4, 2.6, 3.35), (6.0, 4.4, 7.4, 3.35)]:
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.2))
ax.text(5, 0.25, "retention share r = ln LR(level) / ln LR(A)", ha="center", fontsize=9, color=INK, style="italic")

# ---- panel b: nested units on one protein
ax = fig.add_axes([0.46, 0.06, 0.53, 0.88])
L = 400
ax.set_xlim(-95, L + 10); ax.set_ylim(-0.6, 6.2); ax.axis("off")
ax.text(-95, 5.95, "b", fontsize=13, fontweight="bold")
rows = [("A  raw", None), ("B  gene", [(0, L)]), ("C  constraint region", [(0, 140), (140, 260), (260, L)]),
        ("D2  50-residue window", [(i, min(i + 50, L)) for i in range(0, L, 50)]),
        ("E2  window × region", sorted(set([(i, min(i + 50, L)) for i in range(0, L, 50)] +
                                           [(140, 150), (250, 260)])))]
cols = {1: S[0], 2: S[1], 3: S[2], 4: S[4]}
rng = np.random.default_rng(3)
for k, (lab, units) in enumerate(rows):
    y = 5 - k * 1.15
    ax.text(-8, y + 0.2, lab, ha="right", va="center", fontsize=9, color=INK)
    if units is None:
        ax.add_patch(Rectangle((0, y), L, 0.4, fc=GRID, ec="none"))
        ax.text(L / 2, y + 0.2, "all genes pooled", ha="center", va="center", fontsize=8, color=MUTED)
        continue
    cuts = sorted(set([u for ab in units for u in ab]))
    if k == 4:
        cuts = sorted(set(list(range(0, L + 1, 50)) + [140, 260]))
    for j in range(len(cuts) - 1):
        a_, b_ = cuts[j], cuts[j + 1]
        shade = 0.25 + 0.5 * ((j % 2) == 0)
        ax.add_patch(Rectangle((a_, y), b_ - a_, 0.4, fc=cols[k], alpha=shade, ec="white", lw=1))
    if k == 2:
        ax.text(200, y + 0.62, "RMC o/e 0.25 (PM1-like)", ha="center", fontsize=7, color=MUTED)
# de novo markers
for xv in rng.integers(150, 250, 4).tolist() + [40, 330]:
    ax.plot([xv, xv], [-0.35, 5.45], color=S[1], lw=0.6, alpha=0.35, zorder=0)
ax.text(L / 2, -0.55, "protein position (residues); orange lines = observed de novos", ha="center", fontsize=8, color=MUTED)
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig1_design.png"), dpi=200)
print("written fig1_design.png")
