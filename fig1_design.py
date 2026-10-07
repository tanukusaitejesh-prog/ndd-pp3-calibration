"""Figure 1: study design (draft 14): label-free dependence, independent truth, audit of practice."""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7c5cc4"]
INK, MUTED = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9})
fig = plt.figure(figsize=(13, 5.6))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 39); ax.set_ylim(0, 17); ax.axis("off")


def box(x, y, w, h, title, body, col):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1,rounding_size=0.3", fc="white", ec=col, lw=1.5))
    ax.text(x + w / 2, y + h - 0.45, title, ha="center", va="top", fontsize=9.5, fontweight="bold", color=INK)
    ax.text(x + w / 2, y + h - 1.35, body, ha="center", va="top", fontsize=8.2, color=INK, linespacing=1.35)


def arrow(x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.3))


cols = [(0.6, "a", "Label-free dependence", S[0]), (13.6, "b", "Independent truth", S[1]),
        (26.6, "c", "Audit of practice", S[2])]
for x, letter, head, col in cols:
    ax.text(x, 16.3, letter, fontsize=14, fontweight="bold", color=INK)
    ax.text(x + 0.9, 16.3, head, fontsize=11, fontweight="bold", color=col)

# a
box(0.6, 11.0, 11.8, 4.6, "De novo missense, three cohorts",
    "DD 31,058 trios · autism 12,316 probands\npublic 4,256 trios (deduplicated)\n"
    "expected E: Roulette rates scaled\nby synonymous de novos", S[0])
box(0.6, 5.6, 11.8, 4.4, "Likelihood ratio without labels",
    "LR(b) = (X_b / ΣX) / (E_b / ΣE),  X = O − E\nflags: PP3 (ClinGen thresholds),\n"
    "PM1-like (regional o/e, CCR), PP2-like (Z)", S[0])
box(0.6, 0.4, 11.8, 4.2, "Dependence",
    "R = LR(A∧B) / [LR(A)·LR(B)]\nmarginal, within DDG2P genes, within genes\n"
    "PP3 points alone vs beside PM1", S[0])
arrow(6.5, 11.0, 6.5, 10.05); arrow(6.5, 5.6, 6.5, 4.65)

# b
box(13.6, 11.0, 11.8, 4.6, "Truth not built from these criteria",
    "6,177 expert-panel + 87,873 laboratory\nmissense classifications re-scored without\n"
    "PP3 BP4 PM1 PP2 BP1 BP3 PM4 PM5 PS2 PM6", S[1])
box(13.6, 5.6, 11.8, 4.4, "Other truth",
    "ClinVar P/LP vs B/LB (partly circular)\nProteinGym functional assays\n"
    "(82 assays, 70 genes)", S[1])
box(13.6, 0.4, 11.8, 4.2, "Dependence on truth",
    "LR(C) = P(C | P/LP) / P(C | B/LB)\nR for PP3×PM1, PM1×PP2, PP3×PP2\n"
    "gene- and submitter-clustered; replication halves", S[1])
arrow(19.5, 11.0, 19.5, 10.05); arrow(19.5, 5.6, 19.5, 4.65)

# c
box(26.6, 11.0, 11.8, 4.6, "How criteria are combined",
    "co-application of PP3, PM1, PP2\nand PP3 strength beside PM1\n(panels: structured codes; labs: parsed text)", S[2])
box(26.6, 5.6, 11.8, 4.4, "Credit once PM1 applies",
    "points = ln[LR(C∧PM1) / LR(PM1)] / ln 2.406\nearned vs applied, by strength\n"
    "equivalence ±0.5 points", S[2])
box(26.6, 0.4, 11.8, 4.2, "Consequences",
    "re-score under four rules\n(ClinGen limit; triple limit;\nno PP2 beside PM1; PP3 ≤ Supporting)", S[2])
arrow(32.5, 11.0, 32.5, 10.05); arrow(32.5, 5.6, 32.5, 4.65)
arrow(25.5, 13.3, 26.5, 13.3)

fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig1_design.png"), dpi=200)
print("written fig1_design.png")
