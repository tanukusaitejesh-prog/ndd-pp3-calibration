"""Figure 2: points of independent evidence earned by PP3 once PM1 applies, against points applied, by data source and check.
Values from results_fullstudy/{stacking_cv_summary, erepo_eligibility_summary, lab_stacking_summary}.txt.
Run: python fig_credit.py -> fig2_credit.png"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EARN, APPL, INK, INK2, GRID, SURF = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#ebeae6", "#fcfcfb"
rows = [  # label, earned, lo, hi, applied, n (joint), group
    ("All PP3 strengths", 0.60, 0.14, 1.05, 1.24, 304, "Expert panels"),
    ("Supporting PP3", 0.81, 0.24, 1.68, 1.00, 235, "Expert panels"),
    ("Moderate PP3", 0.02, -1.69, 0.89, 2.00, 67, "Expert panels"),
    ("De novo criteria also removed from truth", 0.58, 0.08, 1.04, 1.27, 257, "Expert panels"),
    ("Classifications approved 2025-2026 (held out)", 0.40, -0.25, 1.07, 1.40, 67, "Expert panels"),
    ("PP3 defined by score, not by panel choice", 0.71, 0.35, 1.20, 1.00, 322, "Expert panels"),
    ("All PP3 strengths", 0.85, 0.53, 1.20, 1.09, 1775, "Clinical laboratories"),
    ("Supporting PP3", 0.83, 0.53, 1.24, 1.00, 1681, "Clinical laboratories"),
]
fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=200)
fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
y = []
pos, last = 0, None
for i, r in enumerate(rows):
    if r[6] != last:
        pos += 0.8 if last else 0
        ax.text(-5.3, pos - 0.05, r[6], ha="left", va="center", fontsize=8, fontweight="bold", color=INK)
        pos += 0.9
        last = r[6]
    y.append(pos); pos += 1
for (lab, e, lo, hi, ap, n, g), yy in zip(rows, y):
    ax.plot([max(lo, -1.6), min(hi, 3.0)], [yy, yy], color=EARN, lw=2, solid_capstyle="round", zorder=3)
    if lo < -1.6:
        ax.scatter([-1.6], [yy], marker="<", s=30, color=EARN, zorder=4)
    ax.scatter([e], [yy], s=36, color=EARN, edgecolor=SURF, linewidth=1.5, zorder=4)
    ax.scatter([ap], [yy], s=42, marker="D", facecolor=SURF, edgecolor=APPL, linewidth=2, zorder=5)
    ax.text(-5.3, yy, f"{lab} (n {n})", ha="left", va="center", fontsize=7.2, color=INK2)
ax.axvline(0, color="#c9c8c3", lw=0.8, zorder=1)
ax.axvspan(-5.35, -1.85, color=SURF, zorder=2)
ax.set_xlim(-5.35, 3.0)
ax.set_xticks([-1, 0, 1, 2, 3])
ax.spines["bottom"].set_bounds(-1.7, 3.0)
ax.set_xlabel("Points of independent evidence (ClinGen scale, 2.406 per point)", fontsize=8, color=INK)
ax.set_ylim(pos + 0.2, -0.4)
ax.set_yticks([])
ax.tick_params(axis="x", labelsize=7, colors=INK2)
for sp in ("top", "right", "left"):
    ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_color("#c9c8c3")
ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
ax.scatter([], [], s=36, color=EARN, label="earned beside PM1 (95% CI, gene-clustered)")
ax.scatter([], [], s=42, marker="D", facecolor=SURF, edgecolor=APPL, linewidth=2, label="points applied")
ax.legend(loc="lower center", bbox_to_anchor=(0.6, -0.27), ncol=2, frameon=False, fontsize=7, labelcolor=INK2)
ax.set_title("PP3 applied beside PM1: points earned vs points given", fontsize=9, color=INK, loc="left", x=-0.0)
fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.2)
fig.savefig("fig2_credit.png", facecolor=SURF)
print("saved")
