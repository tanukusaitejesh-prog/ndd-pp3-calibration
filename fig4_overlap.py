"""Figure 4: direct overlap test. R = LR(PP3 and flag) / [LR(PP3) x LR(flag)], with 95% bootstrap CIs."""
import os
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(HERE, "a10_overlap_main.tsv"), sep="\t")
d = d[d["flag"].isin(["PM1-like", "PM1-strict", "PP2-like"])]
cols = {"DD": "#2a78d6", "autism": "#7c5cc4", "public3": "#1baf7a"}
lab = {"DD": "DD (primary)", "autism": "Autism (replication)", "public3": "Public cohort (descriptive)"}
groups = [(t, f) for t in ["AlphaMissense", "REVEL", "VEST4"] for f in ["PM1-like", "PM1-strict", "PP2-like"]]
fig, ax = plt.subplots(figsize=(7.5, 7.2))
y = 0
ticks, tl = [], []
for t, f in groups:
    for k, c in enumerate(["DD", "autism", "public3"]):
        r = d[(d.tool == t) & (d.flag == f) & (d.cohort == c)]
        if not len(r):
            continue
        r = r.iloc[0]
        yy = y - k * 0.25
        ax.plot([r.R_lo95, r.R_hi95], [yy, yy], color=cols[c], lw=1.6)
        ax.plot(r.R, yy, "o", color=cols[c], ms=5, label=lab[c] if (t, f) == groups[0] else None)
    ticks.append(y - 0.25); tl.append(f"{t} × {f}")
    y -= 1
ax.axvline(1, color="#52514e", lw=0.9, ls=(0, (4, 3)))
ax.text(1.01, 0.4, "additive (R = 1)", color="#52514e", fontsize=7.5, va="bottom")
ax.set_yticks(ticks); ax.set_yticklabels(tl, fontsize=8)
ax.set_xlabel("R = joint LR / (LR of PP3 × LR of flag)")
ax.set_xlim(0.2, 1.5)
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig4_overlap.png"), dpi=200)
print("written fig4_overlap.png")
