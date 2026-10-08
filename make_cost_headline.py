"""Draft 14, roadmap section 1: the clinical cost of stacking, in posterior percentage points and in points.
Inputs (saved v13 outputs; no new data): a12_classify_summary.txt (DD class counts O and background E under PS2 + PM2 +
PM1-like + AlphaMissense PP3, uncapped / capped / conditional), a10_overlap.tsv (pairwise R), triple_summary.txt (R3).
Class-level disease-associated fraction = 1 - E/O (de novo excess); 95% interval from an exact Poisson interval on O
with E fixed. Nominal class posteriors use the ClinGen points scale (prior 0.0441, Supporting LR 2.406): Pathogenic
>= 10 points (>= 0.99), Likely Pathogenic 6-9 (0.90-0.99). Overstatement in points from dependence = -ln R / ln 2.406.
Outputs: results_fullstudy/cost_headline.tsv, cost_headline_summary.txt, figures/fig_cost_calibration.png"""
import math
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import chi2

RES, FIG = "results_fullstudy/", "figures/"
PRIOR = 0.0441
PO = PRIOR / (1 - PRIOR)
LN = math.log(2.406)


def post(points):
    o = PO * 2.406 ** points
    return o / (1 + o)


def points_of(p):
    return math.log(p / (1 - p) / PO) / LN


def pois_ci(o, a=0.05):
    lo = chi2.ppf(a / 2, 2 * o) / 2 if o > 0 else 0.0
    hi = chi2.ppf(1 - a / 2, 2 * o + 2) / 2
    return lo, hi


rows = []
for line in open(RES + "a12_classify_summary.txt"):
    m = re.match(r"(\w+)\s+(Pathogenic|Likely pathogenic)\s*:\s+(\d+) de novos, expected background\s+([\d.]+)", line)
    if m and int(m.group(3)) > 0:
        rule, cls, O, E = m.group(1), m.group(2), int(m.group(3)), float(m.group(4))
        lo, hi = pois_ci(O)
        frac, f_lo, f_hi = 1 - E / O, 1 - E / lo, 1 - E / hi
        nom_lo, nom_pts = (0.99, 10) if cls == "Pathogenic" else (0.90, 6)
        rows.append(dict(rule=rule, cls=cls, n=O, background=E, observed=frac, obs_lo=f_lo, obs_hi=f_hi,
                         nominal_min_posterior=nom_lo, shortfall_pp=100 * (nom_lo - frac),
                         points_implied=points_of(frac), points_nominal_min=nom_pts,
                         overstatement_points=nom_pts - points_of(frac)))
C = pd.DataFrame(rows)
C.to_csv(RES + "cost_headline.tsv", sep="\t", index=False)

ov = pd.read_csv(RES + "a10_overlap.tsv", sep="\t")
pairs = ov[(ov.tool == "AlphaMissense") & ov.flag.isin(["PM1-like", "PP2-like"])].copy()
pairs["points_lost"] = -pairs.R.apply(math.log) / LN
trip = []
for line in open(RES + "triple_summary.txt"):
    m = re.match(r"(\w+)\s+AlphaMissense (\w+)\s+x (PM1-\w+)\s+x PP2-like: R3 ([\d.]+) \[([\d.]+), ([\d.]+)\]", line)
    if m:
        trip.append(dict(cohort=m.group(1), pp3=m.group(2), pm1=m.group(3), R3=float(m.group(4)),
                         points_lost=-math.log(float(m.group(4))) / LN,
                         points_lost_lo=-math.log(float(m.group(6))) / LN, points_lost_hi=-math.log(float(m.group(5))) / LN))
T = pd.DataFrame(trip)

out = ["Class-level calibration of PS2 + PM2 + PM1-like + AlphaMissense PP3 in DD (de novo excess):"]
for r in C.itertuples():
    out.append(f"  {r.rule:11s} {r.cls:17s} n {r.n:5d}: observed {r.observed:.2f} [{r.obs_lo:.2f}, {r.obs_hi:.2f}] vs nominal "
               f">= {r.nominal_min_posterior:.2f} -> shortfall {r.shortfall_pp:.0f} percentage points; implied {r.points_implied:.1f} "
               f"points vs >= {r.points_nominal_min} ({r.overstatement_points:.1f} points overstated)")
out.append("Points lost to dependence (-ln R / ln 2.406), AlphaMissense:")
for r in pairs.itertuples():
    out.append(f"  {r.cohort:6s} PP3 {r.pp3_level:20s} x {r.flag:9s}: R {r.R:.2f} -> {r.points_lost:.2f} points over-counted")
for r in T.itertuples():
    out.append(f"  {r.cohort:6s} PP3 {r.pp3:10s} x {r.pm1} x PP2-like: R3 {r.R3:.2f} -> {r.points_lost:.2f} points "
               f"[{r.points_lost_lo:.2f}, {r.points_lost_hi:.2f}] over-counted")
open(RES + "cost_headline_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))

# figure: nominal class posterior band vs observed class-level fraction
SURF, INK, INK2, GRID, BLUE, ORANGE = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(7.2, 2.6))
lab = [f"{r.rule.capitalize()}: {r.cls}\n(n = {r.n:,})" for r in C.itertuples()]
for i, r in enumerate(C.itertuples()):
    hi = 1.0 if r.cls == "Pathogenic" else 0.99
    ax.barh(i, hi - r.nominal_min_posterior, left=r.nominal_min_posterior, height=0.5, color=BLUE, alpha=0.25, zorder=2)
    ax.errorbar(r.observed, i, xerr=[[r.observed - r.obs_lo], [r.obs_hi - r.observed]], fmt="o", color=ORANGE,
                ms=7, capsize=3, zorder=3, mec=SURF, mew=1.5)
    ax.text(r.obs_lo - 0.02, i, f"{r.observed:.0%}", va="center", ha="right", color=INK, fontsize=8)
ax.set_yticks(range(len(C)), lab)
ax.set_xlim(0, 1.0); ax.xaxis.grid(True, color=GRID, lw=0.6, zorder=0)
ax.set_xlabel("Disease-associated fraction of the class (de novo excess, DD)")
ax.scatter([], [], color=ORANGE, label="observed (95% CI)"); ax.barh([], [], color=BLUE, alpha=0.25, label="nominal posterior band")
ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig(FIG + "fig_cost_calibration.png", dpi=300)
