"""Exploratory sensitivity (not pre-registered; addresses limitation 'top 2% only'): retention shares at the top 1%,
2% and 10% of mutational opportunity, rebuilt exactly from the committed coarse tables (no new estimation).
Within-unit LR over a band = sum(X) / sum(X / LR_b); marginal LR = (sum X_top / sum X_all) / (sum E_top / sum E_all)."""
import os
import sys

import numpy as np
import pandas as pd

R = sys.argv[1] if len(sys.argv) > 1 else "/root/ndd_pp3_gate/results_fullstudy"
SUPP = 2.406
LEVELS = {"B": "gene", "C": "rmc", "D2": "win50", "E2": "win50_rmc"}
TOPS = {"top1": 1, "top2": 2, "top10": 3}            # number of highest coarse bins: 99-100, 98-100, 90-100


def lr_top(g, n_top, kind):
    g = g.sort_values("score_lo")
    top = g.iloc[-n_top:]
    if kind == "marginal":
        return (top["X"].sum() / g["X"].sum()) / (top["E"].sum() / g["E"].sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        d = top["X"] / top["LR"]
    return top["X"].sum() / d.sum() if np.all(np.isfinite(d)) and d.sum() > 0 else np.nan


rows = []
for name, path in (("DD", "aim1tools_coarse.tsv"), ("autism", os.path.join("replication_affected", "rep_coarse.tsv")),
                   ("A3", "a1v54_coarse.tsv")):
    c = pd.read_csv(os.path.join(R, path), sep="\t")
    c = c[c["stratum"] == "all_genes"]
    for s, gs in c.groupby("score", sort=False):
        for tname, n in TOPS.items():
            gm = gs[(gs["unit"] == "gene") & (gs["kind"] == "marginal")]
            if len(gm) < 4:
                continue
            A = lr_top(gm, n, "marginal")
            r = {"set": name, "score": s, "top": tname, "A": A}
            for lvl, u in LEVELS.items():
                gw = gs[(gs["unit"] == u) & (gs["kind"] == "within_unit")]
                v = lr_top(gw, n, "within") if len(gw) >= 4 else np.nan
                r[lvl] = v
                r["r_" + lvl] = np.log(v) / np.log(A) if (A > 1 and v > 0) else np.nan
            rows.append(r)
t = pd.DataFrame(rows)
t.to_csv(os.path.join(R, "sens_threshold_retention.tsv"), sep="\t", index=False, float_format="%.4g")
out = []
for (name, tname), g in t.groupby(["set", "top"], sort=False):
    inf = g[(g["score"] != "AlphaMissense") & (g["A"] >= SUPP)]
    am = g[g["score"] == "AlphaMissense"]
    amtxt = (f"; AlphaMissense A {am['A'].iloc[0]:.2f} B {am['B'].iloc[0]:.2f} D2 {am['D2'].iloc[0]:.2f} "
             f"E2 {am['E2'].iloc[0]:.2f}") if len(am) else ""
    out.append(f"{name:6s} {tname:5s} informative {len(inf):2d}; share r(D2) < 0.5 = {(inf['r_D2'] < 0.5).mean():.2f}; "
               f"median r B {inf['r_B'].median():.2f} D2 {inf['r_D2'].median():.2f} E2 {inf['r_E2'].median():.2f}; "
               f"E2 >= 2.406 (point) {int((g['E2'] >= SUPP).sum())}{amtxt}")
txt = "\n".join(out)
open(os.path.join(R, "sens_threshold_summary.txt"), "w").write(
    "Exploratory, not pre-registered: retention by top-score threshold (rebuilt from committed coarse tables)\n" + txt + "\n")
print(txt)
