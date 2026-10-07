"""Within-gene, label-free benchmark of missense predictors from committed de novo outputs (no new analysis of sites).
For each predictor (top 2% of its scores by mutational opportunity): marginal LR (A), within-gene LR (B) and
within-50-residue-window LR (D2) in DD; retention r = ln LR(level) / ln LR(A) with bootstrap intervals (A6);
replication of within-window LR in autism and the public cohort. Outputs supp_table7_benchmark.tsv and a summary.
Run from the project folder: python make_benchmark.py"""
import numpy as np, pandas as pd
from scipy import stats

import os
R = "results_fullstudy/" if os.path.exists("results_fullstudy/a6_retention_ci.tsv") else "results_fullstudy_backup/"
a1 = pd.read_csv(R + "a1_nested.tsv", sep="\t")
dd = a1[a1.stratum == "all_genes"].set_index("score")[["A", "B", "D2", "E2"]]
ddg = a1[a1.stratum == "ddg2p_monoallelic"].set_index("score")[["A", "D2"]].add_prefix("DDG2P_")
v54 = pd.read_csv(R + "a1v54_nested.tsv", sep="\t").set_index("score")[["A", "B", "D2", "E2"]]
dd = pd.concat([dd.assign(set="dbNSFP 4.1a + AlphaMissense"), v54.assign(set="dbNSFP 5.4a (newer)")])
ci = pd.read_csv(R + "a6_retention_ci.tsv", sep="\t")
cdd = ci[ci.set.isin(["aim1tools", "a1v54"])].set_index("score")
for c in ("r_B", "r_B_lo", "r_B_hi", "r_D2", "r_D2_lo", "r_D2_hi"):
    dd[c] = cdd[c].reindex(dd.index)
rep = ci[ci.set == "rep"].set_index("score")
dd["autism_A"] = rep["A"].reindex(dd.index)
dd["autism_D2"] = rep["A"].reindex(dd.index) ** rep["r_D2"].reindex(dd.index)
p3 = pd.read_csv(R + "public3/p3_nested.tsv", sep="\t").set_index("score")
dd["public_A"], dd["public_D2"] = p3["A"].reindex(dd.index), p3["D2"].reindex(dd.index)
dd = dd.join(ddg)
pe = pd.read_csv(R + "a1v54_popeve_ddd_top2pct.tsv", sep="\t")
pe = pe[(pe.score == "popEVE") & (pe.stratum == "all_genes")].set_index("unit")
inf = dd[dd.A >= 2.406].copy()                                   # informative: marginal top-2% LR >= 2.406
inf["D2_lo"] = inf.A ** inf.r_D2_lo
inf["D2_hi"] = inf.A ** inf.r_D2_hi
inf["rank_marginal"] = inf.A.rank(ascending=False).astype(int)
inf["rank_within_window"] = inf.D2.rank(ascending=False).astype(int)
inf["rank_change"] = inf.rank_marginal - inf.rank_within_window
inf["inflation"] = inf.A / inf.D2                                 # marginal LR / within-window LR
inf["points_marginal"] = np.log(inf.A) / np.log(2.406)
inf["points_within_window"] = np.log(inf.D2) / np.log(2.406)
inf = inf.sort_values("D2", ascending=False)
cols = ["set", "A", "B", "D2", "D2_lo", "D2_hi", "r_B", "r_D2", "r_D2_lo", "r_D2_hi", "points_marginal", "points_within_window",
        "inflation", "rank_marginal", "rank_within_window", "rank_change", "DDG2P_A", "DDG2P_D2", "autism_A", "autism_D2",
        "public_A", "public_D2"]
inf[cols].rename(columns={"A": "DD_marginal_LR", "B": "DD_within_gene_LR", "D2": "DD_within_window_LR"}).to_csv(
    "supp_table7_benchmark.tsv", sep="\t", float_format="%.3f")

rho_rank = stats.spearmanr(inf.A, inf.D2)
both_a = inf[["autism_D2", "D2"]].dropna()
both_p = inf[["public_D2", "D2"]].dropna()
both_am = inf[["autism_A", "A"]].dropna()
top5m = inf.sort_values("A", ascending=False).index[:5].tolist()
top5w = inf.index[:5].tolist()
out = [
    f"Informative predictors in DD (marginal top-2% LR >= 2.406): {len(inf)} ({(inf.set.str.startswith('dbNSFP 4')).sum()} dbNSFP 4.1a + AlphaMissense, "
    f"{(inf.set.str.startswith('dbNSFP 5')).sum()} newer)",
    f"Median marginal LR {inf.A.median():.2f}; within gene {inf.B.median():.2f}; within window {inf.D2.median():.2f}",
    f"Median share of log-evidence retained within window {inf.r_D2.median():.2f} (i.e. {1 - inf.r_D2.median():.2f} is gene and position)",
    f"Median inflation (marginal / within-window LR) {inf.inflation.median():.2f}; range {inf.inflation.min():.2f}-{inf.inflation.max():.2f}",
    f"Spearman(marginal LR, within-window LR) across predictors: {rho_rank[0]:.2f} (p = {rho_rank[1]:.2g})",
    f"Top 5 by marginal LR: {', '.join(top5m)}",
    f"Top 5 by within-window LR: {', '.join(top5w)}",
    f"Largest drops in rank: " + ", ".join(f"{k} ({int(v):+d})" for k, v in inf.rank_change.sort_values().head(5).items()),
    f"Largest gains in rank: " + ", ".join(f"{k} ({int(v):+d})" for k, v in inf.rank_change.sort_values(ascending=False).head(5).items()),
    f"Replication of within-window LR across predictors: autism Spearman {stats.spearmanr(both_a.autism_D2, both_a.D2)[0]:.2f} "
    f"(n {len(both_a)}); public Spearman {stats.spearmanr(both_p.public_D2, both_p.D2)[0]:.2f} (n {len(both_p)}); "
    f"marginal LR DD vs autism Spearman {stats.spearmanr(both_am.autism_A, both_am.A)[0]:.2f}",
    f"Predictors whose within-window LR reaches the ClinGen Supporting LR (2.406) with a lower 95% bound above it: "
    f"{', '.join(inf.index[inf.D2_lo > 2.406]) or 'none'}",
]
if len(pe):
    out.append(f"popEVE (DDD only): marginal {pe.LRm_top2.iloc[0]:.2f}; within gene {pe.loc['gene', 'LRw_top2']:.2f}; "
               + "; ".join(f"{u} {pe.loc[u, 'LRw_top2']:.2f}" for u in pe.index if u != "gene"))
open(R + "benchmark_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
print(inf[["A", "D2", "D2_lo", "D2_hi", "r_D2", "rank_marginal", "rank_within_window"]].round(2).to_string())
