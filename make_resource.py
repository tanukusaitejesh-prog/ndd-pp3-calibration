"""Supplementary Table 1: per-predictor conditional-evidence resource, assembled from committed result tables only
(no new estimation). Run in the analysis environment: python make_resource.py [results_dir]."""
import os
import sys

import numpy as np
import pandas as pd

R = sys.argv[1] if len(sys.argv) > 1 else "/root/ndd_pp3_gate/results_fullstudy"
HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = 2.406


def nested(path, stratum=True):
    d = pd.read_csv(os.path.join(R, path), sep="\t")
    return d[d["stratum"] == "all_genes"] if (stratum and "stratum" in d) else d


def popeve():
    t = pd.read_csv(os.path.join(R, "a1v54_popeve_ddd_top2pct.tsv"), sep="\t")
    t = t[t["stratum"] == "all_genes"].set_index("unit")
    A = t.loc["gene", "LRm_top2"]
    row = {"score": "popEVE", "A": A}
    for lvl, u in (("B", "gene"), ("D2", "win50"), ("E2", "win50_rmc")):
        row[lvl] = t.loc[u, "LRw_top2"]
        row[lvl + "_lo_sys"] = t.loc[u, "lo_sys_w"]
        row["r_" + lvl] = np.log(row[lvl]) / np.log(A)
    return pd.DataFrame([row])


cols = ["score", "A", "B", "D2", "E2", "E2_lo_sys", "r_B", "r_D2", "r_E2"]
parts = [nested("a1_nested.tsv")[cols].assign(set="dbNSFP 4.1a + AlphaMissense", cohort_DD="DD (31,058 trios)"),
         nested("a1v54_nested.tsv", False)[cols].assign(set="dbNSFP 5.4a (A3)", cohort_DD="DD (31,058 trios)"),
         popeve()[cols].assign(set="popEVE (DDD only)", cohort_DD="DDD only")]
m4 = os.path.join(R, "m4nested_ddd_nested.tsv")
if os.path.exists(m4):
    parts.append(nested("m4nested_ddd_nested.tsv", False).query("score != 'AlphaMissense'")[cols]
                 .assign(set="this study: null-contrast model", cohort_DD="DDD only"))
m3 = os.path.join(R, "m4v3nested_ddd_nested.tsv")
if os.path.exists(m3):
    parts.append(nested("m4v3nested_ddd_nested.tsv", False).query("score != 'AlphaMissense'")[cols]
                 .assign(set="this study: null-contrast model (gene-matched, post hoc)", cohort_DD="DDD only"))
t = pd.concat(parts, ignore_index=True)
t = t.rename(columns={c: f"DD_{c}" for c in cols if c != "score"})

rep = nested(os.path.join("replication_affected", "rep_nested.tsv"), False)[["score", "A", "r_D2", "r_E2", "E2", "E2_lo_sys"]]
t = t.merge(rep.rename(columns={c: f"autism_{c}" for c in rep.columns if c != "score"}), on="score", how="left")

led = pd.read_csv(os.path.join(HERE, "aim3_training_ledger.tsv"), sep="\t")
lab = dict(zip(led["tool"], led.iloc[:, 4].astype(str).str.split(" ").str[0]))
t["clinical_training_labels"] = t["score"].map(lambda s: "none (this study)" if s.startswith("Aim4") else lab.get(s, "?"))

t["informative_DD"] = t["DD_A"] >= SUPP
t["E2_reaches_Supporting"] = t["DD_E2_lo_sys"] >= SUPP
# share of a predictor's ACMG points that is not explained by gene, position and constraint (points are ~ ln LR)
t["max_PP3_points_with_PM1"] = np.where(t["informative_DD"], np.round(4 * t["DD_r_E2"], 2), np.nan)
t["PP3_plus_PM1_max_points"] = np.round(2 + t["max_PP3_points_with_PM1"], 2)
t["retention_rank_E2"] = t["DD_r_E2"].where(t["informative_DD"]).rank(ascending=False, method="min")
t = t.sort_values(["informative_DD", "DD_r_E2"], ascending=[False, False])
out = os.path.join(R, "supp_table1_predictor_resource.tsv")
t.to_csv(out, sep="\t", index=False, float_format="%.4g")
show = ["score", "set", "DD_A", "DD_B", "DD_E2", "DD_r_B", "DD_r_D2", "DD_r_E2", "autism_r_D2", "clinical_training_labels",
        "max_PP3_points_with_PM1"]
pd.set_option("display.width", 250)
print(t[t["informative_DD"]][show].round(2).to_string(index=False))
print(f"\nwritten {out}: {len(t)} predictors, {int(t['informative_DD'].sum())} informative in DD")
