"""Supplementary Table 4: conditional PP3 calibration. Points PP3 earns by predictor x tier, given a PM1-like / PP2-like flag
and given its absence (2.406 scale; negative values floored at 0), from A10 (results_fullstudy/a10_overlap.tsv)."""
import math, sys
import pandas as pd

t = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else "results_fullstudy/a10_overlap.tsv", sep="\t")
t = t[t["flag"].isin(["PM1-like", "PP2-like"]) & t["pp3_level"].isin(["PP3 Supporting only", "PP3 Moderate only", "PP3 >=Strong(3-4)"])]
pts = lambda x: round(max(0.0, math.log(x) / math.log(2.406)), 2) if x > 0 else 0.0
t["points_given_flag"] = t["LR_PP3_given_F"].map(pts)
t["points_without_flag"] = t["LR_PP3_given_notF"].map(pts)
t["published_points"] = t["pp3_level"].map({"PP3 Supporting only": 1, "PP3 Moderate only": 2, "PP3 >=Strong(3-4)": 4})
cols = ["cohort", "tool", "flag", "pp3_level", "published_points", "points_given_flag", "points_without_flag", "R", "R_lo95", "R_hi95"]
t[cols].to_csv("supp_table4_conditional_pp3.tsv", sep="\t", index=False, float_format="%.3g")
dd = t[t["cohort"] == "DD"].pivot_table(index=["tool", "pp3_level"], columns="flag", values="points_given_flag")
print(dd.round(2).to_string())
