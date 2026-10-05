"""Supplementary Data 1: published vs conditional PP3 + PM1 points for every matched ClinVar missense variant
(AlphaMissense, REVEL, VEST4), from the committed A4.1 match and the A1 retention shares. No new estimation."""
import os
import sys

import numpy as np
import pandas as pd

R = sys.argv[1] if len(sys.argv) > 1 else "/root/ndd_pp3_gate/results_fullstudy"
COL = "/root/ndd_pp3_gate/work/fullstudy/col"
TOOLS = {"AlphaMissense": ("am", [0.792, 0.906, 0.972, 0.990]),
         "REVEL": ("raw_REVEL_score", [0.644, 0.773, 0.879, 0.932]),
         "VEST4": ("raw_VEST4_score", [0.764, 0.861, 0.909, 0.965])}
x = pd.read_csv(os.path.join(R, "a4_clinvar_matched.tsv.gz"), sep="\t", dtype={"chrom": str})
nest = pd.read_csv(os.path.join(R, "a1_nested.tsv"), sep="\t")
r_e2 = nest[nest["stratum"] == "all_genes"].set_index("score")["r_E2"]
parts = []
for c, g in x.groupby("chrom", sort=False):
    key = np.load(os.path.join(COL, f"chr{c}", "key.npy"))
    idx = pd.Series(np.arange(len(key)), index=key).reindex(g["key"].to_numpy()).to_numpy()
    g = g.copy()
    for t, (col, _) in TOOLS.items():
        if col != "am":
            g[col] = np.load(os.path.join(COL, f"chr{c}", col + ".npy"), mmap_mode="r")[idx.astype(np.int64)]
    parts.append(g)
x = pd.concat(parts, ignore_index=True)
out = pd.DataFrame({"clinvar_id": x["clinvar_id"], "chrom": x["chrom"], "pos_grch38": x["key"] // 4,
                    "alt": np.array(list("ACGT"))[x["key"] % 4], "clinvar_class": x["cls"],
                    "pm1_like": x["pm1like"].astype(bool)})
pm1 = np.where(out["pm1_like"], 2.0, 0.0)
for t, (col, cuts) in TOOLS.items():
    s = pd.to_numeric(x[col], errors="coerce").to_numpy(float)
    pub = np.select([s >= cuts[3], s >= cuts[2], s >= cuts[1], s >= cuts[0]], [4, 3, 2, 1], 0).astype(float)
    pub[np.isnan(s)] = np.nan
    out[f"{t}_score"] = s
    out[f"{t}_PP3_published"] = pub
    out[f"{t}_PP3_PM1_capped"] = np.where(out["pm1_like"], np.minimum(pub + pm1, 4.0), pub + pm1)
    out[f"{t}_PP3_PM1_conditional"] = np.round(pub * float(r_e2[t]) + pm1, 3)
dst = os.path.join(R, "supp_data1_clinvar_pp3_pm1_points.tsv.gz")
out.to_csv(dst, sep="\t", index=False, float_format="%.4g")
print(f"written {dst}: {len(out):,} variants; r_E2 used: " + ", ".join(f"{t} {float(r_e2[t]):.3f}" for t in TOOLS))
print(out.groupby("clinvar_class")["pm1_like"].agg(["size", "mean"]).round(3))
