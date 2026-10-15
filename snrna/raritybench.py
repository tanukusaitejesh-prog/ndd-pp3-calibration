#!/usr/bin/env python3
"""RarityBench: one table of rarity-only inflation across every benchmark analysed (A9, A10, A11, A13). Exploratory
summary; reads committed results only and computes nothing new beyond ratios.

    python raritybench.py          # -> results/raritybench.tsv, results/raritybench_report.txt (+ figure if matplotlib)

For each benchmark x unit x tool: raw AUROC, rarity-only null, inflation (null - 0.5), corrected (raw - null), the
share of the above-chance AUROC that rarity alone explains (inflation / (raw - 0.5)), and the rarity-matched AUROC
where the source analysis computed one, deattenuated for an assumed damaging fraction pi of the rare pool (0.05-0.20).
A15's simulation showed corrected AUROC over-corrects, so the deattenuated rarity-matched AUROC is the magnitude to cite. Missing result files are listed and skipped.
"""
import glob
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
SOURCES = [("A9 ClinVar non-coding", "a9_results.tsv"), ("A10 ClinVar missense / SpliceAI", "a10_results.tsv"),
           ("A11 NCBoost (published)", "a11_results.tsv")]
PIS = (0.05, 0.10, 0.20)   # assumed damaging fraction of the rare pool (A15: deattenuated M = (M - pi/2)/(1 - pi))
COLS = ["benchmark", "unit", "tool", "n_P", "n_N", "raw_auroc", "null_mean", "inflation", "corrected", "share_from_rarity",
        "rarity_matched_auroc"] + [f"deattenuated_pi{int(p * 100):02d}" for p in PIS]


def collect(res=RES):
    parts, missing = [], []
    srcs = SOURCES + [(f"A13 {os.path.basename(p)[4:-12]} (published)", os.path.basename(p))
                      for p in sorted(glob.glob(os.path.join(res, "a13_*_results.tsv")))]
    for name, f in srcs:
        p = os.path.join(res, f)
        if not os.path.exists(p):
            missing.append(f); continue
        t = pd.read_csv(p, sep="\t")
        t["benchmark"] = name
        if "unit" not in t:
            t["unit"] = "all"
        for c in ("rarity_matched_auroc",):
            if c not in t:
                t[c] = np.nan
        parts.append(t)
    if not parts:
        return pd.DataFrame(columns=COLS), missing
    t = pd.concat(parts, ignore_index=True)
    above = t.raw_auroc - 0.5
    t["share_from_rarity"] = np.where(above > 0.02, t.inflation / above, np.nan)
    for p in PIS:
        t[f"deattenuated_pi{int(p * 100):02d}"] = (t.rarity_matched_auroc - p / 2) / (1 - p)
    return t[COLS], missing


def main():
    t, missing = collect()
    os.makedirs(RES, exist_ok=True)
    t.to_csv(os.path.join(RES, "raritybench.tsv"), sep="\t", index=False, float_format="%.4f")
    rep = ["RarityBench: rarity-only inflation across benchmarks (exploratory summary of A9-A13)", ""]
    if missing:
        rep.append("missing (skipped): " + ", ".join(missing))
    for b, x in t.groupby("benchmark", sort=False):
        rep.append(f"\n{b}")
        for r in x.itertuples():
            rep.append(f"  {str(r.unit):14s} {r.tool:14s} raw {r.raw_auroc:.3f}  null {r.null_mean:.3f}  inflation {r.inflation:+.3f}"
                       f"  rarity share {r.share_from_rarity:.0%}  rarity-matched {r.rarity_matched_auroc:.3f} (deattenuated {r.deattenuated_pi05:.3f}-{r.deattenuated_pi20:.3f})" if np.isfinite(r.share_from_rarity) else
                       f"  {str(r.unit):14s} {r.tool:14s} raw {r.raw_auroc:.3f}  null {r.null_mean:.3f}  inflation {r.inflation:+.3f}")
    if len(t):
        s = t.groupby("tool").share_from_rarity.median().sort_values(ascending=False)
        rep += ["", "Median share of above-chance AUROC explained by rarity, per tool: "
                + "; ".join(f"{k} {v:.0%}" for k, v in s.dropna().items())]
    open(os.path.join(RES, "raritybench_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        x = t[t.unit.isin(["pooled", "all"]) | t.benchmark.str.startswith("A1")].dropna(subset=["inflation"])
        if len(x):
            fig, ax = plt.subplots(figsize=(7, 0.35 * len(x) + 1))
            lab = x.benchmark.str.split(" ").str[0] + " " + x.unit.astype(str) + " | " + x.tool
            ax.barh(range(len(x)), x.inflation, color="#2a78d6")
            ax.set_yticks(range(len(x)), lab, fontsize=7)
            ax.axvline(0.05, ls="--", c="#52514e", lw=0.8)
            ax.set_xlabel("rarity-only inflation (null AUROC - 0.5)")
            fig.tight_layout()
            os.makedirs(os.path.join(HERE, "paper", "figures"), exist_ok=True)
            fig.savefig(os.path.join(HERE, "paper", "figures", "fig9_raritybench.png"), dpi=200)
    except ImportError:
        pass


if __name__ == "__main__":
    main()
