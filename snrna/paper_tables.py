#!/usr/bin/env python3
"""Main paper tables and the rarity forest plot, generated from committed result files only (no new estimation), so no
number in the manuscript is typed by hand. Missing inputs are listed and their table is skipped.

    python paper_tables.py        # -> paper/tables.md, paper/figures/fig10_rarity_forest.{png,pdf}

Table 1  Blind RNU2-2: sealed (blind_robustness) and literature-only (A8) AUROC, rarity null, excess.
Table 2  Rarity across benchmarks (A9 ClinVar non-coding, A10 ClinVar missense, A11 NCBoost, A13 TraitGym): raw AUROC
         and rarity-matched AUROC with window-bootstrap 95% CIs (A15), deattenuated rarity-matched AUROC (pi 0.05-0.20),
         and the paired tool differences whose intervals exclude 0.
Table 3  PP3 evidence: published thresholds against curated B/LB (A14) and against the rare-variant population (A16),
         with three-way verdicts and break-even pi.
"""
import glob
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
PAPER = os.path.join(HERE, "paper")
NAMES = {"a9": "ClinVar non-coding (A9)", "a10": "ClinVar missense (A10)", "a11": "NCBoost, published (A11)",
         "a13_traitgym_mendelian": "TraitGym Mendelian, published (A13)"}
TOOLS = {"cadd_phred": "CADD", "gpn_msa": "GPN-MSA", "remm": "ReMM", "phylop447": "phyloP447", "alphamissense": "AlphaMissense"}


def md_table(df):
    """Markdown table without the optional tabulate dependency."""
    f = lambda v: f"{v:.3f}" if isinstance(v, float) and not np.isnan(v) else ("" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v))
    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    return "\n".join(lines + ["| " + " | ".join(f(v) for v in row) + " |" for row in df.itertuples(index=False)])


def read(name):
    p = os.path.join(RES, name)
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None


def ci(r, k):
    return f"{r[k + '_point']:.3f} [{r[k + '_lo']:.3f}, {r[k + '_hi']:.3f}]"


def table_rarity():
    rows, diffs, forest = [], [], []
    for p in sorted(glob.glob(os.path.join(RES, "a15_ci_*.tsv"))):
        tag = os.path.basename(p)[7:-4]
        t = pd.read_csv(p, sep="\t")
        m = t[~t.metric.str.startswith("diff_")].pivot(index="tool", columns="metric", values=["point", "lo", "hi"])
        m.columns = [f"{b}_{a}" for a, b in m.columns]
        for tool, r in m.iterrows():
            M = r["rarity_matched_auroc_point"]
            rows.append({"Benchmark": NAMES.get(tag, tag), "Tool": TOOLS.get(tool, tool), "Raw AUROC": ci(r, "raw_auroc"),
                         "Rarity-only null": f"{r['null_mean_point']:.3f}", "Rarity-matched AUROC": ci(r, "rarity_matched_auroc"),
                         "Deattenuated (pi 0.05-0.20)": f"{(M - 0.025) / 0.95:.3f}-{(M - 0.10) / 0.80:.3f}",
                         "Share of above-chance AUROC from rarity": f"{(r['null_mean_point'] - 0.5) / (r['raw_auroc_point'] - 0.5):.0%}"
                         if r["raw_auroc_point"] > 0.52 else "n/a"})
            forest.append((NAMES.get(tag, tag), TOOLS.get(tool, tool), r))
        d = t[t.metric == "diff_rarity_matched_auroc"]
        for r in d.itertuples():
            if r.lo > 0 or r.hi < 0:
                a, b = r.tool.split(" - ")
                diffs.append(f"{NAMES.get(tag, tag)}: {TOOLS.get(a, a)} - {TOOLS.get(b, b)} {r.point:+.3f} [{r.lo:+.3f}, {r.hi:+.3f}]")
    return pd.DataFrame(rows), diffs, forest


def table_blind():
    out = []
    b = read("blind_robustness.tsv")
    if b is not None and {"score", "auroc"} <= set(b.columns):
        x = b[b.get("analysis", pd.Series("", index=b.index)).astype(str).str.contains("null|rarity", case=False, na=False)] if "analysis" in b else b
        out.append(("sealed (22 P/LP)", x))
    a8 = read("a8_results.tsv")
    if a8 is not None:
        out.append(("literature-only (A8)", a8[a8.set.str.contains("primary", na=False)]))
    return out


def table_calibration():
    parts = []
    a14 = read("a14_published.tsv")
    if a14 is not None:
        x = a14[a14.prior.round(4) == 0.0441]
        parts.append(pd.DataFrame({"Tool": x.tool, "Band": x.strength, "Threshold": x.published_threshold,
                                   "Comparison": "curated B/LB (" + x.set + ")", "pi": "-",
                                   "LR (5th pct, min above threshold)": x.min_lr_lo5_above.round(2), "LR needed": x.lr_needed.round(2),
                                   "Verdict": np.where(x.met_at_all_scores_above, "met", "not met (lower bound)"), "Break-even pi": "-"}))
    a16 = read("a16_bands.tsv")
    if a16 is not None:
        parts.append(pd.DataFrame({"Tool": a16.tool, "Band": a16.band, "Threshold": a16.lo, "Comparison": a16.comparison,
                                   "pi": a16.pi, "LR (5th pct, min above threshold)": [f"{l:.1f} [{a:.1f}, {b:.1f}]" for l, a, b in
                                                                                          zip(a16.lr, a16.lr_q5, a16.lr_q95)],
                                   "LR needed": a16.lr_needed.round(2), "Verdict": a16.verdict,
                                   "Break-even pi": a16.break_even_pi.round(3)}))
    return pd.concat(parts, ignore_index=True) if parts else None


def forest_plot(forest):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    if not forest:
        return
    fig, ax = plt.subplots(figsize=(7.5, 0.32 * len(forest) + 1.2))
    for i, (bench, tool, r) in enumerate(forest):
        y = len(forest) - i
        ax.errorbar(r["raw_auroc_point"], y, xerr=[[r["raw_auroc_point"] - r["raw_auroc_lo"]], [r["raw_auroc_hi"] - r["raw_auroc_point"]]],
                    fmt="o", color="#9aa0a6", ms=4, capsize=2)
        ax.errorbar(r["rarity_matched_auroc_point"], y, xerr=[[r["rarity_matched_auroc_point"] - r["rarity_matched_auroc_lo"]],
                                                              [r["rarity_matched_auroc_hi"] - r["rarity_matched_auroc_point"]]],
                    fmt="s", color="#2a78d6", ms=4, capsize=2)
        ax.plot([r["rarity_matched_auroc_point"], r["raw_auroc_point"]], [y, y], color="#c9ccd1", lw=0.8, zorder=0)
    ax.set_yticks(range(len(forest), 0, -1), [f"{b.split(' (')[0]} | {t}" for b, t, _ in forest], fontsize=7)
    ax.axvline(0.5, color="#52514e", lw=0.8, ls="--")
    ax.set_xlabel("AUROC (grey: raw, against common controls; blue: rarity-matched, against rare variants) with 95% CI")
    fig.tight_layout()
    os.makedirs(os.path.join(PAPER, "figures"), exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(PAPER, "figures", f"fig10_rarity_forest.{ext}"), dpi=200)


def main():
    md = ["# Main tables (generated by paper_tables.py from results/; do not edit by hand)", ""]
    blind = table_blind()
    md.append("## Table 1. Blind RNU2-2")
    if blind:
        for name, t in blind:
            md += [f"**{name}**", "", md_table(t), ""]
    else:
        md += ["_inputs missing: results/blind_robustness.tsv, results/a8_results.tsv (DGX; run sync_dgx_outputs.sh)_", ""]
    t2, diffs, forest = table_rarity()
    md.append("## Table 2. Rarity across benchmarks")
    if len(t2):
        md += [md_table(t2), "", "Paired differences in rarity-matched AUROC whose 95% CI excludes 0:", ""]
        md += [f"- {d}" for d in diffs] + [""]
    else:
        md += ["_inputs missing: results/a15_ci_*.tsv_", ""]
    t3 = table_calibration()
    md.append("## Table 3. PP3 evidence against curated benign variants and against the rare-variant population")
    md += [md_table(t3), ""] if t3 is not None else ["_inputs missing: results/a14_published.tsv, results/a16_bands.tsv_", ""]
    os.makedirs(PAPER, exist_ok=True)
    open(os.path.join(PAPER, "tables.md"), "w", encoding="utf-8").write("\n".join(md) + "\n")
    forest_plot(forest)
    print("\n".join(md))


if __name__ == "__main__":
    main()
