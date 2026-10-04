#!/usr/bin/env python3
"""Compare gate.py output on synthetic data with the known truth. Exit code 0 = every check passed."""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

out = sys.argv[1]
biased = "--biased" in sys.argv
cpg = "--cpg" in sys.argv
W, R, T = (os.path.join(out, x) for x in ("work", "results", "truth"))
EDGES = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 98, 99, 100]
fails = []


def check(name, cond, detail=""):
    print(f"  [{'ok' if cond else 'FAIL'}] {name}  {detail}")
    if not cond:
        fails.append(name)


meta = json.load(open(os.path.join(T, "meta.json")))
tr = pd.read_csv(os.path.join(T, "dnv_truth.tsv"), sep="\t", dtype={"chrom": str})
lt = pd.read_csv(os.path.join(T, "lr_truth.tsv"), sep="\t")
print(f"\nSelf-test checks ({'biased detection scenario' if biased else ('CpG-miscalibrated rate scenario' if cpg else 'clean scenario')})")

# 1. liftover
d = pd.read_csv(os.path.join(W, "dnv_grch38.tsv"), sep="\t", dtype={"chrom": str})
m = d.merge(tr, on="dnv_idx", suffixes=("", "_t"))
real = m[m["origin"] != "junk"]
ok = (real["lift"].eq("ok") & real["pos38"].eq(real["pos38_t"]) & real["ref38"].eq(real["ref38_t"])
      & real["alt38"].eq(real["alt38_t"]))
check("liftover recovers GRCh38 position and alleles (incl. reverse-strand chain)", bool(ok.all()),
      f"({int(ok.sum())}/{len(real)})")
junk = sorted(set(m.loc[m["origin"] == "junk", "lift"]))
check("unmappable / indel / chrX rows are flagged", junk == ["not_autosomal", "not_snv", "unmapped"], str(junk))

# 2. transcript statuses
tx = pd.read_pickle(os.path.join(W, "gtf", "transcripts.pkl"))
txq = pd.concat([pd.read_pickle(f) for f in glob.glob(os.path.join(W, "sites", "chr*.txqc.pkl"))])
st = tx.merge(txq[["tx", "reason"]], on="tx", how="left")
st["final"] = np.where(st["status"] != "ok", st["status"], st["reason"])
expect = {"normal": "pass", "AMGAP": "pass", "RGAP": "pass", "DUPA": "pass", "DUPB": "pass", "ISTOP": "internal_stop",
          "NF": "cds_incomplete_tag", "NOGTF": "not_in_gtf", "SHIFT": "low_concordance"}
st["expected"] = st["transcript_id"].map(meta["expected_status"]).map(expect)
bad = st[st["final"] != st["expected"]]
check("transcript QC catches internal stop, NF tag, missing GTF entry, shifted AlphaMissense",
      bad.empty, "" if bad.empty else "\n" + bad[["transcript_id", "final", "expected"]].to_string())

# 3. mutational opportunity (exact)
aggs = [pd.read_pickle(f) for f in sorted(glob.glob(os.path.join(W, "sites", "chr*.agg.pkl")))]
mis = pd.concat([a["mis_tx_am"] for a in aggs])
syn = pd.concat([a["syn_tx"] for a in aggs])
mis, syn = mis[mis["filt_ok"].astype(bool)], syn[syn["filt_ok"].astype(bool)]
cb = np.searchsorted(EDGES, mis["am_bin"].to_numpy(), side="right") - 1
Mb = np.bincount(cb, weights=mis["MR"].to_numpy(), minlength=len(EDGES) - 1)
check("missense mutational mass in every AlphaMissense bin equals the reference enumeration",
      bool(np.allclose(Mb, lt["M_mis"], rtol=1e-9, atol=0)), f"(max rel diff {np.max(np.abs(Mb / lt['M_mis'] - 1)):.2e})")
check("synonymous mutational mass equals the reference enumeration", bool(np.isclose(syn["MR"].sum(), meta["M_syn"], rtol=1e-9)),
      f"({syn['MR'].sum():.6g} vs {meta['M_syn']:.6g})")

# 4. observed counts
lr = pd.read_csv(os.path.join(R, "lr_bins.tsv"), sep="\t")
allm = lr[(lr["analysis"] == "all_genes") & (lr["kind"] == "marginal")].reset_index(drop=True)
tO = np.bincount(tr.loc[tr["universe_cls"] == "mis", "bin"].astype(int), minlength=len(EDGES) - 1)
check("observed de novo missense counts per bin equal the simulated truth", bool(np.array_equal(allm["O"].round().astype(int), tO)),
      f"(pipeline {int(allm['O'].sum())}, truth {int(tO.sum())})")
nc = pd.read_csv(os.path.join(R, "null_checks.tsv"), sep="\t")
osyn = nc.loc[nc["grouping"] == "gene s_het class", "O"].sum()
tsyn = int((tr["universe_cls"] == "syn").sum())
check("observed synonymous de novo count equals the simulated truth", int(round(osyn)) == tsyn, f"({osyn:.0f} vs {tsyn})")

# 5. statistical recovery of the true likelihood ratios
rep = open(os.path.join(R, "gate_report.txt")).read()
verdict = [l for l in rep.splitlines() if "OVERALL:" in l][0].strip()
nv = nc.groupby("grouping")["verdict"].first()
gated = ["Roulette mutability decile", "AlphaMissense context of the codon", "mutation class"]
if not biased:
    # (a) the excess must track the de novos that were simulated as pathogenic, up to background noise
    P_real = np.bincount(tr.loc[(tr["origin"] == "pathogenic") & (tr["universe_cls"] == "mis"), "bin"].astype(int),
                         minlength=len(EDGES) - 1)
    for an in ("all_genes", "ddg2p_monoallelic"):
        x = lr[(lr["analysis"] == an) & (lr["kind"] == "marginal")].reset_index(drop=True)
        tol = 4 * np.sqrt(x["E"] + (0.02 * x["E"]) ** 2) + 3
        dev = np.abs(x["excess"] - P_real)
        check(f"excess de novos track the simulated pathogenic de novos in every bin ({an})", bool((dev <= tol).all()),
              f"(max |excess - pathogenic| / tolerance = {np.max(dev / tol):.2f})")
    # (b) the 95% CIs cover the expected (population) likelihood ratio at about the nominal rate
    for an, kind, col in (("all_genes", "marginal", "LR_marginal"), ("all_genes", "within_gene", "LR_within"),
                          ("ddg2p_monoallelic", "marginal", "LR_ddg2p_marginal")):
        x = lr[(lr["analysis"] == an) & (lr["kind"] == kind)].reset_index(drop=True)
        cover = (x["lo95"] <= lt[col]) & (lt[col] <= x["hi95"])
        print(f"\n    {an} / {kind}:   bin       truth   estimate   95% CI")
        for i in range(len(x)):
            print(f"      {x['bin_lo'][i]:.2f}-{x['bin_hi'][i]:.2f}  {lt[col][i]:8.3f}  {x['LR'][i]:8.3f}   "
                  f"[{x['lo95'][i]:.3f}, {x['hi95'][i]:.3f}] {'' if cover[i] else '  <- outside'}")
        se = (x["hi95"] - x["lo95"]) / 3.92
        z = ((x["LR"] - lt[col]) / se).replace([np.inf, -np.inf], np.nan)
        check(f"true LR recovered ({an}, {kind}): >= 9/13 bins inside the 95% CI, |mean z| < 1, max |z| < 4.5",
              int(cover.sum()) >= len(x) - 4 and abs(np.nanmean(z)) < 1.0 and np.nanmax(np.abs(z)) < 4.5,
              f"({int(cover.sum())}/{len(x)} inside; mean z {np.nanmean(z):+.2f}; max |z| {np.nanmax(np.abs(z)):.2f})")
    if cpg:
        check("pooled mutation-class null check FAILS when CpG rates are miscalibrated", nv.get("mutation class") == "FAIL",
              str(nv.to_dict()))
        check("pre-registered class-specific scaling is triggered and accepted", "fallback accepted" in rep)
        cg = [g + " [class-calibrated]" for g in gated[:2]]
        check("mutability and context null checks pass after class-specific scaling", all(nv.get(g) in ("PASS", "WARN") for g in cg),
              str({g: nv.get(g) for g in cg}))
        check("overall verdict is GO WITH CAUTION", "GO WITH CAUTION" in verdict, verdict)
    else:
        check("no null check FAILS when detection is unbiased (WARN tolerated: it is a 10%-band flag)",
              all(nv.get(g) in ("PASS", "WARN") for g in gated), str({g: nv.get(g) for g in gated}))
        check("overall verdict is GO (or GO WITH CAUTION after a chance WARN)", verdict.endswith("GO") or "GO WITH CAUTION" in verdict,
              verdict)
else:
    check("context null check FAILS when detection depends on protein context",
          nv.get("AlphaMissense context of the codon") == "FAIL", str(nv.to_dict()))
    check("overall verdict is KILL", "KILL" in verdict, verdict)

print("\nSELF-TEST " + ("PASSED" if not fails else f"FAILED: {fails}"))
sys.exit(1 if fails else 0)
