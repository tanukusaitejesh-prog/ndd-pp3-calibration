"""Within-gene decomposition of the points over-counted by summing PP3 with PM1-like / PP2-like evidence (draft 22).
Points over-counted = -ln R / ln 2.406, from the gene-clustered estimates in results_fullstudy/a12_gclust2.tsv:
all sites (R), within DDG2P monoallelic genes (R_within_ddg2p) and within genes across all genes (R_within_all).
Interval endpoints are transformed from the R intervals. Because the ACMG/AMP points framework applies one prior to every
variant, between-gene information enters only through the codes; overlap in it is therefore counted twice as well.
Outputs: results_fullstudy/within_gene_cost.tsv, results_fullstudy/within_gene_cost_summary.txt.
Run: python make_within_gene_cost.py"""
import math, os, re
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results_fullstudy")
d = pd.read_csv(os.path.join(RES, "a12_gclust2.tsv"), sep="\t")
LN = math.log(2.406)


def pts(r):
    return -math.log(r) / LN


def ci(s):
    lo, hi = [float(x) for x in re.findall(r"-?[\d.]+", s)]
    return pts(hi), pts(lo)


rows = []
for r in d.itertuples():
    rec = dict(cohort=r.cohort, predictor=r.tool, flag=r.flag, pp3=r.tier)
    for k, v, c in (("all", r.R, r.R_ci), ("ddg2p", r.R_within_ddg2p, r.R_within_ddg2p_ci), ("within", r.R_within_all, r.R_within_all_ci)):
        lo, hi = ci(c)
        rec.update({f"R_{k}": v, f"pts_{k}": pts(v), f"pts_{k}_lo": lo, f"pts_{k}_hi": hi})
    rec["share_ddg2p"] = rec["pts_ddg2p"] / rec["pts_all"] if rec["pts_all"] > 0.05 else float("nan")
    rows.append(rec)
T = pd.DataFrame(rows)
T.to_csv(os.path.join(RES, "within_gene_cost.tsv"), sep="\t", index=False, float_format="%.3f")

out = ["Points over-counted (-ln R / ln 2.406): all sites | within DDG2P genes | within genes (all genes); 95% intervals"]
for r in T.itertuples():
    out.append(f"  {r.cohort:7s} {r.predictor:13s} {r.flag:10s} PP3 {r.pp3:10s}: {r.pts_all:.2f} [{r.pts_all_lo:.2f}, {r.pts_all_hi:.2f}] | "
               f"{r.pts_ddg2p:.2f} [{r.pts_ddg2p_lo:.2f}, {r.pts_ddg2p_hi:.2f}] | {r.pts_within:.2f} [{r.pts_within_lo:.2f}, {r.pts_within_hi:.2f}]")
am = T[(T.predictor == "AlphaMissense") & T.flag.str.startswith("PM1")]
out.append(f"AlphaMissense x PM1 flags, within DDG2P genes: lower 95% bound > 0 in {(am.pts_ddg2p_lo > 0).sum()}/{len(am)} cohort-flag-strength cells; "
           f"share of the total over-count surviving within DDG2P genes in DD: {am[am.cohort == 'DD'].share_ddg2p.min():.0%}-{am[am.cohort == 'DD'].share_ddg2p.max():.0%}")
pp2 = T[T.flag == "PP2-like"]
out.append(f"PP3 x PP2-like within genes (all genes): points {pp2.pts_within.min():.2f} to {pp2.pts_within.max():.2f}; PP2 is constant within a gene, so its overlap is between-gene")
open(os.path.join(RES, "within_gene_cost_summary.txt"), "w", encoding="utf8").write("\n".join(out) + "\n")
print("\n".join(out[-2:]))
