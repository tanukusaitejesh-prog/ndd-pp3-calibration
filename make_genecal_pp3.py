"""A19 H3 (pre-registered, secondary): calibrate, don't cap.
AlphaMissense PP3 (Supporting) thresholds recalibrated per gene: the lowest score threshold t at which, within the
gene, LR(score >= t) = P(score >= t | P/LP) / P(score >= t | B/LB) reaches 2.406, using the ClinVar classes of the
missense resource (supp_data1) and leaving out every variant in the laboratory truth set. Genes with < 10 P/LP or
< 10 B/LB use the published threshold. On laboratory truth (2026-10): R(flag, PM1) for published vs gene-calibrated
flags, and their ratio (gene-clustered bootstrap)."""
import math

import numpy as np
import pandas as pd

rng = np.random.default_rng(20261008)
r = pd.read_csv("results_fullstudy/supp_data1_clinvar_pp3_pm1_points.tsv.gz", sep="\t",
                usecols=["clinvar_id", "clinvar_class", "AlphaMissense_score", "AlphaMissense_PP3_published"])
vs = pd.read_csv("work/cvt/variant_summary_2026-10.txt.gz", sep="\t", dtype=str, usecols=["VariationID", "GeneSymbol", "Assembly"])
vs = vs[vs.Assembly == "GRCh38"].drop_duplicates("VariationID")
r["VariationID"] = r.clinvar_id.astype(str)
r = r.merge(vs[["VariationID", "GeneSymbol"]], on="VariationID", how="left").dropna(subset=["AlphaMissense_score"])
lab = pd.read_csv("work/cvt/lab_codes_2026-10.tsv.gz", sep="\t", usecols=["VariationID", "gene", "PM1", "indep", "Submitter"])
lab = lab[(lab.indep >= 6) | (lab.indep <= -1)].copy()
lab["y"] = (lab.indep >= 6).astype(int)
lab["VariationID"] = lab.VariationID.astype(str)
held = set(lab.VariationID)
pub_t = r.loc[r.AlphaMissense_PP3_published >= 1, "AlphaMissense_score"].min()
cal = r[r.clinvar_class.isin(["P/LP", "B/LB"]) & ~r.VariationID.isin(held)]
thr = {}
for g, d in cal.groupby("GeneSymbol"):
    p, b = d[d.clinvar_class == "P/LP"].AlphaMissense_score.to_numpy(), d[d.clinvar_class == "B/LB"].AlphaMissense_score.to_numpy()
    if len(p) < 10 or len(b) < 10:
        continue
    for t in np.unique(np.concatenate([p, b])):
        lr = ((p >= t).mean() + 1e-9) / (((b >= t).sum() + 0.5) / (len(b) + 1))
        if lr >= 2.406 and (p >= t).sum() >= 3:
            thr[g] = t
            break
lab = lab.merge(r[["VariationID", "AlphaMissense_score"]].drop_duplicates("VariationID"), on="VariationID")
lab["t_gene"] = lab.gene.map(thr).fillna(pub_t)
lab["PP3pub_on"] = (lab.AlphaMissense_score >= pub_t).astype(int)
lab["PP3gene_on"] = (lab.AlphaMissense_score >= lab.t_gene).astype(int)
lab["PM1_on"] = (lab.PM1 > 0).astype(int)


def lr(d, m):
    p, b = d[d.y == 1], d[d.y == 0]
    return ((m(p).sum() + 0.5) / (len(p) + 1)) / ((m(b).sum() + 0.5) / (len(b) + 1))


def R(d, f):
    return lr(d, lambda x: (x[f] == 1) & (x.PM1_on == 1)) / (lr(d, lambda x: x[f] == 1) * lr(d, lambda x: x.PM1_on == 1))


def ci(d, fn, n=1000):
    g = d.groupby("gene").indices
    k = list(g)
    v = [fn(d.iloc[np.concatenate([g[k[i]] for i in rng.choice(len(k), len(k))])]) for _ in range(n)]
    return np.nanquantile(v, [0.025, 0.975])


out = [f"A19 H3: published AlphaMissense Supporting threshold {pub_t:.3f}; gene-specific thresholds for {len(thr)} genes "
       f"(median {np.median(list(thr.values())):.3f}, IQR {np.quantile(list(thr.values()), .25):.3f}-{np.quantile(list(thr.values()), .75):.3f}); "
       f"laboratory truth SCVs with AlphaMissense {len(lab)} (in calibrated genes {int(lab.gene.isin(thr).sum())})"]
for f, nm in [("PP3pub_on", "published"), ("PP3gene_on", "gene-calibrated")]:
    v, (lo, hi) = R(lab, f), ci(lab, lambda b, f=f: R(b, f))
    out.append(f"  R(AlphaMissense PP3 {nm:15s}, PM1) {v:.2f} [{lo:.2f}, {hi:.2f}]; LR(PP3) {lr(lab, lambda x, f=f: x[f] == 1):.2f}")
rat, (lo, hi) = R(lab, "PP3gene_on") / R(lab, "PP3pub_on"), ci(lab, lambda b: R(b, "PP3gene_on") / R(b, "PP3pub_on"))
out.append(f"  ratio gene-calibrated / published R {rat:.2f} [{lo:.2f}, {hi:.2f}]")
sub = lab[lab.gene.isin(thr)]
for f, nm in [("PP3pub_on", "published"), ("PP3gene_on", "gene-calibrated")]:
    v, (lo, hi) = R(sub, f), ci(sub, lambda b, f=f: R(b, f))
    out.append(f"  calibrated genes only: R({nm:15s}, PM1) {v:.2f} [{lo:.2f}, {hi:.2f}]")
out.append("Interpretation (pre-specified): interval including 1 -> 'calibrate, don't cap'; R stays < 1 -> dependence not an artefact of global thresholds")
open("results_fullstudy/a19_h3_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
