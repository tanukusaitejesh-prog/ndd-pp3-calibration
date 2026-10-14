#!/usr/bin/env python3
"""Post hoc, exploratory (after the single unsealing): would something trivial match the frozen v4 on blind RNU2-2?

Reviewer questions answered here, all with the frozen v4 scores (hash-checked; results_v4/ is read only):
  1. Position-only baselines. Every blind positive was classified with PM1 (critical region), so positives cluster.
     Baselines that use position alone: membership of the guidance Table 2 region (n.1-67, fixed before unsealing in
     table2_regions.tsv), distance to it, 5' proximity (-n), and the position's own per-site gnomAD observation
     (any SNV observed at the position, the crudest constraint). v4 must beat these to be more than "position".
  2. Within-region discrimination: AUROC of v4 (and comparators) among variants inside n.1-67 only, and outside it.
  3. The benign proxy: snRNA loci are hypermutable (~50x de novo rate; D'Souza 2026), so AC >= 9 is a weak benign
     definition. L3 is redefined at AC >= 20, 50 and 100 (PASS), and as AC >= 9 at positions of low Roulette rate.
  4. Paired bootstrap of v4 minus each baseline (2,000 replicates within class).

    python blind_robustness.py --blind sealed/blind_labels.tsv      # -> results/blind_robustness.{tsv,txt}
"""
import argparse
import hashlib
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evaluate as E  # noqa: E402
import snrna as S  # noqa: E402

GENE = "RNU2-2"


def paired(y, a, b, B, rng):
    """AUROC(a) - AUROC(b) with a 95% class-stratified bootstrap interval."""
    d = S.auroc(y, a) - S.auroc(y, b)
    i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
    bs = []
    for _ in range(B):
        idx = np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]
        bs.append(S.auroc(y[idx], a[idx]) - S.auroc(y[idx], b[idx]))
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return d, lo, hi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blind", required=True)
    ap.add_argument("--v4", default=os.path.join(HERE, "results_v4"))
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()
    sp = os.path.join(a.v4, "v4_scores.tsv.gz")
    man = {l.split()[1]: l.split()[0] for l in open(os.path.join(a.v4, "v4_manifest.txt")) if l.strip()}
    if man.get(os.path.relpath(sp, HERE)) != hashlib.sha256(open(sp, "rb").read()).hexdigest():
        sys.exit("v4_scores.tsv.gz does not match the frozen manifest")
    rng = np.random.default_rng(S.SEED)

    var = pd.read_csv(os.path.join(a.inputs, "variants.tsv.gz"), sep="\t")
    v = var[var.gene == GENE].merge(pd.read_csv(sp, sep="\t")[["ID", "v4", "local_constraint"]], on="ID", how="left")
    comp = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    v = v.merge(comp[["ID", "phylop447", "gpn_msa", "remm"]], on="ID", how="left")
    reg = pd.read_csv(os.path.join(HERE, "table2_regions.tsv"), sep="\t", comment="#")
    reg = reg[reg.gene == GENE]
    inreg = np.zeros(len(v), bool)
    dist = np.full(len(v), np.inf)
    for r in reg.itertuples():
        inreg |= v.n.between(r.start, r.end).to_numpy()
        dist = np.minimum(dist, np.where(v.n < r.start, r.start - v.n, np.where(v.n > r.end, v.n - r.end, 0)))
    v["in_table2_region"] = inreg.astype(float)
    v["minus_dist_region"] = -dist
    v["five_prime"] = -v.n.astype(float)
    obs_pos = v[v.AC_pass > 0].groupby("n").size()
    v["pos_unobserved"] = -v.n.map(obs_pos).fillna(0)  # fewer observed alleles at the position = more constrained
    v["rfam_alt_rarity"] = 1 - v.alt_col_freq
    v["rfam_ref_col_freq"] = v.ref_col_freq

    b = pd.read_csv(a.blind, sep="\t")
    b = b[b.ID.str.match(E.SNV_ID)]
    plp = set(b.ID[(b.gene == GENE) & (b.label == "P/LP")])
    pos = v.ID.isin(plp)

    scores = [("v4", "v4 (frozen)"), ("cadd_phred", "CADD"), ("phylop447", "phyloP447"), ("gpn_msa", "GPN-MSA"),
              ("remm", "ReMM"), ("rfam_alt_rarity", "Rfam alt rarity (label-free)"),
              ("rfam_ref_col_freq", "Rfam ref column freq (label-free)"), ("local_constraint", "v4 local constraint"),
              ("in_table2_region", "BASELINE in Table 2 region"), ("minus_dist_region", "BASELINE distance to region"),
              ("five_prime", "BASELINE 5' proximity"), ("pos_unobserved", "BASELINE position unobserved in gnomAD")]
    mr_lo = v.log_mr <= v.log_mr.quantile(0.5)
    negsets = {"L3 AC>=9 (frozen)": v.L3_benign_like.astype(bool),
               "L3 AC>=20": (v["filter"] == "PASS") & (v.AC_pass >= 20),
               "L3 AC>=50": (v["filter"] == "PASS") & (v.AC_pass >= 50),
               "L3 AC>=100": (v["filter"] == "PASS") & (v.AC_pass >= 100),
               "L3 AC>=9, low-mutability half": v.L3_benign_like.astype(bool) & mr_lo}
    subsets = {"all positions": np.ones(len(v), bool), "inside n.1-67": inreg, "outside n.1-67": ~inreg}

    rows = []
    for nname, neg in negsets.items():
        print(f"  negatives: {nname}", flush=True)
        for sname, keep in subsets.items():
            if sname != "all positions" and nname != "L3 AC>=9 (frozen)":
                continue
            m = keep & (pos | (neg.to_numpy() & ~pos.to_numpy()))
            for col, lab in scores:
                ok = m & v[col].notna().to_numpy()
                y = pos[ok].to_numpy()
                row = {"negatives": nname, "positions": sname, "score": lab, "n_P": int(y.sum()), "n_N": int((~y).sum())}
                if y.sum() >= 3 and (~y).sum() >= 3 and v.loc[ok, col].nunique() > 1:
                    x = v.loc[ok, col].to_numpy(float)
                    row["auroc"] = S.auroc(y, x)
                    row["lo"], row["hi"] = S.boot_auroc(y, x, a.boot)
                    if col != "v4":
                        okv = ok & v.v4.notna().to_numpy()
                        yv = pos[okv].to_numpy()
                        row["v4_minus"], row["d_lo"], row["d_hi"] = paired(yv, v.v4[okv].to_numpy(float),
                                                                           v.loc[okv, col].to_numpy(float), a.boot, rng)
                rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.results, "blind_robustness.tsv"), sep="\t", index=False, float_format="%.4f")
    rep = [f"Blind RNU2-2 robustness (post hoc, exploratory; frozen v4 hash verified). Positives: {int(pos.sum())} P/LP SNVs; "
           f"{int((pos & v.in_table2_region.astype(bool)).sum())} inside the Table 2 region n.1-67", ""]
    for (nname, sname), x in out.groupby(["negatives", "positions"], sort=False):
        rep.append(f"{nname} | {sname}: {x.n_P.iloc[0]} P vs {x.n_N.iloc[0]} N")
        for r in x.itertuples():
            au = f"AUROC {r.auroc:.3f} [{r.lo:.3f}, {r.hi:.3f}]" if pd.notna(getattr(r, "auroc", np.nan)) else "n/a"
            dv = f"   v4 - this {r.v4_minus:+.3f} [{r.d_lo:+.3f}, {r.d_hi:+.3f}]" if pd.notna(getattr(r, "v4_minus", np.nan)) else ""
            rep.append(f"  {r.score:40s} {au}{dv}")
        rep.append("")
    rep.append("Reading: v4 is more than position if 'v4 - this' excludes 0 for every BASELINE row, and it ranks within the "
               "region if its 'inside n.1-67' AUROC interval excludes 0.5.")
    open(os.path.join(a.results, "blind_robustness.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
