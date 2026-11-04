#!/usr/bin/env python3
"""Exploratory, post hoc, written 5 Nov 2026 before running: do the v4 PP3 thresholds hold against random RARE variants?

The v4 calibration (v4_calibrate.py) used L3 negatives (gnomAD v4.1 PASS AC >= 9; median AF 1.8e-4). Because scores track
population rarity, a likelihood ratio against population-observed variants can overstate evidence for a rare variant.
This check uses the same out-of-sample v4 scores and the frozen PP3/BP4 thresholds at prior 0.10
(results/v4_calibration_thresholds.tsv). For each unit (pooled; RNU2-2; RNU4-2; RNU4ATAC) it reports:
  - pass_P   : fraction of P/LP SNVs at or above the PP3 Supporting threshold
  - pass_L3  : the same for L3
  - pass_pool: the same for the rare pool (SNVs in the unit's loci that are neither P/LP nor L3, i.e. AC < 9 or absent;
               never labelled, so they include an unknown fraction pi of damaging variants)
  - LR_pool = pass_P / pass_pool. Because the pool contains damaging variants, LR_pool is a LOWER bound on the
    likelihood ratio of passing the threshold for P/LP versus rare benign variants (positive-unlabelled reasoning).
Decision rule fixed in advance: PP3 Supporting is "robust to rarity" in a unit if LR_pool >= 350**(1/8) = 2.08 (the
Supporting LR at prior 0.10); otherwise "not shown". BP4 is mirrored: frac_P / frac_pool at or below the BP4 threshold
must be <= 1/2.08 for "robust". The lower 95% bound of LR_pool (2,000 class-stratified bootstraps) is also reported.

    python v4_rare_passrate.py --blind sealed/blind_labels.tsv   # DGX; needs results_v4/ and results/v4_calibration_thresholds.tsv
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
import v4_calibrate as VC  # noqa: E402

SUP = 350 ** (1 / 8)


def lr_with_ci(s_pos, s_pool, thr, direction, B, rng):
    f = (lambda s: s >= thr) if direction == "PP3" else (lambda s: s <= thr)
    fp, fq = f(s_pos).mean(), f(s_pool).mean()
    lr = fp / fq if fq > 0 else np.inf
    bs = []
    for _ in range(B):
        a = f(rng.choice(s_pos, len(s_pos))).mean()
        b = f(rng.choice(s_pool, len(s_pool))).mean()
        bs.append(a / b if b > 0 else np.inf)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return fp, fq, lr, lo, hi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blind", required=True)
    ap.add_argument("--v4", default=os.path.join(HERE, "results_v4"))
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--labels", default=os.path.join(HERE, "labels"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()
    sp = os.path.join(a.v4, "v4_scores.tsv.gz")
    man = {l.split()[1]: l.split()[0] for l in open(os.path.join(a.v4, "v4_manifest.txt")) if l.strip()}
    if man.get(os.path.relpath(sp, HERE)) != hashlib.sha256(open(sp, "rb").read()).hexdigest():
        sys.exit("v4_scores.tsv.gz does not match the frozen manifest")
    th = pd.read_csv(os.path.join(a.results, "v4_calibration_thresholds.tsv"), sep="\t")
    th = th[(th.score == "v4") & (np.isclose(th.prior, 0.10))].set_index("unit")
    frozen = pd.read_csv(sp, sep="\t")[["ID", "v4"]]
    cf = VC.crossfit_v4(a.inputs, a.evo2, a.labels)  # identical recipe to v4_calibrate.py
    sc = cf.merge(frozen, on="ID", how="left")
    sc["s"] = np.where(sc.gene.isin(S.SEALED_GENES), sc.v4, sc.v4_cf)
    var = pd.read_csv(os.path.join(a.inputs, "variants.tsv.gz"), sep="\t")[["ID", "L3_benign_like"]]
    sc = sc.merge(var, on="ID", how="left")
    L = E.load_l2(a.labels)
    b = pd.read_csv(a.blind, sep="\t")
    b = b[b.ID.str.match(E.SNV_ID)]
    plp = set(L.ID[L.label == "P/LP"]) | set(b.ID[b.label == "P/LP"])
    sc["pos"] = sc.ID.isin(plp)
    sc["l3"] = sc.L3_benign_like.astype(bool) & ~sc.pos
    sc["pool"] = ~sc.pos & ~sc.L3_benign_like.astype(bool)
    sc = sc.dropna(subset=["s"])
    rng = np.random.default_rng(20261105)
    units = {"pooled (all genes)": sc[sc.gene.isin(set(sc.gene[sc.pos]))], "RNU2-2": sc[sc.gene == "RNU2-2"],
             "RNU4-2": sc[sc.gene == "RNU4-2"], "RNU4ATAC": sc[sc.gene == "RNU4ATAC"]}
    rows = []
    for u, x in units.items():
        if u not in th.index:
            continue
        for direction, col in (("PP3", "PP3_Supporting"), ("BP4", "BP4_Supporting")):
            thr = th.loc[u, col]
            if pd.isna(thr):
                continue
            fp, fq, lr, lo, hi = lr_with_ci(x.s[x.pos].to_numpy(), x.s[x.pool].to_numpy(), thr, direction, a.boot, rng)
            fl3 = ((x.s[x.l3] >= thr) if direction == "PP3" else (x.s[x.l3] <= thr)).mean()
            robust = (lr >= SUP) if direction == "PP3" else (lr <= 1 / SUP)
            rows.append({"unit": u, "evidence": f"{direction} Supporting", "threshold": thr, "n_P": int(x.pos.sum()),
                         "n_pool": int(x.pool.sum()), "n_L3": int(x.l3.sum()), "frac_P": fp, "frac_pool": fq, "frac_L3": fl3,
                         "LR_pool": lr, "LR_pool_lo": lo, "LR_pool_hi": hi,
                         "verdict": "robust to rarity" if robust else "not shown"})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.results, "v4_rare_passrate.tsv"), sep="\t", index=False, float_format="%.4g")
    print(f"Supporting LR at prior 0.10 = {SUP:.3f}; LR_pool is a lower bound on LR(P vs rare benign)")
    print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
