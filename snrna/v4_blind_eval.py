#!/usr/bin/env python3
"""Blind evaluation of the frozen v4 scores (amendment A5). Called by final_eval.py after the single unsealing.
Checks results_v4/v4_scores.tsv.gz against results_v4/v4_manifest.txt (frozen before unsealing), then for RNU2-2 and
RNU5B-1: AUROC of L2 P/LP vs L3 SNVs for v4 and CADD (2000-bootstrap CI), delta, and the weighted delta over both.

    python v4_blind_eval.py --blind sealed/blind_labels.tsv
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blind", required=True)
    ap.add_argument("--v4", default=os.path.join(HERE, "results_v4"))
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()
    sp = os.path.join(a.v4, "v4_scores.tsv.gz")
    man = {l.split()[1]: l.split()[0] for l in open(os.path.join(a.v4, "v4_manifest.txt")) if l.strip()}
    key = os.path.relpath(sp, HERE)
    if man.get(key) != hashlib.sha256(open(sp, "rb").read()).hexdigest():
        sys.exit("v4_scores.tsv.gz does not match the frozen manifest; not evaluating")
    sc = pd.read_csv(sp, sep="\t")
    var = pd.read_csv(os.path.join(a.inputs, "variants.tsv.gz"), sep="\t")[["ID", "cadd_phred", "L3_benign_like"]]
    sc = sc.drop(columns=[c for c in ("cadd_phred",) if c in sc.columns]).merge(var, on="ID", how="left")
    b = pd.read_csv(a.blind, sep="\t")
    b = b[b.ID.str.match(E.SNV_ID)]
    plp = set(b.ID[b.label == "P/LP"])
    rng = np.random.default_rng(S.SEED)
    out, sets = [], []
    for gene in S.SEALED_GENES:
        g = sc[sc.gene == gene]
        pos = g.ID.isin(plp); neg = g.L3_benign_like.astype(bool) & ~pos
        ok = (pos | neg) & g.v4.notna() & g.cadd_phred.notna()
        y = pos[ok].to_numpy(); s4 = g.v4[ok].to_numpy(float); c = g.cadd_phred[ok].to_numpy(float)
        row = {"gene": gene, "n_L2": int(y.sum()), "n_L3": int((~y).sum())}
        if y.sum() and (~y).sum():
            for name, x in (("v4", s4), ("cadd", c)):
                row[f"auroc_{name}"] = S.auroc(y, x)
                row[f"{name}_lo"], row[f"{name}_hi"] = S.boot_auroc(y, x, a.boot)
            row["delta"] = row["auroc_v4"] - row["auroc_cadd"]
            sets.append((y, s4, c))
        out.append(row)
    lines = ["v4 blind evaluation (A5; frozen scores, hash verified)", pd.DataFrame(out).to_string(index=False, float_format=lambda x: f"{x:.3f}")]
    if sets:
        w = np.array([y.sum() for y, _, _ in sets], float)
        pt = (w * np.array([S.auroc(y, s) - S.auroc(y, c) for y, s, c in sets])).sum() / w.sum()
        bs = []
        for _ in range(a.boot):
            vals = []
            for y, s, c in sets:
                i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
                idx = np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]
                vals.append(S.auroc(y[idx], s[idx]) - S.auroc(y[idx], c[idx]))
            bs.append((w * np.array(vals)).sum() / w.sum())
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        lines.append(f"Weighted delta AUROC v4 - CADD, blind genes: {pt:+.3f} [{lo:+.3f}, {hi:+.3f}]")
    pd.DataFrame(out).to_csv(os.path.join(a.results, "final_v4_blind.tsv"), sep="\t", index=False, float_format="%.4f")
    open(os.path.join(a.results, "final_v4_blind.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
