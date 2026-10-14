#!/usr/bin/env python3
"""Post hoc, exploratory (after the single unsealing): which v4 features carry the blind RNU2-2 result?
Leave-one-feature-out and single-feature models, refitted exactly as v4.py (non-blind training genes only), scored on
the blind genes. Never replaces the frozen A5 result; reported as an explanation of it.

    python v4_ablation.py --blind sealed/blind_labels.tsv     # -> results/v4_ablation.tsv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evaluate as E  # noqa: E402
import snrna as S  # noqa: E402
import v4 as V  # noqa: E402


def main():
    from sklearn.linear_model import LogisticRegression
    ap = argparse.ArgumentParser()
    ap.add_argument("--blind", required=True)
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--labels", default=os.path.join(HERE, "labels"))
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    a = ap.parse_args()
    ic = os.path.join(HERE, "model_inputs_v3") if os.path.exists(os.path.join(HERE, "model_inputs_v3", "variants.tsv.gz")) else a.inputs
    v = V.build_features(a.inputs, ic, a.evo2)
    genes = pd.read_csv(os.path.join(a.inputs, "genes.tsv.gz"), sep="\t")
    named = set(genes.label[genes.named & ~genes.gene.isin(S.SEALED_GENES)])
    L = E.load_l2(a.labels)
    v["y"] = v.ID.isin(set(L.ID[L.label == "P/LP"]))
    train = v[v.label.isin(named) & (v.y | (v.L3_benign_like.astype(bool) & ~v.y))]
    b = pd.read_csv(a.blind, sep="\t")
    bplp = set(b.ID[(b.label == "P/LP") & b.ID.str.match(E.SNV_ID)])
    t = v[v.gene == "RNU2-2"].copy()
    t["yb"] = t.ID.isin(bplp)
    t = t[t.yb | (t.L3_benign_like.astype(bool) & ~t.yb)]
    yb = t.yb.to_numpy()
    rows = [{"model": "CADD alone (raw)", "auroc_RNU2-2": S.auroc(yb, t.cadd_phred.to_numpy(float))}]
    sets = [("v4 full", V.FEATURES)] + [(f"drop {f}", [g for g in V.FEATURES if g != f]) for f in V.FEATURES] + \
           [(f"only {f}", [f]) for f in V.FEATURES]
    for name, fs in sets:
        X = [f + "_r" for f in fs]
        m = LogisticRegression(C=V.C_REG, class_weight="balanced", max_iter=5000).fit(train[X], train.y)
        rows.append({"model": name, "auroc_RNU2-2": S.auroc(yb, m.decision_function(t[X]))})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.results, "v4_ablation.tsv"), sep="\t", index=False, float_format="%.4f")
    print(f"RNU2-2 blind: {int(yb.sum())} P/LP vs {int((~yb).sum())} L3 (post hoc, exploratory)")
    print(out.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
