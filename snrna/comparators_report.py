#!/usr/bin/env python3
"""Comparator table for the paper (amendment A3; reported, never gating). Same endpoints as evaluate.py:
Aim 2 on the 435 RNU4-2 SGE SNVs, and per-unit AUROC of L2 P/LP SNVs vs L3 SNVs for the G3 units and RNU4-2.
Run after evaluate.py (the model's rows come from results/model_aim2.tsv and model_scores.tsv.gz).

    python comparators_report.py            # -> results/comparators_report.txt, comparators_aim2.tsv, comparators_units.tsv
    python comparators_report.py --shuffle  # code check with permuted labels (numbers meaningless)
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

import evaluate as E
import snrna as S

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = [("cadd_phred", "CADD v1.7"), ("phylop100", "phyloP 100 vertebrates"), ("phylop447", "phyloP 447 mammals"),
        ("remm", "ReMM v0.4"), ("gpn_msa", "GPN-MSA"), ("evo2_gene", "Evo 2 7B (gene window)")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--labels", default=os.path.join(HERE, "labels"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--shuffle", action="store_true")
    args = ap.parse_args()
    if args.shuffle:
        E.SHUFFLE = np.random.default_rng(S.SEED + 3)
    var = pd.read_csv(os.path.join(args.inputs, "variants.tsv.gz"), sep="\t")
    c = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    v = var[["ID", "label", "gene", "L3_benign_like"]].merge(c.drop(columns=["label", "gene"]), on="ID", how="left")
    v = v.merge(var[["ID", "cadd_phred"]], on="ID")
    dp = os.path.join(args.evo2, "deltas.tsv")
    if os.path.exists(dp):
        d = pd.read_csv(dp, sep="\t").drop_duplicates("ID")
        v["evo2_gene"] = -v.ID.map(dict(zip(d.ID, ((d.ll_alt_fwd - d.ll_ref_fwd) + (d.ll_alt_rc - d.ll_ref_rc)) / 2)))
    cols = [(k, n) for k, n in COMP if k in v.columns]
    rep = E.Report(os.path.join(args.results, "comparators_report.txt"),
                   "LABELS SHUFFLED: code check only, numbers meaningless" if args.shuffle else "")

    # Aim 2 endpoints (SGE labels; ReNU vs population)
    sge = pd.read_csv(os.path.join(args.inputs, "sge_eval.tsv.gz"), sep="\t")
    t = sge[sge.is_snv].reset_index(drop=True)
    if args.shuffle:
        t["function_score_mean"] = E.SHUFFLE.permutation(t.function_score_mean.to_numpy())
        t["depleted"] = E.SHUFFLE.permutation(t.depleted.to_numpy())
        t["category"] = E.SHUFFLE.permutation(t.category.to_numpy())
    rs = E.renu_sets(t)
    rows = []
    for k, name in cols:
        x = t.ID.map(dict(zip(v.ID, v[k]))).to_numpy(float)
        rows.append(E.metric_row(name, x, t, args.boot, rs))
    a2 = pd.DataFrame(rows)
    mp = os.path.join(args.results, "model_aim2.tsv")
    if os.path.exists(mp) and not args.shuffle:
        a2 = pd.concat([pd.read_csv(mp, sep="\t"), a2], ignore_index=True).drop_duplicates("predictor")
    a2.to_csv(os.path.join(args.results, "comparators_aim2.tsv"), sep="\t", index=False, float_format="%.4f")
    rep.add("Aim 2 endpoints, 435 RNU4-2 SGE SNVs (model rows from model_aim2.tsv where available)", "")
    for r in a2.to_dict("records"):
        rep.add("  " + E.fmt_row(r))
    rep.add("")

    # per-unit AUROC, L2 P/LP vs L3
    L = E.load_l2(args.labels)
    plp = L[L.label == "P/LP"].drop_duplicates("ID")
    ms_p = os.path.join(args.results, "model_scores.tsv.gz")
    ms = pd.read_csv(ms_p, sep="\t") if os.path.exists(ms_p) else None
    units = {"RNU4-2": ["RNU4-2"], "RNU4ATAC": ["RNU4ATAC"], "RNU12": ["RNU12"], "RNU6ATAC": ["RNU6ATAC"],
             "RNU6 paralogs (pooled)": sorted(g for g in var.gene.unique() if re.match(r"^RNU6-\d+$", str(g)))}
    urows = []
    rep.add("Per unit: AUROC of L2 P/LP SNVs vs L3 SNVs [95% CI]", "")
    for u, gl in units.items():
        x = v[v.gene.isin(gl)].copy()
        if ms is not None:
            x = x.merge(ms[["ID", "f_crossfit", "F_function"]], on="ID", how="left")
        pos = x.ID.isin(plp.ID)
        neg = x.L3_benign_like.astype(bool) & ~pos
        x = x[pos | neg]
        y = E.lab(x.ID.isin(plp.ID).to_numpy())
        if y.sum() < 5 or (~y).sum() < 5:
            rep.add(f"  {u}: {int(y.sum())} P vs {int((~y).sum())} L3 - too few")
            continue
        line = []
        extra = [("f_crossfit", "model f"), ("F_function", "model F")] if ms is not None else []
        for k, name in cols + extra:
            s = x[k].to_numpy(float)
            ok = np.isfinite(s)
            if ok.sum() < len(s) or not ok.any():
                continue
            a = S.auroc(y, s)
            lo, hi = S.boot_auroc(y, s, args.boot)
            urows.append({"unit": u, "score": name, "n_P": int(y.sum()), "n_B": int((~y).sum()), "auroc": a, "lo": lo, "hi": hi})
            line.append(f"{name} {a:.3f} [{lo:.2f}, {hi:.2f}]")
        rep.add(f"  {u} ({int(y.sum())} P vs {int((~y).sum())} L3): " + "; ".join(line))
    pd.DataFrame(urows).to_csv(os.path.join(args.results, "comparators_units.tsv"), sep="\t", index=False, float_format="%.4f")
    rep.write()


if __name__ == "__main__":
    main()
