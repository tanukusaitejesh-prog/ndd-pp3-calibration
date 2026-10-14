#!/usr/bin/env python3
"""Amendment A6: predicting inheritance mode (dominant vs recessive action) of pathogenic snRNA variants from
population constraint. Frozen before the single unsealing; called by final_eval.py.

Rationale: dominant variants are removed from the population in heterozygotes (local depletion); recessive ones are
tolerated in carriers. Non-blind check before freezing (CPU v4 features): local constraint separated dominant from
recessive L2 P/LP SNVs with AUROC 1.0 (14 vs 62; within RNU4-2 14 vs 11), conservation did not (phyloP447 0.17-0.62).

Score: `local_constraint` from the frozen results_v4/v4_scores.tsv.gz (hash checked against v4_manifest.txt).
Threshold: midpoint between the highest non-blind recessive and the lowest non-blind dominant L2 score (if the
classes overlap: the value maximising balanced accuracy on non-blind L2), fixed before any blind label is read.
Blind test: P/LP SNVs of RNU2-2 and RNU5B-1 with an inheritance annotation in the decrypted labels.
Primary: AUROC dominant vs recessive (needs >= 3 of each); secondary: accuracy of the threshold rule, per gene.

    python inheritance_blind_eval.py --blind sealed/blind_labels.tsv
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


def threshold(x, dom):
    lo_dom, hi_rec = x[dom].min(), x[~dom].max()
    if hi_rec < lo_dom:
        return (hi_rec + lo_dom) / 2
    cand = np.unique(x)
    ba = [((x[dom] >= t).mean() + (x[~dom] < t).mean()) / 2 for t in cand]
    return float(cand[int(np.argmax(ba))])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blind", required=True)
    ap.add_argument("--v4", default=os.path.join(HERE, "results_v4"))
    ap.add_argument("--labels", default=os.path.join(HERE, "labels"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    a = ap.parse_args()
    sp = os.path.join(a.v4, "v4_scores.tsv.gz")
    man = {l.split()[1]: l.split()[0] for l in open(os.path.join(a.v4, "v4_manifest.txt")) if l.strip()}
    if man.get(os.path.relpath(sp, HERE)) != hashlib.sha256(open(sp, "rb").read()).hexdigest():
        sys.exit("v4_scores.tsv.gz does not match the frozen manifest; not evaluating")
    sc = pd.read_csv(sp, sep="\t")[["ID", "gene", "local_constraint"]]
    L = E.load_l2(a.labels)
    L = L[(L.label == "P/LP") & L.inheritance.isin(["dominant", "recessive"])].drop_duplicates("ID")
    nb = L[["ID", "inheritance"]].merge(sc, on="ID").dropna(subset=["local_constraint"])
    x, dom = nb.local_constraint.to_numpy(float), (nb.inheritance == "dominant").to_numpy()
    t = threshold(x, dom)
    lines = [f"A6 inheritance-mode prediction. Non-blind: {dom.sum()} dominant vs {(~dom).sum()} recessive, "
             f"AUROC {S.auroc(dom, x):.3f}; threshold {t:.3f} (local_constraint >= t -> dominant)"]
    b = pd.read_csv(a.blind, sep="\t")
    b = b[b.ID.str.match(E.SNV_ID) & (b.label == "P/LP")]
    if "inheritance" not in b.columns:
        lines.append("Blind labels carry no inheritance column: not evaluable")
    else:
        b = b[b.inheritance.isin(["dominant", "recessive"])].drop_duplicates("ID")[["ID", "inheritance"]].merge(sc, on="ID")
        b = b.dropna(subset=["local_constraint"])
        rows = []
        for g, part in list(b.groupby("gene")) + [("both", b)]:
            xb, db = part.local_constraint.to_numpy(float), (part.inheritance == "dominant").to_numpy()
            r = {"gene": g, "n_dominant": int(db.sum()), "n_recessive": int((~db).sum()),
                 "accuracy": float(((xb >= t) == db).mean()) if len(db) else np.nan}
            if db.sum() >= 3 and (~db).sum() >= 3:
                r["auroc"] = S.auroc(db, xb)
                r["auroc_lo"], r["auroc_hi"] = S.boot_auroc(db, xb, 2000)
            rows.append(r)
        t2 = pd.DataFrame(rows)
        t2.to_csv(os.path.join(a.results, "final_inheritance_blind.tsv"), sep="\t", index=False, float_format="%.4f")
        lines.append(t2.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    open(os.path.join(a.results, "final_inheritance_blind.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
