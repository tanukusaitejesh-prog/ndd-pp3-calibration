#!/usr/bin/env python3
"""Post hoc (after the single unsealing): ClinGen PP3/BP4 calibration of v4, and the ClinVar VUS it would move.

Scores are out-of-sample for every labelled variant:
  - blind genes (RNU2-2, RNU5B-1): the A5-frozen v4 scores (hash-checked; results_v4/ is read only);
  - non-blind genes: cross-fitted v4, i.e. the frozen v4 recipe (features, C, class weights) refitted without the
    gene's whole paralog cluster, exactly as v4.py's internal validation. The frozen final model's own scores on its
    training genes are in-sample and are never calibrated.
Calibration is calibrate.py's local-LR method unchanged (PREREG s10): nearest 20% window, >= 10 per class, 1,000
bootstraps, PP3/BP4 at prior 0.10 (sensitivity priors as frozen). Units: pooled over genes (v4 features are
within-gene ranks, so one threshold is meaningful across genes), and each gene with >= 10 P/LP and >= 10 L3.
CADD is calibrated the same way for comparison.
VUS: labels/clinvar_vus_snrna.tsv (ClinVar VUS and conflicting SNVs in the named snRNA genes, fetched 14 Oct 2026),
scored with the frozen v4 (VUS were never labels), against the pooled thresholds.

    python v4_calibrate.py --blind sealed/blind_labels.tsv   # -> results/v4_calibration_*.tsv, results/v4_vus_*.tsv
"""
import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import calibrate as C  # noqa: E402
import evaluate as E  # noqa: E402
import snrna as S  # noqa: E402
import v4 as V  # noqa: E402

STRENGTHS = ("Strong", "Moderate", "Supporting")


def crossfit_v4(inputs, evo2, labels_dir):
    from sklearn.linear_model import LogisticRegression
    ic = os.path.join(HERE, "model_inputs_v3") if os.path.exists(os.path.join(HERE, "model_inputs_v3", "variants.tsv.gz")) else inputs
    v = V.build_features(inputs, ic, evo2)
    genes = pd.read_csv(os.path.join(inputs, "genes.tsv.gz"), sep="\t")
    cl = genes.set_index("label").cluster
    named = set(genes.label[genes.named & ~genes.gene.isin(S.SEALED_GENES)])
    L = E.load_l2(labels_dir)
    v["y"] = v.ID.isin(set(L.ID[L.label == "P/LP"]))
    train = v[v.label.isin(named) & (v.y | (v.L3_benign_like.astype(bool) & ~v.y))]
    X = [f + "_r" for f in V.FEATURES]
    out = pd.Series(np.nan, index=v.index)
    for i, lab in enumerate(sorted(named), 1):
        print(f"  cross-fit {i}/{len(named)}: {lab} (cluster {cl[lab]} held out)", flush=True)
        tr = train[train.label.map(cl) != cl[lab]]
        m = LogisticRegression(C=V.C_REG, class_weight="balanced", max_iter=5000).fit(tr[X], tr.y)
        idx = v.index[v.label == lab]
        out[idx] = m.decision_function(v.loc[idx, X])
    return v[["ID", "label", "gene"]].assign(v4_cf=out)


def calib(name, s, y, B, rows, curves, score):
    if y.sum() < C.MIN_CLASS or (~y).sum() < C.MIN_CLASS:
        return None
    grid, point, lo5, hi95 = C.lr_bootstrap(s, y, B)
    for g, p, a_, b_ in zip(grid, point, lo5, hi95):
        curves.append({"unit": name, "score": score, "s": g, "lr": p, "lr_lo5": a_, "lr_hi95": b_})
    th10 = None
    for prior in C.PRIORS:
        th = C.thresholds(grid, lo5, hi95, prior)
        rows.append({"unit": name, "score": score, "prior": prior, "n_P": int(y.sum()), "n_B": int((~y).sum()),
                     "auroc": S.auroc(y, s), **th})
        if prior == 0.10:
            th10 = th
    return th10


def evidence(s, th):
    if th is None or pd.isna(s):
        return ""
    for k in STRENGTHS:
        t = th.get(f"PP3_{k}")
        if t is not None and not pd.isna(t) and s >= t:
            return f"PP3_{k}"
    for k in STRENGTHS:
        t = th.get(f"BP4_{k}")
        if t is not None and not pd.isna(t) and s <= t:
            return f"BP4_{k}"
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blind", required=True)
    ap.add_argument("--v4", default=os.path.join(HERE, "results_v4"))
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--labels", default=os.path.join(HERE, "labels"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    sp = os.path.join(a.v4, "v4_scores.tsv.gz")
    man = {l.split()[1]: l.split()[0] for l in open(os.path.join(a.v4, "v4_manifest.txt")) if l.strip()}
    if man.get(os.path.relpath(sp, HERE)) != hashlib.sha256(open(sp, "rb").read()).hexdigest():
        sys.exit("v4_scores.tsv.gz does not match the frozen manifest")
    print("building v4 features and cross-fitting (non-blind genes) ...", flush=True)
    frozen = pd.read_csv(sp, sep="\t")[["ID", "gene", "v4"]]
    cf = crossfit_v4(a.inputs, a.evo2, a.labels)
    sc = cf.merge(frozen[["ID", "v4"]], on="ID", how="left")
    sc["v4_oos"] = np.where(sc.gene.isin(S.SEALED_GENES), sc.v4, sc.v4_cf)  # out-of-sample for every gene
    var = pd.read_csv(os.path.join(a.inputs, "variants.tsv.gz"), sep="\t")[["ID", "cadd_phred", "L3_benign_like"]]
    sc = sc.merge(var, on="ID", how="left")

    L = E.load_l2(a.labels)
    b = pd.read_csv(a.blind, sep="\t")
    b = b[b.ID.str.match(E.SNV_ID)]
    plp = set(L.ID[L.label == "P/LP"]) | set(b.ID[b.label == "P/LP"])
    sc["pos"] = sc.ID.isin(plp)
    t = sc[sc.pos | (sc.L3_benign_like.astype(bool) & ~sc.pos)].dropna(subset=["v4_oos", "cadd_phred"])

    rows, curves, th_v4 = [], [], {}
    units = {"pooled (all genes)": t, **{g: x for g, x in t.groupby("gene")}}
    for i, (name, x) in enumerate(units.items(), 1):
        print(f"  calibrating {i}/{len(units)}: {name} ({int(x.pos.sum())} P, {int((~x.pos).sum())} B), {a.boot} bootstraps", flush=True)
        y = x.pos.to_numpy(bool)
        th_v4[name] = calib(name, x.v4_oos.to_numpy(float), y, a.boot, rows, curves, "v4")
        calib(name, x.cadd_phred.to_numpy(float), y, a.boot, rows, curves, "cadd_phred")
    th = pd.DataFrame(rows)
    th.to_csv(os.path.join(a.results, "v4_calibration_thresholds.tsv"), sep="\t", index=False, float_format="%.5g")
    pd.DataFrame(curves).to_csv(os.path.join(a.results, "v4_calibration_lr_curves.tsv.gz"), sep="\t", index=False, float_format="%.5g")

    rep = ["v4 calibration (post hoc; out-of-sample v4: frozen for blind genes, cluster-cross-fitted otherwise). Prior 0.10", ""]
    for r in th[th.prior == 0.10].itertuples():
        pp3 = ", ".join(f"{k[:3]} >= {getattr(r, 'PP3_' + k):.3g}" for k in STRENGTHS[::-1] if getattr(r, "PP3_" + k) is not None and not pd.isna(getattr(r, "PP3_" + k))) or "none"
        bp4 = ", ".join(f"{k[:3]} <= {getattr(r, 'BP4_' + k):.3g}" for k in STRENGTHS[::-1] if getattr(r, "BP4_" + k) is not None and not pd.isna(getattr(r, "BP4_" + k))) or "none"
        rep.append(f"  {r.unit:22s} {r.score:10s} {r.n_P:4d} P {r.n_B:4d} B  AUROC {r.auroc:.3f}  PP3 {pp3} | BP4 {bp4}")

    # VUS: frozen v4 (never trained on VUS) for every gene; thresholds: the gene's own if calibrated, else pooled
    vp = os.path.join(a.labels, "clinvar_vus_snrna.tsv")
    if os.path.exists(vp):
        vus = pd.read_csv(vp, sep="\t")
        vus = vus[vus.ID.str.match(E.SNV_ID) & ~vus.ID.isin(plp)].drop_duplicates("ID")
        vus = vus.merge(sc[["ID", "v4_oos"]], on="ID", how="left").merge(var[["ID", "cadd_phred"]], on="ID", how="left")
        th_cadd = {r.unit: r._asdict() for r in th[(th.prior == 0.10) & (th.score == "cadd_phred")].itertuples()}
        pooled = th_v4["pooled (all genes)"]
        vus["v4_evidence_pooled"] = [evidence(s, pooled) for s in vus.v4_oos]
        vus["v4_evidence_gene"] = [evidence(s, th_v4.get(g)) for s, g in zip(vus.v4_oos, vus.gene)]
        vus["cadd_evidence_pooled"] = [evidence(s, th_cadd.get("pooled (all genes)")) for s in vus.cadd_phred]
        vus.to_csv(os.path.join(a.results, "v4_vus_evidence.tsv"), sep="\t", index=False, float_format="%.5g")
        summ = pd.concat([vus.groupby("gene").v4_evidence_pooled.value_counts().unstack(fill_value=0),
                          vus.groupby("gene").size().rename("n_vus")], axis=1).fillna(0).astype(int)
        summ.to_csv(os.path.join(a.results, "v4_vus_summary.tsv"), sep="\t")
        rep += ["", f"ClinVar VUS/conflicting SNVs in named snRNA genes: {len(vus)} (scored: {int(vus.v4_oos.notna().sum())})",
                "  v4 pooled thresholds: " + ", ".join(f"{k} {v}" for k, v in vus.v4_evidence_pooled.replace("", "none").value_counts().items()),
                "  v4 own-gene thresholds: " + ", ".join(f"{k} {v}" for k, v in vus.v4_evidence_gene.replace("", "none").value_counts().items()),
                "  CADD pooled thresholds: " + ", ".join(f"{k} {v}" for k, v in vus.cadd_evidence_pooled.replace("", "none").value_counts().items()),
                "", summ.to_string()]
    open(os.path.join(a.results, "v4_calibration_report.txt"), "w").write("\n".join(rep) + "\n")
    json.dump({"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "stage": "POST HOC v4 calibration",
               "git_commit": E.git_commit(), "boot": a.boot}, open(os.path.join(a.results, "v4_calibration_run_info.json"), "w"), indent=1)
    print("\n".join(rep))


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
