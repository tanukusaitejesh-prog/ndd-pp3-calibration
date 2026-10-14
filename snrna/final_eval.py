#!/usr/bin/env python3
"""Final evaluation (PREREG s9, amendment A3): the single unsealing of the blind genes' labels. Run once, after
train.py and evaluate.py have finished and the results have been committed. Nothing is refitted afterwards.

    SNRNA_SEAL_PASS=... python final_eval.py --reason "final evaluation, PREREG s9"

1. decrypts sealed/blind_labels.tsv.enc with `snrna.py unseal` (logged in sealed/SEALED_MANIFEST.jsonl)
2. RNU2-2 and RNU5B-1, each scored by the fit that held out its whole paralog cluster: AUROC of L2 P/LP SNVs vs
   L3 SNVs and delta AUROC vs CADD, for the framework score S (RNU5B-1 dominant: f; RNU2-2 mixed: mean rank of f
   and F) and for every component and comparator
3. Aim 3 summary with the blind genes added as a second, blind estimate (weighted delta AUROC vs CADD)
4. calibrate.py with all labels -> final Aim 4
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time

import numpy as np
import pandas as pd

import evaluate as E
import snrna as S

HERE = os.path.dirname(os.path.abspath(__file__))
SEALED = os.path.join(HERE, "sealed")


def enc_path():
    return os.path.join(SEALED, "blind_labels.tsv.enc")


def plain_path():
    return os.path.join(SEALED, "blind_labels.tsv")
BLIND_MODE = {"RNU5B-1": "dominant", "RNU2-2": "mixed"}


def manifest_events():
    p = os.path.join(SEALED, "SEALED_MANIFEST.jsonl")
    return [json.loads(x) for x in open(p) if x.strip()] if os.path.exists(p) else []


def unseal(reason):
    ev = manifest_events()
    if not any(e.get("event") == "SEALED" for e in ev):
        sys.exit("no SEALED event in sealed/SEALED_MANIFEST.jsonl")
    if any(e.get("event") == "UNSEALED" for e in ev):
        if os.path.exists(plain_path()):
            print("already unsealed once; continuing with the decrypted file of that unsealing", flush=True)
            return
        sys.exit("the labels were unsealed before and the decrypted file is gone; decrypt manually and log it")
    if not os.environ.get("SNRNA_SEAL_PASS"):
        sys.exit("set SNRNA_SEAL_PASS (the sealing passphrase)")
    subprocess.run([sys.executable, os.path.join(HERE, "snrna.py"), "unseal", enc_path(), "--reason", reason,
                    "--pass-env", "SNRNA_SEAL_PASS"], check=True)


def blind_gene(args, pred, var, genes, cl, gene, plp, B, hom):
    labs = list(genes.label[genes.gene == gene])
    sc = E.unit_scores(pred, var, genes, labs, lambda lab_: f"lofo_{cl[lab_]}", hom)
    if sc is None:
        return None, None
    sc["F"] = sc.h.where(sc.h.notna(), sc.d)
    sc["S"] = E.framework_score(sc, BLIND_MODE[gene])
    comp = os.path.join(HERE, "comparators", "comparators.tsv.gz")
    if os.path.exists(comp):
        c = pd.read_csv(comp, sep="\t")
        sc = sc.merge(c[["ID", "phylop447", "remm", "gpn_msa"]], on="ID", how="left")
    pos = sc.ID.isin(plp.ID)
    neg = sc.L3_benign_like.astype(bool) & ~pos
    ok = (pos | neg) & sc.S.notna() & sc.cadd_phred.notna()
    y = pos[ok].to_numpy()
    row = {"gene": gene, "mode": BLIND_MODE[gene], "fit": f"lofo_{cl[labs[0]]}", "n_L2_snv": int(y.sum()), "n_L3": int((~y).sum())}
    if not (y.sum() and (~y).sum()):
        return row, None
    for col in ("S", "f", "F", "cadd_phred", "phylop447", "remm", "gpn_msa"):
        x = sc[col][ok]
        if x.notna().all():
            row[f"auroc_{col}"] = S.auroc(y, x.to_numpy(float))
            lo, hi = S.boot_auroc(y, x.to_numpy(float), B)
            row[f"auroc_{col}_lo"], row[f"auroc_{col}_hi"] = lo, hi
    row["delta_S_vs_cadd"] = row["auroc_S"] - row["auroc_cadd_phred"]
    return row, (y, sc.S[ok].to_numpy(float), sc.cadd_phred[ok].to_numpy(float))


def g3_sets(args, pred, var, genes, cl, plp, hom):
    """The pre-registered G3 units (eligible ones) with the gating score, as in evaluate.g3."""
    out = []
    named = genes[genes.named]
    units = {"RNU4ATAC": ["RNU4ATAC"], "RNU6ATAC": ["RNU6ATAC"], "RNU12": ["RNU12"],
             "RNU6 paralogs (pooled)": sorted(x for x in named.gene if re.match(r"^RNU6-\d+$", x))}
    for u, gl in units.items():
        sc = E.unit_scores(pred, var, genes, list(genes.label[genes.gene.isin(gl)]), lambda lab_: f"lofo_{cl[lab_]}", hom)
        if sc is None:
            continue
        sc["S"] = E.framework_score(sc, E.G3_UNITS[u])
        pos = sc.ID.isin(plp.ID)
        neg = sc.L3_benign_like.astype(bool) & ~pos
        ok = (pos | neg) & sc.S.notna() & sc.cadd_phred.notna()
        y = pos[ok].to_numpy()
        if y.sum() >= 10 and (~y).sum():
            out.append((u, y, sc.S[ok].to_numpy(float), sc.cadd_phred[ok].to_numpy(float)))
    return out


def weighted_delta(sets, B):
    w = np.array([y.sum() for _, y, _, _ in sets], float)
    point = float((w * np.array([S.auroc(y, a) - S.auroc(y, c) for _, y, a, c in sets])).sum() / w.sum())
    rng = np.random.default_rng(S.SEED)
    bs = np.empty(B)
    for b in range(B):
        vals = []
        for _, y, a, c in sets:
            i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
            idx = np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))])
            vals.append(S.auroc(y[idx], a[idx]) - S.auroc(y[idx], c[idx]))
        bs[b] = (w * np.array(vals)).sum() / w.sum()
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return point, lo, hi, int(w.sum())


def main():
    global SEALED
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reason", required=True)
    ap.add_argument("--runs", default=os.path.join(HERE, "model_runs"))
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--labels", default=os.path.join(HERE, "labels"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--calib-boot", type=int, default=1000)
    ap.add_argument("--sealed", default=SEALED, help="directory of the sealed labels (tests only)")
    args = ap.parse_args()
    SEALED = args.sealed
    info = json.load(open(os.path.join(args.runs, "data_info.json")))
    if info.get("synthetic_counts") is not None or info.get("quick"):
        sys.exit("final evaluation needs the real, full run (this runs directory is a dry run)")
    if not os.path.exists(os.path.join(args.results, "model_scores.tsv.gz")):
        sys.exit("run evaluate.py first")
    unseal(args.reason)
    blind = pd.read_csv(plain_path(), sep="\t")
    blind = blind[blind.ID.str.match(E.SNV_ID)]
    bplp = blind[blind.label == "P/LP"].drop_duplicates("ID")
    pred, _ = E.load_preds(args.runs)
    var = pd.read_csv(os.path.join(args.inputs, "variants.tsv.gz"), sep="\t")
    genes = pd.read_csv(os.path.join(args.inputs, "genes.tsv.gz"), sep="\t")
    cl = genes.set_index("label").cluster
    hom = E.homology_scores(args.inputs)
    rep = E.Report(os.path.join(args.results, "final_report.txt"), "")
    rep.add("Final evaluation (PREREG s9; single unsealing). Blind genes scored by fits that held out their paralog cluster", "")
    rows, blind_sets = [], []
    for gene in ("RNU2-2", "RNU5B-1"):
        row, st = blind_gene(args, pred, var, genes, cl, gene, bplp, args.boot, hom)
        if row is None:
            rep.add(f"  {gene}: no cross-fitted predictions")
            continue
        rows.append(row)
        if st is None:
            rep.add(f"  {gene}: L2 {row['n_L2_snv']} vs L3 {row['n_L3']} SNVs - not evaluable")
            continue
        blind_sets.append((gene, *st))
        rep.add(f"  {gene} ({row['mode']}; {row['fit']}): {row['n_L2_snv']} P/LP vs {row['n_L3']} L3 SNVs")
        for col, name in (("S", "framework S"), ("f", "selection f"), ("F", "function F"), ("cadd_phred", "CADD"),
                          ("phylop447", "phyloP447"), ("remm", "ReMM"), ("gpn_msa", "GPN-MSA")):
            if f"auroc_{col}" in row:
                rep.add(f"    {name:12s} AUROC {row[f'auroc_{col}']:.3f} [{row[f'auroc_{col}_lo']:.3f}, {row[f'auroc_{col}_hi']:.3f}]")
        rep.add(f"    delta AUROC (S - CADD) {row['delta_S_vs_cadd']:+.3f}", "")
    pd.DataFrame(rows).to_csv(os.path.join(args.results, "final_blind_genes.tsv"), sep="\t", index=False, float_format="%.4f")
    plp = E.load_l2(args.labels)
    plp = plp[plp.label == "P/LP"].drop_duplicates("ID")
    week6 = g3_sets(args, pred, var, genes, cl, plp, hom)
    if week6:
        p, lo, hi, n = weighted_delta(week6, args.boot)
        rep.add(f"Aim 3, week-6 units only ({', '.join(u for u, *_ in week6)}; {n} L2 SNVs): weighted delta AUROC "
                f"{p:+.3f} [{lo:+.3f}, {hi:+.3f}]")
    if blind_sets:
        p, lo, hi, n = weighted_delta(blind_sets, args.boot)
        rep.add(f"Aim 3, blind genes only ({n} L2 SNVs): weighted delta AUROC {p:+.3f} [{lo:+.3f}, {hi:+.3f}]")
        if week6:
            p, lo, hi, n = weighted_delta(week6 + blind_sets, args.boot)
            rep.add(f"Aim 3, week-6 units + blind genes ({n} L2 SNVs): weighted delta AUROC {p:+.3f} [{lo:+.3f}, {hi:+.3f}] "
                    "(second, blind estimate; the G3 decision itself stays as made at evaluate.py)")
    rep.add("")
    rep.write()
    subprocess.run([sys.executable, os.path.join(HERE, "calibrate.py"), "--runs", args.runs, "--inputs", args.inputs,
                    "--labels", args.labels, "--results", args.results, "--evo2", args.evo2, "--blind", plain_path(),
                    "--boot", str(args.calib_boot)], check=True)
    if os.path.exists(os.path.join(HERE, "results_v4", "v4_manifest.txt")):  # amendment A5: frozen v4, same single unsealing
        subprocess.run([sys.executable, os.path.join(HERE, "v4_blind_eval.py"), "--blind", plain_path(),
                        "--results", args.results, "--boot", str(args.boot)], check=False)
    with open(os.path.join(args.results, "final_run_info.json"), "w") as f:
        json.dump({"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "reason": args.reason, "git_commit": E.git_commit(),
                   "boot": args.boot}, f, indent=1)


if __name__ == "__main__":
    main()
