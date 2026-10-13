#!/usr/bin/env python3
"""Aim 4: ClinGen-style calibration of PP3/BP4 evidence (PREREG s10, amendment A3). Runs after evaluate.py.

    python calibrate.py                                  # non-blind genes (preliminary; before the single unsealing)
    python calibrate.py --blind sealed/blind_labels.tsv  # called by final_eval.py after the single unsealing

Scores (higher = more damaging): the framework score S per gene (A2.8b: selection f for dominant genes, function
score F for recessive genes, mean rank of both for mixed genes, f where F is unavailable), its two components, and
the comparators CADD, ViennaRNA dDG (U4/U6 and U4atac/U6atac only), Evo 2, phyloP447, ReMM and GPN-MSA.
Labels: L2 P/LP SNVs (positives) against L3 SNVs (gnomAD PASS AC >= 9, minus L2).

Local LR at every observed score: the nearest 20% of labelled variants by score rank (>= 10 per class), 0.5
pseudocounts; 1,000 class-stratified bootstrap replicates (seed 20261013). PP3 at strength k: the smallest score s*
such that the 5th percentile of the local LR clears the threshold at every score >= s*. BP4 mirrors it with the
95th percentile <= the reciprocal at every score <= s*. Thresholds O_VSt(P)^(1/8, 1/4, 1/2); P = 0.10 primary,
{0.01, 0.0441, 0.2, 0.5} sensitivity. Units: each gene with >= 10 positives and >= 10 negatives, the dominant and
recessive pools, and dominant/recessive positives separately within a gene. Positives whose source classification
used in-silico or conservation evidence are flagged; results are shown with and without them. Ancestry drift: L3
redefined per gnomAD genetic ancestry group (AC_g >= 9 * AN_g / AN), Supporting thresholds per group.
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd

import evaluate as E
import snrna as S
from labels import GENE_MODE

HERE = os.path.dirname(os.path.abspath(__file__))
PRIORS = (0.10, 0.01, 0.0441, 0.2, 0.5)
STRENGTHS = ("Supporting", "Moderate", "Strong")
ANCESTRIES = S.ANCESTRIES
MIN_CLASS = 10
CALIB_MODE = {**GENE_MODE, "RNU4-2": "mixed", "RNU2-2": "mixed", "RNU6-7": "dominant"}
POOLS = {"dominant pool (f)": ("dominant", "f"), "recessive pool (F)": ("recessive", "F")}  # by each positive's inheritance


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ------------------------------------------------------------------------------------------------ scores

def gene_scores(results, inputs, evo2_dir):
    """Per-variant scores: framework S (per gene mode), f, F and comparators. Cross-fitted model scores come from
    evaluate.py's model_scores.tsv.gz."""
    ms = pd.read_csv(os.path.join(results, "model_scores.tsv.gz"), sep="\t")
    var = pd.read_csv(os.path.join(inputs, "variants.tsv.gz"), sep="\t")
    sc = ms.merge(var[["ID", "L3_benign_like", "ddG_bind_sum", "AC_pass"]], on="ID")
    sc["f"] = sc.f_crossfit
    sc["F"] = sc.F_function
    sc["mode"] = sc.gene.map(CALIB_MODE).fillna("unknown")
    parts = []
    for g, x in sc.groupby("gene"):
        x = x.copy()
        mode = x["mode"].iloc[0]
        if mode == "recessive" and x.F.notna().any():
            x["S"] = x.F
        elif mode in ("mixed", "unknown") and x.F.notna().all() and x.f.notna().all():
            x["S"] = (x.f.rank(pct=True) + x.F.rank(pct=True)) / 2
        else:
            x["S"] = x.f
        parts.append(x)
    sc = pd.concat(parts, ignore_index=True)
    fam = sc.gene.map(lambda g: "U4" if g.startswith("RNU4") else "")
    sc["ddG"] = np.where(fam == "U4", sc.ddG_bind_sum, np.nan)
    dp = os.path.join(evo2_dir, "deltas.tsv")
    if os.path.exists(dp):
        d = pd.read_csv(dp, sep="\t").drop_duplicates("ID")
        delta = ((d.ll_alt_fwd - d.ll_ref_fwd) + (d.ll_alt_rc - d.ll_ref_rc)) / 2
        sc["evo2"] = -sc.ID.map(dict(zip(d.ID, delta)))
    cp = os.path.join(HERE, "comparators", "comparators.tsv.gz")
    if os.path.exists(cp):
        c = pd.read_csv(cp, sep="\t")
        sc = sc.merge(c[["ID", "phylop447", "remm", "gpn_msa"]], on="ID", how="left")
    return sc


SCORES = [("S", "Framework score S [model]"), ("f", "Selection f [model component]"), ("F", "Function F [model component]"),
          ("cadd_phred", "CADD v1.7"), ("ddG", "ViennaRNA dDG (U4/U4atac)"), ("evo2", "Evo 2 7B"),
          ("phylop447", "phyloP 447 mammals"), ("remm", "ReMM v0.4"), ("gpn_msa", "GPN-MSA")]


# ---------------------------------------------------------------------------------------------- labels

def labels(labels_dir, blind_path):
    L = E.load_l2(labels_dir)
    if blind_path:
        b = pd.read_csv(blind_path, sep="\t")
        b = b[b.ID.str.match(E.SNV_ID)]
        L = pd.concat([L, b], ignore_index=True)
    if "insilico_or_conservation_used_for_classification" in L.columns:
        L["flagged"] = L.insilico_or_conservation_used_for_classification.astype(str).str.lower().eq("yes")
    else:
        L["flagged"] = False
    plp = L[L.label == "P/LP"].sort_values("flagged").drop_duplicates("ID")
    return plp


def unit_table(sc, plp):
    """Labelled variants with class (1 = L2 positive, 0 = L3 negative), inheritance and flag."""
    t = sc.copy()
    t["pos"] = t.ID.isin(plp.ID)
    t["neg"] = t.L3_benign_like.astype(bool) & ~t.pos
    t = t[t.pos | t.neg].copy()
    inh = dict(zip(plp.ID, plp.inheritance.fillna("unknown")))
    flg = dict(zip(plp.ID, plp.flagged))
    t["inheritance"] = t.ID.map(inh).fillna("")
    t["flagged"] = t.ID.map(flg).fillna(False).astype(bool)
    return t


# ----------------------------------------------------------------------------------------- calibration

def lr_bootstrap(s, y, B, seed=S.SEED):
    grid = np.unique(s)
    point = S.local_lr(s, y, grid)
    rng = np.random.default_rng(seed)
    i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
    boots = np.empty((B, len(grid)))
    for b in range(B):
        idx = np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))])
        boots[b] = S.local_lr(s[idx], y[idx], grid)
    return grid, point, np.percentile(boots, 5, axis=0), np.percentile(boots, 95, axis=0)


def thresholds(grid, lo5, hi95, prior):
    out = {}
    for k in STRENGTHS:
        need = S.lr_thresholds(prior)[k]
        ok = lo5 >= need
        suffix = np.flip(np.cumprod(np.flip(ok)).astype(bool))
        out[f"PP3_{k}"] = float(grid[np.argmax(suffix)]) if suffix.any() else None
        okb = hi95 <= 1 / need
        prefix = np.cumprod(okb).astype(bool)
        out[f"BP4_{k}"] = float(grid[np.flatnonzero(prefix)[-1]]) if prefix.any() else None
    return out


def calibrate_unit(name, t, col, B, curves, rows, extra=None):
    x = t[t[col].notna()]
    y = x.pos.to_numpy(bool)
    if y.sum() < MIN_CLASS or (~y).sum() < MIN_CLASS:
        return None
    s = x[col].to_numpy(float)
    grid, point, lo5, hi95 = lr_bootstrap(s, E.lab(y), B)
    for g, p, a, b in zip(grid, point, lo5, hi95):
        curves.append({"unit": name, "score": col, **(extra or {}), "s": g, "lr": p, "lr_lo5": a, "lr_hi95": b})
    res = {}
    for prior in PRIORS:
        th = thresholds(grid, lo5, hi95, prior)
        rows.append({"unit": name, "score": col, **(extra or {}), "prior": prior, "n_P": int(y.sum()), "n_B": int((~y).sum()),
                     "auroc": S.auroc(E.lab(y), s), **th})
        if prior == 0.10:
            res = th
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=os.path.join(HERE, "model_runs"))
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--labels", default=os.path.join(HERE, "labels"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--blind", default=None, help="decrypted blind labels (final evaluation only)")
    ap.add_argument("--boot", type=int, default=1000)
    args = ap.parse_args()
    info = json.load(open(os.path.join(args.runs, "data_info.json")))
    synthetic = info.get("synthetic_counts") is not None
    E.SHUFFLE = np.random.default_rng(S.SEED + 2) if synthetic else None
    sc = gene_scores(args.results, args.inputs, args.evo2)
    plp = labels(args.labels, args.blind)
    t = unit_table(sc, plp)
    rows, curves, anc_rows = [], [], []
    stage = "FINAL (blind genes included)" if args.blind else "PRELIMINARY (non-blind genes only)"
    rep = E.Report(os.path.join(args.results, "calibration_report.txt"),
                   ("DRY RUN: synthetic counts and SHUFFLED labels; numbers meaningless. " if synthetic else "") + stage)
    rep.add(f"Aim 4 calibration (PREREG s10, A3); bootstrap {args.boot}; positives = L2 P/LP SNVs, negatives = L3", "")

    units = {}
    for g, x in t.groupby("gene"):
        units[g] = x
    for pool, (mode, col) in POOLS.items():
        # positives of this inheritance mode in any gene, against the L3 variants of the genes that contribute them
        genes_m = set(t.gene[t.pos & (t.inheritance == mode)])
        units[pool] = t[(t.pos & (t.inheritance == mode)) | (t.neg & t.gene.isin(genes_m))]
    model_found = []
    for name, x in units.items():
        cols = [c for c, _ in SCORES if c in x.columns]
        if name in POOLS:
            cols = [POOLS[name][1], "cadd_phred"]
        for flag_set in ("all", "unflagged"):
            xx = x if flag_set == "all" else x[~(x.pos & x.flagged)]
            for col in cols:
                th = calibrate_unit(name, xx, col, args.boot, curves, rows, {"positives": flag_set})
                if th is not None and col in ("S", "f", "F") and flag_set == "all" and (name in POOLS or col == "S"):
                    if th["PP3_Supporting"] is not None or th["BP4_Supporting"] is not None:
                        model_found.append(f"{name} ({col})")
        # dominant and recessive positives separately, within a gene
        if name not in POOLS:
            for inh in ("dominant", "recessive"):
                xx = x[~x.pos | (x.inheritance == inh)]
                if (xx.pos.sum() >= MIN_CLASS) and (x.pos & (x.inheritance != inh)).any():
                    calibrate_unit(f"{name} [{inh} positives]", xx, "S", args.boot, curves, rows, {"positives": "all"})
    th = pd.DataFrame(rows)
    th.to_csv(os.path.join(args.results, "calibration_thresholds.tsv"), sep="\t", index=False, float_format="%.5g")
    pd.DataFrame(curves).to_csv(os.path.join(args.results, "calibration_lr_curves.tsv.gz"), sep="\t", index=False, float_format="%.5g")

    # ancestry drift: model score S per gene, Supporting thresholds with L3 redefined per group
    anc = pd.read_csv(os.path.join(HERE, "calibration_inputs", "ancestry_counts.tsv.gz"), sep="\t").set_index("ID")
    for name, x in units.items():
        if name in POOLS:
            continue
        for a in ANCESTRIES:
            ac, an = anc.reindex(x.ID)[f"AC_{a}"].to_numpy(float), anc.reindex(x.ID)[f"AN_{a}"].to_numpy(float)
            tot = anc.reindex(x.ID).AN.to_numpy(float)
            passed = (anc.reindex(x.ID)["filter"].fillna("") == "PASS").to_numpy()
            neg_a = passed & (ac >= np.maximum(1, 9 * np.divide(an, tot, out=np.zeros_like(an), where=tot > 0))) & ~x.pos.to_numpy()
            xa = x[x.pos.to_numpy() | neg_a]
            if xa.pos.sum() >= MIN_CLASS and (~xa.pos).sum() >= MIN_CLASS and xa.S.notna().all():
                y = E.lab(xa.pos.to_numpy())
                grid, point, lo5, hi95 = lr_bootstrap(xa.S.to_numpy(float), y, args.boot)
                tt = thresholds(grid, lo5, hi95, 0.10)
                anc_rows.append({"unit": name, "ancestry": a, "n_P": int(y.sum()), "n_B": int((~y).sum()),
                                 "PP3_Supporting": tt["PP3_Supporting"], "BP4_Supporting": tt["BP4_Supporting"]})
    pd.DataFrame(anc_rows).to_csv(os.path.join(args.results, "calibration_ancestry.tsv"), sep="\t", index=False, float_format="%.5g")

    # report
    prim = th[(th.prior == 0.10) & (th.positives == "all")]
    for name in units:
        sub = prim[prim.unit == name]
        if not len(sub):
            continue
        rep.add(f"{name}: {int(sub.n_P.iloc[0])} P vs {int(sub.n_B.iloc[0])} L3 (prior 0.10)")
        for r in sub.itertuples():
            label = dict(SCORES).get(r.score, r.score)
            pp3 = ", ".join(f"{k[0]} >= {getattr(r, 'PP3_' + k):.3g}" for k in STRENGTHS if getattr(r, "PP3_" + k) is not None) or "none"
            bp4 = ", ".join(f"{k[0]} <= {getattr(r, 'BP4_' + k):.3g}" for k in STRENGTHS if getattr(r, "BP4_" + k) is not None) or "none"
            rep.add(f"  {label:32s} AUROC {r.auroc:.3f}  PP3 {pp3}  |  BP4 {bp4}")
        rep.add("")
    met = bool(model_found)
    rep.add(f"Aim 4 decision ({stage}): {'MET' if met else 'NOT MET'} - model Supporting PP3 or BP4 threshold in: "
            f"{', '.join(model_found) if model_found else 'no gene or pool'}"
            + (" (pooled results are stated as pooled)" if any('pool' in m for m in model_found) else ""), "")
    if anc_rows:
        rep.add("Ancestry drift (framework S, prior 0.10, Supporting): see calibration_ancestry.tsv; reported, not corrected")
    rep.write()

    # clinician table: framework score per gene, primary prior, all positives
    clin = prim[prim.score.isin(["S", "f", "F"]) & ~prim.unit.str.contains(r"\[")].copy()
    clin.to_csv(os.path.join(args.results, "clinician_thresholds.tsv"), sep="\t", index=False, float_format="%.4g")
    with open(os.path.join(args.results, "calibration_run_info.json"), "w") as f:
        json.dump({"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "stage": stage, "boot": args.boot,
                   "synthetic": synthetic, "aim4_met": met, "git_commit": E.git_commit()}, f, indent=1)


if __name__ == "__main__":
    main()
