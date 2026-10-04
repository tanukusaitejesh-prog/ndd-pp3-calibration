#!/usr/bin/env python3
"""fullstudy.py - full-study analyses on top of the frozen feasibility gate (gate.py is imported, never modified).

Plan: FULL_STUDY_PLAN.md (frozen and timestamped before any full-study result).

  python fullstudy.py universe            per-site coding universe, identical to the gate's (validated against its aggregates)
  python fullstudy.py scores              dbNSFP 4.1a rankscores joined to universe missense sites
  python fullstudy.py regions             gnomAD v2.1.1 regional missense constraint + CCR, lifted to GRCh38, per site
  python fullstudy.py aim5|aim2|aim1|aim3|aim4|mavedb
"""
import argparse
import gc
import gzip
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gate  # noqa: E402  (frozen; only its functions and constants are reused)
from gate import P, B2I, AUTOSOMES, RATE_MODELS, PRIMARY_MODEL, SEED, PP3, log, die  # noqa: E402

EXT = "/mnt/x/ndd_gate_data"            # external downloads (dbNSFP, constraint, mechanism labels)
DELTA = 0.091                           # systematic margin, fixed from the gate run (FULL_STUDY_PLAN.md section 2)
N_BOOT = 1000
N_QBINS = 1000                          # fine score bins: quantiles of the score weighted by mutational opportunity
WIN_HALF = 25                           # local window = 51 fine bins (~5% of mutational opportunity)
WIN_STEP = 5


# ======================================================================== universe
def build_universe(work, c, tx, seg):
    """Re-runs gate.process_chrom's universe construction (same code path, same functions) and returns per-site rows."""
    r, n_raw = gate.load_roulette(work, c)
    pos_u, first = np.unique(r["pos"].to_numpy(), return_index=True)
    ref_u = r["rc"].to_numpy()[first].astype(np.int64)
    ok_tx = tx.loc[(tx["chrom"] == c) & (tx["status"] == "ok"), "tx"].to_numpy()
    seg_c = seg[seg["tx"].isin(ok_tx)]
    snv, bad_tx, qc = gate.enumerate_snvs(seg_c, pos_u, ref_u)
    am = gate.load_am(work, c)
    am = am[am["tx"].isin(ok_tx)]
    am = am.assign(key=am["pos"].to_numpy().astype(np.int64) * 4 + am["alt"].to_numpy().astype(np.int64))
    am_core = am[am["prot_pos"] > 1]
    mis = snv[snv["cls"] == 1]
    m = mis[["key", "tx", "gref", "prot_pos", "ref_aa", "alt_aa", "near"]].merge(
        am_core[["key", "tx", "ref", "prot_pos", "ref_aa", "alt_aa", "am"]], on=["key", "tx"], how="inner",
        suffixes=("", "_am"))
    match = ((m["gref"] == m["ref"]) & (m["prot_pos"] == m["prot_pos_am"])
             & (m["ref_aa"] == m["ref_aa_am"]) & (m["alt_aa"] == m["alt_aa_am"]))
    mm = m.loc[match, ["key", "tx", "gref", "prot_pos", "ref_aa", "alt_aa", "near", "am"]]
    q = pd.DataFrame({"n_enum": mis.groupby("tx").size(), "n_am": am_core.groupby("tx").size(),
                      "n_match": mm.groupby("tx").size()}).reindex(ok_tx).fillna(0).astype(np.int64)
    low = (q["n_match"] < gate.CONCORDANCE_MIN * q["n_am"]) | (q["n_match"] < gate.CONCORDANCE_MIN * q["n_enum"])
    bad = low | (q["n_am"] == 0) | q.index.isin(list(bad_tx))
    pass_tx = set(q.index[~bad].tolist())
    del am, am_core, mis, m, match
    in_pass = snv["tx"].isin(pass_tx)
    allk = snv.loc[in_pass & snv["cls"].isin([0, 1, 2]), "key"]
    amb_keys = set(allk[allk.duplicated(keep=False)].tolist())
    syn = snv.loc[in_pass & (snv["cls"] == 0), ["key", "tx", "gref", "prot_pos", "ref_aa", "alt_aa", "near"]].assign(cls=0, am=np.nan)
    misu = mm[mm["tx"].isin(pass_tx)].assign(cls=1)
    u = pd.concat([syn, misu], ignore_index=True)
    u["galt"] = (u["key"].to_numpy() % 4).astype(np.int8)
    u = u[~u["key"].isin(amb_keys)]
    u = u.merge(r[["key", "filt_ok", "pn"] + RATE_MODELS], on="key", how="left")
    del r, syn, misu, snv
    u["filt_ok"] = u["filt_ok"].fillna(False).astype(bool)
    if not isinstance(u["pn"].dtype, pd.CategoricalDtype):
        u["pn"] = u["pn"].astype("category")
    u["mutclass"], _ = gate.mutation_class(u["pn"], u["gref"], u["galt"])
    u["am_bin"] = np.where(u["cls"] == 1, gate.am_fine_bin(u["am"].fillna(0)), -1).astype(np.int16)
    out = pd.DataFrame({
        "key": u["key"].to_numpy(np.int64), "tx": u["tx"].to_numpy(np.int32), "cls": u["cls"].to_numpy(np.int8),
        "am": u["am"].to_numpy(np.float32), "am_bin": u["am_bin"].to_numpy(np.int16),
        "prot_pos": u["prot_pos"].to_numpy(np.int32), "ref_aa": u["ref_aa"].to_numpy(np.uint8),
        "alt_aa": u["alt_aa"].to_numpy(np.uint8), "mutclass": u["mutclass"].to_numpy(np.int8),
        "filt_ok": u["filt_ok"].to_numpy(bool), "near": u["near"].to_numpy(bool)})
    for mdl in RATE_MODELS:
        out[mdl] = u[mdl].to_numpy(np.float64)
    return out


def cmd_universe(a):
    tx = pd.read_pickle(P(a.work, "gtf", "transcripts.pkl"))
    seg = pd.read_pickle(P(a.work, "gtf", "segments.pkl"))
    od = P(a.fs, "universe")
    os.makedirs(od, exist_ok=True)
    for c in gate.parse_chroms(a.chroms):
        f = P(od, f"chr{c}.pkl")
        if os.path.exists(f) and not a.force:
            continue
        t0 = time.time()
        u = build_universe(a.work, c, tx, seg)
        # validation: re-aggregate exactly as the gate did and compare with its saved aggregates
        g = pd.read_pickle(P(a.work, "sites", f"chr{c}.agg.pkl"))["mis_bins"]
        mine = u[u["cls"] == 1].groupby(["am_bin", "mutclass", "filt_ok", "near"], observed=True).agg(
            n=("key", "size"), MR=("MR", "sum")).reset_index()
        cmp = g.merge(mine, on=["am_bin", "mutclass", "filt_ok", "near"], how="outer", suffixes=("_g", "_f"))
        ok = (len(g) == len(mine)) and np.array_equal(cmp["n_g"].to_numpy(), cmp["n_f"].to_numpy()) and \
            np.allclose(cmp["MR_g"].to_numpy(), cmp["MR_f"].to_numpy(), rtol=1e-12, atol=0)
        if not ok:
            die(f"chr{c}: universe does not reproduce the gate's aggregates")
        u.to_pickle(f)
        log(f"chr{c}: {int((u['cls'] == 1).sum()):,} missense + {int((u['cls'] == 0).sum()):,} synonymous sites; "
            f"identical to gate aggregates ({time.time() - t0:.0f}s)")
        del u
        gc.collect()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=P(HERE, "data"))
    ap.add_argument("--work", default=P(HERE, "work"))
    ap.add_argument("--fs", default=P(HERE, "work", "fullstudy"))
    ap.add_argument("--results", default=P(HERE, "results_fullstudy"))
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("universe")
    s.add_argument("--chroms", default="1-22")
    s.add_argument("--force", action="store_true")
    a = ap.parse_args()
    {"universe": cmd_universe}[a.cmd](a)


if __name__ == "__main__":
    main()
