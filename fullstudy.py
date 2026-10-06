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


# ======================================================================== shared loaders
TOOLS = {  # name -> dbNSFP 4.1a rankscore column (higher = more damaging); FULL_STUDY_PLAN.md section 5
    "SIFT": "SIFT_converted_rankscore", "SIFT4G": "SIFT4G_converted_rankscore",
    "PolyPhen2_HDIV": "Polyphen2_HDIV_rankscore", "PolyPhen2_HVAR": "Polyphen2_HVAR_rankscore",
    "LRT": "LRT_converted_rankscore", "MutationTaster": "MutationTaster_converted_rankscore",
    "MutationAssessor": "MutationAssessor_rankscore", "FATHMM": "FATHMM_converted_rankscore",
    "PROVEAN": "PROVEAN_converted_rankscore", "VEST4": "VEST4_rankscore", "MetaSVM": "MetaSVM_rankscore",
    "MetaLR": "MetaLR_rankscore", "M-CAP": "M-CAP_rankscore", "REVEL": "REVEL_rankscore",
    "MutPred": "MutPred_rankscore", "MVP": "MVP_rankscore", "MPC": "MPC_rankscore", "PrimateAI": "PrimateAI_rankscore",
    "DEOGEN2": "DEOGEN2_rankscore", "BayesDel_addAF": "BayesDel_addAF_rankscore",
    "BayesDel_noAF": "BayesDel_noAF_rankscore", "ClinPred": "ClinPred_rankscore", "LIST-S2": "LIST-S2_rankscore",
    "CADD": "CADD_raw_rankscore", "DANN": "DANN_rankscore", "fathmm-MKL": "fathmm-MKL_coding_rankscore",
    "fathmm-XF": "fathmm-XF_coding_rankscore", "Eigen": "Eigen-raw_coding_rankscore",
    "Eigen-PC": "Eigen-PC-raw_coding_rankscore", "GERP++": "GERP++_RS_rankscore",
    "phyloP100way": "phyloP100way_vertebrate_rankscore", "phyloP30way": "phyloP30way_mammalian_rankscore",
    "phastCons100way": "phastCons100way_vertebrate_rankscore", "SiPhy": "SiPhy_29way_logOdds_rankscore",
    "GenoCanyon": "GenoCanyon_rankscore", "fitCons": "integrated_fitCons_rankscore",
}


TOOLS54 = {  # amendment A3: dbNSFP 5.4a tools not analysed from 4.1a; stored in the column store as "v54_<column>"
    "MetaRNN": "MetaRNN_rankscore", "MutPred2": "MutPred2_rankscore", "gMVP": "gMVP_rankscore",
    "MisFit_D": "MisFit_D_rankscore", "MisFit_S": "MisFit_S_rankscore", "VARITY_R": "VARITY_R_rankscore",
    "VARITY_ER": "VARITY_ER_rankscore", "VARITY_R_LOO": "VARITY_R_LOO_rankscore",
    "VARITY_ER_LOO": "VARITY_ER_LOO_rankscore", "ESM1b": "ESM1b_converted_rankscore",
    "PHACTboost": "PHACTboost_rankscore", "MutFormer": "MutFormer_rankscore", "MutScore": "MutScore_rankscore",
    "popEVE": "popEVE_converted_rankscore", "GPN-MSA": "GPN_MSA_converted_rankscore",
    "GERP_92m": "GERP_92_mammals_rankscore", "phyloP470way": "phyloP470way_mammalian_rankscore",
    "phastCons470way": "phastCons470way_mammalian_rankscore",
}
DDD_ONLY = {"popEVE"}


def score_col(s):
    if s in M4_COLS:
        return M4_COLS[s]
    if s in M4V3:
        return M4V3[s]
    if s in M4V4:
        return M4V4[s]
    if s in M4V5:
        return M4V5[s]
    return ("v54_" + TOOLS54[s]) if s in TOOLS54 else TOOLS[s]            # section 6: developed/evaluated on DD de novo cohorts -> DDD de novos only
DB54 = P(EXT, "dbnsfp54", "dbNSFP5.4a_grch38.gz")


def cmd_scores54(a):
    """dbNSFP 5.4a rankscores for the A3 tools, written into the column store in its row order."""
    import subprocess
    if not os.path.exists(DB54 + ".ok"):
        die("dbNSFP 5.4a not downloaded/verified yet")
    with gzip.open(DB54, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    want = sorted(["pos(1-based)", "alt"] + list(TOOLS54.values()), key=header.index)
    idx = ",".join(str(header.index(w) + 1) for w in want)
    for c in gate.parse_chroms(a.chroms):
        od = P(a.fs, "col", f"chr{c}")
        if os.path.exists(P(od, "_v54_done")) and not a.force:
            continue
        t0 = time.time()
        keys = np.load(P(od, "key.npy"))
        tmp = P(a.fs, "col", f"chr{c}.v54.tsv")
        r = subprocess.run(f"tabix '{DB54}' {c} | cut -f{idx} > '{tmp}'", shell=True)
        if r.returncode != 0 or os.path.getsize(tmp) == 0:
            subprocess.run(f"tabix '{DB54}' chr{c} | cut -f{idx} > '{tmp}'", shell=True, check=True)
        got = {}
        for ch in pd.read_csv(tmp, sep="\t", header=None, names=want, na_values=["."], dtype=str, chunksize=1_000_000):
            alt = ch["alt"].map(B2I)
            ok = alt.notna()
            k = ch.loc[ok, "pos(1-based)"].astype(np.int64).to_numpy() * 4 + alt[ok].astype(np.int64).to_numpy()
            sel = np.isin(k, keys)
            sub = ch.loc[ok].loc[sel, list(TOOLS54.values())].apply(pd.to_numeric, errors="coerce").astype(np.float32)
            sub.insert(0, "key", k[sel])
            got.setdefault("parts", []).append(sub)
        os.remove(tmp)
        s = pd.concat(got["parts"], ignore_index=True).drop_duplicates("key").set_index("key").reindex(keys)
        for col in TOOLS54.values():
            np.save(P(od, "v54_" + col + ".npy"), s[col].to_numpy(np.float32))
        open(P(od, "_v54_done"), "w").close()
        log(f"chr{c}: dbNSFP 5.4a A3 tools joined at {s.notna().any(axis=1).mean():.3f} of column-store sites "
            f"({time.time() - t0:.0f}s)")


def cmd_scores54_stream(a):
    """Engineering alternative to `scores54` (deviation log): stream the downloaded byte-range pieces in order
    (their concatenation IS the BGZF file), decompress, keep the A3 columns at column-store sites, and md5 the
    byte stream in the same pass. Columns are committed only if the md5 equals the published checksum."""
    import subprocess
    d54 = P(EXT, "dbnsfp54")
    ranges = pd.read_csv(P(d54, "ranges.tsv"), sep="\t", header=None, names=["seg", "sub", "file", "a", "b"])
    ranges = ranges.sort_values("a")
    for r in ranges.itertuples():
        if os.path.getsize(P(d54, r.file)) != r.b - r.a + 1:
            die(f"piece {r.file} incomplete")
    files = " ".join(f"'{P(d54, f)}'" for f in ranges["file"])
    with gzip.open(P(d54, ranges["file"].iloc[0]), "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    want = sorted(["#chr", "pos(1-based)", "alt"] + list(TOOLS54.values()), key=header.index)
    idx = ",".join(str(header.index(w) + 1) for w in want)
    md5out = P(d54, "stream.md5")
    cmd = f"cat {files} | tee >(md5sum > '{md5out}') | bgzip -dc -@ 4 | cut -f{idx}"
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE)
    t0 = time.time()
    keys_of, buf, done_c = {}, {}, []

    def finish(c):
        if c not in AUTOSOMES or c not in buf:
            return
        od = P(a.fs, "col", f"chr{c}")
        s = pd.concat(buf.pop(c), ignore_index=True).drop_duplicates("key").set_index("key").reindex(keys_of[c])
        for col in TOOLS54.values():
            np.save(P(od, "v54tmp_" + col + ".npy"), s[col].to_numpy(np.float32))
        done_c.append(c)
        log(f"chr{c}: 5.4a A3 columns extracted at {s.notna().any(axis=1).mean():.3f} of sites ({time.time() - t0:.0f}s)")

    cur = None
    only = set(gate.parse_chroms(a.chroms))      # re-extract a subset (the md5 still covers the whole stream)
    for ch in pd.read_csv(proc.stdout, sep="\t", header=0, names=want, dtype=str, na_values=["."], chunksize=500_000):
        ch = ch[ch["#chr"].isin(only)]
        for c, g in ch.groupby("#chr", sort=False):
            if cur is not None and c != cur:
                finish(cur)
            if c in done_c:
                die(f"chromosome {c} is not contiguous in the stream - aborting before anything is committed")
            cur = c
            if c not in keys_of:
                keys_of[c] = np.load(P(a.fs, "col", f"chr{c}", "key.npy"))
            alt = g["alt"].map(B2I)
            ok = alt.notna()
            k = g.loc[ok, "pos(1-based)"].astype(np.int64).to_numpy() * 4 + alt[ok].astype(np.int64).to_numpy()
            sel = np.isin(k, keys_of[c])
            sub = g.loc[ok].loc[sel, list(TOOLS54.values())].apply(pd.to_numeric, errors="coerce").astype(np.float32)
            sub.insert(0, "key", k[sel])
            buf.setdefault(c, []).append(sub)
    finish(cur)
    proc.wait()
    time.sleep(5)
    got = open(md5out).read().split()[0] if os.path.exists(md5out) else ""
    expect = open(P(d54, "dbNSFP5.4a_grch38.gz.md5")).read().split()[0]
    if got != expect or proc.returncode != 0:
        die(f"md5 mismatch or stream error (got {got}, expected {expect}, rc {proc.returncode}) - columns NOT committed")
    for c in done_c:
        od = P(a.fs, "col", f"chr{c}")
        for col in TOOLS54.values():
            os.replace(P(od, "v54tmp_" + col + ".npy"), P(od, "v54_" + col + ".npy"))
        open(P(od, "_v54_done"), "w").close()
    open(P(d54, "dbNSFP5.4a_grch38.gz.ok"), "w").write("verified by streaming md5 " + got + "\n")
    log(f"dbNSFP 5.4a: md5 verified ({got}); A3 columns committed for {len(done_c)} chromosomes ({time.time() - t0:.0f}s)")


def cmd_dom54_stream(a):
    """Amendment A4.2: InterPro domain of every column-store site from the dbNSFP 5.4a stream (md5 re-checked in the
    same pass), then dense unit ids: dom = gene x first-listed domain ("none" pooled per gene); dom_rmc = dom x RMC unit."""
    import subprocess
    d54 = P(EXT, "dbnsfp54")
    ranges = pd.read_csv(P(d54, "ranges.tsv"), sep="\t", header=None, names=["seg", "sub", "file", "a", "b"]).sort_values("a")
    files = " ".join(f"'{P(d54, f)}'" for f in ranges["file"])
    with gzip.open(P(d54, ranges["file"].iloc[0]), "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    want = sorted(["#chr", "pos(1-based)", "alt", "Interpro_domain"], key=header.index)
    idx = ",".join(str(header.index(w) + 1) for w in want)
    md5out = P(d54, "stream_dom.md5")
    proc = subprocess.Popen(["bash", "-c", f"cat {files} | tee >(md5sum > '{md5out}') | bgzip -dc -@ 4 | cut -f{idx}"],
                            stdout=subprocess.PIPE)
    t0 = time.time()
    vocab = {"none": 0}
    keys_of, buf, done_c = {}, {}, []

    def first_dom(s):
        if not isinstance(s, str):
            return 0
        for t in s.split(";"):
            t = t.strip()
            if t and t != ".":
                return vocab.setdefault(t, len(vocab))
        return 0

    def finish(c):
        if c not in AUTOSOMES or c not in buf:
            return
        s = pd.concat(buf.pop(c), ignore_index=True).drop_duplicates("key").set_index("key")["dom"].reindex(keys_of[c])
        np.save(P(a.fs, "col", f"chr{c}", "domtmp_code.npy"), s.fillna(0).to_numpy(np.int32))
        done_c.append(c)
        log(f"chr{c}: InterPro domain at {(s > 0).mean():.3f} of sites ({time.time() - t0:.0f}s)")

    cur = None
    for ch in pd.read_csv(proc.stdout, sep="\t", header=0, names=want, dtype=str, chunksize=500_000):
        for c, g in ch.groupby("#chr", sort=False):
            if cur is not None and c != cur:
                finish(cur)
            if c in done_c:
                die(f"chromosome {c} is not contiguous in the stream")
            cur = c
            if c not in AUTOSOMES:
                continue
            if c not in keys_of:
                keys_of[c] = np.load(P(a.fs, "col", f"chr{c}", "key.npy"))
            alt = g["alt"].map(B2I)
            ok = alt.notna()
            k = g.loc[ok, "pos(1-based)"].astype(np.int64).to_numpy() * 4 + alt[ok].astype(np.int64).to_numpy()
            sel = np.isin(k, keys_of[c])
            buf.setdefault(c, []).append(pd.DataFrame({"key": k[sel], "dom": g.loc[ok, "Interpro_domain"].to_numpy()[sel]})
                                         .assign(dom=lambda x: x["dom"].map(first_dom)))
    finish(cur)
    proc.wait()
    time.sleep(5)
    got = open(md5out).read().split()[0] if os.path.exists(md5out) else ""
    expect = open(P(d54, "dbNSFP5.4a_grch38.gz.md5")).read().split()[0]
    if got != expect or proc.returncode != 0 or len(done_c) != 22:
        die(f"md5 mismatch, stream error or missing chromosomes (got {got}, {len(done_c)} chromosomes) - not committed")
    pd.Series(list(vocab), name="interpro_domain").to_csv(P(a.fs, "regions", "dom_vocab.tsv"), sep="\t")
    n_dom, n_dr = build_dom_units(a, AUTOSOMES)
    log(f"InterPro: md5 verified ({got}); {len(vocab) - 1:,} domain names; {n_dom:,} gene x domain units; "
        f"{n_dr:,} domain x RMC units ({time.time() - t0:.0f}s)")


def build_dom_units(a, chroms):
    """Dense unit ids across chromosomes (transcripts never span chromosomes) from each chromosome's domtmp_code."""
    n_dom = n_dr = 0
    tx_dom, tx_dr = [], []
    for c in chroms:
        od = P(a.fs, "col", f"chr{c}")
        code = np.load(P(od, "domtmp_code.npy")).astype(np.int64)
        tx = np.load(P(od, "tx.npy")).astype(np.int64)
        ru = np.load(P(od, "rmc_unit.npy")).astype(np.int64)
        pair = tx * 1_000_000 + code
        uq, inv = np.unique(pair, return_inverse=True)
        dom_unit = inv + n_dom
        tx_dom.append(uq // 1_000_000)
        n_dom += len(uq)
        ok = ru >= 0
        dr = np.full(len(tx), -1, np.int64)
        if ok.any():
            uq2, inv2 = np.unique(np.c_[dom_unit[ok], ru[ok]], axis=0, return_inverse=True)
            dr[ok] = inv2.ravel() + n_dr
            tx_dr.append(np.asarray(tx_dom[-1])[uq2[:, 0] - (n_dom - len(uq))])
            n_dr += len(uq2)
        np.save(P(od, "dom_unit.npy"), dom_unit)
        np.save(P(od, "domrmc_unit.npy"), dr)
        os.replace(P(od, "domtmp_code.npy"), P(od, "dom_code.npy"))
    np.savez(P(a.fs, "regions", "dom_units.npz"), n_dom=n_dom - 1, n_domrmc=n_dr - 1,
             tx_dom=np.concatenate(tx_dom), tx_domrmc=np.concatenate(tx_dr))
    return n_dom, n_dr


def cmd_scores54_tail(a):
    """Engineering (deviation log): re-extract A3 columns for a few late chromosomes by streaming only the trailing pieces
    of the (already md5-verified) download, starting at the first BGZF block boundary of piece --from-piece. Every BGZF
    block's CRC is still checked by bgzip. Requires that the stream begins inside an earlier chromosome."""
    import subprocess
    d54 = P(EXT, "dbnsfp54")
    ranges = pd.read_csv(P(d54, "ranges.tsv"), sep="\t", header=None, names=["seg", "sub", "file", "a", "b"]).sort_values("a")
    files = list(ranges["file"])
    first = files.index(a.from_piece)
    magic = bytes.fromhex("1f8b08040000000000ff0600424302")
    with open(P(d54, files[first]), "rb") as fh:
        head = fh.read(1 << 20)
    skip = head.find(magic)
    if skip < 0:
        die("no BGZF block start in the first MiB of the piece")
    with gzip.open(P(d54, files[0]), "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    want = sorted(["#chr", "pos(1-based)", "alt"] + list(TOOLS54.values()), key=header.index)
    idx = ",".join(str(header.index(w) + 1) for w in want)
    rest = " ".join(f"'{P(d54, f)}'" for f in files[first + 1:])
    cmd = f"{{ tail -c +{skip + 1} '{P(d54, files[first])}'; cat {rest}; }} | bgzip -dc -@ 4 | tail -n +2 | cut -f{idx}"
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE)
    only = [c for c in gate.parse_chroms(a.chroms)]
    seen_before, buf, keys_of = set(), {c: [] for c in only}, {}
    for ch in pd.read_csv(proc.stdout, sep="\t", header=None, names=want, dtype=str, na_values=["."], chunksize=500_000):
        for c, g in ch.groupby("#chr", sort=False):
            if c not in only:
                if not any(len(v) for v in buf.values()):
                    seen_before.add(c)
                continue
            if not seen_before:
                die(f"stream started inside chr{c}; choose an earlier --from-piece")
            if c not in keys_of:
                keys_of[c] = np.load(P(a.fs, "col", f"chr{c}", "key.npy"))
            alt = g["alt"].map(B2I); ok = alt.notna()
            k = g.loc[ok, "pos(1-based)"].astype(np.int64).to_numpy() * 4 + alt[ok].astype(np.int64).to_numpy()
            sel = np.isin(k, keys_of[c])
            sub = g.loc[ok].loc[sel, list(TOOLS54.values())].apply(pd.to_numeric, errors="coerce").astype(np.float32)
            sub.insert(0, "key", k[sel])
            buf[c].append(sub)
    proc.wait()
    if proc.returncode != 0:
        die(f"stream error rc {proc.returncode}")
    for c in only:
        od = P(a.fs, "col", f"chr{c}")
        s = pd.concat(buf[c], ignore_index=True).drop_duplicates("key").set_index("key").reindex(keys_of[c])
        for col in TOOLS54.values():
            np.save(P(od, "v54_" + col + ".npy"), s[col].to_numpy(np.float32))
        open(P(od, "_v54_done"), "w").close()
        log(f"chr{c}: A3 columns re-extracted from the stream tail at {s.notna().any(axis=1).mean():.3f} of sites "
            f"(stream began in chr{sorted(seen_before)})")


def load_dnv(a):
    """Annotated de novos (gate's per-chromosome tables) + proband index identical to the gate's bootstrap order."""
    meta = pd.read_csv(P(a.work, "dnv_grch38.tsv"), sep="\t", dtype={"chrom": str, "id": str, "cq": str})
    dv = pd.concat([pd.read_pickle(P(a.work, "sites", f"chr{c}.dnv.pkl")) for c in AUTOSOMES], ignore_index=True)
    d = dv.merge(meta[["dnv_idx", "id", "centre", "symbol", "chrom", "pos38", "alt38"]], on="dnv_idx", how="left")
    meta["pkey"] = meta["centre"] + "|" + meta["id"].astype(str)
    pk = pd.Series(np.arange(meta["pkey"].nunique()), index=pd.unique(meta["pkey"]))
    d["p"] = pk.reindex(d["centre"] + "|" + d["id"].astype(str)).to_numpy()
    d = d[d["status"].isin(["missense", "synonymous"])].copy()
    d["key"] = d["pos38"].to_numpy(np.int64) * 4 + d["alt38"].map(B2I).to_numpy(np.int64)
    d["tx"] = d["tx"].astype(np.int64)
    return d, len(pk)


def load_tinfo(a):
    """Per transcript: gene, Ensembl gene id, s_het class, DDG2P monoallelic flag, mechanism labels (Aim 2)."""
    tx = pd.read_pickle(P(a.work, "gtf", "transcripts.pkl")).set_index("tx")
    genes = pd.read_csv(P(a.data, gate.GENE_FILE), sep="\t", low_memory=False, usecols=["symbol", "s_het"])
    shet = genes.dropna(subset=["s_het"]).drop_duplicates("symbol").set_index("symbol")["s_het"]
    dd = pd.read_csv(P(a.data, gate.DDG2P_FILE), dtype=str)
    qual = (dd["DDD category"].str.lower().isin(["confirmed", "probable"])
            & dd["allelic requirement"].str.lower().str.contains("monoallelic", na=False))
    dd_syms = set(dd.loc[qual, "gene symbol"].str.upper())
    for s in dd.loc[qual, "prev symbols"].dropna():
        dd_syms.update(x.strip().upper() for x in s.split(";") if x.strip())
    t = pd.DataFrame(index=tx.index)
    t["gene"] = tx["gene_name"]
    t["ensg"] = tx["gene_id"].str.split(".").str[0]
    t["chrom"] = tx["chrom"]
    t["s_het"] = t["gene"].map(shet)
    t["shet_class"] = np.where(t["s_het"].isna(), 2, (t["s_het"] >= gate.SHET_HIGH).astype(int))
    t["ddg2p"] = t["gene"].str.upper().isin(dd_syms)
    # ---- Aim 2 mechanism labels, precedence 1 (B&M S3 literature) -> 2 (DDG2P mutation consequence) -> 3 (B&M S5 verdict)
    lab = {}
    s3 = pd.read_csv(P(EXT, "aims", "bm_s003_labels.csv"), dtype=str, encoding="latin-1", usecols=["gene", "class"])
    s3["gene"] = s3["gene"].str.upper()
    multi = s3.groupby("gene")["class"].nunique()
    for g, cl in s3.drop_duplicates("gene").set_index("gene")["class"].items():
        lab[g] = ("conflict" if multi[g] > 1 else cl, "BM_S3")
    ddm = dd.loc[qual, ["gene symbol", "mutation consequence"]].dropna()
    for g, mc in zip(ddm["gene symbol"].str.upper(), ddm["mutation consequence"].str.lower()):
        if g in lab:
            continue
        cl = "lof" if mc == "loss of function" else ("nonlof" if mc in ("activating", "dominant negative") else None)
        if cl:
            lab[g] = (cl, "DDG2P")
    t["mech12"] = [lab.get(g.upper(), (None, None))[0] for g in t["gene"]]
    t["mech_src"] = [lab.get(g.upper(), (None, None))[1] for g in t["gene"]]
    s5 = pd.read_excel(P(EXT, "aims", "bm_s005.xlsx"), sheet_name="table_S3")
    s5 = s5[s5["verdict"].isin(["dn", "gof", "lof"])].drop_duplicates("gene").set_index("gene")["verdict"]
    t["mech123"] = t["mech12"]
    miss = t["mech123"].isna()
    t.loc[miss, "mech123"] = t.loc[miss, "gene"].map(s5)
    for col in ("mech12", "mech123"):
        t[col + "_bin"] = t[col].map({"lof": "LoF", "dn": "nonLoF", "gof": "nonLoF", "nonlof": "nonLoF"})
    return t


COL_BASE = ["key", "tx", "am", "am_bin", "near", "ref_aa", "alt_aa", "mutclass", "prot_pos", PRIMARY_MODEL]


def cmd_colstore(a):
    """Engineering: the primary-filter missense sites with their joined scores and regions, one .npy per column per
    chromosome (row order identical to site_chunks' pickle path), so a pass reads only the columns it uses."""
    for c in gate.parse_chroms(a.chroms):
        od = P(a.fs, "col", f"chr{c}")
        if os.path.exists(P(od, "_done")) and not a.force:
            continue
        t0 = time.time()
        os.makedirs(od, exist_ok=True)
        for _, u in _pickle_chunks(a, [c], "exclude", 1, ("scores", "regions")):
            for col in u.columns:
                np.save(P(od, col.replace("/", "_") + ".npy"), u[col].to_numpy())
            n = len(u)
        open(P(od, "_done"), "w").write(str(n))
        log(f"chr{c}: column store written ({n:,} sites, {time.time() - t0:.0f}s)")


def _pickle_chunks(a, chroms, near, cls, extra):
    for c in chroms:
        u = pd.read_pickle(P(a.fs, "universe", f"chr{c}.pkl"))
        m = (u["cls"] == cls) & u["filt_ok"] & (u[PRIMARY_MODEL] > 0)
        if near == "exclude":
            m &= ~u["near"]
        u = u.loc[m, COL_BASE]
        for name in extra:
            u = u.merge(pd.read_pickle(P(a.fs, name, f"chr{c}.pkl")), on="key", how="left")
        yield c, u
        del u
        gc.collect()


def site_chunks(a, chroms=AUTOSOMES, near="exclude", cls=1, extra=(), cols=None):
    """Per chromosome, universe sites of class `cls` passing the primary filters, with optional joined tables.
    Reads the column store (only `cols`) when it exists for the primary missense set; otherwise the pickles."""
    for c in chroms:
        od = P(a.fs, "col", f"chr{c}")
        if cls == 1 and near == "exclude" and os.path.exists(P(od, "_done")):
            names = cols if cols is not None else [f[:-4] for f in sorted(os.listdir(od)) if f.endswith(".npy")]
            u = pd.DataFrame({n: np.load(P(od, n.replace("/", "_") + ".npy")) for n in dict.fromkeys(names)})
            yield c, u
            del u
            gc.collect()
            continue
        u = pd.read_pickle(P(a.fs, "universe", f"chr{c}.pkl"))
        m = (u["cls"] == cls) & u["filt_ok"] & (u[PRIMARY_MODEL] > 0)
        if near == "exclude":
            m &= ~u["near"]
        u = u.loc[m, ["key", "tx", "am", "am_bin", "near", "ref_aa", "alt_aa", "mutclass", "prot_pos", PRIMARY_MODEL]]
        for name in extra:
            u = u.merge(pd.read_pickle(P(a.fs, name, f"chr{c}.pkl")), on="key", how="left")
        yield c, u
        del u
        gc.collect()


def dnv_filter(d, near="exclude"):
    m = d["filt_ok"].fillna(False).astype(bool) & (d[PRIMARY_MODEL].to_numpy(dtype=float) > 0)
    if near == "exclude":
        m &= ~d["near"].fillna(False).astype(bool)
    return d[m].copy()


# ======================================================================== engine
from scipy import sparse  # noqa: E402


class Engine:
    """Gate estimators generalised to arbitrary units (the thing held fixed) and arbitrary score bins.

    run(Mass, syn_mass, syn_p, dnv): Mass = sparse (units x bins) sum of MR over analysed missense sites;
    syn_mass = total synonymous MR; syn_p = proband index of each synonymous de novo; dnv = missense de novos
    with columns p, unit (row of Mass; -1 = outside the analysis) and bin (-1 = unscored)."""

    def __init__(self, n_prob, seed=SEED, n_boot=N_BOOT):
        self.boot = gate.Boot(n_prob, n_boot, seed)

    def counts(self, p, col, ncol, vals=None):
        return self.boot.counts(np.asarray(p, dtype=np.int64), None if col is None else np.asarray(col, dtype=np.int64),
                                ncol, vals)

    def run(self, Mass, syn_mass, syn_p, dnv, unit_mask=None):
        Mass = sparse.csr_matrix(Mass)
        nu, nb = Mass.shape
        if unit_mask is not None:
            Mass = sparse.diags(np.asarray(unit_mask, dtype=float)) @ Mass
            Mass = sparse.csr_matrix(Mass)
        c = self.counts(syn_p, np.zeros(len(syn_p)), 1)[:, 0] / syn_mass          # (B+1,) pooled scale factor
        Mtot = np.asarray(Mass.sum(0)).ravel()
        E = c[:, None] * Mtot[None, :]
        rs = np.asarray(Mass.sum(1)).ravel()
        dd = dnv[(dnv["unit"].to_numpy() >= 0) & (dnv["bin"].to_numpy() >= 0)]
        uu = dd["unit"].to_numpy(np.int64)
        k2 = rs[uu] > 0
        dd, uu = dd[k2], uu[k2]
        O = self.counts(dd["p"], dd["bin"], nb)
        # within-unit reweighting: each de novo contributes its unit's opportunity profile across bins
        Cp = np.zeros((self.boot.n, nb))
        pp = dd["p"].to_numpy(np.int64)
        for s0 in range(0, len(uu), 2000):
            rows = Mass[uu[s0:s0 + 2000]].toarray() / rs[uu[s0:s0 + 2000]][:, None]
            np.add.at(Cp, pp[s0:s0 + 2000], rows)
        Operm = self.boot.W @ Cp
        return {"O": O, "E": E, "Operm": Operm, "c": c, "n_dnv": int(len(dd))}

    @staticmethod
    def lr(res, groups=None):
        """LRs over groups of bins (list of index arrays; None = each bin). Arrays are (B+1, G)."""
        O, E, Op = res["O"], res["E"], res["Operm"]
        Xt = (O - E).sum(1, keepdims=True)
        Et = E.sum(1, keepdims=True)
        if groups is not None:
            O = np.stack([O[:, g].sum(1) for g in groups], 1)
            E = np.stack([E[:, g].sum(1) for g in groups], 1)
            Op = np.stack([Op[:, g].sum(1) for g in groups], 1)
        out = {"O": O, "E": E, "X": O - E}
        with np.errstate(invalid="ignore", divide="ignore"):
            fB = E / Et
            for tag, k in (("", 1.0), ("_dn", 1 + DELTA), ("_up", 1 - DELTA)):
                X = O - k * E
                out["LRm" + tag] = (X / Xt) / fB
                out["LRw" + tag] = X / (Op - E)
        return out


def summ(arr):
    return gate.summary(arr)


def points(lr):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.log(lr) / np.log(PP3[0][1])


def window_groups(nb=N_QBINS, half=WIN_HALF, step=WIN_STEP):
    centres = np.arange(half, nb - half, step)
    return centres, [np.arange(cc - half, cc + half + 1) for cc in centres]


def window_thresholds(lo_sys, hi_sys, centre_scores):
    """Pejaver monotone rule over windows ordered by score (FULL_STUDY_PLAN.md section 2)."""
    out = {}
    n = len(lo_sys)
    for name, cut in PP3:
        thr = None
        for i in range(n):
            if np.all(np.nan_to_num(lo_sys[i:], nan=-np.inf) >= cut):
                thr = float(centre_scores[i])
                break
        out[f"PP3_{name}"] = thr
    for name, cut in PP3:
        thr = None
        for i in range(n - 1, -1, -1):
            if np.all(np.nan_to_num(hi_sys[: i + 1], nan=np.inf) <= 1.0 / cut):
                thr = float(centre_scores[i])
                break
        out[f"BP4_{name}"] = thr
    return out


def weighted_qedges(pairs, nq=N_QBINS, grid=200001):
    """Quantile edges of a [0,1] score weighted by MR, from a fine histogram over (score, weight) chunks."""
    h = np.zeros(grid)
    for s, w in pairs:
        s = np.asarray(s, dtype=np.float64)
        ok = np.isfinite(s)
        idx = np.clip(np.round(s[ok] * (grid - 1)).astype(np.int64), 0, grid - 1)
        h += np.bincount(idx, weights=np.asarray(w)[ok], minlength=grid)
    cdf = np.cumsum(h) / h.sum()
    xs = np.linspace(0, 1, grid)
    qs = np.linspace(0, 1, nq + 1)[1:-1]
    inner = xs[np.minimum(np.searchsorted(cdf, qs), grid - 1)]
    return np.concatenate([[-np.inf], inner, [np.inf]])


def to_bins(s, edges):
    s = np.asarray(s, dtype=np.float64)
    b = np.searchsorted(edges, s, side="right") - 1
    return np.where(np.isfinite(s), np.clip(b, 0, len(edges) - 2), -1).astype(np.int64)


class BinMap:
    """Fine score bins = quantiles of mutational opportunity, robust to tied scores (bug fix, deviation log).

    Each score cell of a 200,001-point grid on [0, 1] is assigned the bin floor(1000 * F(<cell)), where F(<cell) is the
    opportunity strictly below the cell: a block of tied scores always lands, whole, in the bin where it starts, and
    bins it spans beyond are left empty. lo/hi = the score range actually present in each bin."""

    def __init__(self, h, nq=N_QBINS):
        grid = len(h)
        tot = h.sum()
        below = (np.cumsum(h) - h) / tot
        self.grid = grid
        self.map = np.minimum((below * nq).astype(np.int64), nq - 1)
        xs = np.linspace(0, 1, grid)
        filled = h > 0
        self.lo = np.full(nq, np.nan)
        self.hi = np.full(nq, np.nan)
        b = self.map[filled]
        np.fmin.at(self.lo, b, xs[filled])
        np.fmax.at(self.hi, b, xs[filled])
        self.mass = np.bincount(self.map, weights=h, minlength=nq) / tot
        self.top_start = int(self.map[np.flatnonzero(filled)[-1]])     # bin holding the highest score present

    def bin(self, s):
        s = np.asarray(s, dtype=np.float64)
        ok = np.isfinite(s)
        idx = np.clip(np.round(np.where(ok, s, 0) * (self.grid - 1)).astype(np.int64), 0, self.grid - 1)
        return np.where(ok, self.map[idx], -1).astype(np.int64)

    def centres(self):
        c = (self.lo + self.hi) / 2
        return pd.Series(c).ffill().bfill().to_numpy()

    def top_group(self, frac=0.02):
        """Smallest set of top bins holding >= frac of opportunity (= bins 980-999 without ties; with a large tied
        block at the top it extends down to the start of that block). Returns (bins, actual opportunity share)."""
        tail = np.cumsum(self.mass[::-1])[::-1]                 # opportunity in bins >= b
        start = int(np.flatnonzero(tail >= frac - 0.5 / N_QBINS).max())     # half-bin tolerance for discreteness
        g = np.arange(start, N_QBINS)
        return g, float(self.mass[g].sum())


def bin_centre_scores(edges):
    inner = edges[1:-1]
    lo = np.concatenate([[inner[0]], inner])
    hi = np.concatenate([inner, [inner[-1]]])
    return (lo + hi) / 2


# ======================================================================== data joins
def cmd_scores(a):
    """dbNSFP 4.1a rankscores at universe missense sites (key = pos38*4 + alt)."""
    od = P(a.fs, "scores")
    os.makedirs(od, exist_ok=True)
    cols = list(TOOLS.values())
    raw = list(RAW_COLS)
    for c in gate.parse_chroms(a.chroms):
        f = P(od, f"chr{c}.pkl")
        if os.path.exists(f) and not a.force:
            continue
        src = P(EXT, "dbnsfp", f"dbNSFP4.1a_variant.chr{c}.gz")
        if not os.path.exists(src + ".ok"):
            log(f"chr{c}: dbNSFP file not downloaded/verified yet - skipped")
            continue
        t0 = time.time()
        u = pd.read_pickle(P(a.fs, "universe", f"chr{c}.pkl"))
        keys = u.loc[u["cls"] == 1, "key"].to_numpy()
        del u
        # narrow the 180-column file with zcat | cut (fast) before pandas sees it
        with gzip.open(src, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
        want = sorted(["pos(1-based)", "alt"] + cols + raw, key=header.index)     # cut -f emits columns in file order
        idx = [header.index(w) + 1 for w in want]
        narrow = P(a.fs, "scores", f"chr{c}.narrow.tsv")
        import subprocess
        subprocess.run(f"zcat '{src}' | cut -f{','.join(map(str, idx))} > '{narrow}'", shell=True, check=True)
        parts = []
        for ch in pd.read_csv(narrow, sep="\t", header=0, names=want, na_values=["."],
                              dtype={c_: str for c_ in ["alt"] + raw}, chunksize=1_000_000, low_memory=False):
            alt = ch["alt"].map(B2I)
            ok = alt.notna() & ch["pos(1-based)"].notna()
            ch = ch[ok]
            k = ch["pos(1-based)"].to_numpy(np.int64) * 4 + alt[ok].to_numpy(np.int64)
            sel = np.isin(k, keys)
            ch = ch[sel]
            sub = ch[cols].apply(pd.to_numeric, errors="coerce").astype(np.float32)
            for rc, how in RAW_COLS.items():      # raw scores, most damaging across transcripts (as dbNSFP's rankscore)
                sv = ch[rc]
                val = pd.to_numeric(sv, errors="coerce")
                multi = sv.str.contains(";", na=False, regex=False)
                if multi.any():
                    v = sv[multi].str.split(";", expand=True).apply(pd.to_numeric, errors="coerce")
                    val[multi] = v.min(axis=1) if how == "min" else v.max(axis=1)
                sub["raw_" + rc] = val.astype(np.float32).to_numpy()
            sub.insert(0, "key", k[sel])
            parts.append(sub)
        os.remove(narrow)
        s = pd.concat(parts, ignore_index=True).drop_duplicates("key")
        s.to_pickle(f)
        cov = s[cols].notna().mean()
        log(f"chr{c}: dbNSFP rows at {len(s):,}/{len(keys):,} universe missense sites; per-tool coverage "
            f"{cov.min():.2f}-{cov.max():.2f} ({time.time() - t0:.0f}s)")


RAW_COLS = {  # raw scores kept for converting thresholds to Pejaver 2022 units; "min" = smaller is more damaging
    "REVEL_score": "max", "BayesDel_addAF_score": "max", "CADD_phred": "max", "PrimateAI_score": "max",
    "VEST4_score": "max", "MPC_score": "max", "Polyphen2_HVAR_score": "max", "Polyphen2_HDIV_score": "max",
    "SIFT_score": "min", "FATHMM_score": "min", "GERP++_RS": "max", "phyloP100way_vertebrate": "max", "MutPred_score": "max",
}
CCR_EDGES = [0, 50, 80, 90, 95, 99, 100.0001]       # FULL_STUDY_PLAN.md section 3


def cmd_regions(a):
    """gnomAD v2.1.1 RMC regions and CCR intervals lifted to GRCh38, assigned to universe missense sites."""
    from pyliftover import LiftOver
    lo = LiftOver(P(a.data, gate.CHAIN_FILE))
    od = P(a.fs, "regions")
    os.makedirs(od, exist_ok=True)
    tinfo = load_tinfo(a)

    def lift(ch, pos0):
        r = lo.convert_coordinate("chr" + str(ch), int(pos0))
        return (r[0][0].replace("chr", ""), r[0][1]) if r else (None, None)

    # ---- RMC: regions with coordinates (transcripts with RMC) + whole-gene regions (transcripts without RMC)
    rmc_f = P(od, "rmc_regions.pkl")
    if not os.path.exists(rmc_f) or a.force:
        rw = pd.read_csv(P(EXT, "aims", "gnomAD_v2.1.1_transcripts_with_rmc.tsv"), sep="\t")
        rows = []
        for r in rw.itertuples(index=False):
            c1, p1 = str(r.start_coordinate).split(":")
            c2, p2 = str(r.stop_coordinate).split(":")
            a1, b1 = lift(c1, int(p1) - 1)
            a2, b2 = lift(c2, int(p2) - 1)
            ok = a1 is not None and a2 is not None and a1 == a2 == c1 and b2 >= b1 and \
                abs((b2 - b1) - (int(p2) - int(p1))) <= 0.1 * max(1, int(p2) - int(p1))
            rows.append((r.gene_id, r.gene_name, c1, (b1 + 1) if ok else -1, (b2 + 1) if ok else -1, r.oe, ok))
        reg = pd.DataFrame(rows, columns=["ensg", "gene", "chrom", "start38", "end38", "oe", "lifted"])
        wo = pd.read_csv(P(EXT, "aims", "gnomAD_v2.1.1_transcripts_without_rmc.tsv"), sep="\t")
        whole = pd.DataFrame({"ensg": wo["gene_id"], "gene": wo["gene_name"], "chrom": None, "start38": -1, "end38": -1,
                              "oe": wo["oe"], "lifted": True, "whole_gene": True})
        reg["whole_gene"] = False
        reg = pd.concat([reg, whole], ignore_index=True)
        reg["region_id"] = np.arange(len(reg))
        reg.to_pickle(rmc_f)
        log(f"RMC: {int((~reg['whole_gene']).sum()):,} regions in {rw['transcript'].nunique():,} transcripts "
            f"({int((~reg['whole_gene'] & ~reg['lifted']).sum())} failed liftover) + {int(reg['whole_gene'].sum()):,} "
            f"whole-gene units")
    reg = pd.read_pickle(rmc_f)

    # ---- CCR: lift every interval's endpoints
    ccr_f = P(od, "ccr_grch38.pkl")
    if not os.path.exists(ccr_f) or a.force:
        t0 = time.time()
        parts = []
        for ch in pd.read_csv(P(EXT, "aims", "ccrs.autosomes.v2.20180420.bed.gz"), sep="\t",
                              usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "ccr_pct"], header=0,
                              dtype={"chrom": str}, chunksize=1_000_000):
            s38 = [lift(c, s) for c, s in zip(ch["chrom"], ch["start"])]
            e38 = [lift(c, e - 1) for c, e in zip(ch["chrom"], ch["end"])]
            sc = np.array([x[0] for x in s38], dtype=object)
            ec = np.array([x[0] for x in e38], dtype=object)
            sp = np.array([x[1] if x[1] is not None else -1 for x in s38], dtype=np.int64)
            ep = np.array([x[1] if x[1] is not None else -1 for x in e38], dtype=np.int64)
            ok = (sc == ch["chrom"].to_numpy()) & (ec == sc) & (ep - sp == (ch["end"] - 1 - ch["start"]).to_numpy())
            parts.append(pd.DataFrame({"chrom": ch["chrom"].to_numpy()[ok], "s1": sp[ok] + 1, "e1": ep[ok] + 1,
                                       "ccr_pct": ch["ccr_pct"].to_numpy(np.float32)[ok]}))
        ccr = pd.concat(parts, ignore_index=True)
        ccr.to_pickle(ccr_f)
        log(f"CCR: {len(ccr):,} intervals lifted to GRCh38 ({time.time() - t0:.0f}s)")
    ccr = pd.read_pickle(ccr_f)

    for c in gate.parse_chroms(a.chroms):
        f = P(od, f"chr{c}.pkl")
        if os.path.exists(f) and not a.force:
            continue
        u = pd.read_pickle(P(a.fs, "universe", f"chr{c}.pkl"))
        u = u.loc[u["cls"] == 1, ["key", "tx"]].copy()
        u["pos"] = u["key"] // 4
        u["ensg"] = tinfo["ensg"].reindex(u["tx"]).to_numpy()
        u = u.sort_values("pos")
        # RMC regions on this chromosome, by Ensembl gene id
        rc = reg[(~reg["whole_gene"]) & reg["lifted"] & (reg["chrom"] == c)].sort_values("start38")
        m = pd.merge_asof(u[["key", "pos", "ensg"]], rc[["start38", "end38", "ensg", "region_id", "oe"]],
                          left_on="pos", right_on="start38", by="ensg", direction="backward")
        inside = m["end38"].notna() & (m["pos"] <= m["end38"])
        m.loc[~inside, ["region_id", "oe"]] = np.nan
        wg = reg[reg["whole_gene"]].drop_duplicates("ensg").set_index("ensg")
        no_reg = m["region_id"].isna() & m["ensg"].isin(wg.index) & ~m["ensg"].isin(rc["ensg"])
        m.loc[no_reg, "region_id"] = wg["region_id"].reindex(m.loc[no_reg, "ensg"]).to_numpy()
        m.loc[no_reg, "oe"] = wg["oe"].reindex(m.loc[no_reg, "ensg"]).to_numpy()
        # CCR percentile
        cc = ccr[ccr["chrom"] == c].sort_values("s1")
        m2 = pd.merge_asof(u[["key", "pos"]], cc[["s1", "e1", "ccr_pct"]], left_on="pos", right_on="s1", direction="backward")
        ccr_in = m2["e1"].notna() & (m2["pos"] <= m2["e1"])
        out = pd.DataFrame({"key": m["key"].to_numpy(np.int64),
                            "rmc_unit": m["region_id"].fillna(-1).to_numpy(np.int64),
                            "rmc_oe": m["oe"].to_numpy(np.float32)})
        ccr_v = pd.Series(np.where(ccr_in, m2["ccr_pct"], np.nan), index=m2["key"].to_numpy())
        out["ccr_pct"] = ccr_v.reindex(out["key"]).to_numpy(np.float32)
        out.to_pickle(f)
        log(f"chr{c}: RMC region for {(out['rmc_unit'] >= 0).mean():.3f} of missense sites, CCR for "
            f"{out['ccr_pct'].notna().mean():.3f}")


# ======================================================================== generic analysis driver
AM_COARSE = [e / 100 for e in gate.LR_EDGES_PCT]
Q_COARSE = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 980, 990, 1000]   # fine-bin indices: deciles + top 2%/1%


PSEUDO_K = 64
_PSEUDO = {}
_EXCL_TX = None                                              # A12.3: transcripts excluded (top-50 excess genes)


def pseudo_setup(a):
    """Per transcript: number of its gene's RMC regions (1 if none) and protein length, for the pseudo-region control."""
    if "n" in _PSEUDO:
        return
    reg = pd.read_pickle(P(a.fs, "regions", "rmc_regions.pkl"))
    reg = reg[reg["lifted"] | reg["whole_gene"]]
    nreg = reg.groupby("ensg").size()
    t = load_tinfo(a)
    txt = pd.read_pickle(P(a.work, "gtf", "transcripts.pkl")).set_index("tx")
    n = t["ensg"].map(nreg).fillna(1).clip(1, PSEUDO_K).astype(int)
    _PSEUDO["n"] = np.r_[n.reindex(np.arange(18398)).fillna(1).to_numpy(int)]
    _PSEUDO["L"] = np.r_[(txt["cds_len"] // 3).reindex(np.arange(18398)).fillna(1).to_numpy(int)]


WIN = 50                     # amendment A1, D2: fixed 50-residue windows
_POS = {}


def pos_setup(a):
    """Amendment A1 unit numbering, compact per transcript: D2 windows, E1 (chunk x RMC region), E2 (window x region)."""
    if "off_w" in _POS:
        return
    pseudo_setup(a)
    reg = pd.read_pickle(P(a.fs, "regions", "rmc_regions.pkl"))
    t = load_tinfo(a)
    nt = 18398
    L = np.maximum(_PSEUDO["L"], 1)
    nwin = (L + WIN - 1) // WIN
    nreg = _PSEUDO["n"]                                     # RMC regions per gene (>= 1), capped at PSEUDO_K
    base = reg.groupby("ensg")["region_id"].min()           # region ids are contiguous within a gene
    b = t["ensg"].map(base).reindex(np.arange(nt))
    _POS["rmc_base"] = b.fillna(-1).to_numpy(np.int64)
    _POS["nreg"] = nreg
    _POS["off_w"] = np.concatenate([[0], np.cumsum(nwin)])
    _POS["nwin"] = nwin
    _POS["off_e1"] = np.concatenate([[0], np.cumsum(nreg * nreg)])
    _POS["off_e2"] = np.concatenate([[0], np.cumsum(nwin * nreg)])
    nwin_s = (L + WIN // 2 + WIN - 1) // WIN                 # amendment A4.2: windows shifted by 25 residues
    _POS["nwin_s"] = nwin_s
    _POS["off_ws"] = np.concatenate([[0], np.cumsum(nwin_s)])
    for size in SWEEP:                                       # amendment A11.4: unit-size sweep
        nw = (L + size - 1) // size
        _POS[f"nwin{size}"] = nw
        _POS[f"off_w{size}"] = np.concatenate([[0], np.cumsum(nw)])


SWEEP = (25, 100)
POS_OFF = {"win50": "off_w", "pseudo_rmc": "off_e1", "win50_rmc": "off_e2", "win50s": "off_ws", "rand50": "off_w",
           "win25": "off_w25", "rand25": "off_w25", "win100": "off_w100", "rand100": "off_w100"}
_DOM = {}


def dom_setup(a):
    """Amendment A4.2: gene x InterPro-domain units (dom) and dom x RMC units (dom_rmc), built by `dom54_stream`."""
    if "n_dom" not in _DOM:
        z = np.load(P(a.fs, "regions", "dom_units.npz"))
        _DOM.update({k: z[k] for k in z.files})


def _rand50(u, nkey="nwin", okey="off_w"):
    """Amendment A4.2 random-partition control: each residue to one of ceil(L/size) groups of its gene by a fixed hash
    (size 50 by default; A11.4 uses the same hash with 25- and 100-residue group counts)."""
    tx = u["tx"].to_numpy(np.int64)
    z = tx.astype(np.uint64) * np.uint64(1_000_003) + u["prot_pos"].to_numpy(np.int64).astype(np.uint64) + np.uint64(SEED)
    with np.errstate(over="ignore"):                         # splitmix64 finaliser (wrapping arithmetic intended)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        z = z ^ (z >> np.uint64(31))
    g = (z % _POS[nkey][tx].astype(np.uint64)).astype(np.int64)
    return _POS[okey][tx] + g


def _local_rmc(u):
    tx = u["tx"].to_numpy(np.int64)
    r = u["rmc_unit"].to_numpy(np.int64)
    loc = r - _POS["rmc_base"][tx]
    ok = (r >= 0) & (_POS["rmc_base"][tx] >= 0) & (loc >= 0)
    return np.where(ok, np.minimum(loc, _POS["nreg"][tx] - 1), -1), tx


def unit_of(u, kind):
    if kind == "win50s":
        tx = u["tx"].to_numpy(np.int64)
        w = np.minimum((u["prot_pos"].to_numpy(np.int64) - 1 + WIN // 2) // WIN, _POS["nwin_s"][tx] - 1)
        return _POS["off_ws"][tx] + w
    if kind == "rand50":
        return _rand50(u)
    if kind in ("win25", "win100"):
        size = int(kind[3:])
        tx = u["tx"].to_numpy(np.int64)
        w = np.minimum((u["prot_pos"].to_numpy(np.int64) - 1) // size, _POS[f"nwin{size}"][tx] - 1)
        return _POS[f"off_w{size}"][tx] + w
    if kind in ("rand25", "rand100"):
        size = int(kind[4:])
        return _rand50(u, f"nwin{size}", f"off_w{size}")
    if kind == "dom":
        return u["dom_unit"].to_numpy(np.int64)
    if kind == "dom_rmc":
        return u["domrmc_unit"].to_numpy(np.int64)
    if kind == "win50":
        tx = u["tx"].to_numpy(np.int64)
        w = np.minimum((u["prot_pos"].to_numpy(np.int64) - 1) // WIN, _POS["nwin"][tx] - 1)
        return _POS["off_w"][tx] + w
    if kind == "pseudo_rmc":
        loc, tx = _local_rmc(u)
        n = _PSEUDO["n"][tx]
        chunk = np.minimum(((u["prot_pos"].to_numpy(np.int64) - 1) * n) // np.maximum(_PSEUDO["L"][tx], 1), n - 1)
        return np.where(loc >= 0, _POS["off_e1"][tx] + chunk * n + loc, -1)
    if kind == "win50_rmc":
        loc, tx = _local_rmc(u)
        w = np.minimum((u["prot_pos"].to_numpy(np.int64) - 1) // WIN, _POS["nwin"][tx] - 1)
        return np.where(loc >= 0, _POS["off_e2"][tx] + w * _POS["nreg"][tx] + loc, -1)
    if kind == "pseudo":
        tx = u["tx"].to_numpy(np.int64)
        n, L = _PSEUDO["n"][tx], np.maximum(_PSEUDO["L"][tx], 1)
        chunk = np.minimum(((u["prot_pos"].to_numpy(np.int64) - 1) * n) // L, n - 1)
        return tx * PSEUDO_K + chunk
    if kind == "gene":
        return u["tx"].to_numpy(np.int64)
    if kind == "rmc":
        return u["rmc_unit"].to_numpy(np.int64)
    if kind == "ccr":
        cls_ = np.digitize(u["ccr_pct"].to_numpy(np.float64), CCR_EDGES) - 1
        return np.where(u["ccr_pct"].notna().to_numpy(), u["tx"].to_numpy(np.int64) * 6 + cls_, -1)
    raise ValueError(kind)


def n_units(kind, a):
    if kind in POS_OFF:
        pos_setup(a)
        return int(_POS[POS_OFF[kind]][-1]) + 1
    if kind in ("dom", "dom_rmc"):
        dom_setup(a)
        return int(_DOM["n_" + kind.replace("dom_rmc", "domrmc")]) + 1
    if kind == "pseudo":
        return (18397 + 1) * PSEUDO_K
    if kind == "gene":
        return 18397 + 1
    if kind == "ccr":
        return (18397 + 1) * 6
    return len(pd.read_pickle(P(a.fs, "regions", "rmc_regions.pkl"))) + 1


def unit_tx(kind, a, nu):
    """Transcript of each unit (for strata defined on genes)."""
    if kind in POS_OFF:
        pos_setup(a)
        off = _POS[POS_OFF[kind]]
        out = np.repeat(np.arange(len(off) - 1), np.diff(off))
        return np.r_[out, -1][:nu]
    if kind in ("dom", "dom_rmc"):
        dom_setup(a)
        return np.r_[_DOM["tx_" + kind.replace("dom_rmc", "domrmc")], -1][:nu]
    if kind == "pseudo":
        return np.arange(nu) // PSEUDO_K
    if kind == "gene":
        return np.arange(nu)
    if kind == "ccr":
        return np.arange(nu) // 6
    reg = pd.read_pickle(P(a.fs, "regions", "rmc_regions.pkl"))
    tinfo = load_tinfo(a)
    e2t = pd.Series(tinfo.index.to_numpy(), index=tinfo["ensg"].to_numpy())
    e2t = e2t[~e2t.index.duplicated()]
    out = np.full(nu, -1)
    out[reg["region_id"].to_numpy()] = e2t.reindex(reg["ensg"]).fillna(-1).to_numpy(np.int64)
    return out


class Collect:
    """One streaming pass: per score x unit kind, accumulate opportunity per (unit, fine bin) and map de novos."""

    def __init__(self, a, scores, kinds, d_mis, near, restrict=None, fine_edges=None):
        self.a, self.scores, self.kinds, self.near, self.restrict = a, scores, kinds, near, restrict
        self.edges = fine_edges or {}
        self.d_mis = d_mis
        self.dmap = {(s, k): [] for s in scores for k in kinds}

    def cols(self):
        """Columns a pass needs (the column store then reads only these)."""
        c = ["key", "tx", PRIMARY_MODEL]
        for s_ in self.scores:
            c += ["am_bin"] if s_ == "AM_gate" else (["am"] if s_ == "AlphaMissense" else [score_col(s_)])
        if any(k in ("rmc", "ccr") for k in self.kinds) or self.restrict:
            c += ["rmc_unit", "ccr_pct"]
        if any(k in ("pseudo", "win50", "pseudo_rmc", "win50_rmc", "win50s", "rand50", "win25", "rand25", "win100",
                     "rand100") for k in self.kinds):
            c += ["prot_pos"]
        if any(k in ("pseudo_rmc", "win50_rmc") for k in self.kinds):
            c += ["rmc_unit"]
        if "dom" in self.kinds:
            c += ["dom_unit"]
        if "dom_rmc" in self.kinds:
            c += ["domrmc_unit"]
        return list(dict.fromkeys(c))

    def extra(self):
        ex = []
        if any(s not in ("AlphaMissense", "AM_gate") for s in self.scores):
            ex.append("scores")
        if any(k in ("rmc", "ccr") for k in self.kinds) or self.restrict:
            ex.append("regions")
        return tuple(ex)

    def score_of(self, u, s):
        return u["am"].to_numpy(np.float64) if s == "AlphaMissense" else u[score_col(s)].to_numpy(np.float64)

    def nbins(self, s):
        return len(gate.LR_EDGES_PCT) - 1 if s == "AM_gate" else N_QBINS

    def bins_of(self, u, s):
        if s == "AM_gate":
            f2c = np.searchsorted(np.array(gate.LR_EDGES_PCT), np.arange(100), side="right") - 1
            return f2c[u["am_bin"].to_numpy().astype(int)].astype(np.int64)
        return self.edges[s].bin(self.score_of(u, s))

    def edges_pass(self):
        hs = {s: [] for s in self.scores if s not in self.edges and s != "AM_gate"}
        if not hs:
            return
        chunks = {s: [] for s in hs}
        for c, u in site_chunks(self.a, near=self.near, extra=self.extra(), cols=self.cols()):
            u = self._restrict(u)
            for s in hs:
                chunks[s].append((self.score_of(u, s), u[PRIMARY_MODEL].to_numpy()))
            if sum(len(x[0]) for x in chunks[next(iter(hs))]) > 5e6:     # flush to histogram to bound memory
                for s in hs:
                    hs[s].append(self._hist(chunks[s]))
                    chunks[s] = []
        for s in hs:
            hs[s].append(self._hist(chunks[s]))
            h = np.sum(hs[s], 0)
            self.edges[s] = self._edges_from_hist(h)

    @staticmethod
    def _hist(pairs, grid=200001):
        h = np.zeros(grid)
        for sc, w in pairs:
            ok = np.isfinite(sc)
            idx = np.clip(np.round(sc[ok] * (grid - 1)).astype(np.int64), 0, grid - 1)
            h += np.bincount(idx, weights=w[ok], minlength=grid)
        return h

    @staticmethod
    def _edges_from_hist(h, nq=N_QBINS):
        return BinMap(h, nq)

    def _restrict(self, u):
        if self.restrict == "aim1_common":
            u = u[(u["rmc_unit"] >= 0) & u["ccr_pct"].notna()]
        return u

    def mass_pass(self, nunits):
        """nunits: {kind: number of units}. Opportunity is summed into sparse matrices chromosome by chromosome."""
        dm = self.d_mis
        self.M = {(s, k): None for s in self.scores for k in self.kinds}
        for c, u in site_chunks(self.a, near=self.near, extra=self.extra(), cols=self.cols()):
            u = self._restrict(u)
            w = u[PRIMARY_MODEL].to_numpy()
            dc = dm[dm["chrom"] == c][["key", "p", "tx", "centre"]]
            idx_d = np.flatnonzero(np.isin(u["key"].to_numpy(), dc["key"].to_numpy()))
            keys_d = u["key"].to_numpy()[idx_d]
            units = {k: unit_of(u, k) for k in self.kinds}
            for s in self.scores:
                b = self.bins_of(u, s)
                for k in self.kinds:
                    un = units[k]
                    ok = (b >= 0) & (un >= 0)
                    part = sparse.coo_matrix((w[ok], (un[ok], b[ok])), shape=(nunits[k], self.nbins(s))).tocsr()
                    self.M[(s, k)] = part if self.M[(s, k)] is None else self.M[(s, k)] + part
                    okd = ok[idx_d]
                    site = pd.DataFrame({"key": keys_d[okd], "bin": b[idx_d][okd], "unit": un[idx_d][okd]})
                    self.dmap[(s, k)].append(dc.merge(site, on="key", how="inner"))
            del units
            gc.collect()

    def result(self, s, k):
        return self.M[(s, k)], pd.concat(self.dmap[(s, k)], ignore_index=True)


def syn_inputs(a, d, near):
    syn_mass = 0.0
    for c, u in site_chunks(a, near=near, cls=0):
        syn_mass += float(u[PRIMARY_MODEL].sum())
    ds = dnv_filter(d[d["status"] == "synonymous"], near)
    return syn_mass, ds["p"].to_numpy(np.int64), ds


def report(eng, M, syn_mass, syn_p, dmap, edges, unit_mask, tag, explicit_edges=None):
    """Runs the engine; returns (coarse rows, window rows, threshold rows, raw result) for marginal and within-unit LRs.
    explicit_edges: the bins are already the reporting bins (e.g. the gate's AlphaMissense bins) - no windows."""
    res = eng.run(M, syn_mass, syn_p, dmap, unit_mask)
    rows, wrows, trows = [], [], []
    xt = summ((res["O"] - res["E"]).sum(1)[:, None])
    if not xt["lo95"][0] > 0:
        return rows, wrows, trows, res
    if explicit_edges is not None:
        groups = None
        labels = list(zip(explicit_edges[:-1], explicit_edges[1:]))
    else:
        groups = [np.arange(lo_, hi_) for lo_, hi_ in zip(Q_COARSE[:-1], Q_COARSE[1:])]
        with np.errstate(all="ignore"):
            labels = [(float(np.nanmin(edges.lo[g])) if np.isfinite(edges.lo[g]).any() else np.nan,
                       float(np.nanmax(edges.hi[g])) if np.isfinite(edges.hi[g]).any() else np.nan) for g in groups]
    lrc = Engine.lr(res, groups)
    for kind, key in (("marginal", "LRm"), ("within_unit", "LRw")):
        s_ = summ(lrc[key]); lo_s = summ(lrc[key + "_dn"])["lo90"]; hi_s = summ(lrc[key + "_up"])["hi90"]
        for i, (lo_, hi_) in enumerate(labels):
            rows.append(tag + [kind, lo_, hi_, lrc["O"][0, i], lrc["E"][0, i], lrc["X"][0, i], s_["est"][i], s_["lo95"][i],
                               s_["hi95"][i], lo_s[i], hi_s[i], points(np.array([s_["est"][i]]))[0],
                               gate.evidence_label(lo_s[i], hi_s[i])])
        if explicit_edges is not None:
            continue
        centres, wg = window_groups()
        lrw = Engine.lr(res, wg)
        valid = lrw["E"][0] > 0                 # windows wholly inside the empty bins after a tied block carry nothing
        lrw = {k_: v_[:, valid] for k_, v_ in lrw.items()}
        cs = edges.centres()[centres][valid]
        sw = summ(lrw[key]); lo_w = summ(lrw[key + "_dn"])["lo90"]; hi_w = summ(lrw[key + "_up"])["hi90"]
        for i in range(len(cs)):
            wrows.append(tag + [kind, cs[i], lrw["O"][0, i], lrw["E"][0, i], sw["est"][i], sw["lo95"][i], sw["hi95"][i],
                                lo_w[i], hi_w[i]])
        th = window_thresholds(lo_w, hi_w, cs)
        trows.append(tag + [kind] + [th[k] for k in THR_KEYS])
    return rows, wrows, trows, res


COARSE_COLS = ["kind", "score_lo", "score_hi", "O", "E", "X", "LR", "lo95", "hi95", "lo_sys", "hi_sys", "points", "evidence"]
WIN_COLS = ["kind", "centre_score", "O", "E", "LR", "lo95", "hi95", "lo_sys", "hi_sys"]
THR_KEYS = [f"PP3_{n}" for n, _ in PP3] + [f"BP4_{n}" for n, _ in PP3]


# ======================================================================== replication (amendment A2)
REP_FILE = P(EXT, "replication", "zhou2022_dnv_combined.tsv")
REP_COHORTS = ["SPARK", "SSC", "ASC", "MSSNG"]


def load_dnv_replication(a, pheno):
    """Zhou 2022 de novo SNVs (Supplementary Data 1+2) for one phenotype, de-duplicated against Kaplanis 2020,
    lifted to GRCh38 with the gate's chain, REF-checked against Roulette records, annotated from the universe."""
    from pyliftover import LiftOver
    z = pd.read_csv(getattr(a, "rep_file", None) or REP_FILE, sep="\t", dtype=str)
    cohorts = getattr(a, "rep_cohorts", None) or REP_COHORTS
    z = z[(z["Pheno"] == pheno) & z["Cohort"].isin(cohorts) & (z["Ref"].str.len() == 1) & (z["Alt"].str.len() == 1)]
    z = z[z["Chrom"].isin(AUTOSOMES)].copy()
    n0 = len(z)
    # duplicates: an individual sharing >= 2 identical de novo SNVs with one Kaplanis proband (GRCh37 in both)
    k = pd.read_csv(P(a.data, gate.DNV_FILE), sep="\t", dtype=str, usecols=["id", "chrom", "pos", "ref", "alt"])
    k["v"] = k["chrom"] + ":" + k["pos"] + ":" + k["ref"] + ":" + k["alt"]
    z["v"] = z["Chrom"] + ":" + z["Position"] + ":" + z["Ref"] + ":" + z["Alt"]
    shared = z[["IID", "v"]].merge(k[["id", "v"]], on="v")
    pair = shared.groupby(["IID", "id"]).size()
    dup = set(pair[pair >= 2].index.get_level_values(0))
    if getattr(a, "dedup_zhou", False):          # amendment A9: also against Zhou 2022 individuals (GRCh37 in both)
        zz = pd.read_csv(REP_FILE, sep="\t", dtype=str, usecols=["IID", "Chrom", "Position", "Ref", "Alt"])
        zz["v"] = zz["Chrom"] + ":" + zz["Position"] + ":" + zz["Ref"] + ":" + zz["Alt"]
        sh2 = z[["IID", "v"]].merge(zz[["IID", "v"]].rename(columns={"IID": "zid"}), on="v")
        p2 = sh2.groupby(["IID", "zid"]).size()
        dup |= set(p2[p2 >= 2].index.get_level_values(0))
    z = z[~z["IID"].isin(dup)]
    lo = LiftOver(P(a.data, gate.CHAIN_FILE))
    p38, ok = [], []
    for c, pos in zip(z["Chrom"], z["Position"].astype(int)):
        r = lo.convert_coordinate("chr" + c, pos - 1)
        good = bool(r) and r[0][0] == "chr" + c and r[0][2] == "+"
        p38.append(r[0][1] + 1 if good else -1)
        ok.append(good)
    z["pos38"] = p38
    z = z[np.array(ok)]
    z["key"] = z["pos38"].astype(np.int64) * 4 + z["Alt"].map(B2I).astype(np.int64)
    parts = []
    for c in AUTOSOMES:
        zc = z[z["Chrom"] == c]
        if not len(zc):
            continue
        u = pd.read_pickle(P(a.fs, "universe", f"chr{c}.pkl"))
        u = u[u["key"].isin(zc["key"])][["key", "tx", "cls", "filt_ok", "near", PRIMARY_MODEL]]
        rr = pd.read_csv(P(a.work, "roulette", f"chr{c}.tsv.gz"), sep="\t", header=None, usecols=[0, 1], names=["pos", "ref"])
        rr = rr[rr["pos"].isin(zc["pos38"])].drop_duplicates("pos").set_index("pos")["ref"]
        zc = zc.assign(ref38=zc["pos38"].map(rr).to_numpy()).merge(u, on="key", how="inner")
        parts.append(zc[zc["ref38"] == zc["Ref"]])
    d = pd.concat(parts, ignore_index=True)
    d["status"] = np.where(d["cls"] == 1, "missense", np.where(d["cls"] == 0, "synonymous", "other"))
    d = d[d["status"] != "other"].copy()
    d["chrom"] = d["Chrom"]
    d["centre"] = d["Cohort"]
    d["tx"] = d["tx"].astype(np.int64)
    ids = pd.Series(np.arange(z["IID"].nunique()), index=pd.unique(z["IID"]))
    d["p"] = ids.reindex(d["IID"]).to_numpy()
    log(f"replication [{pheno}]: {n0:,} autosomal SNVs; {len(dup)} individuals removed as duplicates of Kaplanis probands; "
        f"{len(d):,} on universe missense/synonymous sites with REF confirmed "
        f"({int((d['status'] == 'missense').sum()):,} missense, {int((d['status'] == 'synonymous').sum()):,} synonymous); "
        f"{len(ids):,} individuals")
    with open(P(a.results, "replication_input.txt"), "a") as f:
        f.write(f"{pheno}\tSNVs {n0}\tduplicates_removed {len(dup)}\tmapped {len(d)}\tmissense "
                f"{int((d['status'] == 'missense').sum())}\tsynonymous {int((d['status'] == 'synonymous').sum())}\t"
                f"individuals {len(ids)}\n")
    return d, len(ids)


# ======================================================================== aims
def setup(a, near="exclude"):
    os.makedirs(a.results, exist_ok=True)
    if getattr(a, "cohort", None):
        d, n_prob = load_dnv_replication(a, a.cohort)
    else:
        d, n_prob = load_dnv(a)
    if getattr(a, "centre_only", None):
        d = d[d["centre"] == a.centre_only].copy()
    eng = Engine(n_prob)
    syn_mass, syn_p, _ = syn_inputs(a, d, near)
    d_mis = dnv_filter(d[d["status"] == "missense"], near)
    tinfo = load_tinfo(a)
    return d, eng, syn_mass, syn_p, d_mis, tinfo


def gene_masks(tinfo, nu, utx):
    """Strata defined on genes, as boolean masks over units (utx = transcript of each unit; -1 = none)."""
    nt = 18398
    t = tinfo.reindex(np.arange(nt))
    ok = utx >= 0
    ix = np.where(ok, utx, 0)

    def per_tx(v):
        return np.asarray(v)[ix] & ok

    dd = per_tx(t["ddg2p"].fillna(False).to_numpy(bool))
    m = {"all_genes": ok.copy(), "ddg2p_monoallelic": dd,
         "shet_high": per_tx((t["shet_class"] == 1).to_numpy()), "shet_low": per_tx((t["shet_class"] == 0).to_numpy())}
    for lab in ("LoF", "nonLoF"):
        m[f"ddg2p_{lab}"] = dd & per_tx((t["mech12_bin"] == lab).to_numpy())
        m[f"allgenes_{lab}"] = per_tx((t["mech123_bin"] == lab).to_numpy())
    for lab in ("lof", "dn", "gof"):
        m[f"ddg2p_{lab}"] = dd & per_tx((t["mech12"] == lab).to_numpy())
    m["a5_cold_genes"] = per_tx(t["gene"].fillna("").map(gene_split).to_numpy(bool))   # A5/A6: Aim 4 cold-test genes
    if "n" in _PSEUDO:                                       # A11.3: genes with >= 2 RMC regions vs one whole-gene region
        m["rmc_multi"] = per_tx(_PSEUDO["n"][:nt] >= 2)
        m["rmc_single"] = per_tx(_PSEUDO["n"][:nt] < 2)
    if _EXCL_TX is not None:                                 # A12.3: leave the top-50 excess genes out
        m["no_top50"] = ok & ~np.isin(utx, _EXCL_TX)
    return m


def write(rows, cols, path, mode="w"):
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(path, sep="\t", index=False, mode=mode, header=(mode == "w" or not os.path.exists(path)))
    return df


def cmd_validate(a):
    """The engine must reproduce the gate's AlphaMissense table (gate filters: exon-edge sites included)."""
    global DELTA
    nc = pd.read_csv(P(HERE, "results", "null_checks.tsv"), sep="\t")
    gated = nc["grouping"].isin(["Roulette mutability decile", "AlphaMissense context of the codon", "mutation class"])
    DELTA = max(gate.DELTA_MIN, float((nc.loc[gated, "OE"] - 1).abs().max()))    # the gate's unrounded delta, for this check only
    log(f"validation uses the gate's unrounded delta {DELTA:.6f} (analyses use the pre-registered 0.091)")
    d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a, near="all")
    col = Collect(a, ["AM_gate"], ["gene"], d_mis, near="all")
    col.mass_pass({"gene": n_units("gene", a)})
    M, dm = col.result("AM_gate", "gene")
    masks = gene_masks(tinfo, M.shape[0], unit_tx("gene", a, M.shape[0]))
    g = pd.read_csv(P(HERE, "results", "lr_bins.tsv"), sep="\t")
    edges = [e / 100 for e in gate.LR_EDGES_PCT]
    worst = 0.0
    for st in ("all_genes", "ddg2p_monoallelic", "shet_high", "shet_low"):
        rows, _, _, _ = report(eng, M, syn_mass, syn_p, dm, None, masks[st], [st], explicit_edges=edges)
        mine = pd.DataFrame(rows, columns=["stratum"] + COARSE_COLS)
        for kind, gk in (("marginal", "marginal"), ("within_unit", "within_gene")):
            x = mine[mine["kind"] == kind]
            y = g[(g["analysis"] == st) & (g["kind"] == gk)]
            for c1, c2 in (("O", "O"), ("E", "E"), ("LR", "LR"), ("lo95", "lo95"), ("hi95", "hi95"), ("lo_sys", "lo_one_sided95_sys")):
                diff = np.nanmax(np.abs(x[c1].to_numpy() - y[c2].to_numpy()) / np.maximum(1e-9, np.abs(y[c2].to_numpy())))
                worst = max(worst, diff)
    log(f"validation against the gate's lr_bins.tsv: max relative difference {worst:.2e}")
    if worst > 1e-6:
        die("engine does not reproduce the gate")


def cmd_aim5(a):
    d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a)
    col = Collect(a, ["AM_gate"], ["gene"], d_mis, near="exclude")
    col.mass_pass({"gene": n_units("gene", a)})
    M, dm = col.result("AM_gate", "gene")
    masks = gene_masks(tinfo, M.shape[0], unit_tx("gene", a, M.shape[0]))
    rows = []
    edges = [e / 100 for e in gate.LR_EDGES_PCT]
    for st, mk in masks.items():
        res = eng.run(M, syn_mass, syn_p, dm, mk)
        O, E = res["O"], res["E"]
        with np.errstate(invalid="ignore", divide="ignore"):
            tot = summ((1 - E.sum(1) / O.sum(1))[:, None])
            byb = summ(1 - E / O)
        ng = int(mk.sum())
        rows.append([st, "all", ng, O[0].sum(), E[0].sum(), tot["est"][0], tot["lo95"][0], tot["hi95"][0]])
        for i in range(len(edges) - 1):
            rows.append([st, f"AM {edges[i]:.2f}-{edges[i + 1]:.2f}", ng, O[0, i], E[0, i], byb["est"][i], byb["lo95"][i], byb["hi95"][i]])
    df = write(rows, ["stratum", "subset", "n_genes", "O", "E", "pathogenic_fraction", "lo95", "hi95"],
               P(a.results, "aim5_priors.tsv"))
    # genes with >= 20 missense de novos (primary filters), confirming the plan's "about 65"
    per_gene = dm.groupby("unit").size()
    n20 = int((per_gene >= 20).sum())
    with open(P(a.results, "aim5_summary.txt"), "w") as f:
        f.write(f"genes with >= 20 missense de novos (primary filters, exon-edge excluded): {n20}\n")
        f.write(f"genes with >= 10: {int((per_gene >= 10).sum())}; >= 5: {int((per_gene >= 5).sum())}\n")
        f.write(f"Pejaver 2022 prior for comparison: 0.0441\n")
        for st in ("all_genes", "ddg2p_monoallelic", "shet_high", "shet_low", "ddg2p_LoF", "ddg2p_nonLoF"):
            r = df[(df["stratum"] == st) & (df["subset"] == "all")].iloc[0]
            f.write(f"{st:22s} pi = {r['pathogenic_fraction']:.3f} [{r['lo95']:.3f}, {r['hi95']:.3f}]   "
                    f"(O {r['O']:.0f}, E {r['E']:.1f}, genes {r['n_genes']})\n")
    log(open(P(a.results, "aim5_summary.txt")).read())


def top_arrays(res, binmap):
    """(LRm, LRw, LRw with E inflated by delta, LRm with delta, opportunity share) over the top-2% group."""
    g, share = binmap.top_group(0.02)
    lr = Engine.lr(res, [g])
    return lr["LRm"][:, 0], lr["LRw"][:, 0], lr["LRw_dn"][:, 0], lr["LRm_dn"][:, 0], share


def run_scores(a, scores, kinds, strata_names, tag, restrict=None, near="exclude", explicit_am=True, top_contrast=None):
    """Shared driver for Aims 1-3: windows + coarse tables + thresholds for each score x unit kind x stratum."""
    d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a, near)
    nunits = {k: n_units(k, a) for k in kinds}
    utx = {k: unit_tx(k, a, nunits[k]) for k in kinds}
    reps = {}
    extra_masks = {}
    if "rmc" in kinds:
        reg = pd.read_pickle(P(a.fs, "regions", "rmc_regions.pkl"))
        oe = np.full(nunits["rmc"], np.nan)
        oe[reg["region_id"].to_numpy()] = reg["oe"].to_numpy()
        extra_masks[("rmc", "pm1like_rmc_oe_lt_0.4")] = oe < 0.4
        extra_masks[("rmc", "not_pm1like_rmc")] = oe >= 0.4
    if "ccr" in kinds:
        cl = np.arange(nunits["ccr"]) % 6
        extra_masks[("ccr", "pm1like_ccr_ge_95")] = cl >= 4
        extra_masks[("ccr", "not_pm1like_ccr")] = cl < 4
    paths = {x: P(a.results, f"{tag}_{x}.tsv") for x in ("coarse", "windows", "thresholds", "top2pct", "contrast")}
    for pth in paths.values():
        if os.path.exists(pth):
            os.remove(pth)
    batch = 8 if list(kinds) == ["gene"] else (4 if len(kinds) <= 3 else 1)    # tools per streaming pass (memory-bound)
    for i0 in range(0, len(scores), batch):
        sc = scores[i0:i0 + batch]
        col = Collect(a, sc, kinds, d_mis, near=near, restrict=restrict)
        col.edges_pass()
        col.mass_pass(nunits)
        for s in sc:
            for k in kinds:
                M, dm = col.result(s, k)
                gm = gene_masks(tinfo, nunits[k], utx[k])
                strata = [(n, gm[n]) for n in strata_names] + [(n, m_ & gm["all_genes"]) for (kk, n), m_ in extra_masks.items() if kk == k]
                tops = {}
                for st, mk in strata:
                    ex = [e / 100 for e in gate.LR_EDGES_PCT] if s == "AM_gate" else None
                    rows, wrows, trows, res = report(eng, M, syn_mass, syn_p, dm, col.edges.get(s), mk, [s, k, st], explicit_edges=ex)
                    write(rows, ["score", "unit", "stratum"] + COARSE_COLS, paths["coarse"], "a")
                    if wrows:
                        write(wrows, ["score", "unit", "stratum"] + WIN_COLS, paths["windows"], "a")
                        write(trows, ["score", "unit", "stratum", "kind"] + THR_KEYS, paths["thresholds"], "a")
                    if s != "AM_gate" and rows:
                        tops[st] = top_arrays(res, col.edges[s])
                        # amendment A6: keep the bootstrap replicates (identical proband weights in every run with the
                        # same proband set and seed), so retention shares get paired confidence intervals
                        reps[f"{s}|{k}|{st}|LRm"] = np.asarray(tops[st][0], np.float64)
                        reps[f"{s}|{k}|{st}|LRw"] = np.asarray(tops[st][1], np.float64)
                        lm, lw = summ(tops[st][0][:, None]), summ(tops[st][1][:, None])
                        lw_sys, lm_sys = summ(tops[st][2][:, None])["lo90"][0], summ(tops[st][3][:, None])["lo90"][0]
                        write([[s, k, st, lm["est"][0], lm["lo95"][0], lm["hi95"][0], lm_sys, lw["est"][0], lw["lo95"][0],
                                lw["hi95"][0], lw_sys, points(np.array([lw["est"][0]]))[0], tops[st][4]]],
                              ["score", "unit", "stratum", "LRm_top2", "lo95", "hi95", "lo_sys", "LRw_top2", "lo95_w", "hi95_w",
                               "lo_sys_w", "points_w", "top_share"],
                              paths["top2pct"], "a")
                if top_contrast and s != "AM_gate":
                    for s1, s2 in top_contrast:
                        if s1 in tops and s2 in tops:
                            with np.errstate(invalid="ignore", divide="ignore"):
                                diff = np.log(tops[s1][1]) - np.log(tops[s2][1])
                            sd_ = summ(diff[:, None])
                            write([[s, k, f"{s1}_vs_{s2}", sd_["est"][0], sd_["lo95"][0], sd_["hi95"][0]]],
                                  ["score", "unit", "contrast", "dlogLRw_top2", "lo95", "hi95"],
                                  P(a.results, f"{tag}_contrast.tsv"), "a")
                log(f"{tag}: {s} x {k} done")
            np.savez(P(a.results, f"{tag}_top2reps.npz"), **reps)        # rewritten after every score (crash-safe)
        del col
        gc.collect()


def cmd_aim1(a):
    run_scores(a, ["AM_gate", "AlphaMissense"], ["gene", "rmc", "ccr"], ["all_genes", "ddg2p_monoallelic"], "aim1",
               restrict="aim1_common")
    c = pd.read_csv(P(a.results, "aim1_coarse.tsv"), sep="\t")
    top = c[(c["score"] == "AM_gate") & (c["score_lo"] == 0.99)]
    with open(P(a.results, "aim1_summary.txt"), "w") as f:
        f.write("Aim 1 readout: AlphaMissense 0.99-1.00, common site set (sites with both an RMC region and a CCR value)\n")
        for st in ("all_genes", "ddg2p_monoallelic"):
            mrg = top[(top["unit"] == "gene") & (top["stratum"] == st) & (top["kind"] == "marginal")].iloc[0]
            f.write(f"\n[{st}]  marginal LR {mrg['LR']:.2f} [{mrg['lo95']:.2f}, {mrg['hi95']:.2f}], one-sided lo (delta) "
                    f"{mrg['lo_sys']:.2f}  -> {mrg['evidence']}\n")
            for k, nm in (("gene", "within gene"), ("rmc", "within gene x RMC region"), ("ccr", "within gene x CCR class")):
                r = top[(top["unit"] == k) & (top["stratum"] == st) & (top["kind"] == "within_unit")].iloc[0]
                frac = np.log(r["LR"]) / np.log(mrg["LR"]) if r["LR"] > 0 else float("nan")
                f.write(f"  {nm:28s} LR {r['LR']:.2f} [{r['lo95']:.2f}, {r['hi95']:.2f}]  one-sided lo (delta) {r['lo_sys']:.2f}"
                        f"  -> {r['evidence']};  share of log-evidence {frac:.2f}\n")
        f.write("\nPre-specified reading: primary readout = within gene x RMC region; PP3 independent of regional constraint "
                "if its one-sided lower bound >= 2.406.\n")
        pm = c[(c["score"] == "AM_gate") & (c["score_lo"] == 0.99) & c["stratum"].str.contains("pm1like")]
        for r in pm.itertuples():
            f.write(f"  {r.stratum:26s} {r.unit:4s} {r.kind:12s} LR {r.LR:.2f} [{r.lo95:.2f}, {r.hi95:.2f}] -> {r.evidence}\n")
    log(open(P(a.results, "aim1_summary.txt")).read())


def cmd_aim1x(a):
    """EXPLORATORY control (added after the Aim 1 result): within-unit LR when each gene is cut into the same number of
    evenly spaced pseudo-regions as its gnomAD RMC regions. Same common site set as Aim 1."""
    pseudo_setup(a)
    run_scores(a, ["AM_gate"], ["gene", "rmc", "pseudo"], ["all_genes", "ddg2p_monoallelic"], "aim1x_control",
               restrict="aim1_common")
    c = pd.read_csv(P(a.results, "aim1x_control_coarse.tsv"), sep="\t")
    t = c[(c["score_lo"] == 0.99) & (c["kind"] == "within_unit")]
    log("EXPLORATORY pseudo-region control, AM 0.99-1.00, within-unit LR:\n"
        + t[["unit", "stratum", "LR", "lo95", "hi95", "lo_sys"]].round(3).to_string(index=False))


def cmd_aim2(a):
    strata = ["ddg2p_monoallelic", "ddg2p_LoF", "ddg2p_nonLoF", "ddg2p_lof", "ddg2p_dn", "ddg2p_gof", "allgenes_LoF", "allgenes_nonLoF"]
    run_scores(a, ["AM_gate", "AlphaMissense"] + (list(TOOLS) if a.tools else []), ["gene"], strata, "aim2",
               top_contrast=[("ddg2p_LoF", "ddg2p_nonLoF"), ("allgenes_LoF", "allgenes_nonLoF")])


def cmd_aim3(a):
    """Every predictor on gene units; strata = Aim 3 (all genes, DDG2P) + Aim 2 mechanism strata, in one set of passes."""
    strata = ["all_genes", "ddg2p_monoallelic", "ddg2p_LoF", "ddg2p_nonLoF", "ddg2p_lof", "ddg2p_dn", "ddg2p_gof",
              "allgenes_LoF", "allgenes_nonLoF"]
    run_scores(a, ["AlphaMissense"] + list(TOOLS), ["gene"], strata, "aim3",
               top_contrast=[("ddg2p_LoF", "ddg2p_nonLoF"), ("allgenes_LoF", "allgenes_nonLoF")])


RAW_OF_TOOL = {"REVEL": "REVEL_score", "BayesDel_addAF": "BayesDel_addAF_score", "CADD": "CADD_phred",
               "PrimateAI": "PrimateAI_score", "VEST4": "VEST4_score", "MPC": "MPC_score",
               "PolyPhen2_HVAR": "Polyphen2_HVAR_score", "PolyPhen2_HDIV": "Polyphen2_HDIV_score", "SIFT": "SIFT_score",
               "FATHMM": "FATHMM_score", "GERP++": "GERP++_RS", "phyloP100way": "phyloP100way_vertebrate",
               "MutPred": "MutPred_score"}


def cmd_aim3_report(a):
    """Ranking table + rankscore thresholds converted to raw-score units (for comparison with Pejaver 2022)."""
    th = pd.read_csv(P(a.results, "aim3_thresholds.tsv"), sep="\t")
    top = pd.read_csv(P(a.results, "aim3_top2pct.tsv"), sep="\t")
    # raw value at each rankscore threshold: the least damaging raw value among sites at or above the rankscore cut
    need = th[th["score"].isin(RAW_OF_TOOL)]
    cuts = {}
    for r in need.itertuples():
        for lvl in THR_KEYS:
            v = getattr(r, lvl)
            if pd.notna(v):
                cuts.setdefault(r.score, set()).add((lvl.startswith("PP3"), float(v)))
    best = {(t, pp, v): None for t, s in cuts.items() for pp, v in s}
    for c, u in site_chunks(a, extra=("scores",)):
        for (t, pp, v) in best:
            rk = u[TOOLS[t]].to_numpy(np.float64)
            raw = u["raw_" + RAW_OF_TOOL[t]].to_numpy(np.float64)
            low_is_bad = RAW_COLS[RAW_OF_TOOL[t]] == "min"
            sel = (rk >= v) if pp else (rk <= v)
            x = raw[sel & np.isfinite(raw)]
            if not len(x):
                continue
            # PP3: boundary = least damaging raw among sites at/above the cut; BP4: most damaging among sites at/below
            if pp:
                val = x.max() if low_is_bad else x.min()
            else:
                val = x.min() if low_is_bad else x.max()
            cur = best[(t, pp, v)]
            if cur is None:
                best[(t, pp, v)] = val
            else:
                best[(t, pp, v)] = (max(cur, val) if low_is_bad else min(cur, val)) if pp else \
                                   (min(cur, val) if low_is_bad else max(cur, val))
    rows = []
    for r in th.itertuples():
        for lvl in THR_KEYS:
            v = getattr(r, lvl)
            raw = best.get((r.score, lvl.startswith("PP3"), float(v))) if pd.notna(v) and r.score in RAW_OF_TOOL else None
            rows.append([r.score, r.stratum, r.kind, lvl, v, raw])
    write(rows, ["score", "stratum", "kind", "level", "rankscore_threshold", "raw_threshold"],
          P(a.results, "aim3_thresholds_long.tsv"))
    # ranking: DDG2P within-gene window PP3 level, then ACMG points of the within-gene LR in the top 2%
    order = {"PP3_VeryStrong": 4, "PP3_Strong": 3, "PP3_Moderate": 2, "PP3_Supporting": 1}
    rk = []
    for s in th["score"].unique():
        t = th[(th["score"] == s) & (th["stratum"] == "ddg2p_monoallelic") & (th["kind"] == "within_unit")]
        lvl, thr = 0, None
        if len(t):
            for k_, o_ in order.items():
                if pd.notna(t.iloc[0][k_]) and o_ > lvl:
                    lvl, thr = o_, t.iloc[0][k_]
        tp = top[(top["score"] == s) & (top["stratum"] == "ddg2p_monoallelic")]
        ta = top[(top["score"] == s) & (top["stratum"] == "all_genes")]
        rk.append([s, lvl, thr, tp["LRw_top2"].iloc[0] if len(tp) else np.nan, tp["lo95_w"].iloc[0] if len(tp) else np.nan,
                   tp["points_w"].iloc[0] if len(tp) else np.nan, ta["LRw_top2"].iloc[0] if len(ta) else np.nan,
                   ta["LRm_top2"].iloc[0] if len(ta) else np.nan])
    df = pd.DataFrame(rk, columns=["score", "ddg2p_within_PP3_level", "at_rankscore", "ddg2p_LRw_top2", "lo95", "points",
                                   "all_LRw_top2", "all_LRm_top2"]).sort_values(["ddg2p_within_PP3_level", "points"],
                                                                                ascending=False)
    df.to_csv(P(a.results, "aim3_ranking.tsv"), sep="\t", index=False)
    log("Aim 3 ranking:\n" + df.round(3).to_string(index=False))


def cmd_a1_report(a):
    """Amendment A1 readouts: nested LRs A-E per tool, log-evidence decomposition, and the three pre-specified tests."""
    t = pd.read_csv(P(a.results, "aim1tools_top2pct.tsv"), sep="\t")
    levels = [("B", "gene"), ("C", "rmc"), ("C_ccr", "ccr"), ("D1", "pseudo"), ("D2", "win50"), ("E1", "pseudo_rmc"),
              ("E2", "win50_rmc")]
    rows = []
    for (s, st), g in t.groupby(["score", "stratum"], sort=False):
        if st not in ("all_genes", "ddg2p_monoallelic"):
            continue
        g = g.set_index("unit")
        if "gene" not in g.index:
            continue
        A = float(g.loc["gene", "LRm_top2"])
        r = {"score": s, "stratum": st, "top_share": float(g.loc["gene", "top_share"]), "A": A,
             "A_lo_sys": float(g.loc["gene", "lo_sys"])}
        for nm, u in levels:
            if u in g.index:
                r[nm] = float(g.loc[u, "LRw_top2"])
                r[nm + "_lo_sys"] = float(g.loc[u, "lo_sys_w"])
        with np.errstate(all="ignore"):
            lnA = np.log(A) if A > 1 else np.nan
            for nm, _ in levels:
                v = r.get(nm, np.nan)
                r["r_" + nm] = (np.log(v) / lnA) if (v > 0 and np.isfinite(lnA)) else np.nan
        r["gene_share"] = 1 - r["r_B"]
        r["position_beyond_gene_D1"] = r["r_B"] - r["r_D1"]
        r["constraint_beyond_position_D1"] = r["r_D1"] - r["r_E1"]
        r["position_beyond_gene_D2"] = r["r_B"] - r["r_D2"]
        r["constraint_beyond_position_D2"] = r["r_D2"] - r["r_E2"]
        r["E1_supporting"] = bool(r.get("E1_lo_sys", np.nan) >= PP3[0][1])
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(P(a.results, "a1_nested.tsv"), sep="\t", index=False)
    out = []
    for st in ("all_genes", "ddg2p_monoallelic"):
        x = df[df["stratum"] == st]
        conf = x[x["score"] != "AlphaMissense"]                    # AlphaMissense D/E are exploratory (amendment A1)
        inf = conf[conf["A"] >= PP3[0][1]]
        frac = float((inf["r_D1"] < 0.5).mean()) if len(inf) else float("nan")
        cbp = inf["constraint_beyond_position_D1"].dropna()
        out.append(f"[{st}] confirmatory tools: {len(conf)}; informative (marginal top-group LR >= 2.406): {len(inf)}")
        out.append(f"  readout 1 - share of informative tools keeping < 50% of log-evidence after gene + position (D1): "
                   f"{frac:.2f} -> {'SYSTEMATIC attenuation' if frac >= 0.75 else 'not systematic'} (pre-specified >= 0.75)")
        out.append(f"  readout 1 (D2 windows, same rule): {float((inf['r_D2'] < 0.5).mean()) if len(inf) else float('nan'):.2f}")
        out.append(f"  readout 2 - constraint beyond position, r(D1)-r(E1): median {cbp.median():.3f} "
                   f"[IQR {cbp.quantile(0.25):.3f}, {cbp.quantile(0.75):.3f}]")
        out.append(f"  readout 3 - tools whose E1 one-sided lower bound (delta) reaches 2.406: "
                   f"{int(conf['E1_supporting'].sum())} of {len(conf)}: {', '.join(conf.loc[conf['E1_supporting'], 'score'])}")
        out.append(f"  median remaining log-evidence share: B {inf['r_B'].median():.2f}, C {inf['r_C'].median():.2f}, "
                   f"D1 {inf['r_D1'].median():.2f}, D2 {inf['r_D2'].median():.2f}, E1 {inf['r_E1'].median():.2f}, "
                   f"E2 {inf['r_E2'].median():.2f}")
    txt = "\n".join(out) + "\n"
    open(P(a.results, "a1_summary.txt"), "w").write(txt)
    log("Amendment A1 readouts:\n" + txt)


def cmd_replicate(a):
    """Amendment A2: replication in autism (Zhou 2022). Affected: power check + nested A-E for every tool.
    Unaffected siblings: negative control (total missense O/E)."""
    base = a.results
    for pheno in ("Unaffected", "Affected"):
        a.results = P(base, f"replication_{pheno.lower()}")
        os.makedirs(a.results, exist_ok=True)
        a.cohort = pheno
        d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a)
        col = Collect(a, ["AM_gate"], ["gene"], d_mis, near="exclude")
        col.mass_pass({"gene": n_units("gene", a)})
        M, dm = col.result("AM_gate", "gene")
        res = eng.run(M, syn_mass, syn_p, dm)
        oe = summ((res["O"].sum(1) / res["E"].sum(1))[:, None])
        with np.errstate(all="ignore"):
            top = summ((res["O"][:, -1] / res["E"][:, -1])[:, None])
        msg = (f"[{pheno}] missense O {res['O'][0].sum():.0f}, E {res['E'][0].sum():.1f}, O/E {oe['est'][0]:.3f} "
               f"[{oe['lo95'][0]:.3f}, {oe['hi95'][0]:.3f}] (one-sided lo {oe['lo90'][0]:.3f}); "
               f"AM 0.99-1.00 O/E {top['est'][0]:.2f} [{top['lo95'][0]:.2f}, {top['hi95'][0]:.2f}]")
        log(msg)
        open(P(a.results, "summary.txt"), "w").write(msg + "\n")
        if pheno == "Unaffected":
            continue
        rows, _, _, _ = report(eng, M, syn_mass, syn_p, dm, None, None, ["AM_gate", "gene", "all_genes"],
                               explicit_edges=[e / 100 for e in gate.LR_EDGES_PCT])
        write(rows, ["score", "unit", "stratum"] + COARSE_COLS, P(a.results, "am_gate_coarse.tsv"))
        pos_setup(a)
        run_scores(a, ["AlphaMissense"] + list(TOOLS), ["gene", "rmc", "ccr", "pseudo", "win50", "pseudo_rmc", "win50_rmc"],
                   ["all_genes"], "rep", restrict="aim1_common")
        rep_report(a, float(oe["lo90"][0]))
    a.results = base
    a.cohort = None


def rep_report(a, oe_lo_one_sided, tag="rep", out_name="rep_summary.txt"):
    """Amendments A2/A3 readouts: power check (A2), D2-primary generalisation, constraint beyond position, E2/E1 Supporting."""
    t = pd.read_csv(P(a.results, f"{tag}_top2pct.tsv"), sep="\t")
    t = t[t["stratum"] == "all_genes"]
    lv = {"B": "gene", "C": "rmc", "D1": "pseudo", "D2": "win50", "E1": "pseudo_rmc", "E2": "win50_rmc"}
    rows = []
    for s, g in t.groupby("score", sort=False):
        g = g.set_index("unit")
        if "gene" not in g.index:
            continue
        r = {"score": s, "A": float(g.loc["gene", "LRm_top2"]), "B_lo95": float(g.loc["gene", "lo95_w"])}
        for nm, u in lv.items():
            if u in g.index:
                r[nm] = float(g.loc[u, "LRw_top2"]); r[nm + "_lo_sys"] = float(g.loc[u, "lo_sys_w"])
        with np.errstate(all="ignore"):
            lnA = np.log(r["A"]) if r["A"] > 1 else np.nan
            for nm in lv:
                v = r.get(nm, np.nan)
                r["r_" + nm] = np.log(v) / lnA if (v > 0 and np.isfinite(lnA)) else np.nan
        rows.append(r)
    df = pd.DataFrame(rows)
    for nm in lv:                                    # a level with no top group in a small cohort is reported as missing
        for col_ in (nm, nm + "_lo_sys", "r_" + nm):
            if col_ not in df:
                df[col_] = np.nan
    df.to_csv(P(a.results, f"{tag}_nested.tsv"), sep="\t", index=False)
    has_am = bool((df["score"] == "AlphaMissense").any())
    am = df[df["score"] == "AlphaMissense"].iloc[0] if has_am else pd.Series({"B_lo95": np.nan})
    powered = has_am and (oe_lo_one_sided > 1) and (am["B_lo95"] > 1)
    conf = df[df["score"] != "AlphaMissense"]
    inf = conf[conf["A"] >= PP3[0][1]]
    f2 = float((inf["r_D2"] < 0.5).mean()) if len(inf) else float("nan")
    f1 = float((inf["r_D1"] < 0.5).mean()) if len(inf) else float("nan")
    cbp2 = (inf["r_D2"] - inf["r_E2"]).dropna()
    cbp1 = (inf["r_D1"] - inf["r_E1"]).dropna()
    is_rep = bool(np.isfinite(oe_lo_one_sided))      # the power check is an A2 (replication) readout only
    out = ([f"POWER CHECK: total missense O/E one-sided lower bound {oe_lo_one_sided:.3f} (> 1 needed); AlphaMissense top-2% "
            f"within-gene LR lower 95% bound {am['B_lo95']:.2f} (> 1 needed) -> {'ADEQUATELY POWERED' if powered else 'UNDERPOWERED (readouts descriptive only)'}"]
           if is_rep else []) + [
           f"tools: {len(conf)}; informative ({'replication ' if is_rep else ''}marginal top-group LR >= 2.406): {len(inf)}",
           f"readout 1 (primary, D2): share with r(D2) < 0.5 = {f2:.2f} -> {'SYSTEMATIC' if f2 >= 0.75 else 'not systematic'} (>= 0.75); D1 secondary: {f1:.2f}",
           f"readout 2: constraint beyond position r(D2)-r(E2) median {cbp2.median():.3f} [IQR {cbp2.quantile(.25):.3f}, {cbp2.quantile(.75):.3f}]; "
           f"r(D1)-r(E1) median {cbp1.median():.3f}",
           f"readout 3: tools with E2 one-sided lower bound >= 2.406: {int((conf['E2_lo_sys'] >= PP3[0][1]).sum())} of {len(conf)} "
           f"({', '.join(conf.loc[conf['E2_lo_sys'] >= PP3[0][1], 'score'])}); E1: {int((conf['E1_lo_sys'] >= PP3[0][1]).sum())}",
           ("AlphaMissense A-E (top 2%): " + ", ".join(f"{k} {am[k]:.2f}" for k in ["A", "B", "C", "D1", "D2", "E1", "E2"]))
           if has_am else "(AlphaMissense not in this run)"]
    txt = "\n".join(out) + "\n"
    open(P(a.results, out_name), "w").write(txt)
    log(("Replication (amendment A2)" if is_rep else f"{tag} (A1 rule, D2 primary)") + " readouts:\n" + txt)


def cmd_v54(a):
    """Amendment A3: dbNSFP 5.4a tools. Aim 3 (+ Aim 2 strata) and A1 levels; popEVE on DDD de novos only (primary),
    every new tool on DDD only as sensitivity. Readouts with D2 as primary position control."""
    new = [t for t in TOOLS54 if t not in DDD_ONLY]
    strata3 = ["all_genes", "ddg2p_monoallelic", "ddg2p_LoF", "ddg2p_nonLoF", "ddg2p_lof", "ddg2p_dn", "ddg2p_gof",
               "allgenes_LoF", "allgenes_nonLoF"]
    kinds = ["gene", "rmc", "ccr", "pseudo", "win50", "pseudo_rmc", "win50_rmc"]
    pos_setup(a)
    # run order = importance (study lead's instruction, 2026-10-05); every analysis below is unchanged from amendment A3
    # 1. the pre-registered A3 confirmatory test (A1 levels on the new tools, D2 primary)
    run_scores(a, new, kinds, ["all_genes", "ddg2p_monoallelic"], "a1v54", restrict="aim1_common")
    rep_report(a, float("nan"), tag="a1v54", out_name="a1v54_summary.txt")
    # 2. popEVE, primary analysis on DDD de novos only (section 6)
    a.centre_only = "DDD"
    run_scores(a, ["popEVE"], kinds, ["all_genes"], "a1v54_popeve_ddd", restrict="aim1_common")
    run_scores(a, ["popEVE"], ["gene"], ["all_genes", "ddg2p_monoallelic"], "aim3v54_popeve_ddd")
    a.centre_only = None
    # 3. Aim 3 table (+ Aim 2 strata) for the new tools
    run_scores(a, new, ["gene"], strata3, "aim3v54", top_contrast=[("ddg2p_LoF", "ddg2p_nonLoF")])
    # 4. DDD-only sensitivity for the other new tools
    a.centre_only = "DDD"
    run_scores(a, new, ["gene"], ["all_genes", "ddg2p_monoallelic"], "aim3v54_ddd")
    run_scores(a, new, kinds, ["all_genes"], "a1v54_ddd", restrict="aim1_common")
    a.centre_only = None


def cmd_aim1tools(a):
    """Aim 1 estimands for every predictor (common site set; gene, gene x RMC region, gene x CCR class)."""
    pos_setup(a)
    run_scores(a, ["AlphaMissense"] + list(TOOLS), ["gene", "rmc", "ccr", "pseudo", "win50", "pseudo_rmc", "win50_rmc"],
               ["all_genes", "ddg2p_monoallelic"], "aim1tools", restrict="aim1_common")


# ======================================================================== Aim 4: null-contrast learning
AA = "ACDEFGHIKLMNPQRSTVWY"
PRED_LR_EDGES = [0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, np.inf]       # calibration bins of predicted LR


AIM4_COLS = ["key", "tx", "am", PRIMARY_MODEL, "ref_aa", "alt_aa", "rmc_oe", "ccr_pct"] + list(TOOLS.values())


def gene_split(sym):
    import hashlib
    return int(hashlib.md5(str(sym).encode()).hexdigest(), 16) % 5 == 0      # True = cold-test gene (20%)


def aim4_features(u, med):
    """Feature matrix (FULL_STUDY_PLAN.md section 7): scores (median-imputed + missing flags), amino acids, RMC o/e, CCR."""
    cols = {"AlphaMissense": u["am"].to_numpy(np.float32)}
    for t, c in TOOLS.items():
        cols[t] = u[c].to_numpy(np.float32)
    cols["rmc_oe"] = u["rmc_oe"].to_numpy(np.float32)
    cols["ccr_pct"] = u["ccr_pct"].to_numpy(np.float32)
    X, names = [], []
    for n, v in cols.items():
        miss = ~np.isfinite(v)
        X.append(np.where(miss, med[n], v)); names.append(n)
        X.append(miss.astype(np.float32)); names.append(n + "_missing")
    ra = u["ref_aa"].to_numpy(); aa = u["alt_aa"].to_numpy()
    for x in AA:
        X.append((ra == ord(x)).astype(np.float32)); names.append("ref_" + x)
        X.append((aa == ord(x)).astype(np.float32)); names.append("alt_" + x)
    return np.column_stack(X), names


def cmd_aim4(a):
    from sklearn.ensemble import HistGradientBoostingClassifier
    os.makedirs(a.results, exist_ok=True)
    d, eng, syn_mass, syn_p_all, d_mis, tinfo = setup(a)
    syn_d = dnv_filter(d[d["status"] == "synonymous"])
    cold_tx = tinfo["gene"].map(gene_split).to_numpy(bool)          # indexed by tx
    rng = np.random.default_rng(SEED)
    extra = ("scores", "regions")

    # ---- pass 1: training sets (GeneDx real de novos; MR-weighted null), medians for imputation
    gdx = d_mis[(d_mis["centre"] == "GeneDx") & ~cold_tx[d_mis["tx"].to_numpy()]]
    n_real = len(gdx)
    n_null = 10 * n_real
    mass_c = {}
    for c, u in site_chunks(a, extra=(), cols=["tx", PRIMARY_MODEL]):
        tr = ~cold_tx[u["tx"].to_numpy()]
        mass_c[c] = float(u.loc[tr, PRIMARY_MODEL].sum())
    tot = sum(mass_c.values())
    draws = dict(zip(mass_c, rng.multinomial(n_null, np.array(list(mass_c.values())) / tot)))
    real_parts, null_parts, med_parts = [], [], []
    for c, u in site_chunks(a, extra=extra, cols=AIM4_COLS):
        tr = ~cold_tx[u["tx"].to_numpy()]
        ut = u[tr]
        w = ut[PRIMARY_MODEL].to_numpy()
        idx = rng.choice(len(ut), size=draws[c], replace=True, p=w / w.sum())
        null_parts.append(ut.iloc[idx])
        real_parts.append(u[u["key"].isin(gdx.loc[gdx["chrom"] == c, "key"])])
        med_parts.append(ut.sample(n=min(len(ut), 20000), random_state=SEED))
    null_df = pd.concat(null_parts, ignore_index=True)
    real_df = pd.concat(real_parts, ignore_index=True)
    medsrc = pd.concat(med_parts, ignore_index=True)
    med = {"AlphaMissense": float(np.nanmedian(medsrc["am"]))}
    for t, cc in TOOLS.items():
        med[t] = float(np.nanmedian(medsrc[cc]))
    med["rmc_oe"] = float(np.nanmedian(medsrc["rmc_oe"]))
    med["ccr_pct"] = float(np.nanmedian(medsrc["ccr_pct"]))
    Xr, names = aim4_features(real_df, med)
    Xn, _ = aim4_features(null_df, med)
    X = np.vstack([Xr, Xn]); y = np.r_[np.ones(len(Xr)), np.zeros(len(Xn))]
    clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, l2_regularization=1.0,
                                         early_stopping=True, validation_fraction=0.1, random_state=SEED)
    clf.fit(X, y)
    # pi on GeneDx training genes with GeneDx's own scale factor
    c_gdx = len(syn_d[syn_d["centre"] == "GeneDx"]) / syn_mass
    E_tr = c_gdx * tot
    pi = 1 - E_tr / len(Xr)          # training de novos matched to analysed sites (= n_real on the full genome)
    ratio = len(Xn) / len(Xr)
    log(f"Aim 4: trained on {len(Xr):,} GeneDx de novos vs {len(Xn):,} null sites ({len(names)} features, "
        f"{clf.n_iter_} iterations); GeneDx pathogenic fraction pi = {pi:.3f}; GeneDx train de novos before matching {n_real:,}")

    def to_lr(Xm):
        p = clf.predict_proba(Xm)[:, 1]
        r = p / np.maximum(1 - p, 1e-12) * ratio
        return np.maximum(0.0, (r - (1 - pi)) / pi)

    # ---- pass 2: predicted LR for every analysed site (features built in slices to bound memory)
    store = {k: [] for k in ("lr", "am", "mr", "tx")}
    dkeys = d_mis[["key"]].assign(row=-1)
    drow = {}
    off = 0
    for c, u in site_chunks(a, extra=extra, cols=AIM4_COLS):
        for s0 in range(0, len(u), 400_000):
            Xm, _ = aim4_features(u.iloc[s0:s0 + 400_000], med)
            store["lr"].append(to_lr(Xm).astype(np.float32))
            del Xm
        store["am"].append(u["am"].to_numpy(np.float32)); store["mr"].append(u[PRIMARY_MODEL].to_numpy(np.float32))
        store["tx"].append(u["tx"].to_numpy(np.int32))
        kk = u["key"].to_numpy()
        hit = np.flatnonzero(np.isin(kk, d_mis.loc[d_mis["chrom"] == c, "key"].to_numpy()))
        drow.update(dict(zip(kk[hit].tolist(), (hit + off).tolist())))
        off += len(u)
        gc.collect()
    S = {k: np.concatenate(v) for k, v in store.items()}
    del store
    gc.collect()
    nu = n_units("gene", a)

    CH = 10_000_000          # rows per chunk (memory: the full arrays are 66M long)

    def mass(bins, nb, mask):
        flat = np.zeros(nu * nb)
        for s0 in range(0, len(bins), CH):
            sl = slice(s0, s0 + CH)
            ok = mask[sl] & (bins[sl] >= 0)
            flat += np.bincount(S["tx"][sl][ok].astype(np.int64) * nb + bins[sl][ok].astype(np.int64),
                                weights=S["mr"][sl][ok].astype(np.float64), minlength=nu * nb)
        return sparse.csr_matrix(flat.reshape(nu, nb))

    def binmap_of(x01):
        """Tie-robust opportunity-quantile bins (BinMap) of a [0,1] score, built and applied in chunks."""
        h = np.zeros(200001)
        for s0 in range(0, len(x01), CH):
            xs, w = x01[s0:s0 + CH], S["mr"][s0:s0 + CH]
            ok = np.isfinite(xs)
            h += np.bincount(np.clip(np.round(xs[ok].astype(np.float64) * 200000).astype(np.int64), 0, 200000),
                             weights=w[ok].astype(np.float64), minlength=200001)
        bm = BinMap(h)
        out = np.empty(len(x01), dtype=np.int16)
        for s0 in range(0, len(x01), CH):
            out[s0:s0 + CH] = bm.bin(x01[s0:s0 + CH])
        return out

    def dmap_of(dd, bins):
        r = dd["key"].map(drow)
        ok = r.notna().to_numpy()
        r = r[ok].to_numpy(np.int64)
        out = dd[ok][["p", "tx", "centre", "key"]].copy()
        out["unit"] = S["tx"][r].astype(np.int64); out["bin"] = bins[r]
        return out

    nb_lr = len(PRED_LR_EDGES) - 1
    lr_bins = np.searchsorted(np.array(PRED_LR_EDGES[1:-1]), S["lr"], side="right").astype(np.int8)
    lmax = float(np.log1p(np.nanmax(S["lr"])))
    lr_q = binmap_of((np.log1p(S["lr"]) / lmax).astype(np.float32))        # monotone in predicted LR
    am_q = binmap_of(S["am"])
    gc.collect()
    cold_site = cold_tx[S["tx"]]
    rows, top_rows = [], []
    for centre in ("DDD", "RUMC", "GeneDx"):
        sp = syn_d.loc[syn_d["centre"] == centre, "p"].to_numpy(np.int64)
        dm_c = d_mis[d_mis["centre"] == centre]
        for gs_name, gmask in (("all_genes", np.ones(len(S["tx"]), bool)), ("cold_genes", cold_site),
                               ("train_genes", ~cold_site)):
            if centre == "GeneDx" and gs_name != "cold_genes":
                continue                    # GeneDx training genes are in-sample; GeneDx is reported on cold genes only
            # calibration by predicted-LR bin
            M = mass(lr_bins, nb_lr, gmask)
            res = eng.run(M, syn_mass, sp, dmap_of(dm_c, lr_bins))
            lrc = Engine.lr(res)
            sm = summ(lrc["LRm"])
            Mt = np.asarray(M.sum(0)).ravel()
            pred_w = np.bincount(lr_bins[gmask], weights=(S["mr"] * S["lr"])[gmask], minlength=nb_lr)
            pred_mean = pred_w / np.maximum(Mt, 1e-300)
            for i in range(nb_lr):
                rows.append([centre, gs_name, PRED_LR_EDGES[i], PRED_LR_EDGES[i + 1], lrc["O"][0, i], lrc["E"][0, i],
                             pred_mean[i], sm["est"][i], sm["lo95"][i], sm["hi95"][i]])
            # top 1/2/5% of the model vs AlphaMissense (MR-weighted quantiles), marginal and within-gene
            for nm, bins in (("model", lr_q), ("AlphaMissense", am_q)):
                M2 = mass(bins, N_QBINS, gmask)
                r2 = eng.run(M2, syn_mass, sp, dmap_of(dm_c, bins))
                for pct in (1, 2, 5):
                    lr2 = Engine.lr(r2, [np.arange(N_QBINS - 10 * pct, N_QBINS)])
                    m_, w_ = summ(lr2["LRm"]), summ(lr2["LRw"])
                    top_rows.append([centre, gs_name, nm, pct, lr2["O"][0, 0], lr2["E"][0, 0], m_["est"][0], m_["lo95"][0],
                                     m_["hi95"][0], w_["est"][0], w_["lo95"][0], w_["hi95"][0]])
    write(rows, ["centre", "genes", "predLR_lo", "predLR_hi", "O", "E", "predicted_LR_mean", "empirical_LR", "lo95", "hi95"],
          P(a.results, "aim4_calibration.tsv"))
    write(top_rows, ["centre", "genes", "score", "top_pct", "O", "E", "LRm", "lo95", "hi95", "LRw", "lo95_w", "hi95_w"],
          P(a.results, "aim4_top.tsv"))
    imp = pd.DataFrame({"feature": names})
    with open(P(a.results, "aim4_model.json"), "w") as f:
        json.dump({"n_real": int(len(Xr)), "n_null": int(len(Xn)), "pi_genedx_train": float(pi), "n_iter": int(clf.n_iter_),
                   "features": names, "medians": med, "cold_gene_fraction": float(cold_tx.mean())}, f, indent=1)
    imp.to_csv(P(a.results, "aim4_features.tsv"), sep="\t", index=False)
    log("Aim 4 written: aim4_calibration.tsv, aim4_top.tsv")


# ======================================================================== amendment A5: Aim 4 v2
AF_COLS = ["af_plddt", "af_nca10", "af_plddtw"]
AF_DIR = "/mnt/c/Users/saite/OneDrive/Desktop/Neuro/G-Test-LOCKIN/ndd_pp3_gate/ndd_pp3_gate/work_af"


def cmd_af_cols(a):
    """AlphaFold structure features into the column store (plan section 3f): transcript -> Swiss-Prot entry by Ensembl
    cross-reference (versions ignored), else by gene symbol; accepted only if >= 95% of the transcript's analysed sites
    have a reference amino acid identical to the AlphaFold residue at the same position."""
    tx = pd.read_pickle(P(a.work, "gtf", "transcripts.pkl"))
    up = pd.read_csv(P(AF_DIR, "uniprot_human_sp.tsv.gz"), sep="\t")
    by_enst, by_gene = {}, {}
    for acc, gene, ens in up[["Entry", "Gene Names (primary)", "Ensembl"]].itertuples(index=False):
        if isinstance(ens, str):
            for item in ens.split(";"):
                e = item.strip().split(" ")[0].split(".")[0]
                if e.startswith("ENST"):
                    by_enst.setdefault(e, acc)
        if isinstance(gene, str):
            by_gene.setdefault(gene, []).append(acc)
    by_gene = {g: v[0] for g, v in by_gene.items() if len(v) == 1}
    af = pd.read_csv(P(AF_DIR, "af_structure_v6.tsv.gz"), sep="\t",
                     dtype={"acc": "category", "pos": np.int32, "aa": "category", "plddt": np.float32,
                            "n_ca10": np.float32, "plddt_w11": np.float32})
    prot = {}
    for acc, g in af.groupby("acc", observed=True, sort=False):
        g = g.sort_values("pos")
        if g["pos"].iloc[0] != 1 or g["pos"].iloc[-1] != len(g):
            continue                                       # positions must be 1..L
        prot[acc] = (np.frombuffer("".join(g["aa"].astype(str)).encode(), np.uint8), g["plddt"].to_numpy(),
                     g["n_ca10"].to_numpy(), g["plddt_w11"].to_numpy())
    del af
    gc.collect()
    tx["enst"] = tx["transcript_id"].str.split(".").str[0]
    cand = {int(t): [x for x in (by_enst.get(e), by_gene.get(gn)) if x in prot]
            for t, e, gn in tx[["tx", "enst", "gene_name"]].itertuples(index=False)}
    stats = {"tx_with_sites": 0, "mapped_enst": 0, "mapped_gene": 0, "rejected": 0, "no_model": 0}
    for c in AUTOSOMES:
        od = P(a.fs, "col", f"chr{c}")
        t_ = np.load(P(od, "tx.npy")).astype(np.int64)
        pp = np.load(P(od, "prot_pos.npy")).astype(np.int64)
        ra = np.load(P(od, "ref_aa.npy")).astype(np.int64)
        out = {k: np.full(len(t_), np.nan, np.float32) for k in AF_COLS}
        order = np.argsort(t_, kind="stable")
        bounds = np.flatnonzero(np.diff(t_[order])) + 1
        for idx in np.split(order, bounds):
            t = int(t_[idx[0]])
            stats["tx_with_sites"] += 1
            if not cand.get(t):
                stats["no_model"] += 1
                continue
            for j, acc in enumerate(dict.fromkeys(cand[t])):
                aa, pl, nc, pw = prot[acc]
                p0 = pp[idx] - 1
                inr = (p0 >= 0) & (p0 < len(aa))
                ident = inr & (aa[np.clip(p0, 0, len(aa) - 1)] == ra[idx])
                if ident.mean() >= 0.95:
                    k = idx[inr]
                    q = p0[inr]
                    out["af_plddt"][k], out["af_nca10"][k], out["af_plddtw"][k] = pl[q], nc[q], pw[q]
                    stats["mapped_enst" if (j == 0 and acc == by_enst.get(str(tx.loc[t, "enst"]))) else "mapped_gene"] += 1
                    break
            else:
                stats["rejected"] += 1
        for k, v in out.items():
            np.save(P(od, k + ".npy"), v)
        log(f"chr{c}: AlphaFold features at {np.isfinite(out['af_plddt']).mean():.3f} of sites")
    log(f"AlphaFold mapping (transcripts with analysed sites): {stats}")


AIM4V2_TOOLS = [t for t in TOOLS54 if t != "popEVE"]       # section 6: popEVE was developed on DD de novo cohorts
M4_COLS = {"Aim4v1": "m4_v1", "Aim4v2_GBDT": "m4_v2_gbdt", "Aim4v2_MLP": "m4_v2_mlp", "Aim4v2_ENS": "m4_v2_ens"}
M4_LRMAX = 1000.0                                           # column-store transform: log1p(min(LR, LRMAX)) / log1p(LRMAX)


def aim4_features_v2(u, med):
    X1, n1 = aim4_features(u, med)
    X, names = [X1], list(n1)
    for t in AIM4V2_TOOLS:
        v = u[score_col(t)].to_numpy(np.float32)
        miss = ~np.isfinite(v)
        X += [np.where(miss, med[t], v)[:, None], miss.astype(np.float32)[:, None]]
        names += [t, t + "_missing"]
    for k in AF_COLS:
        v = u[k].to_numpy(np.float32)
        miss = ~np.isfinite(v)
        X += [np.where(miss, med[k], v)[:, None], miss.astype(np.float32)[:, None]]
        names += [k, k + "_missing"]
    return np.hstack(X).astype(np.float32), names


def cmd_aim4v2(a):
    """Plan section 3f: v1 (section 7, re-fitted identically for the paired comparison), v2-GBDT (primary), v2-MLP,
    v2-ENS. Same training data, LR conversion, gene split and evaluation as section 7; each model's LR is also written
    to the column store (M4_COLS) for the A1 nested levels on DDD de novos."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler
    os.makedirs(a.results, exist_ok=True)
    d, eng, syn_mass, syn_p_all, d_mis, tinfo = setup(a)
    syn_d = dnv_filter(d[d["status"] == "synonymous"])
    cold_tx = tinfo["gene"].map(gene_split).to_numpy(bool)
    rng = np.random.default_rng(SEED)
    cols = list(dict.fromkeys(AIM4_COLS + [score_col(t) for t in AIM4V2_TOOLS] + AF_COLS))

    gdx = d_mis[(d_mis["centre"] == "GeneDx") & ~cold_tx[d_mis["tx"].to_numpy()]]
    n_null = 10 * len(gdx)
    mass_c = {}
    for c, u in site_chunks(a, cols=["tx", PRIMARY_MODEL]):
        mass_c[c] = float(u.loc[~cold_tx[u["tx"].to_numpy()], PRIMARY_MODEL].sum())
    tot = sum(mass_c.values())
    draws = dict(zip(mass_c, rng.multinomial(n_null, np.array(list(mass_c.values())) / tot)))
    real_parts, null_parts, med_parts = [], [], []
    for c, u in site_chunks(a, cols=cols):
        tri = np.flatnonzero(~cold_tx[u["tx"].to_numpy()])          # row positions of training-gene sites (no copy)
        w = u[PRIMARY_MODEL].to_numpy()[tri]
        null_parts.append(u.iloc[tri[rng.choice(len(tri), size=draws[c], replace=True, p=w / w.sum())]])
        real_parts.append(u[u["key"].isin(gdx.loc[gdx["chrom"] == c, "key"])])
        # identical rows to section 7's ut.sample(n, random_state=SEED): positional sample of the training rows
        smp = pd.Series(np.arange(len(tri))).sample(n=min(len(tri), 20000), random_state=SEED).to_numpy()
        med_parts.append(u.iloc[tri[smp]])
        del u, tri, w
        gc.collect()
    null_df, real_df = pd.concat(null_parts, ignore_index=True), pd.concat(real_parts, ignore_index=True)
    medsrc = pd.concat(med_parts, ignore_index=True)
    med = {"AlphaMissense": float(np.nanmedian(medsrc["am"])), "rmc_oe": float(np.nanmedian(medsrc["rmc_oe"])),
           "ccr_pct": float(np.nanmedian(medsrc["ccr_pct"]))}
    for t, cc in TOOLS.items():
        med[t] = float(np.nanmedian(medsrc[cc]))
    for t in AIM4V2_TOOLS:
        med[t] = float(np.nanmedian(medsrc[score_col(t)]))
    for k in AF_COLS:
        med[k] = float(np.nanmedian(medsrc[k]))
    del null_parts, real_parts, med_parts, medsrc
    Xr2, names2 = aim4_features_v2(real_df, med)
    Xn2, _ = aim4_features_v2(null_df, med)
    n1 = len(aim4_features(real_df.iloc[:1], med)[1])          # v1 features are the first n1 columns of v2
    X2 = np.vstack([Xr2, Xn2]); y = np.r_[np.ones(len(Xr2)), np.zeros(len(Xn2))]
    del null_df, real_df
    gc.collect()
    gb_kw = dict(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, l2_regularization=1.0,
                 early_stopping=True, validation_fraction=0.1, random_state=SEED)
    m_v1 = HistGradientBoostingClassifier(**gb_kw).fit(X2[:, :n1], y)
    m_gb = HistGradientBoostingClassifier(**gb_kw).fit(X2, y)
    sc = StandardScaler().fit(X2)
    Xs = sc.transform(X2)
    mlps = [MLPClassifier(hidden_layer_sizes=(256, 128), activation="relu", solver="adam", alpha=1e-4,
                          early_stopping=True, validation_fraction=0.1, max_iter=200, random_state=SEED + k).fit(Xs, y)
            for k in range(5)]
    del Xs
    c_gdx = len(syn_d[syn_d["centre"] == "GeneDx"]) / syn_mass
    pi = 1 - c_gdx * tot / len(Xr2)
    ratio = len(Xn2) / len(Xr2)
    log(f"Aim 4 v2: trained on {len(Xr2):,} GeneDx de novos vs {len(Xn2):,} null sites; v1 {n1} / v2 {len(names2)} "
        f"features; GBDT iterations v1 {m_v1.n_iter_}, v2 {m_gb.n_iter_}; MLP epochs {[m.n_iter_ for m in mlps]}; pi = {pi:.3f}")

    def logodds(p):
        p = np.clip(p, 1e-7, 1 - 1e-7)
        return np.log(p / (1 - p))

    def models_lo(X):
        lo1 = logodds(m_v1.predict_proba(X[:, :n1])[:, 1])
        lo2 = logodds(m_gb.predict_proba(X)[:, 1])
        Xs_ = sc.transform(X)
        lo3 = logodds(np.mean([m.predict_proba(Xs_)[:, 1] for m in mlps], 0))
        return {"Aim4v1": lo1, "Aim4v2_GBDT": lo2, "Aim4v2_MLP": lo3, "Aim4v2_ENS": (lo2 + lo3) / 2}

    def to_lr(lo):
        r = np.exp(lo) * ratio
        return np.maximum(0.0, (r - (1 - pi)) / pi)

    # pass 2: LR of every analysed site for each model; written to the column store and kept for evaluation.
    # Engineering (deviation log): columns are memory-mapped and read in 100k-row slices (the first run, which held a
    # whole chromosome's ~80 columns in memory, was killed for memory on chr7; training is unchanged and deterministic).
    n_real_, n_null_ = len(Xr2), len(Xn2)
    del X2, y, Xr2, Xn2
    gc.collect()
    store = {k: [] for k in list(M4_COLS) + ["am", "mr", "tx"]}
    drow, off = {}, 0
    for c in AUTOSOMES:
        od = P(a.fs, "col", f"chr{c}")
        mm = {n: np.load(P(od, n.replace("/", "_") + ".npy"), mmap_mode="r") for n in cols}
        n_c = len(mm["key"])
        parts = {k: [] for k in M4_COLS}
        for s0 in range(0, n_c, 100_000):
            u = pd.DataFrame({n: np.asarray(v[s0:s0 + 100_000]) for n, v in mm.items()})
            X, _ = aim4_features_v2(u, med)
            for k, lo in models_lo(X).items():
                parts[k].append(to_lr(lo).astype(np.float32))
            del X, u
        for k, col in M4_COLS.items():
            v = np.concatenate(parts[k])
            store[k].append(v)
            np.save(P(od, col + ".npy"), (np.log1p(np.minimum(v, M4_LRMAX)) / np.log1p(M4_LRMAX)).astype(np.float32))
        store["am"].append(np.asarray(mm["am"], np.float32)); store["mr"].append(np.asarray(mm[PRIMARY_MODEL], np.float32))
        store["tx"].append(np.asarray(mm["tx"], np.int32))
        kk = np.asarray(mm["key"])
        hit = np.flatnonzero(np.isin(kk, d_mis.loc[d_mis["chrom"] == c, "key"].to_numpy()))
        drow.update(dict(zip(kk[hit].tolist(), (hit + off).tolist())))
        off += n_c
        log(f"chr{c}: Aim 4 v1/v2 predictions written")
        del mm, parts, kk
        gc.collect()
    S = {k: np.concatenate(v) for k, v in store.items()}
    del store
    gc.collect()
    open(P(a.fs, "col", "_m4_done"), "w").close()
    nu = n_units("gene", a)
    CH = 10_000_000

    def mass(bins, nb, mask):
        flat = np.zeros(nu * nb)
        for s0 in range(0, len(bins), CH):
            sl = slice(s0, s0 + CH)
            ok = mask[sl] & (bins[sl] >= 0)
            flat += np.bincount(S["tx"][sl][ok].astype(np.int64) * nb + bins[sl][ok].astype(np.int64),
                                weights=S["mr"][sl][ok].astype(np.float64), minlength=nu * nb)
        return sparse.csr_matrix(flat.reshape(nu, nb))

    def binmap_of(x01):
        h = np.zeros(200001)
        for s0 in range(0, len(x01), CH):
            xs, w = x01[s0:s0 + CH], S["mr"][s0:s0 + CH]
            ok = np.isfinite(xs)
            h += np.bincount(np.clip(np.round(xs[ok].astype(np.float64) * 200000).astype(np.int64), 0, 200000),
                             weights=w[ok].astype(np.float64), minlength=200001)
        bm = BinMap(h)
        out = np.empty(len(x01), dtype=np.int16)
        for s0 in range(0, len(x01), CH):
            out[s0:s0 + CH] = bm.bin(x01[s0:s0 + CH])
        return out

    def dmap_of(dd, bins):
        r = dd["key"].map(drow)
        ok = r.notna().to_numpy()
        r = r[ok].to_numpy(np.int64)
        out = dd[ok][["p", "tx", "centre", "key"]].copy()
        out["unit"] = S["tx"][r].astype(np.int64); out["bin"] = bins[r]
        return out

    nb_lr = len(PRED_LR_EDGES) - 1
    qbins = {"AlphaMissense": binmap_of(S["am"])}
    lrbins = {}
    for k in M4_COLS:
        qbins[k] = binmap_of((np.log1p(np.minimum(S[k], M4_LRMAX)) / np.log1p(M4_LRMAX)).astype(np.float32))
        lrbins[k] = np.searchsorted(np.array(PRED_LR_EDGES[1:-1]), S[k], side="right").astype(np.int8)
        gc.collect()
    cold_site = cold_tx[S["tx"]]
    rows, top_rows, cmp_rows = [], [], []
    for centre in ("DDD", "RUMC", "GeneDx"):
        sp = syn_d.loc[syn_d["centre"] == centre, "p"].to_numpy(np.int64)
        dm_c = d_mis[d_mis["centre"] == centre]
        for gs_name, gmask in (("all_genes", np.ones(len(S["tx"]), bool)), ("cold_genes", cold_site)):
            if centre == "GeneDx" and gs_name != "cold_genes":
                continue                    # section 7: GeneDx training genes are in-sample; cold genes only
            for k in M4_COLS:
                M = mass(lrbins[k], nb_lr, gmask)
                lrc = Engine.lr(eng.run(M, syn_mass, sp, dmap_of(dm_c, lrbins[k])))
                sm = summ(lrc["LRm"])
                Mt = np.asarray(M.sum(0)).ravel()
                pred_mean = np.bincount(lrbins[k][gmask], weights=(S["mr"] * S[k])[gmask], minlength=nb_lr) / np.maximum(Mt, 1e-300)
                for i in range(nb_lr):
                    rows.append([centre, gs_name, k, PRED_LR_EDGES[i], PRED_LR_EDGES[i + 1], lrc["O"][0, i], lrc["E"][0, i],
                                 pred_mean[i], sm["est"][i], sm["lo95"][i], sm["hi95"][i]])
            reps = {}
            for nm, bins in qbins.items():
                r2 = eng.run(mass(bins, N_QBINS, gmask), syn_mass, sp, dmap_of(dm_c, bins))
                for pct in (1, 2, 5):
                    lr2 = Engine.lr(r2, [np.arange(N_QBINS - 10 * pct, N_QBINS)])
                    m_, w_ = summ(lr2["LRm"]), summ(lr2["LRw"])
                    top_rows.append([centre, gs_name, nm, pct, lr2["O"][0, 0], lr2["E"][0, 0], m_["est"][0], m_["lo95"][0],
                                     m_["hi95"][0], w_["est"][0], w_["lo95"][0], w_["hi95"][0]])
                    if pct == 2:
                        reps[nm] = lr2["LRw"][:, 0]
            for m_name, ref in (("Aim4v2_GBDT", "AlphaMissense"), ("Aim4v2_GBDT", "Aim4v1"), ("Aim4v2_MLP", "Aim4v2_GBDT"),
                                ("Aim4v2_ENS", "Aim4v2_GBDT"), ("Aim4v2_MLP", "AlphaMissense"), ("Aim4v2_ENS", "AlphaMissense"),
                                ("Aim4v1", "AlphaMissense")):
                with np.errstate(invalid="ignore", divide="ignore"):
                    dl = np.log(reps[m_name]) - np.log(reps[ref])
                b = dl[1:][np.isfinite(dl[1:])]
                lo, hi = (np.percentile(b, [2.5, 97.5]) if len(b) else (np.nan, np.nan))
                cmp_rows.append([centre, gs_name, m_name, ref, float(np.exp(dl[0])), float(np.exp(lo)), float(np.exp(hi)),
                                 bool(lo > 0 or hi < 0), len(b)])
    write(rows, ["centre", "genes", "model", "predLR_lo", "predLR_hi", "O", "E", "predicted_LR_mean", "empirical_LR", "lo95",
                 "hi95"], P(a.results, "aim4v2_calibration.tsv"))
    write(top_rows, ["centre", "genes", "score", "top_pct", "O", "E", "LRm", "lo95", "hi95", "LRw", "lo95_w", "hi95_w"],
          P(a.results, "aim4v2_top.tsv"))
    write(cmp_rows, ["centre", "genes", "model", "reference", "ratio_LRw_top2", "lo95", "hi95", "ci_excludes_1", "n_boot_finite"],
          P(a.results, "aim4v2_compare.tsv"))
    with open(P(a.results, "aim4v2_model.json"), "w") as f:
        json.dump({"n_real": int(n_real_), "n_null": int(n_null_), "pi_genedx_train": float(pi),
                   "gbdt_iter": {"v1": int(m_v1.n_iter_), "v2": int(m_gb.n_iter_)}, "mlp_epochs": [int(m.n_iter_) for m in mlps],
                   "features_v2": names2, "n_features_v1": n1, "medians": med}, f, indent=1)
    cmp = pd.DataFrame(cmp_rows, columns=["centre", "genes", "model", "reference", "ratio", "lo", "hi", "sig", "nb"])
    prim = cmp[(cmp.centre == "DDD") & (cmp.genes == "cold_genes") & (cmp.model == "Aim4v2_GBDT") & (cmp.reference == "AlphaMissense")]
    txt = ["PRIMARY (A5): DDD cold genes, within-gene top-2% LR, v2-GBDT / AlphaMissense = "
           f"{prim['ratio'].iloc[0]:.2f} [{prim['lo'].iloc[0]:.2f}, {prim['hi'].iloc[0]:.2f}] -> "
           f"{'v2-GBDT BEATS AlphaMissense' if prim['lo'].iloc[0] > 1 else ('AlphaMissense better' if prim['hi'].iloc[0] < 1 else 'no significant difference')}"]
    for r in cmp.itertuples():
        txt.append(f"{r.centre:5s} {r.genes:10s} {r.model:12s} / {r.reference:13s} {r.ratio:.2f} [{r.lo:.2f}, {r.hi:.2f}]")
    tp = pd.DataFrame(top_rows, columns=["centre", "genes", "score", "pct", "O", "E", "LRm", "lo", "hi", "LRw", "lo_w", "hi_w"])
    for r in tp[(tp.centre == "DDD") & (tp.pct == 2)].itertuples():
        txt.append(f"DDD {r.genes:10s} top2% {r.score:13s} marginal {r.LRm:.2f} within-gene {r.LRw:.2f} [{r.lo_w:.2f}, {r.hi_w:.2f}]")
    open(P(a.results, "aim4v2_summary.txt"), "w").write("\n".join(txt) + "\n")
    log("Aim 4 v2:\n" + "\n".join(txt))


M4V3 = {"Aim4v3_GBDT": "m4_v3_gbdt"}                     # amendment A7: gene-matched null


def cmd_aim4v3(a):
    """Plan section 3h (A7): v2-GBDT with a gene-matched null (10 MR-weighted sites from the same transcript per real
    de novo). Evaluation as A5, paired against AlphaMissense and v2-GBDT; then the A1 nested levels on DDD."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    os.makedirs(a.results, exist_ok=True)
    d, eng, syn_mass, syn_p_all, d_mis, tinfo = setup(a)
    syn_d = dnv_filter(d[d["status"] == "synonymous"])
    cold_tx = tinfo["gene"].map(gene_split).to_numpy(bool)
    rng = np.random.default_rng(SEED)
    cols = list(dict.fromkeys(AIM4_COLS + [score_col(t) for t in AIM4V2_TOOLS] + AF_COLS))
    gdx = d_mis[(d_mis["centre"] == "GeneDx") & ~cold_tx[d_mis["tx"].to_numpy()]]
    real_parts, null_parts, med_parts = [], [], []
    for c in AUTOSOMES:
        od = P(a.fs, "col", f"chr{c}")
        mm = {n: np.load(P(od, n.replace("/", "_") + ".npy"), mmap_mode="r") for n in cols}
        key, tx, mr = np.asarray(mm["key"]), np.asarray(mm["tx"]).astype(np.int64), np.asarray(mm[PRIMARY_MODEL], float)
        rk = gdx.loc[gdx["chrom"] == c, "key"].to_numpy()
        ridx = np.flatnonzero(np.isin(key, rk))
        if len(ridx) == 0:
            continue
        order = np.argsort(tx, kind="stable")
        starts = np.searchsorted(tx[order], np.arange(tx.max() + 2))
        nidx = []
        for t, n_t in zip(*np.unique(tx[ridx], return_counts=True)):
            members = order[starts[t]:starts[t + 1]]
            w = mr[members]
            nidx.append(members[rng.choice(len(members), size=10 * n_t, replace=True, p=w / w.sum())])
        nidx = np.concatenate(nidx)
        tri = np.flatnonzero(~cold_tx[tx])
        smp = tri[pd.Series(np.arange(len(tri))).sample(n=min(len(tri), 20000), random_state=SEED).to_numpy()]
        for parts, ix in ((real_parts, ridx), (null_parts, nidx), (med_parts, smp)):
            parts.append(pd.DataFrame({n: np.asarray(v)[ix] for n, v in mm.items()}))
        del mm, key, tx, mr
        gc.collect()
    real_df, null_df = pd.concat(real_parts, ignore_index=True), pd.concat(null_parts, ignore_index=True)
    medsrc = pd.concat(med_parts, ignore_index=True)
    med = {"AlphaMissense": float(np.nanmedian(medsrc["am"])), "rmc_oe": float(np.nanmedian(medsrc["rmc_oe"])),
           "ccr_pct": float(np.nanmedian(medsrc["ccr_pct"]))}
    for t, cc in TOOLS.items():
        med[t] = float(np.nanmedian(medsrc[cc]))
    for t in AIM4V2_TOOLS:
        med[t] = float(np.nanmedian(medsrc[score_col(t)]))
    for k in AF_COLS:
        med[k] = float(np.nanmedian(medsrc[k]))
    Xr, names = aim4_features_v2(real_df, med)
    Xn, _ = aim4_features_v2(null_df, med)
    n_real_, n_null_ = len(Xr), len(Xn)
    X = np.vstack([Xr, Xn]); y = np.r_[np.ones(n_real_), np.zeros(n_null_)]
    del real_df, null_df, medsrc, Xr, Xn, real_parts, null_parts, med_parts
    gc.collect()
    m3 = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, l2_regularization=1.0,
                                        early_stopping=True, validation_fraction=0.1, random_state=SEED).fit(X, y)
    del X, y
    gc.collect()
    tot = 0.0
    for c in AUTOSOMES:
        tx = np.load(P(a.fs, "col", f"chr{c}", "tx.npy")).astype(np.int64)
        tot += float(np.load(P(a.fs, "col", f"chr{c}", PRIMARY_MODEL + ".npy"))[~cold_tx[tx]].sum())
    pi = 1 - (len(syn_d[syn_d["centre"] == "GeneDx"]) / syn_mass) * tot / n_real_
    ratio = n_null_ / n_real_
    log(f"Aim 4 v3 (gene-matched null): {n_real_:,} GeneDx de novos vs {n_null_:,} same-transcript null sites; "
        f"{len(names)} features; GBDT iterations {m3.n_iter_}; pi = {pi:.3f}")

    def to_lr(X_):
        p = np.clip(m3.predict_proba(X_)[:, 1], 1e-7, 1 - 1e-7)
        return np.maximum(0.0, (p / (1 - p) * ratio - (1 - pi)) / pi)

    S = {k: [] for k in ("v3", "v2", "am", "mr", "tx")}
    drow, off = {}, 0
    for c in AUTOSOMES:
        od = P(a.fs, "col", f"chr{c}")
        mm = {n: np.load(P(od, n.replace("/", "_") + ".npy"), mmap_mode="r") for n in cols}
        n_c = len(mm["key"])
        v = np.concatenate([to_lr(aim4_features_v2(pd.DataFrame({n: np.asarray(x[s0:s0 + 100_000]) for n, x in mm.items()}),
                                                   med)[0]).astype(np.float32) for s0 in range(0, n_c, 100_000)])
        np.save(P(od, M4V3["Aim4v3_GBDT"] + ".npy"), (np.log1p(np.minimum(v, M4_LRMAX)) / np.log1p(M4_LRMAX)).astype(np.float32))
        S["v3"].append(v)
        S["v2"].append(np.load(P(od, M4_COLS["Aim4v2_GBDT"] + ".npy")))
        S["am"].append(np.asarray(mm["am"], np.float32)); S["mr"].append(np.asarray(mm[PRIMARY_MODEL], np.float32))
        S["tx"].append(np.asarray(mm["tx"], np.int32))
        kk = np.asarray(mm["key"])
        hit = np.flatnonzero(np.isin(kk, d_mis.loc[d_mis["chrom"] == c, "key"].to_numpy()))
        drow.update(dict(zip(kk[hit].tolist(), (hit + off).tolist())))
        off += n_c
        log(f"chr{c}: Aim 4 v3 predictions written")
        del mm, kk
        gc.collect()
    S = {k: np.concatenate(v) for k, v in S.items()}
    gc.collect()
    open(P(a.fs, "col", "_m4v3_done"), "w").close()
    nu = n_units("gene", a)
    CH = 10_000_000

    def mass(bins, nb, mask):
        flat = np.zeros(nu * nb)
        for s0 in range(0, len(bins), CH):
            sl = slice(s0, s0 + CH)
            ok = mask[sl] & (bins[sl] >= 0)
            flat += np.bincount(S["tx"][sl][ok].astype(np.int64) * nb + bins[sl][ok].astype(np.int64),
                                weights=S["mr"][sl][ok].astype(np.float64), minlength=nu * nb)
        return sparse.csr_matrix(flat.reshape(nu, nb))

    def binmap_of(x01):
        h = np.zeros(200001)
        for s0 in range(0, len(x01), CH):
            xs, w = x01[s0:s0 + CH], S["mr"][s0:s0 + CH]
            ok = np.isfinite(xs)
            h += np.bincount(np.clip(np.round(xs[ok].astype(np.float64) * 200000).astype(np.int64), 0, 200000),
                             weights=w[ok].astype(np.float64), minlength=200001)
        bm = BinMap(h)
        out = np.empty(len(x01), dtype=np.int16)
        for s0 in range(0, len(x01), CH):
            out[s0:s0 + CH] = bm.bin(x01[s0:s0 + CH])
        return out

    def dmap_of(dd, bins):
        r = dd["key"].map(drow)
        ok = r.notna().to_numpy()
        r = r[ok].to_numpy(np.int64)
        out = dd[ok][["p", "tx", "centre", "key"]].copy()
        out["unit"] = S["tx"][r].astype(np.int64); out["bin"] = bins[r]
        return out

    qb = {"AlphaMissense": binmap_of(S["am"]), "Aim4v2_GBDT": binmap_of(S["v2"]),
          "Aim4v3_GBDT": binmap_of((np.log1p(np.minimum(S["v3"], M4_LRMAX)) / np.log1p(M4_LRMAX)).astype(np.float32))}
    nb_lr = len(PRED_LR_EDGES) - 1
    lrb = np.searchsorted(np.array(PRED_LR_EDGES[1:-1]), S["v3"], side="right").astype(np.int8)
    cold_site = cold_tx[S["tx"]]
    rows, top_rows, cmp_rows = [], [], []
    for centre in ("DDD", "RUMC", "GeneDx"):
        sp = syn_d.loc[syn_d["centre"] == centre, "p"].to_numpy(np.int64)
        dm_c = d_mis[d_mis["centre"] == centre]
        for gs_name, gmask in (("all_genes", np.ones(len(S["tx"]), bool)), ("cold_genes", cold_site)):
            if centre == "GeneDx" and gs_name != "cold_genes":
                continue
            M = mass(lrb, nb_lr, gmask)
            lrc = Engine.lr(eng.run(M, syn_mass, sp, dmap_of(dm_c, lrb)))
            sm = summ(lrc["LRm"])
            Mt = np.asarray(M.sum(0)).ravel()
            pm = np.bincount(lrb[gmask], weights=(S["mr"] * S["v3"])[gmask], minlength=nb_lr) / np.maximum(Mt, 1e-300)
            for i in range(nb_lr):
                rows.append([centre, gs_name, "Aim4v3_GBDT", PRED_LR_EDGES[i], PRED_LR_EDGES[i + 1], lrc["O"][0, i],
                             lrc["E"][0, i], pm[i], sm["est"][i], sm["lo95"][i], sm["hi95"][i]])
            reps = {}
            for nm, bins in qb.items():
                r2 = eng.run(mass(bins, N_QBINS, gmask), syn_mass, sp, dmap_of(dm_c, bins))
                for pct in (1, 2, 5):
                    lr2 = Engine.lr(r2, [np.arange(N_QBINS - 10 * pct, N_QBINS)])
                    m_, w_ = summ(lr2["LRm"]), summ(lr2["LRw"])
                    top_rows.append([centre, gs_name, nm, pct, lr2["O"][0, 0], lr2["E"][0, 0], m_["est"][0], m_["lo95"][0],
                                     m_["hi95"][0], w_["est"][0], w_["lo95"][0], w_["hi95"][0]])
                    if pct == 2:
                        reps[nm] = lr2["LRw"][:, 0]
            for m_name, ref in (("Aim4v3_GBDT", "AlphaMissense"), ("Aim4v3_GBDT", "Aim4v2_GBDT")):
                with np.errstate(invalid="ignore", divide="ignore"):
                    dl = np.log(reps[m_name]) - np.log(reps[ref])
                b = dl[1:][np.isfinite(dl[1:])]
                lo, hi = (np.percentile(b, [2.5, 97.5]) if len(b) else (np.nan, np.nan))
                cmp_rows.append([centre, gs_name, m_name, ref, float(np.exp(dl[0])), float(np.exp(lo)), float(np.exp(hi)),
                                 bool(lo > 0 or hi < 0), len(b)])
    write(rows, ["centre", "genes", "model", "predLR_lo", "predLR_hi", "O", "E", "predicted_LR_mean", "empirical_LR", "lo95",
                 "hi95"], P(a.results, "aim4v3_calibration.tsv"))
    write(top_rows, ["centre", "genes", "score", "top_pct", "O", "E", "LRm", "lo95", "hi95", "LRw", "lo95_w", "hi95_w"],
          P(a.results, "aim4v3_top.tsv"))
    write(cmp_rows, ["centre", "genes", "model", "reference", "ratio_LRw_top2", "lo95", "hi95", "ci_excludes_1", "n_boot_finite"],
          P(a.results, "aim4v3_compare.tsv"))
    with open(P(a.results, "aim4v3_model.json"), "w") as f:
        json.dump({"n_real": int(n_real_), "n_null": int(n_null_), "pi_genedx_train": float(pi), "gbdt_iter": int(m3.n_iter_),
                   "features": names, "medians": med}, f, indent=1)
    cmp = pd.DataFrame(cmp_rows, columns=["centre", "genes", "model", "reference", "ratio", "lo", "hi", "sig", "nb"])
    prim = cmp[(cmp.centre == "DDD") & (cmp.genes == "cold_genes") & (cmp.reference == "AlphaMissense")].iloc[0]
    txt = [f"A7 comparison 1: DDD cold genes, within-gene top-2% LR, v3 / AlphaMissense = {prim.ratio:.2f} [{prim.lo:.2f}, "
           f"{prim.hi:.2f}] -> {'v3 BETTER' if prim.lo > 1 else ('AlphaMissense better' if prim.hi < 1 else 'no significant difference')}"]
    for r in cmp.itertuples():
        txt.append(f"{r.centre:6s} {r.genes:10s} {r.model} / {r.reference:13s} {r.ratio:.2f} [{r.lo:.2f}, {r.hi:.2f}]")
    tp = pd.DataFrame(top_rows, columns=["centre", "genes", "score", "pct", "O", "E", "LRm", "lo", "hi", "LRw", "lo_w", "hi_w"])
    for r in tp[(tp.centre == "DDD") & (tp.pct == 2)].itertuples():
        txt.append(f"DDD {r.genes:10s} top2% {r.score:13s} marginal {r.LRm:.2f} within-gene {r.LRw:.2f} [{r.lo_w:.2f}, {r.hi_w:.2f}]")
    open(P(a.results, "aim4v3_summary.txt"), "w").write("\n".join(txt) + "\n")
    log("Aim 4 v3:\n" + "\n".join(txt))
    del S
    gc.collect()
    # nested levels A-E for v3 on DDD de novos (comparison 4)
    a.centre_only = "DDD"
    run_scores(a, ["AlphaMissense", "Aim4v3_GBDT"], A1_KINDS, ["all_genes"], "m4v3nested_ddd")
    rep_report(a, float("nan"), tag="m4v3nested_ddd", out_name="m4v3nested_ddd_summary.txt")
    a.centre_only = None


M4V4 = {"Aim4v4_GBDT": "m4_v4_gbdt"}                     # amendment A8: gene-matched, all DD centres, autism-only test


def cmd_aim4v4(a):
    """Plan section 3i (A8): v3's model and gene-matched null, trained on all three DD centres and all genes; evaluated
    only in the autism cohort (and siblings), paired against AlphaMissense and v3. One shot, reported regardless."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    os.makedirs(a.results, exist_ok=True)
    d, eng, syn_mass, syn_p_all, d_mis, tinfo = setup(a)
    syn_d = dnv_filter(d[d["status"] == "synonymous"])
    rng = np.random.default_rng(SEED)
    cols = list(dict.fromkeys(AIM4_COLS + [score_col(t) for t in AIM4V2_TOOLS] + AF_COLS))
    real_parts, null_parts, med_parts = [], [], []
    tot = 0.0
    for c in AUTOSOMES:
        od = P(a.fs, "col", f"chr{c}")
        mm = {n: np.load(P(od, n.replace("/", "_") + ".npy"), mmap_mode="r") for n in cols}
        key, tx, mr = np.asarray(mm["key"]), np.asarray(mm["tx"]).astype(np.int64), np.asarray(mm[PRIMARY_MODEL], float)
        tot += float(mr.sum())
        ridx = np.flatnonzero(np.isin(key, d_mis.loc[d_mis["chrom"] == c, "key"].to_numpy()))
        smp = pd.Series(np.arange(len(key))).sample(n=min(len(key), 20000), random_state=SEED).to_numpy()
        med_parts.append(pd.DataFrame({n: np.asarray(v)[smp] for n, v in mm.items()}))
        if len(ridx):
            order = np.argsort(tx, kind="stable")
            starts = np.searchsorted(tx[order], np.arange(tx.max() + 2))
            nidx = []
            for t, n_t in zip(*np.unique(tx[ridx], return_counts=True)):
                members = order[starts[t]:starts[t + 1]]
                w = mr[members]
                nidx.append(members[rng.choice(len(members), size=10 * n_t, replace=True, p=w / w.sum())])
            nidx = np.concatenate(nidx)
            real_parts.append(pd.DataFrame({n: np.asarray(v)[ridx] for n, v in mm.items()}))
            null_parts.append(pd.DataFrame({n: np.asarray(v)[nidx] for n, v in mm.items()}))
        del mm, key, tx, mr
        gc.collect()
    real_df, null_df = pd.concat(real_parts, ignore_index=True), pd.concat(null_parts, ignore_index=True)
    medsrc = pd.concat(med_parts, ignore_index=True)
    med = {"AlphaMissense": float(np.nanmedian(medsrc["am"])), "rmc_oe": float(np.nanmedian(medsrc["rmc_oe"])),
           "ccr_pct": float(np.nanmedian(medsrc["ccr_pct"]))}
    for t, cc in TOOLS.items():
        med[t] = float(np.nanmedian(medsrc[cc]))
    for t in AIM4V2_TOOLS:
        med[t] = float(np.nanmedian(medsrc[score_col(t)]))
    for k in AF_COLS:
        med[k] = float(np.nanmedian(medsrc[k]))
    Xr, names = aim4_features_v2(real_df, med)
    Xn, _ = aim4_features_v2(null_df, med)
    n_real_, n_null_ = len(Xr), len(Xn)
    X = np.vstack([Xr, Xn]); y = np.r_[np.ones(n_real_), np.zeros(n_null_)]
    del real_df, null_df, medsrc, Xr, Xn, real_parts, null_parts, med_parts
    gc.collect()
    m4 = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, l2_regularization=1.0,
                                        early_stopping=True, validation_fraction=0.1, random_state=SEED).fit(X, y)
    del X, y
    gc.collect()
    pi = 1 - (len(syn_d) / syn_mass) * tot / n_real_
    ratio = n_null_ / n_real_
    log(f"Aim 4 v4 (gene-matched, all DD centres): {n_real_:,} de novo sites vs {n_null_:,} same-transcript null sites; "
        f"{len(names)} features; GBDT iterations {m4.n_iter_}; pi = {pi:.3f}")
    for c in AUTOSOMES:
        od = P(a.fs, "col", f"chr{c}")
        mm = {n: np.load(P(od, n.replace("/", "_") + ".npy"), mmap_mode="r") for n in cols}
        n_c = len(mm["key"])
        parts = []
        for s0 in range(0, n_c, 100_000):
            Xs, _ = aim4_features_v2(pd.DataFrame({n: np.asarray(x[s0:s0 + 100_000]) for n, x in mm.items()}), med)
            p = np.clip(m4.predict_proba(Xs)[:, 1], 1e-7, 1 - 1e-7)
            parts.append(np.maximum(0.0, (p / (1 - p) * ratio - (1 - pi)) / pi).astype(np.float32))
        v = np.concatenate(parts)
        np.save(P(od, M4V4["Aim4v4_GBDT"] + ".npy"), (np.log1p(np.minimum(v, M4_LRMAX)) / np.log1p(M4_LRMAX)).astype(np.float32))
        log(f"chr{c}: Aim 4 v4 predictions written")
        del mm, parts
        gc.collect()
    open(P(a.fs, "col", "_m4v4_done"), "w").close()
    with open(P(a.results, "aim4v4_model.json"), "w") as f:
        json.dump({"n_real": int(n_real_), "n_null": int(n_null_), "pi_all_dd": float(pi), "gbdt_iter": int(m4.n_iter_),
                   "features": names, "medians": med}, f, indent=1)
    # evaluation: autism probands and siblings only (A2 / A6.3 specification)
    base = a.results
    pos_setup(a)
    for pheno, sub, tag in (("Affected", "replication_affected", "m4v4rep"), ("Unaffected", "replication_unaffected", "m4v4sib")):
        a.results, a.cohort = P(base, sub), pheno
        os.makedirs(a.results, exist_ok=True)
        run_scores(a, ["AlphaMissense", "Aim4v3_GBDT", "Aim4v4_GBDT"], A1_KINDS, ["all_genes"], tag, restrict="aim1_common")
        if os.path.exists(P(a.results, f"{tag}_top2pct.tsv")):
            rep_report(a, float("nan"), tag=tag, out_name=f"{tag}_summary.txt")
    a.results, a.cohort = base, None
    txt = []
    for tag, sub in (("m4v4rep", "replication_affected"), ("m4v4sib", "replication_unaffected")):
        if not (os.path.exists(P(base, sub, f"{tag}_top2reps.npz")) and os.path.exists(P(base, sub, f"{tag}_top2pct.tsv"))):
            txt.append(f"{tag}: no top-2% table (no de novos in the top group)")
            continue
        z = dict(np.load(P(base, sub, f"{tag}_top2reps.npz")))
        for m_name, ref in (("Aim4v4_GBDT", "AlphaMissense"), ("Aim4v4_GBDT", "Aim4v3_GBDT"), ("Aim4v3_GBDT", "AlphaMissense")):
            x, y_ = z.get(f"{m_name}|gene|all_genes|LRw"), z.get(f"{ref}|gene|all_genes|LRw")
            if x is None or y_ is None:
                continue
            with np.errstate(all="ignore"):
                dl = np.log(x) - np.log(y_)
            b = dl[1:][np.isfinite(dl[1:])]
            lo, hi = np.exp(np.percentile(b, [2.5, 97.5])) if len(b) > 20 else (np.nan, np.nan)
            verdict = "v4 better" if lo > 1 else ("reference better" if hi < 1 else "no significant difference")
            prim = " [A8 PRIMARY]" if (tag == "m4v4rep" and m_name == "Aim4v4_GBDT" and ref == "AlphaMissense") else ""
            txt.append(f"{tag}: within-gene top-2% LR {m_name} / {ref} = {np.exp(dl[0]):.2f} [{lo:.2f}, {hi:.2f}] -> {verdict}{prim}")
        t = pd.read_csv(P(base, sub, f"{tag}_top2pct.tsv"), sep="\t")
        for r in t[t["unit"].isin(["gene", "win50", "win50_rmc"])].itertuples():
            txt.append(f"{tag} {r.score} {r.unit}: marginal {r.LRm_top2:.2f} (lower bound {r.lo_sys:.2f}); within {r.LRw_top2:.2f} "
                       f"[{r.lo95_w:.2f}, {r.hi95_w:.2f}] (lower bound {r.lo_sys_w:.2f})")
    open(P(base, "aim4v4_summary.txt"), "w").write("\n".join(txt) + "\n")
    log("Aim 4 v4:\n" + "\n".join(txt))



M4V5 = {"Aim4v5_GBDT": "m4_v5_gbdt"}                     # amendment A15: label-free inputs only
V5_KEEP = {"SIFT", "SIFT4G", "LRT", "MutationAssessor", "PROVEAN", "PrimateAI", "LIST-S2", "CADD", "DANN", "Eigen", "Eigen-PC",
           "GERP++", "phyloP100way", "phyloP30way", "phastCons100way", "SiPhy", "GenoCanyon", "fitCons", "ESM1b", "MisFit_D",
           "MisFit_S", "GPN-MSA", "GERP_92m", "phyloP470way", "phastCons470way", "AlphaMissense", "rmc_oe", "ccr_pct",
           "af_plddt", "af_nca10", "af_plddtw"}


def aim4_features_v5(u, med):
    X, names = aim4_features_v2(u, med)
    keep = [k for k, n in enumerate(names)
            if n.startswith(("ref_", "alt_")) or (n[:-8] if n.endswith("_missing") else n) in V5_KEEP]
    return X[:, keep], [names[k] for k in keep]


def cmd_aim4v5(a):
    """Amendment A15: v4 with label-free inputs only (V5_KEEP); same null, data, settings and autism evaluation."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    os.makedirs(a.results, exist_ok=True)
    d, eng, syn_mass, syn_p_all, d_mis, tinfo = setup(a)
    syn_d = dnv_filter(d[d["status"] == "synonymous"])
    rng = np.random.default_rng(SEED)
    cols = list(dict.fromkeys(AIM4_COLS + [score_col(t) for t in AIM4V2_TOOLS] + AF_COLS))
    real_parts, null_parts, med_parts = [], [], []
    tot = 0.0
    for c in AUTOSOMES:
        od = P(a.fs, "col", f"chr{c}")
        mm = {n: np.load(P(od, n.replace("/", "_") + ".npy"), mmap_mode="r") for n in cols}
        key, tx, mr = np.asarray(mm["key"]), np.asarray(mm["tx"]).astype(np.int64), np.asarray(mm[PRIMARY_MODEL], float)
        tot += float(mr.sum())
        ridx = np.flatnonzero(np.isin(key, d_mis.loc[d_mis["chrom"] == c, "key"].to_numpy()))
        smp = pd.Series(np.arange(len(key))).sample(n=min(len(key), 20000), random_state=SEED).to_numpy()
        med_parts.append(pd.DataFrame({n: np.asarray(v)[smp] for n, v in mm.items()}))
        if len(ridx):
            order = np.argsort(tx, kind="stable")
            starts = np.searchsorted(tx[order], np.arange(tx.max() + 2))
            nidx = []
            for t, n_t in zip(*np.unique(tx[ridx], return_counts=True)):
                members = order[starts[t]:starts[t + 1]]
                w = mr[members]
                nidx.append(members[rng.choice(len(members), size=10 * n_t, replace=True, p=w / w.sum())])
            nidx = np.concatenate(nidx)
            real_parts.append(pd.DataFrame({n: np.asarray(v)[ridx] for n, v in mm.items()}))
            null_parts.append(pd.DataFrame({n: np.asarray(v)[nidx] for n, v in mm.items()}))
        del mm, key, tx, mr
        gc.collect()
    real_df, null_df = pd.concat(real_parts, ignore_index=True), pd.concat(null_parts, ignore_index=True)
    medsrc = pd.concat(med_parts, ignore_index=True)
    med = {"AlphaMissense": float(np.nanmedian(medsrc["am"])), "rmc_oe": float(np.nanmedian(medsrc["rmc_oe"])),
           "ccr_pct": float(np.nanmedian(medsrc["ccr_pct"]))}
    for t, cc in TOOLS.items():
        med[t] = float(np.nanmedian(medsrc[cc]))
    for t in AIM4V2_TOOLS:
        med[t] = float(np.nanmedian(medsrc[score_col(t)]))
    for k in AF_COLS:
        med[k] = float(np.nanmedian(medsrc[k]))
    Xr, names = aim4_features_v5(real_df, med)
    Xn, _ = aim4_features_v5(null_df, med)
    n_real_, n_null_ = len(Xr), len(Xn)
    X = np.vstack([Xr, Xn]); y = np.r_[np.ones(n_real_), np.zeros(n_null_)]
    del real_df, null_df, medsrc, Xr, Xn, real_parts, null_parts, med_parts
    gc.collect()
    m4 = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, l2_regularization=1.0,
                                        early_stopping=True, validation_fraction=0.1, random_state=SEED).fit(X, y)
    del X, y
    gc.collect()
    pi = 1 - (len(syn_d) / syn_mass) * tot / n_real_
    ratio = n_null_ / n_real_
    log(f"Aim 4 v5 (gene-matched, all DD centres): {n_real_:,} de novo sites vs {n_null_:,} same-transcript null sites; "
        f"{len(names)} features; GBDT iterations {m4.n_iter_}; pi = {pi:.3f}")
    for c in AUTOSOMES:
        od = P(a.fs, "col", f"chr{c}")
        mm = {n: np.load(P(od, n.replace("/", "_") + ".npy"), mmap_mode="r") for n in cols}
        n_c = len(mm["key"])
        parts = []
        for s0 in range(0, n_c, 100_000):
            Xs, _ = aim4_features_v5(pd.DataFrame({n: np.asarray(x[s0:s0 + 100_000]) for n, x in mm.items()}), med)
            p = np.clip(m4.predict_proba(Xs)[:, 1], 1e-7, 1 - 1e-7)
            parts.append(np.maximum(0.0, (p / (1 - p) * ratio - (1 - pi)) / pi).astype(np.float32))
        v = np.concatenate(parts)
        np.save(P(od, M4V5["Aim4v5_GBDT"] + ".npy"), (np.log1p(np.minimum(v, M4_LRMAX)) / np.log1p(M4_LRMAX)).astype(np.float32))
        log(f"chr{c}: Aim 4 v5 predictions written")
        del mm, parts
        gc.collect()
    open(P(a.fs, "col", "_m4v5_done"), "w").close()
    with open(P(a.results, "aim4v5_model.json"), "w") as f:
        json.dump({"n_real": int(n_real_), "n_null": int(n_null_), "pi_all_dd_v5": float(pi), "gbdt_iter": int(m4.n_iter_),
                   "features": names, "medians": med}, f, indent=1)
    # evaluation: autism probands and siblings only (A2 / A6.3 specification)
    base = a.results
    pos_setup(a)
    for pheno, sub, tag in (("Affected", "replication_affected", "m4v5rep"), ("Unaffected", "replication_unaffected", "m4v5sib")):
        a.results, a.cohort = P(base, sub), pheno
        os.makedirs(a.results, exist_ok=True)
        run_scores(a, ["AlphaMissense", "Aim4v4_GBDT", "Aim4v5_GBDT"], A1_KINDS, ["all_genes"], tag, restrict="aim1_common")
        if os.path.exists(P(a.results, f"{tag}_top2pct.tsv")):
            rep_report(a, float("nan"), tag=tag, out_name=f"{tag}_summary.txt")
    a.results, a.cohort = base, None
    txt = []
    for tag, sub in (("m4v5rep", "replication_affected"), ("m4v5sib", "replication_unaffected")):
        if not (os.path.exists(P(base, sub, f"{tag}_top2reps.npz")) and os.path.exists(P(base, sub, f"{tag}_top2pct.tsv"))):
            txt.append(f"{tag}: no top-2% table (no de novos in the top group)")
            continue
        z = dict(np.load(P(base, sub, f"{tag}_top2reps.npz")))
        for m_name, ref in (("Aim4v5_GBDT", "AlphaMissense"), ("Aim4v5_GBDT", "Aim4v4_GBDT"), ("Aim4v4_GBDT", "AlphaMissense")):
            x, y_ = z.get(f"{m_name}|gene|all_genes|LRw"), z.get(f"{ref}|gene|all_genes|LRw")
            if x is None or y_ is None:
                continue
            with np.errstate(all="ignore"):
                dl = np.log(x) - np.log(y_)
            b = dl[1:][np.isfinite(dl[1:])]
            lo, hi = np.exp(np.percentile(b, [2.5, 97.5])) if len(b) > 20 else (np.nan, np.nan)
            verdict = "v5 better" if lo > 1 else ("reference better" if hi < 1 else "no significant difference")
            prim = " [A15.2]" if (tag == "m4v5rep" and m_name == "Aim4v5_GBDT" and ref == "AlphaMissense") else ""
            txt.append(f"{tag}: within-gene top-2% LR {m_name} / {ref} = {np.exp(dl[0]):.2f} [{lo:.2f}, {hi:.2f}] -> {verdict}{prim}")
        t = pd.read_csv(P(base, sub, f"{tag}_top2pct.tsv"), sep="\t")
        for r in t[t["unit"].isin(["gene", "win50", "win50_rmc"])].itertuples():
            txt.append(f"{tag} {r.score} {r.unit}: marginal {r.LRm_top2:.2f} (lower bound {r.lo_sys:.2f}); within {r.LRw_top2:.2f} "
                       f"[{r.lo95_w:.2f}, {r.hi95_w:.2f}] (lower bound {r.lo_sys_w:.2f})")
    open(P(base, "aim4v5_summary.txt"), "w").write("\n".join(txt) + "\n")
    log("Aim 4 v5:\n" + "\n".join(txt))


PUBLIC3 = P(EXT, "public3", "public3_dnv_snv_grch37.tsv")
P3_STUDIES = ["Jin2017", "Fromer2014", "Gulsuner2013", "McCarthy2014", "Willsey2017", "epi4k2013", "Rauch2012"]


def _paired(z, m_name, ref, st="all_genes"):
    x, y_ = z.get(f"{m_name}|gene|{st}|LRw"), z.get(f"{ref}|gene|{st}|LRw")
    if x is None or y_ is None:
        return None
    with np.errstate(all="ignore"):
        dl = np.log(x) - np.log(y_)
    b = dl[1:][np.isfinite(dl[1:])]
    lo, hi = np.exp(np.percentile(b, [2.5, 97.5])) if len(b) > 20 else (np.nan, np.nan)
    return float(np.exp(dl[0])), float(lo), float(hi)


def cmd_public3(a):
    """Plan section 3j (A9): third, public, multi-phenotype cohort. Models first (descriptive), then the A2 readouts for
    AlphaMissense + 36 tools (power-checked), then the Jin 2017 NDD-flagged subset (descriptive)."""
    base = a.results
    src = pd.read_csv(PUBLIC3, sep="\t", dtype=str)
    fmt = pd.DataFrame({"Cohort": src["study"], "IID": src["sample"], "Pheno": "Affected", "Chrom": src["chrom"],
                        "Position": src["pos37"], "Ref": src["ref"], "Alt": src["alt"]})
    ndd = (src["study"] == "Jin2017") & (src["ndd_flag"].str.strip().str.lower() == "yes")
    fmt_ndd = fmt[ndd].assign(Cohort="Jin2017_NDD")
    fp = P(EXT, "public3", "public3_loader_format.tsv")
    pd.concat([fmt, fmt_ndd], ignore_index=True).to_csv(fp, sep="\t", index=False)
    a.rep_file, a.dedup_zhou, a.cohort = fp, True, "Affected"
    a.results = P(base, "public3")
    os.makedirs(a.results, exist_ok=True)
    pos_setup(a)
    out = []
    for subset, cohorts, tag in (("all", P3_STUDIES, "p3"), ("Jin2017 NDD=Yes", ["Jin2017_NDD"], "p3ndd")):
        a.rep_cohorts = cohorts
        d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a)
        col = Collect(a, ["AM_gate"], ["gene"], d_mis, near="exclude")
        col.mass_pass({"gene": n_units("gene", a)})
        M, dm = col.result("AM_gate", "gene")
        res = eng.run(M, syn_mass, syn_p, dm)
        oe = summ((res["O"].sum(1) / res["E"].sum(1))[:, None])
        out.append(f"[{subset}] missense O {res['O'][0].sum():.0f}, E {res['E'][0].sum():.1f}, O/E {oe['est'][0]:.3f} "
                   f"[{oe['lo95'][0]:.3f}, {oe['hi95'][0]:.3f}] (one-sided lo {oe['lo90'][0]:.3f})")
        log(out[-1])
        # models (descriptive): AlphaMissense, v3, v4
        run_scores(a, ["AlphaMissense", "Aim4v3_GBDT", "Aim4v4_GBDT"], A1_KINDS, ["all_genes"], tag + "m",
                   restrict="aim1_common")
        rp = P(a.results, f"{tag}m_top2reps.npz")
        if os.path.exists(rp):
            z = dict(np.load(rp))
            for m_name, ref in (("Aim4v4_GBDT", "AlphaMissense"), ("Aim4v3_GBDT", "AlphaMissense"),
                                ("Aim4v4_GBDT", "Aim4v3_GBDT")):
                r = _paired(z, m_name, ref)
                if r:
                    out.append(f"[{subset}] within-gene top-2% LR {m_name} / {ref} = {r[0]:.2f} [{r[1]:.2f}, {r[2]:.2f}]")
            t = pd.read_csv(P(a.results, f"{tag}m_top2pct.tsv"), sep="\t")
            for r in t[t["unit"].isin(["gene", "win50", "win50_rmc"])].itertuples():
                out.append(f"[{subset}] {r.score} {r.unit}: marginal {r.LRm_top2:.2f}; within {r.LRw_top2:.2f} "
                           f"[{r.lo95_w:.2f}, {r.hi95_w:.2f}]")
        open(P(a.results, "p3_summary_models.txt"), "w").write("\n".join(out) + "\n")
        log("A9 models:\n" + "\n".join(out))
        if tag == "p3":
            # main-finding readouts (A2 specification, D2 primary) for AlphaMissense + 36 tools
            run_scores(a, ["AlphaMissense"] + list(TOOLS), A1_KINDS, ["all_genes"], "p3", restrict="aim1_common")
            rep_report(a, float(oe["lo90"][0]), tag="p3", out_name="p3_summary.txt")
    a.results, a.cohort, a.rep_file, a.rep_cohorts, a.dedup_zhou = base, None, None, None, False


# ======================================================================== amendment A10.1: joint-LR overlap test
GNOMAD_GENE = P(EXT, "gnomad_constraint", "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
A10_TOOLS = {"AlphaMissense": ("am", [0.792, 0.906, 0.972, 0.990]),
             "REVEL": ("raw_REVEL_score", [0.644, 0.773, 0.879, 0.932]),
             "VEST4": ("raw_VEST4_score", [0.764, 0.861, 0.909, 0.965])}


def overlap_cohort(a, label, excl_tx=None, tools=None, fnames=None, lnames=None, site_keep=None, extra_cols=()):
    """R(T, F) = LR(PP3 | F) / LR(PP3) from marginal de novo LRs of the four PP3 x F categories (plan section 3k).
    Streams chromosomes: only the opportunity per category and the category of each de novo site are kept."""
    d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a)
    gz = pd.read_csv(GNOMAD_GENE, sep="\t", compression="gzip", usecols=["gene", "mis_z"])
    misz = gz.groupby("gene")["mis_z"].max()
    pp2_tx = tinfo["gene"].reindex(np.arange(18398)).map(misz).fillna(-np.inf).to_numpy() >= 3.09
    TL = tools if tools is not None else A10_TOOLS
    cols = list(dict.fromkeys(["key", "tx", PRIMARY_MODEL, "rmc_oe", "ccr_pct"] + [c for c, _ in TL.values()] + list(extra_cols)))
    fnames = fnames or ["PM1-like", "PM1-strict", "PP2-like", "PM1-like or PP2-like"]
    lnames = lnames or ["PP3>=Supporting", "PP3 Supporting only", "PP3 Moderate only", "PP3 >=Strong(3-4)"]
    combos = [(t, f, l) for t in TL for f in fnames for l in lnames]
    mass = {k: np.zeros(4) for k in combos}
    dbin = {k: np.full(len(d_mis), -1, np.int64) for k in combos}
    dkey = d_mis["key"].to_numpy()
    dchr = d_mis["chrom"].astype(str).to_numpy()
    for c, u in site_chunks(a, cols=cols):
        tx = u["tx"].to_numpy(np.int64)
        keep_ = ~np.isin(tx, excl_tx) if excl_tx is not None else np.ones(len(tx), bool)
        if site_keep is not None:
            keep_ &= site_keep(u)
        oe, ccr = u["rmc_oe"].to_numpy(float), u["ccr_pct"].to_numpy(float)
        fl = {"PM1-like": (oe < 0.4) | (ccr >= 95), "PM1-strict": oe < 0.2, "PP2-like": pp2_tx[tx]}
        fl["PM1-like or PP2-like"] = fl["PM1-like"] | fl["PP2-like"]
        mr = u[PRIMARY_MODEL].to_numpy(float)
        di = np.flatnonzero(dchr == str(c))
        pos = pd.Series(np.arange(len(u)), index=u["key"].to_numpy()).reindex(dkey[di]).to_numpy()
        okp = np.isfinite(pos)
        di, pos = di[okp], pos[okp].astype(np.int64)
        for tool, (col, cuts) in TL.items():
            s = u[col].to_numpy(float)
            sc = np.isfinite(s) & keep_
            lv = {"PP3>=Supporting": s >= cuts[0], "PP3 Supporting only": (s >= cuts[0]) & (s < cuts[1]),
                  "PP3 Moderate only": (s >= cuts[1]) & (s < cuts[2]), "PP3 >=Strong(3-4)": s >= cuts[2]}
            for f in fnames:
                for l in lnames:
                    cat = np.where(sc, 2 * lv[l].astype(np.int64) + fl[f].astype(np.int64), -1)
                    mass[(tool, f, l)] += np.bincount(cat[sc], weights=mr[sc], minlength=4)
                    dbin[(tool, f, l)][di] = cat[pos]
        del u
        gc.collect()
    rows, txt = [], []
    for (tool, fname, lname) in combos:
        dd = d_mis[["p"]].copy()
        dd["unit"] = 0
        dd["bin"] = dbin[(tool, fname, lname)]
        res = eng.run(sparse.csr_matrix(mass[(tool, fname, lname)].reshape(1, 4)), syn_mass, syn_p, dd)
        lr = Engine.lr(res, [np.array([3]), np.array([1, 3]), np.array([2, 3])])
        with np.errstate(all="ignore"):
            joint, lrF, lrP = (lr["LRm"][:, i] for i in range(3))
            R = joint / (lrF * lrP)
            pp3_given_F = joint / lrF
            Xn = res["O"] - res["E"]
            pp3_given_notF = (Xn[:, 2] / (Xn[:, 0] + Xn[:, 2])) / (res["E"][:, 2] / (res["E"][:, 0] + res["E"][:, 2]))
        b = R[1:][np.isfinite(R[1:])]
        lo, hi = (np.percentile(b, [2.5, 97.5]) if len(b) > 20 else (np.nan, np.nan))
        pts = np.log(pp3_given_F[0]) / np.log(PP3[0][1]) if pp3_given_F[0] > 0 else np.nan
        rows.append({"cohort": label, "tool": tool, "flag": fname, "pp3_level": lname,
                     "n_dnv_PP3_and_F": float(res["O"][0, 3]), "LR_PP3": float(lrP[0]), "LR_F": float(lrF[0]),
                     "LR_joint": float(joint[0]), "LR_PP3_given_F": float(pp3_given_F[0]),
                     "LR_PP3_given_notF": float(pp3_given_notF[0]), "R": float(R[0]), "R_lo95": float(lo),
                     "R_hi95": float(hi), "PP3_points_given_F": float(pts)})
        if lname == "PP3>=Supporting":
            prim = " [A10 PRIMARY]" if (label == "DD" and tool == "AlphaMissense" and fname == "PM1-like") else ""
            txt.append(f"{label:7s} {tool:13s} F={fname:20s} LR(PP3) {lrP[0]:.2f} LR(F) {lrF[0]:.2f} joint {joint[0]:.2f} "
                       f"-> R {R[0]:.2f} [{lo:.2f}, {hi:.2f}]; LR(PP3|F) {pp3_given_F[0]:.2f} (points {pts:.2f}) vs "
                       f"LR(PP3|not F) {pp3_given_notF[0]:.2f}; de novos PP3&F {res['O'][0, 3]:.0f}{prim}")
    return rows, txt


def cmd_a10_post(a):
    """Plan section 3k: A10.2 ClinVar sensitivity (PM1 +1/+2/+4; o/e < 0.4 or < 0.2; all genes / DDG2P monoallelic) and
    A10.3 model equivalence (two one-sided tests, margin [0.80, 1.25], 90% interval) from saved paired replicates."""
    base = a.results
    out = []
    m = pd.read_csv(P(base, "a4_clinvar_matched.tsv.gz"), sep="\t", dtype={"chrom": str})
    s = pd.read_csv(P(base, "supp_data1_clinvar_pp3_pm1_points.tsv.gz"), sep="\t", dtype={"chrom": str})
    x = s.merge(m[["clinvar_id", "chrom", "tx", "rmc_oe", "ccr_pct"]], on=["clinvar_id", "chrom"], how="left")
    tinfo = load_tinfo(a)
    dd = tinfo["ddg2p"].fillna(False).astype(bool)
    x["ddg2p"] = x["tx"].map(dd).fillna(False).astype(bool)
    nest = pd.read_csv(P(base, "a1_nested.tsv"), sep="\t")
    r_e2 = nest[nest["stratum"] == "all_genes"].set_index("score")["r_E2"]
    pm1_defs = {"oe<0.4|CCR>=95": (x["rmc_oe"] < 0.4) | (x["ccr_pct"] >= 95), "oe<0.2": x["rmc_oe"] < 0.2}
    rows = []
    for gset, gm in (("all genes", np.ones(len(x), bool)), ("DDG2P monoallelic", x["ddg2p"].to_numpy())):
        for dname, pmask in pm1_defs.items():
            for pm1_pts in (1, 2, 4):
                for t in ("AlphaMissense", "REVEL", "VEST4"):
                    pub = x[f"{t}_PP3_published"].to_numpy(float)
                    pm = pmask.fillna(False).to_numpy(bool)
                    pm1 = np.where(pm, float(pm1_pts), 0.0)
                    unc = pub + pm1
                    cap = np.where(pm, np.minimum(unc, 4.0), unc)
                    cond = pub * float(r_e2[t]) + pm1
                    for cls in ("P/LP", "VUS"):
                        k = gm & (x["clinvar_class"] == cls).to_numpy() & np.isfinite(pub)
                        kp = k & pm
                        pp3 = k & (pub >= 1)
                        rows.append({"genes": gset, "pm1_def": dname, "pm1_points": pm1_pts, "tool": t, "cls": cls,
                                     "n": int(k.sum()), "share_pp3_at_pm1": float((pp3 & pm).sum() / max(pp3.sum(), 1)),
                                     "ge6_uncapped": int((k & (unc >= 6)).sum()), "ge6_capped": int((k & (cap >= 6)).sum()),
                                     "ge6_conditional": int((k & (cond >= 6)).sum()),
                                     "pm1_mean_capped": float(cap[kp].mean()) if kp.any() else np.nan,
                                     "pm1_mean_conditional": float(cond[kp].mean()) if kp.any() else np.nan})
    cs = pd.DataFrame(rows)
    cs.to_csv(P(base, "a10_clinvar_sensitivity.tsv"), sep="\t", index=False, float_format="%.4g")
    for r in cs[(cs.tool == "AlphaMissense") & (cs.cls == "VUS")].itertuples():
        out.append(f"A10.2 {r.genes:17s} PM1 {r.pm1_def:15s} +{r.pm1_points}: AlphaMissense VUS reaching >=6 from PP3+PM1: "
                   f"uncapped {r.ge6_uncapped:,d}, capped {r.ge6_capped:,d}, conditional {r.ge6_conditional:,d}; "
                   f"PP3-positive at PM1 sites {r.share_pp3_at_pm1:.1%}")
    # A10.3 equivalence
    sets = [("DDD (all genes; cold-gene replicates were not saved)", "m4v3nested_ddd_top2reps.npz", "", ["Aim4v3_GBDT"]),
            ("autism", "m4v4rep_top2reps.npz", "replication_affected", ["Aim4v3_GBDT", "Aim4v4_GBDT"]),
            ("public3", "p3m_top2reps.npz", "public3", ["Aim4v3_GBDT", "Aim4v4_GBDT"])]
    for name, f, sub, models in sets:
        pth = P(base, sub, f)
        if not os.path.exists(pth):
            out.append(f"A10.3 {name}: replicates not found ({f})")
            continue
        z = dict(np.load(pth))
        for mdl in models:
            xx, yy = z.get(f"{mdl}|gene|all_genes|LRw"), z.get(f"AlphaMissense|gene|all_genes|LRw")
            if xx is None or yy is None:
                continue
            with np.errstate(all="ignore"):
                dl = np.log(xx) - np.log(yy)
            b = dl[1:][np.isfinite(dl[1:])]
            lo, hi = np.exp(np.percentile(b, [5, 95]))
            eq = (lo > 0.80) and (hi < 1.25)
            out.append(f"A10.3 {name}: {mdl} / AlphaMissense within-gene top-2% = {np.exp(dl[0]):.2f}, 90% CI [{lo:.2f}, "
                       f"{hi:.2f}] -> {'EQUIVALENT within [0.80, 1.25]' if eq else 'equivalence not shown'}")
    open(P(base, "a10_post_summary.txt"), "w").write("\n".join(out) + "\n")
    log("A10.2/A10.3:\n" + "\n".join(out))


# ======================================================================== amendment A11: robustness (descriptive)
def _a11_tools(a):
    n = pd.read_csv(P(a.results, "a1_nested.tsv"), sep="\t")
    n = n[n["stratum"] == "all_genes"]
    return [s for s in n.loc[n["A"] >= PP3[0][1], "score"] if s != "AlphaMissense"]


def _a11_summary(a, tag, strata, kinds):
    t = pd.read_csv(P(a.results, f"{tag}_top2pct.tsv"), sep="\t")
    out = []
    for st in strata:
        g = t[t["stratum"] == st]
        rows = []
        for s, gs in g.groupby("score", sort=False):
            gs = gs.set_index("unit")
            if "gene" not in gs.index:
                continue
            A = float(gs.loc["gene", "LRm_top2"])
            r = {"score": s, "A": A}
            for k in kinds:
                v = float(gs.loc[k, "LRw_top2"]) if k in gs.index else np.nan
                r[k] = v
                r["r_" + k] = np.log(v) / np.log(A) if (A > 1 and v > 0) else np.nan
            rows.append(r)
        df = pd.DataFrame(rows)
        inf = df[(df["score"] != "AlphaMissense") & (df["A"] >= PP3[0][1])]
        am = df[df["score"] == "AlphaMissense"]
        share = (inf["r_win50"] < 0.5).mean() if "r_win50" in inf else np.nan
        out.append(f"{tag} [{st}] informative {len(inf)}; median r " +
                   ", ".join(f"{k} {inf['r_' + k].median():.3f}" for k in kinds) +
                   (f"; share r(win50) < 0.5 = {share:.2f}" if "win50" in kinds else "") +
                   ("; AlphaMissense A {:.2f} ".format(float(am["A"].iloc[0])) +
                    ", ".join(f"{k} {float(am[k].iloc[0]):.2f}" for k in kinds) if len(am) else ""))
    return out


def _mr_scaled_site_chunks(orig, factor):
    """A11.1: CpG-transition (mutclass 0) rates multiplied by `factor` at every site, missense and synonymous."""
    def sc(a_, chroms=AUTOSOMES, near="exclude", cls=1, extra=(), cols=None):
        cc = None if cols is None else list(dict.fromkeys(list(cols) + ["mutclass"]))
        for c, u in orig(a_, chroms=chroms, near=near, cls=cls, extra=extra, cols=cc):
            u = u.copy()
            u[PRIMARY_MODEL] = np.where(u["mutclass"].to_numpy() == 0, u[PRIMARY_MODEL].to_numpy() * factor,
                                        u[PRIMARY_MODEL].to_numpy())
            yield c, u
    return sc


def _clinvar_excluded_site_chunks(orig, keys_by_chrom):
    """A11.2: sites of ClinVar P/LP missense variants removed from opportunity (and hence from de novo mapping)."""
    def sc(a_, chroms=AUTOSOMES, near="exclude", cls=1, extra=(), cols=None):
        for c, u in orig(a_, chroms=chroms, near=near, cls=cls, extra=extra, cols=cols):
            ks = keys_by_chrom.get(str(c))
            if cls == 1 and ks is not None and "key" in u:
                u = u[~u["key"].isin(ks)]
            yield c, u
    return sc


def cmd_a11(a):
    """Plan section 3l, run in parts (--part): strata, sweep, cpg_down, cpg_up, noclinvar."""
    global site_chunks
    tools = ["AlphaMissense"] + _a11_tools(a)
    pos_setup(a)
    out = []
    if a.part == "strata":
        strata = ["shet_high", "shet_low", "ddg2p_LoF", "ddg2p_nonLoF", "rmc_multi", "rmc_single"]
        kinds = ["gene", "rmc", "win50", "win50_rmc"]
        run_scores(a, tools, kinds, strata, "a11_strata", restrict="aim1_common")
        out = _a11_summary(a, "a11_strata", strata, kinds)
    elif a.part == "sweep":
        kinds = ["gene", "win25", "rand25", "win50", "rand50", "win100", "rand100"]
        run_scores(a, tools, kinds, ["all_genes"], "a11_sweep", restrict="aim1_common")
        out = _a11_summary(a, "a11_sweep", ["all_genes"], kinds)
    elif a.part in ("cpg_down", "cpg_up", "noclinvar"):
        orig = site_chunks
        try:
            if a.part == "noclinvar":
                full = []                                    # A11.2: P/LP of any review status
                with gzip.open(a.clinvar, "rt") as fh:
                    for ln in fh:
                        if ln.startswith("#"):
                            continue
                        f = ln.split("\t", 8)
                        if f[0] not in AUTOSOMES or len(f[3]) != 1 or len(f[4]) != 1 or f[4] not in B2I:
                            continue
                        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
                        if "missense_variant" in info.get("MC", "") and A4_CLASS.get(info.get("CLNSIG", "")) == "P/LP":
                            full.append((f[0], int(f[1]) * 4 + B2I[f[4]]))
                kb = pd.DataFrame(full, columns=["chrom", "key"]).groupby("chrom")["key"].apply(lambda s: set(s.tolist()))
                site_chunks = _clinvar_excluded_site_chunks(orig, {k: np.array(sorted(v)) for k, v in kb.items()})
                log(f"A11.2: {len(full):,} ClinVar P/LP missense sites (any review status) excluded")
            else:
                site_chunks = _mr_scaled_site_chunks(orig, 0.919 if a.part == "cpg_down" else 1.088)
            kinds = ["gene", "rmc", "win50", "win50_rmc"]
            tag = f"a11_{a.part}"
            run_scores(a, tools, kinds, ["all_genes"], tag, restrict="aim1_common")
            out = _a11_summary(a, tag, ["all_genes"], kinds)
        finally:
            site_chunks = orig
    open(P(a.results, f"a11_{a.part}_summary.txt"), "w").write("\n".join(out) + "\n")
    log(f"A11 {a.part}:\n" + "\n".join(out))


def cmd_a11sim(a):
    """Plan section 3l, A11.5: planted within-gene signal at the top 2% of AlphaMissense (lambda in 1, 2, 4), real gene
    excess and real sizes; truth = the estimator applied to expected counts; bias and 95% coverage over replicates."""
    rng = np.random.default_rng(SEED)
    d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a)
    pos_setup(a)
    cols = ["tx", PRIMARY_MODEL, "am", "prot_pos"]
    parts = {k: [] for k in cols}
    for c, u in site_chunks(a, cols=cols):
        for k in cols:
            parts[k].append(u[k].to_numpy())
    tx = np.concatenate(parts["tx"]).astype(np.int64)
    mr = np.concatenate(parts[PRIMARY_MODEL]).astype(np.float64)
    am = np.concatenate(parts["am"]).astype(np.float32)
    pp = np.concatenate(parts["prot_pos"]).astype(np.int64)
    del parts
    gc.collect()
    ok = np.isfinite(am)
    tx, mr, am, pp = tx[ok], mr[ok], am[ok], pp[ok]
    o = np.argsort(-am, kind="stable")
    top = np.zeros(len(am), bool)
    top[o[:np.searchsorted(np.cumsum(mr[o]) / mr.sum(), 0.02) + 1]] = True     # top 2% of opportunity
    units = {"gene": tx, "win50": unit_of(pd.DataFrame({"tx": tx, "prot_pos": pp}), "win50"),
             "rand50": unit_of(pd.DataFrame({"tx": tx, "prot_pos": pp}), "rand50")}
    nb = 2
    mass = {k: sparse.csr_matrix((mr, (u_, top.astype(np.int64))), shape=(int(u_.max()) + 1, nb)) for k, u_ in units.items()}
    c = len(syn_p) / syn_mass
    n_mis = len(d_mis)
    E_tx = np.bincount(tx, weights=mr * c, minlength=tx.max() + 1)
    O_tx = np.bincount(d_mis["tx"].to_numpy(np.int64), minlength=tx.max() + 1)[:len(E_tx)]
    g = np.maximum(O_tx - E_tx, 0.0)                         # real per-gene excess
    g[np.bincount(tx, minlength=len(g)) == 0] = 0.0          # genes without scored sites cannot receive planted de novos
    n_path = int(round(n_mis - c * mr.sum()))
    order = np.argsort(tx, kind="stable")
    starts = np.searchsorted(tx[order], np.arange(tx.max() + 2))
    n_prob = eng.boot.n
    rows = []
    for lam in (1.0, 2.0, 4.0):
        w = mr * np.where(top, lam, 1.0)
        cw = np.cumsum(w[order])
        cw0 = np.r_[0.0, cw]
        # expected pathogenic count per site: P_t * w_i / sum_{j in t} w_j with P_t = n_path * g_t / sum g
        W_t = np.bincount(tx, weights=w, minlength=tx.max() + 1)
        P_t = n_path * g / g.sum()
        exp_path = P_t[tx] * w / np.maximum(W_t[tx], 1e-300)
        truth = {}
        for k, u_ in units.items():
            Xub = sparse.csr_matrix((exp_path, (u_, top.astype(np.int64))), shape=mass[k].shape).toarray()
            Eub = mass[k].toarray() * c
            Xu, Eu = Xub.sum(1), Eub.sum(1)
            with np.errstate(all="ignore"):
                den = np.nansum(np.where(Eu > 0, Xu * Eub[:, 1] / Eu, 0.0))
            truth[k] = Xub[:, 1].sum() / den
        cbg = np.cumsum(mr)
        for rep in range(int(a.reps)):
            n_p = rng.poisson(n_path)
            n_b = rng.poisson(n_mis - n_path)
            bg = np.searchsorted(cbg, rng.random(n_b) * cbg[-1])
            gt = rng.choice(len(g), size=n_p, p=g / g.sum())
            lo_, hi_ = cw0[starts[gt]], cw0[starts[gt + 1]]
            pi = order[np.minimum(np.searchsorted(cw, lo_ + rng.random(n_p) * (hi_ - lo_)), len(order) - 1)]
            sites = np.r_[bg, pi]
            dn = pd.DataFrame({"p": rng.integers(0, n_prob, len(sites)), "bin": top[sites].astype(np.int64)})
            sp = rng.integers(0, n_prob, rng.poisson(c * syn_mass))
            for k, u_ in units.items():
                dk = dn.assign(unit=u_[sites])
                res = eng.run(mass[k], syn_mass, sp, dk)
                lr = Engine.lr(res, [np.array([1])])["LRw"][:, 0]
                est = float(lr[0])
                lo, hi = np.percentile(lr[1:][np.isfinite(lr[1:])], [2.5, 97.5])
                rows.append({"lambda": lam, "unit": k, "rep": rep, "truth": truth[k], "est": est, "lo95": lo, "hi95": hi,
                             "covered": bool(lo <= truth[k] <= hi)})
            if rep % 20 == 0:
                log(f"A11.5 lambda {lam}: replicate {rep}")
        pd.DataFrame(rows).to_csv(P(a.results, "a11_simulation.tsv"), sep="\t", index=False, float_format="%.4g")
    sim = pd.DataFrame(rows)
    out = []
    for (lam, k), gg in sim.groupby(["lambda", "unit"]):
        out.append(f"lambda {lam:.0f} {k:6s}: truth {gg['truth'].iloc[0]:.3f}; mean estimate {gg['est'].mean():.3f} "
                   f"(median {gg['est'].median():.3f}); relative bias {gg['est'].mean() / gg['truth'].iloc[0] - 1:+.3f}; "
                   f"95% coverage {gg['covered'].mean():.2f} over {len(gg)} replicates")
    open(P(a.results, "a11_simulation_summary.txt"), "w").write("\n".join(out) + "\n")
    log("A11.5 simulation:\n" + "\n".join(out))


# ======================================================================== amendment A13: overlap within disease genes
def cmd_a13(a):
    """Plan section 3n: R within DDG2P monoallelic genes (marginal), gene-stratified R (all genes) and gene-stratified R
    within DDG2P genes, for AlphaMissense/REVEL/VEST4 PP3 (Supporting+ and Strong tier) x PM1-like / PM1-strict."""
    base = a.results
    rows, out = [], []
    for label, kw in (("DD", {}), ("autism", {"cohort": "Affected"})):
        a.cohort = kw.get("cohort")
        a.results = P(base, "a13_tmp")
        os.makedirs(a.results, exist_ok=True)
        d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a)
        ntx = 18398
        ddmask = tinfo["ddg2p"].reindex(np.arange(ntx)).fillna(False).to_numpy(bool)
        cols = ["key", "tx", PRIMARY_MODEL, "rmc_oe", "ccr_pct"] + [cc for cc, _ in A10_TOOLS.values()]
        combos = [(t, f, l) for t in A10_TOOLS for f in ("PM1-like", "PM1-strict") for l in ("PP3>=Supporting", "PP3 >=Strong(3-4)")]
        mass = {k: np.zeros(ntx * 4) for k in combos}
        dbin = {k: np.full(len(d_mis), -1, np.int64) for k in combos}
        dk = d_mis["key"].to_numpy(); dch = d_mis["chrom"].astype(str).to_numpy()
        for ch, u in site_chunks(a, cols=cols):
            tx = u["tx"].to_numpy(np.int64)
            oe, ccr = u["rmc_oe"].to_numpy(float), u["ccr_pct"].to_numpy(float)
            fl = {"PM1-like": (oe < 0.4) | (ccr >= 95), "PM1-strict": oe < 0.2}
            mr = u[PRIMARY_MODEL].to_numpy(float)
            di = np.flatnonzero(dch == str(ch))
            pos = pd.Series(np.arange(len(u)), index=u["key"].to_numpy()).reindex(dk[di]).to_numpy()
            okp = np.isfinite(pos)
            di, pos = di[okp], pos[okp].astype(np.int64)
            for t, (col, cuts) in A10_TOOLS.items():
                s_ = u[col].to_numpy(float)
                sc = np.isfinite(s_)
                lv = {"PP3>=Supporting": s_ >= cuts[0], "PP3 >=Strong(3-4)": s_ >= cuts[2]}
                for f in ("PM1-like", "PM1-strict"):
                    for l in lv:
                        cat = np.where(sc, 2 * lv[l].astype(np.int64) + fl[f].astype(np.int64), -1)
                        mass[(t, f, l)] += np.bincount(tx[sc] * 4 + cat[sc], weights=mr[sc], minlength=ntx * 4)
                        dbin[(t, f, l)][di] = cat[pos]
            del u
            gc.collect()
        utx = d_mis["tx"].to_numpy(np.int64)
        for (t, f, l) in combos:
            M = sparse.csr_matrix(mass[(t, f, l)].reshape(ntx, 4))
            dd = d_mis[["p"]].copy(); dd["unit"] = utx; dd["bin"] = dbin[(t, f, l)]
            g3 = [np.array([3]), np.array([1, 3]), np.array([2, 3])]
            res_all = eng.run(M, syn_mass, syn_p, dd)
            res_dd = eng.run(M, syn_mass, syn_p, dd, unit_mask=ddmask)
            for setting, res, kind in (("DDG2P, marginal", res_dd, "LRm"), ("all genes, within gene", res_all, "LRw"),
                                       ("DDG2P, within gene", res_dd, "LRw")):
                lr = Engine.lr(res, g3)[kind]
                with np.errstate(all="ignore"):
                    R = lr[:, 0] / (lr[:, 1] * lr[:, 2])
                    given = lr[:, 0] / lr[:, 1]
                b = R[1:][np.isfinite(R[1:])]
                lo, hi = (np.percentile(b, [2.5, 97.5]) if len(b) > 20 else (np.nan, np.nan))
                rows.append({"cohort": label, "tool": t, "flag": f, "pp3_level": l, "setting": setting, "R": float(R[0]),
                             "R_lo95": float(lo), "R_hi95": float(hi), "LR_PP3": float(lr[0, 2]), "LR_F": float(lr[0, 1]),
                             "LR_joint": float(lr[0, 0]), "LR_PP3_given_F": float(given[0])})
                prim = " [A13 key]" if (setting == "DDG2P, within gene" and t == "AlphaMissense" and f == "PM1-like"
                                        and l == "PP3>=Supporting") else ""
                out.append(f"{label:6s} {t:13s} {f:10s} {l:18s} {setting:22s}: R {R[0]:.2f} [{lo:.2f}, {hi:.2f}]; "
                           f"LR(PP3|F) {given[0]:.2f}; LR(PP3) {lr[0, 2]:.2f}{prim}")
        a.results = base
        pd.DataFrame(rows).to_csv(P(base, "a13_within_genes.tsv"), sep="\t", index=False, float_format="%.4g")
        open(P(base, "a13_summary.txt"), "w").write("\n".join(out) + "\n")
        log(f"A13 {label}:\n" + "\n".join(o for o in out if o.startswith(label)))
    a.cohort = None


# ======================================================================== amendment A12 (descriptive)
def cmd_a12(a):
    """Plan section 3m. part 'decomp' (A12.1): excess-side and null-side P(PP3 | PM1-like) for AlphaMissense, REVEL, VEST4
    in DD and autism. part 'classify' (A12.4): PS2 + PM2 + PM1-like + AlphaMissense PP3 under uncapped / capped /
    conditional rules for DD de novo missense; class counts and de novo-excess fraction pathogenic per class."""
    base = a.results
    if a.part == "lto":                                      # A12.3: leave-top-genes-out (DD)
        global _EXCL_TX
        d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a)
        c = len(syn_p) / syn_mass
        E = np.zeros(18398)
        for ch, u in site_chunks(a, cols=["key", "tx", PRIMARY_MODEL]):
            E += np.bincount(u["tx"].to_numpy(np.int64), weights=u[PRIMARY_MODEL].to_numpy(float) * c, minlength=18398)[:18398]
        O = np.bincount(d_mis["tx"].to_numpy(np.int64), minlength=18398)[:18398]
        top = np.argsort(-(O - E))[:50]
        _EXCL_TX = top
        genes = tinfo["gene"].reindex(top).tolist()
        out = [f"A12.3: excluded 50 transcripts with the largest DD de novo missense excess (O - E): {', '.join(map(str, genes))}",
               f"excluded de novo missense {int(O[top].sum())} of {int(O.sum())}; expected {E[top].sum():.0f}"]
        rows, txt = overlap_cohort(a, "DD_noTop50", excl_tx=top)
        pd.DataFrame(rows).to_csv(P(base, "a12_lto_overlap.tsv"), sep="\t", index=False, float_format="%.4g")
        out += [t for t in txt if "PM1-like " in t or "PM1-strict" in t]
        tools = ["AlphaMissense"] + _a11_tools(a)
        pos_setup(a)
        kinds = ["gene", "rmc", "win50", "win50_rmc", "rand50"]
        run_scores(a, tools, kinds, ["no_top50"], "a12_lto", restrict="aim1_common")
        out += _a11_summary(a, "a12_lto", ["no_top50"], kinds)
        open(P(base, "a12_lto_summary.txt"), "w").write("\n".join(out) + "\n")
        log("\n".join(out))
        return
    if a.part == "plddt":                                    # task 10: overlap by AlphaFold confidence (ordered vs disordered)
        out, allrows = [], []
        for label, kw in (("DD", {}), ("autism", {"cohort": "Affected"})):
            a.cohort = kw.get("cohort")
            for sname, fn in (("pLDDT>=70", lambda u: u["af_plddt"].to_numpy(float) >= 70),
                              ("pLDDT<70", lambda u: u["af_plddt"].to_numpy(float) < 70)):
                rows, txt = overlap_cohort(a, f"{label}|{sname}", fnames=["PM1-like", "PP2-like"],
                                           lnames=["PP3>=Supporting", "PP3 >=Strong(3-4)"], site_keep=fn, extra_cols=("af_plddt",))
                allrows += rows
                out += [t for t in txt if "AlphaMissense" in t or "REVEL" in t or "VEST4" in t]
        a.cohort = None
        pd.DataFrame(allrows).to_csv(P(base, "a12_plddt.tsv"), sep="\t", index=False, float_format="%.4g")
        open(P(base, "a12_plddt_summary.txt"), "w").write("\n".join(out) + "\n")
        log("\n".join(out))
        return
    if a.part == "allr":                                     # A12.2: R for all predictors, top-2% flags
        pop_informed = {"AlphaMissense", "PrimateAI", "MPC", "MVP", "M-CAP", "BayesDel_addAF", "ClinPred", "MetaRNN", "MetaSVM",
                        "MetaLR", "VARITY_R", "VARITY_ER", "VARITY_R_LOO", "VARITY_ER_LOO", "MisFit_D", "MisFit_S", "popEVE"}
        # population constraint or allele frequency among inputs or training labels (training ledger, dbNSFP readme)
        names = ["AlphaMissense"] + list(TOOLS) + list(TOOLS54)
        colof = {t: ("am" if t == "AlphaMissense" else score_col(t)) for t in names}
        out, allrows = [], []
        for label, kw in (("DD", {}), ("autism", {"cohort": "Affected"})):
            a.cohort = kw.get("cohort")
            a.results = P(base, "a12_tmp")
            os.makedirs(a.results, exist_ok=True)
            if label == "DD":                                # thresholds: top 2% of mutational opportunity, all sites
                H = {t: np.zeros(10001) for t in names}
                for ch, u in site_chunks(a, cols=["key", PRIMARY_MODEL] + [colof[t] for t in names]):
                    mr = u[PRIMARY_MODEL].to_numpy(float)
                    for t in names:
                        v = u[colof[t]].to_numpy(float); f = np.isfinite(v)
                        H[t] += np.bincount(np.clip((v[f] * 10000).astype(np.int64), 0, 10000), weights=mr[f], minlength=10001)
                thr = {}
                for t in names:
                    cum = np.cumsum(H[t]) / H[t].sum()
                    thr[t] = (np.searchsorted(cum, 0.98) + 1) / 10000.0
                tl = {t: (colof[t], (thr[t], np.inf, np.inf)) for t in names}
            rows, _ = overlap_cohort(a, label, tools=tl, fnames=["PM1-like", "PP2-like"], lnames=["PP3>=Supporting"])
            allrows += rows
            a.results = base
        df = pd.DataFrame(allrows)
        df["pop_informed"] = df["tool"].isin(pop_informed)
        df["threshold_top2"] = df["tool"].map(lambda t: thr[t])
        df.to_csv(P(base, "a12_allr.tsv"), sep="\t", index=False, float_format="%.4g")
        from scipy.stats import mannwhitneyu
        for (coh, fl), g in df[df["LR_PP3"] >= PP3[0][1]].groupby(["cohort", "flag"]):
            a_, b_ = g.loc[g["pop_informed"], "R"], g.loc[~g["pop_informed"], "R"]
            pv = mannwhitneyu(a_, b_).pvalue if len(a_) > 2 and len(b_) > 2 else np.nan
            out.append(f"A12.2 {coh:6s} F={fl:9s} informative {len(g)}: median R {g['R'].median():.2f} "
                       f"(range {g['R'].min():.2f}-{g['R'].max():.2f}); share R<1 upper bound {(g['R_hi95'] < 1).mean():.2f}; "
                       f"population-informed {a_.median():.2f} (n {len(a_)}) vs other {b_.median():.2f} (n {len(b_)}), "
                       f"Mann-Whitney p {pv:.3g}")
        open(P(base, "a12_allr_summary.txt"), "w").write("\n".join(out) + "\n")
        log("\n".join(out))
        return
    if a.part == "decomp":
        rows = []
        for label, kw in (("DD", {}), ("autism", {"cohort": "Affected"})):
            a.cohort = kw.get("cohort")
            a.results = P(base, "a12_tmp")
            os.makedirs(a.results, exist_ok=True)
            d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a)
            c = len(syn_p) / syn_mass
            cols = ["key", PRIMARY_MODEL, "rmc_oe", "ccr_pct"] + [cc for cc, _ in A10_TOOLS.values()]
            E = {t: np.zeros(4) for t in A10_TOOLS}
            O = {t: np.zeros(4) for t in A10_TOOLS}
            dk = d_mis["key"].to_numpy(); dch = d_mis["chrom"].astype(str).to_numpy()
            for ch, u in site_chunks(a, cols=cols):
                pm = ((u["rmc_oe"].to_numpy(float) < 0.4) | (u["ccr_pct"].to_numpy(float) >= 95)).astype(np.int64)
                mr = u[PRIMARY_MODEL].to_numpy(float)
                di = np.flatnonzero(dch == str(ch))
                pos = pd.Series(np.arange(len(u)), index=u["key"].to_numpy()).reindex(dk[di]).to_numpy()
                pos = pos[np.isfinite(pos)].astype(np.int64)
                for t, (col, cuts) in A10_TOOLS.items():
                    s_ = u[col].to_numpy(float)
                    sc = np.isfinite(s_)
                    cat = np.where(sc, 2 * (s_ >= cuts[0]).astype(np.int64) + pm, -1)
                    E[t] += np.bincount(cat[sc], weights=mr[sc] * c, minlength=4)
                    cd = cat[pos]
                    O[t] += np.bincount(cd[cd >= 0], minlength=4)
            for t in A10_TOOLS:
                X = O[t] - E[t]
                rows.append({"cohort": label, "tool": t, "O": O[t].tolist(), "E": np.round(E[t], 1).tolist(),
                             "P_PP3_given_F_excess": X[3] / (X[1] + X[3]), "P_PP3_given_F_null": E[t][3] / (E[t][1] + E[t][3]),
                             "P_PP3_excess_all": (X[2] + X[3]) / X.sum(), "P_PP3_null_all": (E[t][2] + E[t][3]) / E[t].sum()})
            a.results = base
        df = pd.DataFrame(rows)
        df.to_csv(P(base, "a12_decomp.tsv"), sep="\t", index=False, float_format="%.4g")
        out = [f"{r.cohort:6s} {r.tool:13s} P(PP3|PM1-like): excess {r.P_PP3_given_F_excess:.3f} vs null {r.P_PP3_given_F_null:.3f} "
               f"(ratio = LR(PP3|F) {r.P_PP3_given_F_excess / r.P_PP3_given_F_null:.2f}); overall P(PP3): excess "
               f"{r.P_PP3_excess_all:.3f} vs null {r.P_PP3_null_all:.3f} (LR {r.P_PP3_excess_all / r.P_PP3_null_all:.2f})"
               for r in df.itertuples()]
    else:
        d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a)
        c = len(syn_p) / syn_mass
        ov = pd.read_csv(P(base, "a10_overlap.tsv"), sep="\t")
        ov = ov[(ov.cohort == "DD") & (ov.tool == "AlphaMissense") & (ov.flag == "PM1-like")].set_index("pp3_level")
        lr_tier = {1: ov.loc["PP3 Supporting only", "LR_PP3_given_F"], 2: ov.loc["PP3 Moderate only", "LR_PP3_given_F"],
                   3: ov.loc["PP3 >=Strong(3-4)", "LR_PP3_given_F"], 4: ov.loc["PP3 >=Strong(3-4)", "LR_PP3_given_F"]}
        cond_pts = {k: max(0.0, float(np.log(v) / np.log(PP3[0][1]))) for k, v in lr_tier.items()}
        cp_lookup = np.array([0.0] + [cond_pts[k] for k in (1, 2, 3, 4)])
        cuts = A10_TOOLS["AlphaMissense"][1]
        rules = ("uncapped", "capped", "conditional")
        E = {r: np.zeros(3) for r in rules}
        O = {r: np.zeros(3) for r in rules}
        dk = d_mis["key"].to_numpy(); dch = d_mis["chrom"].astype(str).to_numpy()
        n_change = {"uncapped->capped": 0.0, "uncapped->conditional": 0.0}
        for ch, u in site_chunks(a, cols=["key", PRIMARY_MODEL, "rmc_oe", "ccr_pct", "am"]):
            s_ = u["am"].to_numpy(float)
            sc = np.isfinite(s_)
            pp3 = np.select([s_ >= cuts[3], s_ >= cuts[2], s_ >= cuts[1], s_ >= cuts[0]], [4, 3, 2, 1], 0)
            pm = (u["rmc_oe"].to_numpy(float) < 0.4) | (u["ccr_pct"].to_numpy(float) >= 95)
            base_pts = 5.0                                                    # PS2 (4) + PM2 (1)
            pts = {"uncapped": base_pts + pp3 + 2 * pm,
                   "capped": base_pts + np.where(pm, np.minimum(pp3 + 2, 4), pp3),
                   "conditional": base_pts + 2 * pm + np.where(pm, cp_lookup[pp3], pp3)}
            mr = u[PRIMARY_MODEL].to_numpy(float)
            di = np.flatnonzero(dch == str(ch))
            pos = pd.Series(np.arange(len(u)), index=u["key"].to_numpy()).reindex(dk[di]).to_numpy()
            pos = pos[np.isfinite(pos)].astype(np.int64)
            pos = pos[sc[pos]]
            for r in rules:
                cls = np.select([pts[r] >= 10, pts[r] >= 6], [2, 1], 0)            # 2 = P, 1 = LP, 0 = VUS
                E[r] += np.bincount(cls[sc], weights=mr[sc] * c, minlength=3)
                O[r] += np.bincount(cls[pos], minlength=3)
                if r != "uncapped":
                    cu = np.select([pts["uncapped"] >= 10, pts["uncapped"] >= 6], [2, 1], 0)
                    n_change[f"uncapped->{r}"] += float((cls[pos] != cu[pos]).sum())
        out = [f"conditional PP3 points at PM1-like sites by tier: {cond_pts}"]
        for r in rules:
            for k, nm in ((2, "Pathogenic"), (1, "Likely pathogenic")):
                frac = 1 - E[r][k] / O[r][k] if O[r][k] > 0 else np.nan
                out.append(f"{r:11s} {nm:17s}: {O[r][k]:6.0f} de novos, expected background {E[r][k]:8.1f}, "
                           f"fraction pathogenic {frac:.3f}")
        out += [f"de novos changing class {k}: {int(v):,}" for k, v in n_change.items()]
    open(P(base, f"a12_{a.part}_summary.txt"), "w").write("\n".join(out) + "\n")
    log(f"A12 {a.part}:\n" + "\n".join(out))
    a.cohort = None


def cmd_overlap(a):
    base = a.results
    allrows, alltxt = [], []
    for label, setup_kw in (("DD", {}), ("autism", {"cohort": "Affected"}),
                            ("public3", {"cohort": "Affected", "rep_file": P(EXT, "public3", "public3_loader_format.tsv"),
                                         "rep_cohorts": P3_STUDIES, "dedup_zhou": True})):
        a.cohort, a.rep_file, a.rep_cohorts, a.dedup_zhou = (setup_kw.get("cohort"), setup_kw.get("rep_file"),
                                                             setup_kw.get("rep_cohorts"), setup_kw.get("dedup_zhou", False))
        a.results = P(base, "a10_tmp")
        os.makedirs(a.results, exist_ok=True)
        rows, txt = overlap_cohort(a, label)
        allrows += rows
        alltxt += txt
        a.results = base
        pd.DataFrame(allrows).to_csv(P(base, "a10_overlap.tsv"), sep="\t", index=False, float_format="%.4g")
        open(P(base, "a10_overlap_summary.txt"), "w").write("\n".join(alltxt) + "\n")
        log("A10.1 overlap:\n" + "\n".join(txt))
    a.cohort, a.rep_file, a.rep_cohorts, a.dedup_zhou = None, None, None, False


def cmd_m4_nested(a):
    """Plan section 3f: each Aim 4 model's LR through the A1 nested levels A-E on DDD de novos only."""
    a.centre_only = "DDD"
    kinds = ["gene", "rmc", "ccr", "pseudo", "win50", "pseudo_rmc", "win50_rmc"]
    run_scores(a, ["AlphaMissense"] + list(M4_COLS), kinds, ["all_genes"], "m4nested_ddd")
    rep_report(a, float("nan"), tag="m4nested_ddd", out_name="m4nested_ddd_summary.txt")


# ======================================================================== amendment A6
A1_KINDS = ["gene", "rmc", "ccr", "pseudo", "win50", "pseudo_rmc", "win50_rmc"]
A6_ORIG = [("aim1tools_top2pct.tsv", ""), ("rep_top2pct.tsv", "replication_affected"), ("a1v54_top2pct.tsv", "")]


def cmd_a6(a):
    """Plan section 3g, in order of importance: A6.3 models in autism (+ siblings), A6.2 sibling negative control,
    A6.1 reruns that save bootstrap replicates (point estimates must equal the committed tables)."""
    import shutil
    base = a.results
    pos_setup(a)
    m4 = os.path.exists(P(a.fs, "col", "_m4_done"))
    # A6.3: Aim 4 models, unchanged, in autism probands and siblings
    if m4:
        for pheno, sub, tag in (("Affected", "replication_affected", "m4rep"), ("Unaffected", "replication_unaffected", "m4sib")):
            a.results, a.cohort = P(base, sub), pheno
            os.makedirs(a.results, exist_ok=True)
            models = list(M4_COLS) + (list(M4V3) if os.path.exists(P(a.fs, "col", "_m4v3_done")) else [])
            run_scores(a, ["AlphaMissense"] + models, A1_KINDS, ["all_genes", "a5_cold_genes"], tag,
                       restrict="aim1_common")
    else:
        log("A6.3 skipped: Aim 4 model columns not present (aim4v2 has not completed)")
    # A6.2: sibling negative control at every level (A2 specification)
    a.results, a.cohort = P(base, "replication_unaffected"), "Unaffected"
    run_scores(a, ["AlphaMissense"] + list(TOOLS), A1_KINDS, ["all_genes"], "sib", restrict="aim1_common")
    # A6.1: reruns of the A1, A2 and A3 nested runs, saving replicates; keep the committed tables for comparison
    for fname, sub in A6_ORIG:
        src = P(base, sub, fname)
        if os.path.exists(src) and not os.path.exists(P(base, sub, "a6orig_" + fname)):
            shutil.copy(src, P(base, sub, "a6orig_" + fname))
    a.results, a.cohort = P(base, "replication_affected"), "Affected"
    run_scores(a, ["AlphaMissense"] + list(TOOLS), A1_KINDS, ["all_genes"], "rep", restrict="aim1_common")
    a.results, a.cohort = base, None
    run_scores(a, [t for t in TOOLS54 if t not in DDD_ONLY], A1_KINDS, ["all_genes", "ddg2p_monoallelic"], "a1v54",
               restrict="aim1_common")
    run_scores(a, ["AlphaMissense"] + list(TOOLS), A1_KINDS, ["all_genes", "ddg2p_monoallelic"], "aim1tools",
               restrict="aim1_common")
    cmd_a6_report(a)


def _r_ci(reps, s, k, st):
    """Point estimate and 95% percentile interval of r = ln LRw(k) / ln LRm(A) from paired replicates."""
    A = reps.get(f"{s}|gene|{st}|LRm")
    W = A if k == "A" else reps.get(f"{s}|{k}|{st}|LRw")
    if A is None or W is None:
        return np.nan, np.nan, np.nan, None
    with np.errstate(all="ignore"):
        r = np.where((A > 1) & (W > 0), np.log(W) / np.log(A), np.nan)
    b = r[1:][np.isfinite(r[1:])]
    lo, hi = (np.percentile(b, [2.5, 97.5]) if len(b) > 20 else (np.nan, np.nan))
    return float(r[0]), float(lo), float(hi), r


def cmd_a6_report(a):
    base = a.results
    out = []
    # A6.1: equality with committed tables, then retention CIs and the readout-1 share interval
    for fname, sub in A6_ORIG:
        o, n = P(base, sub, "a6orig_" + fname), P(base, sub, fname)
        if os.path.exists(o) and os.path.exists(n):
            x, y = pd.read_csv(o, sep="\t"), pd.read_csv(n, sep="\t")
            num = x.select_dtypes("number").columns
            same = len(x) == len(y) and np.allclose(x[num].to_numpy(float), y[num].to_numpy(float), rtol=1e-9,
                                                     atol=1e-12, equal_nan=True)
            out.append(f"A6.1 rerun of {sub + '/' if sub else ''}{fname}: point estimates "
                       f"{'IDENTICAL to the committed table' if same else 'DIFFER from the committed table (defect)'}")
    rows = []
    for tag, sub, st in (("aim1tools", "", "all_genes"), ("rep", "replication_affected", "all_genes"),
                         ("a1v54", "", "all_genes"), ("sib", "replication_unaffected", "all_genes")):
        pth = P(base, sub, f"{tag}_top2reps.npz")
        if not os.path.exists(pth):
            continue
        z = dict(np.load(pth))
        scores = sorted({kk.split("|")[0] for kk in z})
        share_reps = []
        for s in scores:
            A = z.get(f"{s}|gene|{st}|LRm")
            if A is None:
                continue
            row = {"set": tag, "score": s, "A": float(A[0])}
            for lvl, k in (("B", "gene"), ("C", "rmc"), ("D1", "pseudo"), ("D2", "win50"), ("E1", "pseudo_rmc"),
                           ("E2", "win50_rmc")):
                e, lo, hi, rr = _r_ci(z, s, k, st)
                row.update({f"r_{lvl}": e, f"r_{lvl}_lo": lo, f"r_{lvl}_hi": hi})
                if lvl == "D2" and s != "AlphaMissense" and A[0] >= PP3[0][1] and rr is not None:
                    share_reps.append(rr < 0.5)
            rows.append(row)
        if share_reps and tag != "sib":
            sh = np.mean(np.vstack(share_reps), 0)
            out.append(f"A6.1 {tag}: readout-1 share r(D2) < 0.5 = {sh[0]:.2f} [95% CI {np.percentile(sh[1:], 2.5):.2f}, "
                       f"{np.percentile(sh[1:], 97.5):.2f}] over {len(share_reps)} informative tools")
    if rows:
        ci = pd.DataFrame(rows)
        ci.to_csv(P(base, "a6_retention_ci.tsv"), sep="\t", index=False)
        for tag in ("aim1tools", "rep", "a1v54"):
            c = ci[(ci["set"] == tag) & (ci["score"] == "AlphaMissense")]
            if len(c):
                c = c.iloc[0]
                out.append(f"A6.1 {tag}: AlphaMissense r(B) {c.r_B:.2f} [{c.r_B_lo:.2f}, {c.r_B_hi:.2f}], r(D2) {c.r_D2:.2f} "
                           f"[{c.r_D2_lo:.2f}, {c.r_D2_hi:.2f}], r(E2) {c.r_E2:.2f} [{c.r_E2_lo:.2f}, {c.r_E2_hi:.2f}]")
    # A4.2 paired difference r(rand50) - r(D2)
    pth = P(base, "a4ctrl_top2reps.npz")
    if os.path.exists(pth):
        z = dict(np.load(pth))
        diffs = []
        for s in sorted({kk.split("|")[0] for kk in z}):
            _, _, _, r1 = _r_ci(z, s, "rand50", "all_genes")
            _, _, _, r2 = _r_ci(z, s, "win50", "all_genes")
            if r1 is not None and r2 is not None and s != "AlphaMissense":
                diffs.append(r1 - r2)
        if diffs:
            med = np.nanmedian(np.vstack(diffs), 0)
            out.append(f"A6.1 A4.2 median over tools of r(rand50) - r(D2): {med[0]:.3f} "
                       f"[95% CI {np.nanpercentile(med[1:], 2.5):.3f}, {np.nanpercentile(med[1:], 97.5):.3f}]")
    # A6.2 sibling pass criterion
    pth = P(base, "replication_unaffected", "sib_top2pct.tsv")
    if os.path.exists(pth):
        t = pd.read_csv(pth, sep="\t")
        t = t[t["stratum"] == "all_genes"]
        g = t[t["unit"] == "gene"]
        n_a = int((g["lo_sys"] >= PP3[0][1]).sum())
        n_b = int((g["lo_sys_w"] >= PP3[0][1]).sum())
        n_d = int((t.loc[t["unit"] == "win50", "lo_sys_w"] >= PP3[0][1]).sum())
        out.append(f"A6.2 sibling negative control: marginal top-2% lower bound >= 2.406 for {n_a} of {g['score'].nunique()} "
                   f"predictors -> {'PASS' if n_a <= 2 else 'FAIL'} (<= 2 pre-specified); level B {n_b}, level D2 {n_d}")
    # A6.3 models in autism: paired comparison v2-GBDT vs AlphaMissense (within gene, top 2%)
    for tag, sub in (("m4rep", "replication_affected"), ("m4sib", "replication_unaffected")):
        pth = P(base, sub, f"{tag}_top2reps.npz")
        if not os.path.exists(pth):
            continue
        z = dict(np.load(pth))
        for st in ("all_genes", "a5_cold_genes"):
            for m_name in list(M4_COLS) + list(M4V3):
                x, y = z.get(f"{m_name}|gene|{st}|LRw"), z.get(f"AlphaMissense|gene|{st}|LRw")
                if x is None or y is None:
                    continue
                with np.errstate(all="ignore"):
                    dl = np.log(x) - np.log(y)
                b = dl[1:][np.isfinite(dl[1:])]
                lo, hi = np.exp(np.percentile(b, [2.5, 97.5])) if len(b) > 20 else (np.nan, np.nan)
                verdict = ("model better" if lo > 1 else ("AlphaMissense better" if hi < 1 else "no significant difference"))
                prim = " [A6 PRIMARY]" if (tag == "m4rep" and st == "all_genes" and m_name == "Aim4v2_GBDT") else ""
                out.append(f"A6.3 {tag} {st}: within-gene top-2% LR {m_name} / AlphaMissense = {np.exp(dl[0]):.2f} "
                           f"[{lo:.2f}, {hi:.2f}] -> {verdict}{prim}")
        t = pd.read_csv(P(base, sub, f"{tag}_top2pct.tsv"), sep="\t")
        g = t[(t["unit"] == "gene") & (t["stratum"] == "all_genes")]
        for r in g.itertuples():
            out.append(f"A6.3 {tag} {r.score}: marginal top-2% {r.LRm_top2:.2f} (lower bound {r.lo_sys:.2f}); within gene "
                       f"{r.LRw_top2:.2f} [{r.lo95_w:.2f}, {r.hi95_w:.2f}]")
    open(P(base, "a6_summary.txt"), "w").write("\n".join(out) + "\n")
    log("A6:\n" + "\n".join(out))


# ======================================================================== amendment A4.1: ClinVar clinical impact
# Published ClinGen intervals (Bergquist et al. 2025 Genet Med Table 1; REVEL/VEST4 as Pejaver et al. 2022):
# lower bounds of the +1, +2, +3, +4 point intervals, on the raw score.
A4_TOOLS = {"AlphaMissense": ("am", [0.792, 0.906, 0.972, 0.990]),
            "REVEL": ("raw_REVEL_score", [0.644, 0.773, 0.879, 0.932]),
            "VEST4": ("raw_VEST4_score", [0.764, 0.861, 0.909, 0.965])}
A4_REVSTAT = {"criteria_provided,_single_submitter", "criteria_provided,_multiple_submitters,_no_conflicts",
              "reviewed_by_expert_panel", "practice_guideline"}
A4_CLASS = {"Pathogenic": "P/LP", "Likely_pathogenic": "P/LP", "Pathogenic/Likely_pathogenic": "P/LP",
            "Benign": "B/LB", "Likely_benign": "B/LB", "Benign/Likely_benign": "B/LB",
            "Uncertain_significance": "VUS"}


def read_clinvar_missense(path):
    """ClinVar GRCh38 VCF -> autosomal missense SNVs, >= 1 star, P/LP, B/LB or VUS (plan section 3e, A4.1)."""
    rows, meta = [], {}
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("##"):
                if ln.startswith("##fileDate"):
                    meta["fileDate"] = ln.strip().split("=", 1)[1]
                continue
            if ln.startswith("#"):
                continue
            f = ln.split("\t", 8)
            if f[0] not in AUTOSOMES or len(f[3]) != 1 or len(f[4]) != 1 or f[4] not in B2I:
                continue
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            if "missense_variant" not in info.get("MC", ""):
                continue
            cls = A4_CLASS.get(info.get("CLNSIG", ""))
            if cls is None or info.get("CLNREVSTAT", "") not in A4_REVSTAT:
                continue
            rows.append((f[0], int(f[1]), f[3], f[4], cls, f[2]))
    d = pd.DataFrame(rows, columns=["chrom", "pos38", "ref", "alt", "cls", "clinvar_id"])
    d["key"] = d["pos38"].to_numpy(np.int64) * 4 + d["alt"].map(B2I).to_numpy(np.int64)
    return d, meta


def cmd_clinvar(a):
    import hashlib
    from sklearn.metrics import roc_auc_score
    os.makedirs(a.results, exist_ok=True)
    cv, meta = read_clinvar_missense(a.clinvar)
    md5 = hashlib.md5(open(a.clinvar, "rb").read()).hexdigest()
    log(f"ClinVar {os.path.basename(a.clinvar)} fileDate {meta.get('fileDate')} md5 {md5}: "
        f"{len(cv):,} autosomal missense SNVs (>= 1 star; P/LP, B/LB, VUS)")
    cols = ["key", "tx", "rmc_oe", "ccr_pct"] + [c for c, _ in A4_TOOLS.values()]
    parts = []
    for c, u in site_chunks(a, cols=cols):
        m = u.merge(cv.loc[cv["chrom"] == c, ["key", "cls", "clinvar_id"]], on="key", how="inner")
        m["chrom"] = c
        parts.append(m)
    x = pd.concat(parts, ignore_index=True)
    log(f"matched to analysed universe sites (primary filters): {len(x):,} "
        f"({len(x) / len(cv):.1%}); by class {x['cls'].value_counts().to_dict()}")
    x["pm1like"] = (x["rmc_oe"] < 0.4) | (x["ccr_pct"] >= 95)

    nest = pd.read_csv(P(a.results, "a1_nested.tsv"), sep="\t")
    nest = nest[nest["stratum"] == "all_genes"].set_index("score")
    thr = pd.read_csv(P(a.results, "aim1tools_thresholds.tsv"), sep="\t")
    rows, lines = [], []
    for tool, (col, cuts) in A4_TOOLS.items():
        s = pd.to_numeric(x[col], errors="coerce").to_numpy(float)
        pub = np.select([s >= cuts[3], s >= cuts[2], s >= cuts[1], s >= cuts[0]], [4, 3, 2, 1], 0).astype(float)
        pub[np.isnan(s)] = np.nan
        r_e2 = float(nest.loc[tool, "r_E2"])
        pm1 = np.where(x["pm1like"], 2.0, 0.0)
        # rule 3: this study's within-unit E2 thresholds (rankscore units for dbNSFP tools; raw for AlphaMissense)
        t3 = thr[(thr["score"] == tool) & (thr["unit"] == "win50_rmc") & (thr["stratum"] == "all_genes") &
                 (thr["kind"] == "within_unit")]
        dn_cols = ["PP3_Supporting", "PP3_Moderate", "PP3_Strong", "PP3_VeryStrong"]
        dn_cuts = [float(t3[k].iloc[0]) if len(t3) and pd.notna(t3[k].iloc[0]) else np.nan for k in dn_cols]
        dn_avail = [k for k, v in zip(dn_cols, dn_cuts) if np.isfinite(v)]
        rule3 = np.zeros(len(x)) if not dn_avail else np.nan * np.zeros(len(x))   # rank-unit conversion only if needed
        if dn_avail:
            die("E2 thresholds exist for " + tool + ": rank-unit conversion not implemented (not expected per A1 results)")
        tot1u = pub + pm1
        tot1c = np.where(x["pm1like"], np.minimum(tot1u, 4.0), tot1u)
        tot2 = pub * r_e2 + pm1
        tot3 = rule3 + pm1
        for cls in ("P/LP", "VUS", "B/LB"):
            k = (x["cls"] == cls).to_numpy() & np.isfinite(pub)
            kp = k & x["pm1like"].to_numpy()
            pp3 = k & (pub >= 1)
            rows.append({
                "tool": tool, "cls": cls, "n": int(k.sum()), "n_pm1like": int(kp.sum()),
                "pp3_any": int(pp3.sum()),
                "share_pp3_in_pm1like": float((pp3 & kp).sum() / max(pp3.sum(), 1)),
                "pm1like_mean_pub_capped": float(tot1c[kp].mean()) if kp.any() else np.nan,
                "pm1like_mean_conditional": float(tot2[kp].mean()) if kp.any() else np.nan,
                "pm1like_mean_denovo_thr": float(tot3[kp].mean()) if kp.any() else np.nan,
                "pm1like_ge4_pub_to_lt4_cond": int((kp & (tot1c >= 4) & (tot2 < 4)).sum()),
                "pm1like_ge4_pub_capped": int((kp & (tot1c >= 4)).sum()),
                "ge6_uncapped": int((k & (tot1u >= 6)).sum()), "ge6_capped": int((k & (tot1c >= 6)).sum()),
                "ge6_conditional": int((k & (tot2 >= 6)).sum()), "r_E2": r_e2})
        ok = np.isfinite(pub) & x["cls"].isin(["P/LP", "B/LB"]).to_numpy()
        yb = (x["cls"] == "P/LP").to_numpy()
        for nm, sub in (("pm1like", x["pm1like"].to_numpy()), ("other", ~x["pm1like"].to_numpy())):
            kk = ok & sub
            if yb[kk].any() and (~yb[kk]).any():
                lines.append(f"{tool} AUC P/LP vs B/LB [{nm}, n={kk.sum():,}]: published points "
                             f"{roc_auc_score(yb[kk], pub[kk]):.3f}; conditional points "
                             f"{roc_auc_score(yb[kk], pub[kk] * r_e2):.3f} (identical ranking; scale only)")
    out = pd.DataFrame(rows)
    out.to_csv(P(a.results, "a4_clinvar_impact.tsv"), sep="\t", index=False)
    x.drop(columns=[c for c, _ in A4_TOOLS.values() if c != "am"]).to_csv(
        P(a.results, "a4_clinvar_matched.tsv.gz"), sep="\t", index=False)
    txt = [f"ClinVar {os.path.basename(a.clinvar)} (fileDate {meta.get('fileDate')}, md5 {md5}); "
           f"{len(x):,} missense variants matched to analysed sites; PM1-like {x['pm1like'].mean():.1%}"]
    for _, r in out.iterrows():
        txt.append(f"{r.tool:13s} {r.cls:5s} n={r.n:6,d} PP3>=Supp {r.pp3_any:6,d} of which PM1-like "
                   f"{r.share_pp3_in_pm1like:.1%} | PM1-like mean PP3+PM1 pts: published(capped) "
                   f"{r.pm1like_mean_pub_capped:.2f} -> conditional {r.pm1like_mean_conditional:.2f} "
                   f"(de novo thr {r.pm1like_mean_denovo_thr:.2f}); >=4 -> <4: {r.pm1like_ge4_pub_to_lt4_cond:,d} of "
                   f"{r.pm1like_ge4_pub_capped:,d} | PP3+PM1 alone >=6 (LP): uncapped {r.ge6_uncapped:,d}, "
                   f"capped {r.ge6_capped:,d}, conditional {r.ge6_conditional:,d}")
    txt += lines
    open(P(a.results, "a4_clinvar_summary.txt"), "w").write("\n".join(txt) + "\n")
    log("A4.1 ClinVar clinical impact:\n" + "\n".join(txt))


# ======================================================================== amendment A4.2: position-unit controls
A4_KINDS = ["gene", "win50", "win50s", "rand50", "dom", "dom_rmc"]


def cmd_a4ctrl(a):
    """AlphaMissense + the informative A1 tools (DD marginal top-2% LR >= 2.406), DD cohort, all genes."""
    nest = pd.read_csv(P(a.results, "a1_nested.tsv"), sep="\t")
    nest = nest[nest["stratum"] == "all_genes"]
    tools = [s for s in nest.loc[nest["A"] >= PP3[0][1], "score"] if s != "AlphaMissense"]
    if a.part in ("nodom", "all"):
        run_scores(a, ["AlphaMissense"] + tools, [k for k in A4_KINDS if not k.startswith("dom")], ["all_genes"], "a4ctrl")
    if a.part in ("dom", "all"):
        run_scores(a, ["AlphaMissense"] + tools, ["gene", "dom", "dom_rmc"], ["all_genes"], "a4ctrl_dom")
    t = pd.read_csv(P(a.results, "a4ctrl_top2pct.tsv"), sep="\t")
    if os.path.exists(P(a.results, "a4ctrl_dom_top2pct.tsv")):
        td = pd.read_csv(P(a.results, "a4ctrl_dom_top2pct.tsv"), sep="\t")
        t = pd.concat([t, td[td["unit"] != "gene"]], ignore_index=True)
    t = t[t["stratum"] == "all_genes"]
    rows = []
    for s, g in t.groupby("score", sort=False):
        g = g.set_index("unit")
        A = float(g.loc["gene", "LRm_top2"])
        r = {"score": s, "A": A}
        for k in A4_KINDS:
            v = float(g.loc[k, "LRw_top2"]) if k in g.index else np.nan
            r[k] = v
            r["r_" + k] = np.log(v) / np.log(A) if (v > 0 and A > 1) else np.nan
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(P(a.results, "a4ctrl_nested.tsv"), sep="\t", index=False)
    inf = df[(df["score"] != "AlphaMissense") & (df["A"] >= PP3[0][1])]
    d = (inf["r_rand50"] - inf["r_win50"]).dropna()
    txt = [f"informative tools: {len(inf)} (+ AlphaMissense reported separately)",
           "median r: " + ", ".join(f"{k} {inf['r_' + k].median():.3f}" for k in A4_KINDS[1:]) +
           f" (B gene {inf['r_gene'].median():.3f})",
           f"paired r(rand50) - r(win50): median {d.median():.3f} [IQR {d.quantile(.25):.3f}, {d.quantile(.75):.3f}]; "
           f"positive for {int((d > 0).sum())} of {len(d)} tools",
           f"r(win50s) - r(win50): median {(inf['r_win50s'] - inf['r_win50']).median():.3f}; "
           f"r(dom) - r(win50): median {(inf['r_dom'] - inf['r_win50']).median():.3f}"]
    am = df[df["score"] == "AlphaMissense"]
    if len(am):
        txt.append("AlphaMissense: " + ", ".join(f"{k} {float(am[k].iloc[0]):.2f}" for k in ["A"] + A4_KINDS))
    open(P(a.results, "a4ctrl_summary.txt"), "w").write("\n".join(txt) + "\n")
    log("A4.2 position controls:\n" + "\n".join(txt))


A14_SCORES = {"Aim4v5_GBDT": "m4_v5_gbdt", "Aim4v4_GBDT": "m4_v4_gbdt", "Aim4v3_GBDT": "m4_v3_gbdt", "AlphaMissense": "am",
              "REVEL": "REVEL_rankscore", "VEST4": "VEST4_rankscore", "popEVE": "v54_popEVE_converted_rankscore",
              "ESM1b": "v54_ESM1b_converted_rankscore", "MisFit_S": "v54_MisFit_S_rankscore",
              "PrimateAI": "PrimateAI_rankscore", "CADD": "CADD_raw_rankscore"}


def _a14_gene_auc(gid, y, s):
    """Per gene: AUROC (average ranks, ties 0.5) and number of P/LP x B/LB pairs; genes need >= 3 of each class."""
    from scipy.stats import rankdata
    out = {}
    df = pd.DataFrame({"g": gid, "y": y, "s": s})
    for g, sub in df.groupby("g", sort=False):
        n1 = int(sub["y"].sum()); n0 = len(sub) - n1
        if n1 < 3 or n0 < 3:
            continue
        r = rankdata(sub["s"].to_numpy())
        out[g] = ((r[sub["y"].to_numpy() == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0), n1 * n0)
    return out


def cmd_a14(a):
    """Amendment A14: ClinVar P/LP vs B/LB, within-gene AUROC of the de novo-trained models vs AlphaMissense and others."""
    m = pd.read_csv(P(a.results, "a4_clinvar_matched.tsv.gz"), sep="\t", dtype={"chrom": str}, usecols=["key", "cls", "chrom"])
    m = m[m["cls"].isin(["P/LP", "B/LB"])]
    meta = pd.read_csv(P(a.work, "dnv_grch38.tsv"), sep="\t", dtype={"chrom": str}, usecols=["chrom", "pos38", "alt38"])
    meta = meta[meta["alt38"].isin(list(gate.B2I))]
    dkey = set(zip(meta["chrom"].str.replace("chr", ""), meta["pos38"].astype(np.int64) * 4 + meta["alt38"].map(gate.B2I)))
    tinfo = load_tinfo(a)
    v5 = os.path.exists(P(a.fs, "col", "_m4v5_done"))
    stem = "a15" if v5 else "a14"
    rows = []
    for c in gate.AUTOSOMES:
        mc = m[m["chrom"].str.replace("chr", "") == c]
        if mc.empty:
            continue
        od = P(a.fs, "col", f"chr{c}")
        keys = np.load(P(od, "key.npy"))
        ix = pd.Index(keys).get_indexer(mc["key"].to_numpy(np.int64))
        ok = ix >= 0
        mc, ix = mc[ok], ix[ok]
        drop = np.array([(c, int(k)) in dkey for k in mc["key"]], dtype=bool)
        mc, ix = mc[~drop], ix[~drop]
        d = {"key": mc["key"].to_numpy(), "y": (mc["cls"] == "P/LP").to_numpy(np.int8),
             "tx": np.load(P(od, "tx.npy"), mmap_mode="r")[ix]}
        for name, col in A14_SCORES.items():
            f = P(od, col + ".npy")
            d[name] = np.asarray(np.load(f, mmap_mode="r")[ix], dtype=np.float64) if os.path.exists(f) else np.full(len(ix), np.nan)
        rows.append(pd.DataFrame(d))
        log(f"A14 chr{c}: {len(ix)} ClinVar sites kept, {int(drop.sum())} DD de novo sites removed")
    df = pd.concat(rows, ignore_index=True)
    df = df[df["tx"] >= 0]
    df["gene"] = tinfo["gene"].reindex(df["tx"].to_numpy()).to_numpy()
    df["ddg2p"] = tinfo["ddg2p"].reindex(df["tx"].to_numpy()).fillna(False).to_numpy(bool)
    rng = np.random.default_rng(20261006)
    out, lines = [], []
    for stratum, smask in (("all_genes", np.ones(len(df), bool)), ("ddg2p_monoallelic", df["ddg2p"].to_numpy())):
        prim = (["Aim4v5_GBDT"] if v5 else []) + ["Aim4v4_GBDT", "Aim4v3_GBDT", "AlphaMissense"]
        sets = [("primary", prim)] + \
               [(f"pair_{t}", [t, "AlphaMissense"]) for t in A14_SCORES if t not in prim and t != "Aim4v5_GBDT"]
        for setname, tools in sets:
            sub = df[smask & df[tools].notna().all(axis=1).to_numpy()]
            aucs = {t: _a14_gene_auc(sub["gene"].to_numpy(), sub["y"].to_numpy(), sub[t].to_numpy()) for t in tools}
            genes = sorted(set.intersection(*[set(v) for v in aucs.values()]))
            if len(genes) < 5:
                continue
            A = {t: np.array([aucs[t][g][0] for g in genes]) for t in tools}
            W = np.array([aucs[tools[0]][g][1] for g in genes], dtype=float)
            C = rng.multinomial(len(genes), np.full(len(genes), 1 / len(genes)), size=1000).astype(float)
            pooled = {t: (A[t] * W).sum() / W.sum() for t in tools}
            boot = {t: (C * (A[t] * W)).sum(1) / (C * W).sum(1) for t in tools}
            for t in tools:
                lo, hi = np.percentile(boot[t], [2.5, 97.5])
                dl = dh = dp = np.nan
                if t != "AlphaMissense":
                    dd_ = boot[t] - boot["AlphaMissense"]; dp = pooled[t] - pooled["AlphaMissense"]
                    dl, dh = np.percentile(dd_, [2.5, 97.5])
                out.append(dict(stratum=stratum, set=setname, tool=t, n_genes=len(genes), n_sites=len(sub),
                                n_plp=int(sub["y"].sum()), auc_pooled=pooled[t], lo95=lo, hi95=hi,
                                auc_gene_mean=A[t].mean(), d_vs_AM=dp, d_lo95=dl, d_hi95=dh))
                if t != "AlphaMissense":
                    tag = " [PRIMARY]" if (stratum == "ddg2p_monoallelic" and setname == "primary" and t == ("Aim4v5_GBDT" if v5 else "Aim4v4_GBDT")) else ""
                    lines.append(f"{stratum:18s} {setname:16s} {t:14s} AUC {pooled[t]:.3f} vs AM {pooled['AlphaMissense']:.3f}  "
                                 f"d {dp:+.3f} [{dl:+.3f}, {dh:+.3f}]  genes {len(genes)} sites {len(sub)}{tag}")
    pd.DataFrame(out).to_csv(P(a.results, f"{stem}_clinvar_withingene.tsv"), sep="\t", index=False, float_format="%.4f")
    with open(P(a.results, f"{stem}_summary.txt"), "w") as f:
        f.write("A14: ClinVar P/LP vs B/LB, within-gene pooled AUROC (genes with >= 3 of each), DD de novo sites removed; "
                "gene bootstrap 1000, seed 20261006\n" + "\n".join(lines) + "\n")
    log("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=P(HERE, "data"))
    ap.add_argument("--work", default=P(HERE, "work"))
    ap.add_argument("--fs", default=P(HERE, "work", "fullstudy"))
    ap.add_argument("--results", default=P(HERE, "results_fullstudy"))
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("universe", "scores", "regions", "colstore", "scores54", "scores54_stream", "scores54_tail"):
        s = sub.add_parser(name)
        s.add_argument("--chroms", default="1-22")
        s.add_argument("--force", action="store_true")
        s.add_argument("--from-piece", default="part.7")
    for name in ("validate", "aim5", "aim1", "aim1x", "aim1tools", "aim3", "aim3_report", "a1_report", "aim4", "replicate", "v54"):
        sub.add_parser(name)
    for name in ("dom54_stream", "af_cols", "aim4v2", "aim4v3", "aim4v4", "public3", "overlap", "a10_post", "a13", "a14", "aim4v5", "m4_nested", "a6", "a6_report"):
        sub.add_parser(name)
    s = sub.add_parser("a4ctrl")
    s.add_argument("--part", default="all", choices=["all", "nodom", "dom", "report"])
    s = sub.add_parser("a11")
    s.add_argument("--part", required=True, choices=["strata", "sweep", "cpg_down", "cpg_up", "noclinvar"])
    s.add_argument("--clinvar", default="/mnt/c/Users/saite/OneDrive/Desktop/Neuro/G-Test-LOCKIN/ndd_pp3_gate/ndd_pp3_gate/work_clinvar/clinvar.vcf.gz")
    s = sub.add_parser("a12")
    s.add_argument("--part", required=True, choices=["decomp", "classify", "lto", "allr", "plddt"])
    s = sub.add_parser("a11sim")
    s.add_argument("--reps", type=int, default=200)
    s = sub.add_parser("clinvar")
    s.add_argument("--clinvar", default="/mnt/c/Users/saite/OneDrive/Desktop/Neuro/G-Test-LOCKIN/ndd_pp3_gate/ndd_pp3_gate/work_clinvar/clinvar.vcf.gz")
    s = sub.add_parser("aim2")
    s.add_argument("--tools", action="store_true", help="also every Aim 3 tool (needs 'scores')")
    a = ap.parse_args()
    {"universe": cmd_universe, "scores": cmd_scores, "colstore": cmd_colstore, "scores54": cmd_scores54, "scores54_stream": cmd_scores54_stream, "scores54_tail": cmd_scores54_tail, "regions": cmd_regions, "validate": cmd_validate,
     "aim5": cmd_aim5, "aim1": cmd_aim1, "aim2": cmd_aim2, "aim3": cmd_aim3, "aim4": cmd_aim4, "aim1x": cmd_aim1x, "aim1tools": cmd_aim1tools, "aim3_report": cmd_aim3_report, "a1_report": cmd_a1_report, "replicate": cmd_replicate, "v54": cmd_v54, "clinvar": cmd_clinvar, "dom54_stream": cmd_dom54_stream, "a4ctrl": cmd_a4ctrl, "af_cols": cmd_af_cols, "aim4v2": cmd_aim4v2, "m4_nested": cmd_m4_nested, "aim4v3": cmd_aim4v3, "aim4v4": cmd_aim4v4, "public3": cmd_public3, "overlap": cmd_overlap, "a10_post": cmd_a10_post, "a11": cmd_a11, "a11sim": cmd_a11sim, "a12": cmd_a12, "a13": cmd_a13, "a14": cmd_a14, "aim4v5": cmd_aim4v5, "a6": cmd_a6, "a6_report": cmd_a6_report}[a.cmd](a)


if __name__ == "__main__":
    main()
