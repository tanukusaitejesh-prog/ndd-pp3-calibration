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


def site_chunks(a, chroms=AUTOSOMES, near="exclude", cls=1, extra=()):
    """Per chromosome, universe sites of class `cls` passing the primary filters, with optional joined tables."""
    for c in chroms:
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


def unit_of(u, kind):
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
    if kind == "pseudo":
        return (18397 + 1) * PSEUDO_K
    if kind == "gene":
        return 18397 + 1
    if kind == "ccr":
        return (18397 + 1) * 6
    return len(pd.read_pickle(P(a.fs, "regions", "rmc_regions.pkl"))) + 1


def unit_tx(kind, a, nu):
    """Transcript of each unit (for strata defined on genes)."""
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

    def extra(self):
        ex = []
        if any(s not in ("AlphaMissense", "AM_gate") for s in self.scores):
            ex.append("scores")
        if any(k in ("rmc", "ccr") for k in self.kinds) or self.restrict:
            ex.append("regions")
        return tuple(ex)

    def score_of(self, u, s):
        return u["am"].to_numpy(np.float64) if s == "AlphaMissense" else u[TOOLS[s]].to_numpy(np.float64)

    def nbins(self, s):
        return len(gate.LR_EDGES_PCT) - 1 if s == "AM_gate" else N_QBINS

    def bins_of(self, u, s):
        if s == "AM_gate":
            f2c = np.searchsorted(np.array(gate.LR_EDGES_PCT), np.arange(100), side="right") - 1
            return f2c[u["am_bin"].to_numpy().astype(int)].astype(np.int64)
        return to_bins(self.score_of(u, s), self.edges[s])

    def edges_pass(self):
        hs = {s: [] for s in self.scores if s not in self.edges and s != "AM_gate"}
        if not hs:
            return
        chunks = {s: [] for s in hs}
        for c, u in site_chunks(self.a, near=self.near, extra=self.extra()):
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
        grid = len(h)
        cdf = np.cumsum(h) / h.sum()
        xs = np.linspace(0, 1, grid)
        inner = xs[np.minimum(np.searchsorted(cdf, np.linspace(0, 1, nq + 1)[1:-1]), grid - 1)]
        return np.concatenate([[-np.inf], inner, [np.inf]])

    def _restrict(self, u):
        if self.restrict == "aim1_common":
            u = u[(u["rmc_unit"] >= 0) & u["ccr_pct"].notna()]
        return u

    def mass_pass(self, nunits):
        """nunits: {kind: number of units}. Opportunity is summed into sparse matrices chromosome by chromosome."""
        dm = self.d_mis
        self.M = {(s, k): None for s in self.scores for k in self.kinds}
        for c, u in site_chunks(self.a, near=self.near, extra=self.extra()):
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
        labels = [(0.0 if lo_ == 0 else float(edges[lo_]), 1.0 if hi_ == N_QBINS else float(edges[hi_]))
                  for lo_, hi_ in zip(Q_COARSE[:-1], Q_COARSE[1:])]
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
        cs = bin_centre_scores(edges)[centres]
        sw = summ(lrw[key]); lo_w = summ(lrw[key + "_dn"])["lo90"]; hi_w = summ(lrw[key + "_up"])["hi90"]
        for i in range(len(centres)):
            wrows.append(tag + [kind, cs[i], lrw["O"][0, i], lrw["E"][0, i], sw["est"][i], sw["lo95"][i], sw["hi95"][i],
                                lo_w[i], hi_w[i]])
        th = window_thresholds(lo_w, hi_w, cs)
        trows.append(tag + [kind] + [th[k] for k in THR_KEYS])
    return rows, wrows, trows, res


COARSE_COLS = ["kind", "score_lo", "score_hi", "O", "E", "X", "LR", "lo95", "hi95", "lo_sys", "hi_sys", "points", "evidence"]
WIN_COLS = ["kind", "centre_score", "O", "E", "LR", "lo95", "hi95", "lo_sys", "hi_sys"]
THR_KEYS = [f"PP3_{n}" for n, _ in PP3] + [f"BP4_{n}" for n, _ in PP3]


# ======================================================================== aims
def setup(a, near="exclude"):
    os.makedirs(a.results, exist_ok=True)
    d, n_prob = load_dnv(a)
    eng = Engine(n_prob)
    syn_mass, syn_p, _ = syn_inputs(a, d, near)
    d_mis = dnv_filter(d[d["status"] == "missense"], near)
    tinfo = load_tinfo(a)
    return d, eng, syn_mass, syn_p, d_mis, tinfo


def gene_masks(tinfo, nu, utx):
    """Strata defined on genes, as boolean masks over units (utx = transcript of each unit)."""
    t = tinfo.reindex(np.where(utx >= 0, utx, -1))
    t.index = np.arange(nu)
    ok = utx >= 0
    dd = t["ddg2p"].fillna(False).to_numpy(bool) & ok
    m = {"all_genes": ok, "ddg2p_monoallelic": dd,
         "shet_high": (t["shet_class"] == 1).to_numpy() & ok, "shet_low": (t["shet_class"] == 0).to_numpy() & ok}
    for lab in ("LoF", "nonLoF"):
        m[f"ddg2p_{lab}"] = dd & (t["mech12_bin"] == lab).to_numpy()
        m[f"allgenes_{lab}"] = ok & (t["mech123_bin"] == lab).to_numpy()
    for lab in ("lof", "dn", "gof"):
        m[f"ddg2p_{lab}"] = dd & (t["mech12"] == lab).to_numpy()
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


def top_arrays(res, lo_bin=980):
    lr = Engine.lr(res, [np.arange(lo_bin, N_QBINS)])
    return lr["LRm"][:, 0], lr["LRw"][:, 0]


def run_scores(a, scores, kinds, strata_names, tag, restrict=None, near="exclude", explicit_am=True, top_contrast=None):
    """Shared driver for Aims 1-3: windows + coarse tables + thresholds for each score x unit kind x stratum."""
    d, eng, syn_mass, syn_p, d_mis, tinfo = setup(a, near)
    nunits = {k: n_units(k, a) for k in kinds}
    utx = {k: unit_tx(k, a, nunits[k]) for k in kinds}
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
    paths = {x: P(a.results, f"{tag}_{x}.tsv") for x in ("coarse", "windows", "thresholds", "top2pct")}
    for pth in paths.values():
        if os.path.exists(pth):
            os.remove(pth)
    batch = 4
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
                        tops[st] = top_arrays(res)
                        lm, lw = summ(tops[st][0][:, None]), summ(tops[st][1][:, None])
                        write([[s, k, st, lm["est"][0], lm["lo95"][0], lm["hi95"][0], lw["est"][0], lw["lo95"][0], lw["hi95"][0],
                                points(np.array([lw["est"][0]]))[0]]],
                              ["score", "unit", "stratum", "LRm_top2", "lo95", "hi95", "LRw_top2", "lo95_w", "hi95_w", "points_w"],
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


def cmd_aim1tools(a):
    """Aim 1 estimands for every predictor (common site set; gene, gene x RMC region, gene x CCR class)."""
    run_scores(a, list(TOOLS), ["gene", "rmc", "ccr"], ["all_genes", "ddg2p_monoallelic"], "aim1tools",
               restrict="aim1_common")


# ======================================================================== Aim 4: null-contrast learning
AA = "ACDEFGHIKLMNPQRSTVWY"
PRED_LR_EDGES = [0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, np.inf]       # calibration bins of predicted LR


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
    for c, u in site_chunks(a, extra=()):
        tr = ~cold_tx[u["tx"].to_numpy()]
        mass_c[c] = float(u.loc[tr, PRIMARY_MODEL].sum())
    tot = sum(mass_c.values())
    draws = dict(zip(mass_c, rng.multinomial(n_null, np.array(list(mass_c.values())) / tot)))
    real_parts, null_parts, med_parts = [], [], []
    for c, u in site_chunks(a, extra=extra):
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
    for c, u in site_chunks(a, extra=extra):
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

    def mass(bins, nb, mask):
        ok = mask & (bins >= 0)
        flat = np.bincount(S["tx"][ok].astype(np.int64) * nb + bins[ok], weights=S["mr"][ok].astype(np.float64),
                           minlength=nu * nb)
        return sparse.csr_matrix(flat.reshape(nu, nb))

    def dmap_of(dd, bins):
        r = dd["key"].map(drow)
        ok = r.notna().to_numpy()
        r = r[ok].to_numpy(np.int64)
        out = dd[ok][["p", "tx", "centre", "key"]].copy()
        out["unit"] = S["tx"][r].astype(np.int64); out["bin"] = bins[r]
        return out

    nb_lr = len(PRED_LR_EDGES) - 1
    lr_bins = np.searchsorted(np.array(PRED_LR_EDGES[1:-1]), S["lr"], side="right").astype(np.int64)

    def qedges(v):
        o = np.argsort(v, kind="stable")
        cw = np.cumsum(S["mr"][o].astype(np.float64)); cw /= cw[-1]
        inner = v[o][np.minimum(np.searchsorted(cw, np.linspace(0, 1, N_QBINS + 1)[1:-1]), len(v) - 1)]
        del o, cw
        return np.concatenate([[-np.inf], inner.astype(np.float64), [np.inf]])
    lv = np.log1p(S["lr"].astype(np.float64))
    lr_q = to_bins(lv, qedges(lv))
    del lv
    am_q = to_bins(S["am"].astype(np.float64), qedges(S["am"].astype(np.float64)))
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=P(HERE, "data"))
    ap.add_argument("--work", default=P(HERE, "work"))
    ap.add_argument("--fs", default=P(HERE, "work", "fullstudy"))
    ap.add_argument("--results", default=P(HERE, "results_fullstudy"))
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("universe", "scores", "regions"):
        s = sub.add_parser(name)
        s.add_argument("--chroms", default="1-22")
        s.add_argument("--force", action="store_true")
    for name in ("validate", "aim5", "aim1", "aim1x", "aim1tools", "aim3", "aim3_report", "aim4"):
        sub.add_parser(name)
    s = sub.add_parser("aim2")
    s.add_argument("--tools", action="store_true", help="also every Aim 3 tool (needs 'scores')")
    a = ap.parse_args()
    {"universe": cmd_universe, "scores": cmd_scores, "regions": cmd_regions, "validate": cmd_validate,
     "aim5": cmd_aim5, "aim1": cmd_aim1, "aim2": cmd_aim2, "aim3": cmd_aim3, "aim4": cmd_aim4, "aim1x": cmd_aim1x, "aim1tools": cmd_aim1tools, "aim3_report": cmd_aim3_report}[a.cmd](a)


if __name__ == "__main__":
    main()
