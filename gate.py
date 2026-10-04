#!/usr/bin/env python3
"""
gate.py - pre-registered feasibility gate for the project

  "De novo excess in developmental-disorder trios as a ClinVar-free truth set for
   calibrating ACMG PP3/BP4 evidence for de novo missense variants"

The question the gate answers: with the published DDD/GeneDx/RUMC de novo list
(31,058 trios), Roulette mutation rates and AlphaMissense, can we (1) reproduce
the synonymous de novo count from the mutation model without bias across
mutability and protein context (null check), and (2) estimate likelihood ratios
for AlphaMissense precisely enough to map them onto the ClinGen PP3/BP4 scale?

Likelihood ratio for a score bin b (de novo missense, NDD cohort):
    E_b = c * sum(mu over all possible missense SNVs in bin b),  c = O_syn / sum(mu_syn)
    X_b = O_b - E_b                      (excess = pathogenic de novos)
    LR_b (marginal)    = [X_b / sum X] / [E_b / sum E]
    LR_b (within gene) = X_b / sum_g O_g * (E_gb / E_g)  - E_b   [gene-composition matched]

Run the subcommands in this order (README.md has the full walk-through):
    python gate.py check-inputs
    python gate.py prep-dnv      lift the published de novo list from GRCh37 to GRCh38
    python gate.py prep-am       split AlphaMissense_hg38 into compact per-chromosome files
    python gate.py prep-gtf      CDS models for AlphaMissense transcripts + region files
    bash get_roulette.sh         download Roulette rates, keep only coding sites
    python gate.py sites         every possible coding SNV: class, AM score, rates
    python gate.py run           null checks, likelihood ratios, thresholds, verdict
"""
import argparse
import glob
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ input files
DNV_FILE = os.path.join("DeNovoWEST", "input",
                        "DDD_RUMC_GDX_denovos_cadd_shet_wweights_minrep_indelfix_2025_09_18.txt.gz")
GENE_FILE = os.path.join("DeNovoWEST", "input", "extended_denovoWEST_results.tab")
DDG2P_FILE = os.path.join("DeNovoWEST", "input", "DDG2P_22_9_2019.csv")
AM_FILE = "AlphaMissense_hg38.tsv.gz"
GTF_FILE = "gencode.v32.annotation.gtf.gz"
CHAIN_FILE = "hg19ToHg38.over.chain.gz"

STUDY_TO_CENTRE = {"DDD": "DDD", "GDX": "GeneDx", "newGDX": "GeneDx", "RUMC": "RUMC"}
CENTRES = ["GeneDx", "DDD", "RUMC"]
AUTOSOMES = [str(i) for i in range(1, 23)]
BASES = "ACGT"
B2I = {b: i for i, b in enumerate(BASES)}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
# standard genetic code; codon index = 16*b1 + 4*b2 + b3 with A=0, C=1, G=2, T=3
CODON_TABLE = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
CT = np.frombuffer(CODON_TABLE.encode("ascii"), dtype=np.uint8)
STOP = ord("*")
RATE_MODELS = ["MR", "AR", "MG", "MC"]
CLS_NAMES = {0: "synonymous", 1: "missense", 2: "stop_gained", 3: "excluded_codon", 4: "internal_stop_tx"}
SPLICE_WIN = 3                      # bases at each internal exon edge flagged as near-splice
MR_LOG_MIN, MR_LOG_STEP, MR_NBINS = -12.0, 0.05, 320   # fine log10(rate) grid, used only for grouping

# ======================= PRE-REGISTERED CONSTANTS (fixed before the first real run) =======================
PRIMARY_MODEL = "MR"                # raw Roulette rate; AR, MG, MC are sensitivity analyses
SEED = 20261004
N_BOOT = 1000                       # Poisson bootstrap over probands
LR_EDGES_PCT = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 98, 99, 100]   # AlphaMissense bins (x100)
CUM_THRESH_PCT = [50, 80, 90, 95, 98, 99]
PP3 = [("Supporting", 2.406), ("Moderate", 5.790), ("Strong", 33.53), ("VeryStrong", 1124.0)]  # Pejaver 2022
SHET_HIGH = 0.15                    # DeNovoWEST s_het split
N_NULL_GROUPS = 10
NULL_PASS = 0.10                    # synonymous O/E tolerance: FAIL if a group's 95% CI lies wholly outside 1 +- 0.10
CONCORDANCE_MIN = 0.99              # enumerated vs AlphaMissense missense per transcript
QC_TX_PASS_MIN = 0.95
QC_DNV_MAP_MIN = 0.85
QC_REF_MISMATCH_MAX = 0.01
PRECISION_MAX_RATIO = 3.0           # hi95 / lo95 of LR in the bins with AM >= 0.90
DELTA_MIN = 0.05                    # systematic margin on expected counts used for evidence bounds:
                                    # delta = max(DELTA_MIN, largest |O/E - 1| in the gated synonymous null groups)
TOP_BINS = 4
# ==========================================================================================================


def die(msg):
    sys.exit("\nERROR: " + msg)


def log(msg):
    print(time.strftime("[%H:%M:%S] ") + msg, flush=True)


def P(*parts):
    return os.path.join(*parts)


def parse_chroms(s):
    out = []
    for part in str(s).replace(" ", "").split(","):
        part = part.replace("chr", "")
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-")
            out += [str(i) for i in range(int(a), int(b) + 1)]
        else:
            out.append(part)
    bad = [c for c in out if c not in AUTOSOMES]
    if bad:
        die(f"--chroms must be autosomes 1-22, got {bad}")
    return out


def ascii_codes(series):
    vals = series.tolist()
    return np.frombuffer("".join(vals).encode("ascii"), dtype=np.uint8).copy()


def cat_lookup(cat_series, mapping, default=-1):
    """Map a categorical Series through a dict, vectorised over its categories."""
    cats = cat_series.cat.categories
    table = np.array([mapping.get(x, default) for x in cats] + [default])
    codes = cat_series.cat.codes.to_numpy()
    return table[np.where(codes >= 0, codes, len(cats))]


def am_header(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for i, line in enumerate(f):
            if line.startswith("#CHROM") or line.startswith("CHROM\t"):
                return i, line.rstrip("\n").lstrip("#").split("\t")
            if i > 500:
                break
    die(f"could not find the '#CHROM' header line in {path}")


def sha256(path, limit=None):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ======================================================================== check-inputs
def cmd_check_inputs(a):
    ok = True

    def report(name, good, detail=""):
        nonlocal ok
        ok = ok and good
        print(f"  [{'OK' if good else 'FAIL'}] {name}  {detail}")

    print(f"Checking inputs in {a.data}")
    f = P(a.data, DNV_FILE)
    if os.path.exists(f):
        d = pd.read_csv(f, sep="\t", nrows=50, dtype=str)
        need = ["chrom", "pos", "ref", "alt", "id", "study", "cq", "symbol"]
        miss = [c for c in need if c not in d.columns]
        report("de novo list", not miss, f"missing columns {miss}" if miss else f"({os.path.basename(f)})")
    else:
        report("de novo list", False, f"not found: {f}  (run get_data.sh)")
    for name, rel, cols in [("gene table", GENE_FILE, ["symbol", "s_het"]),
                            ("DDG2P", DDG2P_FILE, ["gene symbol", "DDD category", "allelic requirement"])]:
        f = P(a.data, rel)
        if os.path.exists(f):
            d = pd.read_csv(f, sep="\t" if rel.endswith(".tab") else ",", nrows=5, dtype=str)
            miss = [c for c in cols if c not in d.columns]
            report(name, not miss, f"missing columns {miss}" if miss else "")
        else:
            report(name, False, f"not found: {f}")
    f = P(a.data, AM_FILE)
    if os.path.exists(f):
        try:
            skip, cols = am_header(f)
            need = ["CHROM", "POS", "REF", "ALT", "transcript_id", "protein_variant", "am_pathogenicity"]
            miss = [c for c in need if c not in cols]
            row = pd.read_csv(f, sep="\t", skiprows=skip + 1, header=None, names=cols, nrows=3, dtype=str)
            detail = f"header at line {skip + 1}; first row: {' '.join(row.iloc[0].astype(str).tolist())}"
            good = (not miss) and str(row.iloc[0]["CHROM"]).startswith("chr")
            report("AlphaMissense", good, f"missing columns {miss}" if miss else detail)
        except (OSError, EOFError) as e:
            report("AlphaMissense", False, f"cannot read ({e}); re-download")
    else:
        report("AlphaMissense", False, f"not found: {f}")
    f = P(a.data, GTF_FILE)
    if os.path.exists(f):
        body = []
        with gzip.open(f, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    body.append(line)
                if len(body) >= 50:
                    break
        good = bool(body) and all(len(l.split("\t")) == 9 for l in body) and any('transcript_id "' in l for l in body)
        report("GENCODE v32 GTF", good, "" if good else f"unexpected content: {body[0][:120] if body else 'empty'}")
    else:
        report("GENCODE v32 GTF", False, f"not found: {f}")
    f = P(a.data, CHAIN_FILE)
    if os.path.exists(f):
        with gzip.open(f, "rt") as fh:
            line = fh.readline()
        report("hg19ToHg38 chain", line.startswith("chain"), "" if line.startswith("chain") else line[:80])
    else:
        report("hg19ToHg38 chain", False, f"not found: {f}")
    try:
        import pyliftover  # noqa: F401
        report("pyliftover", True)
    except ImportError:
        report("pyliftover", False, "pip install pyliftover")
    rfiles = sorted(glob.glob(P(a.work, "roulette", "chr*.tsv.gz")))
    if rfiles:
        r = pd.read_csv(rfiles[0], sep="\t", header=None, nrows=5, dtype=str)
        good = r.shape[1] == 9 and r.iloc[0, 4] is not None and len(str(r.iloc[0, 4])) == 5
        report(f"Roulette extract ({len(rfiles)} chromosomes so far)", good,
               "" if good else f"unexpected first line: {r.iloc[0].tolist()}")
    else:
        print("  [--] Roulette extracts not present yet (made by get_roulette.sh)")
    print("\nAll inputs look right." if ok else "\nFix the FAIL lines above before continuing.")
    if not ok:
        sys.exit(1)


# ======================================================================== prep-dnv
def cmd_prep_dnv(a):
    try:
        from pyliftover import LiftOver
    except ImportError:
        die("pyliftover is not installed: pip install pyliftover")
    d = pd.read_csv(P(a.data, DNV_FILE), sep="\t", low_memory=False,
                    dtype={"chrom": str, "ref": str, "alt": str, "id": str, "study": str, "cq": str, "symbol": str})
    need = ["chrom", "pos", "ref", "alt", "id", "study", "cq", "symbol"]
    miss = [c for c in need if c not in d.columns]
    if miss:
        die(f"de novo file lacks columns {miss}")
    d = d[need].copy()
    d.insert(0, "dnv_idx", np.arange(len(d), dtype=np.int64))
    d["centre"] = d["study"].map(STUDY_TO_CENTRE)
    if d["centre"].isna().any():
        die(f"unknown study labels {sorted(d.loc[d['centre'].isna(), 'study'].unique())}")
    snv = (d["ref"].str.len().eq(1) & d["alt"].str.len().eq(1)
           & d["ref"].isin(list(BASES)) & d["alt"].isin(list(BASES))).to_numpy()
    auto = d["chrom"].isin(AUTOSOMES).to_numpy()
    lo = LiftOver(P(a.data, CHAIN_FILE))
    status = np.where(~snv, "not_snv", np.where(~auto, "not_autosomal", "")).astype(object)
    pos38 = np.full(len(d), -1, dtype=np.int64)
    ref38 = d["ref"].to_numpy(dtype=object).copy()
    alt38 = d["alt"].to_numpy(dtype=object).copy()
    chroms = d["chrom"].to_numpy(dtype=object)
    poss = d["pos"].to_numpy(dtype=np.int64)
    for i in np.flatnonzero(status == ""):
        res = lo.convert_coordinate("chr" + chroms[i], int(poss[i]) - 1)
        if not res:
            status[i] = "unmapped"
            continue
        if len(res) > 1:
            status[i] = "multi_mapped"
            continue
        c38, p0, strand, _ = res[0]
        if c38.replace("chr", "") != chroms[i]:
            status[i] = "other_chrom"
            continue
        pos38[i] = p0 + 1
        if strand == "-":
            ref38[i] = COMP[ref38[i]]
            alt38[i] = COMP[alt38[i]]
        status[i] = "ok"
    d["lift"] = status
    d["pos38"] = pos38
    d["ref38"] = ref38
    d["alt38"] = alt38
    os.makedirs(a.work, exist_ok=True)
    d.to_csv(P(a.work, "dnv_grch38.tsv"), sep="\t", index=False)
    log(f"de novo list: {len(d)} rows, {d['id'].nunique()} probands with >=1 DNV")
    print(d["lift"].value_counts().to_string())
    print(d.groupby("centre")["id"].nunique().rename("probands").to_string())
    log(f"wrote {P(a.work, 'dnv_grch38.tsv')}")


# ======================================================================== prep-am
def cmd_prep_am(a):
    path = P(a.data, AM_FILE)
    out = P(a.work, "am")
    if os.path.exists(P(out, "transcripts.pkl")) and not a.force:
        log("prep-am already done (use --force to redo)")
        return
    skip, cols = am_header(path)
    need = ["CHROM", "POS", "REF", "ALT", "transcript_id", "protein_variant", "am_pathogenicity"]
    miss = [c for c in need if c not in cols]
    if miss:
        die(f"AlphaMissense file lacks columns {miss}; header was {cols}")
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    tx_ids, tx_chrom = {}, {}
    n_in = n_auto = n_bad = 0
    part = 0
    reader = pd.read_csv(path, sep="\t", skiprows=skip + 1, header=None, names=cols, usecols=need,
                         dtype={"CHROM": str, "POS": np.int64, "REF": str, "ALT": str, "transcript_id": str,
                                "protein_variant": str, "am_pathogenicity": np.float64},
                         chunksize=a.chunksize)
    for ch in reader:
        n_in += len(ch)
        chrom = ch["CHROM"].str.replace("chr", "", regex=False)
        keep = chrom.isin(AUTOSOMES).to_numpy()
        ch, chrom = ch.loc[keep], chrom.loc[keep]
        n_auto += len(ch)
        pv = ch["protein_variant"].fillna("")
        ppos = pd.to_numeric(pv.str[1:-1], errors="coerce")
        ref = ch["REF"].map(B2I)
        alt = ch["ALT"].map(B2I)
        good = (ppos.notna() & ref.notna() & alt.notna() & pv.str.len().ge(3)
                & ch["am_pathogenicity"].notna()).to_numpy()
        n_bad += int((~good).sum())
        ch, chrom, pv, ppos, ref, alt = (x.loc[good] for x in (ch, chrom, pv, ppos, ref, alt))
        if len(ch) == 0:
            continue
        codes, uniq = pd.factorize(ch["transcript_id"])
        gid = np.array([tx_ids.setdefault(t, len(tx_ids)) for t in uniq], dtype=np.int32)
        cvals = chrom.to_numpy(dtype=object)
        _, first = np.unique(codes, return_index=True)
        for k, i in enumerate(first):
            tx_chrom.setdefault(int(gid[k]), cvals[i])
        df = pd.DataFrame({"pos": ch["POS"].to_numpy(dtype=np.int64).astype(np.int32),
                           "ref": ref.to_numpy().astype(np.int8),
                           "alt": alt.to_numpy().astype(np.int8),
                           "tx": gid[codes],
                           "prot_pos": ppos.to_numpy().astype(np.int32),
                           "ref_aa": ascii_codes(pv.str[0]),
                           "alt_aa": ascii_codes(pv.str[-1]),
                           "am": ch["am_pathogenicity"].to_numpy(dtype=np.float64)})
        for c in pd.unique(cvals):
            m = cvals == c
            df.loc[m].reset_index(drop=True).to_pickle(P(out, f"chr{c}.part{part:04d}.pkl"))
        part += 1
        log(f"  AlphaMissense rows read: {n_in:,}")
    if n_auto == 0:
        die("no autosomal rows read from AlphaMissense")
    tx = pd.DataFrame({"tx": np.arange(len(tx_ids), dtype=np.int32), "transcript_id": list(tx_ids.keys())})
    tx["am_chrom"] = tx["tx"].map(tx_chrom)
    tx.to_pickle(P(out, "transcripts.pkl"))
    log(f"AlphaMissense: {n_in:,} rows, {n_auto:,} autosomal, {n_bad:,} unparsable, {len(tx):,} transcripts")


def load_am(work, c):
    files = sorted(glob.glob(P(work, "am", f"chr{c}.part*.pkl")))
    if not files:
        die(f"no AlphaMissense rows for chr{c}; run prep-am")
    return pd.concat([pd.read_pickle(f) for f in files], ignore_index=True)


# ======================================================================== prep-gtf
def cmd_prep_gtf(a):
    am = pd.read_pickle(P(a.work, "am", "transcripts.pkl"))
    wanted = set(am["transcript_id"])
    path = P(a.data, GTF_FILE)
    n_head = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                n_head += 1
            else:
                break
    cols = ["chrom", "src", "feature", "start", "end", "score", "strand", "phase", "attr"]
    parts = []
    for ch in pd.read_csv(path, sep="\t", header=None, names=cols, skiprows=n_head, quoting=3,
                          usecols=["chrom", "feature", "start", "end", "strand", "phase", "attr"],
                          dtype={"chrom": str, "feature": str, "start": np.int64, "end": np.int64,
                                 "strand": str, "phase": str, "attr": str}, chunksize=500_000):
        ch = ch[ch["feature"] == "CDS"]
        tid = ch["attr"].str.extract(r'transcript_id "([^"]+)"', expand=False)
        keep = tid.isin(wanted)
        ch = ch.loc[keep].copy()
        ch["transcript_id"] = tid.loc[keep]
        parts.append(ch)
    cds = pd.concat(parts, ignore_index=True)
    if len(cds) == 0:
        die("no CDS lines matched AlphaMissense transcript IDs - is this the GENCODE v32 GTF?")
    attr = cds["attr"]
    cds["gene_name"] = attr.str.extract(r'gene_name "([^"]+)"', expand=False)
    cds["gene_id"] = attr.str.extract(r'gene_id "([^"]+)"', expand=False)
    cds["nf"] = attr.str.contains(r'tag "(?:cds_start_NF|cds_end_NF)"', regex=True)
    cds["chrom"] = cds["chrom"].str.replace("chr", "", regex=False)
    cds["len"] = cds["end"] - cds["start"] + 1
    cds = cds.drop(columns=["attr", "feature"]).sort_values(["transcript_id", "start"]).reset_index(drop=True)
    g = cds.groupby("transcript_id", sort=False)
    t = g.agg(chrom=("chrom", "first"), n_chrom=("chrom", "nunique"), strand=("strand", "first"),
              n_strand=("strand", "nunique"), cds_len=("len", "sum"), n_seg=("len", "size"), nf=("nf", "any"),
              gene_name=("gene_name", "first"), gene_id=("gene_id", "first"),
              phase_lo=("phase", "first"), phase_hi=("phase", "last")).reset_index()
    t["first_phase"] = np.where(t["strand"] == "+", t["phase_lo"], t["phase_hi"])
    prev_end = g["end"].shift(1)
    overlapping = set(cds.loc[cds["start"] <= prev_end, "transcript_id"])
    tx = am.merge(t, on="transcript_id", how="left")
    status = np.full(len(tx), "ok", dtype=object)
    conds = [("overlapping_segments", tx["transcript_id"].isin(overlapping)),
             ("len_not_multiple_of_3", tx["cds_len"].fillna(0) % 3 != 0),
             ("first_phase_not_0", tx["first_phase"].astype(str) != "0"),
             ("cds_incomplete_tag", tx["nf"].fillna(False).astype(bool)),
             ("chrom_or_strand_mismatch", (tx["chrom"] != tx["am_chrom"]) | (tx["n_chrom"] > 1) | (tx["n_strand"] > 1)),
             ("not_in_gtf", tx["chrom"].isna())]
    for name, cond in conds:              # later entries override earlier ones
        status[cond.fillna(False).to_numpy(dtype=bool)] = name
    tx["status"] = status
    tx = tx[["tx", "transcript_id", "am_chrom", "chrom", "strand", "cds_len", "n_seg", "gene_name", "gene_id", "status"]]
    ok = tx[tx["status"] == "ok"]
    seg = cds.merge(ok[["transcript_id", "tx"]], on="transcript_id").sort_values(["tx", "start"])
    seg["rank"] = seg.groupby("tx").cumcount()
    nseg = seg.groupby("tx")["start"].transform("size")
    minus = (seg["strand"] == "-").to_numpy()
    rank = seg["rank"].to_numpy().copy()
    rank[minus] = nseg.to_numpy()[minus] - 1 - rank[minus]
    seg["rank"] = rank
    seg = seg.sort_values(["tx", "rank"])[["tx", "chrom", "strand", "start", "end", "rank"]].reset_index(drop=True)
    os.makedirs(P(a.work, "gtf"), exist_ok=True)
    os.makedirs(P(a.work, "regions"), exist_ok=True)
    tx.to_pickle(P(a.work, "gtf", "transcripts.pkl"))
    seg.to_pickle(P(a.work, "gtf", "segments.pkl"))
    total = 0
    for c in AUTOSOMES:
        s = seg.loc[seg["chrom"] == c, ["start", "end"]].sort_values("start").to_numpy()
        merged = []
        for st, en in s:
            if merged and st <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], en)
            else:
                merged.append([int(st), int(en)])
        with open(P(a.work, "regions", f"chr{c}.regions.tsv"), "w") as f:
            for st, en in merged:
                f.write(f"{c}\t{st}\t{en}\n")
                total += en - st + 1
    log(f"AlphaMissense transcripts: {len(tx):,}")
    print(tx["status"].value_counts().to_string())
    log(f"coding bases to extract from Roulette: {total:,}; region files in {P(a.work, 'regions')}")


# ======================================================================== sites
def load_roulette(work, c):
    f = P(work, "roulette", f"chr{c}.tsv.gz")
    r = pd.read_csv(f, sep="\t", header=None, names=["pos", "ref", "alt", "filt", "pn"] + RATE_MODELS,
                    dtype={"pos": np.int64, "ref": "category", "alt": "category", "filt": "category",
                           "pn": "category", "MR": np.float64, "AR": np.float64, "MG": np.float64, "MC": np.float64},
                    na_values=["."], keep_default_na=False)
    # Roulette ReadMe: AR is MR corrected for Pol III / TF-binding hypermutability, written only where that
    # correction applies, so a missing AR means "no adjustment" (AR = MR), not "no rate"
    r["AR"] = r["AR"].fillna(r["MR"])
    r["rc"] = cat_lookup(r["ref"], B2I).astype(np.int8)
    r["ac"] = cat_lookup(r["alt"], B2I).astype(np.int8)
    filt_ok = {x: ("low" not in str(x).split(";")) for x in r["filt"].cat.categories}
    r["filt_ok"] = cat_lookup(r["filt"], filt_ok, default=True).astype(bool)
    n0 = len(r)
    r = r[(r["rc"] >= 0) & (r["ac"] >= 0) & (r["rc"] != r["ac"])]
    r = r.assign(key=r["pos"].to_numpy() * 4 + r["ac"].to_numpy().astype(np.int64))
    r = r.drop_duplicates("key")
    return r.drop(columns=["ref", "alt", "filt"]).reset_index(drop=True), n0


def build_positions(seg_c):
    """CDS positions of every transcript in transcript (5'->3') order."""
    pos_l, tx_l, off_l, near_l, minus_l = [], [], [], [], []
    for t, g in seg_c.groupby("tx", sort=False):
        strand = g["strand"].iat[0]
        starts, ends = g["start"].to_numpy(), g["end"].to_numpy()
        k = len(starts)
        arrs, nears = [], []
        for i in range(k):
            arr = np.arange(starts[i], ends[i] + 1) if strand == "+" else np.arange(ends[i], starts[i] - 1, -1)
            fl = np.zeros(len(arr), dtype=bool)
            if i > 0:
                fl[:SPLICE_WIN] = True
            if i < k - 1:
                fl[-SPLICE_WIN:] = True
            arrs.append(arr)
            nears.append(fl)
        pos = np.concatenate(arrs)
        n = len(pos)
        pos_l.append(pos)
        tx_l.append(np.full(n, t, dtype=np.int32))
        off_l.append(np.arange(n, dtype=np.int64))
        near_l.append(np.concatenate(nears))
        minus_l.append(np.full(n, strand == "-"))
    return (np.concatenate(pos_l).astype(np.int64), np.concatenate(tx_l), np.concatenate(off_l),
            np.concatenate(near_l), np.concatenate(minus_l))


def enumerate_snvs(seg_c, pos_u, ref_u):
    gpos, txp, off, near, minus = build_positions(seg_c)
    n = len(gpos)
    if n % 3:
        die("internal error: concatenated CDS length not a multiple of 3")
    k = np.minimum(np.searchsorted(pos_u, gpos), len(pos_u) - 1)
    found = pos_u[k] == gpos
    gb = np.where(found, ref_u[k], -1).astype(np.int64)
    tb = np.where(gb < 0, -1, np.where(minus, 3 - gb, gb))
    cod = tb.reshape(-1, 3)
    cvalid = (cod >= 0).all(axis=1)
    codc = np.where(cod < 0, 0, cod)
    cint = 16 * codc[:, 0] + 4 * codc[:, 1] + codc[:, 2]
    caa = CT[cint].copy()
    caa[~cvalid] = ord("X")
    ctx = txp[::3]
    cidx = off[::3] // 3
    is_last = np.r_[ctx[1:] != ctx[:-1], True]
    stop = caa == STOP
    bad_tx = np.unique(ctx[stop & ~is_last])
    cexcl = (~cvalid) | (cidx == 0) | (is_last & stop)
    cint_p = np.repeat(cint, 3)
    caa_p = np.repeat(caa, 3)
    cexcl_p = np.repeat(cexcl, 3)
    w = 4 ** (2 - (off % 3))
    tbc = np.where(tb < 0, 0, tb)
    alt_t = (tbc[:, None] + np.array([1, 2, 3])[None, :]) % 4
    alt_cint = cint_p[:, None] + (alt_t - tbc[:, None]) * w[:, None]
    alt_aa = CT[alt_cint]
    cls = np.where(alt_aa == caa_p[:, None], 0, np.where(alt_aa == STOP, 2, 1)).astype(np.int8)
    cls[cexcl_p] = 3
    galt = np.where(minus[:, None], 3 - alt_t, alt_t)
    sel = np.repeat(gb >= 0, 3)
    snv = pd.DataFrame({
        "pos": np.repeat(gpos, 3)[sel],
        "gref": np.repeat(gb, 3)[sel].astype(np.int8),
        "galt": galt.ravel()[sel].astype(np.int8),
        "tx": np.repeat(txp, 3)[sel],
        "prot_pos": (np.repeat(off // 3, 3) + 1)[sel].astype(np.int32),
        "cls": cls.ravel()[sel],
        "ref_aa": np.repeat(caa_p, 3)[sel],
        "alt_aa": alt_aa.ravel()[sel],
        "near": np.repeat(near, 3)[sel],
    })
    snv["key"] = snv["pos"].to_numpy() * 4 + snv["galt"].to_numpy().astype(np.int64)
    snv.loc[snv["tx"].isin(bad_tx), "cls"] = 4
    qc = {"cds_positions": int(n), "positions_without_reference_base": int((gb < 0).sum()),
          "codons": int(len(cint)), "codons_excluded": int(cexcl.sum())}
    return snv, set(int(x) for x in bad_tx), qc


def mutation_class(pn_series, gref, galt):
    """0 = CpG transition, 1 = other transition, 2 = transversion, 3 = context unavailable/inconsistent."""
    cats = list(pn_series.cat.categories)
    info = np.full((len(cats) + 1, 3), -1, dtype=np.int64)
    for i, s in enumerate(cats):
        s = str(s).upper()
        if len(s) == 5:
            info[i] = [B2I.get(s[1], -1), B2I.get(s[2], -1), B2I.get(s[3], -1)]
    codes = pn_series.cat.codes.to_numpy()
    rows = info[np.where(codes >= 0, codes, len(cats))]
    prv, ctr, nxt = rows[:, 0], rows[:, 1], rows[:, 2]
    gref = np.asarray(gref, dtype=np.int64)
    galt = np.asarray(galt, dtype=np.int64)
    ok = ctr == gref
    ti = np.abs(gref - galt) == 2
    cpg = ((gref == 1) & (nxt == 2)) | ((gref == 2) & (prv == 1))
    return np.where(~ok, 3, np.where(ti & cpg, 0, np.where(ti, 1, 2))).astype(np.int8), float(ok.mean()) if len(ok) else 1.0


def am_fine_bin(am):
    return np.clip(np.floor(np.round(np.asarray(am, dtype=np.float64) * 100, 6)), 0, 99)


def agg_sums(df, keys):
    g = df.groupby(keys, observed=True, sort=False)
    spec = {"n": ("key", "size")}
    for m in RATE_MODELS:
        spec[m] = (m, "sum")
    return g.agg(**spec).reset_index()


def process_chrom(a, c, tx, seg, dnv):
    t0 = time.time()
    r, n_raw = load_roulette(a.work, c)
    pos_u, first = np.unique(r["pos"].to_numpy(), return_index=True)
    ref_u = r["rc"].to_numpy()[first].astype(np.int64)
    ok_tx = tx.loc[(tx["chrom"] == c) & (tx["status"] == "ok"), "tx"].to_numpy()
    seg_c = seg[seg["tx"].isin(ok_tx)]
    snv, bad_tx, qc = enumerate_snvs(seg_c, pos_u, ref_u)
    qc.update({"chrom": c, "roulette_records": int(n_raw)})

    # ---- AlphaMissense join and per-transcript concordance
    am = load_am(a.work, c)
    am = am[am["tx"].isin(ok_tx)]
    am = am.assign(key=am["pos"].to_numpy().astype(np.int64) * 4 + am["alt"].to_numpy().astype(np.int64))
    am_core = am[am["prot_pos"] > 1]
    mis = snv[snv["cls"] == 1]
    m = mis[["key", "tx", "gref", "prot_pos", "ref_aa", "alt_aa", "near"]].merge(
        am_core[["key", "tx", "ref", "prot_pos", "ref_aa", "alt_aa", "am"]], on=["key", "tx"], how="inner",
        suffixes=("", "_am"))
    match = ((m["gref"] == m["ref"]) & (m["prot_pos"] == m["prot_pos_am"])
             & (m["ref_aa"] == m["ref_aa_am"]) & (m["alt_aa"] == m["alt_aa_am"]))
    mm = m.loc[match, ["key", "tx", "gref", "prot_pos", "near", "am"]]
    n_am_core = len(am_core)
    q = pd.DataFrame({"n_enum": mis.groupby("tx").size(), "n_am": am_core.groupby("tx").size(),
                      "n_match": mm.groupby("tx").size()}).reindex(ok_tx).fillna(0).astype(np.int64)
    q.index.name = "tx"
    q = q.reset_index()
    reason = np.full(len(q), "pass", dtype=object)
    low = (q["n_match"] < CONCORDANCE_MIN * q["n_am"]) | (q["n_match"] < CONCORDANCE_MIN * q["n_enum"])
    reason[low.to_numpy()] = "low_concordance"
    reason[(q["n_am"] == 0).to_numpy()] = "no_alphamissense_rows"
    reason[q["tx"].isin(bad_tx).to_numpy()] = "internal_stop"
    q["reason"] = reason
    q["chrom"] = c
    pass_tx = set(q.loc[q["reason"] == "pass", "tx"].tolist())
    del am, am_core, mis, m, match          # free memory: chr1 otherwise peaks near 7.5 GB

    # ---- universe: synonymous + AlphaMissense-matched missense SNVs in passing transcripts, unambiguous
    in_pass = snv["tx"].isin(pass_tx)
    allk = snv.loc[in_pass & snv["cls"].isin([0, 1, 2]), "key"]
    amb_keys = set(allk[allk.duplicated(keep=False)].tolist())
    syn = snv.loc[in_pass & (snv["cls"] == 0), ["key", "tx", "gref", "prot_pos", "near"]].assign(cls=0, am=np.nan)
    misu = mm[mm["tx"].isin(pass_tx)].assign(cls=1)
    u = pd.concat([syn, misu], ignore_index=True)
    u["galt"] = (u["key"].to_numpy() % 4).astype(np.int8)
    u["pos"] = u["key"].to_numpy() // 4
    ctx_mean = misu.groupby(["tx", "prot_pos"])["am"].mean().rename("ctx_am").reset_index()
    u = u.merge(ctx_mean, on=["tx", "prot_pos"], how="left")
    n_before = len(u)
    u = u[~u["key"].isin(amb_keys)]
    u = u.merge(r[["key", "filt_ok", "pn"] + RATE_MODELS], on="key", how="left")
    del r, syn, misu
    u["filt_ok"] = u["filt_ok"].fillna(False).astype(bool)
    if not isinstance(u["pn"].dtype, pd.CategoricalDtype):
        u["pn"] = u["pn"].astype("category")
    u["mutclass"], pn_ok = mutation_class(u["pn"], u["gref"], u["galt"])
    u["am_bin"] = np.where(u["cls"] == 1, am_fine_bin(u["am"].fillna(0)), -1).astype(np.int16)
    mr = u[PRIMARY_MODEL].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        mb = np.floor((np.log10(mr) - MR_LOG_MIN) / MR_LOG_STEP)
    u["mr_bin"] = np.where(np.isfinite(mb) & (mr > 0), np.clip(mb, 0, MR_NBINS - 1), -1).astype(np.int16)
    ctxv = u["ctx_am"].to_numpy(dtype=np.float64)
    u["ctx_bin"] = np.where(np.isfinite(ctxv), np.clip(np.floor(np.round(np.nan_to_num(ctxv) * 100, 6)), 0, 99),
                            -1).astype(np.int16)
    qc.update({"tx_ok_in_gtf": int(len(ok_tx)), "tx_pass": int(len(pass_tx)),
               "universe_missense": int((u["cls"] == 1).sum()), "universe_synonymous": int((u["cls"] == 0).sum()),
               "ambiguous_snvs_dropped": int(n_before - len(u)),
               "missing_rate_MR": int(u["MR"].isna().sum()), "filter_low": int((~u["filt_ok"]).sum()),
               "pentamer_matches_ref": round(pn_ok, 5),
               "am_rows_matched": int(len(mm)), "am_rows_core": int(n_am_core)})
    qc.update({f"rate_present_{m_}": int(u[m_].notna().sum()) for m_ in RATE_MODELS})

    # ---- aggregates of mutational opportunity (the expected-count denominators)
    musu = u[u["cls"] == 1]
    sysu = u[u["cls"] == 0]
    agg = {"mis_tx_am": agg_sums(musu, ["tx", "am_bin", "mutclass", "filt_ok", "near"]),
           "syn_tx": agg_sums(sysu, ["tx", "mutclass", "filt_ok", "near"]),
           "syn_bins": agg_sums(sysu, ["mr_bin", "ctx_bin", "mutclass", "filt_ok", "near"]),
           "mis_bins": agg_sums(musu, ["am_bin", "mutclass", "filt_ok", "near"])}
    del musu, sysu, mm

    # ---- annotate the observed de novos on this chromosome
    dv = dnv[(dnv["chrom"] == c) & (dnv["lift"] == "ok")][["dnv_idx", "pos38", "ref38", "alt38"]].copy()
    dv["gref"] = dv["ref38"].map(B2I).astype(np.int64)
    dv["key"] = dv["pos38"].to_numpy(dtype=np.int64) * 4 + dv["alt38"].map(B2I).to_numpy(dtype=np.int64)
    if len(pos_u):
        k = np.minimum(np.searchsorted(pos_u, dv["pos38"].to_numpy()), len(pos_u) - 1)
        has = pos_u[k] == dv["pos38"].to_numpy()
        rref = np.where(has, ref_u[k], -1)
    else:
        has = np.zeros(len(dv), dtype=bool)
        rref = np.full(len(dv), -1)
    st = np.full(len(dv), "not_in_cds", dtype=object)
    st[has & (rref != dv["gref"].to_numpy())] = "ref_mismatch"
    snv_first = snv.drop_duplicates("key").set_index("key")
    del snv
    hit = dv["key"].isin(snv_first.index).to_numpy() & (st != "ref_mismatch")
    cls_hit = snv_first["cls"].reindex(dv["key"]).to_numpy()
    tx_hit = snv_first["tx"].reindex(dv["key"]).to_numpy()
    st[hit & (cls_hit == 2)] = "stop_gained"
    st[hit & (cls_hit == 3)] = "excluded_codon"
    st[hit & (cls_hit == 4)] = "internal_stop_tx"
    fail_tx = hit & np.isin(cls_hit, [0, 1]) & ~np.isin(tx_hit, list(pass_tx))
    st[fail_tx] = "tx_failed_concordance"
    st[hit & np.isin(cls_hit, [0, 1]) & np.isin(tx_hit, list(pass_tx))] = "missense_without_am"
    st[dv["key"].isin(amb_keys).to_numpy() & (st != "ref_mismatch")] = "ambiguous_multi_tx"
    uu = u.set_index("key")
    del u
    inu = dv["key"].isin(uu.index).to_numpy() & (st != "ref_mismatch")
    st[inu] = np.where(uu["cls"].reindex(dv["key"][inu]).to_numpy() == 1, "missense", "synonymous")
    dv["status"] = st
    cols = ["tx", "cls", "am", "am_bin", "ctx_am", "ctx_bin", "mr_bin", "mutclass", "filt_ok", "near"] + RATE_MODELS
    ann = uu[cols].reindex(dv["key"]).reset_index(drop=True)
    dv = pd.concat([dv[["dnv_idx", "status"]].reset_index(drop=True), ann], axis=1)
    outdir = P(a.work, "sites")
    os.makedirs(outdir, exist_ok=True)
    pd.to_pickle(agg, P(outdir, f"chr{c}.agg.pkl"))
    dv.to_pickle(P(outdir, f"chr{c}.dnv.pkl"))
    q.to_pickle(P(outdir, f"chr{c}.txqc.pkl"))
    with open(P(outdir, f"chr{c}.qc.json"), "w") as f:
        json.dump(qc, f, indent=1)
    open(P(outdir, f"chr{c}.done"), "w").close()
    log(f"chr{c}: {len(pass_tx)}/{len(ok_tx)} transcripts pass, {qc['universe_missense']:,} missense + "
        f"{qc['universe_synonymous']:,} synonymous sites, {int((dv['status'] == 'missense').sum())} missense / "
        f"{int((dv['status'] == 'synonymous').sum())} synonymous de novos ({time.time() - t0:.0f}s)")


def cmd_sites(a):
    tx = pd.read_pickle(P(a.work, "gtf", "transcripts.pkl"))
    seg = pd.read_pickle(P(a.work, "gtf", "segments.pkl"))
    dnv = pd.read_csv(P(a.work, "dnv_grch38.tsv"), sep="\t", dtype={"chrom": str, "ref38": str, "alt38": str, "id": str})
    todo = parse_chroms(a.chroms)
    for c in todo:
        if os.path.exists(P(a.work, "sites", f"chr{c}.done")) and not a.force:
            log(f"chr{c}: already done")
            continue
        if not os.path.exists(P(a.work, "roulette", f"chr{c}.tsv.gz.ok")) and not a.no_ok_check:
            log(f"chr{c}: Roulette extract not ready yet - skipped (run get_roulette.sh)")
            continue
        process_chrom(a, c, tx, seg, dnv)


# ======================================================================== run
def load_sites(work, chroms):
    agg = {k: [] for k in ["mis_tx_am", "syn_tx", "syn_bins", "mis_bins"]}
    dvs, txq, qcs = [], [], []
    for c in chroms:
        d = pd.read_pickle(P(work, "sites", f"chr{c}.agg.pkl"))
        for k in agg:
            agg[k].append(d[k])
        dvs.append(pd.read_pickle(P(work, "sites", f"chr{c}.dnv.pkl")))
        txq.append(pd.read_pickle(P(work, "sites", f"chr{c}.txqc.pkl")))
        with open(P(work, "sites", f"chr{c}.qc.json")) as f:
            qcs.append(json.load(f))
    return ({k: pd.concat(v, ignore_index=True) for k, v in agg.items()}, pd.concat(dvs, ignore_index=True),
            pd.concat(txq, ignore_index=True), pd.DataFrame(qcs))


def pct(bs, q):
    with np.errstate(all="ignore"):
        return np.nanpercentile(np.where(np.isfinite(bs), bs, np.nan), q, axis=0)


def summary(arr):
    """arr: (B+1) x k; row 0 = point estimate, rows 1.. = bootstrap replicates."""
    bs = arr[1:]
    return {"est": arr[0], "lo95": pct(bs, 2.5), "hi95": pct(bs, 97.5), "lo90": pct(bs, 5), "hi90": pct(bs, 95)}


class Boot:
    def __init__(self, n_prob, n_boot, seed):
        rng = np.random.default_rng(seed)
        self.W = np.vstack([np.ones((1, n_prob)), rng.poisson(1.0, size=(n_boot, n_prob))]).astype(np.float64)
        self.n = n_prob

    def counts(self, pidx, col, ncol, vals=None):
        C = np.zeros((self.n, ncol))
        if vals is None:
            np.add.at(C, (np.asarray(pidx), np.asarray(col)), 1.0)
        else:
            np.add.at(C, np.asarray(pidx), vals)
        return self.W @ C


def evidence_label(lo90, hi90):
    for name, cut in reversed(PP3):
        if lo90 >= cut:
            return f"PP3_{name}"
    for name, cut in reversed(PP3):
        if hi90 <= 1.0 / cut:
            return f"BP4_{name}"
    return "-"


def thresholds(lo90, hi90, edges):
    out = []
    nb = len(lo90)
    for name, cut in PP3:
        thr = None
        for i in range(nb):
            if np.all(np.nan_to_num(lo90[i:], nan=-np.inf) >= cut):
                thr = edges[i]
                break
        out.append(("PP3", name, cut, thr))
    for name, cut in PP3:
        thr = None
        for i in range(nb - 1, -1, -1):
            if np.all(np.nan_to_num(hi90[: i + 1], nan=np.inf) <= 1.0 / cut):
                thr = edges[i + 1]
                break
        out.append(("BP4", name, round(1.0 / cut, 5), thr))
    return out


def equal_mass_groups(bins, mass, k):
    """Merge ordered fine bins into ~k groups of equal mutational mass; returns {bin: group}."""
    order = np.argsort(bins)
    b, mss = np.asarray(bins)[order], np.asarray(mass, dtype=float)[order]
    cum = np.cumsum(mss) - mss / 2
    grp = np.minimum((cum / mss.sum() * k).astype(int), k - 1)
    grp = np.unique(grp, return_inverse=True)[1].ravel()   # renumber densely: a heavy fine bin can skip a group index
    return dict(zip(b.tolist(), grp.tolist()))


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{x:.{nd}f}"


def cmd_run(a):
    t0 = time.time()
    os.makedirs(a.results, exist_ok=True)
    done = [c for c in AUTOSOMES if os.path.exists(P(a.work, "sites", f"chr{c}.done"))]
    if len(done) < 22 and not a.allow_partial:
        die(f"only {len(done)}/22 chromosomes processed ({','.join(done)}); finish 'sites' or use --allow-partial")
    if not done:
        die("no chromosomes processed yet")
    tx = pd.read_pickle(P(a.work, "gtf", "transcripts.pkl"))
    agg, dv, txq, qcc = load_sites(a.work, done)
    meta = pd.read_csv(P(a.work, "dnv_grch38.tsv"), sep="\t", dtype={"chrom": str, "id": str, "cq": str})
    rep = []

    def out(line=""):
        rep.append(line)
        print(line)

    out("=" * 100)
    out("NDD de novo PP3/BP4 calibration - feasibility gate (pre-registered; see GATE_PLAN.md)")
    out(f"run: {time.strftime('%Y-%m-%d %H:%M')}   chromosomes: {','.join(done)}   rate model (primary): "
        f"{PRIMARY_MODEL}   bootstrap: {a.n_boot} Poisson replicates over probands, seed {SEED}")
    out("=" * 100)

    # ------------------------------------------------------------ [0] pipeline QC
    out("\n[0] PIPELINE QC")
    auto_meta = meta[meta["chrom"].isin(done)]
    lift = auto_meta["lift"].value_counts()
    snv_auto = auto_meta[auto_meta["lift"] != "not_snv"]
    n_ok = int((snv_auto["lift"] == "ok").sum())
    other = {k: int(v) for k, v in lift.drop("ok", errors="ignore").items()}
    out(f"  autosomal SNV de novos: {len(snv_auto):,}; lifted to GRCh38: {n_ok:,} "
        f"({n_ok / max(len(snv_auto), 1):.2%}); other: {other}")
    am_tx_auto = tx[tx["am_chrom"].isin(done)]
    n_tx_pass = int((txq["reason"] == "pass").sum())
    tx_rate = n_tx_pass / max(len(am_tx_auto), 1)
    out(f"  AlphaMissense transcripts on these chromosomes: {len(am_tx_auto):,}; pass all checks: {n_tx_pass:,} ({tx_rate:.2%})")
    reasons = pd.concat([am_tx_auto.loc[am_tx_auto["status"] != "ok", "status"], txq.loc[txq["reason"] != "pass", "reason"]])
    if len(reasons):
        out("    excluded: " + ", ".join(f"{k} {v}" for k, v in reasons.value_counts().items()))
    out(f"  coding universe: {int(qcc['universe_missense'].sum()):,} missense + {int(qcc['universe_synonymous'].sum()):,} "
        f"synonymous possible SNVs; ambiguous (>1 transcript) dropped {int(qcc['ambiguous_snvs_dropped'].sum()):,}; "
        f"positions lacking a Roulette record {int(qcc['positions_without_reference_base'].sum()):,}")
    out(f"  Roulette FILTER=low among universe sites: {int(qcc['filter_low'].sum()):,}; missing MR: "
        f"{int(qcc['missing_rate_MR'].sum()):,}; pentamer centre matches REF: {qcc['pentamer_matches_ref'].min():.4f} (min over chromosomes)")
    dmeta = dv.merge(meta[["dnv_idx", "id", "centre", "study", "cq", "symbol"]], on="dnv_idx", how="left")
    out("  lifted de novos by status: " + ", ".join(f"{k} {v}" for k, v in dmeta["status"].value_counts().items()))
    vep = dmeta[dmeta["cq"].isin(["missense_variant", "synonymous_variant"])]
    # G0 mapping rate: VEP missense/synonymous autosomal SNV de novos (denominator includes those that failed
    # liftover) that land on a universe site, counted before the Roulette FILTER/rate exclusion (a data exclusion,
    # not a mapping failure; the post-exclusion rate is reported alongside)
    vep_all = snv_auto[snv_auto["cq"].isin(["missense_variant", "synonymous_variant"])]
    in_u = vep["status"].isin(["missense", "synonymous"])
    map_rate = in_u.sum() / len(vep_all) if len(vep_all) else 0.0
    in_u_kept = in_u & vep["filt_ok"].fillna(False).astype(bool) & (vep[PRIMARY_MODEL].fillna(0) > 0)
    map_rate_kept = in_u_kept.sum() / len(vep_all) if len(vep_all) else 0.0
    # REF mismatch: among lifted de novos whose GRCh38 position has a Roulette record, i.e. whose REF can be checked
    checkable = dmeta["status"] != "not_in_cds"
    ref_mm = (dmeta["status"] == "ref_mismatch").sum() / checkable.sum() if checkable.sum() else 0.0
    xt = pd.crosstab(vep["cq"], vep["status"])
    out("  source VEP consequence x our class (canonical AlphaMissense transcript):")
    for line in xt.to_string().splitlines():
        out("    " + line)
    out(f"  VEP missense/synonymous SNV de novos: {len(vep_all):,}; on a universe site {int(in_u.sum()):,} "
        f"({map_rate:.2%}); of these at a FILTER-passing site with an MR rate {int(in_u_kept.sum()):,} "
        f"({map_rate_kept:.2%} of all, informational); REF checkable after liftover {int(checkable.sum()):,}")
    g0 = (tx_rate >= QC_TX_PASS_MIN) and (map_rate >= QC_DNV_MAP_MIN) and (ref_mm <= QC_REF_MISMATCH_MAX)
    out(f"  G0 pipeline QC: transcripts pass {tx_rate:.2%} (>= {QC_TX_PASS_MIN:.0%}), VEP missense/synonymous "
        f"de novos mapped {map_rate:.2%} (>= {QC_DNV_MAP_MIN:.0%}), REF mismatch {ref_mm:.3%} "
        f"(<= {QC_REF_MISMATCH_MAX:.0%})  ->  {'PASS' if g0 else 'FAIL'}")
    qc_rows = [("transcripts_pass_rate", tx_rate), ("dnv_map_rate", map_rate),
               ("dnv_map_rate_after_filter_informational", map_rate_kept), ("ref_mismatch_rate", ref_mm),
               ("G0", "PASS" if g0 else "FAIL")]
    pd.DataFrame(qc_rows, columns=["metric", "value"]).to_csv(P(a.results, "qc.tsv"), sep="\t", index=False)
    qcc.to_csv(P(a.results, "qc_by_chrom.tsv"), sep="\t", index=False)
    if a.qc_only:
        out("\n--qc-only: stopping before any calibration result is computed.")
        with open(P(a.results, "gate_report.txt"), "w") as f:
            f.write("\n".join(rep) + "\n")
        return

    # ------------------------------------------------------------ gene attributes
    genes = pd.read_csv(P(a.data, GENE_FILE), sep="\t", low_memory=False, usecols=["symbol", "s_het"])
    shet = genes.dropna(subset=["s_het"]).drop_duplicates("symbol").set_index("symbol")["s_het"]
    dd = pd.read_csv(P(a.data, DDG2P_FILE), dtype=str)
    qual = (dd["DDD category"].str.lower().isin(["confirmed", "probable"])
            & dd["allelic requirement"].str.lower().str.contains("monoallelic", na=False))
    dd_syms = set(dd.loc[qual, "gene symbol"].str.upper())
    if "prev symbols" in dd.columns:
        for s in dd.loc[qual, "prev symbols"].dropna():
            dd_syms.update(x.strip().upper() for x in s.split(";") if x.strip())
    utx = sorted(set(agg["mis_tx_am"]["tx"]) | set(agg["syn_tx"]["tx"]))
    tinfo = tx.set_index("tx").loc[utx, ["gene_name"]].copy()
    tinfo["s_het"] = tinfo["gene_name"].map(shet)
    tinfo["shet_class"] = np.where(tinfo["s_het"].isna(), 2, (tinfo["s_het"] >= SHET_HIGH).astype(int))
    tinfo["ddg2p"] = tinfo["gene_name"].str.upper().isin(dd_syms)
    row_of = pd.Series(np.arange(len(utx)), index=utx)
    out(f"\n  genes in universe {len(utx):,}; with s_het {int(tinfo['s_het'].notna().sum()):,} "
        f"(high >= {SHET_HIGH}: {int((tinfo['shet_class'] == 1).sum()):,}); DDG2P monoallelic confirmed/probable: "
        f"{int(tinfo['ddg2p'].sum()):,}")

    # ------------------------------------------------------------ probands, bootstrap
    meta["pkey"] = meta["centre"] + "|" + meta["id"].astype(str)
    pk = pd.Series(np.arange(meta["pkey"].nunique()), index=pd.unique(meta["pkey"]))
    dmeta["p"] = pk.reindex(dmeta["centre"] + "|" + dmeta["id"].astype(str)).to_numpy()
    boot = Boot(len(pk), a.n_boot, SEED)
    nb = len(LR_EDGES_PCT) - 1
    fine_to_coarse = np.searchsorted(np.array(LR_EDGES_PCT), np.arange(100), side="right") - 1
    edges = [e / 100 for e in LR_EDGES_PCT]

    NCLS = 4   # mutation classes: 0 CpG transition, 1 other transition, 2 transversion, 3 context unavailable
    delta = [DELTA_MIN]   # systematic margin, finalised after the null checks

    def variant_data(model, filt="ok", near="all"):
        def f(df):
            m = np.ones(len(df), dtype=bool)
            if filt == "ok":
                m &= df["filt_ok"].astype(bool).to_numpy()
            if near == "exclude":
                m &= ~df["near"].astype(bool).to_numpy()
            return m
        mt = agg["mis_tx_am"][f(agg["mis_tx_am"])]
        M = np.zeros((len(utx), 100, NCLS))
        np.add.at(M, (row_of.reindex(mt["tx"]).to_numpy(), mt["am_bin"].to_numpy().astype(int),
                      mt["mutclass"].to_numpy().astype(int)), mt[model].to_numpy())
        st = agg["syn_tx"][f(agg["syn_tx"])]
        Msyn = np.zeros((len(utx), NCLS))
        np.add.at(Msyn, (row_of.reindex(st["tx"]).to_numpy(), st["mutclass"].to_numpy().astype(int)), st[model].to_numpy())
        d = dmeta[dmeta["status"].isin(["missense", "synonymous"])].copy()
        d = d[f(d) & (d[model].to_numpy(dtype=float) > 0)].copy()
        d["row"] = row_of.reindex(d["tx"].astype(int)).to_numpy().astype(int)
        d["mc"] = d["mutclass"].to_numpy().astype(int)
        return M, Msyn, d

    def syn_counts(d, centre=None):
        """(B+1) x NCLS synonymous de novo counts by mutation class."""
        s_ = d[d["cls"] == 0]
        if centre:
            s_ = s_[s_["centre"] == centre]
        return boot.counts(s_["p"], s_["mc"], NCLS)

    def scaling(Os_k, Msyn_k, calib):
        """Scale factors, (B+1) x NCLS. pooled: one c for all classes; class: c per class (unknown class pooled)."""
        pooled = Os_k.sum(1) / Msyn_k.sum()
        if calib == "pooled":
            return np.repeat(pooled[:, None], NCLS, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            c = Os_k / Msyn_k[None, :]
        c[:, 3] = pooled
        return np.where(np.isfinite(c), c, pooled[:, None])

    def lr_curve(M, Msyn, d, rows_mask, centre=None, binmap=fine_to_coarse, nbins=nb, calib="pooled"):
        Mc = np.zeros((M.shape[0], nbins, NCLS))
        for fb in range(100):
            Mc[:, binmap[fb], :] += M[:, fb, :]
        Mc = Mc[rows_mask]
        c = scaling(syn_counts(d, centre), Msyn.sum(0), calib)
        E = c @ Mc.sum(0).T
        Eg = np.einsum("gbk,k->gb", Mc, c[0])          # expected per gene and bin at the point estimate of c
        Egs = Eg.sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            pi = np.where(Egs[:, None] > 0, Eg / Egs[:, None], 0.0)
        local = np.full(M.shape[0], -1)
        local[np.flatnonzero(rows_mask)] = np.arange(rows_mask.sum())
        mis = d[(d["cls"] == 1)]
        if centre:
            mis = mis[mis["centre"] == centre]
        mis = mis[local[mis["row"].to_numpy()] >= 0]
        bins = binmap[mis["am_bin"].to_numpy().astype(int)]
        O = boot.counts(mis["p"], bins, nbins)
        Operm = boot.counts(mis["p"].to_numpy(), None, nbins, vals=pi[local[mis["row"].to_numpy()]])
        X = O - E
        Xup = O - (1 - delta[0]) * E      # background over-estimated by delta: more excess (benign-direction bound)
        Xdn = O - (1 + delta[0]) * E      # background under-estimated by delta: less excess (pathogenic-direction bound)
        with np.errstate(invalid="ignore", divide="ignore"):
            Xt = X.sum(1, keepdims=True)
            fB = E / E.sum(1, keepdims=True)
            den = Operm - E
            res = {"O": O, "E": E, "X": X, "c": c,
                   "LRm": (X / Xt) / fB, "LRm_up": (Xup / Xt) / fB, "LRm_dn": (Xdn / Xt) / fB,
                   "LRw": X / den, "LRw_up": Xup / den, "LRw_dn": Xdn / den}
        return res

    # ------------------------------------------------------------ [1] null calibration
    M, Msyn, d = variant_data(PRIMARY_MODEL)
    Os_k = syn_counts(d)
    out("\n[1] NULL CALIBRATION - synonymous de novos; expected = c x sum(Roulette MR); c fitted on all synonymous de novos")
    for cen in CENTRES:
        oc = syn_counts(d, cen)[0].sum()
        out(f"  {cen:7s} synonymous de novos in universe: {int(oc):,}  (scale c = {oc / Msyn.sum():.4g})")
    sd = d[d["cls"] == 0]
    sb = agg["syn_bins"]
    sb = sb[sb["filt_ok"].astype(bool)]
    null_rows = []
    from scipy.stats import chi2, norm

    def null_table(name, group_of_dnv, Mk, labels, c, gate, verdicts, tag, ordered=False):
        """Mk: K x NCLS synonymous mutational mass; c: (B+1) x NCLS scale factors."""
        K = Mk.shape[0]
        keep = group_of_dnv >= 0
        Ok = boot.counts(sd["p"].to_numpy()[keep], group_of_dnv[keep], K)
        Ek = c @ Mk.T
        with np.errstate(invalid="ignore", divide="ignore"):
            s = summary(Ok / Ek)
        o0 = Ok[0]
        e0 = o0.sum() * Ek[0] / Ek[0].sum()
        okk = e0 > 0
        stat = float((((o0 - e0) ** 2)[okk] / e0[okk]).sum())
        pval = float(chi2.sf(stat, okk.sum() - 1)) if okk.sum() > 1 else float("nan")
        dev = float(np.nanmax(np.abs(s["est"] - 1))) if K else float("nan")
        # FAIL: a group is significantly miscalibrated by more than NULL_PASS (95% CI wholly outside 1 +- NULL_PASS)
        # WARN: a point estimate outside 1 +- NULL_PASS, omnibus chi2 p < 0.001, or (ordered groups) trend p < 0.001
        sig_off = bool(np.any((s["lo95"] > 1 + NULL_PASS) | (s["hi95"] < 1 - NULL_PASS)))
        trend_txt, trend_p = "", float("nan")
        use = okk & (o0 > 0)
        if ordered and use.sum() > 2:
            x = np.linspace(0, 1, K)[use]
            y = np.log(o0[use] / Ek[0][use])
            w = o0[use]
            xb = (w * x).sum() / w.sum()
            sxx = (w * (x - xb) ** 2).sum()
            beta = (w * (x - xb) * (y - (w * y).sum() / w.sum())).sum() / sxx
            trend_p = float(2 * norm.sf(abs(beta) * np.sqrt(sxx)))
            trend_txt = f"   trend: last/first group ratio {np.exp(beta):.3f}, p = {trend_p:.3g}"
        if sig_off:
            verdict = "FAIL"
        elif dev > NULL_PASS or pval < 0.001 or trend_p < 0.001:
            verdict = "WARN"
        else:
            verdict = "PASS"
        if gate:
            verdicts[name] = verdict
        out(f"\n  by {name}{tag}:  max |O/E-1| = {dev:.3f}   chi2 p = {pval:.3g}{trend_txt}   ->  "
            f"{verdict}{'' if gate else ' (informational)'}")
        out(f"    {'group':34s} {'O':>7s} {'E':>9s} {'O/E':>7s}  95% CI")
        for k in range(K):
            out(f"    {labels[k]:34s} {o0[k]:7.0f} {Ek[0, k]:9.1f} {s['est'][k]:7.3f}  "
                f"[{fmt(s['lo95'][k])}, {fmt(s['hi95'][k])}]")
            null_rows.append([name + tag, labels[k], o0[k], Ek[0, k], s["est"][k], s["lo95"][k], s["hi95"][k], pval, verdict])

    def grouped_mass(col, kgroups):
        m = sb.groupby([col, "mutclass"])[PRIMARY_MODEL].sum().reset_index()
        m = m[m[col] >= 0]
        tot = m.groupby(col)[PRIMARY_MODEL].sum()
        gmap = equal_mass_groups(tot.index.to_numpy(), tot.to_numpy(), kgroups)
        K = max(gmap.values()) + 1
        Mk = np.zeros((K, NCLS))
        np.add.at(Mk, (m[col].map(gmap).to_numpy().astype(int), m["mutclass"].to_numpy().astype(int)), m[PRIMARY_MODEL].to_numpy())
        members = {k: [b for b, g in gmap.items() if g == k] for k in range(K)}
        return gmap, Mk, members

    mr_gmap, mr_Mk, mr_mem = grouped_mass("mr_bin", N_NULL_GROUPS)
    mr_lab = [f"MR {10 ** (MR_LOG_MIN + min(b) * MR_LOG_STEP):.3g}-{10 ** (MR_LOG_MIN + (max(b) + 1) * MR_LOG_STEP):.3g}"
              for b in mr_mem.values()]
    mr_grp = np.array([mr_gmap.get(int(b), -1) for b in sd["mr_bin"].to_numpy()])
    cx_gmap, cx_Mk, cx_mem = grouped_mass("ctx_bin", N_NULL_GROUPS)
    cx_lab = [f"codon mean AM {min(b) / 100:.2f}-{(max(b) + 1) / 100:.2f}" for b in cx_mem.values()]
    cx_grp = np.array([cx_gmap.get(int(b), -1) for b in sd["ctx_bin"].to_numpy()])
    mc_tot = sb.groupby("mutclass")[PRIMARY_MODEL].sum().reindex(range(NCLS)).fillna(0).to_numpy()
    mc_Mk = np.zeros((3, NCLS))
    for j in range(3):
        mc_Mk[j, j] = mc_tot[j]
    mc_grp = np.where(sd["mc"].to_numpy() <= 2, sd["mc"].to_numpy(), -1)

    def run_nulls(calib, tag=""):
        c = scaling(Os_k, Msyn.sum(0), calib)
        v = {}
        null_table("Roulette mutability decile", mr_grp, mr_Mk, mr_lab, c, True, v, tag, ordered=True)
        null_table("AlphaMissense context of the codon", cx_grp, cx_Mk, cx_lab, c, True, v, tag, ordered=True)
        if calib == "pooled":
            null_table("mutation class", mc_grp, mc_Mk, ["CpG transition", "non-CpG transition", "transversion"], c, True, v, tag)
        return v, c

    verdicts, c_pool = run_nulls("pooled")
    shc = tinfo["shet_class"].to_numpy()
    null_table("gene s_het class", shc[sd["row"].to_numpy()], np.array([Msyn[shc == k].sum(0) for k in (0, 1, 2)]),
               ["s_het < 0.15", "s_het >= 0.15", "s_het unknown"], c_pool, False, {}, "")
    ddv = tinfo["ddg2p"].to_numpy().astype(int)
    null_table("DDG2P gene", ddv[sd["row"].to_numpy()], np.array([Msyn[ddv == k].sum(0) for k in (0, 1)]),
               ["not DDG2P", "DDG2P monoallelic"], c_pool, False, {}, "")
    st_all = agg["syn_tx"][agg["syn_tx"]["filt_ok"].astype(bool)]
    near_Mk = np.zeros((2, NCLS))
    np.add.at(near_Mk, (st_all["near"].astype(bool).to_numpy().astype(int), st_all["mutclass"].to_numpy().astype(int)),
              st_all[PRIMARY_MODEL].to_numpy())
    null_table("exon-edge position", sd["near"].astype(bool).to_numpy().astype(int), near_Mk,
               ["interior", f"within {SPLICE_WIN} bp of a splice junction"], c_pool, False, {}, "")
    g1 = "FAIL" if "FAIL" in verdicts.values() else ("WARN" if "WARN" in verdicts.values() else "PASS")
    calib = "pooled"
    if g1 == "FAIL" and verdicts.get("mutation class") == "FAIL":
        out("\n  The mutation-class check failed with one pooled scale factor (a class miscalibration also distorts the "
            "mutability groups). Pre-registered fallback: a separate scale factor per mutation class, fitted on synonymous "
            "de novos, then mutability and protein context are re-checked:")
        v2, _ = run_nulls("class", " [class-calibrated]")
        if all(v != "FAIL" for v in v2.values()):
            calib, g1 = "class", "WARN"
            out("  fallback accepted: all downstream expectations use class-specific scaling (G1 -> WARN)")
        else:
            out("  fallback rejected: miscalibration persists after class-specific scaling")
    pd.DataFrame(null_rows, columns=["grouping", "group", "O", "E", "OE", "OE_lo95", "OE_hi95", "chi2_p", "verdict"]).to_csv(
        P(a.results, "null_checks.tsv"), sep="\t", index=False)
    out(f"\n  G1 null calibration (mutability, protein context, mutation class): {g1}   [expectation: {calib} scaling]")
    tagsel = "" if calib == "pooled" else " [class-calibrated]"
    gnames = {"Roulette mutability decile", "AlphaMissense context of the codon"} | ({"mutation class"} if calib == "pooled" else set())
    devs = [abs(r[4] - 1) for r in null_rows if r[0] in {g + tagsel for g in gnames} and np.isfinite(r[4])]
    delta[0] = max([DELTA_MIN] + devs)
    out(f"  systematic margin for evidence bounds: delta = {delta[0]:.3f} (largest |O/E - 1| above, floor {DELTA_MIN})")

    # ------------------------------------------------------------ [2] missense excess
    allrows = np.ones(len(utx), dtype=bool)
    base = lr_curve(M, Msyn, d, allrows, calib=calib)
    O_t = base["O"].sum(1)
    E_t = base["E"].sum(1)
    s_oe = summary((O_t / E_t)[:, None])
    s_x = summary((O_t - E_t)[:, None])
    out("\n[2] MISSENSE EXCESS (all genes)")
    out(f"  observed {O_t[0]:.0f}, expected {E_t[0]:.1f}, O/E {s_oe['est'][0]:.3f} [{s_oe['lo95'][0]:.3f}, "
        f"{s_oe['hi95'][0]:.3f}], excess {s_x['est'][0]:.0f} [{s_x['lo95'][0]:.0f}, {s_x['hi95'][0]:.0f}] "
        f"= pathogenic fraction {1 - E_t[0] / O_t[0]:.3f}")
    nc_res = lr_curve(M, Msyn, d, shc == 0, binmap=(np.arange(100) >= 10).astype(int), nbins=2, calib=calib)
    with np.errstate(invalid="ignore", divide="ignore"):
        nc_s = summary(nc_res["O"] / nc_res["E"])
    nc_flag = bool(nc_s["lo95"][0] > 1 + NULL_PASS or nc_s["hi95"][0] < 1 - NULL_PASS)
    out(f"  negative control - missense with AlphaMissense < 0.10 in genes with s_het < {SHET_HIGH}: "
        f"O {nc_res['O'][0, 0]:.0f}, E {nc_res['E'][0, 0]:.1f}, O/E {nc_s['est'][0]:.3f} "
        f"[{fmt(nc_s['lo95'][0])}, {fmt(nc_s['hi95'][0])}]  ->  {'WARN (synonymous scaling may not transfer to missense)' if nc_flag else 'PASS'}")
    if nc_flag and g1 == "PASS":
        g1 = "WARN"

    # ------------------------------------------------------------ [3] LR curves
    lr_rows, thr_rows, curves = [], [], {}

    def record(analysis, kind, res, key):
        s = summary(res[key])
        lo_sys = summary(res[key + "_dn"])["lo90"]     # one-sided 95% bound with E inflated by delta
        hi_sys = summary(res[key + "_up"])["hi90"]     # one-sided 95% bound with E deflated by delta
        curves[(analysis, kind)] = s
        for i in range(nb):
            lr_rows.append([analysis, kind, edges[i], edges[i + 1], res["O"][0, i], res["E"][0, i], res["X"][0, i],
                            s["est"][i], s["lo95"][i], s["hi95"][i], s["lo90"][i], s["hi90"][i], lo_sys[i], hi_sys[i],
                            evidence_label(lo_sys[i], hi_sys[i])])
        for lvl in thresholds(lo_sys, hi_sys, edges):
            thr_rows.append([analysis, kind] + list(lvl))
        return s

    no_excess = set()

    def print_curve(title, res, s_m, s_w=None, name=None):
        out(f"\n  {title}")
        xt = summary(res["X"].sum(1)[:, None])
        if not (xt["lo95"][0] > 0):
            out(f"    no detectable missense excess in this stratum (excess {xt['est'][0]:.0f} "
                f"[{xt['lo95'][0]:.0f}, {xt['hi95'][0]:.0f}]); likelihood ratios are undefined here")
            if name:
                no_excess.add(name)
            return
        hdr = f"    {'AM bin':12s} {'O':>6s} {'E':>8s} {'excess':>7s}   {'LR marginal [95% CI]':28s}"
        if s_w is not None:
            hdr += f" {'LR within-gene [95% CI]':28s}"
        out(hdr)
        for i in range(nb):
            line = (f"    {edges[i]:.2f}-{edges[i + 1]:.2f}    {res['O'][0, i]:6.0f} {res['E'][0, i]:8.1f} "
                    f"{res['X'][0, i]:7.1f}   {fmt(s_m['est'][i], 2):>6s} [{fmt(s_m['lo95'][i], 2)}, {fmt(s_m['hi95'][i], 2)}]"
                    f"{'':8s}")
            if s_w is not None:
                line += f" {fmt(s_w['est'][i], 2):>6s} [{fmt(s_w['lo95'][i], 2)}, {fmt(s_w['hi95'][i], 2)}]"
            out(line)

    out("\n[3] LIKELIHOOD RATIOS BY ALPHAMISSENSE BIN (de novo missense; pathogenic = excess, benign = background)")
    out("    statistical 95% CIs; a negative value means fewer de novos than expected (read as ~0)")
    s_m = record("all_genes", "marginal", base, "LRm")
    s_w = record("all_genes", "within_gene", base, "LRw")
    print_curve("all genes (marginal = gene + variant information; within-gene = variant information only)", base, s_m, s_w,
                "all_genes")
    ddrows = tinfo["ddg2p"].to_numpy()
    res_dd = lr_curve(M, Msyn, d, ddrows, calib=calib)
    s1 = record("ddg2p_monoallelic", "marginal", res_dd, "LRm")
    s2 = record("ddg2p_monoallelic", "within_gene", res_dd, "LRw")
    print_curve("DDG2P monoallelic confirmed/probable genes (the clinical use case: known NDD gene)", res_dd, s1, s2,
                "ddg2p_monoallelic")
    for k, nm in ((1, "shet_high"), (0, "shet_low")):
        res_s = lr_curve(M, Msyn, d, tinfo["shet_class"].to_numpy() == k, calib=calib)
        a1 = record(nm, "marginal", res_s, "LRm")
        a2 = record(nm, "within_gene", res_s, "LRw")
        print_curve(f"genes with s_het {'>=' if k else '<'} {SHET_HIGH}", res_s, a1, a2, nm)

    # cumulative LR(score >= t)
    out("\n  cumulative: LR for 'AlphaMissense >= t' (all genes)")
    cum_rows = []
    for tp in CUM_THRESH_PCT:
        bm = (np.arange(100) >= tp).astype(int)
        rc = lr_curve(M, Msyn, d, allrows, binmap=bm, nbins=2, calib=calib)
        smm, sww = summary(rc["LRm"]), summary(rc["LRw"])
        out(f"    AM >= {tp / 100:.2f}: marginal {fmt(smm['est'][1], 2)} [{fmt(smm['lo95'][1], 2)}, {fmt(smm['hi95'][1], 2)}]"
            f"   within-gene {fmt(sww['est'][1], 2)} [{fmt(sww['lo95'][1], 2)}, {fmt(sww['hi95'][1], 2)}]")
        cum_rows.append([tp / 100, smm["est"][1], smm["lo95"][1], smm["hi95"][1], sww["est"][1], sww["lo95"][1], sww["hi95"][1]])
    pd.DataFrame(cum_rows, columns=["am_at_least", "LR_marginal", "lo95", "hi95", "LR_within", "lo95_w", "hi95_w"]).to_csv(
        P(a.results, "cumulative_lr.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------ [4] ClinGen thresholds
    out("\n[4] CLINGEN THRESHOLDS (Pejaver et al. 2022 rule: smallest score whose bins all have one-sided 95% bound past the LR cut-off)")
    out(f"    bounds include a +-{delta[0]:.3f} systematic margin on expected counts; bin resolution only (the full study uses "
        "sliding local windows)")
    out("    BP4 levels rest on the accuracy of the expected counts in low-score bins - read them with the null checks")
    thr = pd.DataFrame(thr_rows, columns=["analysis", "kind", "direction", "level", "LR_cutoff", "AM_threshold"])
    for (an, kind), g in thr.groupby(["analysis", "kind"], sort=False):
        if an.startswith("centre_") or an.startswith("sens:"):
            continue
        if an in no_excess:
            out(f"  {an:18s} {kind:12s} (no detectable excess - not estimable)")
            continue
        parts = [f"{r.direction}_{r.level} " + ("not reached" if pd.isna(r.AM_threshold) else
                 f"{'>=' if r.direction == 'PP3' else '<='} {r.AM_threshold:.2f}") for r in g.itertuples()]
        out(f"  {an:18s} {kind:12s} " + "; ".join(parts))
    thr.to_csv(P(a.results, "thresholds.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------ [5] centre consistency
    out("\n[5] CENTRE CONSISTENCY (marginal LR, all genes; centre-specific synonymous scaling)")
    from scipy.stats import chi2
    per = {cen: lr_curve(M, Msyn, d, allrows, centre=cen, calib=calib) for cen in CENTRES}
    for cen in CENTRES:
        record(f"centre_{cen}", "marginal", per[cen], "LRm")
        record(f"centre_{cen}", "within_gene", per[cen], "LRw")
    het_p = []
    out(f"    {'AM bin':12s}" + "".join(f" {cen + ' LR [95% CI]':26s}" for cen in CENTRES) + "  Q p-value")
    for i in range(nb):
        ests, vars_ = [], []
        line = f"    {edges[i]:.2f}-{edges[i + 1]:.2f}  "
        for cen in CENTRES:
            arr = per[cen]["LRm"][:, i]
            bs = arr[1:][np.isfinite(arr[1:])]
            ests.append(arr[0])
            vars_.append(np.var(bs) if len(bs) > 10 else np.nan)
            ss = summary(per[cen]["LRm"][:, [i]])
            line += f" {fmt(arr[0], 2):>6s} [{fmt(ss['lo95'][0], 2)}, {fmt(ss['hi95'][0], 2)}]{'':6s}"
        ests, vars_ = np.array(ests), np.array(vars_)
        if np.all(np.isfinite(ests)) and np.all(vars_ > 0):
            wv = 1 / vars_
            mu = (wv * ests).sum() / wv.sum()
            pq = float(chi2.sf(((ests - mu) ** 2 * wv).sum(), len(ests) - 1))
        else:
            pq = float("nan")
        het_p.append(pq)
        out(line + f"  {pq:.3g}")
    top_het = np.array(het_p[-TOP_BINS:])
    g3 = "PASS" if np.all(np.isfinite(top_het) & (np.nan_to_num(top_het, nan=0.0) >= 0.01)) else "WARN"
    out(f"  G3 centre consistency (Q p >= 0.01 in the top {TOP_BINS} bins; a p-value that cannot be computed counts as WARN): {g3}")

    # ------------------------------------------------------------ [6] sensitivity analyses
    out("\n[6] SENSITIVITY (marginal / within-gene LR in the top bins; all genes)")
    other = "class" if calib == "pooled" else "pooled"
    sens = [("rate model AR (adjusted Roulette)", "AR", "ok", "all", calib), ("rate model MG (gnomAD)", "MG", "ok", "all", calib),
            ("rate model MC (Carlson)", "MC", "ok", "all", calib),
            ("exclude sites near splice junctions", PRIMARY_MODEL, "ok", "exclude", calib),
            ("include Roulette FILTER=low sites", PRIMARY_MODEL, "all", "all", calib),
            (f"{other}-calibrated expectation", PRIMARY_MODEL, "ok", "all", other)]
    out(f"    {'analysis':40s}" + "".join(f"  {edges[i]:.2f}-{edges[i + 1]:.2f} (marg / within)" for i in range(nb - TOP_BINS, nb)))
    out(f"    {'primary (MR, ' + calib + ' scaling)':40s}" + "".join(f"  {fmt(s_m['est'][i], 2):>6s} / {fmt(s_w['est'][i], 2):<6s}{'':12s}" for i in range(nb - TOP_BINS, nb)))
    for label, model, filt, near, cal in sens:
        Ms, Mss, ds = variant_data(model, filt, near)
        if ds.empty or Mss.sum() <= 0:
            out(f"    {label:40s}  (no valid rates)")
            continue
        rs = lr_curve(Ms, Mss, ds, allrows, calib=cal)
        a1 = record("sens:" + label, "marginal", rs, "LRm")
        a2 = record("sens:" + label, "within_gene", rs, "LRw")
        cov = ""
        if f"rate_present_{model}" in qcc.columns:
            n_u = qcc["universe_missense"].sum() + qcc["universe_synonymous"].sum()
            cov = f"  [{model} present at {qcc[f'rate_present_{model}'].sum() / max(n_u, 1):.1%} of universe sites]"
        out(f"    {label:40s}" + "".join(f"  {fmt(a1['est'][i], 2):>6s} / {fmt(a2['est'][i], 2):<6s}{'':12s}" for i in range(nb - TOP_BINS, nb)) + cov)

    pd.DataFrame(lr_rows, columns=["analysis", "kind", "bin_lo", "bin_hi", "O", "E", "excess", "LR", "lo95", "hi95",
                                   "lo_one_sided95", "hi_one_sided95", "lo_one_sided95_sys", "hi_one_sided95_sys",
                                   "evidence"]).to_csv(
        P(a.results, "lr_bins.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------ [7] verdict
    top = slice(nb - TOP_BINS, nb)
    ratio = s_m["hi95"][top] / s_m["lo95"][top]
    g2 = bool(np.all(np.isfinite(ratio)) and np.all(s_m["lo95"][top] > 0) and np.all(ratio <= PRECISION_MAX_RATIO)
              and s_m["lo95"][nb - 1] > 1)
    out("\n[7] GATE VERDICT")
    out(f"  G0 pipeline QC ................................ {'PASS' if g0 else 'FAIL'}")
    out(f"  G1 null calibration (kill criterion) ........... {g1}   [expectation: {calib} scaling]")
    out(f"  G2 precision: top {TOP_BINS} bins hi/lo <= {PRECISION_MAX_RATIO}, top-bin lower bound > 1 ... "
        f"{'PASS' if g2 else 'FAIL'}   (hi/lo = {', '.join(fmt(x, 2) for x in ratio)})")
    out(f"  G3 centre consistency .......................... {g3}")
    if not g0:
        overall = "FIX PIPELINE - QC failed; results are not interpretable yet"
    elif g1 == "FAIL":
        overall = "KILL - the mutation model does not reproduce synonymous de novos; expected counts cannot be trusted"
    elif not g2:
        overall = "UNDERPOWERED - signal too imprecise for ACMG-scale calibration"
    elif g1 == "WARN" or g3 == "WARN":
        overall = "GO WITH CAUTION - investigate the WARN items before the full study"
    else:
        overall = "GO"
    out(f"\n  OVERALL: {overall}")
    out(f"\n  key readout (not a gate): within-gene LR in AM {edges[nb - 1]:.2f}-1.00 = {fmt(s_w['est'][nb - 1], 2)} "
        f"[{fmt(s_w['lo95'][nb - 1], 2)}, {fmt(s_w['hi95'][nb - 1], 2)}]; marginal = {fmt(s_m['est'][nb - 1], 2)}")
    out(f"\n  finished in {time.time() - t0:.0f}s; tables and figures in {a.results}")
    with open(P(a.results, "gate_report.txt"), "w") as f:
        f.write("\n".join(rep) + "\n")
    write_run_info(a, done, overall, calib)
    try:
        plot(a, curves, null_rows, edges)
    except Exception as e:  # plotting must never hide the verdict
        print(f"(figures skipped: {e})")


def write_run_info(a, chroms, verdict, calib):
    info = {"verdict": verdict, "expectation_scaling": calib, "chromosomes": chroms, "python": platform.python_version(),
            "numpy": np.__version__, "pandas": pd.__version__, "seed": SEED, "n_boot": a.n_boot,
            "primary_model": PRIMARY_MODEL, "inputs": {}}
    for rel in [DNV_FILE, GENE_FILE, DDG2P_FILE, AM_FILE, GTF_FILE, CHAIN_FILE]:
        f = P(a.data, rel)
        if os.path.exists(f):
            info["inputs"][rel] = {"bytes": os.path.getsize(f), "sha256": sha256(f)}
    try:
        info["DeNovoWEST_commit"] = subprocess.run(["git", "-C", P(a.data, "DeNovoWEST"), "rev-parse", "HEAD"],
                                                   capture_output=True, text=True).stdout.strip()
    except OSError:
        pass
    info["gate_py_sha256"] = sha256(os.path.abspath(__file__))
    for mod in ("scipy", "matplotlib", "pyliftover"):
        try:
            info[mod] = __import__(mod).__version__
        except Exception:
            pass
    for tool in ("bcftools", "bgzip"):
        try:
            info[tool] = subprocess.run([tool, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
        except (OSError, IndexError):
            pass
    info["work_files"] = {}
    for f in sorted(glob.glob(P(a.work, "roulette", "chr*.tsv.gz")) + glob.glob(P(a.work, "regions", "chr*.regions.tsv"))):
        info["work_files"][os.path.relpath(f, a.work)] = {"bytes": os.path.getsize(f), "sha256": sha256(f)}
    with open(P(a.results, "run_info.json"), "w") as f:
        json.dump(info, f, indent=1)


def plot(a, curves, null_rows, edges):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mids = np.array([(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)])
    xs = np.arange(len(mids))
    labels = [f"{edges[i]:.2f}" for i in range(len(edges) - 1)]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True)
    for ax, kind, title in ((axes[0], "marginal", "Marginal LR (gene + variant)"),
                            (axes[1], "within_gene", "Within-gene LR (variant only)")):
        for j, (an, col) in enumerate((("all_genes", "C0"), ("ddg2p_monoallelic", "C3"), ("shet_high", "C2"))):
            s = curves.get((an, kind))
            if s is None:
                continue
            est = np.clip(s["est"], 0.02, None)
            lo = np.clip(s["lo95"], 0.02, None)
            hi = np.clip(s["hi95"], 0.02, None)
            off = (j - 1) * 0.18
            ax.errorbar(xs + off, est, yerr=[np.maximum(est - lo, 0), np.maximum(hi - est, 0)], fmt="o-", ms=4, lw=1,
                        capsize=2, color=col, label=an)
        for _, cut in PP3[:3]:
            ax.axhline(cut, color="0.5", ls="--", lw=0.8)
            ax.axhline(1 / cut, color="0.5", ls=":", lw=0.8)
        ax.axhline(1, color="k", lw=0.6)
        ax.set_yscale("log")
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, rotation=45)
        ax.set_xlabel("AlphaMissense bin (lower edge)")
        ax.set_title(title)
    axes[0].set_ylabel("likelihood ratio (95% CI); dashed = PP3 2.406/5.79/33.5, dotted = BP4")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(P(a.results, "lr_curves.png"), dpi=150)
    plt.close(fig)
    nr = pd.DataFrame(null_rows, columns=["grouping", "group", "O", "E", "OE", "lo", "hi", "p", "verdict"])
    groups = list(dict.fromkeys(nr["grouping"]))
    fig, axes = plt.subplots(1, len(groups), figsize=(3.2 * len(groups), 4.2), sharey=True)
    for ax, gname in zip(np.atleast_1d(axes), groups):
        g = nr[nr["grouping"] == gname]
        x = np.arange(len(g))
        ax.axhspan(1 - NULL_PASS, 1 + NULL_PASS, color="0.9")
        ax.errorbar(x, g["OE"], yerr=[np.maximum(g["OE"] - g["lo"], 0), np.maximum(g["hi"] - g["OE"], 0)], fmt="o", ms=4,
                    capsize=2)
        ax.axhline(1, color="k", lw=0.6)
        ax.set_title(gname, fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels([str(i + 1) for i in x], fontsize=7)
    np.atleast_1d(axes)[0].set_ylabel("synonymous O/E (95% CI)")
    fig.tight_layout()
    fig.savefig(P(a.results, "null_checks.png"), dpi=150)
    plt.close(fig)


# ======================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.environ.get("DATA", P(HERE, "data")))
    ap.add_argument("--work", default=os.environ.get("WORK", P(HERE, "work")))
    ap.add_argument("--results", default=os.environ.get("RESULTS", P(HERE, "results")))
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check-inputs")
    sub.add_parser("prep-dnv")
    s = sub.add_parser("prep-am")
    s.add_argument("--chunksize", type=int, default=2_000_000)
    s.add_argument("--force", action="store_true")
    sub.add_parser("prep-gtf")
    s = sub.add_parser("sites")
    s.add_argument("--chroms", default="1-22")
    s.add_argument("--force", action="store_true")
    s.add_argument("--no-ok-check", action="store_true", help=argparse.SUPPRESS)
    s = sub.add_parser("run")
    s.add_argument("--n-boot", type=int, default=N_BOOT)
    s.add_argument("--allow-partial", action="store_true", help="run on the chromosomes processed so far")
    s.add_argument("--qc-only", action="store_true", help="pipeline QC only; computes no calibration result")
    a = ap.parse_args()
    {"check-inputs": cmd_check_inputs, "prep-dnv": cmd_prep_dnv, "prep-am": cmd_prep_am, "prep-gtf": cmd_prep_gtf,
     "sites": cmd_sites, "run": cmd_run}[a.cmd](a)


if __name__ == "__main__":
    main()
