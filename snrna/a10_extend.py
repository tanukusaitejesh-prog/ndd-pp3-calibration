#!/usr/bin/env python3
"""Amendment A10: does the rarity-only inflation (A9) reach missense benchmarks and SpliceAI, and what removes it?

    python a10_extend.py select                 # local: ClinVar VCF -> a10/missense_positives.tsv, a10/blb_*.tsv (frozen)
    python a10_extend.py missense_fetch         # DGX: windows around missense positives: CADD, gnomAD, GPN-MSA, ReMM,
                                                #      phyloP447 (UCSC API) and AlphaMissense (bulk file, filtered)
    python a10_extend.py spliceai_prep          # DGX: VCF of A9 positives, common negatives and a fixed pool sample
    (run SpliceAI in its own env, see run_next.sh part5)
    python a10_extend.py spliceai_parse         # DGX: SpliceAI VCF -> max delta score per SNV
    python a10_extend.py eval                   # DGX: -> results/a10_report.txt, results/a10_results.tsv

Design (pre-registered, PREREG s15 A10; same null and metrics as A9):
  Missense: ClinVar P/LP SNVs (>= 1 star) whose molecular consequence includes missense_variant and no splice-site,
    nonsense, stop or start term; at most 1,500 (seed 20261019). Windows +-100 bp. Analysis restricted to SNVs with an
    AlphaMissense score (i.e. missense in AlphaMissense's canonical transcripts): positives, negatives (gnomAD v4.1 genome
    PASS AC >= 9, not P/LP) and the rarity pool (neither) are all missense. Tools: AlphaMissense, CADD, GPN-MSA, phyloP447.
  SpliceAI: the A9 positives and negatives plus, per positive, up to 10 pool SNVs from its own window (seed 20261019),
    scored with SpliceAI (max delta score, distance 50; unscored = 0). Null drawn from that sample.
  Fix metrics, for every tool in both sets: (a) rarity-matched AUROC = positives vs ALL pool SNVs in their windows;
    (b) AUROC vs ClinVar Benign/Likely benign SNVs (>= 1 star) inside the same windows (same consequence class).
  H1: AlphaMissense pooled inflation >= 0.05. H2: SpliceAI inflation >= 0.05 on intronic positives.
  H3 (descriptive): for each tool, the rarity-matched AUROC and the B/LB AUROC lie closer to (null-)corrected AUROC + 0.5
    than the raw AUROC does.
"""
import argparse
import gzip
import json
import os
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import a9_inflation as A  # noqa: E402
import snrna as S  # noqa: E402

OUT = os.path.join(HERE, "a10")
RES = os.path.join(HERE, "results")
SEED, CAP, FLANK, POOL_PER_POS = 20261019, 1500, 100, 10
AM_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
MISSENSE_EXCLUDE = A.EXCLUDE - {"missense_variant"} | {"nonsense", "stop_gained", "stop_lost", "start_lost", "initiator_codon_variant"}
BLB = {"Benign", "Likely_benign", "Benign/Likely_benign"}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def clinvar_rows(vcf):
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t", 8)
            if len(c[3]) != 1 or len(c[4]) != 1 or c[4] not in "ACGT":
                continue
            info = dict(kv.split("=", 1) for kv in c[7].strip().split(";") if "=" in kv)
            if set(info.get("CLNREVSTAT", "").split(",")) & A.NO_STAR:
                continue
            mc = {t.split("|")[1] for t in info.get("MC", "").split(",") if "|" in t}
            yield {"ID": f"chr{c[0]}-{c[1]}-{c[3]}-{c[4]}", "chrom": "chr" + c[0], "pos": int(c[1]), "ref": c[3], "alt": c[4],
                   "clnsig": info.get("CLNSIG", ""), "mc": mc, "gene": info.get("GENEINFO", "").split(":")[0]}


def nc_class(mc):
    if not mc or mc & A.EXCLUDE:
        return None
    return next((lab for term, lab in A.CLASSES if term in mc), None)


def in_windows(df, centers):
    """Rows of df within +-FLANK of any centre on the same chromosome."""
    keep = np.zeros(len(df), bool)
    for c, x in centers.groupby("chrom"):
        pos = np.sort(x.pos.to_numpy())
        m = (df.chrom == c).to_numpy()
        if not m.any():
            continue
        p = df.pos.to_numpy()[m]
        i = np.searchsorted(pos, p)
        lo = np.abs(p - pos[np.clip(i - 1, 0, len(pos) - 1)]) <= FLANK
        hi = np.abs(p - pos[np.clip(i, 0, len(pos) - 1)]) <= FLANK
        keep[np.flatnonzero(m)] = lo | hi
    return df[keep]


def cmd_select(a):
    os.makedirs(OUT, exist_ok=True)
    vcf = os.path.join(S.DATA, "clinvar_GRCh38.vcf.gz")
    rows = list(clinvar_rows(vcf))
    cv = pd.DataFrame(rows)
    cv = cv[cv.chrom.isin([f"chr{i}" for i in list(range(1, 23)) + ["X"]])]
    plp = cv.clnsig.isin(A.PLP)
    mis = plp & cv.mc.map(lambda m: "missense_variant" in m and not (m & MISSENSE_EXCLUDE))
    mp = cv[mis].drop_duplicates("ID")
    n_avail = len(mp)
    mp = mp.sample(min(len(mp), CAP), random_state=SEED).assign(mc=lambda d: d.mc.map(lambda m: ",".join(sorted(m))))
    mp.to_csv(os.path.join(OUT, "missense_positives.tsv"), sep="\t", index=False)
    blb = cv[cv.clnsig.isin(BLB)].drop_duplicates("ID")
    a9p = pd.read_csv(os.path.join(HERE, "a9", "positives.tsv"), sep="\t")
    nc = blb[blb.mc.map(nc_class).notna()].assign(**{"class": lambda d: d.mc.map(nc_class)})
    nc = in_windows(nc, a9p).assign(mc=lambda d: d.mc.map(lambda m: ",".join(sorted(m))))
    nc.to_csv(os.path.join(OUT, "blb_noncoding.tsv"), sep="\t", index=False)
    mb = blb[blb.mc.map(lambda m: "missense_variant" in m and not (m & MISSENSE_EXCLUDE))]
    mb = in_windows(mb, mp).assign(mc=lambda d: d.mc.map(lambda m: ",".join(sorted(m))))
    mb.to_csv(os.path.join(OUT, "blb_missense.tsv"), sep="\t", index=False)
    info = {"missense_PLP_available": int(n_avail), "missense_PLP_kept": len(mp), "blb_noncoding_in_A9_windows": len(nc),
            "blb_noncoding_by_class": nc["class"].value_counts().to_dict(), "blb_missense_in_windows": len(mb),
            "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(info, open(os.path.join(OUT, "select_info.json"), "w"), indent=1)
    log(json.dumps(info))


# ---------------------------------------------------------------------------------------------------- fetch

def fetch_windows(wins, workers, timeout):
    import multiprocessing as mp
    parts, pending = [], wins
    for rnd in range(1, 4):
        if not pending:
            break
        failed = []
        pool = mp.get_context("fork").Pool(workers)
        jobs = [(w, pool.apply_async(A.fetch_window, (w,))) for w in pending]
        for k, (w, j) in enumerate(jobs, 1):
            try:
                d = j.get(timeout=timeout)
                if len(d):
                    parts.append(d)
            except Exception as ex:
                failed.append(w)
                log(f"  round {rnd}: {w[0]}:{w[1]}-{w[2]} not done ({type(ex).__name__}); will retry")
            if k % 25 == 0 or k == len(jobs):
                log(f"  round {rnd}: {k}/{len(jobs)} windows checked, {len(failed)} pending retry")
        pool.terminate()
        pending = failed
    return pd.concat(parts, ignore_index=True).drop_duplicates("ID"), pending


def add_phylop(d):
    import comparators as CP
    wins = sorted(set(d.win))
    vals = {}
    for k, w in enumerate(wins, 1):
        chrom, r = w.split(":")
        s, e = map(int, r.split("-"))
        try:
            vals[w] = CP.ucsc_track("phyloP447way", chrom, s, e)
        except Exception:
            vals[w] = {}
        time.sleep(0.25)
        if k % 50 == 0 or k == len(wins):
            log(f"  phyloP {k}/{len(wins)} windows")
    d["phylop447"] = [vals[w].get(p, np.nan) for w, p in zip(d.win, d.pos)]
    return d


def add_alphamissense(d):
    path = os.path.join(S.DATA, "AlphaMissense_hg38.tsv.gz")
    if not os.path.exists(path):
        log(f"downloading {AM_URL} (~640 MB)")
        urllib.request.urlretrieve(AM_URL, path)
    want = set(d.ID)
    am = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            k = f"{c[0]}-{c[1]}-{c[2]}-{c[3]}"
            if k in want:
                v = float(c[8])
                am[k] = max(am.get(k, v), v)  # several transcripts: keep the maximum
    d["alphamissense"] = d.ID.map(am)
    log(f"AlphaMissense scores for {d.alphamissense.notna().sum()} of {len(d)} window SNVs")
    return d


def cmd_missense_fetch(a):
    p = pd.read_csv(os.path.join(OUT, "missense_positives.tsv"), sep="\t")
    wins = [(c, s, e) for c, x in p.groupby("chrom") for s, e in A.merged_windows(x.pos.tolist())]
    log(f"{len(p)} missense positives in {len(wins)} windows")
    d, skipped = fetch_windows(wins, a.workers, a.timeout)
    d = add_alphamissense(add_phylop(d))
    d.to_csv(os.path.join(OUT, "missense_window_snvs.tsv.gz"), sep="\t", index=False, float_format="%.5g", compression="gzip")
    json.dump({"windows": len(wins), "skipped": [f"{c}:{s}-{e}" for c, s, e in skipped]},
              open(os.path.join(OUT, "missense_fetch_info.json"), "w"), indent=1)
    log(f"{len(d)} SNVs; coverage " + ", ".join(f"{c} {d[c].notna().mean():.2f}" for c in ("cadd_phred", "gpn_msa", "phylop447", "alphamissense")))


# -------------------------------------------------------------------------------------------------- SpliceAI

def cmd_spliceai_prep(a):
    rng = np.random.default_rng(SEED)
    p = pd.read_csv(os.path.join(HERE, "a9", "positives.tsv"), sep="\t")
    d = pd.read_csv(os.path.join(HERE, "a9", "window_snvs.tsv.gz"), sep="\t")
    d = d.sort_values(["chrom", "pos"]).reset_index(drop=True)
    pos = set(p.ID)
    neg = (d["filter"] == "PASS") & (d.AC >= 9) & ~d.ID.isin(pos)
    pool = ~neg & ~d.ID.isin(pos)
    keep = set(d.ID[d.ID.isin(pos) | neg])
    by_chrom = {c: x for c, x in d[pool].groupby("chrom")}
    for r in p.itertuples():
        x = by_chrom.get(r.chrom)
        if x is None:
            continue
        w = x[(x.pos >= r.pos - FLANK) & (x.pos <= r.pos + FLANK)]
        if len(w):
            keep.update(rng.choice(w.ID.to_numpy(), min(POOL_PER_POS, len(w)), replace=False))
    v = d[d.ID.isin(keep)][["chrom", "pos", "ID", "ref", "alt"]]
    order = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23}
    v = v.assign(o=v.chrom.map(order)).sort_values(["o", "pos"])
    path = os.path.join(OUT, "spliceai_in.vcf")
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in sorted(order, key=order.get):
            f.write(f"##contig=<ID={c.replace('chr', '')}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in v.itertuples():
            f.write(f"{r.chrom.replace('chr', '')}\t{r.pos}\t{r.ID}\t{r.ref}\t{r.alt}\t.\t.\t.\n")
    log(f"wrote {path}: {len(v)} SNVs ({len(pos)} positives, {int(neg.sum())} common negatives, pool sample)")


def cmd_spliceai_parse(a):
    rows = []
    with open(os.path.join(OUT, "spliceai_out.vcf")) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            best = 0.0
            for kv in c[7].split(";"):
                if kv.startswith("SpliceAI="):
                    for ann in kv[9:].split(","):
                        x = ann.split("|")
                        best = max([best] + [float(v) for v in x[2:6] if v not in (".", "")])
            rows.append({"ID": c[2], "spliceai": best})
    s = pd.DataFrame(rows)
    s.to_csv(os.path.join(OUT, "spliceai_scores.tsv"), sep="\t", index=False)
    log(f"{len(s)} SNVs parsed; {int((s.spliceai > 0).sum())} with a delta score > 0")


# ------------------------------------------------------------------------------------------------------ eval

def inflation(p, d, tools, unit_col, rng, boot, extra_neg=None):
    """A9 metrics per unit and pooled, plus the fix metrics: rarity-matched AUROC (vs all pool SNVs in the
    positives' windows) and AUROC vs ClinVar B/LB in the same windows (extra_neg: set of IDs)."""
    d = d.sort_values(["chrom", "pos"]).reset_index(drop=True)
    plp = set(p.ID)
    d["pos_flag"] = d.ID.isin(plp)
    d["neg_flag"] = (d["filter"] == "PASS") & (d.AC >= 9) & ~d.pos_flag
    d["pool_flag"] = ~d.neg_flag & ~d.pos_flag & ~d.ID.isin(extra_neg or set())
    by_chrom = {c: x for c, x in d.groupby("chrom")}
    rows, rep = [], []
    units = [(u, x) for u, x in p.groupby(unit_col) if len(x) >= 50] + [("pooled", p)] if unit_col else [("pooled", p)]
    for unit, pu in units:
        pos_idx, neg_idx, pool_all, pools, blb_set = [], set(), set(), [], set()
        for r in pu.itertuples():
            x = by_chrom.get(r.chrom)
            if x is None:
                continue
            w = x[(x.pos >= r.pos - FLANK) & (x.pos <= r.pos + FLANK)]
            hit = w.index[w.ID == r.ID]
            pl = w.index[w.pool_flag].to_numpy()
            if not len(hit) or not len(pl):
                continue
            pos_idx.append(hit[0]); neg_idx.update(w.index[w.neg_flag]); pool_all.update(pl); pools.append(pl)
            blb_set.update(w.index[w.ID.isin(extra_neg or set())])
        pos_idx, neg_idx, pool_all = np.array(pos_idx), np.array(sorted(neg_idx)), np.array(sorted(pool_all))
        if len(pos_idx) < 30 or len(neg_idx) < 30:
            continue
        blb_idx = np.array(sorted(blb_set), dtype=int)
        y = np.r_[np.ones(len(pos_idx), bool), np.zeros(len(neg_idx), bool)]
        rep.append(f"{unit}: {len(pos_idx)} P/LP vs {len(neg_idx)} common SNVs; pool {len(pool_all)}; ClinVar B/LB {len(blb_idx)}")
        for col, lab in tools:
            s = d[col].to_numpy(float)
            s = np.where(np.isnan(s), np.nanmedian(s), s)
            raw = S.auroc(y, s[np.r_[pos_idx, neg_idx]])
            nulls = np.array([S.auroc(y, s[np.r_[np.array([rng.choice(pl) for pl in pools]), neg_idx]]) for _ in range(A.NULL_DRAWS)])
            nm = float(nulls.mean())
            matched = S.auroc(np.r_[np.ones(len(pos_idx), bool), np.zeros(len(pool_all), bool)], s[np.r_[pos_idx, pool_all]])
            vs_blb = S.auroc(np.r_[np.ones(len(pos_idx), bool), np.zeros(len(blb_idx), bool)], s[np.r_[pos_idx, blb_idx]]) \
                if len(blb_idx) >= 20 else np.nan
            rows.append({"unit": unit, "tool": lab, "n_P": len(pos_idx), "n_N": len(neg_idx), "raw_auroc": raw, "null_mean": nm,
                         "inflation": nm - 0.5, "corrected": raw - nm, "rarity_matched_auroc": matched, "auroc_vs_clinvar_blb": vs_blb,
                         "n_blb": len(blb_idx)})
            rep.append(f"  {lab:14s} raw {raw:.3f}  null {nm:.3f}  inflation {nm - 0.5:+.3f}  corrected {raw - nm:+.3f}  | "
                       f"rarity-matched {matched:.3f}  vs ClinVar B/LB {vs_blb:.3f}")
        rep.append("")
    return pd.DataFrame(rows), rep


def cmd_eval(a):
    rng = np.random.default_rng(SEED)
    tabs, rep = [], ["A10: rarity-only inflation in missense and SpliceAI benchmarks, and two fixes", ""]
    # missense
    mp = os.path.join(OUT, "missense_window_snvs.tsv.gz")
    if os.path.exists(mp):
        d = pd.read_csv(mp, sep="\t")
        d = d[d.alphamissense.notna()]
        p = pd.read_csv(os.path.join(OUT, "missense_positives.tsv"), sep="\t")
        p = p[p.ID.isin(set(d.ID))]
        blb = set(pd.read_csv(os.path.join(OUT, "blb_missense.tsv"), sep="\t").ID)
        rep.append(f"== MISSENSE ({len(p)} ClinVar P/LP with an AlphaMissense score)")
        t, r = inflation(p, d, [("alphamissense", "AlphaMissense"), ("cadd_phred", "CADD v1.7"), ("gpn_msa", "GPN-MSA"),
                                ("phylop447", "phyloP447")], None, rng, a.boot, blb)
        tabs.append(t.assign(set="missense")); rep += r
        am = t[(t.unit == "pooled") & (t.tool == "AlphaMissense")]
        if len(am):
            rep.append(f"H1 (AlphaMissense pooled inflation >= 0.05): {'MET' if am.inflation.iloc[0] >= 0.05 else 'NOT MET'} "
                       f"({am.inflation.iloc[0]:+.3f})\n")
    # SpliceAI on the A9 set
    sp = os.path.join(OUT, "spliceai_scores.tsv")
    if os.path.exists(sp):
        s = pd.read_csv(sp, sep="\t")
        d = pd.read_csv(os.path.join(HERE, "a9_phylop", "window_snvs.tsv.gz"), sep="\t") if \
            os.path.exists(os.path.join(HERE, "a9_phylop", "window_snvs.tsv.gz")) else pd.read_csv(os.path.join(HERE, "a9", "window_snvs.tsv.gz"), sep="\t")
        d = d.merge(s, on="ID", how="inner")
        p = pd.read_csv(os.path.join(HERE, "a9", "positives.tsv"), sep="\t")
        blb = set(pd.read_csv(os.path.join(OUT, "blb_noncoding.tsv"), sep="\t").ID)
        rep.append(f"== SPLICEAI on the A9 non-coding set (scored subset: {len(d)} SNVs)")
        tools = [("spliceai", "SpliceAI"), ("cadd_phred", "CADD v1.7"), ("gpn_msa", "GPN-MSA"), ("remm", "ReMM v0.4")] + \
                ([("phylop447", "phyloP447")] if "phylop447" in d.columns and d.phylop447.notna().mean() > 0.5 else [])
        t, r = inflation(p, d, tools, "class", rng, a.boot, blb)
        tabs.append(t.assign(set="noncoding_spliceai")); rep += r
        sa = t[(t.unit == "intronic") & (t.tool == "SpliceAI")]
        if len(sa):
            rep.append(f"H2 (SpliceAI intronic inflation >= 0.05): {'MET' if sa.inflation.iloc[0] >= 0.05 else 'NOT MET'} "
                       f"({sa.inflation.iloc[0]:+.3f})\n")
    if tabs:
        t = pd.concat(tabs)
        t.to_csv(os.path.join(RES, "a10_results.tsv"), sep="\t", index=False, float_format="%.4f")
        h3 = t.dropna(subset=["auroc_vs_clinvar_blb"])
        closer = ((h3.rarity_matched_auroc - (h3.corrected + 0.5)).abs() < (h3.raw_auroc - (h3.corrected + 0.5)).abs()).mean() if len(h3) else np.nan
        rep.append(f"H3 (descriptive): rarity-matched AUROC closer to corrected+0.5 than raw is, in {closer:.0%} of tool x unit rows")
    open(os.path.join(RES, "a10_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["select", "missense_fetch", "spliceai_prep", "spliceai_parse", "eval"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True); os.makedirs(RES, exist_ok=True)
    {"select": cmd_select, "missense_fetch": cmd_missense_fetch, "spliceai_prep": cmd_spliceai_prep,
     "spliceai_parse": cmd_spliceai_parse, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
