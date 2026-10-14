#!/usr/bin/env python3
"""Amendment A9: how much of a non-coding predictor's AUROC against population "benign" variants is rarity alone?

Benchmarks of non-coding variant-effect predictors commonly use ClinVar P/LP variants as positives and common gnomAD
variants as negatives. Any score that tracks selection or conservation then separates a RANDOM rare variant from the
common ones, so the null AUROC of each score is above 0.5 by a score-specific amount (seen first on RNU2-2, A8 rows of
PREREG s14). A9 measures that inflation across ClinVar's non-coding P/LP SNVs, genome-wide, for widely used scores.

    python a9_inflation.py select      # ClinVar VCF -> a9/positives.tsv (non-coding P/LP SNVs, >= 1 star, capped per class)
    python a9_inflation.py fetch       # gnomAD, CADD, GPN-MSA, ReMM (+ phyloP447 if pyBigWig) in +-100 bp windows
    python a9_inflation.py eval        # -> results/a9_report.txt, results/a9_results.tsv

Design (pre-registered, PREREG s15 A9):
  Positives: ClinVar GRCh38 VCF (fileDate recorded) SNVs with CLNSIG Pathogenic / Likely_pathogenic /
    Pathogenic/Likely_pathogenic, review >= 1 star, and molecular consequences only non-coding (any coding or canonical
    splice-site term excludes). Class by first match: 5'UTR, 3'UTR, ncRNA exon (non-coding_transcript_variant), intronic,
    upstream. At most 1,500 per class (random, seed 20261017).
  Negatives: gnomAD v4.1 genome PASS SNVs with AC >= 9 within +-100 bp of a positive, not ClinVar P/LP.
  Rarity-only null: for each positive, one SNV drawn from its own window that is neither a negative nor ClinVar
    P/LP (rare or absent in gnomAD); 200 draws. Null AUROC = AUROC of the score for drawn SNVs vs the same negatives.
  Metrics per class (>= 50 positives) and pooled: raw AUROC [1,000 bootstraps], null mean [2.5, 97.5 pct],
    inflation = null mean - 0.5, corrected = raw - null mean.
  H1 (primary): pooled inflation >= 0.05 for at least 2 of CADD, GPN-MSA, ReMM.
  H2: in at least one class with >= 200 positives, the ranking of tools by raw AUROC differs from that by corrected.
"""
import argparse
import gzip
import json
import os
import re
import sys
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snrna as S  # noqa: E402

OUT = os.path.join(HERE, "a9")
RES = os.path.join(HERE, "results")
CLINVAR = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
PHYLOP_BW = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP447way/hg38.phyloP447way.bw"
SEED, CAP, FLANK, NULL_DRAWS = 20261017, 1500, 100, 200
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
NO_STAR = {"no_assertion_criteria_provided", "no_classification_provided", "no_interpretation_for_the_single_variant",
           "no_classifications_from_unflagged_records", "no_assertion_provided"}
EXCLUDE = {"missense_variant", "synonymous_variant", "nonsense", "stop_lost", "initiator_codon_variant", "frameshift_variant",
           "inframe_insertion", "inframe_deletion", "inframe_indel", "splice_donor_variant", "splice_acceptor_variant",
           "stop_gained", "start_lost", "coding_sequence_variant"}
CLASSES = [("5_prime_UTR_variant", "5'UTR"), ("3_prime_UTR_variant", "3'UTR"), ("non-coding_transcript_variant", "ncRNA exon"),
           ("intron_variant", "intronic"), ("genic_upstream_transcript_variant", "upstream"), ("upstream_gene_variant", "upstream")]
TOOLS = [("cadd_phred", "CADD v1.7"), ("gpn_msa", "GPN-MSA"), ("remm", "ReMM v0.4"), ("phylop447", "phyloP447")]
CORE = ("cadd_phred", "gpn_msa", "remm")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ---------------------------------------------------------------------------------------------------- select

def cmd_select(a):
    os.makedirs(OUT, exist_ok=True)
    vcf = os.path.join(S.DATA, "clinvar_GRCh38.vcf.gz")
    if not os.path.exists(vcf) or a.refresh:
        log(f"downloading {CLINVAR}")
        urllib.request.urlretrieve(CLINVAR, vcf)
    rows, fdate = [], ""
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##fileDate"):
                fdate = line.strip().split("=")[1]
            if line.startswith("#"):
                continue
            c = line.split("\t", 8)
            ref, alt = c[3], c[4]
            if len(ref) != 1 or len(alt) != 1 or alt not in "ACGT":
                continue
            info = dict(kv.split("=", 1) for kv in c[7].strip().split(";") if "=" in kv)
            if info.get("CLNSIG") not in PLP or set(info.get("CLNREVSTAT", "").split(",")) & NO_STAR:
                continue
            mc = {t.split("|")[1] for t in info.get("MC", "").split(",") if "|" in t}
            if not mc or mc & EXCLUDE:
                continue
            cls = next((lab for term, lab in CLASSES if term in mc), None)
            if cls is None:
                continue
            chrom = "chr" + c[0]
            rows.append({"ID": f"{chrom}-{c[1]}-{ref}-{alt}", "chrom": chrom, "pos": int(c[1]), "ref": ref, "alt": alt,
                         "class": cls, "gene": info.get("GENEINFO", "").split(":")[0], "clnsig": info["CLNSIG"],
                         "revstat": info.get("CLNREVSTAT", ""), "mc": ",".join(sorted(mc))})
    p = pd.DataFrame(rows).drop_duplicates("ID")
    p = p[p.chrom.isin([f"chr{i}" for i in list(range(1, 23)) + ["X"]])]
    counts = p["class"].value_counts().to_dict()
    p = pd.concat([x.sample(min(len(x), CAP), random_state=SEED) for _, x in p.groupby("class")])
    p.to_csv(os.path.join(OUT, "positives.tsv"), sep="\t", index=False)
    json.dump({"clinvar_fileDate": fdate, "available_by_class": counts, "kept_by_class": p["class"].value_counts().to_dict(),
               "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, open(os.path.join(OUT, "select_info.json"), "w"), indent=1)
    log(f"ClinVar fileDate {fdate}; available {counts}; kept {p['class'].value_counts().to_dict()}")


# ----------------------------------------------------------------------------------------------------- fetch

def merged_windows(pos):
    pos = sorted(pos)
    out = [[pos[0] - FLANK, pos[0] + FLANK]]
    for p in pos[1:]:
        if p - FLANK <= out[-1][1] + 1:
            out[-1][1] = p + FLANK
        else:
            out.append([p - FLANK, p + FLANK])
    return out


def fetch_chrom(args):
    import comparators as CP
    chrom, wins = args
    parts = []
    bw = None
    try:
        import pyBigWig
        bw = pyBigWig.open(PHYLOP_BW)
    except Exception:
        bw = None
    for k, (s, e) in enumerate(wins):
        cadd = S.cadd_records(chrom, s, e)
        if not len(cadd):
            continue
        d = cadd.assign(ID=cadd.chrom + "-" + cadd.pos.astype(str) + "-" + cadd.ref + "-" + cadd.alt)[["ID", "chrom", "pos", "ref", "alt", "cadd_phred"]]
        gn = S.gnomad_records(chrom, s, e)
        gn = gn.assign(ID=gn.chrom + "-" + gn.pos.astype(str) + "-" + gn.ref + "-" + gn.alt).drop_duplicates("ID").set_index("ID")
        d["filter"] = d.ID.map(gn["filter"])
        d["AC"] = d.ID.map(gn.AC).fillna(0).astype(int)
        gpn = {(int(r[1]), r[2], r[3]): float(r[4]) for r in CP.tabix_rows(CP.GPN, chrom.replace("chr", ""), s, e)}
        remm = {int(r[1]): float(r[2]) for r in CP.tabix_rows(CP.REMM, chrom, s, e)}
        d["gpn_msa"] = [-gpn[x] if x in gpn else np.nan for x in zip(d.pos, d.ref, d.alt)]
        d["remm"] = d.pos.map(remm)
        if bw is not None:
            vals = bw.values(chrom, s - 1, e)
            d["phylop447"] = [vals[p - s] if 0 <= p - s < len(vals) else np.nan for p in d.pos]
        else:
            d["phylop447"] = np.nan
        d["win"] = f"{chrom}:{s}-{e}"
        parts.append(d)
        if k % 50 == 0:
            print(f"  {chrom}: window {k + 1}/{len(wins)}", flush=True)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def cmd_fetch(a):
    p = pd.read_csv(os.path.join(OUT, "positives.tsv"), sep="\t")
    jobs = [(c, merged_windows(x.pos.tolist())) for c, x in p.groupby("chrom")]
    log(f"{len(p)} positives in {sum(len(w) for _, w in jobs)} merged windows on {len(jobs)} chromosomes; {a.workers} workers")
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        parts = [d for d in ex.map(fetch_chrom, jobs) if len(d)]
    d = pd.concat(parts, ignore_index=True).drop_duplicates("ID")
    d.to_csv(os.path.join(OUT, "window_snvs.tsv.gz"), sep="\t", index=False, float_format="%.5g", compression="gzip")
    log(f"{len(d)} SNVs in windows; coverage " + ", ".join(f"{c} {d[c].notna().mean():.2f}" for c, _ in TOOLS))


# ------------------------------------------------------------------------------------------------------ eval

def cmd_eval(a):
    p = pd.read_csv(os.path.join(OUT, "positives.tsv"), sep="\t")
    d = pd.read_csv(os.path.join(OUT, "window_snvs.tsv.gz"), sep="\t")
    plp_all = set(p.ID)
    d["pos_flag"] = d.ID.isin(plp_all)
    d["neg_flag"] = (d["filter"] == "PASS") & (d.AC >= 9) & ~d.pos_flag
    d["pool_flag"] = ~d.neg_flag & ~d.pos_flag
    # each positive's own +-100 bp neighbourhood (for the window-matched null)
    d = d.sort_values(["chrom", "pos"]).reset_index(drop=True)
    by_chrom = {c: x for c, x in d.groupby("chrom")}
    rng = np.random.default_rng(SEED)
    tools = [(c, l) for c, l in TOOLS if d[c].notna().mean() > 0.5]
    rows, rep = [], [f"A9: rarity-only inflation of non-coding predictor AUROCs (ClinVar P/LP vs gnomAD PASS AC >= 9, +-{FLANK} bp)", ""]
    units = [(cl, x) for cl, x in p.groupby("class") if len(x) >= 50] + [("pooled", p)]
    for unit, pu in units:
        pos_idx, neg_idx, pools = [], set(), []
        for r in pu.itertuples():
            x = by_chrom.get(r.chrom)
            if x is None:
                continue
            w = x[(x.pos >= r.pos - FLANK) & (x.pos <= r.pos + FLANK)]
            hit = w.index[w.ID == r.ID]
            pool = w.index[w.pool_flag].to_numpy()
            if not len(hit) or not len(pool):
                continue
            pos_idx.append(hit[0])
            neg_idx.update(w.index[w.neg_flag])
            pools.append(pool)
        neg_idx = np.array(sorted(neg_idx))
        if len(pos_idx) < 50 or len(neg_idx) < 50:
            continue
        pos_idx = np.array(pos_idx)
        y = np.r_[np.ones(len(pos_idx), bool), np.zeros(len(neg_idx), bool)]
        rep.append(f"{unit}: {len(pos_idx)} P/LP vs {len(neg_idx)} common SNVs")
        res = {}
        for col, lab in tools:
            s = d[col].to_numpy(float)
            s = np.where(np.isnan(s), np.nanmedian(s), s)
            raw = S.auroc(y, s[np.r_[pos_idx, neg_idx]])
            lo, hi = S.boot_auroc(y, s[np.r_[pos_idx, neg_idx]], a.boot)
            nulls = np.array([S.auroc(y, s[np.r_[np.array([rng.choice(pl) for pl in pools]), neg_idx]]) for _ in range(NULL_DRAWS)])
            nm, nlo, nhi = float(nulls.mean()), *np.percentile(nulls, [2.5, 97.5])
            res[col] = (raw, raw - nm)
            rows.append({"unit": unit, "tool": lab, "n_P": len(pos_idx), "n_N": len(neg_idx), "raw_auroc": raw, "raw_lo": lo,
                         "raw_hi": hi, "null_mean": nm, "null_lo": nlo, "null_hi": nhi, "inflation": nm - 0.5, "corrected": raw - nm})
            rep.append(f"  {lab:12s} raw {raw:.3f} [{lo:.3f}, {hi:.3f}]  rarity null {nm:.3f} [{nlo:.3f}, {nhi:.3f}]  "
                       f"inflation {nm - 0.5:+.3f}  corrected {raw - nm:+.3f}")
        rr = sorted(res, key=lambda c: -res[c][0])
        rc = sorted(res, key=lambda c: -res[c][1])
        rep.append(f"  ranking raw: {' > '.join(rr)} | corrected: {' > '.join(rc)}{'  (CHANGED)' if rr != rc else ''}")
        rep.append("")
        for r in rows:
            if r["unit"] == unit:
                r["rank_changed"] = rr != rc
    t = pd.DataFrame(rows)
    pooled = t[(t.unit == "pooled") & t.tool.isin([l for c, l in TOOLS if c in CORE])]
    h1 = int((pooled.inflation >= 0.05).sum()) >= 2
    big = t[(t.n_P >= 200) & (t.unit != "pooled")]
    h2 = bool(big.rank_changed.any()) if len(big) else False
    rep.append(f"H1 (pooled inflation >= 0.05 for >= 2 of CADD, GPN-MSA, ReMM): {'MET' if h1 else 'NOT MET'} "
               f"({', '.join(f'{r.tool} {r.inflation:+.3f}' for r in pooled.itertuples())})")
    rep.append(f"H2 (tool ranking changes after correction in a class with >= 200 positives): {'MET' if h2 else 'NOT MET'}")
    os.makedirs(RES, exist_ok=True)
    t.to_csv(os.path.join(RES, "a9_results.tsv"), sep="\t", index=False, float_format="%.4f")
    open(os.path.join(RES, "a9_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["select", "fetch", "eval"])
    ap.add_argument("--refresh", action="store_true", help="select: re-download the ClinVar VCF")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    {"select": cmd_select, "fetch": cmd_fetch, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
