#!/usr/bin/env python3
"""Amendment A11: re-analysis of a published non-coding benchmark that used common variants as negatives.

NCBoost (Caron, Luo & Rausell, Genome Biology 2019, doi 10.1186/s13059-019-1634-2) compared non-coding predictors on
737 curated pathogenic non-coding SNVs (655 public, Table S1) against 10x region-matched common dbSNP SNVs (MAF > 5%,
no clinical assertion; Table S3, 6,550), and ranked ReMM above Eigen, DeepSEA and CADD. Tables are CC BY 4.0 and are
copied to a11/ as TSV (GRCh37). We re-score this exact benchmark and add the rarity-only null.

    python a11_ncboost.py fetch     # DGX: liftover to GRCh38 (pyliftover), CADD/gnomAD/GPN-MSA/ReMM (+-100 bp around
                                    #      positives for the pool; the published negatives themselves), phyloP447 (UCSC API)
    python a11_ncboost.py eval      # -> results/a11_report.txt, results/a11_results.tsv

Metrics as A9: raw AUROC (published positives vs published negatives), window-matched rarity-only null (per positive,
one SNV from its own +-100 bp window that is neither positive nor gnomAD-common, 200 draws, against the same published
negatives), inflation, corrected, and rarity-matched AUROC (positives vs all pool SNVs in their windows). Variants that
fail liftover or whose GRCh38 reference base disagrees are dropped and counted. NCBoost's own score is not openly
downloadable and is not re-scored.
H1: inflation >= 0.05 for at least 2 of CADD, ReMM, GPN-MSA, phyloP447 on the published benchmark.
H2 (descriptive): whether the published order ReMM > CADD holds for raw, corrected and rarity-matched AUROC.
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import a9_inflation as A  # noqa: E402
import a10_extend as X  # noqa: E402
import snrna as S  # noqa: E402

OUT = os.path.join(HERE, "a11")
RES = os.path.join(HERE, "results")
SEED = 20261020
TOOLS = [("cadd_phred", "CADD v1.7"), ("remm", "ReMM v0.4"), ("gpn_msa", "GPN-MSA"), ("phylop447", "phyloP447")]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def lift(df):
    from pyliftover import LiftOver
    lo = LiftOver("hg19", "hg38")
    out = []
    for r in df.itertuples():
        c = f"chr{r.chr}"
        m = lo.convert_coordinate(c, int(r.pos) - 1)
        if not m or m[0][0] != c or m[0][2] != "+":
            continue
        out.append({"chrom": c, "pos": int(m[0][1]) + 1, "ref": r.ref, "alt": r.alt, "region": r.annovar_genomic_region,
                    "pos37": int(r.pos)})
    o = pd.DataFrame(out)
    o["ID"] = o.chrom + "-" + o.pos.astype(str) + "-" + o.ref + "-" + o.alt
    return o


def cmd_fetch(a):
    p = lift(pd.read_csv(os.path.join(OUT, "ncboost_positives_GRCh37.tsv"), sep="\t", dtype={"chr": str}))
    n = lift(pd.read_csv(os.path.join(OUT, "ncboost_negatives_GRCh37.tsv"), sep="\t", dtype={"chr": str}))
    log(f"lifted: {len(p)} of 655 positives, {len(n)} of 6,550 negatives")
    wins = [(c, s, e) for c, x in p.groupby("chrom") for s, e in A.merged_windows(x.pos.tolist())]
    nwins = [(r.chrom, r.pos, r.pos) for r in n.drop_duplicates(["chrom", "pos"]).itertuples()]
    log(f"{len(wins)} positive windows (+-100 bp) and {len(nwins)} negative positions")
    d_all, skipped = X.fetch_windows(wins + nwins, a.workers, a.timeout)
    poswin = {f"{c}:{s}-{e}" for c, s, e in wins}
    d = X.add_phylop(d_all[d_all.win.isin(poswin)].copy())
    neg_rows = d_all[d_all.ID.isin(set(n.ID))].drop_duplicates("ID").copy()
    import comparators as CP
    blocks = []  # phyloP for negative positions: one UCSC request per merged block (gap <= 20 kb, span <= 50 kb)
    for c, x in neg_rows.groupby("chrom"):
        ps = np.sort(x.pos.unique())
        s = e = ps[0]
        for q in ps[1:]:
            if q - e <= 20000 and q - s <= 50000:
                e = q
            else:
                blocks.append((c, s, e)); s = e = q
        blocks.append((c, s, e))
    vals = {}
    for k, (c, s, e) in enumerate(blocks, 1):
        try:
            vals.update({(c, p): v for p, v in CP.ucsc_track("phyloP447way", c, s, e).items()})
        except Exception:
            pass
        time.sleep(0.25)
        if k % 100 == 0 or k == len(blocks):
            log(f"  phyloP negatives {k}/{len(blocks)} blocks")
    neg_rows["phylop447"] = [vals.get((c, p), np.nan) for c, p in zip(neg_rows.chrom, neg_rows.pos)]
    neg_rows["win"] = "negative"
    t = pd.concat([d, neg_rows], ignore_index=True).drop_duplicates("ID")
    ref_ok = set(t.ID)
    info = {"positives_lifted": len(p), "negatives_lifted": len(n), "positives_scored": int(p.ID.isin(ref_ok).sum()),
            "negatives_scored": int(n.ID.isin(ref_ok).sum()), "skipped_windows": len(skipped)}
    p.to_csv(os.path.join(OUT, "positives_GRCh38.tsv"), sep="\t", index=False)
    n.to_csv(os.path.join(OUT, "negatives_GRCh38.tsv"), sep="\t", index=False)
    t.to_csv(os.path.join(OUT, "window_snvs.tsv.gz"), sep="\t", index=False, float_format="%.5g", compression="gzip")
    json.dump(info, open(os.path.join(OUT, "fetch_info.json"), "w"), indent=1)
    log(json.dumps(info))


def cmd_eval(a):
    rng = np.random.default_rng(SEED)
    d = pd.read_csv(os.path.join(OUT, "window_snvs.tsv.gz"), sep="\t").drop_duplicates("ID").reset_index(drop=True)
    p = pd.read_csv(os.path.join(OUT, "positives_GRCh38.tsv"), sep="\t")
    n = pd.read_csv(os.path.join(OUT, "negatives_GRCh38.tsv"), sep="\t")
    P, N = set(p.ID) & set(d.ID), set(n.ID) & set(d.ID)
    common = (d["filter"] == "PASS") & (d.AC >= 9)
    pool_flag = (~d.ID.isin(P | N) & ~common & (d.win != "negative")).to_numpy()
    by_chrom = {c: x for c, x in d[pool_flag].groupby("chrom")}
    pos_idx = d.index[d.ID.isin(P)].to_numpy()
    neg_idx = d.index[d.ID.isin(N)].to_numpy()
    pools, keep = [], []
    for i in pos_idx:
        x = by_chrom.get(d.chrom[i])
        w = x.index[(x.pos >= d.pos[i] - A.FLANK) & (x.pos <= d.pos[i] + A.FLANK)].to_numpy() if x is not None else []
        if len(w):
            pools.append(w); keep.append(i)
    pos_idx = np.array(keep)
    pool_all = np.unique(np.concatenate(pools))
    y = np.r_[np.ones(len(pos_idx), bool), np.zeros(len(neg_idx), bool)]
    rep = [f"A11: NCBoost published benchmark re-scored (GRCh38): {len(pos_idx)} pathogenic vs {len(neg_idx)} common "
           f"(published negatives); pool {len(pool_all)} rare/absent SNVs in +-100 bp of positives", ""]
    rows = []
    for col, lab in TOOLS:
        s = d[col].to_numpy(float)
        cover = float(np.isfinite(s[np.r_[pos_idx, neg_idx]]).mean())
        s = np.where(np.isnan(s), np.nanmedian(s), s)
        raw = S.auroc(y, s[np.r_[pos_idx, neg_idx]])
        lo, hi = S.boot_auroc(y, s[np.r_[pos_idx, neg_idx]], a.boot)
        nulls = np.array([S.auroc(y, s[np.r_[np.array([rng.choice(pl) for pl in pools]), neg_idx]]) for _ in range(A.NULL_DRAWS)])
        nm = float(nulls.mean())
        matched = S.auroc(np.r_[np.ones(len(pos_idx), bool), np.zeros(len(pool_all), bool)], s[np.r_[pos_idx, pool_all]])
        rows.append({"tool": lab, "n_P": len(pos_idx), "n_N": len(neg_idx), "coverage": cover, "raw_auroc": raw, "raw_lo": lo,
                     "raw_hi": hi, "null_mean": nm, "inflation": nm - 0.5, "corrected": raw - nm, "rarity_matched_auroc": matched})
        rep.append(f"  {lab:11s} raw {raw:.3f} [{lo:.3f}, {hi:.3f}]  null {nm:.3f}  inflation {nm - 0.5:+.3f}  corrected "
                   f"{raw - nm:+.3f}  rarity-matched {matched:.3f}  (coverage {cover:.2f})")
    t = pd.DataFrame(rows).set_index("tool")
    h1 = int((t.inflation >= 0.05).sum()) >= 2
    rep += ["", f"H1 (inflation >= 0.05 for >= 2 tools on the published benchmark): {'MET' if h1 else 'NOT MET'}"]
    for m in ("raw_auroc", "corrected", "rarity_matched_auroc"):
        rep.append(f"H2 ({m}): ReMM {'>' if t.loc['ReMM v0.4', m] > t.loc['CADD v1.7', m] else '<='} CADD "
                   f"({t.loc['ReMM v0.4', m]:.3f} vs {t.loc['CADD v1.7', m]:.3f}); order: "
                   + " > ".join(t[m].sort_values(ascending=False).index))
    t.reset_index().to_csv(os.path.join(RES, "a11_results.tsv"), sep="\t", index=False, float_format="%.4f")
    open(os.path.join(RES, "a11_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["fetch", "eval"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    {"fetch": cmd_fetch, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
