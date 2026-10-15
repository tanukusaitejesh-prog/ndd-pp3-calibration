#!/usr/bin/env python3
"""Amendment A13: the A11 re-analysis generalised to any published variant benchmark, starting with TraitGym
(PREREG s15, A13). Written before any A13 score was fetched.

A11 re-scored NCBoost and found its common-variant negatives inflate every tool. A13 applies the same pipeline to
other published benchmarks so the effect is shown on more than one benchmark design, and first checks whether the
benchmark's negatives are in fact common (the rarity-only null only bites when they are).

    python a13_benchmarks.py freeze traitgym_mendelian      # download + commit a13/<name>/benchmark.tsv (md5 recorded)
    python a13_benchmarks.py freeze mybench --tsv my.tsv --build GRCh37   # any TSV: chrom, pos, ref, alt, label (1/0)
    python a13_benchmarks.py fetch traitgym_mendelian       # DGX/any: CADD, gnomAD, GPN-MSA, ReMM in +-100 bp, phyloP447
    python a13_benchmarks.py eval traitgym_mendelian        # -> results/a13_<name>_{report.txt,results.tsv}

Presets:
  traitgym_mendelian  TraitGym (Benegas et al. 2025) Mendelian non-coding causal variants vs 9 matched controls each
                      (songlab/TraitGym, mendelian_traits_matched_9/test.parquet; GRCh38). Seen before this amendment:
                      column names and class counts only (338 positives, 3,042 controls); no score of any variant.
Metrics, as A11: raw AUROC (published positives vs published negatives), control commonness (share of negatives that
are gnomAD v4.1 genome PASS AC >= 9, and median AF), window-matched rarity-only null (one SNV per positive from its own
+-100 bp window that is neither positive, negative nor gnomAD-common; 200 draws), inflation = null - 0.5,
corrected = raw - null, rarity-matched AUROC (positives vs every pool SNV in their windows) and, where the benchmark
has match groups, the within-group AUROC (pairs only within a match group, which is how the benchmark was built).
H1: inflation >= 0.05 for at least 2 of CADD, ReMM, GPN-MSA, phyloP447 - evaluated only if >= 50% of negatives are
gnomAD-common; otherwise H1 is "not applicable" and the null is reported descriptively.
H2 (descriptive): whether the tool order changes between raw, corrected and rarity-matched AUROC (ties within
overlapping intervals are reported as ties).
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

OUT = os.path.join(HERE, "a13")
RES = os.path.join(HERE, "results")
SEED = 20261022
NULL_DRAWS, FLANK = 200, 100
TOOLS = [("cadd_phred", "CADD v1.7"), ("remm", "ReMM v0.4"), ("gpn_msa", "GPN-MSA"), ("phylop447", "phyloP447")]
PRESETS = {
    "traitgym_mendelian": {"url": "https://huggingface.co/datasets/songlab/TraitGym/resolve/main/mendelian_traits_matched_9/test.parquet",
                           "build": "GRCh38", "group": "match_group",
                           "cite": "Benegas et al. 2025, TraitGym (bioRxiv 10.1101/2025.02.11.637758)"},
}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def md5(p):  # line endings normalised, so a Windows and a Linux checkout give the same hash
    return hashlib.md5(open(p, "rb").read().replace(b"\r\n", b"\n")).hexdigest()


def normalise(d, group=None):
    d = d.rename(columns={c: c.lower() for c in d.columns})
    d["chrom"] = "chr" + d.chrom.astype(str).str.replace("chr", "", regex=False)
    d["label"] = d.label.astype(str).str.lower().isin(["1", "true", "pathogenic", "positive", "yes"]).astype(int)
    keep = ["chrom", "pos", "ref", "alt", "label"] + ([group.lower()] if group and group.lower() in d else [])
    d = d[keep].rename(columns={group.lower(): "group"} if group else {})
    d = d[(d.ref.str.len() == 1) & (d.alt.str.len() == 1)].copy()
    d["pos"] = d.pos.astype(int)
    return d.drop_duplicates(["chrom", "pos", "ref", "alt"]).reset_index(drop=True)


def cmd_freeze(a):
    od = os.path.join(OUT, a.name)
    os.makedirs(od, exist_ok=True)
    pre = PRESETS.get(a.name, {})
    if a.tsv:
        raw, build, src = pd.read_csv(a.tsv, sep="\t"), a.build, os.path.abspath(a.tsv)
    elif pre:
        rp = os.path.join(od, "source" + os.path.splitext(pre["url"])[1])
        if not os.path.exists(rp):
            log(f"downloading {pre['url']}")
            urllib.request.urlretrieve(pre["url"], rp)
        raw = pd.read_parquet(rp) if rp.endswith(".parquet") else pd.read_csv(rp, sep="\t")
        build, src = pre["build"], pre["url"]
    else:
        sys.exit("unknown preset: give --tsv and --build")
    d = normalise(raw, pre.get("group") or a.group)
    if build == "GRCh37":
        from pyliftover import LiftOver
        lo, keep = LiftOver("hg19", "hg38"), []
        for r in d.itertuples():
            m = lo.convert_coordinate(r.chrom, r.pos - 1)
            keep.append(int(m[0][1]) + 1 if m and m[0][0] == r.chrom and m[0][2] == "+" else -1)
        d["pos37"], d["pos"] = d.pos, keep
        d = d[d.pos > 0]
    d["ID"] = d.chrom + "-" + d.pos.astype(str) + "-" + d.ref + "-" + d.alt
    bp = os.path.join(od, "benchmark.tsv")
    d.to_csv(bp, sep="\t", index=False, lineterminator="\n")
    info = {"name": a.name, "source": src, "build_in": build, "n_pos": int(d.label.sum()), "n_neg": int((d.label == 0).sum()),
            "has_groups": "group" in d, "benchmark_md5": md5(bp), "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "cite": pre.get("cite", "")}
    json.dump(info, open(os.path.join(od, "freeze_info.json"), "w"), indent=1)
    log(json.dumps(info))
    print(f"NEXT: git add a13/{a.name}/benchmark.tsv a13/{a.name}/freeze_info.json && git commit (before fetch)")


def cmd_fetch(a):
    import a9_inflation as A
    import a10_extend as X
    od = os.path.join(OUT, a.name)
    info = json.load(open(os.path.join(od, "freeze_info.json")))
    bp = os.path.join(od, "benchmark.tsv")
    if md5(bp) != info["benchmark_md5"]:
        sys.exit("benchmark.tsv changed since freeze")
    d = pd.read_csv(bp, sep="\t")
    p, n = d[d.label == 1], d[d.label == 0]
    wins = [(c, s, e) for c, x in p.groupby("chrom") for s, e in A.merged_windows(x.pos.tolist())]
    nwins = [(r.chrom, r.pos, r.pos) for r in n.drop_duplicates(["chrom", "pos"]).itertuples()]
    log(f"{len(p)} positives -> {len(wins)} windows (+-{FLANK} bp); {len(nwins)} negative positions")
    t, skipped = X.fetch_windows(wins + nwins, a.workers, a.timeout)
    poswin = {f"{c}:{s}-{e}" for c, s, e in wins}
    tp = X.add_phylop(t[t.win.isin(poswin)].copy())
    tn = t[t.ID.isin(set(n.ID)) & ~t.ID.isin(set(tp.ID))].drop_duplicates("ID").copy()
    if len(tn):
        import comparators as CP
        vals = {}
        for c, x in tn.groupby("chrom"):  # one UCSC request per merged block (gap <= 20 kb, span <= 50 kb), as A11
            ps = np.sort(x.pos.unique())
            s = e = ps[0]
            blocks = []
            for q in ps[1:]:
                if q - e <= 20000 and q - s <= 50000:
                    e = q
                else:
                    blocks.append((s, e)); s = e = q
            blocks.append((s, e))
            for s, e in blocks:
                try:
                    vals.update({(c, k): v for k, v in CP.ucsc_track("phyloP447way", c, s, e).items()})
                except Exception:
                    pass
                time.sleep(0.25)
        tn["phylop447"] = [vals.get((c, q), np.nan) for c, q in zip(tn.chrom, tn.pos)]
        tn["win"] = "negative"
    out = pd.concat([tp, tn], ignore_index=True).drop_duplicates("ID")
    out.to_csv(os.path.join(od, "window_snvs.tsv.gz"), sep="\t", index=False, float_format="%.5g", compression="gzip")
    fi = {"positives_scored": int(p.ID.isin(set(out.ID)).sum()), "negatives_scored": int(n.ID.isin(set(out.ID)).sum()),
          "skipped_windows": len(skipped)}
    json.dump(fi, open(os.path.join(od, "fetch_info.json"), "w"), indent=1)
    log(json.dumps(fi))


def evaluate_benchmark(d, bench, auroc, boot_auroc, boot, rng, tools=TOOLS):
    """Core metrics (tested on synthetic data). d: window SNV table (ID, chrom, pos, filter, AC, AF?, tool columns);
    bench: frozen benchmark (ID, label, optional group)."""
    d = d.drop_duplicates("ID").reset_index(drop=True)
    P = set(bench.ID[bench.label == 1]) & set(d.ID)
    N = set(bench.ID[bench.label == 0]) & set(d.ID)
    common = ((d["filter"] == "PASS") & (d.AC >= 9)).to_numpy()
    neg_idx = d.index[d.ID.isin(N)].to_numpy()
    frac_common = float(common[neg_idx].mean()) if len(neg_idx) else float("nan")
    af = d.AF.to_numpy(float)[neg_idx] if "AF" in d else np.array([np.nan])
    pool_flag = (~d.ID.isin(P | N) & ~common & (d.get("win", "") != "negative")).to_numpy()
    by_chrom = {c: x for c, x in d[pool_flag].groupby("chrom")}
    pools, keep = [], []
    for i in d.index[d.ID.isin(P)].to_numpy():
        x = by_chrom.get(d.chrom[i])
        w = x.index[(x.pos >= d.pos[i] - FLANK) & (x.pos <= d.pos[i] + FLANK)].to_numpy() if x is not None else []
        if len(w):
            pools.append(w); keep.append(i)
    pos_idx = np.array(keep, int)
    pool_all = np.unique(np.concatenate(pools)) if pools else np.array([], int)
    y = np.r_[np.ones(len(pos_idx), bool), np.zeros(len(neg_idx), bool)]
    grp = bench.set_index("ID").get("group")
    rows = []
    for col, lab in tools:
        if col not in d or not np.isfinite(d[col].to_numpy(float)).any():
            continue
        s = d[col].to_numpy(float)
        cover = float(np.isfinite(s[np.r_[pos_idx, neg_idx]]).mean())
        s = np.where(np.isnan(s), np.nanmedian(s), s)
        raw = auroc(y, s[np.r_[pos_idx, neg_idx]])
        lo, hi = boot_auroc(y, s[np.r_[pos_idx, neg_idx]], boot)
        nulls = np.array([auroc(y, s[np.r_[np.array([rng.choice(pl) for pl in pools], int), neg_idx]]) for _ in range(NULL_DRAWS)])
        nm = float(nulls.mean())
        matched = auroc(np.r_[np.ones(len(pos_idx), bool), np.zeros(len(pool_all), bool)], s[np.r_[pos_idx, pool_all]]) if len(pool_all) else np.nan
        within = np.nan
        if grp is not None:  # pairs only inside match groups
            g = d.ID.map(grp).to_numpy()
            num = den = 0.0
            for gv in pd.unique(g[pos_idx]):
                sp, sn = s[pos_idx[g[pos_idx] == gv]], s[neg_idx[g[neg_idx] == gv]]
                if len(sp) and len(sn):
                    num += ((sp[:, None] > sn[None, :]).sum() + 0.5 * (sp[:, None] == sn[None, :]).sum())
                    den += len(sp) * len(sn)
            within = num / den if den else np.nan
        rows.append({"tool": lab, "n_P": len(pos_idx), "n_N": len(neg_idx), "coverage": cover, "raw_auroc": raw, "raw_lo": lo,
                     "raw_hi": hi, "null_mean": nm, "null_lo": float(np.percentile(nulls, 2.5)), "null_hi": float(np.percentile(nulls, 97.5)),
                     "inflation": nm - 0.5, "corrected": raw - nm, "rarity_matched_auroc": matched, "within_group_auroc": within})
    return pd.DataFrame(rows), {"frac_neg_common": frac_common, "median_neg_af": float(np.nanmedian(af)) if np.isfinite(af).any() else None,
                                "n_pool": int(len(pool_all))}


def cmd_eval(a):
    import snrna as S
    od = os.path.join(OUT, a.name)
    for f in ("benchmark.tsv", "freeze_info.json"):
        tracked = subprocess.run(["git", "-C", HERE, "ls-files", "--error-unmatch", os.path.join(od, f)], capture_output=True).returncode == 0
        if not tracked:
            sys.exit(f"a13/{a.name}/{f} must be committed before eval")
    bench = pd.read_csv(os.path.join(od, "benchmark.tsv"), sep="\t")
    d = pd.read_csv(os.path.join(od, "window_snvs.tsv.gz"), sep="\t")
    t, q = evaluate_benchmark(d, bench, S.auroc, S.boot_auroc, a.boot, np.random.default_rng(SEED))
    applicable = q["frac_neg_common"] >= 0.5
    rep = [f"A13 {a.name}: {int(t.n_P.iloc[0]) if len(t) else 0} positives vs {int(t.n_N.iloc[0]) if len(t) else 0} published negatives; "
           f"pool {q['n_pool']} rare/absent SNVs", f"Negatives gnomAD-common (PASS AC >= 9): {q['frac_neg_common']:.2f}; "
           f"median AF {q['median_neg_af']}", ""]
    for r in t.itertuples():
        rep.append(f"  {r.tool:11s} raw {r.raw_auroc:.3f} [{r.raw_lo:.3f}, {r.raw_hi:.3f}]  null {r.null_mean:.3f}  inflation "
                   f"{r.inflation:+.3f}  corrected {r.corrected:+.3f}  rarity-matched {r.rarity_matched_auroc:.3f}  "
                   f"within-group {r.within_group_auroc:.3f}  (coverage {r.coverage:.2f})")
    h1 = int((t.inflation >= 0.05).sum()) >= 2
    rep += ["", f"H1: {'NOT APPLICABLE (negatives mostly not common)' if not applicable else 'MET' if h1 else 'NOT MET'}"]
    for m in ("raw_auroc", "corrected", "rarity_matched_auroc"):
        rep.append(f"H2 order by {m}: " + " > ".join(t.set_index("tool")[m].sort_values(ascending=False).index))
    os.makedirs(RES, exist_ok=True)
    t.assign(benchmark=a.name).to_csv(os.path.join(RES, f"a13_{a.name}_results.tsv"), sep="\t", index=False, float_format="%.4f")
    open(os.path.join(RES, f"a13_{a.name}_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["freeze", "fetch", "eval"])
    ap.add_argument("name", help="preset name or your own benchmark name")
    ap.add_argument("--tsv", help="freeze: your own benchmark (chrom, pos, ref, alt, label)")
    ap.add_argument("--build", default="GRCh38", choices=["GRCh37", "GRCh38"])
    ap.add_argument("--group", help="freeze: optional match-group column")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    {"freeze": cmd_freeze, "fetch": cmd_fetch, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
