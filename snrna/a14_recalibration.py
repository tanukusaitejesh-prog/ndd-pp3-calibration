#!/usr/bin/env python3
"""Amendment A14: do ClinGen PP3/BP4 thresholds survive when the benign set is as rare as the variants clinicians see?
(PREREG s15, A14). Written before any A14 variant was selected or scored.

A10 showed that ClinVar B/LB negatives do not remove the rarity-only inflation (AlphaMissense 0.916 against B/LB vs 0.920
against gnomAD-common SNVs): most ClinVar benign classifications rest on population frequency (BA1/BS1). ClinGen's
PP3/BP4 calibrations (Pejaver 2022; Bergquist 2025) use such B/LB sets. In practice PP3 is applied to variants that are
already rare (PM2), so the evidence a score carries is its separation of pathogenic from RARE benign variants.
A14 recalibrates missense predictors with B/LB restricted to rare variants and counts the evidence that changes.

    python a14_recalibration.py select            # ClinVar missense P/LP + B/LB (>= 1 star), capped -> a14/selected.tsv (commit)
    python a14_recalibration.py fetch             # gnomAD v4.1 AC, CADD, GPN-MSA (remote tabix) + AlphaMissense (+ REVEL)
    python a14_recalibration.py fetch --revel revel_with_transcript_ids.csv   # optional REVEL v1.3 file (GRCh38 column used)
    python a14_recalibration.py eval              # -> results/a14_{report.txt,thresholds.tsv,published.tsv,curves.tsv}

Selection: ClinVar GRCh38 VCF (fileDate recorded), SNVs whose consequences include missense_variant and no splice,
nonsense, stop or start term (A10 rule), P/LP or B/LB with >= 1 review star; at most 30 per gene per class (so a few
genes cannot dominate) and at most 6,000 per class (seed 20261023).
Sets (gnomAD v4.1 genomes, PASS AC): RARE = absent or AC <= 5 (the "absent or extremely rare" range where PM2 applies).
  ALL:  P/LP vs all B/LB                (the ClinGen calibration design)
  RARE: rare P/LP vs rare B/LB          (the clinically relevant contrast)
Method: calibrate.py's local-LR procedure unchanged (nearest 20%, >= 10 per class, 1,000 bootstraps; PP3 at strength k
where the 5th percentile local LR clears the threshold at every higher score; BP4 mirrored). Prior 0.0441 (Pejaver 2022)
primary, 0.10 secondary.
Tools: AlphaMissense (primary), CADD v1.7, GPN-MSA, REVEL if supplied. Published PP3 thresholds checked:
  AlphaMissense >= 0.792 / 0.906 / 0.972 (Supporting / Moderate / Strong; Bergquist et al. 2025)
  REVEL >= 0.644 / 0.773 / 0.932 (Pejaver et al. 2022); CADD >= 25.3 / 28.1 (Supporting / Moderate; Pejaver et al. 2022)
H1 (primary): at prior 0.0441, AlphaMissense's published PP3 Strong threshold (0.972) does NOT reach Strong in the RARE
  set (the 5th-percentile local LR falls below the Strong LR somewhere at or above 0.972), while the ALL set is reported
  alongside. Either outcome is reported.
H2 (descriptive): for every tool x strength, the threshold under ALL vs RARE (stricter, unchanged, or unreachable), and
  the share of rare P/LP variants whose PP3 strength drops when RARE thresholds replace ALL thresholds.
Validity check (reported first): under ALL, our AlphaMissense thresholds are compared with the published ones.
Caveat stated in advance: rare B/LB were classified without frequency evidence and may rest partly on BP4 (low
in-silico scores), which makes them EASIER to separate; any shift towards stricter thresholds is therefore conservative.
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

OUT = os.path.join(HERE, "a14")
RES = os.path.join(HERE, "results")
SEED = 20261023
CAP, PER_GENE, RARE_AC = 6000, 30, 5
PRIORS = (0.0441, 0.10)
STRENGTHS = ("Supporting", "Moderate", "Strong")
TOOLS = [("alphamissense", "AlphaMissense"), ("revel", "REVEL"), ("cadd_phred", "CADD v1.7"), ("gpn_msa", "GPN-MSA")]
PUBLISHED = {"alphamissense": {"Supporting": 0.792, "Moderate": 0.906, "Strong": 0.972},
             "revel": {"Supporting": 0.644, "Moderate": 0.773, "Strong": 0.932},
             "cadd_phred": {"Supporting": 25.3, "Moderate": 28.1}}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def cap_select(df, seed=SEED):
    """At most PER_GENE per gene, then at most CAP overall (tested)."""
    df = df.sample(frac=1, random_state=seed)
    df = df.groupby(df.gene.fillna(""), sort=False).head(PER_GENE)
    return df.head(CAP)


def cmd_select(a):
    import a10_extend as X
    import snrna as S
    os.makedirs(OUT, exist_ok=True)
    vcf = a.vcf or os.path.join(S.DATA, "clinvar_GRCh38.vcf.gz")
    if not os.path.exists(vcf):
        sys.exit(f"{vcf} missing: pass --vcf (e.g. the A12 download data/a12_clinvar_<date>.vcf.gz)")
    import gzip
    fdate = next((l.strip().split("=")[1] for l in gzip.open(vcf, "rt") if l.startswith("##fileDate")), "")
    cv = pd.DataFrame(list(X.clinvar_rows(vcf)))
    cv = cv[cv.chrom.isin([f"chr{i}" for i in list(range(1, 23)) + ["X"]])]
    mis = cv.mc.map(lambda m: "missense_variant" in m and not (m & X.MISSENSE_EXCLUDE))
    cv = cv[mis].drop_duplicates("ID")
    plp, blb = cv[cv.clnsig.isin(X.A.PLP)], cv[cv.clnsig.isin(X.BLB)]
    sel = pd.concat([cap_select(plp).assign(cls="P/LP"), cap_select(blb).assign(cls="B/LB")])
    sel = sel.assign(mc=sel.mc.map(lambda m: ",".join(sorted(m))))
    sel.to_csv(os.path.join(OUT, "selected.tsv"), sep="\t", index=False, lineterminator="\n")
    info = {"clinvar_fileDate": fdate, "plp_available": len(plp), "blb_available": len(blb),
            "kept": sel.cls.value_counts().to_dict(), "genes": int(sel.gene.nunique()),
            "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(info, open(os.path.join(OUT, "select_info.json"), "w"), indent=1)
    log(json.dumps(info))
    print("NEXT: git add a14/selected.tsv a14/select_info.json && git commit (before fetch)")


def site_windows(sel, gap=500):
    wins = []
    for c, x in sel.groupby("chrom"):
        ps = np.sort(x.pos.unique())
        s = e = ps[0]
        for p in ps[1:]:
            if p - e <= gap:
                e = p
            else:
                wins.append((c, int(s), int(e))); s = e = p
        wins.append((c, int(s), int(e)))
    return wins


def add_revel(d, path):
    want = set(d.ID)
    rv = {}
    for ch in pd.read_csv(path, usecols=["chr", "grch38_pos", "ref", "alt", "REVEL"], dtype=str, chunksize=2_000_000):
        ch = ch[ch.grch38_pos != "."]
        ids = "chr" + ch.chr + "-" + ch.grch38_pos + "-" + ch.ref + "-" + ch.alt
        m = ids.isin(want)
        for k, v in zip(ids[m], ch.REVEL[m].astype(float)):
            rv[k] = max(rv.get(k, v), v)
    d["revel"] = d.ID.map(rv)
    return d


def cmd_fetch(a):
    import a10_extend as X
    sel = pd.read_csv(os.path.join(OUT, "selected.tsv"), sep="\t")
    wins = site_windows(sel)
    log(f"{len(sel)} selected SNVs in {len(wins)} fetch windows")
    d, skipped = X.fetch_windows(wins, a.workers, a.timeout)
    d = d[d.ID.isin(set(sel.ID))].drop_duplicates("ID")
    d = X.add_alphamissense(d)
    if a.revel:
        d = add_revel(d, a.revel)
    d = sel[["ID", "gene", "cls"]].merge(d, on="ID", how="left")
    d.to_csv(os.path.join(OUT, "scores.tsv.gz"), sep="\t", index=False, float_format="%.5g", compression="gzip")
    json.dump({"windows": len(wins), "skipped": len(skipped), "scored": int(d.cadd_phred.notna().sum())},
              open(os.path.join(OUT, "fetch_info.json"), "w"), indent=1)
    log("coverage " + ", ".join(f"{c} {d[c].notna().mean():.2f}" for c, _ in TOOLS if c in d))


def define_sets(d):
    """ALL and RARE label sets (tested). AC is gnomAD v4.1 genome PASS allele count (NaN/absent -> 0)."""
    fetched = d.cadd_phred.notna().to_numpy()          # CADD covers every SNV, so NaN means the site was not fetched
    ac = d.AC.fillna(0).where((d["filter"].fillna("PASS") == "PASS") | (d.AC.fillna(0) == 0), np.inf)
    rare = (ac <= RARE_AC).to_numpy()
    y = (d.cls == "P/LP").to_numpy()
    return {"ALL": (fetched, y), "RARE": (fetched & rare, y)}


def suffix_met(grid, lo5, need, t):
    """Is the strength met (5th-pct local LR >= need) at every grid score >= t?"""
    m = grid >= t
    return bool(m.any() and (lo5[m] >= need).all())


def cmd_eval(a):
    import calibrate as C
    import snrna as S
    d = pd.read_csv(os.path.join(OUT, "scores.tsv.gz"), sep="\t")
    sets = define_sets(d)
    rows, pub, curves = [], [], []
    rep = ["A14: ClinGen PP3/BP4 recalibration with rare benign variants (missense)", ""]
    for sname, (mask, y_all) in sets.items():
        rep.append(f"{sname}: {int((mask & y_all).sum())} P/LP vs {int((mask & ~y_all).sum())} B/LB")
    rep.append("")
    store = {}
    for col, lab in TOOLS:
        if col not in d or d[col].notna().sum() < 50:
            continue
        for sname, (mask, y_all) in sets.items():
            m = mask & d[col].notna().to_numpy()
            s, y = d[col].to_numpy(float)[m], y_all[m]
            if y.sum() < 10 or (~y).sum() < 10:
                continue
            grid, point, lo5, hi95 = C.lr_bootstrap(s, y, a.boot, seed=SEED)
            store[(col, sname)] = (grid, lo5, s, y)
            curves += [{"tool": lab, "set": sname, "s": g, "lr": p, "lr_lo5": l, "lr_hi95": h} for g, p, l, h in zip(grid, point, lo5, hi95)]
            for prior in PRIORS:
                th = C.thresholds(grid, lo5, hi95, prior)
                rows.append({"tool": lab, "set": sname, "prior": prior, "n_P": int(y.sum()), "n_B": int((~y).sum()),
                             "auroc": S.auroc(y, s), **th})
                for k, t in PUBLISHED.get(col, {}).items():
                    need = S.lr_thresholds(prior)[k]
                    above = s >= t
                    pub.append({"tool": lab, "set": sname, "prior": prior, "strength": k, "published_threshold": t,
                                "met_at_all_scores_above": suffix_met(grid, lo5, need, t), "lr_needed": need,
                                "min_lr_lo5_above": float(lo5[grid >= t].min()) if (grid >= t).any() else np.nan,
                                "share_P_above": float(above[y].mean()), "share_B_above": float(above[~y].mean())})
    T, Pb = pd.DataFrame(rows), pd.DataFrame(pub)
    for prior in PRIORS:
        rep.append(f"--- prior {prior}")
        for lab in T.tool.unique():
            for sname in ("ALL", "RARE"):
                r = T[(T.tool == lab) & (T.set == sname) & (T.prior == prior)]
                if len(r):
                    r = r.iloc[0]
                    f = lambda v: "none" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.3g}"
                    rep.append(f"  {lab:13s} {sname:4s} AUROC {r.auroc:.3f} | PP3 Sup {f(r.PP3_Supporting)} Mod {f(r.PP3_Moderate)} "
                               f"Str {f(r.PP3_Strong)} | BP4 Sup {f(r.BP4_Supporting)} Mod {f(r.BP4_Moderate)} Str {f(r.BP4_Strong)}"
                               f"  ({r.n_P} P / {r.n_B} B)")
        rep.append("")
    rep.append("Published thresholds: is the strength met at every score above the threshold? (5th-pct local LR)")
    for r in Pb[Pb.prior == PRIORS[0]].itertuples():
        rep.append(f"  {r.tool:13s} {r.set:4s} {r.strength:10s} >= {r.published_threshold}: {'MET' if r.met_at_all_scores_above else 'NOT met'}"
                   f" (min LR lo5 {r.min_lr_lo5_above:.2f}, needed {r.lr_needed:.2f}; P above {r.share_P_above:.2f}, B above {r.share_B_above:.2f})")
    am = Pb[(Pb.tool == "AlphaMissense") & (Pb.strength == "Strong") & (Pb.prior == PRIORS[0])].set_index("set")
    if "RARE" in am.index:
        h1 = not bool(am.loc["RARE", "met_at_all_scores_above"])
        rep += ["", f"H1 (AlphaMissense published PP3 Strong >= 0.972 NOT met in the RARE set, prior 0.0441): "
                    f"{'SUPPORTED' if h1 else 'NOT SUPPORTED'} (ALL set: {'met' if am.loc['ALL', 'met_at_all_scores_above'] else 'not met'})"]
    # H2: strength lost by rare P/LP when RARE thresholds replace ALL thresholds
    rep.append("")
    for lab in T.tool.unique():
        a_ = T[(T.tool == lab) & (T.set == "ALL") & (T.prior == PRIORS[0])]
        r_ = T[(T.tool == lab) & (T.set == "RARE") & (T.prior == PRIORS[0])]
        col = dict((l, c) for c, l in TOOLS)[lab]
        if not len(a_) or not len(r_) or (col, "RARE") not in store:
            continue
        s, y = store[(col, "RARE")][2], store[(col, "RARE")][3]
        def level(v, row):
            lv = np.zeros(len(v), int)
            for i, k in enumerate(STRENGTHS, 1):
                t = row.iloc[0][f"PP3_{k}"]
                if t is not None and np.isfinite(t):
                    lv = np.where(v >= t, i, lv)
            return lv
        la, lr = level(s[y], a_), level(s[y], r_)
        rep.append(f"H2 {lab}: rare P/LP whose PP3 strength drops under RARE thresholds: {np.mean(lr < la):.1%} "
                   f"(rises {np.mean(lr > la):.1%}; n = {int(y.sum())})")
    os.makedirs(RES, exist_ok=True)
    T.to_csv(os.path.join(RES, "a14_thresholds.tsv"), sep="\t", index=False, float_format="%.4g", lineterminator="\n")
    Pb.to_csv(os.path.join(RES, "a14_published.tsv"), sep="\t", index=False, float_format="%.4g", lineterminator="\n")
    pd.DataFrame(curves).to_csv(os.path.join(RES, "a14_curves.tsv.gz"), sep="\t", index=False, float_format="%.4g", compression="gzip")
    open(os.path.join(RES, "a14_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["select", "fetch", "eval"])
    ap.add_argument("--vcf", help="select: ClinVar GRCh38 VCF (default data/clinvar_GRCh38.vcf.gz)")
    ap.add_argument("--revel", help="fetch: REVEL v1.3 CSV (columns chr, grch38_pos, ref, alt, REVEL)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    {"select": cmd_select, "fetch": cmd_fetch, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
