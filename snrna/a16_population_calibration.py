#!/usr/bin/env python3
"""Amendment A16: PP3 evidence against the rare-variant population that clinicians actually see (PREREG s15, A16).
Written after A14 (results seen) and before any A16 window was sampled or scored.

A14 found that curated rare ClinVar B/LB variants are far easier to separate from P/LP (AlphaMissense AUROC 0.967) than
random rare variants are (A10/A15: 0.694; deattenuated 0.70-0.74). Clinical PP3 is applied to a rare variant drawn
from that population, not from ClinVar's benign set. A16 measures the evidence each published PP3 band carries when
the comparison set is random rare missense variants observed in healthy people (gnomAD v4.1 genomes, AC 1-5), which
is an unlabelled set containing an unknown fraction pi of damaging variants (positive-unlabelled calibration).

    python a16_population_calibration.py sample            # 3,000 random 1-kb coding windows from the AlphaMissense file -> a16/windows.tsv (commit)
    python a16_population_calibration.py fetch --vcf data/a12_clinvar_<date>.vcf.gz [--revel <REVEL csv>]   # gnomAD + CADD + GPN-MSA + AM
    python a16_population_calibration.py eval              # -> results/a16_{report.txt,bands.tsv,thresholds.tsv,curves.tsv.gz}

Sets. P: rare ClinVar P/LP missense from A14 (gnomAD genome AC <= 5). U: every missense SNV observed in gnomAD v4.1
genomes (PASS, AC 1-5) in the sampled windows, minus any ClinVar P/LP, B/LB or VUS. B_cur: A14's curated B/LB (ALL and
RARE), for the side-by-side contrast.
Positive-unlabelled correction. With U = (1 - pi) B + pi P, the observed ratio LR_obs = f_P / f_U and the true
LR = (1 - pi) / (1 / LR_obs - pi) (infinite when 1 / LR_obs <= pi). pi is not identified; its upper bound is
pi_max = min over scores of 1 / LR_obs (reported with a bootstrap interval). pi grid {0, 0.01, 0.02, 0.05}; pi = 0
treats U as benign (the most conservative case), larger pi raises the LR.
Band LRs (primary readout, robust to sparse high-score tails). For each published band of a tool (AlphaMissense
Supporting [0.792, 0.906), Moderate [0.906, 0.972), Strong >= 0.972; REVEL 0.644 / 0.773 / 0.932; CADD 25.3 / 28.1),
LR_band = P(band | P) / P(band | comparison set), 0.5 pseudocounts, 1,000 bootstrap replicates, corrected for pi.
Verdict per band, three-way: DEMONSTRABLY MET if the 5th percentile >= the ClinGen LR for that strength (prior 0.0441),
DEMONSTRABLY NOT MET if the 95th percentile < it, otherwise UNDETERMINED (precision; A14's lesson).
H1 (primary): AlphaMissense's published Strong band (>= 0.972) is DEMONSTRABLY NOT MET against U at pi = 0.02 (the
assumed share of observed rare missense variants that would be classified P/LP); pi = 0 (U benign) and pi = 0.05
(stress test) are reported alongside. The correction grows steeply as pi approaches pi_max, so verdicts near pi_max
are not interpreted.
Break-even pi (assumption-free readout): for each band, the pi at which its corrected LR would reach the ClinGen LR,
pi* = (need / LR_obs - 1) / (need - 1), computed from the point LR and from the 95th-percentile LR (the most favourable);
if pi* exceeds pi_max, the band cannot reach that strength for any pi the data allow.
H2 (descriptive): the same verdicts for every tool x band x pi, against U and against curated B/LB, and the
population false-PP3 burden (share of U variants in each band, which needs no pi).
Also: local-LR thresholds (calibrate.py procedure) against U at each pi, as A14.
"""
import argparse
import gzip
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

OUT = os.path.join(HERE, "a16")
RES = os.path.join(HERE, "results")
SEED = 20261025
N_WIN, HALF = 3000, 500
RARE_AC = 5
PIS = (0.0, 0.01, 0.02, 0.05)
PRIOR = 0.0441
BANDS = {"alphamissense": [("Supporting", 0.792, 0.906), ("Moderate", 0.906, 0.972), ("Strong", 0.972, np.inf)],
         "revel": [("Supporting", 0.644, 0.773), ("Moderate", 0.773, 0.932), ("Strong", 0.932, np.inf)],
         "cadd_phred": [("Supporting", 25.3, 28.1), ("Moderate", 28.1, np.inf)]}
LABELS = {"alphamissense": "AlphaMissense", "revel": "REVEL", "cadd_phred": "CADD v1.7", "gpn_msa": "GPN-MSA"}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ---------------------------------------------------------------------------------------------------- sample

def reservoir(lines, k, rng):
    """Uniform sample of k items from an iterable of unknown length (tested)."""
    res = []
    for i, x in enumerate(lines):
        if i < k:
            res.append(x)
        else:
            j = rng.integers(0, i + 1)
            if j < k:
                res[j] = x
    return res


def cmd_sample(a):
    import snrna as S
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(S.DATA, "AlphaMissense_hg38.tsv.gz")
    if not os.path.exists(path):
        sys.exit(f"{path} missing (a10_extend.py missense_fetch downloads it)")
    rng = np.random.default_rng(SEED)
    log("streaming the AlphaMissense file (~70 M lines, several minutes)")
    with gzip.open(path, "rt") as f:
        pick = reservoir((l.split("\t", 2)[:2] for l in f if not l.startswith("#")), N_WIN, rng)
    w = pd.DataFrame(pick, columns=["chrom", "pos"]).assign(pos=lambda d: d.pos.astype(int))
    w = w[w.chrom.isin([f"chr{i}" for i in list(range(1, 23)) + ["X"]])]
    w = w.assign(start=w.pos - HALF, end=w.pos + HALF).sort_values(["chrom", "start"]).drop_duplicates(["chrom", "start"])
    w[["chrom", "start", "end"]].to_csv(os.path.join(OUT, "windows.tsv"), sep="\t", index=False, lineterminator="\n")
    json.dump({"windows": len(w), "seed": SEED, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
              open(os.path.join(OUT, "sample_info.json"), "w"), indent=1)
    log(f"{len(w)} windows written")
    print("NEXT: git add a16/windows.tsv a16/sample_info.json && git commit (before fetch)")


# ----------------------------------------------------------------------------------------------------- fetch

def clinvar_labelled(vcf):
    import a10_extend as X
    ids = set()
    for r in X.clinvar_rows(vcf):
        if r["clnsig"]:
            ids.add(r["ID"])
    return ids


def cmd_fetch(a):
    import a10_extend as X
    w = pd.read_csv(os.path.join(OUT, "windows.tsv"), sep="\t")
    wins = [(r.chrom, int(r.start), int(r.end)) for r in w.itertuples()]
    log(f"{len(wins)} windows; {a.workers} workers")
    d, skipped = X.fetch_windows(wins, a.workers, a.timeout)
    d = d.drop_duplicates("ID")
    obs = d[(d["filter"] == "PASS") & (d.AC >= 1) & (d.AC <= RARE_AC)].copy()
    log(f"{len(obs)} observed rare SNVs (PASS, AC 1-{RARE_AC}) before the missense filter")
    obs = X.add_alphamissense(obs)
    obs = obs[obs.alphamissense.notna()]
    n0 = len(obs)
    if a.vcf:
        lab = clinvar_labelled(a.vcf)
        obs = obs[~obs.ID.isin(lab)]
    if a.revel:
        import a14_recalibration as A14
        obs = A14.add_revel(obs, a.revel)
    cols = [c for c in ["ID", "chrom", "pos", "AC", "cadd_phred", "gpn_msa", "alphamissense", "revel"] if c in obs]
    obs[cols].to_csv(os.path.join(OUT, "unlabelled.tsv.gz"), sep="\t", index=False, float_format="%.5g", compression="gzip")
    info = {"windows": len(wins), "skipped": len(skipped), "observed_rare_missense": n0, "after_clinvar_exclusion": len(obs),
            "clinvar_excluded": bool(a.vcf)}
    json.dump(info, open(os.path.join(OUT, "fetch_info.json"), "w"), indent=1)
    log(json.dumps(info))


# ------------------------------------------------------------------------------------------------------ eval

def pu_correct(lr_obs, pi):
    """True LR from LR_obs = f_P / f_U when U = (1 - pi) B + pi P (tested). inf where 1 / LR_obs <= pi."""
    lr_obs = np.asarray(lr_obs, float)
    if pi == 0:
        return lr_obs
    inv = 1 / np.maximum(lr_obs, 1e-12) - pi
    return np.where(inv > 0, (1 - pi) / np.where(inv > 0, inv, 1), np.inf)


def band_lr(sp, sn, lo, hi):
    """P(band | P) / P(band | N) with 0.5 pseudocounts."""
    a = ((sp >= lo) & (sp < hi)).sum() + 0.5
    b = ((sn >= lo) & (sn < hi)).sum() + 0.5
    return (a / (len(sp) + 1)) / (b / (len(sn) + 1))


def band_table(sp, sn, bands, pis, B, rng, need, comparison, tool):
    rows = []
    idx_p, idx_n = np.arange(len(sp)), np.arange(len(sn))
    boots = {b: [] for b, _, _ in bands}
    for _ in range(B):
        bp, bn = sp[rng.choice(idx_p, len(idx_p))], sn[rng.choice(idx_n, len(idx_n))]
        for b, lo, hi in bands:
            boots[b].append(band_lr(bp, bn, lo, hi))
    for b, lo, hi in bands:
        pt, bs = band_lr(sp, sn, lo, hi), np.array(boots[b])
        for pi in pis:
            c_pt, c_bs = float(pu_correct(pt, pi)), pu_correct(bs, pi)
            q5, q95 = np.percentile(c_bs, [5, 95])
            n = need[b]
            verdict = "DEMONSTRABLY MET" if q5 >= n else "DEMONSTRABLY NOT MET" if q95 < n else "UNDETERMINED"
            be = lambda L: float((n / L - 1) / (n - 1)) if L > 0 else np.inf
            rows.append({"tool": tool, "comparison": comparison, "band": b, "lo": lo, "hi": hi, "pi": pi, "lr": c_pt, "lr_q5": q5,
                         "break_even_pi": be(pt), "break_even_pi_from_q95": be(float(np.percentile(bs, 95))),
                         "lr_q95": q95, "lr_needed": n, "verdict": verdict, "share_P_in_band": float(((sp >= lo) & (sp < hi)).mean()),
                         "share_comparison_in_band": float(((sn >= lo) & (sn < hi)).mean()), "n_P": len(sp), "n_comparison": len(sn)})
    return rows


def pi_max(sp, su, B, rng):
    """Upper bound on pi: min over the score grid of 1 / LR_obs (local LR, calibrate.py's estimator)."""
    import snrna as S
    s, y = np.r_[sp, su], np.r_[np.ones(len(sp), bool), np.zeros(len(su), bool)]
    grid = np.quantile(s, np.linspace(0.02, 0.98, 49))
    pt = float(np.min(1 / np.maximum(S.local_lr(s, y, grid), 1e-12)))
    bs = []
    for _ in range(B):
        i = np.r_[rng.choice(np.flatnonzero(y), y.sum()), rng.choice(np.flatnonzero(~y), (~y).sum())]
        bs.append(float(np.min(1 / np.maximum(S.local_lr(s[i], y[i], grid), 1e-12))))
    return pt, *np.percentile(bs, [5, 95])


def cmd_eval(a):
    import a14_recalibration as A14
    import calibrate as C
    import snrna as S
    rng = np.random.default_rng(SEED)
    sc = pd.read_csv(os.path.join(HERE, "a14", "scores.tsv.gz"), sep="\t")
    sets = A14.define_sets(sc)
    rare_mask, y = sets["RARE"]
    all_mask, _ = sets["ALL"]
    P = sc[rare_mask & y]
    curated = {"curated B/LB (ALL)": sc[all_mask & ~y], "curated B/LB (RARE)": sc[rare_mask & ~y]}
    U = pd.read_csv(os.path.join(OUT, "unlabelled.tsv.gz"), sep="\t")
    need = {k: v for k, v in S.lr_thresholds(PRIOR).items()}
    rep = [f"A16: PP3 evidence against the rare-variant population (positive-unlabelled calibration)",
           f"P: {len(P)} rare ClinVar P/LP missense (A14); U: {len(U)} random observed rare missense (gnomAD AC 1-{RARE_AC}); "
           f"curated B/LB ALL {len(curated['curated B/LB (ALL)'])}, RARE {len(curated['curated B/LB (RARE)'])}", ""]
    bands_rows, th_rows, curves = [], [], []
    for col, lab in LABELS.items():
        if col not in U or col not in P or U[col].notna().sum() < 100:
            continue
        sp, su = P[col].dropna().to_numpy(float), U[col].dropna().to_numpy(float)
        au = S.auroc(np.r_[np.ones(len(sp), bool), np.zeros(len(su), bool)], np.r_[sp, su])
        pm = pi_max(sp, su, 200, rng)
        rep.append(f"{lab}: AUROC P vs U {au:.3f} (vs curated RARE B/LB "
                   f"{S.auroc(np.r_[np.ones(len(sp), bool), np.zeros(curated['curated B/LB (RARE)'][col].notna().sum(), bool)], np.r_[sp, curated['curated B/LB (RARE)'][col].dropna().to_numpy(float)]):.3f}); "
                   f"pi_max {pm[0]:.3f} [{pm[1]:.3f}, {pm[2]:.3f}]")
        if col in BANDS:
            bands_rows += band_table(sp, su, BANDS[col], PIS, a.boot, rng, need, "population U", lab)
            for cname, cd in curated.items():
                bands_rows += band_table(sp, cd[col].dropna().to_numpy(float), BANDS[col], (0.0,), a.boot, rng, need, cname, lab)
        s, yy = np.r_[sp, su], np.r_[np.ones(len(sp), bool), np.zeros(len(su), bool)]
        grid, point, lo5, hi95 = C.lr_bootstrap(s, yy, a.boot, seed=SEED)
        for pi in PIS:
            th = C.thresholds(grid, pu_correct(lo5, pi), pu_correct(hi95, pi), PRIOR)
            th_rows.append({"tool": lab, "pi": pi, "auroc": au, "pi_max": pm[0], **th})
        curves += [{"tool": lab, "s": g, "lr_obs": p, "lr_lo5": l, "lr_hi95": h} for g, p, l, h in zip(grid, point, lo5, hi95)]
    Bd, Th = pd.DataFrame(bands_rows), pd.DataFrame(th_rows)
    rep += ["", f"Published PP3 bands (prior {PRIOR}; LR [5th, 95th pct]; verdict):"]
    for r in Bd.itertuples():
        rep.append(f"  {r.tool:13s} {r.band:10s} vs {r.comparison:20s} pi {r.pi:<4} LR {r.lr:7.2f} [{r.lr_q5:7.2f}, {r.lr_q95:7.2f}] "
                   f"need {r.lr_needed:5.2f}  {r.verdict:21s} P in band {r.share_P_in_band:.2f}, comparison in band {r.share_comparison_in_band:.3f}")
    rep += ["", "Local-LR thresholds against U (calibrate.py procedure, pi-corrected):"]
    for r in Th.itertuples():
        f = lambda v: "none" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.3g}"
        rep.append(f"  {r.tool:13s} pi {r.pi:<4} PP3 Sup {f(r.PP3_Supporting)} Mod {f(r.PP3_Moderate)} Str {f(r.PP3_Strong)} | "
                   f"BP4 Sup {f(r.BP4_Supporting)} Mod {f(r.BP4_Moderate)} Str {f(r.BP4_Strong)}")
    rep += ["", "Break-even pi (damaging share of U needed for the band to reach its strength; compare with pi_max above):"]
    for r in Bd[(Bd.comparison == "population U") & (Bd.pi == 0.0)].itertuples():
        rep.append(f"  {r.tool:13s} {r.band:10s} pi* {r.break_even_pi:.3f} (most favourable {r.break_even_pi_from_q95:.3f})")
    h = Bd[(Bd.tool == "AlphaMissense") & (Bd.band == "Strong") & (Bd.comparison == "population U") & (Bd.pi == 0.02)]
    if len(h):
        rep += ["", f"H1 (AlphaMissense Strong band DEMONSTRABLY NOT MET against U at pi 0.02): "
                    f"{'SUPPORTED' if h.verdict.iloc[0] == 'DEMONSTRABLY NOT MET' else 'NOT SUPPORTED'} ({h.verdict.iloc[0]}; "
                    f"LR {h.lr.iloc[0]:.1f} [{h.lr_q5.iloc[0]:.1f}, {h.lr_q95.iloc[0]:.1f}], need {h.lr_needed.iloc[0]:.1f})"]
    os.makedirs(RES, exist_ok=True)
    Bd.to_csv(os.path.join(RES, "a16_bands.tsv"), sep="\t", index=False, float_format="%.4g", lineterminator="\n")
    Th.to_csv(os.path.join(RES, "a16_thresholds.tsv"), sep="\t", index=False, float_format="%.4g", lineterminator="\n")
    pd.DataFrame(curves).to_csv(os.path.join(RES, "a16_curves.tsv.gz"), sep="\t", index=False, float_format="%.4g", compression="gzip")
    open(os.path.join(RES, "a16_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["sample", "fetch", "eval"])
    ap.add_argument("--vcf", help="fetch: ClinVar VCF; its labelled variants are removed from U")
    ap.add_argument("--revel", help="fetch: REVEL v1.3 CSV")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    {"sample": cmd_sample, "fetch": cmd_fetch, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
