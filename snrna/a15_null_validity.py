#!/usr/bin/env python3
"""Amendment A15: validity of the rarity-only null (ground-truth simulation) and confidence intervals for every
corrected and rarity-matched AUROC (PREREG s15, A15). Written before the simulation or any interval was computed.

Reviewers will object that rare variants are themselves enriched for damaging ones, so "inflation" is partly real
signal and "corrected AUROC" has no clear meaning. A15 answers with a simulation where the truth is known, and states
the estimand: in clinical use the variant is already rare (PM2), so the target is T = AUROC of pathogenic vs RARE
BENIGN variants. The rarity-matched AUROC M (pathogenic vs all rare variants) estimates T diluted by the damaging
fraction pi of the rare pool: M = (1 - pi) T + pi A_pd, A_pd ~ 0.5, so T ~ (M - pi/2) / (1 - pi) (deattenuated).

    python a15_null_validity.py simulate                 # CPU, no data -> results/a15_simulation.{tsv,txt}
    python a15_null_validity.py ci a9|a10|a11|a13:<name> # DGX: window bootstrap -> results/a15_ci_<x>.tsv

Simulation (pre-registered grid): windows each with 1 pathogenic positive, 60 rare SNVs (a fraction pi damaging) and
6 common SNVs (benign). Score = beta * damaging + gamma * rare + N(0, 1); beta in {0, 0.5, 1, 2} (pathogenicity
signal), gamma in {0, 0.5, 1, 2} (rarity/constraint tracking), pi in {0.05, 0.2}; 1,000 windows, 100 replicates.
Estimators: raw AUROC (positives vs common), corrected + 0.5, rarity-matched M, deattenuated M (pi known), against T.
Expectations fixed in advance: (E1) whenever gamma > 0, raw AUROC overstates T and |M - T| < |raw - T|; (E2) when
beta = 0, M and corrected + 0.5 are near 0.5 while raw is not; (E3) deattenuated M is within 0.02 of T on average.
CIs: cluster (window) bootstrap over positives (with their windows' negatives and pool), 200 replicates, 50 null
draws each, for raw, null mean, corrected and rarity-matched AUROC; and the paired difference in corrected AUROC
between every pair of tools (does a reordering exceed its interval?).
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
SEED = 20261024


def auroc(y, s):
    """Mann-Whitney AUROC with average ranks for ties."""
    y = np.asarray(y, bool)
    r = pd.Series(np.asarray(s, float)).rank().to_numpy()
    n1, n0 = y.sum(), (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)) if n1 and n0 else np.nan


def simulate_once(beta, gamma, pi, rng, n_win=1000, n_rare=60, n_common=6):
    pos = beta * 1 + gamma * 1 + rng.normal(size=n_win)
    dmg = rng.random((n_win, n_rare)) < pi
    rare = beta * dmg + gamma * 1 + rng.normal(size=(n_win, n_rare))
    common = rng.normal(size=(n_win, n_common))
    T = auroc(np.r_[np.ones(n_win, bool), np.zeros((~dmg).sum(), bool)], np.r_[pos, rare[~dmg]])
    raw = auroc(np.r_[np.ones(n_win, bool), np.zeros(common.size, bool)], np.r_[pos, common.ravel()])
    pick = rare[np.arange(n_win), rng.integers(0, n_rare, n_win)]       # one null draw per window (A9 design)
    null = auroc(np.r_[np.ones(n_win, bool), np.zeros(common.size, bool)], np.r_[pick, common.ravel()])
    M = auroc(np.r_[np.ones(n_win, bool), np.zeros(rare.size, bool)], np.r_[pos, rare.ravel()])
    return {"T": T, "raw": raw, "null": null, "corrected_plus_half": raw - null + 0.5, "matched": M,
            "matched_deattenuated": (M - pi / 2) / (1 - pi)}


def cmd_simulate(a):
    rng = np.random.default_rng(SEED)
    rows = []
    for beta, gamma, pi in itertools.product((0, 0.5, 1, 2), (0, 0.5, 1, 2), (0.05, 0.2)):
        reps = pd.DataFrame([simulate_once(beta, gamma, pi, rng) for _ in range(a.reps)])
        m = reps.mean()
        rows.append({"beta": beta, "gamma": gamma, "pi": pi, **{k: m[k] for k in reps},
                     "raw_err": m.raw - m.T, "matched_err": m.matched - m["T"], "deatt_err": m.matched_deattenuated - m["T"],
                     "corr_err": m.corrected_plus_half - m["T"]})
    t = pd.DataFrame(rows)
    g = t[t.gamma > 0]
    e1 = bool(((g.raw_err > 0) & (g.matched_err.abs() < g.raw_err.abs())).all())
    b0 = t[t.beta == 0]
    e2 = bool((b0.matched.sub(0.5).abs() < 0.03).all() and (b0.corrected_plus_half.sub(0.5).abs() < 0.03).all()
              and (b0[b0.gamma > 0].raw > 0.6).all())
    e3 = float(t.deatt_err.abs().mean())
    rep = ["A15 simulation: rarity-only null with known truth (T = AUROC pathogenic vs rare benign)", "",
           t.round(3).to_string(index=False), "",
           f"E1 (gamma > 0: raw overstates T and matched is closer to T than raw): {'HOLDS' if e1 else 'FAILS'}",
           f"E2 (beta = 0: matched and corrected + 0.5 near 0.5, raw not): {'HOLDS' if e2 else 'FAILS'}",
           f"E3 (deattenuated matched within 0.02 of T on average): {'HOLDS' if e3 <= 0.02 else 'FAILS'} (mean |error| {e3:.4f})",
           f"Corrected + 0.5 vs T: mean error {t.corr_err.mean():+.3f} (range {t.corr_err.min():+.3f} to {t.corr_err.max():+.3f})"]
    os.makedirs(RES, exist_ok=True)
    t.to_csv(os.path.join(RES, "a15_simulation.tsv"), sep="\t", index=False, float_format="%.4f", lineterminator="\n")
    open(os.path.join(RES, "a15_simulation.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


# --------------------------------------------------------------------------------------------- intervals

def window_index(d, pos_ids, neg_ids=None, flank=100):
    """Per positive: its row, its window's pool rows (rare/absent, not positive/negative) and, if neg_ids is None, its
    window's gnomAD-common rows as negatives (A9/A10 design); else the given published negatives (A11/A13 design)."""
    d = d.drop_duplicates("ID").sort_values(["chrom", "pos"]).reset_index(drop=True)
    common = ((d["filter"] == "PASS") & (d.AC >= 9)).to_numpy()
    isP = d.ID.isin(pos_ids).to_numpy()
    isN = d.ID.isin(neg_ids).to_numpy() if neg_ids is not None else common & ~isP
    pubneg = d.ID.isin(neg_ids).to_numpy() if neg_ids is not None else np.zeros(len(d), bool)
    pool = ~isP & ~common & ~pubneg
    if "win" in d:
        pool &= (d.win != "negative").to_numpy()
    units = []
    by = {c: x for c, x in d.groupby("chrom")}
    for i in np.flatnonzero(isP):
        x = by[d.chrom[i]]
        w = x.index[(x.pos >= d.pos[i] - flank) & (x.pos <= d.pos[i] + flank)].to_numpy()
        pl = w[pool[w]]
        if len(pl):
            units.append((i, pl, w[isN[w]] if neg_ids is None else None))
    return d, units, (np.flatnonzero(isN) if neg_ids is not None else None)


def metrics(s, units, sel, global_neg, rng, draws):
    pos = np.array([units[k][0] for k in sel])
    pools = [units[k][1] for k in sel]
    neg = np.unique(np.concatenate([units[k][2] for k in sel])) if global_neg is None else global_neg
    y = np.r_[np.ones(len(pos), bool), np.zeros(len(neg), bool)]
    raw = auroc(y, s[np.r_[pos, neg]])
    null = np.mean([auroc(y, s[np.r_[np.array([rng.choice(p) for p in pools]), neg]]) for _ in range(draws)])
    allpool = np.unique(np.concatenate(pools))
    M = auroc(np.r_[np.ones(len(pos), bool), np.zeros(len(allpool), bool)], s[np.r_[pos, allpool]])
    return raw, null, raw - null, M


def cmd_ci(a):
    which, _, name = a.analysis.partition(":")
    if which == "a9":
        d = pd.read_csv(os.path.join(HERE, "a9", "window_snvs.tsv.gz"), sep="\t")
        P, N = set(pd.read_csv(os.path.join(HERE, "a9", "positives.tsv"), sep="\t").ID), None
        tools = ["cadd_phred", "gpn_msa", "remm", "phylop447"]
    elif which == "a10":
        d = pd.read_csv(os.path.join(HERE, "a10", "missense_window_snvs.tsv.gz"), sep="\t")
        P, N = set(pd.read_csv(os.path.join(HERE, "a10", "missense_positives.tsv"), sep="\t").ID), None
        d = d[d.alphamissense.notna()]
        tools = ["alphamissense", "cadd_phred", "gpn_msa", "phylop447"]
    elif which in ("a11", "a13"):
        base = os.path.join(HERE, "a11") if which == "a11" else os.path.join(HERE, "a13", name)
        d = pd.read_csv(os.path.join(base, "window_snvs.tsv.gz"), sep="\t")
        if which == "a11":
            P = set(pd.read_csv(os.path.join(base, "positives_GRCh38.tsv"), sep="\t").ID)
            N = set(pd.read_csv(os.path.join(base, "negatives_GRCh38.tsv"), sep="\t").ID)
        else:
            b = pd.read_csv(os.path.join(base, "benchmark.tsv"), sep="\t")
            P, N = set(b.ID[b.label == 1]), set(b.ID[b.label == 0])
        tools = ["cadd_phred", "remm", "gpn_msa", "phylop447"]
    else:
        sys.exit("analysis must be a9, a10, a11 or a13:<name>")
    if which == "a9" and "phylop447" in d and d.phylop447.isna().all():
        pp = os.path.join(HERE, "a9_phylop", "window_snvs.tsv.gz")   # extras.py a9_phylop output, if present
        if os.path.exists(pp):
            d = d.drop(columns="phylop447").merge(pd.read_csv(pp, sep="\t", usecols=["ID", "phylop447"]).drop_duplicates("ID"), on="ID", how="left")
    d, units, gneg = window_index(d, P, N)
    tools = [t for t in tools if t in d and np.isfinite(d[t].to_numpy(float)).any()]
    S = {t: np.where(np.isnan(d[t].to_numpy(float)), np.nanmedian(d[t].to_numpy(float)), d[t].to_numpy(float)) for t in tools}
    rng = np.random.default_rng(SEED)
    full = np.arange(len(units))
    point = {t: metrics(S[t], units, full, gneg, rng, a.draws) for t in tools}
    boots = {t: [] for t in tools}
    for b in range(a.boot):
        sel = rng.integers(0, len(units), len(units))
        for t in tools:
            boots[t].append(metrics(S[t], units, sel, gneg, rng, a.draws))
        if (b + 1) % 20 == 0:
            print(f"  bootstrap {b + 1}/{a.boot}", flush=True)
    rows = []
    for t in tools:
        B = np.array(boots[t])
        lo, hi = np.percentile(B, [2.5, 97.5], axis=0)
        for k, nm in enumerate(["raw_auroc", "null_mean", "corrected", "rarity_matched_auroc"]):
            rows.append({"analysis": a.analysis, "tool": t, "metric": nm, "point": point[t][k], "lo": lo[k], "hi": hi[k]})
    for t1, t2 in itertools.combinations(tools, 2):     # paired differences in corrected and rarity-matched AUROC
        for k, nm in ((2, "corrected"), (3, "rarity_matched_auroc")):
            dd = np.array(boots[t1])[:, k] - np.array(boots[t2])[:, k]
            lo, hi = np.percentile(dd, [2.5, 97.5])
            rows.append({"analysis": a.analysis, "tool": f"{t1} - {t2}", "metric": f"diff_{nm}",
                         "point": point[t1][k] - point[t2][k], "lo": lo, "hi": hi})
    out = pd.DataFrame(rows)
    tag = a.analysis.replace(":", "_")
    out.to_csv(os.path.join(RES, f"a15_ci_{tag}.tsv"), sep="\t", index=False, float_format="%.4f", lineterminator="\n")
    print(f"A15 CIs for {a.analysis}: {len(units)} positives with windows; {a.boot} window-bootstrap replicates")
    print(out.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("simulate")
    s.add_argument("--reps", type=int, default=100)
    c = sub.add_parser("ci")
    c.add_argument("analysis", help="a9, a10, a11 or a13:<name>")
    c.add_argument("--boot", type=int, default=200)
    c.add_argument("--draws", type=int, default=50)
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    {"simulate": cmd_simulate, "ci": cmd_ci}[a.cmd](a)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
