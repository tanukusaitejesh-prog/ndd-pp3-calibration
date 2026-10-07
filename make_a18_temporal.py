"""A18: does stacked computational + location evidence predict later downgrade? (ClinVar 2019 -> 2026)

Specification: prereg_fullstudy/A18_TEMPORAL_RECLASSIFICATION_PREREG.md (pushed in 8a57850 before this ran).
Inputs: work/clinvar_20190923.vcf.gz and work/clinvar_20250504.vcf.gz (gnomAD public bucket,
gs://gcp-public-data--gnomad/resources/grch38/clinvar/), results_fullstudy/supp_data1_clinvar_pp3_pm1_points.tsv.gz.
Output: results_fullstudy/a18_temporal_summary.txt, results_fullstudy/a18_temporal_cohort.tsv.gz
"""
import gzip
import hashlib
import os

import numpy as np
import pandas as pd

R = "results_fullstudy"
W = "work"
AUTOSOMES = {str(i) for i in range(1, 23)}
REVSTAT = {"criteria_provided,_single_submitter", "criteria_provided,_multiple_submitters,_no_conflicts",
           "reviewed_by_expert_panel", "practice_guideline"}
STARS = {"criteria_provided,_single_submitter": 1, "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "reviewed_by_expert_panel": 3, "practice_guideline": 4}
PP3 = {"REVEL": ("REVEL_score", 0.644), "AlphaMissense": ("AlphaMissense_score", 0.792),
       "VEST4": ("VEST4_score", 0.764)}
SEED, NBOOT = 20261018, 2000
out_lines = []


def say(s=""):
    print(s)
    out_lines.append(s)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_vcf(path, missense_only=True):
    rows = []
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln[0] == "#":
                continue
            f = ln.split("\t", 8)
            if f[0] not in AUTOSOMES or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            if missense_only and "missense_variant" not in info.get("MC", ""):
                continue
            gene = info.get("GENEINFO", "").split("|")[0].split(":")[0]
            rows.append((f[0], int(f[1]), f[4], info.get("CLNSIG", ""), info.get("CLNREVSTAT", ""), gene))
    return pd.DataFrame(rows, columns=["chrom", "pos", "alt", "clnsig", "revstat", "gene"])


def mh_or(d, expo="grp", a_lab="Both", b_lab="One", y="down", strata="gene"):
    """Mantel-Haenszel OR (exposed a_lab vs b_lab), stratified by gene."""
    x = d[d[expo].isin([a_lab, b_lab])]
    num = den = 0.0
    for _, s in x.groupby(strata):
        e = s[expo] == a_lab
        if e.all() or (~e).all():
            continue
        n = len(s)
        a = (e & s[y]).sum(); b = (e & ~s[y]).sum()
        c = (~e & s[y]).sum(); dd = (~e & ~s[y]).sum()
        num += a * dd / n
        den += b * c / n
    return num / den if den > 0 else np.nan


def _mh_terms(d, expo="grp", a_lab="Both", b_lab="One", y="down", strata="gene"):
    """Per-gene Mantel-Haenszel numerator and denominator terms (0 for uninformative genes)."""
    x = d[d[expo].isin([a_lab, b_lab])]
    e = (x[expo] == a_lab).to_numpy()
    yy = x[y].to_numpy(bool)
    t = pd.DataFrame({"g": x[strata].to_numpy(), "a": e & yy, "b": e & ~yy, "c": ~e & yy, "d": ~e & ~yy})
    s = t.groupby("g")[["a", "b", "c", "d"]].sum()
    n = s.sum(axis=1)
    num = (s.a * s.d / n).where(n > 0, 0.0)
    den = (s.b * s.c / n).where(n > 0, 0.0)
    return num, den


def boot_mh(d, **kw):
    """Gene-clustered bootstrap: resample all cohort genes; uninformative genes contribute zero terms."""
    rng = np.random.default_rng(SEED)
    num, den = _mh_terms(d, **kw)
    genes = d.gene.unique()
    num = num.reindex(genes, fill_value=0.0).to_numpy()
    den = den.reindex(genes, fill_value=0.0).to_numpy()
    idx = rng.integers(0, len(genes), size=(NBOOT, len(genes)))
    sn, sd = num[idx].sum(1), den[idx].sum(1)
    v = np.where(sd > 0, sn / np.where(sd > 0, sd, 1), np.nan)
    v = v[np.isfinite(v) & (v > 0)]
    lo, hi = np.percentile(v, [2.5, 97.5])
    return lo, hi, float(np.std(np.log(v))), len(v)


def rate_ci(d, mask):
    rng = np.random.default_rng(SEED)
    x = d[mask]
    genes = x.gene.unique()
    by = {g: s.down.to_numpy() for g, s in x.groupby("gene")}
    vals = [np.concatenate([by[g] for g in rng.choice(genes, len(genes))]).mean() for _ in range(1000)]
    return x.down.mean(), *np.percentile(vals, [2.5, 97.5]), len(x)


def groups(d, tool):
    col, thr = PP3[tool]
    d = d[d[col].notna()].copy()
    pp3 = d[col] >= thr
    pm1 = d.pm1_like.astype(bool)
    d["grp"] = np.select([pp3 & pm1, pp3 ^ pm1], ["Both", "One"], "Neither")
    d["grp_detail"] = np.select([pp3 & pm1, pp3 & ~pm1, ~pp3 & pm1], ["Both", "PP3 only", "PM1-like only"], "Neither")
    return d


def main():
    p19, p25 = f"{W}/clinvar_20190923.vcf.gz", f"{W}/clinvar_20250504.vcf.gz"
    say("A18 temporal reclassification (pre-registered: prereg_fullstudy/A18_TEMPORAL_RECLASSIFICATION_PREREG.md, 8a57850)")
    say(f"md5 clinvar_20190923.vcf.gz {md5(p19)}; clinvar_20250504.vcf.gz {md5(p25)}")

    v19 = read_vcf(p19)
    v19 = v19[v19.revstat.isin(REVSTAT)]
    v19["stars"] = v19.revstat.map(STARS)
    s26 = pd.read_csv(f"{R}/supp_data1_clinvar_pp3_pm1_points.tsv.gz", sep="\t", dtype={"chrom": str})
    s26 = s26.rename(columns={"pos_grch38": "pos", "clinvar_class": "cls26"}).drop_duplicates(["chrom", "pos", "alt"])

    cohorts = {}
    for base, sig in [("LP", "Likely_pathogenic"), ("P", "Pathogenic")]:
        b = v19[v19.clnsig == sig].drop_duplicates(["chrom", "pos", "alt"])
        m = b.merge(s26, on=["chrom", "pos", "alt"], how="left")
        present = m.cls26.notna()
        say(f"\n2019 {sig} (>=1 star, autosomal missense SNV): {len(b)}; in 2026 table: {present.sum()}; "
            f"absent: {(~present).sum()}")
        m = m[present].copy()
        m["down"] = m.cls26.isin(["VUS", "B/LB"])
        cohorts[base] = (m, b)

    # ---- primary: 2019 LP, REVEL, MH OR Both vs One
    lp, lp_all = cohorts["LP"]
    d = groups(lp, "REVEL")
    say("\n== PRIMARY: 2019 LP, PP3 = REVEL >= 0.644, outcome 2026 VUS or B/LB ==")
    for g in ["Both", "One", "Neither"]:
        r, lo, hi, n = rate_ci(d, d.grp == g)
        say(f"  {g:8s} n {n:5d}  downgraded {r:.3f} [{lo:.3f}, {hi:.3f}]")
    for g in ["PP3 only", "PM1-like only"]:
        r, lo, hi, n = rate_ci(d, d.grp_detail == g)
        say(f"    {g:14s} n {n:5d}  downgraded {r:.3f} [{lo:.3f}, {hi:.3f}]")
    orr = mh_or(d)
    lo, hi, se, nb = boot_mh(d)
    inf = d[d.grp.isin(["Both", "One"])].groupby("gene").grp.nunique()
    n_inf = int((inf == 2).sum())
    mde = float(np.exp(2.8 * se))
    say(f"  MH OR Both vs One (gene-stratified): {orr:.3f} [{lo:.3f}, {hi:.3f}] "
        f"({n_inf} informative genes; {nb} valid bootstrap replicates; min detectable OR ~{mde:.2f})")
    verdict = "MAIN TEXT (interval excludes 1)" if (lo > 1 or hi < 1) else "SUPPLEMENT (null)"
    say(f"  Pre-specified rule -> {verdict}")

    # ---- secondary 2: other predictors
    say("\n== Secondary 2: other PP3 proxies (2019 LP) ==")
    for tool in ["AlphaMissense", "VEST4"]:
        dt = groups(lp, tool)
        o = mh_or(dt)
        l2, h2, _, _ = boot_mh(dt)
        rates = " ".join(f"{g} {dt[dt.grp == g].down.mean():.3f} (n {int((dt.grp == g).sum())})"
                         for g in ["Both", "One", "Neither"])
        say(f"  {tool:13s} MH OR {o:.3f} [{l2:.3f}, {h2:.3f}]; {rates}")

    # ---- secondary 3: 2019 Pathogenic cohort
    p, _ = cohorts["P"]
    dp = groups(p, "REVEL")
    o = mh_or(dp)
    l3, h3, _, _ = boot_mh(dp)
    rates = " ".join(f"{g} {dp[dp.grp == g].down.mean():.3f} (n {int((dp.grp == g).sum())})"
                     for g in ["Both", "One", "Neither"])
    say(f"\n== Secondary 3: 2019 Pathogenic cohort (REVEL) ==\n  MH OR {o:.3f} [{l3:.3f}, {h3:.3f}]; {rates}")

    # ---- secondary 4: logistic regression with stars, gene-clustered SEs
    import statsmodels.formula.api as smf
    dl = d.copy()
    dl["grp"] = pd.Categorical(dl.grp, ["One", "Both", "Neither"])
    dl["y"] = dl.down.astype(int)
    fit = smf.logit("y ~ C(grp) + C(stars)", data=dl).fit(
        disp=0, cov_type="cluster", cov_kwds={"groups": pd.factorize(dl.gene)[0]})
    say("\n== Secondary 4: logistic, downgrade ~ group + 2019 stars (gene-clustered SE; ref = One) ==")
    for k in fit.params.index:
        if k.startswith("C(grp)"):
            ci = fit.conf_int().loc[k]
            say(f"  {k:28s} OR {np.exp(fit.params[k]):.3f} [{np.exp(ci[0]):.3f}, {np.exp(ci[1]):.3f}] p {fit.pvalues[k]:.3g}")

    # ---- secondary 5: drop-out bound using ClinVar 2025-05
    v25 = read_vcf(p25)
    v25 = v25.drop_duplicates(["chrom", "pos", "alt"])
    absent = lp_all.merge(s26[["chrom", "pos", "alt"]], on=["chrom", "pos", "alt"], how="left", indicator=True)
    absent = absent[absent._merge == "left_only"].drop(columns="_merge")
    a = absent.merge(v25[["chrom", "pos", "alt", "clnsig", "revstat"]], on=["chrom", "pos", "alt"],
                     how="left", suffixes=("", "_25"))

    def bucket(s, r):
        if pd.isna(s):
            return "not in 2025 release"
        if "Conflicting" in s:
            return "conflicting"
        if s in ("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"):
            return "P/LP (0 star or other)" if r not in REVSTAT else "P/LP"
        if s == "Uncertain_significance":
            return "VUS"
        if s in ("Benign", "Likely_benign", "Benign/Likely_benign"):
            return "B/LB"
        return "other"
    a["status25"] = [bucket(s, r) for s, r in zip(a.clnsig_25, a.revstat_25)]
    say(f"\n== Secondary 5: 2019 LP variants absent from the 2026 table (n {len(a)}), status in ClinVar 2025-05 ==")
    say("  " + "; ".join(f"{k} {v}" for k, v in a.status25.value_counts().items()))
    k_down = int(a.status25.isin(["conflicting", "VUS", "B/LB"]).sum())
    # bound: add all downgraded drop-outs to Both (upper) or to One (lower), in a single pseudo-stratum per gene
    dd = d[["gene", "grp", "down"]].copy()
    extra = a[a.status25.isin(["conflicting", "VUS", "B/LB"])][["gene"]].assign(down=True)
    up = mh_or(pd.concat([dd, extra.assign(grp="Both")], ignore_index=True))
    low = mh_or(pd.concat([dd, extra.assign(grp="One")], ignore_index=True))
    say(f"  {k_down} drop-outs were downgraded or conflicting in 2025. Primary MH OR if all were Both: {up:.3f}; "
        f"if all were One: {low:.3f}")

    d.assign(cohort="2019_LP")[["chrom", "pos", "alt", "gene", "stars", "cls26", "down", "grp", "grp_detail",
                                "REVEL_score", "AlphaMissense_score", "pm1_like"]] \
        .to_csv(f"{R}/a18_temporal_cohort.tsv.gz", sep="\t", index=False)
    with open(f"{R}/a18_temporal_summary.txt", "w") as fh:
        fh.write("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    main()
