#!/usr/bin/env python3
"""Amendment A7: zero-shot test of RNA-family covariance profiles on non-snRNA disease RNA genes (PREREG s15, A7).

Hypothesis (from the post hoc v4 ablation on blind RNU2-2): at structured non-coding RNA genes, the Rfam family
profile ranks patient variants above population variants better than genome-wide scores. Tested on genes that no
model, threshold or analysis in this project has touched. Two stages, with a commit between them:

    python a7_ncrna.py score     # label-free: features and scores for every SNV -> a7/a7_scores.tsv.gz + a7_manifest.txt
    (commit + push a7/a7_manifest.txt and a7/a7_scores.tsv.gz; ots stamp a7/a7_manifest.txt)
    python a7_ncrna.py eval      # refuses unless the manifest is committed and the scores match it -> results/a7_report.txt

Scores (higher = more damaging):
  rfam_alt_rarity   PRIMARY. 1 - frequency of the alternative base in the variant's Rfam consensus column (cmemit
                    profile, 20,000 samples, as model_features.py); insert-state positions get 1 (as frozen in v4).
                    Parameter-free; chosen as the best single feature of the post hoc v4 ablation on blind RNU2-2.
  v4_frozen         SECONDARY. The A5-frozen v4 coefficients (results_v4/v4_model.json) applied to the same feature
                    definitions; Evo 2 and lab function are missing here and take the frozen missing value (rank 0.5);
                    ddG_bind and alt_inter_broken are 0 (no partner RNAs defined for these families).
  rfam_ref_col_freq, cadd_phred, phylop447, gpn_msa, remm   comparators.
Labels (eval only): ClinVar P/LP SNVs with >= 1 star (positives) vs gnomAD v4.1 genome PASS AC >= 9 not P/LP.
Units: RMRP; TERC; SNORD118 + RNU7-1 pooled (each has < 10 P/LP SNVs). Feasibility counts seen before freezing
(14 Oct 2026, counts only): RMRP 28, TERC 15, SNORD118 6, RNU7-1 5 P/LP SNVs.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snrna as S  # noqa: E402

OUT = os.path.join(HERE, "a7")          # frozen, label-free scores (committed before eval)
RES = os.path.join(HERE, "results")     # eval outputs
SEED = 20261015
# GRCh38 (Ensembl lookup, 14 Oct 2026); Rfam 15.1 family; established inheritance
GENES = pd.DataFrame([
    ("RMRP", "ENSG00000277027", "chr9", 35657750, 35658019, "-", "RF00030", "recessive", "RMRP"),
    ("TERC", "ENSG00000270141", "chr3", 169764610, 169765060, "-", "RF00024", "dominant", "TERC"),
    ("SNORD118", "ENSG00000200463", "chr17", 8173452, 8173588, "-", "RF00096", "recessive", "SNORD118 + RNU7-1"),
    ("RNU7-1", "ENSG00000238923", "chr12", 6943816, 6943878, "+", "RF00066", "recessive", "SNORD118 + RNU7-1"),
], columns=["gene", "ensembl", "chrom", "start", "end", "strand", "rfam", "inheritance", "unit"])
SCORES = [("rfam_alt_rarity", "Rfam alt rarity [PRIMARY]"), ("v4_frozen", "v4 frozen coefficients [secondary]"),
          ("rfam_ref_col_freq", "Rfam ref column freq"), ("cadd_phred", "CADD v1.7"), ("phylop447", "phyloP447"),
          ("gpn_msa", "GPN-MSA"), ("remm", "ReMM v0.4")]
PRIMARY, COMPARATORS = "rfam_alt_rarity", ("cadd_phred", "phylop447")  # co-primary comparisons, Holm-adjusted
BASES = "ACGU"


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


# ------------------------------------------------------------------------------------------------------ score

def gene_snvs(g):
    seq = S.ucsc_seq(g.chrom, g.start - 1, g.end)
    rows = []
    L = len(seq)
    for i, ref in enumerate(seq):
        pos = g.start + i
        n = i + 1 if g.strand == "+" else L - i
        for alt in "ACGT":
            if alt != ref:
                rt, at = (ref, alt) if g.strand == "+" else (ref.translate(S.COMP), alt.translate(S.COMP))
                rows.append({"ID": f"{g.chrom}-{pos}-{ref}-{alt}", "gene": g.gene, "label": g.gene, "chrom": g.chrom,
                             "pos": pos, "ref": ref, "alt": alt, "n": n, "ref_t": rt.replace("T", "U"),
                             "alt_t": at.replace("T", "U")})
    tseq = seq if g.strand == "+" else S.revcomp(seq)
    return pd.DataFrame(rows), tseq.replace("T", "U")


def family_features(g, v, tseq):
    import model_features as MF
    cm = os.path.join(S.DATA, "rfam", f"{g.rfam}.cm")
    if not os.path.exists(cm):
        os.makedirs(os.path.dirname(cm), exist_ok=True)
        S.http_get(f"https://rfam.org/family/{g.rfam}/cm", cm)
    freq, _, _ = MF.consensus_profile(cm)
    cmap, _, bits = MF.cmalign_columns(cm, {g.gene: tseq})
    cmap = cmap[g.gene]
    with tempfile.TemporaryDirectory() as td:
        pairs, _ = S.cmalign_pairs(cm, {g.gene: tseq}, td)
    paired = {n: p for n, p, _, _ in pairs[g.gene]}
    k = v.n.map(lambda n: cmap[n - 1])
    v["cons_col"] = k
    v["alt_col_freq"] = [freq[kk][BASES.index(a)] if kk >= 0 else 0.0 for kk, a in zip(k, v.alt_t)]
    v["ref_col_freq"] = [freq[kk][BASES.index(r)] if kk >= 0 else 0.0 for kk, r in zip(k, v.ref_t)]
    v["ref_paired"] = v.n.map(lambda n: float(paired.get(n, False)))
    g0 = MF.mfe(tseq)
    dd = {}
    for n, a in set(zip(v.n, v.alt_t)):
        dd[(n, a)] = MF.mfe(tseq[:n - 1] + a + tseq[n:]) - g0
    v["abs_ddG_fold"] = [abs(dd[(n, a)]) for n, a in zip(v.n, v.alt_t)]
    log(f"  {g.gene}: Rfam {g.rfam} bits {bits[g.gene]:.1f}; {np.mean([c >= 0 for c in cmap]):.0%} of positions at consensus columns")
    return v, float(bits[g.gene]), float(np.mean([c >= 0 for c in cmap]))


def population(g, v):
    gn = S.gnomad_records(g.chrom, g.start, g.end)
    gn = gn.assign(ID=gn.chrom + "-" + gn.pos.astype(str) + "-" + gn.ref + "-" + gn.alt).drop_duplicates("ID").set_index("ID")
    ro = S.roulette_records(g.chrom, g.start, g.end)
    ro = ro.assign(ID=ro.chrom + "-" + ro.pos.astype(str) + "-" + ro.ref + "-" + ro.alt).drop_duplicates("ID").set_index("ID")
    flt = v.ID.map(gn["filter"])
    ac = v.ID.map(gn.AC).fillna(0)
    af = v.ID.map(gn.AF).fillna(0.0)
    passrec = flt == "PASS"
    v["filter"] = flt
    v["AC_pass"] = np.where(passrec, ac, 0).astype(int)
    v["MR"] = v.ID.map(ro.MR)
    v["log_mr"] = np.log(v.MR)
    # inclusion as in model_features.py (rare-variant sites with a usable mutation rate)
    v["include"] = (~(passrec & (af >= 0.001)) & ~(flt.notna() & ~passrec) & v.MR.notna()
                    & (v.ID.map(ro.roulette_filter).fillna("") != "low")).astype(int)
    return v


def comparators(g, v):
    import comparators as CP
    cadd = S.cadd_records(g.chrom, g.start, g.end)
    cadd = dict(zip(cadd.chrom + "-" + cadd.pos.astype(str) + "-" + cadd.ref + "-" + cadd.alt, cadd.cadd_phred))
    v["cadd_phred"] = v.ID.map(cadd)
    v["phylop447"] = v.pos.map(CP.ucsc_track("phyloP447way", g.chrom, g.start, g.end))
    remm = {int(r[1]): float(r[2]) for r in CP.tabix_rows(CP.REMM, g.chrom, g.start, g.end)}
    gpn = {(int(r[1]), r[2], r[3]): float(r[4]) for r in CP.tabix_rows(CP.GPN, g.chrom.replace("chr", ""), g.start, g.end)}
    v["remm"] = v.pos.map(remm)
    v["gpn_msa"] = [-gpn[k] if k in gpn else np.nan for k in zip(v.pos, v.ref, v.alt)]
    return v


def cmd_score(a):
    import v4 as V
    os.makedirs(OUT, exist_ok=True)
    parts, info = [], {}
    for g in GENES.itertuples():
        log(f"{g.gene} ({g.chrom}:{g.start}-{g.end} {g.strand})")
        v, tseq = gene_snvs(g)
        v, bits, cov = family_features(g, v, tseq)
        v = population(g, v)
        v = comparators(g, v)
        info[g.gene] = {"length": len(tseq), "rfam": g.rfam, "cm_bits": bits, "consensus_fraction": cov}
        parts.append(v)
    v = pd.concat(parts, ignore_index=True)
    v["rfam_alt_rarity"] = 1 - v.alt_col_freq
    v["rfam_ref_col_freq"] = v.ref_col_freq
    v["alt_rarity"] = v.rfam_alt_rarity
    v["local_constraint"] = V.local_constraint(v)
    v["evo2"] = np.nan
    v["lab_function"] = np.nan
    v["ddG_bind_max"] = 0.0
    v["alt_inter_broken"] = 0.0
    mdl = json.load(open(os.path.join(HERE, "results_v4", "v4_model.json")))
    x = mdl["intercept"]
    for f in mdl["features"]:
        x = x + mdl["coef"][f] * v.groupby("label")[f].rank(pct=True).fillna(0.5)
    v["v4_frozen"] = x
    keep = ["ID", "gene", "chrom", "pos", "ref", "alt", "n", "ref_t", "alt_t", "cons_col", "filter", "AC_pass", "include",
            "log_mr", "alt_col_freq", "ref_col_freq", "ref_paired", "abs_ddG_fold", "local_constraint"] + [c for c, _ in SCORES]
    sp = os.path.join(OUT, "a7_scores.tsv.gz")
    v[list(dict.fromkeys(keep))].to_csv(sp, sep="\t", index=False, float_format="%.6g", compression="gzip")
    json.dump({"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "genes": info, "n_snvs": len(v),
               "coverage": {c: float(v[c].notna().mean()) for c, _ in SCORES}}, open(os.path.join(OUT, "a7_score_info.json"), "w"), indent=1)
    files = [sp, os.path.join(OUT, "a7_score_info.json"), os.path.join(HERE, "a7_ncrna.py"),
             os.path.join(HERE, "results_v4", "v4_model.json")]
    with open(os.path.join(OUT, "a7_manifest.txt"), "w") as f:
        for p in files:
            f.write(f"{sha(p)}  {os.path.relpath(p, HERE)}\n")
    log(f"{len(v)} SNVs scored; no labels read. NEXT: commit and push a7/ (manifest + scores + info), "
        f"`ots stamp a7/a7_manifest.txt`, then run `python a7_ncrna.py eval`.")


# ------------------------------------------------------------------------------------------------------- eval

def check_frozen():
    man = os.path.join(OUT, "a7_manifest.txt")
    tracked = subprocess.run(["git", "-C", HERE, "ls-files", "--error-unmatch", man], capture_output=True).returncode == 0
    dirty = subprocess.run(["git", "-C", HERE, "status", "--porcelain", "--", OUT], capture_output=True, text=True).stdout.strip()
    if not tracked or dirty:
        sys.exit("a7/ must be committed (and unchanged) before eval: commit + push a7/, stamp the manifest, then rerun")
    for line in open(man):
        h, p = line.split()
        if p.endswith("v4_model.json") and not os.path.exists(os.path.join(HERE, p)):
            continue
        if sha(os.path.join(HERE, p)) != h:
            sys.exit(f"{p} does not match the committed A7 manifest")


def fetch_labels(acmg):
    import labels as LB
    rows = []
    for g in GENES.itertuples():
        for r in LB.clinvar(g.gene, g.start, g.end, g.chrom):
            if r["classification"] in LB.PLP and r["review_status"] not in LB.NO_STAR:
                rows.append({"ID": r["ID"], "gene": g.gene, "source": r["source"], "review_status": r["review_status"]})
    lab = pd.DataFrame(rows).drop_duplicates("ID")
    lab["pp3_bp4_applied"] = np.nan
    if acmg:
        import circularity_check as CC
        flags = {}
        for r in lab.itertuples():
            subs = []
            for acc in {s.split(":", 1)[1] for s in lab.source[lab.ID == r.ID]}:
                time.sleep(0.4)
                subs += CC.vcv_texts(acc)
            codes = {m.group(1) for s in subs for m in CC.ACMG.finditer(s["text"])}
            flags[r.ID] = float(bool({"PP3", "BP4"} & codes)) if codes else np.nan
        lab["pp3_bp4_applied"] = lab.ID.map(flags)
    return lab


def boot_delta(units, a_col, b_col, B, rng):
    """Weighted (by positives) mean of per-unit AUROC(a) - AUROC(b), class-stratified bootstrap within unit."""
    w = np.array([y.sum() for y, _ in units], float)
    pt = float((w * np.array([S.auroc(y, d[a_col]) - S.auroc(y, d[b_col]) for y, d in units])).sum() / w.sum())
    bs = []
    for _ in range(B):
        vals = []
        for y, d in units:
            i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
            idx = np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]
            vals.append(S.auroc(y[idx], d[a_col][idx]) - S.auroc(y[idx], d[b_col][idx]))
        bs.append((w * np.array(vals)).sum() / w.sum())
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    p = float(min(1.0, 2 * min(np.mean(np.array(bs) <= 0), np.mean(np.array(bs) >= 0))))
    return pt, lo, hi, p


def evaluate(v, lab, name, B, rng, keep=None):
    v = v if keep is None else v[keep]
    pos = v.ID.isin(set(lab.ID))
    neg = (v["filter"] == "PASS") & (v.AC_pass >= 9) & ~pos
    rows, units = [], []
    for unit, x in v[pos | neg].groupby(v.gene.map(dict(zip(GENES.gene, GENES.unit)))):
        x = x.dropna(subset=[c for c, _ in SCORES if c not in ("gpn_msa", "remm", "phylop447")])
        y = x.ID.isin(set(lab.ID)).to_numpy()
        row = {"set": name, "unit": unit, "n_P": int(y.sum()), "n_N": int((~y).sum())}
        if y.sum() >= 5 and (~y).sum() >= 5:
            d = {}
            for c, _ in SCORES:
                s = x[c].to_numpy(float)
                d[c] = np.where(np.isnan(s), np.nanmedian(s), s)  # comparator gaps -> median (coverage reported)
                row[f"auroc_{c}"] = S.auroc(y, d[c])
                row[f"{c}_lo"], row[f"{c}_hi"] = S.boot_auroc(y, d[c], B)
            units.append((y, d))
        rows.append(row)
    res = {}
    if units:
        for comp in COMPARATORS + ("gpn_msa",):
            res[comp] = boot_delta(units, PRIMARY, comp, B, rng)
        res["v4_vs_cadd"] = boot_delta(units, "v4_frozen", "cadd_phred", B, rng)
    return pd.DataFrame(rows), res


def cmd_eval(a):
    check_frozen()
    v = pd.read_csv(os.path.join(OUT, "a7_scores.tsv.gz"), sep="\t")
    lab = fetch_labels(a.acmg)
    lab.to_csv(os.path.join(RES, "a7_labels.tsv"), sep="\t", index=False)
    rng = np.random.default_rng(SEED)
    sets = {"PRIMARY: ClinVar P/LP >= 1 star vs PASS AC >= 9": (lab, None),
            "consensus-column positions only": (lab, v.cons_col >= 0),
            "negatives PASS AC >= 20": (lab, ~((v["filter"] == "PASS") & v.AC_pass.between(9, 19)))}
    if a.acmg:
        sets["positives without PP3/BP4"] = (lab[lab.pp3_bp4_applied != 1], None)
    tabs, rep = [], ["Amendment A7: zero-shot RNA-family profile on non-snRNA disease RNA genes (scores frozen before labels)", ""]
    decision = None
    for name, (L, keep) in sets.items():
        t, res = evaluate(v, L, name, a.boot, rng, keep)
        tabs.append(t)
        rep.append(f"== {name}")
        for r in t.itertuples():
            rep.append(f"  {r.unit}: {r.n_P} P vs {r.n_N} N")
            for c, lbl in SCORES:
                if hasattr(r, f"auroc_{c}") and pd.notna(getattr(r, f"auroc_{c}")):
                    rep.append(f"    {lbl:36s} AUROC {getattr(r, f'auroc_{c}'):.3f} [{getattr(r, f'{c}_lo'):.3f}, {getattr(r, f'{c}_hi'):.3f}]")
        if res:
            ps = sorted((res[c][3], c) for c in COMPARATORS)
            holm = {c: min(1.0, p * (len(ps) - i)) for i, (p, c) in enumerate(ps)}
            for c in COMPARATORS + ("gpn_msa",):
                pt, lo, hi, p = res[c]
                rep.append(f"  weighted delta AUROC primary - {c}: {pt:+.3f} [{lo:+.3f}, {hi:+.3f}]  p {p:.4f}"
                           + (f"  Holm p {holm[c]:.4f}" if c in holm else "  (reported)"))
            pt, lo, hi, p = res["v4_vs_cadd"]
            rep.append(f"  weighted delta AUROC v4_frozen - CADD: {pt:+.3f} [{lo:+.3f}, {hi:+.3f}] (secondary)")
            if name.startswith("PRIMARY"):
                t_ok = t.dropna(subset=[f"auroc_{PRIMARY}"])
                n70 = int((t_ok[f"auroc_{PRIMARY}"] >= 0.70).sum())
                beats = all(res[c][1] > 0 and holm[c] < 0.05 for c in COMPARATORS)
                decision = (f"A7 decision: {'MET' if beats and n70 >= 2 else 'NOT MET'} (pre-registered: weighted delta vs CADD and "
                            f"vs phyloP447 both > 0 with interval excluding 0 and Holm p < 0.05, and primary AUROC >= 0.70 in "
                            f">= 2 of {len(t_ok)} units; here {n70} units >= 0.70, comparisons {'pass' if beats else 'fail'})")
        rep.append("")
    if decision:
        rep.append(decision)
    pd.concat(tabs).to_csv(os.path.join(RES, "a7_results.tsv"), sep="\t", index=False, float_format="%.4f")
    open(os.path.join(RES, "a7_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["score", "eval"])
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--acmg", action="store_true", help="eval: also read ACMG codes from ClinVar submissions (slow)")
    a = ap.parse_args()
    {"score": cmd_score, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
