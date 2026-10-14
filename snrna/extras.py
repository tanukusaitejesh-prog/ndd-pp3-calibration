#!/usr/bin/env python3
"""Post hoc robustness and explanation analyses after A8 and A9 (exploratory; logged in PREREG s14 before running).

    python extras.py a8_sensitivity --blind sealed/blind_labels.tsv   # A8 under alternative literature label rules
    python extras.py a9_phylop                                        # add phyloP447 to A9 via the UCSC API, re-evaluate
    python extras.py mechanism --blind sealed/blind_labels.tsv        # per-gene hypermutability vs score performance

a8_sensitivity: the frozen wording rule rejects a variant if any source calls it VUS. Two alternatives: the curator's
  [STRICT] rows (explicit P/LP or disease-causing in at least one source), and the union of both rules. For each, new
  RNU2-2 positives (not sealed, not training) vs L3: AUROC, rarity-only null, excess and v4 - CADD, as a8_literature.py.
a9_phylop: phyloP447 per window from the UCSC REST API (as comparators.py), merged into the A9 window table; the
  frozen A9 evaluation is re-run on the augmented table into results/a9_phylop/.
mechanism: for every gene with >= 5 positives and >= 5 L3 (named snRNA genes incl. blind + A8 labels; A7 genes),
  the fraction of possible SNVs observed in gnomAD v4.1 genomes (PASS), mean Roulette rate, and AUROC of CADD,
  phyloP447, GPN-MSA, Rfam alt rarity and ref column freq. Descriptive (few genes): does genome-wide score performance
  fall where loci are hypermutable?
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snrna as S  # noqa: E402

SEED = 20261018
RES = os.path.join(HERE, "results")
HGVS = re.compile(r"^n\.(\d+)([ACGTU])>([ACGTU])$")
COMP = str.maketrans("ACGT", "TGCA")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def frozen_v4():
    sp = os.path.join(HERE, "results_v4", "v4_scores.tsv.gz")
    man = {l.split()[1]: l.split()[0] for l in open(os.path.join(HERE, "results_v4", "v4_manifest.txt")) if l.strip()}
    if man.get(os.path.relpath(sp, HERE)) != hashlib.sha256(open(sp, "rb").read()).hexdigest():
        sys.exit("v4_scores.tsv.gz does not match the frozen manifest")
    return pd.read_csv(sp, sep="\t")[["ID", "v4"]]


def null_excess(y_pos, neg_scores, pool_scores, x_pos, rng, draws=1000):
    """AUROC of positives vs negatives, its rarity-only null (random pool SNVs as positives) and one-sided p."""
    au = S.auroc(np.r_[np.ones(len(x_pos), bool), np.zeros(len(neg_scores), bool)], np.r_[x_pos, neg_scores])
    yy = np.r_[np.ones(len(x_pos), bool), np.zeros(len(neg_scores), bool)]
    nulls = np.array([S.auroc(yy, np.r_[rng.choice(pool_scores, len(x_pos), replace=False), neg_scores]) for _ in range(draws)])
    return au, float(nulls.mean()), float((np.sum(nulls >= au) + 1) / (draws + 1))


# ------------------------------------------------------------------------------------------- A8 sensitivity

def map_rows(t, genes, tx):
    out = []
    for r in t.itertuples():
        m = HGVS.match(str(r.hgvs_n).strip())
        if not m or r.gene not in genes.index:
            continue
        g = genes.loc[r.gene]
        n, ref, alt = int(m.group(1)), m.group(2).replace("T", "U"), m.group(3).replace("T", "U")
        if tx[g.label].get(n) != ref:
            continue
        pos = g.start + n - 1 if g.strand == "+" else g.end - n + 1
        rd, ad = ref.replace("U", "T"), alt.replace("U", "T")
        if g.strand == "-":
            rd, ad = rd.translate(COMP), ad.translate(COMP)
        out.append(f"{g.chrom}-{pos}-{rd}-{ad}")
    return set(out)


def cmd_a8_sensitivity(a):
    import evaluate as E
    rng = np.random.default_rng(SEED)
    genes = pd.read_csv(os.path.join(HERE, "model_inputs", "genes.tsv.gz"), sep="\t").set_index("gene")
    nodes = pd.read_csv(os.path.join(HERE, "model_inputs", "nodes.tsv.gz"), sep="\t")
    tx = {lab: dict(zip(x.n, x.base)) for lab, x in nodes[nodes.molecule == "target"].groupby("label")}
    lit = pd.read_csv(os.path.join(HERE, "lit", "a8_RNU2-2_2026-10-14.tsv"), sep="\t", dtype=str).fillna("")
    strict = map_rows(lit[lit.notes.str.startswith("[STRICT]")], genes, tx)
    frozen = set(pd.read_csv(os.path.join(HERE, "labels", "A8_literature.tsv"), sep="\t").query("gene == 'RNU2-2'").ID)
    blind = pd.read_csv(a.blind, sep="\t")
    sealed = set(blind.ID[(blind.label == "P/LP") & blind.ID.str.match(E.SNV_ID)])
    trained = set(E.load_l2(os.path.join(HERE, "labels")).ID)
    var = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    v = var[var.gene == "RNU2-2"].merge(frozen_v4(), on="ID", how="left")
    v = v.merge(pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")[["ID", "phylop447", "gpn_msa"]], on="ID", how="left")
    v["rfam_alt_rarity"], v["rfam_ref_col_freq"] = 1 - v.alt_col_freq, v.ref_col_freq
    allpos = sealed | frozen | strict
    neg = v.L3_benign_like.astype(bool) & ~v.ID.isin(allpos)
    pool = ~v.L3_benign_like.astype(bool) & ~v.ID.isin(allpos)
    rep = ["A8 label-rule sensitivity (post hoc): new RNU2-2 literature positives under three rules", ""]
    rows = []
    for name, P in (("frozen wording rule (as A8)", frozen), ("curator STRICT", strict), ("union", frozen | strict)):
        new = (P - sealed - trained) & set(v.ID)
        pos = v.ID.isin(new)
        rep.append(f"{name}: {int(pos.sum())} new positives vs {int(neg.sum())} L3")
        for c in ("v4", "rfam_alt_rarity", "rfam_ref_col_freq", "cadd_phred", "phylop447", "gpn_msa"):
            x = v[c].to_numpy(float)
            x = np.where(np.isnan(x), np.nanmedian(x), x)
            au, nm, p = null_excess(None, x[neg.to_numpy()], x[pool.to_numpy()], x[pos.to_numpy()], rng)
            rows.append({"rule": name, "score": c, "n_new": int(pos.sum()), "auroc": au, "null": nm, "excess": au - nm, "p": p})
            rep.append(f"  {c:18s} AUROC {au:.3f}  null {nm:.3f}  excess {au - nm:+.3f}  p {p:.4f}")
        y = np.r_[np.ones(int(pos.sum()), bool), np.zeros(int(neg.sum()), bool)]
        s4 = np.r_[v.v4[pos].to_numpy(float), v.v4[neg].to_numpy(float)]
        sc = np.r_[v.cadd_phred[pos].to_numpy(float), v.cadd_phred[neg].to_numpy(float)]
        i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
        d = [S.auroc(y[ix], s4[ix]) - S.auroc(y[ix], sc[ix]) for ix in
             (np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))] for _ in range(2000))]
        lo, hi = np.nanpercentile(d, [2.5, 97.5])
        rep += [f"  v4 - CADD {S.auroc(y, s4) - S.auroc(y, sc):+.3f} [{lo:+.3f}, {hi:+.3f}]", ""]
    pd.DataFrame(rows).to_csv(os.path.join(RES, "a8_sensitivity.tsv"), sep="\t", index=False, float_format="%.4f")
    open(os.path.join(RES, "a8_sensitivity.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


# ------------------------------------------------------------------------------------------------ A9 phyloP

def cmd_a9_phylop(a):
    import a9_inflation as A
    import comparators as CP
    d = pd.read_csv(os.path.join(HERE, "a9", "window_snvs.tsv.gz"), sep="\t")
    wins = sorted(set(d.win))
    log(f"phyloP447 for {len(wins)} windows via the UCSC API (cached; ~1 request/s)")
    vals = {}
    for k, w in enumerate(wins, 1):
        chrom, rng_ = w.split(":")
        s, e = map(int, rng_.split("-"))
        try:
            vals[w] = CP.ucsc_track("phyloP447way", chrom, s, e)
        except Exception as ex:
            log(f"  {w}: failed ({ex}); left missing")
            vals[w] = {}
        time.sleep(0.25)
        if k % 50 == 0 or k == len(wins):
            log(f"  {k}/{len(wins)} windows")
    d["phylop447"] = [vals[w].get(p, np.nan) for w, p in zip(d.win, d.pos)]
    out = os.path.join(HERE, "a9_phylop")
    os.makedirs(out, exist_ok=True)
    shutil.copy(os.path.join(HERE, "a9", "positives.tsv"), os.path.join(out, "positives.tsv"))
    d.to_csv(os.path.join(out, "window_snvs.tsv.gz"), sep="\t", index=False, float_format="%.5g", compression="gzip")
    log(f"phyloP447 coverage {d.phylop447.notna().mean():.2f}; re-running the frozen A9 evaluation")
    A.OUT, A.RES = out, os.path.join(RES, "a9_phylop")
    os.makedirs(A.RES, exist_ok=True)
    A.cmd_eval(argparse.Namespace(boot=1000))


# ------------------------------------------------------------------------------------------------- mechanism

def cmd_mechanism(a):
    import evaluate as E
    var = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    var = var.merge(pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")[["ID", "phylop447", "gpn_msa"]], on="ID", how="left")
    L = E.load_l2(os.path.join(HERE, "labels"))
    b = pd.read_csv(a.blind, sep="\t")
    lit = pd.read_csv(os.path.join(HERE, "labels", "A8_literature.tsv"), sep="\t")
    plp = set(L.ID[L.label == "P/LP"]) | set(b.ID[b.label == "P/LP"]) | set(lit.ID)
    var["pos"] = var.ID.isin(plp)
    var["neg"] = var.L3_benign_like.astype(bool) & ~var.pos
    var["observed"] = (var.AC_pass > 0) & (var.include == 1)
    var["MR"] = np.exp(var.log_mr)
    a7 = pd.read_csv(os.path.join(HERE, "a7", "a7_scores.tsv.gz"), sep="\t")
    a7lab = pd.read_csv(os.path.join(RES, "a7_labels.tsv"), sep="\t")
    a7["pos"] = a7.ID.isin(set(a7lab.ID))
    a7["neg"] = (a7["filter"] == "PASS") & (a7.AC_pass >= 9) & ~a7.pos
    a7["observed"] = (a7.AC_pass > 0) & (a7.include == 1)
    a7["MR"] = np.exp(a7.log_mr)
    a7["family"] = "non-snRNA (A7)"
    var["family"] = "spliceosomal snRNA"
    cols = ["ID", "gene", "family", "pos", "neg", "observed", "include", "MR", "cadd_phred", "phylop447", "gpn_msa", "alt_col_freq", "ref_col_freq"]
    t = pd.concat([var[cols], a7[cols]], ignore_index=True)
    t["rfam_alt_rarity"] = 1 - t.alt_col_freq
    t["rfam_ref_col_freq"] = t.ref_col_freq
    rows = []
    for g, x in t.groupby("gene"):
        if x.pos.sum() < 5 or x.neg.sum() < 5:
            continue
        inc = x[x.include == 1]
        r = {"gene": g, "family": x.family.iloc[0], "n_P": int(x.pos.sum()), "n_L3": int(x.neg.sum()),
             "frac_observed": float(inc.observed.mean()), "mean_roulette_rate": float(inc.MR.mean()),
             "obs_per_rate": float(inc.observed.mean() / inc.MR.mean()) if inc.MR.mean() > 0 else np.nan}
        xx = x[x.pos | x.neg]
        y = xx.pos.to_numpy()
        for c in ("cadd_phred", "phylop447", "gpn_msa", "rfam_alt_rarity", "rfam_ref_col_freq"):
            s = xx[c].to_numpy(float)
            r[f"auroc_{c}"] = S.auroc(y, np.where(np.isnan(s), np.nanmedian(s), s))
        rows.append(r)
    m = pd.DataFrame(rows).sort_values("frac_observed")
    m.to_csv(os.path.join(RES, "mechanism_by_gene.tsv"), sep="\t", index=False, float_format="%.4f")
    rep = ["Mechanism (exploratory, descriptive): per-gene hypermutability vs score performance (all positives incl. blind + A8)", "",
           m.to_string(index=False, float_format=lambda z: f"{z:.3f}"), ""]
    if len(m) >= 4:
        for c in ("cadd_phred", "phylop447", "gpn_msa", "rfam_alt_rarity", "rfam_ref_col_freq"):
            rho = S.spearman(m.frac_observed.to_numpy(float), m[f"auroc_{c}"].to_numpy(float))
            rep.append(f"  Spearman(fraction observed, AUROC {c}) = {rho:+.2f}  (n = {len(m)} genes)")
    open(os.path.join(RES, "mechanism_by_gene.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["a8_sensitivity", "a9_phylop", "mechanism"])
    ap.add_argument("--blind", default=os.path.join(HERE, "sealed", "blind_labels.tsv"))
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    {"a8_sensitivity": cmd_a8_sensitivity, "a9_phylop": cmd_a9_phylop, "mechanism": cmd_mechanism}[a.cmd](a)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
