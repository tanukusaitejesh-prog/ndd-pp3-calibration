#!/usr/bin/env python3
"""Amendment A12: prospective blind tests of the frozen v4 on patient variants that appear after the freeze, including
genes that never had a single training label (PREREG s15, A12). Written before any A12 label was built or scored.

v4 and its scores for all 17,385 SNVs were hash-frozen on 14 Oct 2026 (A5). Any P/LP variant that is not in a label file
committed by then (L2_nonblind, L2_curated_additions, A8_literature, the sealed blind set) is unseen by v4. Several
snRNA genes never had any label at all (RNU5A-1, RNU6-2, RNU6-7, RNU11, RNU4-1, RNU5D/E/F-1, RNVU1-*), so they are
blind genes by construction. A12 collects such variants from ClinVar (re-scanned monthly) and from the literature
(curated with lit/A12_TALLY_PROMPT.md), and evaluates each pre-registered unit once, when it reaches its minimum.

    python a12_new_variants.py scan                       # ClinVar VCF -> a12/clinvar_scan_<fileDate>.tsv + a12/scan_log.tsv
    python a12_new_variants.py map lit/a12_*.tsv          # literature tallies -> labels/A12_literature{,_rejected,_indels}.tsv
    python a12_new_variants.py build [--blind sealed/blind_labels.tsv]   # -> labels/A12_labels.tsv (commit it)
    python a12_new_variants.py status                     # positives per unit vs the pre-registered minimum (no scores)
    python a12_new_variants.py eval --unit U1_RNU6_adRP [--blind ...]    # one look per unit; refuses if labels uncommitted

Units (pre-registered; a gene belongs to one unit):
  U1_RNU6_adRP   RNU6-1, RNU6-2, RNU6-7, RNU6-8 (adRP; dominant)         PRIMARY  >= 8 new positives
  U2_U5          RNU5A-1, RNU5B-1 (U5 NDD; dominant)                      PRIMARY  >= 5   (RNU5B-1 needs --blind)
  U3_RNU2-2      RNU2-2 beyond the sealed and A8 sets                     PRIMARY  >= 8   (needs --blind)
  U4_other       every other locus without training labels               secondary >= 5
  U5_trained     RNU4-2, RNU4ATAC, RNU12, RNU6ATAC (new variants only)    secondary >= 5   (genes seen in training)
Positives: ClinVar P/LP with >= 1 review star, or literature P/LP / disease-causing under the frozen wording rule of
curate_literature.py, SNVs, not in any earlier label file; a literature positive that ClinVar classifies B/LB (>= 1 star)
is dropped and listed. Negatives: L3 (gnomAD v4.1 PASS AC >= 9) in the unit's loci, minus every positive ever used.
Rarity-only null: 1,000 draws of as many SNVs from the unit's loci that are neither L3 nor any positive.
Scores: frozen v4 (PRIMARY; hash-checked against results_v4/v4_manifest.txt), Rfam ref column freq (key secondary,
the Rfam feature that held up in A7), Rfam alt rarity, CADD, phyloP447, GPN-MSA.
Decision per primary unit: MET if (i) v4 excess over its rarity null has one-sided p < 0.05/3 (Bonferroni over the three
primary units, fixed in advance so that units evaluated at different times share one alpha) and (ii) v4 - CADD AUROC
has a 95% class-stratified bootstrap interval above 0. Below the minimum: descriptive only. Each unit is evaluated once
(a12/eval_log.tsv); a repeat needs --again and is labelled exploratory.
Reclassification arm (descriptive): committed ClinVar VUS (labels/clinvar_vus_snrna.tsv) that later become P/LP or B/LB;
AUROC of frozen v4 for P/LP vs B/LB once each class has >= 5.
"""
import argparse
import bisect
import gzip
import hashlib
import os
import re
import subprocess
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import curate_literature as CL  # noqa: E402

OUT = os.path.join(HERE, "a12")
LABELS = os.path.join(HERE, "labels")
LAB = os.path.join(LABELS, "A12_labels.tsv")
LIT = os.path.join(LABELS, "A12_literature.tsv")
CLINVAR = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
SEED = 20261021
NULL_DRAWS = 1000
ALPHA = 0.05 / 3
SEALED = ("RNU2-2", "RNU5B-1")
UNITS = {
    "U1_RNU6_adRP": {"genes": ("RNU6-1", "RNU6-2", "RNU6-7", "RNU6-8"), "primary": True, "min": 8},
    "U2_U5": {"genes": ("RNU5A-1", "RNU5B-1"), "primary": True, "min": 5},
    "U3_RNU2-2": {"genes": ("RNU2-2",), "primary": True, "min": 8},
    "U5_trained": {"genes": ("RNU4-2", "RNU4ATAC", "RNU12", "RNU6ATAC"), "primary": False, "min": 5},
    "U4_other": {"genes": None, "primary": False, "min": 5},   # everything else
}
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BLB = {"Benign", "Likely_benign", "Benign/Likely_benign"}
VUS = {"Uncertain_significance", "Conflicting_classifications_of_pathogenicity", "Conflicting_interpretations_of_pathogenicity"}
NO_STAR = {"no_assertion_criteria_provided", "no_classification_provided", "no_interpretation_for_the_single_variant",
           "no_classifications_from_unflagged_records", "no_assertion_provided"}
PRIOR_FILES = ("L2_nonblind.tsv", "L2_curated_additions.tsv", "A8_literature.tsv")
SCORES = [("v4", "v4 (frozen) [PRIMARY]"), ("rfam_ref_col_freq", "Rfam ref column freq [key secondary]"),
          ("rfam_alt_rarity", "Rfam alt rarity"), ("cadd_phred", "CADD"), ("phylop447", "phyloP447"), ("gpn_msa", "GPN-MSA")]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def unit_of(gene):
    for u, d in UNITS.items():
        if d["genes"] and gene in d["genes"]:
            return u
    return "U4_other"


def genes_table():
    return pd.read_csv(os.path.join(HERE, "model_inputs", "genes.tsv.gz"), sep="\t")


def prior_ids(blind=None):
    """Every SNV already used as a label anywhere in the study (any class)."""
    ids = set()
    for f in PRIOR_FILES:
        p = os.path.join(LABELS, f)
        if os.path.exists(p):
            ids |= set(pd.read_csv(p, sep="\t").ID)
    if blind:
        ids |= set(pd.read_csv(blind, sep="\t").ID)
    return ids


def committed(path):
    tracked = subprocess.run(["git", "-C", HERE, "ls-files", "--error-unmatch", path], capture_output=True).returncode == 0
    dirty = subprocess.run(["git", "-C", HERE, "status", "--porcelain", "--", path], capture_output=True, text=True).stdout.strip()
    return tracked and not dirty


# ------------------------------------------------------------------------------------------------------------ scan

def parse_clinvar(vcf, genes):
    """Every ClinVar record inside an analysis-set locus: SNVs with classification and star flag; indels counted."""
    iv = {}
    for r in genes.itertuples():
        iv.setdefault(r.chrom.replace("chr", ""), []).append((r.start, r.end, r.label, r.gene))
    for c in iv:
        iv[c].sort()
    starts = {c: [x[0] for x in v] for c, v in iv.items()}
    rows, fdate, n_indel = [], "", 0
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##fileDate"):
                fdate = line.strip().split("=")[1]
            if line.startswith("#"):
                continue
            c = line.split("\t", 8)
            if c[0] not in iv:
                continue
            pos = int(c[1])
            k = bisect.bisect_right(starts[c[0]], pos) - 1
            if k < 0 or pos > iv[c[0]][k][1]:
                continue
            s, e, label, gene = iv[c[0]][k]
            if len(c[3]) != 1 or len(c[4]) != 1 or c[4] not in "ACGT":
                n_indel += 1
                continue
            info = dict(kv.split("=", 1) for kv in c[7].strip().split(";") if "=" in kv)
            sig = info.get("CLNSIG", "")
            cls = "P/LP" if sig in PLP else "B/LB" if sig in BLB else "VUS" if sig in VUS else "other"
            rev = info.get("CLNREVSTAT", "")
            rows.append({"ID": f"chr{c[0]}-{pos}-{c[3]}-{c[4]}", "gene": gene, "label": label, "clinvar_class": cls,
                         "clnsig": sig, "review_status": rev, "star1": not (set(rev.split(",")) & NO_STAR),
                         "variation_id": c[2]})
    return pd.DataFrame(rows, columns=["ID", "gene", "label", "clinvar_class", "clnsig", "review_status", "star1",
                                       "variation_id"]), fdate, n_indel


def cmd_scan(a):
    os.makedirs(OUT, exist_ok=True)
    vcf = a.vcf or os.path.join(HERE, "data", f"a12_clinvar_{time.strftime('%Y%m%d')}.vcf.gz")  # never the frozen A9 file
    if not os.path.exists(vcf):
        os.makedirs(os.path.dirname(vcf), exist_ok=True)
        log(f"downloading {CLINVAR} -> {vcf}")
        urllib.request.urlretrieve(CLINVAR, vcf)
    d, fdate, n_indel = parse_clinvar(vcf, genes_table())
    prior = prior_ids()
    d["prior_label"] = d.ID.isin(prior)
    d["unit"] = d.gene.map(unit_of)
    d.to_csv(os.path.join(OUT, f"clinvar_scan_{fdate}.tsv"), sep="\t", index=False)
    vus = pd.read_csv(os.path.join(LABELS, "clinvar_vus_snrna.tsv"), sep="\t")
    now = d.set_index("ID").clinvar_class
    recl = vus.assign(now=vus.ID.map(now)).query("now in ['P/LP', 'B/LB']")
    recl.to_csv(os.path.join(OUT, f"vus_reclassified_{fdate}.tsv"), sep="\t", index=False)
    new = d[(d.clinvar_class == "P/LP") & d.star1 & ~d.prior_label]
    rows = []
    for g, x in d.groupby("gene"):
        n = new[new.gene == g]
        rows.append({"clinvar_fileDate": fdate, "gene": g, "unit": unit_of(g), "snvs_in_clinvar": len(x),
                     "plp_star1": int(((x.clinvar_class == "P/LP") & x.star1).sum()), "plp_any": int((x.clinvar_class == "P/LP").sum()),
                     "new_plp_star1": len(n), "new_needs_blind_check": len(n) if g in SEALED else 0,
                     "blb_star1": int(((x.clinvar_class == "B/LB") & x.star1).sum()), "vus": int((x.clinvar_class == "VUS").sum()),
                     "vus_to_plp": int(((recl.gene == g) & (recl.now == "P/LP")).sum()),
                     "vus_to_blb": int(((recl.gene == g) & (recl.now == "B/LB")).sum()),
                     "scan_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
    s = pd.DataFrame(rows)
    lp = os.path.join(OUT, "scan_log.tsv")
    old = pd.read_csv(lp, sep="\t") if os.path.exists(lp) else pd.DataFrame()
    old = old[old.get("clinvar_fileDate", pd.Series(dtype=str)).astype(str) != str(fdate)] if len(old) else old
    pd.concat([old, s], ignore_index=True).to_csv(lp, sep="\t", index=False)
    log(f"ClinVar fileDate {fdate}: {len(d)} SNV records in {d.gene.nunique()} loci ({n_indel} indels skipped)")
    print(s[["gene", "unit", "plp_star1", "new_plp_star1", "blb_star1", "vus", "vus_to_plp", "vus_to_blb"]].to_string(index=False))
    if (s.new_needs_blind_check > 0).any():
        print("NOTE: RNU2-2/RNU5B-1 'new' counts still include sealed labels; `build --blind` on the DGX removes them")
    print(f"NEXT: commit a12/clinvar_scan_{fdate}.tsv a12/scan_log.tsv a12/vus_reclassified_{fdate}.tsv")


# ------------------------------------------------------------------------------------------------------------- map

def cmd_map(a):
    genes = genes_table()
    named = genes[genes.named.astype(bool)].drop_duplicates("gene").set_index("gene")
    nodes = pd.read_csv(os.path.join(HERE, "model_inputs", "nodes.tsv.gz"), sep="\t")
    t = nodes[nodes.molecule == "target"]
    tx = {lab: dict(zip(x.n, x.base)) for lab, x in t.groupby("label")}
    rows, indels, rejected = [], [], []
    for p in a.files:
        for r in pd.read_csv(p, sep="\t", dtype=str).fillna("").to_dict("records"):
            gene = r["gene"].strip()
            cls = r.get("classification_as_stated", "").lower()
            if gene not in named.index:
                rejected.append({**r, "reason": "not a named analysis-set gene"}); continue
            if not any(w in cls for w in CL.POSITIVE_WORDS) or any(w in cls for w in CL.NOT_PER_VARIANT):
                rejected.append({**r, "reason": "not reported as P/LP or disease-causing"}); continue
            g = named.loc[gene]
            hg = r.get("hgvs_n", "").strip().replace(" ", "")
            m = CL.HGVS_SNV.match(hg)
            if not m:
                indels.append({**r, "file": os.path.basename(p)}); continue
            n, ref, alt = int(m.group(1)), m.group(2).replace("T", "U"), m.group(3).replace("T", "U")
            if tx[g.label].get(n) != ref:
                rejected.append({**r, "reason": f"reference {ref} at n.{n} but the transcript has {tx[g.label].get(n)}"}); continue
            pos = g.start + n - 1 if g.strand == "+" else g.end - n + 1
            gref, galt = (CL.to_dna(ref), CL.to_dna(alt)) if g.strand == "+" else \
                (CL.to_dna(ref).translate(CL.COMP), CL.to_dna(alt).translate(CL.COMP))
            vid = f"{g.chrom}-{pos}-{gref}-{galt}"
            gm = re.match(r"^(?:chr)?(\w+)[:\-](?:g\.)?(\d+)", r.get("genomic_grch38", "").strip())
            if gm and [gm.group(1), gm.group(2)] != [g.chrom.replace("chr", ""), str(pos)]:
                rejected.append({**r, "reason": f"source GRCh38 {r['genomic_grch38']} disagrees with the mapping {vid}"}); continue
            inh = r.get("inheritance", "").strip().lower()
            rows.append({"ID": vid, "gene": gene, "source": f"literature:{r.get('source', '')}", "hgvs_n": hg,
                         "classification": r.get("classification_as_stated", ""), "label": "P/LP",
                         "inheritance": inh if inh in ("dominant", "recessive") else "unknown",
                         "first_publication_date": r.get("first_publication_date", "") or r.get("first_publication_year", ""),
                         "insilico_or_conservation_used_for_classification": r.get("insilico_or_conservation_used_for_classification", "unknown"),
                         "tally_date": time.strftime("%Y-%m-%d", time.gmtime())})
    out = pd.DataFrame(rows, columns=["ID", "gene", "source", "hgvs_n", "classification", "label", "inheritance",
                                      "first_publication_date", "insilico_or_conservation_used_for_classification", "tally_date"])
    if len(out):
        src = out.groupby("ID").source.agg(lambda s: ";".join(sorted(set(s))))
        out = out.drop_duplicates("ID").assign(source=lambda d: d.ID.map(src))
    out.to_csv(LIT, sep="\t", index=False)
    pd.DataFrame(rejected).to_csv(LIT.replace(".tsv", "_rejected.tsv"), sep="\t", index=False)
    pd.DataFrame(indels).to_csv(LIT.replace(".tsv", "_indels.tsv"), sep="\t", index=False)
    print(f"A12 literature SNVs {len(out)}; indels {len(indels)}; rejected {len(rejected)}")
    if len(out):
        print(out.groupby("gene").ID.nunique().to_string())


# ----------------------------------------------------------------------------------------------------------- build

def build_labels(scans, lit, prior, have_blind):
    """Pure function (tested): new positives from ClinVar scans + literature, one row per SNV."""
    rows, dropped = [], []
    cv = pd.concat(scans, ignore_index=True) if scans else pd.DataFrame(columns=["ID", "gene", "clinvar_class", "star1"])
    if len(cv):
        cv = cv.sort_values("clinvar_fileDate").drop_duplicates("ID", keep="last")   # latest classification wins
    blb = set(cv.ID[(cv.clinvar_class == "B/LB") & cv.star1.astype(bool)]) if len(cv) else set()
    for r in cv[(cv.clinvar_class == "P/LP") & cv.star1.astype(bool)].itertuples() if len(cv) else []:
        rows.append({"ID": r.ID, "gene": r.gene, "source": f"clinvar:{r.variation_id}", "classification": r.clnsig,
                     "review_status": r.review_status, "first_seen": r.clinvar_fileDate})
    for r in lit.itertuples() if len(lit) else []:
        if r.ID in blb:
            dropped.append({"ID": r.ID, "gene": r.gene, "reason": "literature P/LP but ClinVar B/LB (>= 1 star)"}); continue
        rows.append({"ID": r.ID, "gene": r.gene, "source": r.source, "classification": r.classification,
                     "review_status": "literature", "first_seen": r.tally_date})
    out = pd.DataFrame(rows, columns=["ID", "gene", "source", "classification", "review_status", "first_seen"])
    if len(out):
        src = out.groupby("ID").source.agg(lambda s: ";".join(sorted(set(s))))
        out = out.drop_duplicates("ID").assign(source=lambda d: d.ID.map(src))
    keep = ~out.ID.isin(prior)
    if not have_blind:
        dropped += [{"ID": i, "gene": g, "reason": "sealed gene: rebuild with --blind on the DGX"}
                    for i, g in zip(out.ID[out.gene.isin(SEALED) & keep], out.gene[out.gene.isin(SEALED) & keep])]
        keep &= ~out.gene.isin(SEALED)
    out = out[keep].assign(unit=lambda d: d.gene.map(unit_of), label="P/LP")
    return out.sort_values(["unit", "gene", "ID"]).reset_index(drop=True), pd.DataFrame(dropped)


def cmd_build(a):
    import glob
    scans = [pd.read_csv(p, sep="\t").assign(clinvar_fileDate=os.path.basename(p)[13:-4])
             for p in sorted(glob.glob(os.path.join(OUT, "clinvar_scan_*.tsv")))]
    lit = pd.read_csv(LIT, sep="\t") if os.path.exists(LIT) else pd.DataFrame()
    out, dropped = build_labels(scans, lit, prior_ids(a.blind), bool(a.blind))
    out.assign(build_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())).to_csv(LAB, sep="\t", index=False)
    dropped.to_csv(LAB.replace(".tsv", "_dropped.tsv"), sep="\t", index=False)
    print(f"A12 labels: {len(out)} new positives; dropped {len(dropped)}")
    status(out)
    print("NEXT: commit labels/A12_labels.tsv (+ _dropped); `eval --unit <U>` only for a unit at or above its minimum")


def status(lab):
    for u, d in UNITS.items():
        n = int((lab.unit == u).sum()) if len(lab) else 0
        tag = "PRIMARY" if d["primary"] else "secondary"
        print(f"  {u:14s} {tag:9s} {n:3d} new positives (minimum {d['min']}){'  -> evaluable' if n >= d['min'] else ''}")


def cmd_status(a):
    status(pd.read_csv(LAB, sep="\t") if os.path.exists(LAB) else pd.DataFrame(columns=["unit"]))


# ------------------------------------------------------------------------------------------------------------ eval

def load_scores(v4_dir):
    var = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    var = var.merge(pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")[["ID", "phylop447", "gpn_msa"]],
                    on="ID", how="left")
    var["rfam_alt_rarity"], var["rfam_ref_col_freq"] = 1 - var.alt_col_freq, var.ref_col_freq
    sp, mp = os.path.join(v4_dir, "v4_scores.tsv.gz"), os.path.join(v4_dir, "v4_manifest.txt")
    note = "frozen v4 hash verified"
    if os.path.exists(sp) and os.path.exists(mp):
        man = {l.split()[1]: l.split()[0] for l in open(mp) if l.strip() and not l.startswith("#")}
        if man.get(os.path.relpath(sp, HERE).replace("\\", "/")) != hashlib.sha256(open(sp, "rb").read()).hexdigest():
            sys.exit("results_v4/v4_scores.tsv.gz does not match the frozen manifest")
        var = var.merge(pd.read_csv(sp, sep="\t")[["ID", "v4"]], on="ID", how="left")
    else:
        var["v4"], note = np.nan, "results_v4/ absent: v4 not scored, decision pending (run on the DGX)"
    return var, note


def unit_metrics(g, P, allpos, rng, boot, auroc, boot_auroc, scores=SCORES):
    """Metrics for one unit (tested): g = variant table of the unit's loci, P = positive IDs."""
    pos = g.ID.isin(P).to_numpy()
    l3 = g.L3_benign_like.astype(bool).to_numpy()
    neg = l3 & ~g.ID.isin(allpos).to_numpy()
    pool = np.flatnonzero(~l3 & ~g.ID.isin(allpos).to_numpy())
    idx = np.flatnonzero(pos | neg)
    y = pos[idx]
    rows, res = [], {}
    negi = np.flatnonzero(neg)
    for c, lab in scores:
        full = g[c].to_numpy(float)
        if not np.isfinite(full).any():
            continue
        full = np.where(np.isnan(full), np.nanmedian(full), full)
        s = full[idx]
        au = auroc(y, s)
        lo, hi = boot_auroc(y, s, boot)
        k = int(y.sum())
        nulls = np.array([auroc(np.r_[np.ones(k, bool), np.zeros(len(negi), bool)],
                                full[np.r_[rng.choice(pool, k, replace=len(pool) < k), negi]]) for _ in range(NULL_DRAWS)])
        p = float((np.sum(nulls >= au) + 1) / (len(nulls) + 1))
        res[c] = s
        rows.append({"score": lab, "n_P": k, "n_N": int((~y).sum()), "auroc": au, "lo": lo, "hi": hi,
                     "null_mean": float(nulls.mean()), "excess": au - float(nulls.mean()), "p_vs_null": p})
    delta = None
    if "v4" in res and "cadd_phred" in res:
        i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
        d = [auroc(y[ix], res["v4"][ix]) - auroc(y[ix], res["cadd_phred"][ix])
             for ix in (np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))] for _ in range(boot))]
        delta = (auroc(y, res["v4"]) - auroc(y, res["cadd_phred"]), *np.nanpercentile(d, [2.5, 97.5]))
    return rows, delta


def decide(unit, rows, delta, n):
    d = UNITS[unit]
    if not d["primary"]:
        return f"{unit}: secondary, descriptive"
    if n < d["min"]:
        return f"{unit}: DESCRIPTIVE ONLY ({n} new positives < {d['min']})"
    v4 = [r for r in rows if r["score"].startswith("v4")]
    if not v4 or delta is None:
        return f"{unit}: PENDING (v4 not scored here)"
    ok = v4[0]["p_vs_null"] < ALPHA and delta[1] > 0
    return (f"{unit}: {'MET' if ok else 'NOT MET'} (v4 excess p {v4[0]['p_vs_null']:.4f} [< {ALPHA:.4f} needed]; "
            f"v4 - CADD {delta[0]:+.3f} [{delta[1]:+.3f}, {delta[2]:+.3f}] [interval > 0 needed]; {n} new positives)")


def cmd_eval(a):
    import snrna as S
    if not committed(LAB):
        sys.exit("labels/A12_labels.tsv must be committed and unchanged before eval")
    lp = os.path.join(OUT, "eval_log.tsv")
    elog = pd.read_csv(lp, sep="\t") if os.path.exists(lp) else pd.DataFrame(columns=["unit"])
    if a.unit in set(elog.unit) and not a.again:
        sys.exit(f"{a.unit} was already evaluated (a12/eval_log.tsv); --again runs an exploratory repeat")
    lab = pd.read_csv(LAB, sep="\t")
    if lab.gene.isin(SEALED).any() and not a.blind:
        sys.exit("labels include sealed-gene rows: pass --blind sealed/blind_labels.tsv")
    var, note = load_scores(a.v4)
    allpos = prior_ids(a.blind) | set(lab.ID)
    genes = UNITS[a.unit]["genes"] or sorted(set(var.gene) - {g for d in UNITS.values() if d["genes"] for g in d["genes"]})
    g = var[var.gene.isin(genes)].reset_index(drop=True)
    P = set(lab.ID[lab.unit == a.unit]) & set(g.ID)
    rng = np.random.default_rng(SEED)
    rep = [f"A12 {a.unit}{' (EXPLORATORY repeat)' if a.again else ''}: {len(P)} new positives in {', '.join(genes)}; {note}", ""]
    rows, delta = unit_metrics(g, P, allpos, rng, a.boot, S.auroc, S.boot_auroc) if len(P) >= 3 else ([], None)
    for r in rows:
        rep.append(f"  {r['score']:38s} AUROC {r['auroc']:.3f} [{r['lo']:.3f}, {r['hi']:.3f}]  rarity null {r['null_mean']:.3f}"
                   f"  excess {r['excess']:+.3f}  p {r['p_vs_null']:.4f}")
    if delta:
        rep.append(f"  v4 - CADD {delta[0]:+.3f} [{delta[1]:+.3f}, {delta[2]:+.3f}]")
    dec = decide(a.unit, rows, delta, len(P)) if len(P) >= 3 else f"{a.unit}: too few positives ({len(P)})"
    rep += ["", dec]
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    tag = a.unit + ("_again" if a.again else "")
    pd.DataFrame(rows).assign(unit=a.unit).to_csv(os.path.join(HERE, "results", f"a12_{tag}.tsv"), sep="\t", index=False, float_format="%.4f")
    open(os.path.join(HERE, "results", f"a12_{tag}_report.txt"), "w").write("\n".join(rep) + "\n")
    commit = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    pd.concat([elog, pd.DataFrame([{"unit": a.unit, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "n_P": len(P),
                                    "decision": dec, "exploratory_repeat": bool(a.again), "labels_commit": commit}])],
              ignore_index=True).to_csv(lp, sep="\t", index=False)
    print("\n".join(rep))
    # reclassification arm (descriptive)
    import glob
    rec = [pd.read_csv(p, sep="\t") for p in sorted(glob.glob(os.path.join(OUT, "vus_reclassified_*.tsv")))]
    if rec and np.isfinite(var.v4).any():
        r = pd.concat(rec).drop_duplicates("ID", keep="last").merge(var[["ID", "v4", "cadd_phred"]], on="ID")
        if (r.now == "P/LP").sum() >= 5 and (r.now == "B/LB").sum() >= 5:
            y = (r.now == "P/LP").to_numpy()
            print(f"VUS reclassified: {int(y.sum())} -> P/LP, {int((~y).sum())} -> B/LB; v4 AUROC {S.auroc(y, r.v4.to_numpy(float)):.3f}, "
                  f"CADD {S.auroc(y, r.cadd_phred.to_numpy(float)):.3f} (descriptive)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("scan")
    s.add_argument("--vcf", help="local ClinVar GRCh38 VCF (default: download today's into data/)")
    m = sub.add_parser("map")
    m.add_argument("files", nargs="+")
    b = sub.add_parser("build")
    b.add_argument("--blind", help="decrypted sealed labels (DGX): needed to include RNU2-2 and RNU5B-1")
    sub.add_parser("status")
    e = sub.add_parser("eval")
    e.add_argument("--unit", required=True, choices=list(UNITS))
    e.add_argument("--blind")
    e.add_argument("--v4", default=os.path.join(HERE, "results_v4"))
    e.add_argument("--boot", type=int, default=2000)
    e.add_argument("--again", action="store_true", help="exploratory repeat of an already evaluated unit (logged)")
    a = ap.parse_args()
    {"scan": cmd_scan, "map": cmd_map, "build": cmd_build, "status": cmd_status, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
