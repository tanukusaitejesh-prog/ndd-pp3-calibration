#!/usr/bin/env python3
"""Post hoc, exploratory (after the single unsealing): is the blind RNU2-2 result circular?

Two routes by which the labels could already contain the predictors' information:

1. Positive labels. A ClinVar submitter may have used in-silico or conservation evidence (PP3, CADD, phyloP, ...)
   to call a variant P/LP; v4 uses conservation features, so such positives would be easy for it. For every blind
   P/LP SNV, all ClinVar records at that variant are fetched (VCV XML, every submission) and the ACMG criteria the
   submitters applied are read from the submission text. Classes:
     insilico   - PP3 or BP4 applied by at least one submitter       (excluded from both filtered sets)
     clean      - criteria listed, neither PP3 nor BP4
     unknown    - no criteria stated
   Filtered sets: "strict" = clean only; "lenient" = clean + unknown. Also: strict without PM2 (population rarity,
   which overlaps the L3 definition and the constraint feature), and PS2 (de novo) only. In-silico keywords are kept
   as an informational column; they also match mechanism prose ("predicted to disrupt"), so they do not classify.
2. Negative labels. L3 is gnomAD AC >= 9, and v4's local-constraint feature counts observed SNVs in a 21-nt window
   that includes the variant's own position, so each L3 variant raises its own window's count. Check: recompute the
   feature leaving the focal position out, re-rank within gene, and score with the FROZEN v4 coefficients (no refit).

Scores: frozen v4 (hash-checked, results_v4/ is read only), the v2 framework S with its components f and F (scored
exactly as final_eval.py does), CADD and the comparators.

    python circularity_check.py --blind none --fetch-only                  # any machine: sealed/blind_evidence.tsv for
                                                                           # every ClinVar/RNUdb P/LP SNV (not committed)
    python circularity_check.py --blind sealed/blind_labels.tsv --fetch   # or fetch for the blind P/LP set itself
    python circularity_check.py --blind sealed/blind_labels.tsv           # reuses it -> results/circularity_check.{txt,tsv}
"""
import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evaluate as E  # noqa: E402
import final_eval as FE  # noqa: E402
import labels as LB  # noqa: E402
import snrna as S  # noqa: E402
import v4 as V  # noqa: E402
from gene_discovery import fit_k  # noqa: E402

GENE = "RNU2-2"
INSILICO = re.compile(r"\bPP3\b|\bBP4\b|in[ -]?silico|computational|predict|conserv|\bCADD\b|phylop|\bGERP|"
                      r"spliceai|\bREVEL\b|\bremm\b", re.I)  # informational only: matches mechanism prose too
ACMG = re.compile(r"\b(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(?:_[A-Za-z-]+)?\b")


# --------------------------------------------------------------------------------------- 1. label evidence

def vcv_texts(acc):
    """Submission-level text (comments, method descriptions, evidence descriptions) of one ClinVar VCV record."""
    num = int(re.sub(r"\D", "", acc.split(".")[0]))
    url = f"{LB.EUTILS}/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id={num}&from_esearch=true"
    for i in range(6):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "snrna-vep/2 (research)"})
            with urllib.request.urlopen(req, timeout=120) as r:
                root = ET.fromstring(r.read())
            break
        except Exception as e:
            time.sleep(1.5 * (i + 1))
            last = e
    else:
        raise RuntimeError(f"{acc}: {last}")
    out = []
    for ca in root.iter("ClinicalAssertion"):
        sub = ca.find("ClinVarAccession")
        out.append({"who": sub.get("SubmitterName", "") if sub is not None else "",
                    "comment": " ".join((c.text or "").strip() for c in ca.iter("Comment")).strip(),
                    "text": " ".join((el.text or "").strip() for el in ca.iter()
                                     if (el.text or "").strip() and not (el.tag == "Attribute" and el.get("Type") == "HGVS")),
                    "origin": {(o.text or "").strip() for o in ca.iter("Origin")},
                    "pmids": {(i.text or "").strip() for i in ca.iter("ID") if i.get("Source") == "PubMed"}})
    return out


def fetch_evidence(plp_ids, inputs, path):
    """Evidence class per P/LP SNV. plp_ids None (no blind file at hand): every ClinVar P/LP (>= 1 star) and RNUdb P/LP
    SNV of the gene, a superset of the blind positives, which are a subset of exactly these sources (labels.py)."""
    genes = pd.read_csv(os.path.join(inputs, "genes.tsv.gz"), sep="\t")
    g = genes[genes.gene == GENE].iloc[0]
    accs, cv_plp = {}, set()
    for r in LB.clinvar(GENE, g.start, g.end, g.chrom):  # every ClinVar record at these positions, not only the kept one
        accs.setdefault(r["ID"], set()).add(r["source"].split(":", 1)[1])
        if r["classification"] in LB.PLP and r["review_status"] not in LB.NO_STAR:
            cv_plp.add(r["ID"])
    in_rnudb = {r["ID"] for r in LB.rnudb(GENE)}
    if plp_ids is None:
        plp_ids = {i for i in cv_plp | in_rnudb if E.SNV_ID.match(i)}
    rows = []
    for i, vid in enumerate(sorted(plp_ids), 1):
        print(f"  ClinVar evidence {i}/{len(plp_ids)}: {vid}", flush=True)
        subs = []
        for acc in sorted(accs.get(vid, ())):
            time.sleep(0.4)
            subs += vcv_texts(acc)
        hits = sorted({m.group(0).lower() for s in subs for m in INSILICO.finditer(s["text"])})
        codes = sorted({m.group(1) for s in subs for m in ACMG.finditer(s["text"])})
        rnudb = vid in in_rnudb
        cls = "insilico" if {"PP3", "BP4"} & set(codes) else "clean" if codes else "unknown"
        origin = set().union(*[s["origin"] for s in subs]) if subs else set()
        rows.append({"ID": vid, "n_submissions": len(subs), "submitters": "; ".join(sorted({s["who"] for s in subs if s["who"]})),
                     "rnudb": rnudb, "clinvar_plp": vid in cv_plp, "origin": ",".join(sorted(origin)),
                     "de_novo": any("de novo" in o.lower() for o in origin),
                     "pmids": ",".join(sorted(set().union(*[s["pmids"] for s in subs]) - {""})) if subs else "",
                     "acmg_codes": ",".join(codes), "insilico_terms": ",".join(hits), "evidence_class": cls,
                     "comments": " || ".join(s["comment"] for s in subs if s["comment"])[:2000]})
    ev = pd.DataFrame(rows)
    ev.to_csv(path, sep="\t", index=False)
    return ev


# ------------------------------------------------------------------------- 2. leave-position-out constraint

def local_constraint_lpo(v, w=21):
    """V.local_constraint with the focal position's own observed/expected/variance removed from its window."""
    out = pd.Series(np.nan, index=v.index)
    ok = (v.include == 1) & v.log_mr.notna()
    for lab, g in v[ok].groupby("label"):
        m = np.exp(g.log_mr.to_numpy()); obs = (g.AC_pass > 0).to_numpy(float)
        p = 1 - np.exp(-fit_k(obs, m) * m)
        n = g.n.to_numpy(); n0 = n.min(); L = n.max() - n0 + 1
        O = np.bincount(n - n0, obs, L); Ex = np.bincount(n - n0, p, L); Vr = np.bincount(n - n0, p * (1 - p), L)
        h = w // 2
        cs = lambda x: np.concatenate([[0], np.cumsum(x)])
        cO, cE, cV = cs(O), cs(Ex), cs(Vr)
        lo = np.clip(np.arange(L) - h, 0, L); hi = np.clip(np.arange(L) + h + 1, 0, L)
        o = cO[hi] - cO[lo] - O; e = cE[hi] - cE[lo] - Ex; var = cV[hi] - cV[lo] - Vr
        z = (o - e) / np.sqrt(np.maximum(var, 1e-9))
        allg = v.label == lab
        out[allg] = -z[(v.n[allg] - n0).clip(0, L - 1).to_numpy()]
    return out


def v4_rescored(inputs, evo2, v4dir):
    """Frozen v4 coefficients applied to (a) the original features, as a reproduction check, and (b) the features with
    leave-position-out local constraint."""
    ic = os.path.join(HERE, "model_inputs_v3") if os.path.exists(os.path.join(HERE, "model_inputs_v3", "variants.tsv.gz")) else inputs
    v = V.build_features(inputs, ic, evo2)
    vc = pd.read_csv(os.path.join(ic, "variants.tsv.gz"), sep="\t")
    v["local_constraint_lpo"] = local_constraint_lpo(v.assign(AC_pass=v.ID.map(dict(zip(vc.ID, vc.AC_pass)))))
    v["local_constraint_lpo_r"] = v.groupby("label").local_constraint_lpo.rank(pct=True).fillna(0.5)
    mdl = json.load(open(os.path.join(v4dir, "v4_model.json")))
    lin = lambda cols: mdl["intercept"] + sum(mdl["coef"][f] * v[c] for f, c in zip(mdl["features"], cols))
    X = [f + "_r" for f in mdl["features"]]
    v["v4_repro"] = lin(X)
    v["v4_lpo"] = lin([("local_constraint_lpo_r" if c == "local_constraint_r" else c) for c in X])
    return v[["ID", "v4_repro", "v4_lpo"]]


# ------------------------------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blind", required=True)
    ap.add_argument("--fetch", action="store_true", help="(re)fetch ClinVar submission evidence")
    ap.add_argument("--evidence", default=os.path.join(HERE, "sealed", "blind_evidence.tsv"))
    ap.add_argument("--v4", default=os.path.join(HERE, "results_v4"))
    ap.add_argument("--runs", default=os.path.join(HERE, "model_runs"))
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--fetch-only", action="store_true",
                    help="no blind file needed: evidence for every ClinVar/RNUdb P/LP SNV of the gene, then stop")
    a = ap.parse_args()
    if a.fetch_only:
        os.makedirs(os.path.dirname(a.evidence), exist_ok=True)
        ev = fetch_evidence(None, a.inputs, a.evidence)
        print(ev.drop(columns="comments").to_string(index=False))
        print(ev.evidence_class.value_counts().to_string())
        return

    sp = os.path.join(a.v4, "v4_scores.tsv.gz")
    man = {l.split()[1]: l.split()[0] for l in open(os.path.join(a.v4, "v4_manifest.txt")) if l.strip()}
    if man.get(os.path.relpath(sp, HERE)) != hashlib.sha256(open(sp, "rb").read()).hexdigest():
        sys.exit("v4_scores.tsv.gz does not match the frozen manifest")
    blind = pd.read_csv(a.blind, sep="\t")
    blind = blind[blind.ID.str.match(E.SNV_ID)]
    bplp = set(blind.ID[(blind.gene == GENE) & (blind.label == "P/LP")])
    ev = fetch_evidence(bplp, a.inputs, a.evidence) if a.fetch or not os.path.exists(a.evidence) else \
        pd.read_csv(a.evidence, sep="\t")
    if missing := bplp - set(ev.ID):
        sys.exit(f"{len(missing)} blind P/LP SNVs have no evidence row; rerun with --fetch")

    # v2 framework scores for the blind gene, exactly as final_eval.blind_gene builds them
    pred, _ = E.load_preds(a.runs)
    var = pd.read_csv(os.path.join(a.inputs, "variants.tsv.gz"), sep="\t")
    genes = pd.read_csv(os.path.join(a.inputs, "genes.tsv.gz"), sep="\t")
    cl = genes.set_index("label").cluster
    hom = E.homology_scores(a.inputs)
    sc = E.unit_scores(pred, var, genes, list(genes.label[genes.gene == GENE]), lambda l: f"lofo_{cl[l]}", hom)
    sc["F"] = sc.h.where(sc.h.notna(), sc.d)
    sc["S"] = E.framework_score(sc, FE.BLIND_MODE[GENE])
    comp = pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")
    sc = sc.merge(comp[["ID", "phylop447", "remm", "gpn_msa"]], on="ID", how="left")
    print("scoring: v2 framework, comparators, frozen v4, leave-position-out v4 ...", flush=True)
    v4 = pd.read_csv(sp, sep="\t")[["ID", "v4"]]
    sc = sc.merge(v4, on="ID", how="left").merge(v4_rescored(a.inputs, a.evo2, a.v4), on="ID", how="left")

    plp_all = bplp
    cls = dict(zip(ev.ID, ev.evidence_class))
    dn = set(ev.ID[ev.de_novo.astype(str).str.lower() == "true"])
    codes = {r.ID: set(str(r.acmg_codes).split(",")) if pd.notna(r.acmg_codes) else set() for r in ev.itertuples()}
    subsets = {"all (as frozen)": plp_all,
               "lenient: no PP3/BP4 applied": {i for i in plp_all if cls.get(i) in ("clean", "unknown")},
               "strict: criteria listed, no PP3/BP4": {i for i in plp_all if cls.get(i) == "clean"},
               "strict, and PM2 (population rarity) not applied": {i for i in plp_all if cls.get(i) == "clean" and "PM2" not in codes[i]},
               "PS2 applied (de novo)": {i for i in plp_all if "PS2" in codes[i]},
               "de novo origin in a submission": plp_all & dn,
               "PP3/BP4 applied": {i for i in plp_all if cls.get(i) == "insilico"}}
    cols = [("v4", "v4 (frozen)"), ("v4_lpo", "v4, leave-position-out constraint"), ("S", "v2 framework S"),
            ("f", "v2 selection f"), ("F", "v2 function F"), ("cadd_phred", "CADD"), ("phylop447", "phyloP447"),
            ("remm", "ReMM"), ("gpn_msa", "GPN-MSA")]
    rows = []
    for name, plp in subsets.items():
        pos = sc.ID.isin(plp)
        neg = sc.L3_benign_like.astype(bool) & ~sc.ID.isin(plp_all)  # L3 fixed; only positives are filtered
        for col, lab in cols:
            ok = (pos | neg) & sc[col].notna()
            y = pos[ok].to_numpy(); x = sc[col][ok].to_numpy(float)
            row = {"label_set": name, "score": lab, "n_P": int(y.sum()), "n_L3": int((~y).sum())}
            if y.sum() >= 3 and (~y).sum():
                row["auroc"] = S.auroc(y, x)
                row["lo"], row["hi"] = S.boot_auroc(y, x, a.boot)
            rows.append(row)
    out = pd.DataFrame(rows)
    rep = ["Circularity check, blind RNU2-2 (post hoc, exploratory; frozen v4 hash verified, nothing refitted)", "",
           "Positive-label evidence classes: " + ", ".join(f"{k} {v}" for k, v in
                                                           ev[ev.ID.isin(plp_all)].evidence_class.value_counts().items()),
           "Classes from the ACMG criteria applied in the ClinVar submissions (PP3/BP4 = in silico); L3 negatives unchanged.", ""]
    for name, x in out.groupby("label_set", sort=False):
        rep.append(f"{name}: {x.n_P.iloc[0]} P/LP vs {x.n_L3.iloc[0]} L3")
        for r in x.itertuples():
            rep.append(f"  {r.score:36s} " + (f"AUROC {r.auroc:.3f} [{r.lo:.3f}, {r.hi:.3f}]" if pd.notna(getattr(r, 'auroc', np.nan)) else "n/a (< 3 positives)"))
        rep.append("")
    d = sc[["v4", "v4_repro"]].dropna()
    rep.append(f"Reproduction of the frozen v4 from v4_model.json: max |diff| {np.abs(d.v4 - d.v4_repro).max():.2e} "
               f"over {len(d)} SNVs (should be ~0; if not, the v4_lpo row is not comparable)")
    out.to_csv(os.path.join(a.results, "circularity_check.tsv"), sep="\t", index=False, float_format="%.4f")
    open(os.path.join(a.results, "circularity_check.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
