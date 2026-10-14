#!/usr/bin/env python3
"""Post hoc, exploratory (after the single unsealing): is the blind RNU2-2 result circular?

Two routes by which the labels could already contain the predictors' information:

1. Positive labels. A ClinVar submitter may have used in-silico or conservation evidence (PP3, CADD, phyloP, ...)
   to call a variant P/LP; v4 uses conservation features, so such positives would be easy for it. For every blind
   P/LP SNV, all ClinVar records at that variant are fetched (VCV XML, every submission) and each submission's text is
   searched for in-silico/conservation wording. Classes:
     insilico   - at least one submission mentions it               (excluded from both filtered sets)
     clean      - submissions carry text, none mentions it, or RNUdb literature curation only
     unknown    - no submission text to judge
   Filtered sets: "strict" = clean only; "lenient" = clean + unknown. RNUdb entries are literature curations (de novo
   occurrence, recurrence, segregation); treating them as clean is an assumption stated in the report.
2. Negative labels. L3 is gnomAD AC >= 9, and v4's local-constraint feature counts observed SNVs in a 21-nt window
   that includes the variant's own position, so each L3 variant raises its own window's count. Check: recompute the
   feature leaving the focal position out, re-rank within gene, and score with the FROZEN v4 coefficients (no refit).

Scores: frozen v4 (hash-checked, results_v4/ is read only), the v2 framework S with its components f and F (scored
exactly as final_eval.py does), CADD and the comparators.

    python circularity_check.py --blind sealed/blind_labels.tsv --fetch   # writes sealed/blind_evidence.tsv (not committed)
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
                      r"spliceai|\bREVEL\b|\bremm\b", re.I)


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
        who = sub.get("SubmitterName", "") if sub is not None else ""
        txt = " ".join(t.strip() for el in ca.iter() if el.tag in ("Comment", "Description", "Attribute")
                       for t in [el.text or ""] if t.strip())
        out.append((who, txt))
    return out


def fetch_evidence(blind, inputs, path):
    genes = pd.read_csv(os.path.join(inputs, "genes.tsv.gz"), sep="\t")
    g = genes[genes.gene == GENE].iloc[0]
    accs = {}
    for r in LB.clinvar(GENE, g.start, g.end, g.chrom):  # every ClinVar record at these positions, not only the kept one
        accs.setdefault(r["ID"], set()).add(r["source"].split(":", 1)[1])
    plp = blind[(blind.gene == GENE) & (blind.label == "P/LP")]
    for r in plp.itertuples():
        if r.source.startswith("ClinVar:"):
            accs.setdefault(r.ID, set()).add(r.source.split(":", 1)[1])
    in_rnudb = {r["ID"] for r in LB.rnudb(GENE)} | set(plp.ID[plp.source == "RNUdb"])  # labels.collect keeps one source
    rows = []
    for vid in sorted(set(plp.ID)):
        subs = []
        for acc in sorted(accs.get(vid, ())):
            time.sleep(0.4)
            subs += vcv_texts(acc)
        hits = sorted({m.group(0).lower() for _, t in subs for m in INSILICO.finditer(t)})
        rnudb = vid in in_rnudb
        has_text = any(t for _, t in subs)
        cls = "insilico" if hits else "clean" if (has_text or (rnudb and not subs)) else "unknown"
        rows.append({"ID": vid, "n_submissions": len(subs), "submitters": "; ".join(sorted({w for w, _ in subs if w})),
                     "rnudb": bool(rnudb), "insilico_terms": ",".join(hits), "evidence_class": cls})
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
    a = ap.parse_args()

    sp = os.path.join(a.v4, "v4_scores.tsv.gz")
    man = {l.split()[1]: l.split()[0] for l in open(os.path.join(a.v4, "v4_manifest.txt")) if l.strip()}
    if man.get(os.path.relpath(sp, HERE)) != hashlib.sha256(open(sp, "rb").read()).hexdigest():
        sys.exit("v4_scores.tsv.gz does not match the frozen manifest")
    blind = pd.read_csv(a.blind, sep="\t")
    blind = blind[blind.ID.str.match(E.SNV_ID)]
    ev = fetch_evidence(blind, a.inputs, a.evidence) if a.fetch or not os.path.exists(a.evidence) else \
        pd.read_csv(a.evidence, sep="\t")

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
    v4 = pd.read_csv(sp, sep="\t")[["ID", "v4"]]
    sc = sc.merge(v4, on="ID", how="left").merge(v4_rescored(a.inputs, a.evo2, a.v4), on="ID", how="left")

    plp_all = set(blind.ID[(blind.gene == GENE) & (blind.label == "P/LP")])
    cls = dict(zip(ev.ID, ev.evidence_class))
    subsets = {"all (as frozen)": plp_all,
               "lenient: no in-silico mention": {i for i in plp_all if cls.get(i) in ("clean", "unknown")},
               "strict: documented non-in-silico": {i for i in plp_all if cls.get(i) == "clean"},
               "in-silico mentioned only": {i for i in plp_all if cls.get(i) == "insilico"}}
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
           "RNUdb-only positives are counted as clean (literature curation; assumption).", ""]
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
