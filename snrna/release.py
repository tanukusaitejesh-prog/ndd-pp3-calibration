#!/usr/bin/env python3
"""Package the scores for a dated public release (Zenodo), PREREG s11 pseudo-prospective test (amendment A3).

    python release.py --stage scores    # right after the run: scores only (no thresholds), starts the prospective clock
    python release.py --stage final     # after final_eval.py: scores + calibrated PP3/BP4 thresholds and evidence labels

Writes release/snrna-vep-v2-<stage>-<date>/ and a .zip next to it: scores.tsv.gz (every possible SNV in the 55
analysis-set snRNA loci), thresholds.tsv (final only), README.md, zenodo.json, SHA256SUMS. Upload the zip to Zenodo
(new upload -> files -> metadata from zenodo.json) and record the DOI and date in PREREG s14.
"""
import argparse
import hashlib
import json
import os
import shutil
import time

import numpy as np
import pandas as pd

import calibrate as C

HERE = os.path.dirname(os.path.abspath(__file__))
STRENGTH_ORDER = ("Strong", "Moderate", "Supporting")


def evidence(row, th):
    """Highest PP3 or BP4 strength whose threshold the framework score clears for this gene (prior 0.10)."""
    if th is None or pd.isna(row.S):
        return ""
    for k in STRENGTH_ORDER:
        t = th.get(f"PP3_{k}")
        if t is not None and not pd.isna(t) and row.S >= t:
            return f"PP3_{k}"
    for k in STRENGTH_ORDER:
        t = th.get(f"BP4_{k}")
        if t is not None and not pd.isna(t) and row.S <= t:
            return f"BP4_{k}"
    return ""


README = """# snRNA-VEP v2 scores ({stage}, {date})

Variant-effect scores for every possible single-nucleotide variant in {n_genes} spliceosomal snRNA loci
({n_snvs} SNVs, GRCh38), from a selection-supervised model trained on gnomAD v4.1 allele counts with each gene's
whole paralog cluster held out (cross-fitted), plus a function score transferred from RNU4-2 saturation genome
editing. Pre-registered analysis plan, code and timestamps: {repo}

## Columns of scores.tsv.gz
- `ID` chr-pos-ref-alt (GRCh38, + strand); `locus` (gene, or family|Ensembl ID for unnamed copies); `gene`; `n` transcript position; `ref_t`/`alt_t` RNA bases
- `S` framework score (higher = more damaging): selection f for dominant genes, function F for recessive genes,
  mean within-gene percentile rank of both for mixed genes; `mode` the gene's inheritance mode used
- `f` cross-fitted selection score (log scale, >= 0); `F` function score on the SGE scale (homology transfer where
  mapped to RNU4-2, otherwise the model's function head); `gamma` the gene gate
- comparators: `cadd_phred`, `phylop447`, `remm`, `gpn_msa` (minus the log-likelihood ratio)
{final_cols}
## Intended use and caveats
- Research use. Evidence labels follow the ClinGen local-LR calibration (Pejaver et al. 2022) at prior 0.10 and
  apply only to the gene they were calibrated in; they are one line of evidence, never a classification.
- RNU4-2 scores are from a fit that never saw the U4 family; the measured SGE scores should be used for RNU4-2 itself.
- Scores were frozen on {date}; variants first published after this date form the pseudo-prospective test set.

License: CC BY 4.0.
"""


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=["scores", "final"], required=True)
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--out", default=os.path.join(HERE, "release"))
    ap.add_argument("--repo", default="https://github.com/tanukusaitejesh-prog/ndd-pp3-calibration (snrna/)")
    args = ap.parse_args()
    if args.stage == "final":
        ri = json.load(open(os.path.join(args.results, "calibration_run_info.json")))
        if not ri.get("stage", "").startswith("FINAL") or ri.get("synthetic"):
            raise SystemExit("final release needs the FINAL calibration (run final_eval.py first)")
    date = time.strftime("%Y-%m-%d", time.gmtime())
    sc = C.gene_scores(args.results, args.inputs, args.evo2)
    sc["locus"] = sc.label
    keep = ["ID", "locus", "gene", "n", "ref_t", "alt_t", "S", "mode", "f", "F", "gamma", "cadd_phred", "phylop447", "remm", "gpn_msa"]
    out = sc[[c for c in keep if c in sc.columns]].copy()
    final_cols = ""
    th_out = None
    if args.stage == "final":
        th = pd.read_csv(os.path.join(args.results, "calibration_thresholds.tsv"), sep="\t")
        th = th[(th.prior == 0.10) & (th.positives == "all") & (th.score == "S") & ~th.unit.str.contains(r"\[|pool")]
        by_gene = {r.unit: r._asdict() for r in th.itertuples()}
        out["evidence_prior_0.10"] = [evidence(r, by_gene.get(r.gene)) for r in out.itertuples()]
        th_out = th
        final_cols = ("- `evidence_prior_0.10` PP3/BP4 strength reached by `S` under that gene's calibrated thresholds "
                      "(empty if none or the gene could not be calibrated)\n")
    # amendment A5: the frozen v4 score (hash-checked), and its post hoc PP3/BP4 label from v4_calibrate.py if run
    v4p, v4m = os.path.join(HERE, "results_v4", "v4_scores.tsv.gz"), os.path.join(HERE, "results_v4", "v4_manifest.txt")
    if os.path.exists(v4p):
        man = {l.split()[1]: l.split()[0] for l in open(v4m) if l.strip()}
        if man.get(os.path.relpath(v4p, HERE)) != sha256(v4p):
            raise SystemExit("results_v4/v4_scores.tsv.gz does not match its frozen manifest")
        out = out.merge(pd.read_csv(v4p, sep="\t")[["ID", "v4"]], on="ID", how="left")
        final_cols += "- `v4` frozen stacked score (amendment A5; higher = more damaging; the paper's headline model)\n"
        v4t = os.path.join(args.results, "v4_calibration_thresholds.tsv")
        if args.stage == "final" and os.path.exists(v4t):
            t4 = pd.read_csv(v4t, sep="\t")
            t4 = t4[(t4.prior == 0.10) & (t4.score == "v4") & (t4.unit == "pooled (all genes)")]
            if len(t4):
                pooled = t4.iloc[0].to_dict()
                out["v4_evidence_prior_0.10"] = [evidence(pd.Series({"S": s}), pooled) for s in out.v4]
                th_out = pd.concat([th_out, t4], ignore_index=True) if th_out is not None else t4
                final_cols += ("- `v4_evidence_prior_0.10` PP3/BP4 strength reached by `v4` under the pooled post hoc "
                               "calibration (v4_calibrate.py)\n")
    name = f"snrna-vep-v2-{args.stage}-{date}"
    d = os.path.join(args.out, name)
    os.makedirs(d, exist_ok=True)
    out.to_csv(os.path.join(d, "scores.tsv.gz"), sep="\t", index=False, float_format="%.5g")
    if th_out is not None:
        th_out.to_csv(os.path.join(d, "thresholds.tsv"), sep="\t", index=False, float_format="%.5g")
    with open(os.path.join(d, "README.md"), "w") as f:
        f.write(README.format(stage=args.stage, date=date, n_genes=out.locus.nunique(), n_snvs=len(out), repo=args.repo,
                              final_cols=final_cols))
    meta = {"metadata": {
        "upload_type": "dataset", "title": f"snRNA-VEP v2: variant-effect scores for spliceosomal snRNA genes ({args.stage}, {date})",
        "description": "Cross-fitted selection-supervised and function-transfer scores for every possible SNV in 55 "
                       "spliceosomal snRNA loci, with ClinGen-calibrated PP3/BP4 thresholds (final release). "
                       "Pre-registered; see the linked repository.",
        "creators": [{"name": "AUTHOR, Name (edit before upload)"}], "license": "cc-by-4.0",
        "keywords": ["snRNA", "spliceosome", "variant effect prediction", "ReNU syndrome", "RNU4-2", "RNU4ATAC",
                     "ClinGen", "PP3", "BP4", "neurodevelopmental disorders"],
        "related_identifiers": [{"identifier": args.repo.split()[0], "relation": "isSupplementTo"}],
        "publication_date": date}}
    with open(os.path.join(d, "zenodo.json"), "w") as f:
        json.dump(meta, f, indent=1)
    with open(os.path.join(d, "SHA256SUMS"), "w") as f:
        for fn in sorted(os.listdir(d)):
            if fn != "SHA256SUMS":
                f.write(f"{sha256(os.path.join(d, fn))}  {fn}\n")
    z = shutil.make_archive(d, "zip", args.out, name)
    print(f"release written: {z} ({os.path.getsize(z) / 1e6:.1f} MB); {len(out)} SNVs, {out.locus.nunique()} loci"
          + (f"; evidence labels for {int((out['evidence_prior_0.10'] != '').sum())} SNVs" if args.stage == "final" else ""))


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
