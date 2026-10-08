"""Cached per-submission laboratory dataset for the v14 analyses (review 20 roadmap).

For each ClinVar vintage: germline, non-expert-panel SCVs whose comment lists any ACMG/AMP code, on missense SNVs
(variant_summary 2026-10 annotation, GRCh38). Codes are parsed with the published parser of make_lab_stacking.py
(unchanged; strengths, negations skipped). Stores every code's points, the reported class, submitter, gene, SCV,
evaluation date and the comment text (needed for PM1 rationale and templating analyses).
Usage: python make_lab_dataset.py 2026-10 [2024-01 2022-01 2020-01]  -> work/cvt/lab_codes_<vintage>.tsv.gz"""
import json
import re
import sys

import pandas as pd

src = open("make_lab_stacking.py").read().split('path = CVT + "submission_summary_2026-10.txt.gz"')[0]
G = {"__name__": "lab_dataset"}
exec(compile(src, "make_lab_stacking_head", "exec"), G)
parse, CVT, EXCL = G["parse"], G["CVT"], G["EXCL"]

vs = pd.read_csv(CVT + "variant_summary_2026-10.txt.gz", sep="\t", dtype=str, compression="gzip",
                 usecols=["VariationID", "Name", "GeneSymbol", "Type", "Assembly"])
vs = vs[(vs.Assembly == "GRCh38") & (vs.Type == "single nucleotide variant")].drop_duplicates("VariationID")
MIS = re.compile(r"\(p\.(?:[A-Z][a-z]{2})\d+(?!Ter|=|del|ins|dup|fs|ext)[A-Z][a-z]{2}\)")
vs = vs[vs.Name.str.contains(MIS, na=False)][["VariationID", "GeneSymbol", "Name"]]

for v in sys.argv[1:]:
    path = CVT + f"submission_summary_{v}.txt.gz"
    import gzip
    with gzip.open(path, "rt") as fh:
        n = next(i for i, line in enumerate(fh) if line.startswith("#VariationID\t"))
    s = pd.read_csv(path, sep="\t", skiprows=n, header=0, dtype=str, quoting=3, compression="gzip", on_bad_lines="skip",
                    usecols=["#VariationID", "ClinicalSignificance", "Description", "OriginCounts", "ReviewStatus", "Submitter",
                             "SubmittedGeneSymbol", "DateLastEvaluated", "SCV"]).rename(columns={"#VariationID": "VariationID"})
    s = s[s.OriginCounts.fillna("").str.contains("germline", case=False)]
    s = s[~s.ReviewStatus.fillna("").str.contains("expert panel|practice guideline", case=False)]
    s = s[s.Description.fillna("").str.contains(r"\b(?:PVS|PS|PM|PP|BA|BS|BP)[1-9]\b", regex=True)]
    s = s.merge(vs, on="VariationID")
    ev = s.Description.map(parse)
    s["codes"] = ev.map(json.dumps)
    for c in ["PP3", "PM1", "PP2", "PS2", "PM2", "PM6", "PP1", "PP4", "PS3", "PM5", "PS1", "BP4"]:
        s[c] = ev.map(lambda e, c=c: e.get(c, 0))
    s["indep"] = ev.map(lambda e: sum(x for k, x in e.items() if k not in EXCL))
    s["total"] = ev.map(lambda e: sum(e.values()))
    s["gene"] = s.GeneSymbol
    s["vintage"] = v
    keep = ["VariationID", "SCV", "Submitter", "gene", "Name", "ClinicalSignificance", "ReviewStatus", "DateLastEvaluated",
            "vintage", "codes", "PP3", "PM1", "PP2", "PS2", "PM2", "PM6", "PP1", "PP4", "PS3", "PM5", "PS1", "BP4", "indep",
            "total", "Description"]
    s[keep].to_csv(CVT + f"lab_codes_{v}.tsv.gz", sep="\t", index=False, compression="gzip")
    print(f"{v}: {len(s)} SCVs, {s.Submitter.nunique()} submitters, PP3&PM1 {int(((s.PP3 > 0) & (s.PM1 > 0)).sum())}", flush=True)
