#!/usr/bin/env python3
"""v3 (exploratory, NOT the pre-registered v2): add allele counts from further public cohorts to the training counts.

    python v3_counts.py cohorts/ukb.tsv [cohorts/bravo.tsv ...]      # -> model_inputs_v3/
    python v3_counts.py cohorts/ukb.csv     # all UK Biobank browser exports concatenated into one CSV (see ukb_regions.txt)

Each cohort file is a TSV with columns chrom, pos, ref, alt, AC, AN (GRCh38, + strand, PASS sites only), one file per
cohort, from cohorts that do not overlap gnomAD v4.1 genomes or each other (state overlap handling in the file name /
README). The training count becomes AC_pass = gnomAD PASS AC + sum of cohort ACs. Everything else (sites, `include`
flags, L3 labels, Evo 2 features) is unchanged, so v3 differs from v2 only in the counts the model learns from.
"""
import json
import os
import shutil
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC, DST = os.path.join(HERE, "model_inputs"), os.path.join(HERE, "model_inputs_v3")


def main(paths):
    if not paths:
        raise SystemExit(__doc__)
    var = pd.read_csv(os.path.join(SRC, "variants.tsv.gz"), sep="\t")
    var["AC_gnomad"] = var.AC_pass
    summ = {"cohorts": {}}
    for p in paths:
        if p.endswith(".csv"):  # UK Biobank Allele Frequency Browser export ("Variant ID", "Allele Count", ...)
            c = pd.read_csv(p).rename(columns={"Variant ID": "ID", "Allele Count": "AC", "Allele Number": "AN"})
        else:
            c = pd.read_csv(p, sep="\t", dtype={"chrom": str})
            c["chrom"] = "chr" + c.chrom.str.replace("chr", "", regex=False)
            c["ID"] = c.chrom + "-" + c.pos.astype(str) + "-" + c.ref + "-" + c.alt
        c = c.drop_duplicates("ID")
        name = os.path.splitext(os.path.basename(p))[0]
        ac = var.ID.map(dict(zip(c.ID, c.AC))).fillna(0).astype(int)
        var[f"AC_{name}"] = ac
        var["AC_pass"] = var.AC_pass + ac
        summ["cohorts"][name] = {"rows": len(c), "matched_snvs": int(var.ID.isin(c.ID).sum()), "total_AC": int(ac.sum()),
                                 "median_AN": float(c.AN.median()) if len(c) else None}
    os.makedirs(DST, exist_ok=True)
    for f in os.listdir(SRC):
        if f != "variants.tsv.gz":
            shutil.copy(os.path.join(SRC, f), DST)
    var.to_csv(os.path.join(DST, "variants.tsv.gz"), sep="\t", index=False)
    summ["observed_snvs_gnomad"] = int((var.AC_gnomad > 0).sum())
    summ["observed_snvs_combined"] = int((var.AC_pass > 0).sum())
    with open(os.path.join(DST, "v3_counts.json"), "w") as f:
        json.dump(summ, f, indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
