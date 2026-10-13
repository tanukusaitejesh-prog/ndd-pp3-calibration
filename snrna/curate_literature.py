#!/usr/bin/env python3
"""Literature tally -> labels/L2_curated_additions.tsv (PREREG s4 L2: variants reported as P/LP or disease-causing by
the reporting study; amendment A3). Non-blind genes only: any RNU2-2 or RNU5B-1 row stops the script.

    python curate_literature.py lit/*.tsv

Each SNV given as transcript HGVS (n.<pos><ref>><alt>) is mapped to GRCh38 with the analysis-set transcript
coordinates (model_inputs/genes.tsv.gz); the reference base must match the transcript base at that position, else
the row is rejected and listed. A genomic GRCh38 coordinate given by the source must agree with the mapping.
Indels are kept in labels/L2_curated_indels.tsv (not scored by the model). Output columns follow L2_nonblind.tsv,
plus the source's in-silico/conservation flag, zygosity, phenotype and first publication date.
"""
import os
import re
import sys
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SEALED = ("RNU2-2", "RNU5B-1")
HGVS_SNV = re.compile(r"^n\.(\d+)([ACGTU])>([ACGTU])$")
COMP = str.maketrans("ACGT", "TGCA")
POSITIVE_WORDS = ("pathogenic", "likely pathogenic", "disease-causing", "disease causing", "causal", "causative")


def to_dna(b):
    return b.replace("U", "T")


def main(paths):
    genes = pd.read_csv(os.path.join(HERE, "model_inputs", "genes.tsv.gz"), sep="\t")
    nodes = pd.read_csv(os.path.join(HERE, "model_inputs", "nodes.tsv.gz"), sep="\t")
    t = nodes[nodes.molecule == "target"]
    tx = {lab: dict(zip(x.n, x.base)) for lab, x in t.groupby("label")}
    named = genes[genes.named].set_index("gene")
    rows, indels, rejected = [], [], []
    for p in paths:
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        for r in df.to_dict("records"):
            gene = r["gene"].strip()
            if gene in SEALED:
                raise SystemExit(f"{p}: contains a sealed-gene row; remove the file and re-run the tally")
            cls = r.get("classification_as_stated", "").lower()
            if not any(w in cls for w in POSITIVE_WORDS) or "uncertain" in cls or "benign" in cls:
                rejected.append({**r, "reason": "not reported as P/LP or disease-causing"})
                continue
            if gene not in named.index:
                rejected.append({**r, "reason": "gene not in the analysis set"})
                continue
            g = named.loc[gene]
            hg = r.get("hgvs_n", "").strip().replace(" ", "")
            m = HGVS_SNV.match(hg)
            if not m:
                indels.append({**r, "file": os.path.basename(p)})
                continue
            n, ref, alt = int(m.group(1)), m.group(2).replace("T", "U"), m.group(3).replace("T", "U")
            base = tx[g.label].get(n)
            if base != ref:
                rejected.append({**r, "reason": f"reference {ref} at n.{n} but the transcript has {base}"})
                continue
            pos = g.start + n - 1 if g.strand == "+" else g.end - n + 1
            gref, galt = (to_dna(ref), to_dna(alt)) if g.strand == "+" else (to_dna(ref).translate(COMP), to_dna(alt).translate(COMP))
            vid = f"{g.chrom}-{pos}-{gref}-{galt}"
            given = r.get("genomic_grch38", "").replace("chr", "").replace(":", "-").strip()
            if given and given.split("-")[:2] != [g.chrom.replace("chr", ""), str(pos)]:
                rejected.append({**r, "reason": f"source GRCh38 {r['genomic_grch38']} disagrees with the mapping {vid}"})
                continue
            inh = r.get("inheritance", "").strip().lower()
            rows.append({"ID": vid, "gene": gene, "source": f"literature:{r.get('source', '')}",
                         "classification": r.get("classification_as_stated", ""),
                         "review_status": "literature (ACMG given)" if r.get("acmg_classification_given", "").lower() == "yes"
                         else "literature", "inheritance": inh if inh in ("dominant", "recessive") else "unknown",
                         "title": f"{r.get('transcript', '')}({gene}):{hg}", "label": "P/LP",
                         "insilico_or_conservation_used_for_classification": r.get("insilico_or_conservation_used_for_classification", "unknown"),
                         "zygosity": r.get("zygosity_observed", ""), "phenotype": r.get("phenotype", ""),
                         "first_publication_date": r.get("first_publication_date", "") or r.get("first_publication_year", ""),
                         "tally_date": time.strftime("%Y-%m-%d", time.gmtime())})
    out = pd.DataFrame(rows).drop_duplicates("ID") if rows else pd.DataFrame()
    os.makedirs(os.path.join(HERE, "labels"), exist_ok=True)
    out.to_csv(os.path.join(HERE, "labels", "L2_curated_additions.tsv"), sep="\t", index=False)
    pd.DataFrame(indels).to_csv(os.path.join(HERE, "labels", "L2_curated_indels.tsv"), sep="\t", index=False)
    pd.DataFrame(rejected).to_csv(os.path.join(HERE, "labels", "L2_curated_rejected.tsv"), sep="\t", index=False)
    old = pd.read_csv(os.path.join(HERE, "labels", "L2_nonblind.tsv"), sep="\t")
    new = out[~out.ID.isin(old[old.label == "P/LP"].ID)] if len(out) else out
    print(f"SNVs accepted {len(out)} ({len(new)} not already in L2_nonblind); indels {len(indels)}; rejected {len(rejected)}")
    if len(out):
        print(out.groupby("gene").ID.nunique().to_string())
        print("new by gene:\n" + (new.groupby("gene").ID.nunique().to_string() if len(new) else "none"))


if __name__ == "__main__":
    main(sys.argv[1:])
