#!/usr/bin/env python3
"""Homology transfer of RNU4-2 SGE function scores to RNU4ATAC (amendment A2, upgrade U3). Label-free for RNU4ATAC.

U4atac is the minor-spliceosome counterpart of U4: both form the same 5' stem-loop with a k-turn (bound by
Snu13/15.5K and Prp31 in both spliceosomes) and pair with U6 / U6atac through stems I and II. RNU4-2 and RNU4ATAC
are aligned to both Rfam covariance models (U4 RF00015, U4atac RF00618; Infernal cmalign, match states only). A
RNU4ATAC position is mapped to a RNU4-2 position only where both alignments agree.

Transfer score h(v) for a RNU4ATAC SNV at a mapped position (higher = more damaging, the SGE function-score scale):
  - same reference base in RNU4-2: minus the SGE score of the same alternative base
  - different reference base: minus the mean SGE score of the three SNVs at that RNU4-2 position
Unmapped positions get no h. Writes model_inputs/homology_rnu4atac.tsv.gz and a summary (no patient label is read).

    INFERNAL_BIN=/path/to/infernal/bin python homology_transfer.py
"""
import json
import os
import subprocess
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
INPUTS = os.path.join(HERE, "model_inputs")
CMS = ("RF00015", "RF00618")


def residues(sto):
    seqs = {}
    for line in open(sto):
        if line.startswith(("RNU4-2 ", "RNU4ATAC ")):
            k, v = line.split()
            seqs[k] = seqs.get(k, "") + v
    return seqs


def column_mapping(sto):
    """RNU4ATAC position -> RNU4-2 position for columns where both are match-state residues (upper case)."""
    s = residues(sto)
    i = j = 0
    m = {}
    for x, y in zip(s["RNU4-2"], s["RNU4ATAC"]):
        xa, yb = x not in ".-", y not in ".-"
        i += xa
        j += yb
        if xa and yb and x.isupper() and y.isupper():
            m[j] = i
    return m


def main():
    nodes = pd.read_csv(os.path.join(INPUTS, "nodes.tsv.gz"), sep="\t")
    nodes = nodes[nodes.molecule == "target"]
    seq = {lab: "".join(nodes[nodes.label == lab].sort_values("n").base) for lab in ("RNU4-2", "RNU4ATAC")}
    binp = os.environ.get("INFERNAL_BIN", "")
    maps = {}
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, "u4.fa")
        with open(fa, "w") as f:
            for k, v in seq.items():
                f.write(f">{k}\n{v}\n")
        for cm in CMS:
            sto = os.path.join(tmp, cm + ".sto")
            subprocess.run([os.path.join(binp, "cmalign"), "--notrunc", "-o", sto, os.path.join(HERE, "data", "rfam", cm + ".cm"), fa],
                           check=True, capture_output=True)
            maps[cm] = column_mapping(sto)
    agree = {p: q for p, q in maps[CMS[0]].items() if maps[CMS[1]].get(p) == q}

    var = pd.read_csv(os.path.join(INPUTS, "variants.tsv.gz"), sep="\t")
    sge = pd.read_csv(os.path.join(INPUTS, "sge_eval.tsv.gz"), sep="\t")
    s42 = sge[sge.is_snv].merge(var[["ID", "n", "ref_t", "alt_t"]], on=["ID", "n"])
    by_alt = {(r.n, r.alt_t): r.function_score_mean for r in s42.itertuples()}
    pos_mean = s42.groupby("n").function_score_mean.mean().to_dict()
    ref42 = dict(zip(s42.n, s42.ref_t))
    v4 = var[var.label == "RNU4ATAC"][["ID", "n", "ref_t", "alt_t"]].copy()
    rows = []
    for r in v4.itertuples():
        q = agree.get(r.n)
        if q is None or q not in pos_mean:
            rows.append({"ID": r.ID, "n": r.n, "rnu42_n": None, "same_ref": None, "h": None})
            continue
        same = ref42.get(q) == r.ref_t
        h = -by_alt[(q, r.alt_t)] if same and (q, r.alt_t) in by_alt else -pos_mean[q]
        rows.append({"ID": r.ID, "n": r.n, "rnu42_n": q, "same_ref": same, "h": h})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(INPUTS, "homology_rnu4atac.tsv.gz"), sep="\t", index=False, float_format="%.5g")
    summ = {"cms": list(CMS), "pairs_per_cm": {k: len(v) for k, v in maps.items()}, "agreed_positions": len(agree),
            "agreed_rnu4atac_span": [min(agree), max(agree)] if agree else None,
            "variants_with_h": int(out.h.notna().sum()), "variants": len(out),
            "same_ref_fraction": float(out.same_ref.dropna().astype(bool).mean()) if out.h.notna().any() else None}
    with open(os.path.join(INPUTS, "homology_rnu4atac.json"), "w") as f:
        json.dump(summ, f, indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
