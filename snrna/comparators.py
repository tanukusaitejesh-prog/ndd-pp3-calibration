#!/usr/bin/env python3
"""Additional comparators for every analysis-set SNV (amendment A3; reported only, never gating: the frozen bars
compare against CADD). Higher = more damaging for every column.

  phylop100   phyloP, 100 vertebrates (UCSC track phyloP100way)             position-level conservation
  phylop447   phyloP, 447 mammals (UCSC track phyloP447way, Zoonomia)        position-level conservation
  remm        ReMM v0.4 hg38 (regulatory Mendelian mutation score)           position-level, non-coding
  gpn_msa     minus the GPN-MSA log-likelihood ratio (alt vs ref; hg38)       allele-level, zero-shot

    python comparators.py          # -> comparators/comparators.tsv.gz (cached under data/cache/comparators)
"""
import json
import os
import time
import urllib.request

import pandas as pd

import snrna as S

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "comparators")
CACHE = os.path.join(S.DATA, "cache", "comparators")
REMM = "https://kircherlab.bihealth.org/download/ReMM/ReMM.v0.4.hg38.tsv.gz"
GPN = "https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz"
UCSC = "https://api.genome.ucsc.edu/getData/track?genome=hg38;track={track};chrom={chrom};start={start};end={end}"


def cached(key, fn):
    os.makedirs(CACHE, exist_ok=True)
    p = os.path.join(CACHE, key + ".json")
    if os.path.exists(p):
        return json.load(open(p))
    x = fn()
    with open(p, "w") as f:
        json.dump(x, f)
    return x


def ucsc_track(track, chrom, start, end):
    """{1-based position: value} from the UCSC REST API (start/end 1-based inclusive)."""
    def get():
        url = UCSC.format(track=track, chrom=chrom, start=start - 1, end=end)
        for i in range(6):
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    d = json.load(r)
                break
            except Exception:
                time.sleep(2 * (i + 1))
        else:
            raise RuntimeError(url)
        vals = d.get(track, [])
        if isinstance(vals, dict):
            vals = vals.get(chrom, [])
        out = {}
        for v in vals:
            for p in range(v["start"] + 1, v["end"] + 1):
                out[str(p)] = v["value"]
        return out
    return {int(k): v for k, v in cached(f"{track}_{chrom}_{start}_{end}", get).items()}


def tabix_rows(url, chrom, start, end):
    import pysam

    def get():
        with S.remote_index_dir():
            tb = pysam.TabixFile(url, index=url + ".tbi")
            return [line.split("\t") for line in tb.fetch(chrom, start - 1, end)]
    return cached(f"{os.path.basename(url)}_{chrom}_{start}_{end}", get)


def main():
    var = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    genes = pd.read_csv(os.path.join(HERE, "model_inputs", "genes.tsv.gz"), sep="\t")
    parts = []
    for g in genes.itertuples():
        v = var[var.label == g.label][["ID", "label", "gene", "chrom", "pos", "ref", "alt"]].copy()
        lo, hi = int(v.pos.min()), int(v.pos.max())
        p100 = ucsc_track("phyloP100way", g.chrom, lo, hi)
        p447 = ucsc_track("phyloP447way", g.chrom, lo, hi)
        remm = {int(r[1]): float(r[2]) for r in tabix_rows(REMM, g.chrom, lo, hi)}
        gpn = {(int(r[1]), r[2], r[3]): float(r[4]) for r in tabix_rows(GPN, g.chrom.replace("chr", ""), lo, hi)}
        v["phylop100"] = v.pos.map(p100)
        v["phylop447"] = v.pos.map(p447)
        v["remm"] = v.pos.map(remm)
        v["gpn_msa"] = [-gpn[k] if k in gpn else float("nan") for k in zip(v.pos, v.ref, v.alt)]
        parts.append(v)
        print(f"{g.label}: {len(v)} SNVs; coverage phyloP100 {v.phylop100.notna().mean():.2f}, phyloP447 "
              f"{v.phylop447.notna().mean():.2f}, ReMM {v.remm.notna().mean():.2f}, GPN-MSA {v.gpn_msa.notna().mean():.2f}", flush=True)
    out = pd.concat(parts, ignore_index=True)
    os.makedirs(OUT, exist_ok=True)
    out.drop(columns=["chrom", "pos", "ref", "alt"]).to_csv(os.path.join(OUT, "comparators.tsv.gz"), sep="\t", index=False,
                                                           float_format="%.5g")
    cov = {c: float(out[c].notna().mean()) for c in ("phylop100", "phylop447", "remm", "gpn_msa")}
    with open(os.path.join(OUT, "comparators.json"), "w") as f:
        json.dump({"n_snvs": len(out), "coverage": cov, "sources": {"phyloP100way": "UCSC API", "phyloP447way": "UCSC API",
                                                                    "ReMM": REMM, "GPN-MSA": GPN},
                   "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, f, indent=1)
    print("coverage:", cov)


if __name__ == "__main__":
    main()
