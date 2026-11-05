#!/usr/bin/env python3
"""Post hoc mechanism check, written before running (5 Nov 2026; PREREG s14).

alignment_depth.py found phyloP447 flat and low across RNU4-2 (max 0.90) but not at its paralog RNU4-1 (8.03), and
proposed that the 447-way whole-genome alignment captures few species at RNU4-2. That was never checked against the
alignment itself. This script reads the UCSC Cactus 447-way alignment (hg38.cactus447way.bb, bigMaf) for each named
disease gene and its paralog RNU4-1, using HTTP range requests (no full download), and reports per gene:
  - aligned species per human position (rows with a non-gap base in that column): median and minimum;
  - number of distinct species rows with >= 50% of the gene's positions aligned;
  - mean identity to the human base among aligned rows.
Prediction stated in advance: RNU4-2 has far fewer aligned species per position than RNU4-1 and the other genes.
If not, the alignment-depth explanation is wrong and must be withdrawn. Reported whatever the result.

    python alignment_species.py   # -> results/alignment_species.tsv
"""
import os
import struct
import zlib

import numpy as np
import pandas as pd
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus447way/hg38.cactus447way.bb"
CA = "/root/.ccr/ca-bundle.crt" if os.path.exists("/root/.ccr/ca-bundle.crt") else True
GENES = ["RNU4-2", "RNU4-1", "RNU2-2", "RNU12", "RNU4ATAC", "RNU6ATAC", "RNU5B-1"]
S = requests.Session()


def rng(off, n):
    r = S.get(URL, headers={"Range": f"bytes={off}-{off + n - 1}"}, verify=CA, timeout=120)
    r.raise_for_status()
    return r.content


class BigBed:
    def __init__(self):
        h = rng(0, 64)
        (magic, ver, zl, self.ctree, self.data, self.index, fc, dfc, aso, tso, self.ubuf, _) = struct.unpack("<IHHQQQHHQQIQ", h)
        assert magic == 0x8789F2EB
        self.chroms = self._chroms()

    def _chroms(self):
        h = rng(self.ctree, 32)
        magic, bs, ks, vs, n, _ = struct.unpack("<IIIIQQ", h)
        out = {}

        def node(off):
            hdr = rng(off, 4)
            leaf, _, cnt = struct.unpack("<BBH", hdr)
            item = ks + (8 if leaf else 8)
            body = rng(off + 4, cnt * item)
            for i in range(cnt):
                b = body[i * item:(i + 1) * item]
                key = b[:ks].rstrip(b"\0").decode()
                if leaf:
                    cid, size = struct.unpack("<II", b[ks:ks + 8]); out[key] = cid
                else:
                    node(struct.unpack("<Q", b[ks:ks + 8])[0])
        node(self.ctree + 32)
        return out

    def entries(self, chrom, start, end):
        cid = self.chroms[chrom]
        blocks = []

        def node(off):
            leaf, _, cnt = struct.unpack("<BBH", rng(off, 4))
            item = 32 if leaf else 24
            body = rng(off + 4, cnt * item)
            for i in range(cnt):
                b = body[i * item:(i + 1) * item]
                sc, sb, ec, eb = struct.unpack("<IIII", b[:16])
                if (sc, sb) < (cid, end) and (ec, eb) > (cid, start):
                    if leaf:
                        blocks.append(struct.unpack("<QQ", b[16:32]))
                    else:
                        node(struct.unpack("<Q", b[16:24])[0])
        node(self.index + 48)
        out = []
        for off, size in blocks:
            raw = rng(off, size)
            raw = zlib.decompress(raw) if self.ubuf else raw
            p = 0
            while p < len(raw):
                c, s, e = struct.unpack("<III", raw[p:p + 12]); p += 12
                z = raw.index(b"\0", p); rest = raw[p:z].decode(); p = z + 1
                if c == cid and s < end and e > start:
                    out.append((s, e, rest))
        return out


def gene_stats(bb, chrom, start, end):
    """start/end 0-based half-open human interval."""
    cover = {pos: [] for pos in range(start, end)}   # pos -> list of (src, base)
    human = {}
    for s, e, txt in bb.entries(chrom, start, end):
        rows = [l.split() for l in txt.split(";") if l.startswith("s ")]
        h = [r for r in rows if r[1].startswith("hg38.")]
        if not h:
            continue
        h = h[0]
        hstart, htext = int(h[2]), h[6]
        cols = []
        p = hstart
        for j, ch in enumerate(htext):
            if ch != "-":
                if start <= p < end:
                    cols.append((j, p, ch.upper()))
                p += 1
        for r in rows:
            if r is h:
                continue
            sp = r[1].split(".")[0]
            for j, p, hb in cols:
                b = r[6][j].upper()
                if b != "-":
                    cover[p].append((sp, b == hb))
        for j, p, hb in cols:
            human[p] = hb
    n_aligned = np.array([len({sp for sp, _ in cover[p]}) for p in range(start, end)])
    ident = [np.mean([ok for _, ok in cover[p]]) for p in range(start, end) if cover[p]]
    per_sp = {}
    for p in range(start, end):
        for sp in {sp for sp, _ in cover[p]}:
            per_sp[sp] = per_sp.get(sp, 0) + 1
    L = end - start
    return {"positions": L, "human_positions_in_alignment": len(human),
            "aligned_species_median": float(np.median(n_aligned)), "aligned_species_min": int(n_aligned.min()),
            "aligned_species_max": int(n_aligned.max()),
            "species_with_half_aligned": int(sum(1 for v in per_sp.values() if v >= 0.5 * L)),
            "mean_identity_to_human": float(np.mean(ident)) if ident else np.nan}


def main():
    v = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t", usecols=["gene", "chrom", "pos"])
    bb = BigBed()
    rows = []
    for g in GENES:
        x = v[v.gene == g]
        chrom = str(x.chrom.iloc[0]); chrom = chrom if chrom.startswith("chr") else "chr" + chrom
        st = gene_stats(bb, chrom, int(x.pos.min()) - 1, int(x.pos.max()))
        st.update({"gene": g, "chrom": chrom, "start1": int(x.pos.min()), "end1": int(x.pos.max())})
        rows.append(st)
        print(st, flush=True)
    out = pd.DataFrame(rows)[["gene", "chrom", "start1", "end1", "positions", "human_positions_in_alignment",
                              "aligned_species_median", "aligned_species_min", "aligned_species_max",
                              "species_with_half_aligned", "mean_identity_to_human"]]
    out.to_csv(os.path.join(HERE, "results", "alignment_species.tsv"), sep="\t", index=False, float_format="%.3f")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
