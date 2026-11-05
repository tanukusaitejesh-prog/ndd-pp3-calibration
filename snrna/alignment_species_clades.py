#!/usr/bin/env python3
"""Descriptive extension of alignment_species.py (written 5 Nov 2026 AFTER its counts were seen; PREREG s14): which
species are aligned at each gene? Splits species with >= 50% of the gene's positions aligned into primates (the clade
of hg38 and the mouse lemur, Microcebus murinus, in UCSC's hg38.447way.nh.txt) and non-primates.

    python alignment_species_clades.py   # -> results/alignment_species_clades.tsv
"""
import os
import re

import pandas as pd

import alignment_species as A

HERE = os.path.dirname(os.path.abspath(__file__))


def parse_newick(s):
    s = re.sub(r"\s+", "", s).rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        kids = []
        if s[pos] == "(":
            pos += 1
            while True:
                kids.append(node())
                if s[pos] == ",":
                    pos += 1; continue
                pos += 1; break
        m = re.match(r"([^:,()]*)(:[0-9.eE-]+)?", s[pos:])
        pos += m.end()
        return (m.group(1), kids)
    return node()


def leaves(n):
    return [n[0]] if not n[1] else [x for k in n[1] for x in leaves(k)]


def primate_set(tree, a="hg38", b="Microcebus_murinus"):
    best = None

    def walk(n):
        nonlocal best
        L = leaves(n)
        if a in L and b in L and (best is None or len(L) < len(best)):
            best = L
        for k in n[1]:
            walk(k)
    walk(tree)
    return set(best)


def main():
    tree = parse_newick(open(os.path.join(HERE, "data_ext", "hg38.447way.nh.txt")).read())
    allsp = set(leaves(tree)); prim = primate_set(tree)
    v = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t", usecols=["gene", "chrom", "pos"])
    bb = A.BigBed(); rows = []
    for g in A.GENES:
        x = v[v.gene == g]; chrom = str(x.chrom.iloc[0]); chrom = chrom if chrom.startswith("chr") else "chr" + chrom
        start, end = int(x.pos.min()) - 1, int(x.pos.max())
        cnt = {}
        for s, e, txt in bb.entries(chrom, start, end):
            rws = [l.split() for l in txt.split(";") if l.startswith("s ")]
            h = [r for r in rws if r[1].startswith("hg38.")]
            if not h: continue
            h = h[0]; p = int(h[2]); cols = []
            for j, ch in enumerate(h[6]):
                if ch != "-":
                    if start <= p < end: cols.append(j)
                    p += 1
            for r in rws:
                if r is h: continue
                sp = r[1].split(".")[0]
                cnt[sp] = cnt.get(sp, 0) + sum(1 for j in cols if r[6][j] != "-")
        L = end - start
        al = {sp for sp, c in cnt.items() if c >= 0.5 * L}
        rows.append({"gene": g, "aligned_species": len(al), "primates": len(al & prim),
                     "non_primates": len(al - prim), "primates_in_tree": len(prim & allsp) - 1,
                     "non_primates_in_tree": len(allsp - prim)})
        print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(HERE, "results", "alignment_species_clades.tsv"), sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
