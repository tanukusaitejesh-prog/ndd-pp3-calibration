#!/usr/bin/env python3
"""In-vivo readout of snRNA variants in healthy carriers: 1000 Genomes 30x genotypes + Geuvadis LCL RNA-seq junctions.

Exploratory (added after the v2 endpoints were seen; hypotheses and analysis fixed in this file and committed before
the test was run, see PREREG s14 2026-10-14). Open data only.

Data
  * Geuvadis GD462.JunctionQuantCount (E-GEUV-1; hg19 junctions chr_start_end, normalised counts, 462 LCLs)
  * hg19.2bit (UCSC) to classify introns: U12-type if the 5' splice site matches [AG]TATCC[TC] at intron positions
    +1..+7 AND a U12 branch-point motif TCCTT[AG]A[CT] lies within the last 40 nt; strand from the motifs
  * 1000 Genomes 30x panel (3202 samples, GRCh38) at the 55 analysis-set snRNA loci (model_inputs/genes.tsv.gz)

Per-sample readout (MI score): for each U12 junction with mean count >= 10, log((count + 1) / (total junction count
of the same gene + 1)), z-scored within population (YRI vs the rest); MI = mean z over U12 junctions. U2 control
score: the same over all other junctions of the same genes. Primary readout: MI - U2 (minor-intron-specific usage).

Pre-specified tests (one-sided Mann-Whitney, carriers vs samples carrying no rare variant at any snRNA locus):
  H1  carriers of a rare SNV (1000G panel AF < 0.005) in a minor-spliceosome snRNA (RNU4ATAC, RNU6ATAC, RNU12, RNU11)
      have lower MI - U2 than non-carriers.
  H1b same, restricted to carriers whose variant has CADD >= 15 or a homology-transfer score h (RNU4ATAC) in the
      top half (variants predicted damaging).
  N1  negative control: carriers of rare SNVs in major-spliceosome named snRNAs only show no shift (two-sided).
Reported whatever the outcome; one-sided p for H1 and H1b, two-sided for N1, no multiplicity correction claimed beyond
listing all three. Heterozygous recessive-gene variants are expected to have small effects: a null is plausible.

    python geuvadis_carriers.py --junctions geuv/junc.txt.gz --twobit geuv/hg19.2bit [--out results/geuvadis_carriers]
"""
import argparse
import json
import os
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MINOR = ("RNU4ATAC", "RNU6ATAC", "RNU12", "RNU11")
VCF = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
       "20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{c}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
U12_5SS = re.compile(r"^[AG]TATCC[TC]")
U12_BP = re.compile(r"TCCTT[AG]A[CT]")
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def classify_introns(junc, twobit):
    """Return a Series: 'U12', 'U2' or '' per junction ID (hg19, 1-based intron start/end inferred from GT-AG)."""
    import py2bit
    tb = py2bit.open(twobit)
    chroms = set(tb.chroms())
    out = {}
    for jid, ch, s, e in zip(junc.TargetID, junc.Chr.astype(str), junc.s, junc.e):
        c = "chr" + ch
        if c not in chroms:
            out[jid] = ""
            continue
        seq = tb.sequence(c, s - 1, e).upper()  # intron s..e (1-based, inclusive) if coordinates are intron bounds
        if len(seq) < 60:
            out[jid] = ""
            continue
        plus, minus = seq, rc(seq)
        lab = ""
        for t in (plus, minus):
            if U12_5SS.match(t) and U12_BP.search(t[-40:]):
                lab = "U12"
                break
            if t.startswith("GT") and t.endswith("AG"):
                lab = "U2"
        out[jid] = lab
    return pd.Series(out)


def infer_offsets(junc, twobit, n=3000):
    """Coordinates convention check: which (start, end) offsets make GT..AG (either strand) most frequent."""
    import py2bit
    tb = py2bit.open(twobit)
    sub = junc.sample(min(n, len(junc)), random_state=0)
    best = None
    for ds in (-1, 0, 1):
        for de in (-1, 0, 1):
            ok = 0
            for ch, s, e in zip(sub.Chr.astype(str), sub.s, sub.e):
                try:
                    seq = tb.sequence("chr" + ch, s + ds - 1, e + de).upper()
                except Exception:
                    continue
                ok += (seq[:2] == "GT" and seq[-2:] == "AG") or (seq[:2] == "CT" and seq[-2:] == "AC")
            if best is None or ok > best[0]:
                best = (ok, ds, de)
    return best


def sample_scores(junc, cls, samples, pop):
    x = junc.set_index("TargetID")[samples].astype(float)
    gene_tot = junc.groupby("Gene_Symbol")[samples].transform("sum").set_index(junc.TargetID)
    keep_gene = junc.groupby("Gene_Symbol").TargetID.transform(lambda t: (cls.reindex(t) == "U12").any()).to_numpy()
    j = junc[keep_gene & (x.mean(axis=1).to_numpy() >= 10)]
    lr = np.log((x.loc[j.TargetID] + 1) / (gene_tot.loc[j.TargetID] + 1))
    z = lr.copy()
    for p in np.unique(pop):
        cols = [s for s, q in zip(samples, pop) if q == p]
        sub = lr[cols]
        z[cols] = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    c = cls.reindex(j.TargetID).to_numpy()
    mi, u2 = z[c == "U12"].mean(axis=0), z[c == "U2"].mean(axis=0)
    return pd.DataFrame({"MI": mi, "U2": u2, "MI_minus_U2": mi - u2}), int((c == "U12").sum()), int((c == "U2").sum())


def carriers(samples, genes_tab, variants):
    import pysam
    rows = []
    for _, g in genes_tab[genes_tab.named].iterrows():
        vf = pysam.VariantFile(VCF.format(c=g.chrom))
        panel = list(vf.header.samples)
        idx = [panel.index(s) for s in samples if s in panel]
        for rec in vf.fetch(g.chrom, int(g.start) - 1, int(g.end)):
            if len(rec.ref) != 1 or any(len(a) != 1 for a in rec.alts or ()):
                continue
            gts = [rec.samples[i]["GT"] for i in range(len(panel))]
            for ai, alt in enumerate(rec.alts, start=1):
                ac = sum(sum(1 for a in gt if a == ai) for gt in gts if gt)
                an = sum(sum(1 for a in gt if a is not None) for gt in gts if gt)
                af = ac / an if an else 1
                if af >= 0.005 or ac == 0:
                    continue
                vid = f"{rec.chrom}-{rec.pos}-{rec.ref}-{alt}"
                for i in idx:
                    gt = gts[i]
                    if gt and ai in gt:
                        rows.append({"sample": panel[i], "gene": g.gene, "ID": vid, "panel_AF": af,
                                     "dosage": sum(1 for a in gt if a == ai)})
    car = pd.DataFrame(rows)
    if len(car):
        car = car.merge(variants[["ID", "cadd_phred", "ddG_bind_max"]], on="ID", how="left")
    return car


def mw(a, b, alternative):
    from scipy.stats import mannwhitneyu
    if len(a) == 0 or len(b) == 0:
        return {"n": int(len(a)), "p": None, "median_diff": None}
    r = mannwhitneyu(a, b, alternative=alternative)
    return {"n": int(len(a)), "n_ref": int(len(b)), "p": float(r.pvalue), "median_diff": float(np.median(a) - np.median(b))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--twobit", required=True)
    ap.add_argument("--out", default=os.path.join(HERE, "results", "geuvadis_carriers"))
    a = ap.parse_args()
    junc = pd.read_csv(a.junctions, sep="\t")
    se = junc.TargetID.str.split("_", expand=True)
    junc["s"], junc["e"] = se[1].astype(int), se[2].astype(int)
    ok, ds, de = infer_offsets(junc, a.twobit)
    junc["s"], junc["e"] = junc.s + ds, junc.e + de
    cls = classify_introns(junc, a.twobit)
    samples = [c for c in junc.columns if re.match(r"^(HG|NA)\d+$", c)]
    yri = set(s for s in samples if s.startswith("NA18") or s.startswith("NA19"))  # Geuvadis YRI are NA18xxx/NA19xxx
    pop = np.array(["YRI" if s in yri else "EUR" for s in samples])
    sc, n_u12, n_u2 = sample_scores(junc, cls, samples, pop)
    genes = pd.read_csv(os.path.join(HERE, "model_inputs", "genes.tsv.gz"), sep="\t")
    variants = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    car = carriers(samples, genes, variants)
    hom = os.path.join(HERE, "model_inputs", "homology_rnu4atac.tsv.gz")
    if len(car) and os.path.exists(hom):
        h = pd.read_csv(hom, sep="\t")
        hcol = [c for c in h.columns if c.startswith("h")][:1]
        if hcol and "ID" in h.columns:
            car = car.merge(h[["ID", hcol[0]]].rename(columns={hcol[0]: "h"}), on="ID", how="left")
    if "h" not in car:
        car["h"] = np.nan
    any_car = set(car["sample"]) if len(car) else set()
    ref = sc.loc[[s for s in samples if s not in any_car], "MI_minus_U2"].to_numpy()
    minor_c = car[car.gene.isin(MINOR)] if len(car) else car
    hmed = car.h.median()
    dmg = minor_c[(minor_c.cadd_phred >= 15) | (minor_c.h >= hmed)] if len(minor_c) else minor_c
    major_only = sorted(any_car - set(minor_c["sample"])) if len(car) else []
    res = {"coordinate_offsets": {"start": ds, "end": de, "gt_ag_in_sample": ok},
           "n_junctions_U12": n_u12, "n_junctions_U2_same_genes": n_u2, "n_samples": len(samples),
           "n_carrier_rows": int(len(car)),
           "H1_minor_carriers": mw(sc.loc[sorted(set(minor_c["sample"])), "MI_minus_U2"].to_numpy(), ref, "less"),
           "H1b_minor_predicted_damaging": mw(sc.loc[sorted(set(dmg["sample"])), "MI_minus_U2"].to_numpy(), ref, "less"),
           "N1_major_only_carriers": mw(sc.loc[major_only, "MI_minus_U2"].to_numpy(), ref, "two-sided")}
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    sc.to_csv(a.out + "_scores.tsv", sep="\t", float_format="%.4g")
    car.to_csv(a.out + "_carriers.tsv", sep="\t", index=False)
    if len(minor_c):
        per = minor_c.merge(sc, left_on="sample", right_index=True)
        per.to_csv(a.out + "_minor_carriers.tsv", sep="\t", index=False, float_format="%.4g")
    with open(a.out + ".json", "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
