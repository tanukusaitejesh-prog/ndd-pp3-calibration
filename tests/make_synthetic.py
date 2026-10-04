#!/usr/bin/env python3
"""
Synthetic inputs with a known answer, for the gate self-test.

Writes under OUTDIR:
  data/            files with the same names and formats as the real inputs:
                   DeNovoWEST de novo list (synthetic 'GRCh37'), gene table, DDG2P,
                   AlphaMissense_hg38-like table, GENCODE-like GTF, UCSC chain file
  roulette_full/   Roulette-like rates for every position (bcftools-query output format)
  roulette_vcf/    VCF header used by selftest.sh to build Roulette-style VCFs for the bcftools path
  truth/           ground truth: GRCh38 coordinates of every de novo, the coding universe,
                   expected counts and the true likelihood ratios

This reference implementation is written independently of gate.py (string codons,
TCAG-ordered genetic code, per-transcript loops), so agreement is evidence both are right.

Ground truth: background de novos occur at rate k * mu (centre-specific detection);
pathogenic de novo missense occur only in 'DD genes' at rate A * mu * w_gene * AM^3.
"""
import argparse
import gzip
import json
import os
import shutil

import numpy as np
import pandas as pd

TCAG = "TCAG"
AA_TCAG = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {x + y + z: AA_TCAG[16 * i + 4 * j + k]
        for i, x in enumerate(TCAG) for j, y in enumerate(TCAG) for k, z in enumerate(TCAG)}
SENSE = sorted(c for c, aa in CODE.items() if aa != "*")
CMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
L = 700_000
CHROMS = ["21", "22"]
EDGES = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 98, 99, 100]
CENTRES = {"GeneDx": (12000, 1.00), "DDD": (7000, 0.90), "RUMC": (1500, 0.80)}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def sig4(x):
    """Round to 4 significant digits (what '%.4g' writes)."""
    x = np.asarray(x, dtype=float)
    out = np.full_like(x, np.nan)
    ok = np.isfinite(x) & (x > 0)
    e = np.floor(np.log10(x[ok]))
    out[ok] = np.round(x[ok] / 10 ** e, 3) * 10 ** e
    return out


def coarse(am):
    fb = np.clip(np.floor(np.round(np.asarray(am, float) * 100, 6)), 0, 99).astype(int)
    return np.searchsorted(EDGES, fb, side="right") - 1


def layout(rng, start, strand, n_codons, n_exons):
    u5, u3 = int(rng.integers(30, 121)), int(rng.integers(30, 151))
    lc = 3 * n_codons
    mlen = u5 + lc + 3 + u3
    s_lo, s_hi = u5 + lc, u5 + lc + 3
    while True:
        cuts = np.sort(rng.choice(np.arange(25, mlen - 25), size=n_exons - 1, replace=False))
        b = np.r_[0, cuts, mlen]
        if np.all(np.diff(b) >= 25) and not np.any((cuts > s_lo) & (cuts < s_hi)):
            break
    exl = np.diff(b)
    introns = rng.integers(150, 1500, size=n_exons - 1)
    asc_l = list(exl) if strand == "+" else list(exl[::-1])
    asc_i = list(introns) if strand == "+" else list(introns[::-1])
    blocks, g = [], start
    for k, el in enumerate(asc_l):
        blocks.append((g, g + el - 1))
        g += el + (asc_i[k] if k < len(asc_i) else 0)
    if strand == "-":
        blocks = blocks[::-1]
    m2g = np.empty(mlen, dtype=np.int64)
    a = 0
    for (gs, ge), el in zip(blocks, exl):
        m2g[a:a + el] = np.arange(gs, ge + 1) if strand == "+" else np.arange(ge, gs - 1, -1)
        a += el
    return dict(u5=u5, lc=lc, mlen=mlen, blocks=blocks, exl=list(exl), m2g=m2g, end=g)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--biased", action="store_true",
                    help="miss 40%% of de novos at codons with mean AM > 0.75 (null check must fail)")
    ap.add_argument("--sim-seed", type=int, default=777, help="seed for the de novo simulation only")
    ap.add_argument("--cpg-factor", type=float, default=1.0,
                    help="true de novo rate of CpG transitions relative to the rate model (mutation-class miscalibration)")
    a = ap.parse_args()
    rng = np.random.default_rng(12345)
    out = a.out
    for sub in ["data/DeNovoWEST/input", "data/roulette", "roulette_full", "roulette_vcf", "truth"]:
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    D = os.path.join(out, "data")

    genome = {c: rng.choice(np.array(list("ACGT")), size=L, p=[0.29, 0.21, 0.21, 0.29]) for c in CHROMS}

    # -------------------------------------------------------------- transcripts
    specials = {"ISTOP": "internal stop codon", "NF": "cds_start_NF tag", "NOGTF": "absent from GTF",
                "SHIFT": "AlphaMissense protein positions shifted", "DUPA": "duplicate annotation A",
                "DUPB": "duplicate annotation B", "AMGAP": "AlphaMissense lacks 0.5% of rows",
                "RGAP": "one CDS base missing from Roulette"}
    txs = []
    n_per = {"21": 48, "22": 48}
    special_slots = {"22": ["ISTOP", "NF", "NOGTF", "SHIFT", "DUPA", "AMGAP", "RGAP"]}
    gi = 0
    for c in CHROMS:
        pos = 5000
        kinds = ["normal"] * n_per[c] + special_slots.get(c, [])
        rng.shuffle(kinds)
        for kind in kinds:
            gi += 1
            strand = "+" if rng.random() < 0.5 else "-"
            n_codons = int(rng.integers(200, 700)) if kind != "NF" else 300
            n_exons = int(rng.integers(2, 9))
            lay = layout(rng, pos, strand, n_codons, n_exons)
            if lay["end"] > L - 12000:
                break
            name = f"SYN{gi:02d}" if kind == "normal" else f"SYN{gi:02d}_{kind}"
            t = dict(chrom=c, strand=strand, kind=kind, gene_name=name, n_codons=n_codons,
                     gene_id=f"ENSGS{gi:011d}.1", transcript_id=f"ENSTS{gi:011d}.{1 + gi % 3}", **lay)
            txs.append(t)
            if kind == "DUPA":
                gi += 1
                t2 = dict(t)
                t2.update(kind="DUPB", gene_name=f"SYN{gi:02d}_DUPB", gene_id=f"ENSGS{gi:011d}.1",
                          transcript_id=f"ENSTS{gi:011d}.1")
                txs.append(t2)
            pos = lay["end"] + int(rng.integers(1000, 3000))

    # write sequences
    for t in txs:
        if t["kind"] == "DUPB":
            continue
        cds = ["ATG"] + [SENSE[i] for i in rng.integers(0, len(SENSE), size=t["n_codons"] - 1)]
        if t["kind"] == "ISTOP":
            cds[40] = "TGA"
        stop = ["TAA", "TAG", "TGA"][int(rng.integers(0, 3))]
        u5 = "".join(rng.choice(list("ACGT"), size=t["u5"]))
        u3 = "".join(rng.choice(list("ACGT"), size=t["mlen"] - t["u5"] - t["lc"] - 3))
        mrna = u5 + "".join(cds) + stop + u3
        g = genome[t["chrom"]]
        for m, base in enumerate(mrna):
            g[t["m2g"][m]] = base if t["strand"] == "+" else CMP[base]

    # -------------------------------------------------------------- GTF
    lines = ["##description: synthetic GENCODE-like annotation for the gate self-test",
             "##provider: GENCODE", "##contact: none", "##format: gtf", "##date: 2026-10-04"]

    def seg_of(t, lo, hi):
        res, a = [], 0
        for k, ((gs, ge), el) in enumerate(zip(t["blocks"], t["exl"])):
            s, e = max(a, lo), min(a + el, hi)
            if s < e:
                gg = t["m2g"][s:e]
                res.append((k, int(gg.min()) + 1, int(gg.max()) + 1, s - lo))
            a += el
        return res

    for t in txs:
        if t["kind"] == "NOGTF":
            continue
        c = "chr" + t["chrom"]
        tags = 'tag "basic";' + (' tag "cds_start_NF";' if t["kind"] == "NF" else "")
        base = (f'gene_id "{t["gene_id"]}"; transcript_id "{t["transcript_id"]}"; gene_type "protein_coding"; '
                f'gene_name "{t["gene_name"]}"; transcript_type "protein_coding"; transcript_name "{t["gene_name"]}-201";')
        lo = min(b[0] for b in t["blocks"]) + 1
        hi = max(b[1] for b in t["blocks"]) + 1
        st = t["strand"]
        lines.append(f'{c}\tSYN\tgene\t{lo}\t{hi}\t.\t{st}\t.\tgene_id "{t["gene_id"]}"; gene_type "protein_coding"; gene_name "{t["gene_name"]}"; level 2;')
        lines.append(f"{c}\tSYN\ttranscript\t{lo}\t{hi}\t.\t{st}\t.\t{base} level 2; {tags}")
        for k, (gs, ge) in enumerate(t["blocks"]):
            lines.append(f"{c}\tSYN\texon\t{gs + 1}\t{ge + 1}\t.\t{st}\t.\t{base} exon_number {k + 1}; level 2; {tags}")
        u5, lc = t["u5"], t["lc"]
        for k, s, e, before in seg_of(t, u5, u5 + lc):
            phase = (3 - before % 3) % 3
            lines.append(f'{c}\tSYN\tCDS\t{s}\t{e}\t.\t{st}\t{phase}\t{base} exon_number {k + 1}; level 2; '
                         f'protein_id "ENSPS{t["gene_id"][5:]}"; {tags}')
        for k, s, e, _ in seg_of(t, u5, u5 + 3):
            lines.append(f"{c}\tSYN\tstart_codon\t{s}\t{e}\t.\t{st}\t0\t{base} exon_number {k + 1}; {tags}")
        for k, s, e, _ in seg_of(t, u5 + lc, u5 + lc + 3):
            lines.append(f"{c}\tSYN\tstop_codon\t{s}\t{e}\t.\t{st}\t0\t{base} exon_number {k + 1}; {tags}")
        for lo_, hi_ in ((0, u5), (u5 + lc + 3, t["mlen"])):
            for k, s, e, _ in seg_of(t, lo_, hi_):
                lines.append(f"{c}\tSYN\tUTR\t{s}\t{e}\t.\t{st}\t.\t{base} exon_number {k + 1}; {tags}")
    with gzip.open(os.path.join(D, "gencode.v32.annotation.gtf.gz"), "wt") as f:
        f.write("\n".join(lines) + "\n")

    # -------------------------------------------------------------- Roulette-like rates
    rates, filt, mr_missing = {}, {}, {}
    for c in CHROMS:
        g = genome[c]
        prev = np.r_[["N"], g[:-1]]
        nxt = np.r_[g[1:], ["N"]]
        region = np.repeat(rng.lognormal(0, 0.25, size=L // 5000 + 1), 5000)[:L]
        R = np.full((L, 4), np.nan)
        for ai, ab in enumerate("ACGT"):
            ti = np.array([(r, ab) in TRANSITIONS for r in "ACGT"])[np.searchsorted(np.array(list("ACGT")), g)]
            cpg = ((g == "C") & (nxt == "G")) | ((g == "G") & (prev == "C"))
            base_rate = np.where(ti & cpg, 12.0, np.where(ti, 1.6, 0.5))
            r = 0.05 * base_rate * region * rng.lognormal(0, 0.4, size=L)
            r = sig4(r)
            r[g == ab] = np.nan
            R[:, ai] = r
        rates[c] = R
        fl = np.full(L, "high", dtype=object)
        fl[rng.random(L) < 0.02] = "SFS_bump"
        fl[rng.random(L) < 0.01] = "TFBS"
        for s in rng.integers(0, L - 2000, size=6):
            fl[s:s + 2000] = "low"
        filt[c] = fl
        mr_missing[c] = rng.random((L, 4)) < 0.0005

    rgap = {}
    for t in txs:
        if t["kind"] == "RGAP":
            gpos = int(t["m2g"][t["u5"] + 3 * 150 + 1])
            rgap[t["chrom"]] = gpos

    vcf_head = ["##fileformat=VCFv4.3"] + [f"##contig=<ID={c},length={L},assembly=GRCh38>" for c in CHROMS] + [
        '##FILTER=<ID=low,Description="Low quality regions">', '##FILTER=<ID=SFS_bump,Description="SFS bump">',
        '##FILTER=<ID=TFBS,Description="TF binding site">', '##FILTER=<ID=high,Description="Nothing suspicious">',
        '##INFO=<ID=PN,Number=1,Type=String,Description="Pentanucleotide context">',
        '##INFO=<ID=MR,Number=1,Type=Float,Description="Roulette mutation rate estimate">',
        '##INFO=<ID=AR,Number=1,Type=Float,Description="Adjusted Roulette mutation rate estimate">',
        '##INFO=<ID=MG,Number=1,Type=Float,Description="gnomAD mutation rate estimate">',
        '##INFO=<ID=MC,Number=1,Type=Float,Description="Carlson mutation rate estimate">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    with open(os.path.join(out, "roulette_vcf", "header.txt"), "w") as f:
        f.write("\n".join(vcf_head) + "\n")
    bases = np.array(list("ACGT"))
    for c in CHROMS:
        g = genome[c]
        gp = np.r_[["N", "N"], g, ["N", "N"]]
        pn = np.char.add(np.char.add(np.char.add(np.char.add(gp[:-4], gp[1:-3]), gp[2:-2]), gp[3:-1]), gp[4:])
        AR = sig4(rates[c] * rng.lognormal(0, 0.1, size=(L, 4)))
        MG = sig4(rates[c] * 1e-8 * rng.lognormal(0, 0.3, size=(L, 4)))
        MC = sig4(rates[c] * 0.02 * rng.lognormal(0, 0.3, size=(L, 4)))
        MR = np.where(mr_missing[c], np.nan, rates[c])
        pos = np.repeat(np.arange(L), 4)
        alt = np.tile(np.arange(4), L)
        keep = bases[alt] != g[pos]
        if c in rgap:
            keep &= pos != rgap[c]
        pos, alt = pos[keep], alt[keep]
        def txt(v):  # shortest round-trip text, '.' for missing (as bcftools prints)
            t = v.astype(str)
            t[~np.isfinite(v)] = "."
            return t
        df = pd.DataFrame({"POS": (pos + 1).astype(str), "REF": g[pos], "ALT": bases[alt], "FILTER": filt[c][pos],
                           "PN": pn[pos], "MR": txt(MR[pos, alt]), "AR": txt(AR[pos, alt]), "MG": txt(MG[pos, alt]),
                           "MC": txt(MC[pos, alt])})
        df.to_csv(os.path.join(out, "roulette_full", f"chr{c}.tsv.gz"), sep="\t", header=False, index=False,
                  compression={"method": "gzip", "compresslevel": 1})

    # -------------------------------------------------------------- reference enumeration + AlphaMissense
    am_rows, snv_rows = [], []
    for ti_, t in enumerate(txs):
        g = genome[t["chrom"]]
        cds_g = t["m2g"][t["u5"]:t["u5"] + t["lc"]]
        seq = "".join(g[x] if t["strand"] == "+" else CMP[g[x]] for x in cds_g)
        imp = rng.beta(0.5, 0.5, size=t["n_codons"])
        for ci in range(t["n_codons"]):
            codon = seq[3 * ci:3 * ci + 3]
            aa = CODE[codon]
            for j in range(3):
                gp_ = int(cds_g[3 * ci + j])
                for ab in "ACGT":
                    if ab == codon[j]:
                        continue
                    aa2 = CODE[codon[:j] + ab + codon[j + 1:]]
                    cls = "syn" if aa2 == aa else ("stop" if aa2 == "*" else ("mis" if aa != "*" else "stoplost"))
                    galt = ab if t["strand"] == "+" else CMP[ab]
                    snv_rows.append((ti_, t["chrom"], gp_ + 1, g[gp_], galt, ci + 1, cls, ci))
                    if cls == "mis":
                        lg = np.log(imp[ci] / (1 - imp[ci] + 1e-12) + 1e-12) + rng.normal(0, 0.8)
                        am = float(np.round(1 / (1 + np.exp(-lg)), 4))
                        am_rows.append((ti_, t["chrom"], gp_ + 1, g[gp_], galt, ci + 1, aa, aa2, am))
    snv = pd.DataFrame(snv_rows, columns=["ti", "chrom", "pos", "ref", "alt", "prot_pos", "cls", "codon"])
    amt = pd.DataFrame(am_rows, columns=["ti", "chrom", "pos", "ref", "alt", "prot_pos", "aa", "aa2", "am"])
    tdf = pd.DataFrame([{k: t[k] for k in ("chrom", "strand", "kind", "gene_name", "gene_id", "transcript_id", "n_codons")}
                        for t in txs])
    amt = amt.merge(tdf[["transcript_id"]].reset_index().rename(columns={"index": "ti"}), on="ti")
    amgap_ti = tdf.index[tdf["kind"] == "AMGAP"].tolist()
    drop = amt["ti"].isin(amgap_ti) & (amt["prot_pos"] > 1) & (rng.random(len(amt)) < 0.005)
    amt_out = amt[~drop].copy()
    shift = amt_out["ti"].isin(tdf.index[tdf["kind"] == "SHIFT"])
    amt_out["pp"] = np.where(shift, amt_out["prot_pos"] + 1, amt_out["prot_pos"])
    amt_out["protein_variant"] = amt_out["aa"] + amt_out["pp"].astype(str) + amt_out["aa2"]
    amt_out["am_class"] = np.where(amt_out["am"] > 0.564, "likely_pathogenic",
                                   np.where(amt_out["am"] < 0.34, "likely_benign", "ambiguous"))
    amt_out["CHROM"] = "chr" + amt_out["chrom"]
    amt_out = amt_out.sort_values(["CHROM", "pos"])
    with gzip.open(os.path.join(D, "AlphaMissense_hg38.tsv.gz"), "wt") as f:
        f.write("# Copyright 2023 DeepMind Technologies Limited (synthetic stand-in for the self-test)\n#\n"
                "# Licensed under CC BY 4.0\n")
        f.write("#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n")
        pd.DataFrame({"CHROM": amt_out["CHROM"], "POS": amt_out["pos"], "REF": amt_out["ref"], "ALT": amt_out["alt"],
                      "genome": "hg38", "uniprot_id": "P" + amt_out["ti"].astype(str).str.zfill(5),
                      "transcript_id": amt_out["transcript_id"], "protein_variant": amt_out["protein_variant"],
                      "am_pathogenicity": amt_out["am"].map(lambda x: f"{x:.4f}"), "am_class": amt_out["am_class"]}).to_csv(
            f, sep="\t", header=False, index=False)

    # -------------------------------------------------------------- the reference universe (what gate.py must find)
    good_kinds = {"normal", "AMGAP", "RGAP"}
    snv = snv.merge(tdf[["kind", "gene_name"]].reset_index().rename(columns={"index": "ti"}), on="ti")
    u = snv[snv["kind"].isin(good_kinds) & snv["cls"].isin(["syn", "mis"]) & (snv["prot_pos"] > 1)].copy()
    dropped = amt[drop][["ti", "pos", "alt"]].assign(gone=True)
    u = u.merge(dropped, on=["ti", "pos", "alt"], how="left")
    u = u[u["gone"].isna()].drop(columns="gone")
    for c, gp_ in rgap.items():
        bad = snv[(snv["chrom"] == c) & (snv["pos"] == gp_ + 1)][["ti", "codon"]].drop_duplicates()
        u = u.merge(bad.assign(badc=True), on=["ti", "codon"], how="left")
        u = u[u["badc"].isna()].drop(columns="badc")
    u = u.merge(amt[["ti", "pos", "alt", "am"]], on=["ti", "pos", "alt"], how="left")
    ai = {"A": 0, "C": 1, "G": 2, "T": 3}
    u["MR"] = [np.nan if mr_missing[c][p - 1, ai[al]] else rates[c][p - 1, ai[al]]
               for c, p, al in zip(u["chrom"], u["pos"], u["alt"])]
    u["filt_ok"] = [filt[c][p - 1] != "low" for c, p in zip(u["chrom"], u["pos"])]
    ctx = u[u["cls"] == "mis"].groupby(["ti", "codon"])["am"].mean().rename("ctx").reset_index()
    u = u.merge(ctx, on=["ti", "codon"], how="left")
    uu = u[u["filt_ok"] & u["MR"].notna()].copy()
    uu["bin"] = np.where(uu["cls"] == "mis", coarse(uu["am"].fillna(0)), -1)

    # -------------------------------------------------------------- simulate de novos
    sim = np.random.default_rng(a.sim_seed)
    dd_ti = sorted(sim.choice(tdf.index[tdf["kind"] == "normal"].to_numpy(), size=8, replace=False).tolist())
    wgene = {t: float(sim.lognormal(0, 0.6)) for t in dd_ti}
    # every SNV in every CDS can mutate (also excluded sites), plus some non-coding background
    allsnv = snv.copy()
    allsnv["rate"] = [rates[c][p - 1, ai[al]] for c, p, al in zip(allsnv["chrom"], allsnv["pos"], allsnv["alt"])]

    def cpg_ti(df):
        out_ = np.zeros(len(df), dtype=bool)
        for c in CHROMS:
            m = (df["chrom"] == c).to_numpy()
            p0 = df["pos"].to_numpy()[m] - 1
            g = genome[c]
            prv = np.where(p0 > 0, g[np.maximum(p0 - 1, 0)], "N")
            nxt = np.where(p0 < L - 1, g[np.minimum(p0 + 1, L - 1)], "N")
            ref, alt = df["ref"].to_numpy()[m], df["alt"].to_numpy()[m]
            ti = np.array([(r, x) in TRANSITIONS for r, x in zip(ref, alt)], dtype=bool)
            out_[m] = ti & (((ref == "C") & (nxt == "G")) | ((ref == "G") & (prv == "C")))
        return out_

    allsnv["rate"] = allsnv["rate"] * np.where(cpg_ti(allsnv), a.cpg_factor, 1.0)
    uu["rate_true"] = uu["MR"] * np.where(cpg_ti(uu), a.cpg_factor, 1.0)
    m_mis_total = uu.loc[uu["cls"] == "mis", "rate_true"].sum()
    k0 = 20000.0 / (m_mis_total * sum(n * s for n, s in CENTRES.values()))
    pat = uu[(uu["cls"] == "mis") & uu["ti"].isin(dd_ti)].copy()
    pat["lam"] = pat["rate_true"] * pat["ti"].map(wgene) * pat["am"] ** 3
    A0 = 6000.0 / (pat["lam"].sum() * sum(n * s for n, s in CENTRES.values()))
    ctx_all = snv.merge(amt[["ti", "pos", "alt", "am"]], on=["ti", "pos", "alt"], how="left").groupby(["ti", "codon"])["am"].mean()
    cols = ["chrom", "pos", "ref", "alt", "ti", "codon"]
    parts = []
    for cen, (n_prob, sens) in CENTRES.items():
        cnt = sim.poisson(k0 * n_prob * sens * allsnv["rate"].to_numpy())
        parts.append(allsnv.iloc[np.repeat(np.arange(len(allsnv)), cnt)][cols].assign(centre=cen, origin="background"))
        cntp = sim.poisson(A0 * n_prob * sens * pat["lam"].to_numpy())
        parts.append(pat.iloc[np.repeat(np.arange(len(pat)), cntp)][cols].assign(centre=cen, origin="pathogenic"))
    dn = pd.concat(parts, ignore_index=True)[["centre"] + cols + ["origin"]]
    dn["pos"] = dn["pos"].astype(int)
    if a.biased:
        cm = ctx_all.reindex(list(zip(dn["ti"], dn["codon"]))).to_numpy()
        lose = (np.nan_to_num(cm) > 0.75) & (sim.random(len(dn)) < 0.4)
        dn = dn[~lose]
    # junk rows that gate.py must drop: non-coding SNVs, indels, chrX, a deleted-in-hg38 region
    junk = []
    for c in CHROMS:
        for p in sim.integers(2000, L - 2000, size=40):
            alt = [b for b in "ACGT" if b != genome[c][p]][int(sim.integers(0, 3))]
            junk.append(("GeneDx", c, int(p) + 1, genome[c][p], alt, -1, -1, "noncoding"))
    dn = pd.concat([dn, pd.DataFrame(junk, columns=dn.columns)], ignore_index=True)
    dn = dn.sample(frac=1.0, random_state=5).reset_index(drop=True)
    pid = []
    for cen, (n_prob, _) in CENTRES.items():
        m = (dn["centre"] == cen).to_numpy()
        pid_c = sim.integers(0, n_prob, size=m.sum())
        pid.append(pd.Series(pid_c, index=dn.index[m]))
    dn["pnum"] = pd.concat(pid).sort_index()
    study = np.where(dn["centre"] == "GeneDx", np.where(dn["pnum"] < 6000, "GDX", "newGDX"), dn["centre"])
    dn["study"] = study
    dn["id"] = dn["centre"].str[:3].str.upper() + "_" + dn["pnum"].astype(str).str.zfill(5)

    # synthetic GRCh37 coordinates: chr22 offset chain with a deletion, chr21 reverse-complement chain
    def to_hg19(c, p38, ref, alt):
        q0 = p38 - 1
        if c == "22":
            if 1000 <= q0 < 151000:
                return q0 - 1000 + 1, ref, alt
            if 151000 <= q0 < L - 100:
                return q0 - 900 + 1, ref, alt
            return None
        return L - 1 - q0 + 1, CMP[ref], CMP[alt]

    rows = []
    truth = []
    for i, r in dn.iterrows():
        h = to_hg19(r["chrom"], r["pos"], r["ref"], r["alt"])
        if h is None:
            continue
        rows.append((r["chrom"], h[0], h[1], h[2], r["id"], r["study"]))
        truth.append((len(rows) - 1, r["chrom"], r["pos"], r["ref"], r["alt"], r["ti"], r["codon"], r["origin"]))
    extra = [("22", 150050, "A", "G", "GDX_00001", "GDX"),              # hg19 bases deleted in hg38: unmapped
             ("22", 2000, "AT", "A", "GDX_00002", "GDX"),               # indel
             ("X", 1000000, "C", "T", "DDD_00003", "DDD")]              # chrX
    for e in extra:
        rows.append(e)
        truth.append((len(rows) - 1, e[0], -1, e[2], e[3], -1, -1, "junk"))
    dnv = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "id", "study"])
    tr = pd.DataFrame(truth, columns=["dnv_idx", "chrom", "pos38", "ref38", "alt38", "ti", "codon", "origin"])
    ukey = uu[["chrom", "pos", "alt", "cls", "bin"]].rename(columns={"pos": "pos38", "alt": "alt38"})
    ukey = ukey.drop_duplicates(["chrom", "pos38", "alt38"], keep=False)
    key = tr.merge(ukey, on=["chrom", "pos38", "alt38"], how="left")
    tr["universe_cls"] = key["cls"].fillna("none").to_numpy()
    tr["bin"] = key["bin"].fillna(-1).astype(int).to_numpy()
    snvcls = snv[["ti", "pos", "alt", "cls"]].rename(columns={"pos": "pos38", "alt": "alt38", "cls": "ref_cls"})
    tr = tr.merge(snvcls, on=["ti", "pos38", "alt38"], how="left")
    cqmap = {"mis": "missense_variant", "syn": "synonymous_variant", "stop": "stop_gained", "stoplost": "stop_lost"}
    dnv["cq"] = tr["ref_cls"].map(cqmap).fillna("intron_variant").to_numpy()
    dnv["consequence"] = dnv["cq"]
    sym = tdf["gene_name"].to_dict()
    dnv["symbol"] = tr["ti"].map(sym).fillna("").to_numpy()
    for col, val in (("altprop_child", 0.5), ("hgnc_id", 0), ("maf", 0.0), ("prob", 0.99), ("raw", 1.0), ("score", 20.0),
                     ("constrained", "FALSE"), ("shethigh", "FALSE")):
        dnv[col] = val
    dnv = dnv[["alt", "altprop_child", "chrom", "consequence", "cq", "hgnc_id", "id", "maf", "pos", "prob", "raw",
               "ref", "score", "study", "symbol", "constrained", "shethigh"]]
    dnv.to_csv(os.path.join(D, "DeNovoWEST/input/DDD_RUMC_GDX_denovos_cadd_shet_wweights_minrep_indelfix_2025_09_18.txt.gz"),
               sep="\t", index=False)
    tr.to_csv(os.path.join(out, "truth", "dnv_truth.tsv"), sep="\t", index=False)

    # chain file (hg19 = target/source of liftOver, hg38 = query)
    with gzip.open(os.path.join(D, "hg19ToHg38.over.chain.gz"), "wt") as f:
        f.write(f"chain 1000 chr22 {L - 1000} + 0 {L - 1000} chr22 {L} + 1000 {L - 100} 1\n150000\t100\t0\n{L - 1000 - 150100}\n\n")
        f.write(f"chain 1000 chr21 {L} + 0 {L} chr21 {L} - 0 {L} 2\n{L}\n\n")

    # gene table and DDG2P
    gt = tdf[["gene_name", "chrom"]].rename(columns={"gene_name": "symbol", "chrom": "chr"})
    gt["s_het"] = [0.3 if i in dd_ti else float(sim.uniform(0.005, 0.25)) for i in tdf.index]
    gt.loc[gt.index[:2], "s_het"] = np.nan
    gt.to_csv(os.path.join(D, "DeNovoWEST/input/extended_denovoWEST_results.tab"), sep="\t", index=False)
    ddr = []
    for i in dd_ti:
        nm = tdf.loc[i, "gene_name"]
        if i == dd_ti[0]:     # listed under a newer symbol; matched through 'prev symbols'
            ddr.append([nm + "X", "", "SYN DISORDER", "", "confirmed", "monoallelic", "all missense/in frame", "", "", "", "DD", nm, "1", ""])
        else:
            ddr.append([nm, "", "SYN DISORDER", "", "confirmed", "monoallelic", "all missense/in frame", "", "", "", "DD", "", "1", ""])
    for i in tdf.index[tdf["kind"] == "normal"][:3]:
        if i not in dd_ti:
            ddr.append([tdf.loc[i, "gene_name"], "", "SYN RECESSIVE", "", "confirmed", "biallelic", "loss of function", "", "", "", "DD", "", "1", ""])
    pd.DataFrame(ddr, columns=["gene symbol", "gene mim", "disease name", "disease mim", "DDD category",
                               "allelic requirement", "mutation consequence", "phenotypes", "organ specificity list",
                               "pmids", "panel", "prev symbols", "hgnc id", "gene disease pair entry date"]).to_csv(
        os.path.join(D, "DeNovoWEST/input/DDG2P_22_9_2019.csv"), index=False)

    # -------------------------------------------------------------- true likelihood ratios
    nb = len(EDGES) - 1
    mis_u = uu[uu["cls"] == "mis"]
    M_mis = np.bincount(mis_u["bin"], weights=mis_u["MR"], minlength=nb)
    fB = np.bincount(mis_u["bin"], weights=mis_u["rate_true"], minlength=nb)
    lamP = np.bincount(coarse(pat["am"]), weights=pat["lam"], minlength=nb)
    lr_marg = (lamP / lamP.sum()) / (fB / fB.sum())
    Xg = pat.groupby("ti")["lam"].sum()
    den = np.zeros(nb)
    for t_, xg in Xg.items():
        mg = mis_u[mis_u["ti"] == t_]
        den += xg * np.bincount(mg["bin"], weights=mg["rate_true"], minlength=nb) / mg["rate_true"].sum()
    lr_within = lamP / den
    mis_dd = mis_u[mis_u["ti"].isin(dd_ti)]
    fBd = np.bincount(mis_dd["bin"], weights=mis_dd["rate_true"], minlength=nb)
    lr_dd = (lamP / lamP.sum()) / (fBd / fBd.sum())
    syn_u = uu[uu["cls"] == "syn"]
    pd.DataFrame({"bin_lo": [e / 100 for e in EDGES[:-1]], "bin_hi": [e / 100 for e in EDGES[1:]],
                  "M_mis": M_mis, "LR_marginal": lr_marg, "LR_within": lr_within, "LR_ddg2p_marginal": lr_dd}).to_csv(
        os.path.join(out, "truth", "lr_truth.tsv"), sep="\t", index=False)
    meta = {"M_syn": float(syn_u["MR"].sum()), "M_mis": float(mis_u["MR"].sum()), "biased": bool(a.biased),
            "cpg_factor": a.cpg_factor,
            "dd_genes": [tdf.loc[i, "gene_name"] for i in dd_ti],
            "expected_status": {r["transcript_id"]: r["kind"] for _, r in tdf.iterrows()},
            "n_dnv_rows": int(len(dnv))}
    with open(os.path.join(out, "truth", "meta.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print(f"synthetic data written to {out}: {len(txs)} transcripts, {len(dnv)} de novo rows "
          f"({int((tr['origin'] == 'pathogenic').sum())} pathogenic), biased={a.biased}, cpg_factor={a.cpg_factor}")


if __name__ == "__main__":
    main()
