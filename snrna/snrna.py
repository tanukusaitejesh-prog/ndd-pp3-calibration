#!/usr/bin/env python3
"""snRNA-VEP v2 pipeline. Every analysis is defined in PREREG.md; section numbers below refer to it.

Commands (run in order; each caches its downloads under data/ and is safe to re-run):
  fetch        pinned small inputs: SGE table, GENCODE v49, GIAB stratification
  genes        111-gene universe, gnomAD coverage and paralog QC -> analysis set        (s2)
  sites        every possible SNV in the analysis set + RNU4-2: Roulette, CADD, gnomAD   (s3, s7.1)
  sge          SGE table QC, GRCh38 REF check, gnomAD/CADD joins                         (s3, s4)
  vienna       recompute U4/U6 ddG_bind for SGE variants (QC against the paper)          (s6.2)
  evo2-prep    8,192-nt ref/alt windows for evo2_score.py (run that on a GPU)            (s6.1)
  g1           gate G1 metrics, bootstrap intervals, reproduction check, decision        (s6)
  seal/unseal  sealed blind-gene labels                                                  (s5)

Remote data are read by tabix over HTTPS (pysam), so nothing genome-wide is downloaded.
"""
import argparse
import datetime
import getpass
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("SNRNA_DATA", os.path.join(HERE, "data"))
WORK = os.environ.get("SNRNA_WORK", os.path.join(HERE, "work"))
RESULTS = os.environ.get("SNRNA_RESULTS", os.path.join(HERE, "results"))
SEED = 20261013

URL = {
    "sge": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10334-9/MediaObjects/41586_2026_10334_MOESM2_ESM.xlsx",
    "gencode": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz",
    "giab": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/Union/GRCh38_alllowmapandsegdupregions.bed.gz",
    "gnomad_vcf": "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz",
    "gnomad_api": "https://gnomad.broadinstitute.org/api",
    "roulette": "https://genetics.bwh.harvard.edu/downloads/Vova/Roulette/{n}_rate_v5.2_TFBS_correction_all.vcf.bgz",
    "cadd": "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz",
    "ucsc_seq": "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}",
}
FILES = {"sge": "sge_ST1.xlsx", "gencode": "gencode.v49.annotation.gtf.gz", "giab": "GRCh38_alllowmapandsegdupregions.bed.gz"}

# s2: GENCODE v49 spliceosomal snRNA families (pseudogenes, i.e. names ending in P, are excluded separately)
SPLICEOSOMAL = re.compile(r"^(RNU1-\d+|RNVU1-\w+|RNU2-\d+|RNU4-\d+|RNU4ATAC|RNU5[A-F]-\d+|RNU6-\d+|RNU6ATAC|RNU11|RNU12"
                          r"|U1|U2|U4|U6|U4atac|U6atac|U11|U12)$")
SEALED_GENES = ("RNU2-2", "RNU5B-1")
U4_FAMILY_NAMES = ("RNU4-1", "RNU4-2", "RNU4ATAC", "U4")

# s4: RNU4-2 regions (Crick code assign_region, integer positions on NR_003137.3) and the critical region
RNU4_2_REGIONS = [("Stem II", [(3, 16)]), ("k-turn", [(27, 35), (41, 46)]), ("Stem I", [(56, 62)]),
                  ("T-loop", [(63, 70)]), ("Stem III", [(75, 79)]), ("3' stem loop", [(85, 117)]),
                  ("Sm site", [(118, 126)]), ("terminal stem loop", [(127, 144)])]
CRITICAL = list(range(62, 71)) + list(range(75, 79))
SGE_DEPLETED, SGE_STRONG = -0.302, -0.9
BLIND_HOTSPOTS = {"RNU2-2": (1, 67), "RNU5B-1": (37, 44)}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def log(*a):
    print(*a, flush=True)


def ensure_dirs():
    for d in (DATA, WORK, RESULTS, os.path.join(DATA, "cache")):
        os.makedirs(d, exist_ok=True)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def http_get(url, dest, tries=6):
    if os.path.exists(dest + ".ok"):
        return dest
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "snrna-vep/2 (research pipeline)"})
            with urllib.request.urlopen(req, timeout=300) as r, open(dest + ".tmp", "wb") as f:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b)
            os.replace(dest + ".tmp", dest)
            open(dest + ".ok", "w").close()
            return dest
        except Exception as e:  # network hiccups: back off and retry
            log(f"  retry {i + 1}/{tries} for {url}: {e}")
            time.sleep(2 ** (i + 1))
    raise RuntimeError(f"download failed: {url}")


def cached_json(key, fn):
    path = os.path.join(DATA, "cache", key + ".json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    out = fn()
    with open(path + ".tmp", "w") as f:
        json.dump(out, f)
    os.replace(path + ".tmp", path)
    return out


def cached_tsv(key, fn):
    path = os.path.join(DATA, "cache", key + ".tsv.gz")
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t", dtype={"chrom": str})
    df = fn()
    df.to_csv(path + ".tmp.gz", sep="\t", index=False)
    os.replace(path + ".tmp.gz", path)
    return df


def run_info(name, extra=None):
    info = {"command": " ".join(sys.argv), "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__}
    try:
        info["git_commit"] = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        info["git_dirty"] = bool(subprocess.run(["git", "-C", HERE, "status", "--porcelain", "--", "."], capture_output=True, text=True).stdout.strip())
    except FileNotFoundError:
        pass
    mf = os.path.join(DATA, "manifest.json")
    if os.path.exists(mf):
        with open(mf) as f:
            info["inputs_md5"] = json.load(f)
    info.update(extra or {})
    with open(os.path.join(RESULTS, f"run_info_{name}.json"), "w") as f:
        json.dump(info, f, indent=1)


# ----------------------------------------------------------------------------------------------- remote sources

def ucsc_seq(chrom, start0, end):
    """GRCh38 bases [start0, end) (0-based, half-open), upper case, cached."""
    def get():
        url = URL["ucsc_seq"].format(chrom=chrom, start=start0, end=end)
        for i in range(6):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    return json.load(r)["dna"].upper()
            except Exception as e:
                log(f"  UCSC retry {i + 1}: {e}")
                time.sleep(2 ** (i + 1))
        raise RuntimeError(url)
    return cached_json(f"seq_{chrom}_{start0}_{end}", get)


def gnomad_api(query, key):
    def get():
        body = json.dumps({"query": query}).encode()
        for i in range(8):
            try:
                req = urllib.request.Request(URL["gnomad_api"], data=body, headers={"Content-Type": "application/json"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    out = json.load(r)
                if "errors" in out:
                    raise RuntimeError(out["errors"])
                time.sleep(1.0)  # be polite: the public API is rate limited
                return out["data"]
            except Exception as e:
                log(f"  gnomAD API retry {i + 1}: {e}")
                time.sleep(5 * (i + 1))
        raise RuntimeError(key)
    return cached_json(key, get)


def gnomad_coverage(chrom, start, end):
    q = ('{ region(chrom:"%s", start:%d, stop:%d, reference_genome:GRCh38){ coverage(dataset:gnomad_r4){ genome{ pos median } } } }'
         % (chrom.replace("chr", ""), start, end))
    d = gnomad_api(q, f"cov_{chrom}_{start}_{end}")
    cov = {c["pos"]: c["median"] for c in d["region"]["coverage"]["genome"]}
    return [cov.get(p, 0) for p in range(start, end + 1)]


ANCESTRIES = ("afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "remaining", "sas")


def gnomad_records(chrom, start, end):
    """All gnomAD v4.1 genome records overlapping [start, end] (1-based), one row per alt allele."""
    def get():
        import pysam
        rows = []
        vf = pysam.VariantFile(URL["gnomad_vcf"].format(chrom=chrom))
        for r in vf.fetch(chrom, start - 1, end):
            flt = ";".join(r.filter.keys()) or "PASS"
            for i, alt in enumerate(r.alts):
                info = r.info
                ab = info.get("ab_hist_alt_bin_freq")
                ab = [int(x) for x in ab[i].split("|")] if ab else [0] * 20
                row = {"chrom": chrom, "pos": r.pos, "ref": r.ref, "alt": alt, "filter": flt,
                       "AC": info["AC"][i], "AN": info["AN"], "AF": info["AF"][i] if info.get("AF") else 0.0,
                       "nhomalt": info["nhomalt"][i], "ab_het_lt025": sum(ab[:5]), "ab_het_all": sum(ab)}
                for a in ANCESTRIES:
                    row[f"AC_{a}"] = info.get(f"AC_{a}", (0,))[i] if info.get(f"AC_{a}") is not None else 0
                    row[f"AN_{a}"] = info.get(f"AN_{a}", 0)
                rows.append(row)
        return pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "filter", "AC", "AN", "AF", "nhomalt", "ab_het_lt025",
                                           "ab_het_all"] + [f"{p}_{a}" for a in ANCESTRIES for p in ("AC", "AN")])
    return cached_tsv(f"gnomad_{chrom}_{start}_{end}", get)


def roulette_records(chrom, start, end):
    def get():
        import pysam
        n = chrom.replace("chr", "")
        url = URL["roulette"].format(n=n)
        tb = pysam.TabixFile(url, index=url + ".csi")
        rows = []
        for line in tb.fetch(n, start - 1, end):
            c = line.split("\t")
            info = dict(kv.split("=") for kv in c[7].split(";"))
            rows.append({"chrom": chrom, "pos": int(c[1]), "ref": c[3], "alt": c[4], "roulette_filter": c[6],
                         "MR": float(info.get("MR", "nan")), "AR": float(info.get("AR", "nan")), "PN": info.get("PN", "")})
        return pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "roulette_filter", "MR", "AR", "PN"])
    return cached_tsv(f"roulette_{chrom}_{start}_{end}", get)


def cadd_records(chrom, start, end):
    def get():
        import pysam
        tb = pysam.TabixFile(URL["cadd"], index=URL["cadd"] + ".tbi")
        rows = [line.split("\t") for line in tb.fetch(chrom.replace("chr", ""), start - 1, end)]
        df = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "cadd_raw", "cadd_phred"])
        df["chrom"] = chrom
        return df.astype({"pos": int, "cadd_raw": float, "cadd_phred": float})
    return cached_tsv(f"cadd_{chrom}_{start}_{end}", get)


# ----------------------------------------------------------------------------------------------------- s2 genes

def cmd_fetch(args):
    ensure_dirs()
    man = {}
    for k, fn in FILES.items():
        dest = os.path.join(DATA, fn)
        log(f"{k}: {fn}")
        http_get(URL[k], dest)
        man[fn] = md5(dest)
    with open(os.path.join(DATA, "manifest.json"), "w") as f:
        json.dump(man, f, indent=1)
    log(json.dumps(man, indent=1))


def load_genes():
    """GENCODE v49 spliceosomal snRNA genes (s2.1), one row per gene with its transcript span and strand."""
    path = os.path.join(WORK, "genes_v49.tsv")
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t")
    rows, exons = [], {}
    with gzip.open(os.path.join(DATA, FILES["gencode"]), "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if 'gene_type "snRNA"' not in c[8]:
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', c[8]))
            if c[2] == "gene":
                rows.append({"gene_id": a["gene_id"], "gene": a["gene_name"], "chrom": c[0], "start": int(c[3]),
                             "end": int(c[4]), "strand": c[6]})
            elif c[2] == "exon":
                exons.setdefault(a["gene_id"], []).append((a["transcript_id"], int(c[3]), int(c[4])))
    g = pd.DataFrame(rows)
    keep = g.gene.str.match(SPLICEOSOMAL) & ~g.gene.str.endswith("P") & g.chrom.str.match(r"^chr([0-9]+|X|Y)$")
    g = g[keep].copy()
    n_exons = {gid: len(v) for gid, v in exons.items()}
    n_tx = {gid: len({t for t, _, _ in v}) for gid, v in exons.items()}
    g["n_transcripts"] = g.gene_id.map(n_tx)
    g["n_exons"] = g.gene_id.map(n_exons)
    g["length"] = g.end - g.start + 1
    # unnamed loci (U1, U2, ...) get a unique label from their gene_id
    g["label"] = np.where(g.gene.str.match(r"^U\d|^U\d+atac"), g.gene + "|" + g.gene_id.str.split(".").str[0], g.gene)
    g = g.sort_values(["chrom", "start"]).reset_index(drop=True)
    g.to_csv(path, sep="\t", index=False)
    return g


def tx_positions(row):
    """Genomic positions of transcript n.1..n.L for a single-exon gene."""
    pos = np.arange(row.start, row.end + 1)
    return pos if row.strand == "+" else pos[::-1]


def read_bed_overlap(bed_gz, regions):
    """For each (chrom, start, end) 1-based inclusive region, number of bases inside the BED intervals."""
    want = {}
    for i, (c, s, e) in enumerate(regions):
        want.setdefault(c, []).append((i, s, e))
    out = np.zeros(len(regions), dtype=int)
    with gzip.open(bed_gz, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            c, s0, e = line.split("\t")[:3]
            if c not in want:
                continue
            s0, e = int(s0), int(e)
            for i, rs, re_ in want[c]:
                lo, hi = max(rs, s0 + 1), min(re_, e)
                if hi >= lo:
                    out[i] += hi - lo + 1
    return out


def cmd_genes(args):
    ensure_dirs()
    g = load_genes()
    log(f"spliceosomal snRNA genes in GENCODE v49: {len(g)} (guidance: 111)")
    med, mn = [], []
    for r in g.itertuples():
        cov = gnomad_coverage(r.chrom, r.start, r.end)
        med.append(float(np.median(cov)))
        mn.append(float(np.min(cov)))
    g["cov_median"], g["cov_min"] = med, mn
    g["pass_coverage"] = (g.cov_median >= 20) & (g.cov_min >= 10)
    g["giab_lowmap_segdup_bp"] = read_bed_overlap(os.path.join(DATA, FILES["giab"]), list(zip(g.chrom, g.start, g.end)))
    g["giab_frac"] = g.giab_lowmap_segdup_bp / g.length
    lt, tot, npass = [], [], []
    for r in g.itertuples():
        v = gnomad_records(r.chrom, r.start, r.end)
        v = v[v["filter"] == "PASS"]
        lt.append(int(v.ab_het_lt025.sum()))
        tot.append(int(v.ab_het_all.sum()))
        npass.append(len(v))
    g["het_ab_lt025"], g["het_ab_all"], g["n_pass_records"] = lt, tot, npass
    g["low_vaf_frac"] = np.where(g.het_ab_all > 0, g.het_ab_lt025 / g.het_ab_all.clip(lower=1), np.nan)
    g["pass_giab"] = g.giab_frac <= 0.5
    g["pass_low_vaf"] = ~(g.low_vaf_frac > 0.5)
    g["analysis_set"] = g.pass_coverage & g.pass_giab & g.pass_low_vaf
    g["sealed"] = g.gene.isin(SEALED_GENES)
    g["u4_family"] = g.gene.isin(U4_FAMILY_NAMES)
    reasons = []
    for r in g.itertuples():
        why = []
        if not r.pass_coverage:
            why.append(f"coverage (median {r.cov_median:g}, min {r.cov_min:g})")
        if not r.pass_giab:
            why.append(f"GIAB lowmap/segdup {r.giab_frac:.0%}")
        if not r.pass_low_vaf:
            why.append(f"low-VAF hets {r.low_vaf_frac:.0%}")
        reasons.append("; ".join(why))
    g["excluded_because"] = reasons
    g.to_csv(os.path.join(RESULTS, "week0_gene_qc.tsv"), sep="\t", index=False)
    g[g.analysis_set].to_csv(os.path.join(WORK, "analysis_genes.tsv"), sep="\t", index=False)
    lines = [f"Spliceosomal snRNA genes (GENCODE v49): {len(g)}",
             f"Pass gnomAD genome coverage (median >= 20, min >= 10): {int(g.pass_coverage.sum())} (guidance: 87)",
             f"Fail GIAB low-mappability/segdup > 50%: {int((~g.pass_giab).sum())}",
             f"Fail mostly-low-VAF heterozygotes: {int((~g.pass_low_vaf).sum())}",
             f"Analysis set: {int(g.analysis_set.sum())}"]
    for name in ("RNU4-2", "RNU4-1", "RNU4ATAC", "RNU2-2", "RNU5B-1", "RNU6ATAC", "RNU12", "RNU11", "RNU2-1", "RNU5A-1"):
        r = g[g.gene == name]
        if len(r):
            r = r.iloc[0]
            lines.append(f"  {name:9s} cov median {r.cov_median:5.1f} min {r.cov_min:5.1f}  GIAB {r.giab_frac:4.0%}  "
                         f"low-VAF {r.low_vaf_frac:5.1%}  -> {'IN' if r.analysis_set else 'OUT: ' + r.excluded_because}")
    fams = g[g.analysis_set].gene.str.replace(r"-\w+$", "", regex=True).value_counts()
    lines.append("Analysis set by family: " + ", ".join(f"{k} {v}" for k, v in fams.items()))
    txt = "\n".join(lines)
    with open(os.path.join(RESULTS, "week0_gene_qc.txt"), "w") as f:
        f.write(txt + "\n")
    log(txt)
    run_info("genes")


# ----------------------------------------------------------------------------------------------------- s3 sites

def gene_sites(r, flank=0):
    """Every possible SNV at each transcript base of gene row r, with transcript position n."""
    seq = ucsc_seq(r.chrom, r.start - 1 - flank, r.end + flank)
    pos = np.arange(r.start - flank, r.end + flank + 1)
    rows = []
    for p, b in zip(pos, seq):
        n = (p - r.start + 1) if r.strand == "+" else (r.end - p + 1)
        for alt in "ACGT":
            if alt != b:
                rows.append((r.label, r.gene, r.chrom, int(p), b, alt, int(n), r.strand))
    return pd.DataFrame(rows, columns=["label", "gene", "chrom", "pos", "ref", "alt", "n", "strand"])


def cmd_sites(args):
    ensure_dirs()
    g = pd.read_csv(os.path.join(RESULTS, "week0_gene_qc.tsv"), sep="\t")
    use = g[g.analysis_set | (g.gene == "RNU4-2")]
    out = []
    for r in use.itertuples():
        s = gene_sites(r)
        rl = roulette_records(r.chrom, r.start, r.end)
        cd = cadd_records(r.chrom, r.start, r.end)
        gn = gnomad_records(r.chrom, r.start, r.end)
        gn = gn[(gn.ref.str.len() == 1) & (gn.alt.str.len() == 1)]
        s = s.merge(rl[["pos", "ref", "alt", "roulette_filter", "MR", "AR", "PN"]], on=["pos", "ref", "alt"], how="left")
        s = s.merge(cd[["pos", "ref", "alt", "cadd_phred"]], on=["pos", "ref", "alt"], how="left")
        s = s.merge(gn.drop(columns=["chrom"]), on=["pos", "ref", "alt"], how="left")
        s["in_analysis_set"] = bool(r.analysis_set)
        # REF check against every source that reports a REF at this base
        refs = {}
        for src in (rl, cd, gn):
            for p, b in zip(src.pos, src.ref.str[0]):
                refs.setdefault(p, set()).add(b)
        s["ref_conflict"] = [bool(refs.get(p)) and (b not in refs[p]) for p, b in zip(s.pos, s.ref)]
        out.append(s)
        log(f"{r.label}: {len(s)} SNVs, Roulette {s.MR.notna().mean():.0%}, CADD {s.cadd_phred.notna().mean():.0%}, "
            f"gnomAD PASS {((s['filter'] == 'PASS') & (s.AC > 0)).sum()}, REF conflicts {int(s.ref_conflict.sum())}")
    sites = pd.concat(out, ignore_index=True)
    sites.to_csv(os.path.join(WORK, "sites.tsv.gz"), sep="\t", index=False)
    n_conf = int(sites.ref_conflict.sum())
    summ = (sites.groupby("label").agg(n_snv=("pos", "size"), roulette=("MR", lambda x: x.notna().mean()),
                                       roulette_low=("roulette_filter", lambda x: (x == "low").mean()),
                                       cadd=("cadd_phred", lambda x: x.notna().mean()),
                                       gnomad_pass=("filter", lambda x: (x == "PASS").sum()))
            .reset_index())
    summ.to_csv(os.path.join(RESULTS, "week0_sites_summary.tsv"), sep="\t", index=False)
    txt = (f"sites: {len(sites)} SNVs in {sites.label.nunique()} genes; REF conflicts {n_conf} "
           f"({n_conf / len(sites):.3%}); Roulette coverage {sites.MR.notna().mean():.1%} "
           f"(FILTER low {(sites.roulette_filter == 'low').mean():.1%}); CADD coverage {sites.cadd_phred.notna().mean():.1%}")
    log(txt)
    with open(os.path.join(RESULTS, "week0_sites.txt"), "w") as f:
        f.write(txt + "\n")
    if n_conf / len(sites) > 0.01:
        sys.exit("REF conflict rate above 1%: pipeline QC failure (PREREG s3)")
    run_info("sites")


# ------------------------------------------------------------------------------------------------------- s4 SGE

def normalise(chrom, pos, ref, alt, seqfn):
    """Left-align and trim a VCF-style allele pair against the reference (small indels only)."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    if len(ref) != len(alt):
        for _ in range(100):  # shift left while the last bases match
            if ref[-1] != alt[-1]:
                break
            prev = seqfn(chrom, pos - 1)
            ref, alt, pos = prev + ref[:-1], prev + alt[:-1], pos - 1
    return pos, ref, alt


def rnu4_2_region(n):
    if pd.isna(n):
        return "outside transcript"
    for name, spans in RNU4_2_REGIONS:
        if any(a <= n <= b for a, b in spans):
            return name
    return "Other"


def load_sge():
    d = pd.read_excel(os.path.join(DATA, FILES["sge"]), sheet_name="ST1_20251212_RNU4_2_merged_gf", header=1)
    d = d[d.ID.astype(str).str.startswith("chr12-")].copy()
    p = d.ID.str.split("-", expand=True)
    d["chrom"], d["pos"], d["ref"], d["alt"] = p[0], p[1].astype(int), p[2], p[3]
    d["n"] = pd.to_numeric(d.nucleotide, errors="coerce")
    d["in_transcript"] = d.n.between(1, 145)
    d["is_snv"] = d.Type == "SNV"
    d["depleted"] = d.function_score_mean < SGE_DEPLETED
    d["strong"] = d.function_score_mean < SGE_STRONG
    d["region"] = d.n.map(rnu4_2_region)
    d["critical"] = d.n.isin(CRITICAL)
    return d.reset_index(drop=True)


def cmd_sge(args):
    ensure_dirs()
    d = load_sge()
    lo, hi = d.pos.min(), d.pos.max()
    ref = ucsc_seq("chr12", lo - 101, hi + 100)  # [lo-101, hi+100) 0-based covers lo-100..hi+100

    def base(chrom, p):
        return ref[p - (lo - 100)]

    d["ref_ok"] = [ref[p - (lo - 100): p - (lo - 100) + len(r)] == r for p, r in zip(d.pos, d.ref)]
    norm = [normalise("chr12", p, r, a, base) for p, r, a in zip(d.pos, d.ref, d.alt)]
    d["npos"], d["nref"], d["nalt"] = zip(*norm)
    gn = gnomad_records("chr12", lo - 100, hi + 100)
    gn_norm = {normalise("chr12", p, r, a, base): (f, ac, af) for p, r, a, f, ac, af in
               zip(gn.pos, gn.ref, gn.alt, gn["filter"], gn.AC, gn.AF)}
    hit = [gn_norm.get(k) for k in zip(d.npos, d.nref, d.nalt)]
    d["gnomad_filter"] = [h[0] if h else "" for h in hit]
    d["gnomad_AC_raw"] = [h[1] if h else 0 for h in hit]
    d["gnomad_AC"] = [h[1] if (h and h[0] == "PASS") else 0 for h in hit]
    d["gnomad_AF"] = [h[2] if (h and h[0] == "PASS") else 0.0 for h in hit]
    cd = cadd_records("chr12", lo, hi)
    cadd = {(p, r, a): s for p, r, a, s in zip(cd.pos, cd.ref, cd.alt, cd.cadd_phred)}
    d["cadd_v17_file"] = [cadd.get((p, r, a), np.nan) if s else np.nan for p, r, a, s in zip(d.pos, d.ref, d.alt, d.is_snv)]
    d.to_csv(os.path.join(WORK, "sge.tsv"), sep="\t", index=False)

    t = d[d.in_transcript]
    rep = t[["function_score_r1", "function_score_r2", "function_score_r3"]].corr().values[np.triu_indices(3, 1)]
    snv_cadd = d[d.is_snv & d.CADD_score.notna() & d.cadd_v17_file.notna()]
    cadd_diff = (snv_cadd.CADD_score - snv_cadd.cadd_v17_file).abs()
    checks = [
        ("variants in SGE ST1", len(d), 539),
        ("depleted (< -0.302)", int(d.depleted.sum()), 151),
        ("in-transcript (n.1-145)", int(d.in_transcript.sum()), 521),
        ("ReNU SNVs", int(((d.category == "ReNU_syndrome") & d.is_snv).sum()), 12),
        ("population SNVs (UKBB/AllofUs)", int(((d.category == "UKBB/AllofUs") & d.is_snv).sum()), 346),
        ("library start chr12", int(lo), 120291753),
        ("library end chr12", int(hi), 120291903),
        ("REF mismatches vs GRCh38", int((~d.ref_ok).sum()), 0),
    ]
    lines = ["Week-0 SGE QC (PREREG s3 known-fact checks)"]
    ok = True
    for name, got, want in checks:
        flag = "OK" if got == want else "MISMATCH"
        ok &= got == want
        lines.append(f"  {name:34s} {got:>10}  expected {want:>10}  {flag}")
    lines += [f"  strong (< -0.9): {int(d.strong.sum())}; SNVs {int(d.is_snv.sum())}; in-transcript SNVs {int((d.in_transcript & d.is_snv).sum())}, indels {int((d.in_transcript & ~d.is_snv).sum())}",
              f"  replicate Pearson r (in-transcript): {', '.join(f'{x:.3f}' for x in rep)} (paper: 0.83-0.86)",
              f"  CADD column vs CADD v1.7 file, SNVs: n={len(snv_cadd)}, max |diff| {cadd_diff.max():.3f}, "
              f"share within 0.01: {(cadd_diff <= 0.01).mean():.1%}",
              f"  in-transcript variants found in gnomAD v4.1 genomes: PASS {int((t.gnomad_AC > 0).sum())}, "
              f"non-PASS {int(((t.gnomad_filter != '') & (t.gnomad_filter != 'PASS')).sum())}",
              f"  in-transcript non-ReNU SNVs with gnomAD PASS AC >= 1 (primary G1 control set): "
              f"{int((t.is_snv & (t.category != 'ReNU_syndrome') & (t.gnomad_AC > 0)).sum())}",
              f"  ReNU SNVs present in gnomAD PASS: {int((t.is_snv & (t.category == 'ReNU_syndrome') & (t.gnomad_AC > 0)).sum())}"]
    txt = "\n".join(lines)
    with open(os.path.join(RESULTS, "week0_sge_qc.txt"), "w") as f:
        f.write(txt + "\n")
    log(txt)
    run_info("sge")
    if not ok:
        sys.exit("SGE known-fact check failed (PREREG s3): fix the pipeline before any endpoint")


# --------------------------------------------------------------------------------------------------- s6.2 Vienna

RNU4_2_RNA = "AGCUUUGCGCAGUGGCAGUAUCGUAGCCAAUGAGGUUUAUCCGAGGCGCGAUUAUUGCUAAUUGAAAACUUUUCCCAAUACCCCGCCAUGACGACUUGAAAUAUAGUCGGCAUUGGCAAUUUUUGACAGUCUCUACGGAGACUGA"
RNU6_RNA = "GUGCUCGCUUCGGCAGCACAUAUACUAAAAUUGGAACGAUACAGAGAAGAUUAGCAUGGCCCCUGCGCAAGGAUGACACGCAAAUUCGUGAAGCGUUCCAUAUUUU"


def apply_hgvs_n(seq, hgvs):
    """Apply an n. HGVS edit (substitution, del, dup, ins) to an RNA sequence; DNA letters are converted."""
    h = hgvs.split(":")[-1]
    body = h[2:].replace("T", "U")
    m = re.fullmatch(r"(\d+)([ACGU])>([ACGU])", body)
    if m:
        p = int(m.group(1))
        assert seq[p - 1] == m.group(2), (hgvs, seq[p - 1])
        return seq[:p - 1] + m.group(3) + seq[p:]
    m = re.fullmatch(r"(\d+)(?:_(\d+))?del[ACGU]*", body)
    if m:
        a = int(m.group(1)); b = int(m.group(2) or a)
        return seq[:a - 1] + seq[b:]
    m = re.fullmatch(r"(\d+)(?:_(\d+))?dup[ACGU]*", body)
    if m:
        a = int(m.group(1)); b = int(m.group(2) or a)
        return seq[:b] + seq[a - 1:b] + seq[b:]
    m = re.fullmatch(r"(\d+)_(\d+)ins([ACGU]+)", body)
    if m:
        a = int(m.group(1))
        return seq[:a] + m.group(3) + seq[a:]
    raise ValueError(hgvs)


def ddg_bind(u4_mut, u4=RNU4_2_RNA, u6=RNU6_RNA):
    """ddG of U4/U6 binding (kcal/mol), Crick definition: dG_bind = dG(U4&U6) - dG(U4) - dG(U6)."""
    import RNA
    md = RNA.md()
    md.temperature = 37.0

    def single(s):
        return RNA.fold_compound(s, md).mfe()[1]

    def dimer(a, b):
        return RNA.fold_compound(a + "&" + b, md).mfe_dimer()[1]
    g6 = single(u6)
    wt = dimer(u4, u6) - single(u4) - g6
    mt = dimer(u4_mut, u6) - single(u4_mut) - g6
    return mt - wt


def cmd_vienna(args):
    import RNA
    d = pd.read_csv(os.path.join(WORK, "sge.tsv"), sep="\t")
    t = d[d.in_transcript].copy()
    vals = []
    for h in t.HGVS_transcript:
        try:
            vals.append(ddg_bind(apply_hgvs_n(RNU4_2_RNA, h)))
        except (ValueError, AssertionError) as e:
            log(f"  unparsed {h}: {e}")
            vals.append(np.nan)
    t["ddG_bind_ours"] = vals
    both = t[t.ddG_bind.notna() & t.ddG_bind_ours.notna()]
    diff = (both.ddG_bind - both.ddG_bind_ours).abs()
    r = np.corrcoef(both.ddG_bind, both.ddG_bind_ours)[0, 1]
    txt = (f"ViennaRNA {RNA.__version__} ddG_bind recomputation vs SGE table column: n={len(both)}, Pearson r {r:.4f}, "
           f"max |diff| {diff.max():.3f} kcal/mol, share within 0.05: {(diff <= 0.05).mean():.1%}")
    log(txt)
    t[["ID", "HGVS_transcript", "ddG_bind", "ddG_bind_ours"]].to_csv(os.path.join(WORK, "sge_vienna.tsv"), sep="\t", index=False)
    with open(os.path.join(RESULTS, "week1_vienna_qc.txt"), "w") as f:
        f.write(txt + "\n")
    run_info("vienna")


# ------------------------------------------------------------------------------------------------- s6.1 Evo 2

EVO2_WINDOW = 8192


def cmd_evo2_prep(args):
    d = pd.read_csv(os.path.join(WORK, "sge.tsv"), sep="\t")
    t = d[d.in_transcript]
    half = EVO2_WINDOW // 2
    lo, hi = int(t.pos.min()) - half - 10, int(t.pos.max()) + half + 10
    ref = ucsc_seq("chr12", lo - 1, hi)  # 1-based lo..hi

    def sub(a, b):  # 1-based inclusive
        return ref[a - lo: b - lo + 1]
    rows = []
    for r in t.itertuples():
        c = r.pos  # first reference base of the VCF-style allele
        a, b = c - half, c + half - 1
        rseq = sub(a, b)
        off = c - a
        assert rseq[off:off + len(r.ref)] == r.ref, r.ID
        aseq = rseq[:off] + r.alt + rseq[off + len(r.ref):]
        rows.append({"ID": r.ID, "ref_seq": rseq, "alt_seq": aseq})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(WORK, "evo2_inputs.tsv.gz"), sep="\t", index=False)
    log(f"evo2 inputs: {len(out)} variants, window {EVO2_WINDOW} nt -> {os.path.join(WORK, 'evo2_inputs.tsv.gz')}")
    log("Run on a GPU (>= 24 GB): python evo2_score.py work/evo2_inputs.tsv.gz work/evo2_scores.tsv")


# ------------------------------------------------------------------------------------------------ metrics (s6.3)

def auroc(y, s):
    """AUROC of score s for binary y (Mann-Whitney, ties count 1/2)."""
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    n1, n0 = y.sum(), (~y).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = stats.rankdata(s)
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def spearman(a, b):
    return stats.spearmanr(a, b).statistic


def boot_auroc(y, s, B=2000, seed=SEED):
    """Class-stratified bootstrap 95% interval for AUROC."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
    out = np.empty(B)
    for b in range(B):
        idx = np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))])
        out[b] = auroc(y[idx], s[idx])
    return np.nanpercentile(out, [2.5, 97.5])


def boot_spearman(a, b, B=2000, seed=SEED, clusters=None):
    """Bootstrap 95% interval for Spearman rho, resampling variants or (if given) clusters of variants."""
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    out = np.empty(B)
    if clusters is None:
        for k in range(B):
            idx = rng.integers(0, len(a), len(a))
            out[k] = spearman(a[idx], b[idx])
    else:
        clusters = np.asarray(clusters)
        uc = np.unique(clusters)
        members = {c: np.flatnonzero(clusters == c) for c in uc}
        for k in range(B):
            idx = np.concatenate([members[c] for c in rng.choice(uc, len(uc))])
            out[k] = spearman(a[idx], b[idx])
    return np.nanpercentile(out, [2.5, 97.5])


def position_floor(t):
    """Leave-one-position-out OLS prediction of function score from region one-hot + log1p(distance to C)."""
    regions = [name for name, _ in RNU4_2_REGIONS] + ["Other"]
    X = np.column_stack([(t.region == r).astype(float) for r in regions]
                        + [np.log1p(np.array([min(abs(n - c) for c in CRITICAL) for n in t.n]))])
    y = t.function_score_mean.to_numpy(float)
    pos = t.n.to_numpy()
    pred = np.empty(len(t))
    for p in np.unique(pos):
        test = pos == p
        coef, *_ = np.linalg.lstsq(X[~test], y[~test], rcond=None)
        pred[test] = X[test] @ coef
    return pred


def damage_scores(t, evo2=None):
    """Higher = more damaging, for every comparator in PREREG s6.2."""
    sc = {"CADD v1.7 (paper column)": t.CADD_score.to_numpy(float),
          "ViennaRNA U4/U6 ddG_bind (paper column)": t.ddG_bind.to_numpy(float),
          "Allele-count floor": -np.log10(t.gnomad_AC.to_numpy(float) + 1),
          "Position-only floor (LOPO)": -position_floor(t)}
    if evo2 is not None:
        sc["Evo 2 7B zero-shot"] = evo2
    return sc


def g1_metrics(t, scores, B=2000):
    rows = []
    fs = t.function_score_mean.to_numpy(float)
    dep = t.depleted.to_numpy(bool)
    for name, d in scores.items():
        ok = np.isfinite(d)
        rho = spearman(-d[ok], fs[ok])
        ci = boot_spearman(-d[ok], fs[ok], B)
        ci_pos = boot_spearman(-d[ok], fs[ok], B, clusters=t.n.to_numpy()[ok])
        auc = auroc(dep[ok], d[ok])
        aci = boot_auroc(dep[ok], d[ok], B)
        snv = ok & t.is_snv.to_numpy(bool)
        rows.append({"predictor": name, "n": int(ok.sum()), "spearman": rho, "spearman_lo": ci[0], "spearman_hi": ci[1],
                     "spearman_poscluster_lo": ci_pos[0], "spearman_poscluster_hi": ci_pos[1],
                     "auroc_depleted": auc, "auroc_lo": aci[0], "auroc_hi": aci[1],
                     "spearman_snv": spearman(-d[snv], fs[snv]), "auroc_depleted_snv": auroc(dep[snv], d[snv]),
                     "spearman_diploid": spearman(-d[ok & t.function_score_dip.notna().to_numpy()],
                                                  t.function_score_dip.to_numpy(float)[ok & t.function_score_dip.notna().to_numpy()])})
    return pd.DataFrame(rows)


def g1_decision(m, evo2_name="Evo 2 7B zero-shot"):
    lines = []
    floor = m[m.predictor == "Position-only floor (LOPO)"].spearman.iloc[0]
    best = m[~m.predictor.str.contains("floor")].sort_values("spearman", ascending=False).iloc[0]
    flag = floor >= best.spearman - 0.05
    lines.append(f"Position-only floor rho {floor:.3f} vs best non-floor predictor ({best.predictor}) {best.spearman:.3f}: "
                 + ("WITHIN 0.05 -> deep-model framing inside RNU4-2 ends (PLAN.md)" if flag else "not within 0.05"))
    if evo2_name in set(m.predictor):
        e = m[m.predictor == evo2_name].iloc[0]
        fires = e.spearman >= 0.60 and e.auroc_depleted >= 0.85
        verdict = "G1 FIRES (Evo 2 is enough): drop the new-model framing" if fires else "G1 PASSES: continue"
        lines.append(f"Evo 2: rho {e.spearman:.3f} (need >= 0.60), AUROC {e.auroc_depleted:.3f} (need >= 0.85) -> {verdict}")
    else:
        verdict = "G1 PENDING: Evo 2 scores not available (run evo2_score.py on a GPU)"
        lines.append(verdict)
    return verdict, lines


def cmd_g1(args):
    ensure_dirs()
    d = pd.read_csv(os.path.join(WORK, "sge.tsv"), sep="\t")
    t = d[d.in_transcript].reset_index(drop=True)
    evo2 = None
    ep = os.path.join(WORK, "evo2_scores.tsv")
    if os.path.exists(ep):
        e = pd.read_csv(ep, sep="\t")
        delta = ((e.ll_alt_fwd - e.ll_ref_fwd) + (e.ll_alt_rc - e.ll_ref_rc)) / 2
        evo2 = -t.ID.map(dict(zip(e.ID, delta))).to_numpy(float)
    scores = damage_scores(t, evo2)
    m = g1_metrics(t, scores, args.boot)
    m.to_csv(os.path.join(RESULTS, "g1_metrics.tsv"), sep="\t", index=False, float_format="%.4f")

    # per-region rho
    reg = []
    for name, s in scores.items():
        for r, sub in t.groupby("region"):
            ok = np.isfinite(s[sub.index]) & (len(sub) >= 10)
            if ok.sum() >= 10:
                reg.append({"predictor": name, "region": r, "n": int(ok.sum()),
                            "spearman": spearman(-s[sub.index][ok], sub.function_score_mean.to_numpy(float)[ok])})
    pd.DataFrame(reg).to_csv(os.path.join(RESULTS, "g1_by_region.tsv"), sep="\t", index=False, float_format="%.4f")

    # s6.4 ReNU vs population
    snv = t[t.is_snv].copy()
    renu = snv.category == "ReNU_syndrome"
    prim = renu | ((~renu) & (snv.gnomad_AC > 0))
    paper = renu | (snv.category == "UKBB/AllofUs")
    pos_rows = []
    sge_score = -t.function_score_mean.to_numpy(float)
    allsc = dict(scores)
    allsc["SGE function score (truth reference)"] = sge_score
    for name, s in allsc.items():
        s_snv = s[snv.index]
        for setname, mask in (("primary: gnomAD v4.1 PASS AC>=1 controls", prim), ("reproduction: paper's 346 UKB/AoU SNVs", paper)):
            if name == "Allele-count floor" and setname.startswith("primary"):
                continue
            ok = mask.to_numpy() & np.isfinite(s_snv)
            y = renu.to_numpy()[ok]
            a = auroc(y, s_snv[ok])
            ci = boot_auroc(y, s_snv[ok], args.boot)
            pos_rows.append({"predictor": name, "control_set": setname, "n_renu": int(y.sum()), "n_controls": int((~y).sum()),
                             "auroc": a, "lo": ci[0], "hi": ci[1]})
    pr = pd.DataFrame(pos_rows)
    pr.to_csv(os.path.join(RESULTS, "g1_renu_vs_population.tsv"), sep="\t", index=False, float_format="%.4f")

    rep = pr[pr.control_set.str.startswith("reproduction")].set_index("predictor").auroc
    repro = [("CADD v1.7 (paper column)", 0.65), ("ViennaRNA U4/U6 ddG_bind (paper column)", 0.72),
             ("SGE function score (truth reference)", 0.93)]
    repro_ok = all(abs(round(rep[k], 2) - v) <= 0.01 for k, v in repro)
    verdict, dec = g1_decision(m)
    noise = t[["function_score_r1", "function_score_r2", "function_score_r3"]].corr().values[np.triu_indices(3, 1)]

    L = ["Gate G1 (PREREG s6), RNU4-2 SGE, 521 in-transcript variants (435 SNVs, 86 indels)",
         f"Labels: haploid function_score_mean; depleted < {SGE_DEPLETED} (n={int(t.depleted.sum())}); "
         f"replicate noise ceiling r = {noise.min():.2f}-{noise.max():.2f}",
         f"Bootstrap: {args.boot} replicates, seed {SEED}; rho interval resamples variants (position-cluster interval in brackets)",
         "",
         f"{'predictor':42s} {'rho [95% CI]':>24s} {'pos-cluster CI':>17s} {'AUROC depleted [95% CI]':>27s} {'rho SNV':>8s} {'AUC SNV':>8s}"]
    for r in m.itertuples():
        L.append(f"{r.predictor:42s} {r.spearman:7.3f} [{r.spearman_lo:6.3f},{r.spearman_hi:6.3f}] "
                 f"[{r.spearman_poscluster_lo:6.3f},{r.spearman_poscluster_hi:6.3f}] "
                 f"{r.auroc_depleted:9.3f} [{r.auroc_lo:5.3f},{r.auroc_hi:5.3f}]   {r.spearman_snv:8.3f} {r.auroc_depleted_snv:8.3f}")
    L += ["", "ReNU SNVs vs population SNVs (s6.4)"]
    for r in pr.itertuples():
        L.append(f"  {r.predictor:42s} {r.control_set:44s} {r.n_renu:2d} vs {r.n_controls:3d}: AUROC {r.auroc:.3f} [{r.lo:.3f}, {r.hi:.3f}]")
    L += ["", "Reproduction of published AUROCs (paper's control set): " + ", ".join(f"{k.split(' (')[0]} {rep[k]:.3f} (published {v})" for k, v in repro)
          + (" -> OK" if repro_ok else " -> MISMATCH: pipeline error, stop"), ""]
    L += ["Decision"] + ["  " + x for x in dec]
    txt = "\n".join(L)
    with open(os.path.join(RESULTS, "g1_report.txt"), "w") as f:
        f.write(txt + "\n")
    log(txt)
    run_info("g1", {"verdict": verdict, "reproduction_ok": bool(repro_ok)})


# ------------------------------------------------------------------------------------------ s7 G2 statistics core

def fit_kappa(obs, m, lo=1e-6, hi=1e6):
    """ML scale k for P(observed) = 1 - exp(-k m) (Bernoulli), by bisection on the score equation."""
    obs = np.asarray(obs, float)
    m = np.asarray(m, float)

    def score(k):
        q = np.exp(-k * m)
        p = np.clip(1 - q, 1e-12, 1)
        return np.sum(obs * m * q / p - (1 - obs) * m)
    a, b = np.log(lo), np.log(hi)
    for _ in range(200):
        c = (a + b) / 2
        if score(np.exp(c)) > 0:
            a = c
        else:
            b = c
    return float(np.exp((a + b) / 2))


def poisson_binomial_cdf(k, p):
    """P(X <= k) for X = sum of independent Bernoulli(p_i), by dynamic programming."""
    dist = np.zeros(len(p) + 1)
    dist[0] = 1.0
    for i, pi in enumerate(p):
        dist[1:i + 2] = dist[1:i + 2] * (1 - pi) + dist[:i + 1] * pi
        dist[0] *= 1 - pi
    return float(dist[:int(k) + 1].sum())


def window_table(n, obs, p, width=5):
    """Per-window O, E, O/E and Poisson-binomial lower tail D_w; n = transcript position of each site."""
    n = np.asarray(n)
    L = int(n.max())
    rows = []
    for s in range(1, L - width + 2):
        sel = (n >= s) & (n < s + width)
        o, e = obs[sel].sum(), p[sel].sum()
        rows.append({"start": s, "end": s + width - 1, "O": o, "E": e, "oe": o / e if e > 0 else np.nan,
                     "D": poisson_binomial_cdf(o, p[sel])})
    w = pd.DataFrame(rows)
    w["pct"] = w.D.rank(method="average") / len(w)
    w["top_decile"] = w.pct <= 0.10
    return w


def critical_region_test(n, obs, p, region, width=5):
    """PREREG s7.4(a): O/E of the region vs the gene's median window O/E, exact one-sided test, plus circular-shift p."""
    n = np.asarray(n)
    w = window_table(n, obs, p, width)
    m_g = float(np.nanmedian(w.oe))
    sel = np.isin(n, region)
    oe = obs[sel].sum() / p[sel].sum()
    pval = poisson_binomial_cdf(obs[sel].sum(), np.minimum(1.0, m_g * p[sel]))
    L = int(n.max())
    perm = []
    for s in range(1, L):
        shifted = ((np.asarray(region) - 1 + s) % L) + 1
        ss = np.isin(n, shifted)
        perm.append(obs[ss].sum() / p[ss].sum())
    perm_p = (1 + np.sum(np.array(perm) <= oe)) / L
    return {"oe": oe, "median_window_oe": m_g, "p_exact": pval, "p_circular": perm_p,
            "recovered": bool(oe < m_g and pval < 0.01)}


def hotspot_test(n, obs, p, hotspot, width=5):
    """PREREG s7.4(b): top-decile window inside H and circular-shift permutation p of mean depletion percentile."""
    n = np.asarray(n)
    w = window_table(n, obs, p, width)
    L = int(n.max())
    a, b = hotspot

    def stat(lo, hi):
        if lo <= hi:
            inside = (w.start >= lo) & (w.end <= hi)
        else:  # wraps around the end of the transcript
            inside = (w.start >= lo) | (w.end <= hi)
        return w.pct[inside].mean(), bool(w.top_decile[inside].any())
    s_obs, top = stat(a, b)
    perm = []
    for s in range(1, L):
        lo, hi = (a - 1 + s) % L + 1, (b - 1 + s) % L + 1
        perm.append(stat(lo, hi)[0])
    perm = np.array(perm)
    p_perm = (1 + np.sum(perm[np.isfinite(perm)] <= s_obs)) / L
    return {"mean_pct": s_obs, "top_decile_window": top, "p_circular": p_perm, "recovered": bool(top and p_perm < 0.01)}


# ------------------------------------------------------------------------------------------- s10 calibration core

def o_vst(prior):
    """Odds of pathogenicity for Very Strong evidence at a given prior (Tavtigian 2018 framework)."""
    return (9 * (1 - prior) / prior) ** (4 / 3)


def lr_thresholds(prior):
    o = o_vst(prior)
    return {"Supporting": o ** (1 / 8), "Moderate": o ** (1 / 4), "Strong": o ** (1 / 2), "Very Strong": o}


def local_lr(scores, y, grid, frac=0.2, min_class=10):
    """Local LR at each grid point: nearest `frac` of labelled variants by score rank (>= min_class per class)."""
    scores = np.asarray(scores, float)
    y = np.asarray(y, bool)
    order = np.argsort(scores)
    s_sorted, y_sorted = scores[order], y[order]
    N, NP, NB = len(scores), y.sum(), (~y).sum()
    k0 = min(N, max(int(np.ceil(frac * N)), 2 * min_class))
    cp = np.concatenate([[0], np.cumsum(y_sorted)])  # positives in [lo, hi) = cp[hi] - cp[lo]
    out = np.empty(len(grid))
    for i, g in enumerate(grid):
        c = np.searchsorted(s_sorted, g)
        k = k0
        while True:  # widen until both classes have min_class members (or the window is everything)
            lo = max(0, min(c - k // 2, N - k))
            npos = cp[lo + k] - cp[lo]
            if (npos >= min(min_class, NP) and k - npos >= min(min_class, NB)) or k >= N:
                break
            k = min(N, k + 2)
        out[i] = ((npos + 0.5) / NP) / ((k - npos + 0.5) / NB)
    return out


def pp3_threshold(scores, y, prior=0.10, B=1000, seed=SEED, strength="Supporting"):
    """Smallest grid score above which the 5th bootstrap percentile of local LR clears the strength (PREREG s10)."""
    rng = np.random.default_rng(seed)
    scores = np.asarray(scores, float)
    y = np.asarray(y, bool)
    grid = np.unique(scores)
    i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
    boots = np.empty((B, len(grid)))
    for b in range(B):
        idx = np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))])
        boots[b] = local_lr(scores[idx], y[idx], grid)
    lo5 = np.percentile(boots, 5, axis=0)
    need = lr_thresholds(prior)[strength]
    ok = lo5 >= need
    # all grid points at or above s* must clear it
    suffix_ok = np.flip(np.cumprod(np.flip(ok)).astype(bool))
    return float(grid[np.argmax(suffix_ok)]) if suffix_ok.any() else None


# ---------------------------------------------------------------------------------------------------- s5 sealing

def cmd_seal(args):
    src = args.file
    if not os.path.exists(src):
        sys.exit(f"{src} not found")
    enc = src + ".enc"
    pw = ["-pass", f"env:{args.pass_env}"] if args.pass_env else []
    subprocess.run(["openssl", "enc", "-aes-256-cbc", "-pbkdf2", "-iter", "200000", "-salt", "-in", src, "-out", enc] + pw, check=True)
    rec = {"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), "plaintext": os.path.basename(src),
           "sha256_plaintext": sha256(src), "sha256_ciphertext": sha256(enc), "user": getpass.getuser()}
    with open(os.path.join(os.path.dirname(os.path.abspath(src)), "SEALED_MANIFEST.jsonl"), "a") as f:
        f.write(json.dumps(rec) + "\n")
    log(json.dumps(rec, indent=1))
    log("Record both SHA-256 values in a dated PREREG amendment, then delete the plaintext (shred -u) before any model is trained.")


def cmd_unseal(args):
    if not args.reason:
        sys.exit("--reason is required: unsealing is logged and happens once, at the final evaluation (PREREG s9)")
    enc = args.file
    out = enc[:-4] if enc.endswith(".enc") else enc + ".dec"
    pw = ["-pass", f"env:{args.pass_env}"] if args.pass_env else []
    subprocess.run(["openssl", "enc", "-d", "-aes-256-cbc", "-pbkdf2", "-iter", "200000", "-in", enc, "-out", out] + pw, check=True)
    rec = {"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), "event": "UNSEALED",
           "file": os.path.basename(enc), "sha256_plaintext": sha256(out), "user": getpass.getuser(), "reason": args.reason}
    with open(os.path.join(os.path.dirname(os.path.abspath(enc)), "SEALED_MANIFEST.jsonl"), "a") as f:
        f.write(json.dumps(rec) + "\n")
    log(json.dumps(rec, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("fetch", "genes", "sites", "sge", "vienna", "evo2-prep"):
        sub.add_parser(name)
    p = sub.add_parser("g1")
    p.add_argument("--boot", type=int, default=2000)
    p = sub.add_parser("seal")
    p.add_argument("file")
    p.add_argument("--pass-env", help="read the passphrase from this environment variable instead of prompting")
    p = sub.add_parser("unseal")
    p.add_argument("file")
    p.add_argument("--reason", required=False)
    p.add_argument("--pass-env", help="read the passphrase from this environment variable instead of prompting")
    args = ap.parse_args(argv)
    ensure_dirs()
    {"fetch": cmd_fetch, "genes": cmd_genes, "sites": cmd_sites, "sge": cmd_sge, "vienna": cmd_vienna,
     "evo2-prep": cmd_evo2_prep, "g1": cmd_g1, "seal": cmd_seal, "unseal": cmd_unseal}[args.cmd](args)


if __name__ == "__main__":
    main()
