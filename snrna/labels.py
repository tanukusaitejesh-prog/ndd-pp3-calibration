#!/usr/bin/env python3
"""Label sets L2 (patient P/LP) for evaluation and calibration, and sealing of the blind genes (PREREG s4-s5).

    python labels.py nonblind                         # labels/L2_nonblind.tsv for every non-sealed analysis-set gene
    SNRNA_SEAL_PASS=... python labels.py seal          # RNU2-2 + RNU5B-1 -> sealed/blind_labels.tsv.enc (never printed)

Sources (public): ClinVar via NCBI E-utilities (P/LP or B/LB with >= 1 review star), RNUdb curated variants
(RNU4-2, RNU2-2), and the ReNU SNVs of the SGE paper. Literature-only variants not yet in these resources are added
later by a curator (labels/L2_curated_additions.tsv, same columns) and are merged by evaluate.py.
"""
import argparse
import datetime
import getpass
import hashlib
import io
import json
import os
import subprocess
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
INPUTS = os.path.join(HERE, "model_inputs")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RNUDB = "https://rnudb.rarediseasegenomics.org/api/genes/{gene}/variants"
SEALED = ("RNU2-2", "RNU5B-1")
PLP = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BLB = {"Benign", "Likely benign", "Benign/Likely benign"}
NO_STAR = {"no assertion criteria provided", "no classification provided", "no assertion provided",
           "no classifications from unflagged records"}
# gene-level inheritance where the gene has a single established mode (guidance Table 1); otherwise per variant
GENE_MODE = {"RNU4ATAC": "recessive", "RNU6ATAC": "recessive", "RNU12": "recessive", "RNU6-1": "dominant",
             "RNU6-2": "dominant", "RNU6-8": "dominant", "RNU6-9": "dominant", "RNU5B-1": "dominant", "RNU5A-1": "dominant"}
COLUMNS = ["ID", "gene", "source", "classification", "review_status", "inheritance", "title"]


def get_json(url, tries=6):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "snrna-vep/2 (research)"})
            with urllib.request.urlopen(req, timeout=120) as r:
                d = json.load(r)
            if isinstance(d, dict) and "error" in d and "rate" in str(d["error"]).lower():
                raise RuntimeError(d["error"])
            return d
        except Exception as e:
            time.sleep(1.5 * (i + 1))
            last = e
    raise RuntimeError(f"{url}: {last}")


def clinvar(gene, start, end, chrom):
    """ClinVar records annotated to `gene` whose GRCh38 position lies in [start, end]."""
    term = urllib.parse.quote(f"{gene}[gene]")
    ids = get_json(f"{EUTILS}/esearch.fcgi?db=clinvar&term={term}&retmax=5000&retmode=json")["esearchresult"]["idlist"]
    rows = []
    for k in range(0, len(ids), 150):
        time.sleep(0.4)
        d = get_json(f"{EUTILS}/esummary.fcgi?db=clinvar&id={','.join(ids[k:k + 150])}&retmode=json")
        d = d.get("result", d)
        for uid in d.get("uids", []):
            r = d[uid]
            gc = r.get("germline_classification") or {}
            for vs in r.get("variation_set", []):
                spdi = vs.get("canonical_spdi") or ""
                if not spdi.startswith("NC_"):
                    continue
                acc, pos0, ref, alt = spdi.split(":")
                c = acc.split(".")[0].replace("NC_0000", "").lstrip("0")
                c = {"23": "X", "24": "Y"}.get(c, c)
                if f"chr{c}" != chrom:
                    continue
                pos = int(pos0) + 1
                if not (start <= pos <= end):
                    continue
                if len(ref) == 1 and len(alt) == 1:
                    vid = f"chr{c}-{pos}-{ref}-{alt}"
                else:
                    vid = f"chr{c}-{pos}-{ref or '-'}-{alt or '-'}"  # indel, SPDI-style (not scored by the model)
                rows.append({"ID": vid, "gene": gene, "source": f"ClinVar:{r.get('accession', uid)}",
                             "classification": gc.get("description", ""), "review_status": gc.get("review_status", ""),
                             "title": r.get("title", "")})
    return rows


def rnudb(gene):
    try:
        x = get_json(RNUDB.format(gene=gene))
    except Exception:
        return []
    rows = []
    for v in x if isinstance(x, list) else []:
        cs = v.get("clinical_significance")
        if cs in ("Pathogenic", "Likely Pathogenic"):
            dt = (v.get("disease_type") or "").lower()
            mode = "dominant" if "dominant" in dt or "retinitis" in dt else "recessive" if "recessive" in dt else "unknown"
            rows.append({"ID": v["id"], "gene": gene, "source": "RNUdb", "classification": cs.replace("Likely Pathogenic", "Likely pathogenic"),
                         "review_status": "RNUdb curated", "inheritance": mode, "title": v.get("hgvs") or ""})
    return rows


def collect(genes):
    """L2 (P/LP, >= 1 star or curated) and ClinVar B/LB (>= 1 star) for the given gene rows."""
    out = []
    for g in genes.itertuples():
        for r in clinvar(g.gene, g.start, g.end, g.chrom):
            stars_ok = r["review_status"] not in NO_STAR
            if r["classification"] in PLP and stars_ok:
                r["label"] = "P/LP"
            elif r["classification"] in BLB and stars_ok:
                r["label"] = "B/LB"
            else:
                continue
            r["inheritance"] = GENE_MODE.get(g.gene, "unknown") if r["label"] == "P/LP" else ""
            out.append(r)
        if g.gene in ("RNU4-2", "RNU2-2", "RNU5B-1"):  # RNUdb (13 Oct 2026: no entries for any other non-blind gene)
            for r in rnudb(g.gene):
                r["label"] = "P/LP"
                out.append(r)
    df = pd.DataFrame(out, columns=COLUMNS + ["label"])
    if len(df):
        # one row per variant and label; prefer curated inheritance over gene-level defaults
        df["pri"] = df.inheritance.map({"dominant": 0, "recessive": 0, "unknown": 1, "": 2}).fillna(2)
        df = df.sort_values(["ID", "label", "pri"]).drop_duplicates(["ID", "label"]).drop(columns="pri")
    return df


def cmd_nonblind(args):
    genes = pd.read_csv(os.path.join(INPUTS, "genes.tsv.gz"), sep="\t")
    named = genes[genes.named & ~genes.gene.isin(SEALED)]
    df = collect(named)
    # ReNU SNVs of the SGE paper (dominant, published) if not already present
    sge = pd.read_csv(os.path.join(INPUTS, "sge_eval.tsv.gz"), sep="\t")
    renu = sge[(sge.category == "ReNU_syndrome")]
    extra = pd.DataFrame({"ID": renu.ID, "gene": "RNU4-2", "source": "SGE paper (ReNU)", "classification": "Pathogenic",
                          "review_status": "published", "inheritance": "dominant", "title": renu.HGVS_transcript, "label": "P/LP"})
    df = pd.concat([df, extra[~extra.ID.isin(df[df.label == "P/LP"].ID)]], ignore_index=True)
    os.makedirs(os.path.join(HERE, "labels"), exist_ok=True)
    path = os.path.join(HERE, "labels", "L2_nonblind.tsv")
    df.to_csv(path, sep="\t", index=False)
    summ = df.groupby(["gene", "label"]).size().unstack(fill_value=0)
    snv = df[df.ID.str.count("-") == 3]
    print(f"wrote {path}: {len(df)} rows")
    print(summ.to_string())
    print("P/LP SNVs by gene and inheritance:")
    print(snv[snv.label == "P/LP"].groupby(["gene", "inheritance"]).ID.nunique().to_string())


def cmd_seal(args):
    """Fetch blind-gene labels straight into an encrypted file. Contents are never printed or written in clear."""
    pw_env = "SNRNA_SEAL_PASS"
    if not os.environ.get(pw_env):
        sys.exit(f"set {pw_env} (a passphrase the analyst does not know) before sealing")
    out_dir = os.path.join(HERE, "sealed")
    enc = os.path.join(out_dir, "blind_labels.tsv.enc")
    if os.path.exists(enc) and not args.force:
        sys.exit(f"{enc} exists; sealing happens once (use --force only before any model is trained)")
    genes = pd.read_csv(os.path.join(INPUTS, "genes.tsv.gz"), sep="\t")
    blind = genes[genes.gene.isin(SEALED)]
    df = collect(blind)
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False)
    plain = buf.getvalue().encode()
    os.makedirs(out_dir, exist_ok=True)
    subprocess.run(["openssl", "enc", "-aes-256-cbc", "-pbkdf2", "-iter", "200000", "-salt", "-out", enc, "-pass", f"env:{pw_env}"],
                   input=plain, check=True)
    rec = {"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), "event": "SEALED",
           "genes": list(SEALED), "n_rows": int(len(df)), "sha256_plaintext": hashlib.sha256(plain).hexdigest(),
           "sha256_ciphertext": hashlib.sha256(open(enc, "rb").read()).hexdigest(), "user": getpass.getuser(),
           "sources": "ClinVar E-utilities (P/LP, B/LB >= 1 star) + RNUdb curated P/LP"}
    with open(os.path.join(out_dir, "SEALED_MANIFEST.jsonl"), "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["nonblind", "seal"])
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    {"nonblind": cmd_nonblind, "seal": cmd_seal}[args.cmd](args)


if __name__ == "__main__":
    main()
