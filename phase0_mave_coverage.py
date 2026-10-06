"""Phase 0.1 (feasibility, counts only, no outcomes): MaveDB score sets mapped to DDG2P monoallelic genes.
Run in the WSL repo: python phase0_mave_coverage.py -> results_fullstudy/phase0_mave_coverage.tsv"""
import json, os, sys, urllib.request
import pandas as pd

sys.path.insert(0, ".")
import gate

API = "https://api.mavedb.org/api/v1"


def get(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)


dd = pd.read_csv(os.path.join("data", gate.DDG2P_FILE), dtype=str)
qual = (dd["DDD category"].str.lower().isin(["confirmed", "probable"])
        & dd["allelic requirement"].str.lower().str.contains("monoallelic", na=False))
syms = set(dd.loc[qual, "gene symbol"].str.upper())
m = get(f"{API}/score-sets/mapped-genes")
rows = []
for urn, genes in m.items():
    hit = [g for g in genes if g.upper() in syms]
    if not hit:
        continue
    s = get(f"{API}/score-sets/{urn}")
    tg = (s.get("targetGenes") or [{}])[0]
    seq = ((tg.get("targetSequence") or {}).get("sequence") or "")
    rows.append({"gene": hit[0], "urn": urn, "title": s.get("title", "")[:120], "n_variants": s.get("numVariants"),
                 "target_len": len(seq), "target_type": (tg.get("targetSequence") or {}).get("sequenceType"),
                 "published": s.get("publishedDate"), "superseded": s.get("supersededScoreSet") is not None})
x = pd.DataFrame(rows).sort_values(["gene", "n_variants"], ascending=[True, False])
x.to_csv("results_fullstudy/phase0_mave_coverage.tsv", sep="\t", index=False)
cur = x[~x.superseded]
print(f"DDG2P monoallelic genes: {len(syms)}; with any MaveDB score set: {cur.gene.nunique()}; score sets {len(cur)}")
g = cur.groupby("gene").agg(n_sets=("urn", "size"), max_variants=("n_variants", "max"), max_len=("target_len", "max"))
print(g.sort_values("max_variants", ascending=False).to_string())
