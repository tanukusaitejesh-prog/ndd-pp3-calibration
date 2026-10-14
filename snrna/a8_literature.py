#!/usr/bin/env python3
"""Amendment A8: a second blind test of the frozen v4 on literature-reported variants of the blind genes (PREREG s15).

RNU2-2 and RNU5B-1 were never in any v4 training set and the v4 scores are hash-frozen (A5). Their first blind test
used ClinVar (>= 1 star) + RNUdb labels sealed on 13 Oct 2026. Patient variants reported only in the literature
(2025-2026 cohort papers) form a second, independent positive set. Written before any literature variant was curated.

    python a8_literature.py map lit/a8_*.tsv     # any machine: HGVS n. -> GRCh38 with reference check
                                                 #   -> labels/A8_literature.tsv (+ _rejected, _indels); commit it
    python a8_literature.py eval --blind sealed/blind_labels.tsv    # DGX: refuses unless the A8 labels are committed

Positives (primary): literature P/LP or disease-causing SNVs (curate_literature.py wording rule) in RNU2-2 that are NOT
in the sealed blind P/LP set and NOT in any training label file. Negatives: L3 (gnomAD PASS AC >= 9), minus every
positive of either set. Scores: frozen v4 (primary), Rfam alt rarity, Rfam ref column freq, CADD, phyloP447, GPN-MSA.
Decision (pre-registered; needs >= 8 new positives, else descriptive only):
  MET if (i) v4 AUROC exceeds its rarity-only null (1,000 draws of as many random rare RNU2-2 SNVs) with one-sided
  p < 0.05, and (ii) v4 - CADD AUROC > 0 with the 95% class-stratified bootstrap interval excluding 0.
Secondary: the same for RNU5B-1 (>= 5 new positives); sealed + new positives pooled; A6 inheritance prediction run
unchanged (inheritance_blind_eval.py, frozen score and threshold) on every RNU2-2/RNU5B-1 P/LP SNV whose
inheritance the literature states (A6 was not evaluable at unsealing for lack of annotation).
"""
import argparse
import hashlib
import os
import re
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import curate_literature as CL  # noqa: E402

BLIND = ("RNU2-2", "RNU5B-1")
LAB = os.path.join(HERE, "labels", "A8_literature.tsv")
SEED = 20261016


def cmd_map(a):
    genes = pd.read_csv(os.path.join(HERE, "model_inputs", "genes.tsv.gz"), sep="\t")
    nodes = pd.read_csv(os.path.join(HERE, "model_inputs", "nodes.tsv.gz"), sep="\t")
    t = nodes[nodes.molecule == "target"]
    tx = {lab: dict(zip(x.n, x.base)) for lab, x in t.groupby("label")}
    g_of = genes.set_index("gene")
    rows, indels, rejected = [], [], []
    for p in a.files:
        for r in pd.read_csv(p, sep="\t", dtype=str).fillna("").to_dict("records"):
            gene = r["gene"].strip()
            cls = r.get("classification_as_stated", "").lower()
            if gene not in BLIND:
                rejected.append({**r, "reason": "A8 covers RNU2-2 and RNU5B-1 only"})
                continue
            if not any(w in cls for w in CL.POSITIVE_WORDS) or any(w in cls for w in CL.NOT_PER_VARIANT):
                rejected.append({**r, "reason": "not reported as P/LP or disease-causing"})
                continue
            g = g_of.loc[gene]
            hg = r.get("hgvs_n", "").strip().replace(" ", "")
            m = CL.HGVS_SNV.match(hg)
            if not m:
                indels.append({**r, "file": os.path.basename(p)})
                continue
            n, ref, alt = int(m.group(1)), m.group(2).replace("T", "U"), m.group(3).replace("T", "U")
            if tx[g.label].get(n) != ref:
                rejected.append({**r, "reason": f"reference {ref} at n.{n} but the transcript has {tx[g.label].get(n)}"})
                continue
            pos = g.start + n - 1 if g.strand == "+" else g.end - n + 1
            gref, galt = (CL.to_dna(ref), CL.to_dna(alt)) if g.strand == "+" else \
                (CL.to_dna(ref).translate(CL.COMP), CL.to_dna(alt).translate(CL.COMP))
            vid = f"{g.chrom}-{pos}-{gref}-{galt}"
            gm = re.match(r"^(?:chr)?(\w+)[:\-](?:g\.)?(\d+)", r.get("genomic_grch38", "").strip())  # chr11:62841806 C>T, chr15:g.65304713G>C, 11-62841806
            if gm and [gm.group(1), gm.group(2)] != [g.chrom.replace("chr", ""), str(pos)]:
                rejected.append({**r, "reason": f"source GRCh38 {r['genomic_grch38']} disagrees with the mapping {vid}"})
                continue
            inh = r.get("inheritance", "").strip().lower()
            rows.append({"ID": vid, "gene": gene, "source": f"literature:{r.get('source', '')}", "hgvs_n": hg,
                         "classification": r.get("classification_as_stated", ""), "label": "P/LP",
                         "inheritance": inh if inh in ("dominant", "recessive") else "unknown",
                         "zygosity": r.get("zygosity_observed", ""), "n_families": r.get("n_families", ""),
                         "insilico_or_conservation_used_for_classification": r.get("insilico_or_conservation_used_for_classification", "unknown"),
                         "first_publication_date": r.get("first_publication_date", "") or r.get("first_publication_year", ""),
                         "tally_date": time.strftime("%Y-%m-%d", time.gmtime())})
    out = pd.DataFrame(rows)
    if len(out):  # one row per variant; a variant reported with both modes keeps "unknown"
        inh = out.groupby("ID").inheritance.agg(lambda s: s.iloc[0] if s.nunique() == 1 else "unknown")
        src = out.groupby("ID").source.agg(lambda s: ";".join(sorted(set(s))))
        out = out.drop_duplicates("ID").assign(inheritance=lambda d: d.ID.map(inh), source=lambda d: d.ID.map(src))
    out.to_csv(LAB, sep="\t", index=False)
    pd.DataFrame(rejected).to_csv(LAB.replace(".tsv", "_rejected.tsv"), sep="\t", index=False)
    pd.DataFrame(indels).to_csv(LAB.replace(".tsv", "_indels.tsv"), sep="\t", index=False)
    print(f"A8 literature SNVs {len(out)}; indels {len(indels)}; rejected {len(rejected)}")
    if len(out):
        print(out.groupby(["gene", "inheritance"]).ID.nunique().to_string())
    print("NEXT: commit labels/A8_literature*.tsv and lit/a8_*.tsv, push, then `python a8_literature.py eval` on the DGX")


def auroc_ci(S, y, x, B):
    return S.auroc(y, x), *S.boot_auroc(y, x, B)


def cmd_eval(a):
    import evaluate as E
    import snrna as S
    tracked = subprocess.run(["git", "-C", HERE, "ls-files", "--error-unmatch", LAB], capture_output=True).returncode == 0
    dirty = subprocess.run(["git", "-C", HERE, "status", "--porcelain", "--", LAB], capture_output=True, text=True).stdout.strip()
    if not tracked or dirty:
        sys.exit("labels/A8_literature.tsv must be committed and unchanged before eval")
    sp = os.path.join(HERE, "results_v4", "v4_scores.tsv.gz")
    man = {l.split()[1]: l.split()[0] for l in open(os.path.join(HERE, "results_v4", "v4_manifest.txt")) if l.strip()}
    if man.get(os.path.relpath(sp, HERE)) != hashlib.sha256(open(sp, "rb").read()).hexdigest():
        sys.exit("v4_scores.tsv.gz does not match the frozen manifest")
    rng = np.random.default_rng(SEED)
    lit = pd.read_csv(LAB, sep="\t")
    blind = pd.read_csv(a.blind, sep="\t")
    blind = blind[blind.ID.str.match(E.SNV_ID)]
    sealed = set(blind.ID[blind.label == "P/LP"])
    train = E.load_l2(os.path.join(HERE, "labels"))
    trained = set(train.ID)
    var = pd.read_csv(os.path.join(HERE, "model_inputs", "variants.tsv.gz"), sep="\t")
    v = var[var.gene.isin(BLIND)].merge(pd.read_csv(sp, sep="\t")[["ID", "v4"]], on="ID", how="left")
    v = v.merge(pd.read_csv(os.path.join(HERE, "comparators", "comparators.tsv.gz"), sep="\t")[["ID", "phylop447", "gpn_msa"]],
                on="ID", how="left")
    v["rfam_alt_rarity"], v["rfam_ref_col_freq"] = 1 - v.alt_col_freq, v.ref_col_freq
    scores = [("v4", "v4 (frozen) [PRIMARY]"), ("rfam_alt_rarity", "Rfam alt rarity"), ("rfam_ref_col_freq", "Rfam ref column freq"),
              ("cadd_phred", "CADD"), ("phylop447", "phyloP447"), ("gpn_msa", "GPN-MSA")]
    allpos = sealed | set(lit.ID)
    rep, rows, decision = ["A8: second blind test of frozen v4 on literature-only variants of the blind genes", ""], [], None
    for gene in BLIND:
        g = v[v.gene == gene]
        new = set(lit.ID[(lit.gene == gene)]) - sealed - trained
        sets = {"NEW literature positives [primary]" if gene == "RNU2-2" else "NEW literature positives": new,
                "sealed + new pooled": (set(lit.ID[lit.gene == gene]) | sealed) & set(g.ID)}
        for sname, P in sets.items():
            pos = g.ID.isin(P)
            neg = g.L3_benign_like.astype(bool) & ~g.ID.isin(allpos)
            x = g[pos | neg]
            y = x.ID.isin(P).to_numpy()
            rep.append(f"{gene} | {sname}: {int(y.sum())} P vs {int((~y).sum())} L3")
            need = 8 if gene == "RNU2-2" else 5
            if y.sum() < 3:
                rep.append("  too few positives")
                continue
            pool = np.flatnonzero((~g.L3_benign_like.astype(bool) & ~g.ID.isin(allpos)).to_numpy())
            res = {}
            for c, lab in scores:
                s = x[c].to_numpy(float)
                s = np.where(np.isnan(s), np.nanmedian(s), s)
                au, lo, hi = auroc_ci(S, y, s, a.boot)
                full = g[c].to_numpy(float)
                full = np.where(np.isnan(full), np.nanmedian(full), full)
                negi = np.flatnonzero(neg.to_numpy())
                nulls = np.array([S.auroc(np.r_[np.ones(y.sum(), bool), np.zeros(len(negi), bool)],
                                          full[np.r_[rng.choice(pool, int(y.sum()), replace=False), negi]]) for _ in range(1000)])
                p_null = float((np.sum(nulls >= au) + 1) / (len(nulls) + 1))
                res[c] = (au, s)
                rows.append({"gene": gene, "set": sname, "score": lab, "n_P": int(y.sum()), "n_N": int((~y).sum()), "auroc": au,
                             "lo": lo, "hi": hi, "null_mean": float(nulls.mean()), "excess": au - float(nulls.mean()), "p_vs_null": p_null})
                rep.append(f"  {lab:28s} AUROC {au:.3f} [{lo:.3f}, {hi:.3f}]  rarity null {nulls.mean():.3f}  excess "
                           f"{au - nulls.mean():+.3f}  p {p_null:.4f}")
            i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
            d = [S.auroc(y[ix], res["v4"][1][ix]) - S.auroc(y[ix], res["cadd_phred"][1][ix])
                 for ix in (np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))] for _ in range(a.boot))]
            dlo, dhi = np.nanpercentile(d, [2.5, 97.5])
            dpt = res["v4"][0] - res["cadd_phred"][0]
            rep.append(f"  v4 - CADD {dpt:+.3f} [{dlo:+.3f}, {dhi:+.3f}]")
            if sname.endswith("[primary]"):
                pv = [r for r in rows if r["gene"] == gene and r["set"] == sname and r["score"].startswith("v4")][0]["p_vs_null"]
                if y.sum() < need:
                    decision = f"A8 decision: DESCRIPTIVE ONLY ({int(y.sum())} new RNU2-2 positives < {need})"
                else:
                    ok = pv < 0.05 and dlo > 0
                    decision = (f"A8 decision: {'MET' if ok else 'NOT MET'} (v4 excess over rarity null p {pv:.4f} [< 0.05 needed]; "
                                f"v4 - CADD {dpt:+.3f} [{dlo:+.3f}, {dhi:+.3f}] [interval > 0 needed]; {int(y.sum())} new positives)")
            rep.append("")
    if decision:
        rep.append(decision)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "results", "a8_results.tsv"), sep="\t", index=False, float_format="%.4f")
    # A6, unchanged code: every blind-gene P/LP SNV whose inheritance the literature states
    inh = lit[lit.inheritance.isin(["dominant", "recessive"])][["ID", "gene", "inheritance"]].assign(label="P/LP")
    ip = os.path.join(HERE, "results", "a8_inheritance_labels.tsv")
    inh.to_csv(ip, sep="\t", index=False)
    rep.append(f"\nA6 on literature-annotated blind-gene P/LP SNVs ({len(inh)}: "
               f"{int((inh.inheritance == 'dominant').sum())} dominant, {int((inh.inheritance == 'recessive').sum())} recessive):")
    os.makedirs(os.path.join(HERE, "results", "a8_a6"), exist_ok=True)
    a6 = subprocess.run([sys.executable, os.path.join(HERE, "inheritance_blind_eval.py"), "--blind", ip,
                         "--results", os.path.join(HERE, "results", "a8_a6")], capture_output=True, text=True)
    rep.append(a6.stdout.strip() or a6.stderr.strip()[-800:])
    open(os.path.join(HERE, "results", "a8_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("map")
    m.add_argument("files", nargs="+")
    e = sub.add_parser("eval")
    e.add_argument("--blind", required=True)
    e.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()
    {"map": cmd_map, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
