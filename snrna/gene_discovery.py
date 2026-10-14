#!/usr/bin/env python3
"""Gene-level constraint as a disease-gene finder (exploratory, post hoc; see PREREG s14 2026-10-14).

Question: would population constraint alone have ranked the snRNA genes found to cause disease since 2024 at the top
of the 55 analysis-set loci, and which loci are ranked next?

Per locus, on the `include` sites of variants.tsv.gz (P(observed) = 1 - exp(-k m), m = exp(log_mr), as in G2):
  * oe_gene : observed / expected SNVs over the whole locus, E from the pooled k of all loci
  * oe_min_window : lowest O/E over sliding windows of W transcript positions (hotspot constraint; the disease
    regions are short, e.g. RNU4-2 n.62-79)
  * rel_* : the same window search with k refitted within the locus (primary score). Between-locus O/E is confounded
    by locus mutability: the expressed snRNA genes are germline mutation hotspots (whole-gene O/E ~2-3 under the
    pooled k) while pseudogene copies look 'constrained' (O/E ~0.2-0.4); the first run, which ranked on the pooled
    scores, put every disease gene near the bottom (AUROC 0.08-0.16) and is kept in the output as gene_* / win_*.
  * z = (O - E) / sqrt(sum p(1-p)) for each, and an upper 95% bound of O/E
No patient label is read; disease status is gene-level public knowledge (DISEASE_GENES below, with first report year).
Evaluation: rank of each disease gene, AUROC (disease vs other loci) with a permutation p, for all disease genes,
NDD genes only, and genes first reported from 2024 on (the 'retrospective prediction' set).

    python gene_discovery.py [--inputs model_inputs] [--window 20] [--out results/gene_discovery]
    python gene_discovery.py --inputs model_inputs_v3 --out results/gene_discovery_v3     # with UK Biobank counts
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

# gene: (first report year, phenotype class). Public gene-level knowledge only.
DISEASE_GENES = {
    "RNU4ATAC": (2011, "NDD"),   # MOPD1/TALS, Roifman, Lowry-Wood
    "RNU12": (2017, "NDD"),      # SCAR33, CDAGS
    "RNU4-2": (2024, "NDD"),     # ReNU syndrome
    "RNU5B-1": (2025, "NDD"),
    "RNU2-2": (2025, "NDD"),
    "RNU6ATAC": (2025, "other"),  # autoimmune diabetes / growth-microcephaly
    "RNU6-1": (2025, "RP"),      # adRP (Quinodoz 2026)
    "RNU6-2": (2025, "RP"),
    "RNU6-8": (2025, "RP"),
    "RNU6-9": (2025, "RP"),
}


def fit_k(obs, m, lo=1e-8, hi=1e8):
    """ML scale k for P(observed) = 1 - exp(-k m), bisection on the score equation (as snrna.fit_kappa)."""
    a, b = np.log(lo), np.log(hi)
    for _ in range(200):
        c = (a + b) / 2
        q = np.exp(-np.exp(c) * m)
        p = np.clip(1 - q, 1e-12, 1)
        if np.sum(obs * m * q / p - (1 - obs) * m) > 0:
            a = c
        else:
            b = c
    return float(np.exp((a + b) / 2))


def oe_stats(obs, p):
    o, e = float(obs.sum()), float(p.sum())
    sd = float(np.sqrt(np.sum(p * (1 - p))))
    return {"O": o, "E": e, "oe": o / e if e > 0 else np.nan, "z": (o - e) / sd if sd > 0 else np.nan,
            "oe_upper95": (o + 1.645 * sd) / e if e > 0 else np.nan}


def min_window(g, w):
    """Lowest-O/E window of w consecutive transcript positions with E >= 2 (all sites at those positions)."""
    best = None
    ns = np.sort(g.n.unique())
    for start in ns:
        sub = g[(g.n >= start) & (g.n < start + w)]
        if sub.p.sum() < 2:
            continue
        st = oe_stats(sub.obs.to_numpy(), sub.p.to_numpy())
        if best is None or st["z"] < best["z"]:
            best = {**st, "start": int(start), "end": int(start + w - 1)}
    return best or {"O": np.nan, "E": np.nan, "oe": np.nan, "z": np.nan, "oe_upper95": np.nan, "start": -1, "end": -1}


def auroc(score_pos, score_neg):
    """P(score of a disease gene < score of another locus): lower constraint score = more constrained."""
    sp, sn = np.asarray(score_pos)[:, None], np.asarray(score_neg)[None, :]
    return float(np.mean((sp < sn) + 0.5 * (sp == sn)))


def evaluate(tab, score, genes, n_perm=20000, seed=0):
    rng = np.random.default_rng(seed)
    t = tab.dropna(subset=[score])
    pos = t[t.label.isin(genes)][score].to_numpy()
    neg = t[~t.label.isin(DISEASE_GENES)][score].to_numpy()  # other disease genes are left out of the negatives
    if len(pos) == 0 or len(neg) == 0:
        return {"n_pos": int(len(pos)), "n_neg": int(len(neg)), "auroc": None, "p_perm": None}
    a = auroc(pos, neg)
    pool = np.concatenate([pos, neg])
    perm = np.array([auroc(*np.split(rng.permutation(pool), [len(pos)])) for _ in range(n_perm)])
    return {"n_pos": int(len(pos)), "n_neg": int(len(neg)), "auroc": a, "p_perm": float((np.sum(perm >= a) + 1) / (n_perm + 1))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--window", type=int, default=20)
    ap.add_argument("--out", default=os.path.join(HERE, "results", "gene_discovery"))
    a = ap.parse_args()
    v = pd.read_csv(os.path.join(a.inputs, "variants.tsv.gz"), sep="\t")
    genes = pd.read_csv(os.path.join(a.inputs, "genes.tsv.gz"), sep="\t")
    v = v[(v.include == 1) & v.log_mr.notna()].copy()
    v["m"] = np.exp(v.log_mr)
    v["obs"] = (v.AC_pass > 0).astype(float)
    k = fit_k(v.obs.to_numpy(), v.m.to_numpy())
    v["p"] = 1 - np.exp(-k * v.m)
    rows = []
    for lab, g in v.groupby("label"):
        whole = oe_stats(g.obs.to_numpy(), g.p.to_numpy())
        win = min_window(g, a.window)
        # within-gene: k refitted on the locus itself, so locus-wide mutability (snRNA genes are transcription-
        # associated mutation hotspots) cancels; the score is the most depleted window relative to its own gene
        g = g.assign(p=1 - np.exp(-fit_k(g.obs.to_numpy(), g.m.to_numpy()) * g.m))
        rel = min_window(g, a.window)
        rows.append({"label": lab, "gene": g.gene.iloc[0], "n_sites": len(g), **{f"gene_{x}": y for x, y in whole.items()},
                     **{f"win_{x}": y for x, y in win.items()}, **{f"rel_{x}": y for x, y in rel.items()}})
    tab = pd.DataFrame(rows).merge(genes[["label", "family", "named", "sealed"]], on="label", how="left")
    tab["disease"] = tab.gene.map(lambda x: DISEASE_GENES.get(x, ("", ""))[1])
    tab["first_report"] = tab.gene.map(lambda x: DISEASE_GENES.get(x, ("", ""))[0])
    tab.loc[~tab.label.isin(DISEASE_GENES) & ~tab.gene.isin(DISEASE_GENES), "disease"] = ""
    tab["label"] = np.where(tab.gene.isin(DISEASE_GENES), tab.gene, tab.label)  # disease genes keyed by gene name
    for s in ("gene_z", "win_z", "rel_z"):
        tab[f"rank_{s}"] = tab[s].rank(method="min")  # 1 = most constrained
    tab = tab.sort_values("rel_z")
    sets = {"all_disease": list(DISEASE_GENES), "ndd": [g for g, (_, c) in DISEASE_GENES.items() if c == "NDD"],
            "reported_2024_on": [g for g, (y, _) in DISEASE_GENES.items() if y >= 2024],
            "ndd_2024_on": [g for g, (y, c) in DISEASE_GENES.items() if y >= 2024 and c == "NDD"]}
    present = set(tab.label)
    ev = {s: {sc: evaluate(tab, sc, [g for g in gs if g in present]) for sc in ("rel_z", "rel_oe", "gene_z", "win_z")}
          for s, gs in sets.items()}
    # power control: rel_z grows with the number of observed variants, which is higher at the expressed (named) genes
    # than at the pseudogene copies, so repeat the test among named loci only
    named_tab = tab[tab.named.fillna(False).astype(bool)]
    ev["named_only"] = {s: evaluate(named_tab, "rel_z", [g for g in gs if g in present]) for s, gs in sets.items()}
    candidates = tab[tab.disease == ""].head(10)[["label", "family", "rel_oe", "rel_z", "rel_start", "rel_end", "gene_oe"]]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    tab.to_csv(a.out + ".tsv", sep="\t", index=False, float_format="%.4g")
    info = {"inputs": a.inputs, "window": a.window, "pooled_k": k, "n_loci": int(len(tab)),
            "disease_genes_missing": sorted(set(DISEASE_GENES) - present), "evaluation": ev,
            "top10_candidates_not_yet_disease": candidates.to_dict("records"),
            "note": "Exploratory and post hoc: gene-level O/E was seen during G2 (results/g2_qc_by_gene.tsv). "
                    "No patient variant label is used."}
    with open(a.out + ".json", "w") as f:
        json.dump(info, f, indent=1, default=float)
    show = ["label", "disease", "first_report", "n_sites", "gene_oe", "rel_oe", "rel_z", "rel_start", "rel_end"]
    print(tab[show].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("\nAUROC (disease gene more constrained than other loci) and permutation p:")
    for s, d in ((s, d) for s, d in ev.items() if s != "named_only"):
        print(f"  {s:18s} " + "  ".join(f"{sc} {r['auroc']:.3f} (p {r['p_perm']:.3g}, n {r['n_pos']})" if r["auroc"] is not None
                                         else f"{sc} n/a" for sc, r in d.items()))
    print("  named loci only, rel_z: " + "  ".join(f"{s} {r['auroc']:.3f} (p {r['p_perm']:.3g}, n {r['n_pos']} vs {r['n_neg']})"
                                                for s, r in ev["named_only"].items()))
    print("missing from analysis set:", info["disease_genes_missing"])


if __name__ == "__main__":
    main()
