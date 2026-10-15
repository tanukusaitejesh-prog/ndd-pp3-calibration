"""Known-answer tests for a12_new_variants.py, a13_benchmarks.py and raritybench.py (no network, no torch)."""
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import a12_new_variants as A12  # noqa: E402
import a13_benchmarks as A13  # noqa: E402
import raritybench as RB  # noqa: E402
import snrna as S  # noqa: E402


def test_units_cover_every_gene_once():
    named = [g for d in A12.UNITS.values() if d["genes"] for g in d["genes"]]
    assert len(named) == len(set(named))
    assert A12.unit_of("RNU6-2") == "U1_RNU6_adRP" and A12.unit_of("RNU11") == "U4_other" and A12.unit_of("U6") == "U4_other"
    assert sum(d["primary"] for d in A12.UNITS.values()) == 3


def test_parse_clinvar_in_loci_only(tmp_path):
    genes = pd.DataFrame({"label": ["RNU6-2"], "gene": ["RNU6-2"], "chrom": ["chr19"], "start": [1021522], "end": [1021628]})
    lines = ["##fileDate=2026-11-01", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
             "19\t1021530\t111\tC\tT\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter",
             "19\t1021531\t112\tC\tT\t.\t.\tCLNSIG=Likely_pathogenic;CLNREVSTAT=no_assertion_criteria_provided",
             "19\t1021540\t113\tG\tA\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter",
             "19\t1021541\t114\tGA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter",
             "19\t2000000\t115\tC\tT\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter"]
    p = tmp_path / "cv.vcf.gz"
    with gzip.open(p, "wt") as f:
        f.write("\n".join(lines) + "\n")
    d, fdate, n_indel = A12.parse_clinvar(str(p), genes)
    assert fdate == "2026-11-01" and n_indel == 1 and len(d) == 3
    r = d.set_index("ID")
    assert r.loc["chr19-1021530-C-T", "clinvar_class"] == "P/LP" and r.loc["chr19-1021530-C-T", "star1"]
    assert not r.loc["chr19-1021531-C-T", "star1"]
    assert r.loc["chr19-1021540-G-A", "clinvar_class"] == "B/LB"


def test_build_excludes_prior_and_sealed_and_conflicts():
    scan = pd.DataFrame({"ID": ["a", "b", "c", "d"], "gene": ["RNU6-2", "RNU6-7", "RNU2-2", "RNU6-2"],
                         "clinvar_class": ["P/LP", "P/LP", "P/LP", "B/LB"], "star1": [True, True, True, True],
                         "clnsig": ["Pathogenic"] * 4, "review_status": ["x"] * 4, "variation_id": [1, 2, 3, 4],
                         "clinvar_fileDate": ["2026-11-01"] * 4})
    lit = pd.DataFrame({"ID": ["e", "d"], "gene": ["RNU5A-1", "RNU6-2"], "source": ["literature:1", "literature:2"],
                        "classification": ["pathogenic", "pathogenic"], "tally_date": ["2026-11-02"] * 2})
    out, dropped = A12.build_labels([scan], lit, prior={"b"}, have_blind=False)
    assert set(out.ID) == {"a", "e"}                      # b prior, c sealed without --blind, d ClinVar B/LB
    assert set(dropped.ID) == {"c", "d"}
    assert dict(zip(out.ID, out.unit)) == {"a": "U1_RNU6_adRP", "e": "U2_U5"}
    out2, _ = A12.build_labels([scan], lit, prior={"b"}, have_blind=True)
    assert "c" in set(out2.ID)


def _unit_table(n=300, seed=0, signal=0.0):
    rng = np.random.default_rng(seed)
    g = pd.DataFrame({"ID": [f"v{i}" for i in range(n)], "L3_benign_like": rng.random(n) < 0.3})
    for c, _ in A12.SCORES:
        g[c] = rng.normal(size=n)
    return g, rng


def test_unit_metrics_null_and_signal():
    g, rng = _unit_table()
    P = set(g.ID[~g.L3_benign_like].iloc[:20])
    rows, delta = A12.unit_metrics(g, P, P, np.random.default_rng(1), 200, S.auroc, S.boot_auroc)
    v4 = rows[0]
    assert abs(v4["excess"]) < 0.12 and v4["p_vs_null"] > 0.01   # random scores: no excess
    g.loc[g.ID.isin(P), "v4"] += 3                                 # planted signal in v4 only
    rows, delta = A12.unit_metrics(g, P, P, np.random.default_rng(1), 200, S.auroc, S.boot_auroc)
    assert rows[0]["excess"] > 0.3 and rows[0]["p_vs_null"] < 0.01 and delta[1] > 0
    assert A12.decide("U1_RNU6_adRP", rows, delta, 20).startswith("U1_RNU6_adRP: MET")
    assert "DESCRIPTIVE" in A12.decide("U1_RNU6_adRP", rows, delta, 5)


def test_a13_random_scores_give_no_inflation_and_rarity_scores_do():
    rng = np.random.default_rng(3)
    rows = []
    for k in range(60):  # 60 windows: 1 positive, 40 rare pool SNVs, 4 common negatives each
        c, base = "chr1", 10_000 * (k + 1)
        rows.append({"ID": f"p{k}", "chrom": c, "pos": base, "filter": "PASS", "AC": 0, "lab": 1})
        rows += [{"ID": f"r{k}_{j}", "chrom": c, "pos": base + j - 20, "filter": "PASS", "AC": 1, "lab": -1} for j in range(40) if j != 20]
        rows += [{"ID": f"n{k}_{j}", "chrom": c, "pos": base + 50 + j, "filter": "PASS", "AC": 500, "lab": 0} for j in range(4)]
    d = pd.DataFrame(rows)
    d["cadd_phred"] = rng.normal(size=len(d))                                    # random: inflation ~ 0
    d["remm"] = np.where(d.AC >= 9, 0.0, 1.0) + rng.normal(0, 0.3, len(d))       # tracks rarity: inflated
    bench = d[d.lab >= 0][["ID"]].assign(label=d.lab[d.lab >= 0].to_numpy())
    t, q = A13.evaluate_benchmark(d, bench, S.auroc, S.boot_auroc, 100, np.random.default_rng(0),
                                  tools=[("cadd_phred", "CADD"), ("remm", "ReMM")])
    t = t.set_index("tool")
    assert q["frac_neg_common"] == 1.0
    assert abs(t.loc["CADD", "inflation"]) < 0.08
    assert t.loc["ReMM", "inflation"] > 0.3 and abs(t.loc["ReMM", "corrected"]) < 0.1


def test_raritybench_collects(tmp_path):
    pd.DataFrame({"unit": ["pooled"], "tool": ["CADD v1.7"], "n_P": [10], "n_N": [20], "raw_auroc": [0.8], "null_mean": [0.65],
                  "inflation": [0.15], "corrected": [0.15]}).to_csv(tmp_path / "a9_results.tsv", sep="\t", index=False)
    pd.DataFrame({"tool": ["ReMM v0.4"], "n_P": [5], "n_N": [9], "raw_auroc": [0.7], "null_mean": [0.6], "inflation": [0.1],
                  "corrected": [0.1], "rarity_matched_auroc": [0.6]}).to_csv(tmp_path / "a13_traitgym_mendelian_results.tsv", sep="\t", index=False)
    t, missing = RB.collect(str(tmp_path))
    assert len(t) == 2 and "a10_results.tsv" in missing
    assert abs(t.share_from_rarity.iloc[0] - 0.5) < 1e-9
    assert t.benchmark.iloc[1].startswith("A13 traitgym_mendelian")
