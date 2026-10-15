"""Known-answer tests for a16_population_calibration.py and paper_tables.py (no network, no torch)."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import a16_population_calibration as A16  # noqa: E402
import paper_tables as PT  # noqa: E402


def test_reservoir_is_uniform():
    rng = np.random.default_rng(0)
    counts = np.zeros(100)
    for _ in range(2000):
        for x in A16.reservoir(range(100), 5, rng):
            counts[x] += 1
    assert abs(counts.mean() - 100) < 1e-9 and counts.std() < 20


def test_pu_correction_and_break_even():
    assert A16.pu_correct(10.0, 0.0) == 10.0
    # U = (1 - pi) B + pi P  with true LR 20 gives LR_obs = 20 / ((1 - pi) + pi * 20)
    pi, true = 0.02, 20.0
    obs = true / ((1 - pi) + pi * true)
    assert abs(float(A16.pu_correct(obs, pi)) - true) < 1e-9
    assert np.isinf(A16.pu_correct(60.0, 0.02))                       # 1 / 60 < 0.02: unbounded
    need, L = 33.6, 15.0
    pstar = (need / L - 1) / (need - 1)
    assert abs(float(A16.pu_correct(L, pstar)) - need) < 1e-6


def test_band_table_verdicts():
    rng = np.random.default_rng(1)
    sp = np.r_[np.full(500, 0.99), np.full(500, 0.5)]            # half of P in the band
    sn_low = np.r_[np.full(10, 0.99), np.full(9990, 0.1)]         # 0.1% of comparison in band: LR ~ 500
    sn_high = np.r_[np.full(2000, 0.99), np.full(8000, 0.1)]      # 20% in band: LR ~ 2.5
    need = {"Strong": 33.6}
    met = A16.band_table(sp, sn_low, [("Strong", 0.972, np.inf)], (0.0,), 200, rng, need, "c", "t")[0]
    notmet = A16.band_table(sp, sn_high, [("Strong", 0.972, np.inf)], (0.0,), 200, rng, need, "c", "t")[0]
    assert met["verdict"] == "DEMONSTRABLY MET" and notmet["verdict"] == "DEMONSTRABLY NOT MET"
    assert notmet["break_even_pi"] > 0 and met["break_even_pi"] < 0


def test_paper_tables_from_synthetic_results(tmp_path, monkeypatch):
    rows = []
    for tool, raw, m in (("cadd_phred", 0.83, 0.70), ("gpn_msa", 0.84, 0.57)):
        for metric, pt in (("raw_auroc", raw), ("null_mean", 0.7), ("corrected", raw - 0.7), ("rarity_matched_auroc", m)):
            rows.append({"analysis": "a11", "tool": tool, "metric": metric, "point": pt, "lo": pt - 0.01, "hi": pt + 0.01})
    rows.append({"analysis": "a11", "tool": "cadd_phred - gpn_msa", "metric": "diff_rarity_matched_auroc", "point": 0.13, "lo": 0.1, "hi": 0.16})
    pd.DataFrame(rows).to_csv(tmp_path / "a15_ci_a11.tsv", sep="\t", index=False)
    monkeypatch.setattr(PT, "RES", str(tmp_path))
    t2, diffs, forest = PT.table_rarity()
    assert len(t2) == 2 and len(forest) == 2
    assert diffs and "CADD - GPN-MSA +0.130" in diffs[0]
    assert "|" in PT.md_table(t2)
