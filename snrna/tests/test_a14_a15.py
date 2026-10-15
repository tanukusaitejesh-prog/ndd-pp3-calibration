"""Known-answer tests for a14_recalibration.py and a15_null_validity.py (no network, no torch)."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import a14_recalibration as A14  # noqa: E402
import a15_null_validity as A15  # noqa: E402


def test_cap_select_limits_genes_and_total():
    df = pd.DataFrame({"ID": range(500), "gene": ["TTN"] * 400 + [f"G{i}" for i in range(100)]})
    s = A14.cap_select(df)
    assert (s.gene == "TTN").sum() == A14.PER_GENE and len(s) == A14.PER_GENE + 100


def test_define_sets_rare_and_unfetched():
    d = pd.DataFrame({"cls": ["P/LP", "B/LB", "B/LB", "B/LB", "P/LP"], "cadd_phred": [30, 10, 12, np.nan, 25],
                      "AC": [np.nan, 3, 400, np.nan, 7], "filter": [np.nan, "PASS", "PASS", np.nan, "AC0"]})
    s = A14.define_sets(d)
    assert s["ALL"][0].tolist() == [True, True, True, False, True]          # row 3 not fetched
    assert s["RARE"][0].tolist() == [True, True, False, False, False]       # absent, AC 3 rare; AC 400 and non-PASS not
    assert s["ALL"][1].tolist() == [True, False, False, False, True]


def test_suffix_met():
    grid = np.array([0.1, 0.5, 0.8, 0.9, 0.99])
    lo5 = np.array([0.2, 1.0, 3.0, 2.0, 5.0])
    assert A14.suffix_met(grid, lo5, 1.5, 0.8) and not A14.suffix_met(grid, lo5, 2.5, 0.8)
    assert not A14.suffix_met(grid, lo5, 1.0, 1.5)


def test_calibration_shifts_when_benign_are_rarity_matched():
    import calibrate as C
    import snrna as S
    rng = np.random.default_rng(0)
    p = rng.normal(2, 1, 800)
    b_common, b_rare = rng.normal(-1, 1, 800), rng.normal(1, 1, 800)   # rare benign sit closer to pathogenic
    th = {}
    for name, b in (("ALL", b_common), ("RARE", b_rare)):
        s, y = np.r_[p, b], np.r_[np.ones(800, bool), np.zeros(800, bool)]
        grid, point, lo5, hi95 = C.lr_bootstrap(s, y, 100, seed=1)
        th[name] = C.thresholds(grid, lo5, hi95, 0.0441)["PP3_Supporting"]
    assert th["RARE"] is None or th["RARE"] > th["ALL"]


def test_simulation_known_answers():
    rng = np.random.default_rng(5)
    r = pd.DataFrame([A15.simulate_once(0, 1, 0.1, rng, n_win=400) for _ in range(5)]).mean()
    assert abs(r["T"] - 0.5) < 0.03 and r.raw > 0.7 and abs(r.matched - 0.5) < 0.03 and abs(r.corrected_plus_half - 0.5) < 0.04
    r = pd.DataFrame([A15.simulate_once(1, 1, 0.2, rng, n_win=400) for _ in range(5)]).mean()
    assert r.raw > r["T"] > r.matched and abs(r.matched_deattenuated - r["T"]) < 0.02


def test_window_bootstrap_metrics():
    rows = []
    for k in range(30):
        b = 10_000 * (k + 1)
        rows.append({"ID": f"p{k}", "chrom": "chr1", "pos": b, "filter": "PASS", "AC": 0, "x": 2.0})
        rows += [{"ID": f"r{k}_{j}", "chrom": "chr1", "pos": b + j - 10, "filter": "PASS", "AC": 1, "x": 0.0} for j in range(20) if j != 10]
        rows += [{"ID": f"n{k}_{j}", "chrom": "chr1", "pos": b + 50 + j, "filter": "PASS", "AC": 99, "x": -1.0} for j in range(3)]
    d = pd.DataFrame(rows)
    d, units, gneg = A15.window_index(d, {f"p{k}" for k in range(30)})
    assert len(units) == 30 and gneg is None
    assert all(len(u[1]) == 19 and len(u[2]) == 3 for u in units)                   # pool and negatives per window
    assert all(d.ID[u[1]].str.startswith("r").all() and d.ID[u[2]].str.startswith("n").all() for u in units)
    raw, null, corr, M = A15.metrics(d.x.to_numpy(float), units, np.arange(30), None, np.random.default_rng(0), 10)
    assert raw == 1.0 and null == 1.0 and corr == 0.0 and M == 1.0
