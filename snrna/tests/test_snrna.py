"""Self-tests with known answers for snrna.py (run: python -m pytest -q tests). No network access needed."""
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import snrna  # noqa: E402

RNG = np.random.default_rng(7)


# ---------------------------------------------------------------- metrics

def test_auroc_matches_definition():
    y = np.array([1, 1, 0, 0, 0, 1, 0], bool)
    s = np.array([0.9, 0.4, 0.4, 0.1, 0.3, 0.8, 0.95])
    # pairs (pos, neg): count pos > neg, ties 1/2
    pos, neg = s[y], s[~y]
    want = np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg])
    assert snrna.auroc(y, s) == pytest.approx(want)
    assert snrna.auroc(y, -s) == pytest.approx(1 - want)


def test_bootstrap_intervals_cover_truth():
    n, rho = 500, 0.6
    z = RNG.multivariate_normal([0, 0], [[1, rho], [rho, 1]], n)
    lo, hi = snrna.boot_spearman(z[:, 0], z[:, 1], B=500)
    true_spearman = 6 / np.pi * np.arcsin(rho / 2)
    assert lo < true_spearman < hi
    y = RNG.random(n) < 0.3
    s = y * 1.0 + RNG.normal(0, 1, n)
    lo, hi = snrna.boot_auroc(y, s, B=500)
    from scipy.stats import norm
    assert lo < norm.cdf(1 / np.sqrt(2)) < hi


def test_position_floor_has_no_leakage():
    # scores that depend only on region are predictable; per-position noise must not be
    n = np.repeat(np.arange(1, 146), 3)
    t = pd.DataFrame({"n": n})
    t["region"] = t.n.map(snrna.rnu4_2_region)
    effect = {name: RNG.normal(0, 1) for name, _ in snrna.RNU4_2_REGIONS} | {"Other": 0.0}
    t["function_score_mean"] = t.region.map(effect) + RNG.normal(0, 0.1, len(t))
    pred = snrna.position_floor(t)
    assert snrna.spearman(pred, t.function_score_mean) > 0.8
    pos_noise = dict(zip(range(1, 146), RNG.normal(0, 1, 145)))
    t["function_score_mean"] = t.n.map(pos_noise)
    pred = snrna.position_floor(t)
    assert abs(snrna.spearman(pred, t.function_score_mean)) < 0.25


# ---------------------------------------------------------------- G2 statistics

def simulate_gene(L=145, kappa=3.0, depleted=(), factor=1.0, seed=0):
    rng = np.random.default_rng(seed)
    n = np.repeat(np.arange(1, L + 1), 3)
    m = rng.lognormal(np.log(0.1), 0.8, len(n))  # heterogeneous per-site rates, Roulette-like spread
    sel = np.where(np.isin(n, depleted), factor, 1.0)
    obs = (rng.random(len(n)) < 1 - np.exp(-kappa * m * sel)).astype(float)
    return n, m, obs


def test_poisson_binomial_matches_binomial():
    from scipy.stats import binom
    p = np.full(30, 0.2)
    for k in (0, 3, 6, 12, 30):
        assert snrna.poisson_binomial_cdf(k, p) == pytest.approx(binom.cdf(k, 30, 0.2), rel=1e-9)


def test_fit_kappa_recovers_scale():
    n, m, obs = simulate_gene(L=2000, kappa=3.0, seed=1)
    assert snrna.fit_kappa(obs, m) == pytest.approx(3.0, rel=0.1)


def test_critical_region_test_null_and_power():
    region = snrna.CRITICAL
    null_hits, power_hits = 0, 0
    sims = 300
    for s in range(sims):
        n, m, obs = simulate_gene(seed=100 + s)
        p = 1 - np.exp(-snrna.fit_kappa(obs, m) * m)
        null_hits += snrna.critical_region_test(n, obs, p, region)["recovered"]
        n, m, obs = simulate_gene(depleted=region, factor=0.05, seed=1000 + s)
        p = 1 - np.exp(-snrna.fit_kappa(obs, m) * m)
        power_hits += snrna.critical_region_test(n, obs, p, region)["recovered"]
    # nominal one-sided alpha 0.01: allow Monte Carlo slack, but a badly anti-conservative test fails
    assert null_hits / sims <= 0.03
    assert power_hits / sims >= 0.8


def test_hotspot_test_null_and_power():
    hot = (37, 44)
    null_hits, power_hits = 0, 0
    sims = 200
    for s in range(sims):
        n, m, obs = simulate_gene(L=117, seed=5000 + s)
        p = 1 - np.exp(-snrna.fit_kappa(obs, m) * m)
        null_hits += snrna.hotspot_test(n, obs, p, hot)["recovered"]
        n, m, obs = simulate_gene(L=117, depleted=range(37, 45), factor=0.02, seed=9000 + s)
        p = 1 - np.exp(-snrna.fit_kappa(obs, m) * m)
        power_hits += snrna.hotspot_test(n, obs, p, hot)["recovered"]
    assert null_hits / sims <= 0.03
    assert power_hits / sims >= 0.8  # amended rule A1; the frozen rule had ~0.57 here


# ---------------------------------------------------------------- calibration

def test_evidence_thresholds_match_published():
    assert snrna.o_vst(0.10) == pytest.approx(350, rel=0.01)  # Tavtigian 2018
    t = snrna.lr_thresholds(0.0441)  # Pejaver 2022: 2.406, 5.790, 33.53, 1124
    assert t["Supporting"] == pytest.approx(2.406, abs=0.01)
    assert t["Moderate"] == pytest.approx(5.790, abs=0.05)
    assert t["Strong"] == pytest.approx(33.53, rel=0.01)
    assert t["Very Strong"] == pytest.approx(1124, rel=0.01)


def test_pp3_threshold_is_conservative_and_finite():
    rng = np.random.default_rng(3)
    y = np.r_[np.ones(400, bool), np.zeros(4000, bool)]
    s = np.r_[rng.normal(2, 1, 400), rng.normal(0, 1, 4000)]
    # true LR(s) = exp(2s - 2); Supporting at prior 0.1 needs LR >= 2.08 -> s >= 1.37
    thr = snrna.pp3_threshold(s, y, prior=0.10, B=200)
    assert thr is not None and 1.2 <= thr <= 3.0
    # no signal -> no threshold
    s0 = rng.normal(0, 1, len(y))
    assert snrna.pp3_threshold(s0, y, prior=0.10, B=200) is None


# ---------------------------------------------------------------- variant helpers

def test_hgvs_edits():
    seq = "AGCUUUGCG"
    assert snrna.apply_hgvs_n(seq, "n.2G>A") == "AACUUUGCG"
    assert snrna.apply_hgvs_n(seq, "NR_003137.3:n.3C>T") == "AGUUUUGCG"
    assert snrna.apply_hgvs_n(seq, "n.4del") == "AGCUUGCG"
    assert snrna.apply_hgvs_n(seq, "n.3dup") == "AGCCUUUGCG"
    assert snrna.apply_hgvs_n(seq, "n.2_3insUU") == "AGUUCUUUGCG"
    with pytest.raises(AssertionError):
        snrna.apply_hgvs_n(seq, "n.1C>A")


def test_normalise_left_aligns_indels():
    ref = "GATTTTC"  # positions 1..7

    def base(c, p):
        return ref[p - 1]
    # deletion of one T written at the right end of the run -> leftmost representation
    assert snrna.normalise("x", 5, "TT", "T", base) == (2, "AT", "A")
    assert snrna.normalise("x", 1, "G", "A", base) == (1, "G", "A")


def test_ddg_bind_is_zero_for_wild_type():
    pytest.importorskip("RNA")
    assert snrna.ddg_bind(snrna.RNU4_2_RNA) == pytest.approx(0.0, abs=1e-9)


def test_seal_roundtrip(tmp_path):
    f = tmp_path / "blind_labels.tsv"
    f.write_text("chrom\tpos\tref\talt\nchr11\t1\tA\tG\n")
    env = dict(os.environ, SEALPW="test-passphrase")
    py = [sys.executable, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "snrna.py")]
    subprocess.run(py + ["seal", str(f), "--pass-env", "SEALPW"], check=True, env=env, capture_output=True)
    plain = f.read_text()
    f.unlink()
    r = subprocess.run(py + ["unseal", str(f) + ".enc", "--pass-env", "SEALPW"], env=env, capture_output=True, text=True)
    assert r.returncode != 0  # --reason is mandatory
    subprocess.run(py + ["unseal", str(f) + ".enc", "--pass-env", "SEALPW", "--reason", "test"], check=True, env=env,
                   capture_output=True)
    assert f.read_text() == plain
    log = (tmp_path / "SEALED_MANIFEST.jsonl").read_text().strip().splitlines()
    assert len(log) == 2 and "UNSEALED" in log[1]


# ---------------------------------------------------------------- Evo 2 scoring loop (fake scorer, no GPU)

def test_evo2_run_scores_four_strands_and_resumes(tmp_path):
    import evo2_score
    inp = tmp_path / "in.tsv.gz"
    pd.DataFrame({"ID": ["v1", "v2", "v3"], "ref_seq": ["ACGTAC", "GGGTTT", "AAAACC"],
                  "alt_seq": ["ACTTAC", "GGGGTTT", "AAACC"]}).to_csv(inp, sep="\t", index=False)
    calls = []

    def fake(seqs):
        calls.append(seqs)
        return [float(len(s)) + 0.1 * s.count("G") for s in seqs]
    out = tmp_path / "scores.tsv"
    assert evo2_score.run(str(inp), str(out), fake, "fake", limit=2, log=lambda *a: None) == 2
    assert evo2_score.run(str(inp), str(out), fake, "fake", log=lambda *a: None) == 1  # resumes
    s = pd.read_csv(out, sep="\t")
    assert list(s.ID) == ["v1", "v2", "v3"] and list(s.columns) == evo2_score.COLUMNS
    assert calls[0] == ["ACGTAC", "ACTTAC", "GTACGT", "GTAAGT"]  # ref, alt, rc(ref), rc(alt)
    assert s.n_alt.tolist() == [6, 7, 5]
