"""Model tests on the real feature graphs with SYNTHETIC counts (no real gnomAD counts are fitted here).
Run: python -m pytest -q tests/test_model.py   (about 5-10 min on CPU)"""
import os
import sys

import numpy as np
import pandas as pd
import pytest
import torch

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import model as M  # noqa: E402
import train as T  # noqa: E402

INPUTS = os.path.join(HERE, "model_inputs")
pytestmark = pytest.mark.skipif(not os.path.exists(os.path.join(INPUTS, "variants.tsv.gz")), reason="model_inputs not built")


def spearman(a, b):
    ra, rb = pd.Series(a).rank().to_numpy(), pd.Series(b).rank().to_numpy()
    return float(np.corrcoef(ra, rb)[0, 1])


@pytest.fixture(scope="module")
def data():
    return M.Data(INPUTS, None)


def planted(data, seed=3):
    nodes = pd.read_csv(os.path.join(INPUTS, "nodes.tsv.gz"), sep="\t")
    return M.planted_counts(data.var, nodes, data.genes, seed)


def test_planted_selection_is_recovered_on_held_out_genes(data):
    f_true, gamma, y = planted(data)
    data.v_y = torch.tensor(y, dtype=torch.float32)
    g = data.genes
    clusters = sorted(g.cluster.unique())
    rng = np.random.default_rng(1)
    rng.shuffle(clusters)
    test_cl = set(clusters[:10])
    val_cl = set(clusters[10:16])
    train = list(g.label[~g.cluster.isin(test_cl | val_cl)])
    val = list(g.label[g.cluster.isin(val_cl)])
    test = set(g.label[g.cluster.isin(test_cl)])
    cfg = {"d": 32, "layers": 2, "dropout": 0.1, "lr": 5e-3, "wd": 1e-4, "epochs": 300, "patience": 8, "eval_every": 5,
           "attention": False}
    model, info = M.train_one(data, train, val, cfg, seed=0)
    f, gam, _ = M.predict(model, data)
    v = data.var
    sel = v.label.isin(test).to_numpy() & (gamma > 0)
    rho = spearman(f[sel], f_true[sel])
    named = g.set_index("label").named.astype(bool)
    gam_named = np.mean([gam[data.gidx[x]] for x in g.label if named[x]])
    gam_unnamed = np.mean([gam[data.gidx[x]] for x in g.label if not named[x]])
    print(f"held-out rho {rho:.3f}; gate named {gam_named:.2f} vs unnamed {gam_unnamed:.2f}; best epoch {info['best_epoch']}")
    assert rho > 0.6
    assert gam_named > gam_unnamed


def test_attention_path_runs_and_learns(data):
    f_true, gamma, y = planted(data, seed=4)
    data.v_y = torch.tensor(y, dtype=torch.float32)
    labs = list(data.genes.label)
    cfg = {"d": 16, "layers": 1, "dropout": 0.0, "lr": 5e-3, "wd": 0.0, "epochs": 15, "patience": 100, "eval_every": 5}
    model, info = M.train_one(data, labs[:40], labs[40:], cfg, seed=0)
    hist_ll = [h[1] for h in info.get("history", [])] if "history" in info else None
    f, _, _ = M.predict(model, data)
    assert np.isfinite(f).all() and (f >= 0).all()


def test_nb_loglik_matches_reference():
    from math import lgamma, log
    y, mu, r = 4.0, 2.5, 3.0
    ref = lgamma(y + r) - lgamma(r) - lgamma(y + 1) + r * log(r / (r + mu)) + y * log(mu / (r + mu))
    got = float(M.nb_loglik(torch.tensor(y), torch.tensor(mu), torch.tensor(r)))
    assert got == pytest.approx(ref, abs=1e-5)


def test_alpha_is_profiled_out(data):
    """Adding a constant to every log rate in a gene must not change the likelihood (alpha absorbs it)."""
    f_true, gamma, y = planted(data, seed=5)
    data.v_y = torch.tensor(y, dtype=torch.float32)
    m = M.SelectionModel(data.d_node, data.d_var, data.d_gene, d=16, layers=1)
    tm = torch.ones(len(data.labels), dtype=torch.bool)
    with torch.no_grad():
        f, g = m(data)
        a, _ = M.gene_loglik(m, data, f, g, tm)
        old = data.v_logm.clone()
        data.v_logm = data.v_logm + 3.0
        b, _ = M.gene_loglik(m, data, f, g, tm)
        data.v_logm = old
    assert torch.allclose(a, b, atol=1e-4)


def test_function_head_recovers_linear_signal():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(300, 10))
    y = X[:, 0] * 2 - X[:, 3] + rng.normal(scale=0.3, size=300)
    groups = np.repeat(np.arange(6), 50)
    pred, info = M.function_head(X, y, groups, X)
    assert info["cv_r2"] > 0.8 and np.corrcoef(pred, y)[0, 1] > 0.9


def test_no_held_out_or_sealed_gene_is_ever_trained_on():
    genes = pd.read_csv(os.path.join(INPUTS, "genes.tsv.gz"), sep="\t")
    g, sealed, u4, rnu42 = T.gene_sets(genes)
    allg = list(g.index)
    strict = [x for x in allg if x not in sealed and x not in u4]
    loose = [x for x in allg if x not in sealed and x not in rnu42]
    assert not set(strict) & (sealed | u4)
    assert not set(loose) & (sealed | rnu42)
    assert {"RNU4-2", "RNU4-1", "RNU4ATAC"} <= set(g.loc[list(u4)].gene)
    for c in g.cluster.unique():
        held = set(g.index[g.cluster == c])
        train = [x for x in allg if x not in sealed and g.cluster[x] != c]
        assert not set(train) & held and not set(train) & sealed
    # RNU5B-1 is held out together with its near-identical paralogs
    c5 = g.loc[g.gene == "RNU5B-1", "cluster"].iloc[0]
    assert {"RNU5A-1", "RNU5B-1"} <= set(g.gene[g.cluster == c5])


def test_ablation_switches_change_inputs():
    full = M.Data(INPUTS, None)
    no_inter = M.Data(INPUTS, None, use_inter=False)
    no_cov = M.Data(INPUTS, None, use_covariation=False)
    assert no_inter.N < full.N and int((no_inter.etype == 2).sum()) == 0
    assert no_cov.d_node < full.d_node and no_cov.d_var < full.d_var
    rnd = M.Data(INPUTS, None, random_graph_seed=7)
    assert not torch.equal(rnd.etm, full.etm)
    shuf = M.Data(INPUTS, None, shuffle_counts_seed=7)
    by_gene = lambda d: d.var.groupby("label").AC_pass.sum()  # noqa: E731
    assert (by_gene(shuf) == by_gene(full)).all() and not torch.equal(shuf.v_y, full.v_y)
