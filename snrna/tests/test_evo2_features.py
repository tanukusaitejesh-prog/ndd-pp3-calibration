"""Evo 2 feature script on the mock model: coordinates, batch invariance, OOM fallback, resume.
Run: python -m pytest -q tests/test_evo2_features.py"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import torch

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "tests"))
import evo2_features as E  # noqa: E402
from mock_evo2 import Evo2  # noqa: E402

INPUTS = os.path.join(HERE, "model_inputs")
pytestmark = pytest.mark.skipif(not os.path.exists(os.path.join(INPUTS, "windows.tsv.gz")), reason="model_inputs not built")
COMP = str.maketrans("ACGU", "UGCA")


def test_transcript_positions_map_to_the_right_window_bases():
    genes = pd.read_csv(os.path.join(INPUTS, "genes.tsv.gz"), sep="\t")
    wins = pd.read_csv(os.path.join(INPUTS, "windows.tsv.gz"), sep="\t").set_index("label")
    nodes = pd.read_csv(os.path.join(INPUTS, "nodes.tsv.gz"), sep="\t")
    nodes = nodes[nodes.molecule == "target"]
    W = len(wins.iloc[0]["seq"])
    for g in genes.itertuples():
        seq = wins.loc[g.label, "seq"].upper().replace("T", "U")
        tx = nodes[nodes.label == g.label].sort_values("n")
        assert len(tx) == g.length
        idx = np.array([E.transcript_window_index(g, n, W) for n in tx.n])
        got = "".join(seq[i] for i in idx)
        if g.strand == "-":
            got = got.translate(COMP)
        assert got == "".join(tx.base), g.label
        # the reverse-complement window index used for embeddings points at the same genomic base
        rc = E.revcomp(wins.loc[g.label, "seq"].upper())
        assert all(rc[W - 1 - i] == E.revcomp(wins.loc[g.label, "seq"].upper()[i]) for i in idx[:5])


def test_batched_scores_do_not_depend_on_batch_size_and_survive_oom():
    m = Evo2()
    seqs = ["ACGT" * 64, "TTGA" * 64, "GGCA" * 64, "CATG" * 64, "AAAA" * 64]
    ref = E.Scorer(m, "cpu", 1).ll_batch(seqs)
    assert np.allclose(E.Scorer(m, "cpu", 3).ll_batch(seqs), ref, atol=1e-3)

    class OOM:  # raises CUDA OOM for batches larger than 2, like a full GPU
        def __init__(self, net):
            self.net, self.calls = net, []

        def __call__(self, ids):
            self.calls.append(ids.shape[0])
            if ids.shape[0] > 2:
                raise torch.cuda.OutOfMemoryError("mock")
            return self.net(ids)

    m2 = Evo2()
    m2.model = OOM(m2.model)
    sc = E.Scorer(m2, "cpu", 8)
    got = sc.ll_batch(seqs)
    assert sc.batch == 2 and np.allclose(got, ref, atol=1e-3)
    assert set(m2.model.calls[-3:]) == {2}  # whole list rescored at the final batch size, last chunk padded


def test_full_script_mock_run_and_resume(tmp_path):
    env = {**os.environ, "SNRNA_EVO2_MOCK": "1"}
    cmd = [sys.executable, os.path.join(HERE, "evo2_features.py"), "--out", str(tmp_path), "--limit-genes", "1", "--batch", "4"]
    subprocess.run(cmd, check=True, env=env, capture_output=True)
    d = pd.read_csv(tmp_path / "deltas.tsv", sep="\t")
    var = pd.read_csv(os.path.join(INPUTS, "variants.tsv.gz"), sep="\t")
    g0 = pd.read_csv(os.path.join(INPUTS, "genes.tsv.gz"), sep="\t").iloc[0]
    assert set(d.ID) == set(var.ID[var.label == g0.label]) and d.ID.is_unique
    assert d.ll_ref_fwd.nunique() == 1 and np.isfinite(d.drop(columns="ID").to_numpy()).all()
    emb = np.load(tmp_path / "emb" / f"{g0.label.replace('|', '_')}.npy")
    assert emb.shape == (g0.length, 4096) and emb.dtype == np.float16
    # a re-run scores nothing new and leaves the file unchanged
    before = (tmp_path / "deltas.tsv").read_text()
    out = subprocess.run(cmd, check=True, env=env, capture_output=True, text=True).stdout
    assert "0 variants to score" in out and (tmp_path / "deltas.tsv").read_text() == before
