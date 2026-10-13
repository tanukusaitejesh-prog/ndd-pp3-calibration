"""Homology transfer RNU4-2 -> RNU4ATAC: consistent mapping and correct use of the SGE table (no patient labels).
Run: python -m pytest -q tests/test_homology.py"""
import os

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS = os.path.join(HERE, "model_inputs")
pytestmark = pytest.mark.skipif(not os.path.exists(os.path.join(INPUTS, "homology_rnu4atac.tsv.gz")),
                                reason="homology table not built")


def test_mapping_is_collinear_and_bases_are_right():
    h = pd.read_csv(os.path.join(INPUTS, "homology_rnu4atac.tsv.gz"), sep="\t").dropna(subset=["h"])
    nodes = pd.read_csv(os.path.join(INPUTS, "nodes.tsv.gz"), sep="\t")
    t = nodes[nodes.molecule == "target"]
    b42 = dict(zip(t[t.label == "RNU4-2"].n, t[t.label == "RNU4-2"].base))
    b4a = dict(zip(t[t.label == "RNU4ATAC"].n, t[t.label == "RNU4ATAC"].base))
    pos = h.drop_duplicates("n").sort_values("n")
    assert (np.diff(pos.rnu42_n.to_numpy()) > 0).all()            # no crossing alignment
    assert len(pos) == 38 and (h.groupby("n").size() == 3).all()   # three SNVs per mapped position
    for r in h.itertuples():
        assert bool(r.same_ref) == (b42[int(r.rnu42_n)] == b4a[r.n])


def test_transfer_values_come_from_the_sge_table():
    h = pd.read_csv(os.path.join(INPUTS, "homology_rnu4atac.tsv.gz"), sep="\t").dropna(subset=["h"])
    var = pd.read_csv(os.path.join(INPUTS, "variants.tsv.gz"), sep="\t")
    sge = pd.read_csv(os.path.join(INPUTS, "sge_eval.tsv.gz"), sep="\t")
    s = sge[sge.is_snv].merge(var[["ID", "n", "alt_t"]], on=["ID", "n"])
    alt4a = dict(zip(var.ID, var.alt_t))
    for r in h.itertuples():
        at = s[s.n == int(r.rnu42_n)]
        if r.same_ref:
            assert r.h == pytest.approx(-at[at.alt_t == alt4a[r.ID]].function_score_mean.iloc[0], abs=1e-4)
        else:
            assert r.h == pytest.approx(-at.function_score_mean.mean(), abs=1e-4)
