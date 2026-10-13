"""End-to-end check of evaluate.py -> final_eval.py -> calibrate.py on random predictions and FAKE blind labels
(randomly chosen variant IDs, encrypted with a test passphrase). The real sealed file is never touched.
Run: python -m pytest -q tests/test_final_eval.py   (about 2 min)"""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS = os.path.join(HERE, "model_inputs")
pytestmark = pytest.mark.skipif(not os.path.exists(os.path.join(INPUTS, "variants.tsv.gz")), reason="model_inputs not built")


def fake_run(tmp):
    genes = pd.read_csv(os.path.join(INPUTS, "genes.tsv.gz"), sep="\t")
    var = pd.read_csv(os.path.join(INPUTS, "variants.tsv.gz"), sep="\t")
    runs = tmp / "runs"
    (runs / "preds").mkdir(parents=True)
    json.dump({"has_evo2_emb": False, "has_evo2_delta": False, "device": "cpu", "quick": False, "synthetic_counts": None},
              open(runs / "data_info.json", "w"))
    rng = np.random.default_rng(0)
    scen = {"strict": list(genes.label)}
    for g in ("RNU2-2", "RNU5B-1", "RNU4ATAC", "RNU12", "RNU6ATAC", "RNU6-1", "RNU4-2"):
        c = genes.cluster[genes.gene == g].iloc[0]
        scen[f"lofo_{c}"] = list(genes.label[genes.cluster == c])
    for name, labs in scen.items():
        for seed in (11, 22):
            v = var[var.label.isin(labs)]
            df = pd.DataFrame({"scenario": name, "seed": seed, "ID": v.ID, "label": v.label, "f": rng.gamma(1, 1, len(v)),
                               "gamma": 0.8, "d": rng.normal(size=len(v))})
            df.to_csv(runs / "preds" / f"{name}__seed{seed}.tsv.gz", sep="\t", index=False)
            json.dump({"scenario": name, "seed": seed, "beta": 0.1, "log_r": 1.0, "cfg": {}, "best_epoch": 1,
                       "best_val_ll": -1.0}, open(runs / "preds" / f"{name}__seed{seed}.json", "w"))
    return runs, var


def test_final_evaluation_end_to_end(tmp_path):
    runs, var = fake_run(tmp_path)
    res = tmp_path / "results"
    env = {**os.environ, "SNRNA_RESULTS": str(tmp_path / "r"), "SNRNA_WORK": str(tmp_path / "w")}
    subprocess.run([sys.executable, os.path.join(HERE, "evaluate.py"), "--runs", str(runs), "--results", str(res),
                    "--boot", "20", "--evo2", str(tmp_path / "none")], check=True, capture_output=True, env=env)
    # fake blind labels, sealed with a test passphrase
    rng = np.random.default_rng(1)
    rows = []
    for g in ("RNU2-2", "RNU5B-1"):
        ids = var.ID[(var.gene == g) & ~var.L3_benign_like.astype(bool)].to_numpy()
        for i in rng.choice(ids, 15, replace=False):
            rows.append({"ID": i, "gene": g, "source": "FAKE", "classification": "Pathogenic", "review_status": "test",
                         "inheritance": "dominant", "title": "", "label": "P/LP"})
    sealed = tmp_path / "sealed"
    sealed.mkdir()
    plain = pd.DataFrame(rows).to_csv(sep="\t", index=False).encode()
    env["SNRNA_SEAL_PASS"] = "test-pass"
    subprocess.run(["openssl", "enc", "-aes-256-cbc", "-pbkdf2", "-iter", "200000", "-salt", "-out",
                    str(sealed / "blind_labels.tsv.enc"), "-pass", "env:SNRNA_SEAL_PASS"], input=plain, check=True, env=env)
    (sealed / "SEALED_MANIFEST.jsonl").write_text(json.dumps({"event": "SEALED"}) + "\n")
    cmd = [sys.executable, os.path.join(HERE, "final_eval.py"), "--reason", "test", "--runs", str(runs), "--results", str(res),
           "--sealed", str(sealed), "--boot", "20", "--calib-boot", "20", "--evo2", str(tmp_path / "none")]
    out = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert out.returncode == 0, out.stderr[-2000:]
    rep = (res / "final_report.txt").read_text()
    assert "RNU2-2" in rep and "RNU5B-1" in rep and "framework S" in rep
    assert "FINAL" in (res / "calibration_report.txt").read_text()
    events = [json.loads(x)["event"] for x in (sealed / "SEALED_MANIFEST.jsonl").read_text().splitlines()]
    assert events == ["SEALED", "UNSEALED"]
    # a second call continues with the same decrypted file and does not unseal again
    out2 = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert out2.returncode == 0 and "already unsealed once" in out2.stdout
    assert [json.loads(x)["event"] for x in (sealed / "SEALED_MANIFEST.jsonl").read_text().splitlines()] == events


def test_final_evaluation_refuses_dry_runs(tmp_path):
    runs = tmp_path / "runs"
    runs.mkdir()
    json.dump({"synthetic_counts": 1, "quick": True}, open(runs / "data_info.json", "w"))
    out = subprocess.run([sys.executable, os.path.join(HERE, "final_eval.py"), "--reason", "x", "--runs", str(runs),
                          "--sealed", str(tmp_path)], capture_output=True, text=True)
    assert out.returncode != 0 and "dry run" in out.stderr
