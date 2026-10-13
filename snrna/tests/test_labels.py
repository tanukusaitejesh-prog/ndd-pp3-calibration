"""Label handling: the blind genes are never queried by the non-blind builder, and sealing never exposes the labels.
No network: the ClinVar and RNUdb fetchers are replaced by fakes. Run: python -m pytest -q tests/test_labels.py"""
import hashlib
import json
import os
import subprocess
import sys
import types

import pandas as pd
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import labels as LB  # noqa: E402

pytestmark = pytest.mark.skipif(not os.path.exists(os.path.join(HERE, "model_inputs", "genes.tsv.gz")),
                                reason="model_inputs not built")
SECRET = "chr17-999-A-G"


def fake_fetchers(monkeypatch, calls):
    def clinvar(gene, start, end, chrom):
        calls.append(gene)
        return [{"ID": SECRET if gene in LB.SEALED else f"{chrom}-{start}-A-G", "gene": gene, "source": "ClinVar:VCVx",
                 "classification": "Pathogenic", "review_status": "criteria provided, single submitter", "title": "t"}]

    def rnudb(gene):
        calls.append(gene + ":rnudb")
        return []
    monkeypatch.setattr(LB, "clinvar", clinvar)
    monkeypatch.setattr(LB, "rnudb", rnudb)


def test_nonblind_never_queries_sealed_genes(monkeypatch, tmp_path):
    calls = []
    fake_fetchers(monkeypatch, calls)
    monkeypatch.setattr(LB, "HERE", str(tmp_path))
    LB.cmd_nonblind(types.SimpleNamespace())
    assert calls and not any(c.split(":")[0] in LB.SEALED for c in calls)
    out = pd.read_csv(tmp_path / "labels" / "L2_nonblind.tsv", sep="\t")
    assert not out.gene.isin(LB.SEALED).any() and SECRET not in set(out.ID)


def test_seal_encrypts_without_printing_and_decrypts_with_the_passphrase(monkeypatch, tmp_path, capsys):
    calls = []
    fake_fetchers(monkeypatch, calls)
    monkeypatch.setattr(LB, "HERE", str(tmp_path))
    monkeypatch.setenv("SNRNA_SEAL_PASS", "test-passphrase")
    LB.cmd_seal(types.SimpleNamespace(force=False))
    printed = capsys.readouterr().out
    assert set(c.split(":")[0] for c in calls) == set(LB.SEALED)
    assert SECRET not in printed and "chr17" not in printed
    enc = tmp_path / "sealed" / "blind_labels.tsv.enc"
    assert SECRET.encode() not in enc.read_bytes()
    plain = subprocess.run(["openssl", "enc", "-d", "-aes-256-cbc", "-pbkdf2", "-iter", "200000", "-in", str(enc),
                            "-pass", "env:SNRNA_SEAL_PASS"], capture_output=True, check=True).stdout
    assert SECRET.encode() in plain
    rec = json.loads((tmp_path / "sealed" / "SEALED_MANIFEST.jsonl").read_text().splitlines()[-1])
    assert rec["event"] == "SEALED" and rec["sha256_plaintext"] == hashlib.sha256(plain).hexdigest()
    assert rec["sha256_ciphertext"] == hashlib.sha256(enc.read_bytes()).hexdigest()
    # the existing unseal command reads the same format
    out = subprocess.run([sys.executable, os.path.join(HERE, "snrna.py"), "unseal", str(enc), "--reason", "test",
                          "--pass-env", "SNRNA_SEAL_PASS"], capture_output=True, text=True,
                         env={**os.environ, "SNRNA_RESULTS": str(tmp_path / "r"), "SNRNA_WORK": str(tmp_path / "w")})
    assert out.returncode == 0, out.stderr
    assert (tmp_path / "sealed" / "blind_labels.tsv").read_bytes() == plain
    # sealing happens once
    with pytest.raises(SystemExit):
        LB.cmd_seal(types.SimpleNamespace(force=False))


def test_seal_refuses_without_passphrase(monkeypatch, tmp_path):
    monkeypatch.setattr(LB, "HERE", str(tmp_path))
    monkeypatch.delenv("SNRNA_SEAL_PASS", raising=False)
    with pytest.raises(SystemExit):
        LB.cmd_seal(types.SimpleNamespace(force=False))


def test_evaluator_refuses_label_files_with_sealed_genes(tmp_path):
    import evaluate as EV
    pd.DataFrame({"ID": ["chr1-1-A-G"], "gene": ["RNU2-2"], "label": ["P/LP"]}).to_csv(tmp_path / "L2_nonblind.tsv", sep="\t", index=False)
    with pytest.raises(SystemExit):
        EV.load_l2(str(tmp_path))
