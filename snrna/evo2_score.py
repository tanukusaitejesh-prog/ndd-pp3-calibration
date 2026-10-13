#!/usr/bin/env python3
"""Evo 2 zero-shot scores for gate G1 (PREREG s6.1). Needs a CUDA GPU with 24 GB (e.g. RTX 4090) and the `evo2` package.

    python evo2_score.py evo2_run/evo2_inputs.tsv.gz evo2_run/evo2_scores.tsv --limit 5   # smoke test + timing
    python evo2_score.py evo2_run/evo2_inputs.tsv.gz evo2_run/evo2_scores.tsv             # all 521 variants (resumes)

Each variant has an 8,192-nt GRCh38 reference window centred on it and the same window with the variant applied
(made by `snrna.py evo2-prep`). The script scores four sequences per variant with Evo 2's own `score_sequences`
(sum of per-token log-likelihoods, no BOS token): reference and alternative, forward strand and reverse complement.
`snrna.py g1` then forms delta = mean over strands of (alt - ref), and the damage score is -delta.
A JSON sidecar (<output>.json) records the model, package versions and GPU.
"""
import argparse
import json
import os
import sys
import time

import pandas as pd

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
COLUMNS = ["ID", "ll_ref_fwd", "ll_alt_fwd", "ll_ref_rc", "ll_alt_rc", "n_ref", "n_alt", "model"]


def revcomp(s):
    return s.translate(COMP)[::-1]


def run(inputs, output, score, model_name, limit=None, log=print):
    """Score every variant in `inputs` not already in `output`; score(list_of_seqs) -> list of summed log-likelihoods."""
    d = pd.read_csv(inputs, sep="\t")
    for col in ("ref_seq", "alt_seq"):
        bad = d[~d[col].str.fullmatch(r"[ACGT]+")]
        if len(bad):
            sys.exit(f"{len(bad)} {col} values contain non-ACGT bases (first: {bad.ID.iloc[0]})")
    if (d.ref_seq == d.alt_seq).any():
        sys.exit("some alternative windows equal their reference window: re-run snrna.py evo2-prep")
    done = set(pd.read_csv(output, sep="\t").ID) if os.path.exists(output) else set()
    todo = d[~d.ID.isin(done)]
    if limit:
        todo = todo.head(limit)
    log(f"{len(d)} variants, {len(done)} already scored, scoring {len(todo)} now")
    new = not os.path.exists(output)
    t0 = time.time()
    with open(output, "a") as f:
        if new:
            f.write("\t".join(COLUMNS) + "\n")
        for i, r in enumerate(todo.itertuples(), 1):
            rf, af, rr, ar = score([r.ref_seq, r.alt_seq, revcomp(r.ref_seq), revcomp(r.alt_seq)])
            f.write(f"{r.ID}\t{rf:.6f}\t{af:.6f}\t{rr:.6f}\t{ar:.6f}\t{len(r.ref_seq)}\t{len(r.alt_seq)}\t{model_name}\n")
            f.flush()
            if i <= 3 or i % 25 == 0 or i == len(todo):
                per = (time.time() - t0) / i
                log(f"  {i}/{len(todo)} {r.ID}: delta fwd {af - rf:+.3f}, rc {ar - rr:+.3f} "
                    f"({per:.1f} s/variant, ~{per * (len(todo) - i) / 60:.0f} min left)")
    return len(todo)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs")
    ap.add_argument("output")
    ap.add_argument("--model", default="evo2_7b", help="frozen in PREREG s6.1; change only with a logged deviation")
    ap.add_argument("--limit", type=int, help="score at most this many new variants (smoke test)")
    args = ap.parse_args()

    import torch
    from evo2 import Evo2
    try:
        from importlib.metadata import version
        evo2_version = version("evo2")
    except Exception:  # installed from source without metadata
        from evo2.version import __version__ as evo2_version

    if not torch.cuda.is_available():
        sys.exit("No CUDA device visible. In WSL2, check that `nvidia-smi` works first.")
    gpu = torch.cuda.get_device_properties(0)
    print(f"GPU: {gpu.name}, {gpu.total_memory / 2**30:.1f} GB; torch {torch.__version__}; evo2 {evo2_version}")
    model = Evo2(args.model)

    def score(seqs):
        # batch_size=1: windows differ in length for indels, and one 8-kb sequence at a time fits in 24 GB
        return model.score_sequences(seqs, batch_size=1, prepend_bos=False, reduce_method="sum")

    n = run(args.inputs, args.output, score, args.model, args.limit)
    meta = {"model": args.model, "evo2": evo2_version, "torch": torch.__version__, "cuda": torch.version.cuda,
            "gpu": gpu.name, "reduce_method": "sum", "prepend_bos": False, "use_kernels": False,
            "inputs": os.path.basename(args.inputs), "scored_this_run": n,
            "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    with open(args.output + ".json", "w") as f:
        json.dump(meta, f, indent=1)
    print(f"done: {args.output}  (peak GPU memory {torch.cuda.max_memory_allocated() / 2**30:.1f} GB)")


if __name__ == "__main__":
    main()
