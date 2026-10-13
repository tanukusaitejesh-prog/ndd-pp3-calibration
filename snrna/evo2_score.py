#!/usr/bin/env python3
"""Evo 2 zero-shot scores for gate G1 (PREREG s6.1). Needs a CUDA GPU with >= 24 GB and the `evo2` package.

    pip install evo2            # see github.com/ArcInstitute/evo2 for CUDA / flash-attn requirements
    python snrna.py evo2-prep   # on any machine: writes work/evo2_inputs.tsv.gz (521 variants, 8,192-nt windows)
    python evo2_score.py work/evo2_inputs.tsv.gz work/evo2_scores.tsv [--model evo2_7b]
    python snrna.py g1          # picks up work/evo2_scores.tsv and applies the G1 rule

For each variant the summed per-token log-likelihood of the reference and alternative windows is computed on the
forward strand and on the reverse complement. snrna.py g1 forms delta = mean over strands of (alt - ref).
The script resumes: variants already in the output file are skipped.
"""
import argparse
import os
import sys

import pandas as pd

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs")
    ap.add_argument("output")
    ap.add_argument("--model", default="evo2_7b")
    args = ap.parse_args()

    import torch
    from evo2 import Evo2

    if not torch.cuda.is_available():
        sys.exit("No CUDA device: Evo 2 7B needs a GPU with >= 24 GB")
    model = Evo2(args.model)
    tok = model.tokenizer

    @torch.inference_mode()
    def loglik(seq):
        ids = torch.tensor(tok.tokenize(seq), dtype=torch.int).unsqueeze(0).to("cuda:0")
        logits, _ = model(ids)
        logp = torch.log_softmax(logits[0, :-1].float(), dim=-1)
        nxt = ids[0, 1:].long()
        return float(logp.gather(1, nxt.unsqueeze(1)).sum()), int(ids.shape[1])

    d = pd.read_csv(args.inputs, sep="\t")
    done = set()
    if os.path.exists(args.output):
        done = set(pd.read_csv(args.output, sep="\t").ID)
    new = not os.path.exists(args.output)
    with open(args.output, "a") as f:
        if new:
            f.write("ID\tll_ref_fwd\tll_alt_fwd\tll_ref_rc\tll_alt_rc\tn_ref\tn_alt\tmodel\n")
        for i, r in enumerate(d.itertuples()):
            if r.ID in done:
                continue
            rf, nr = loglik(r.ref_seq)
            af, na = loglik(r.alt_seq)
            rr, _ = loglik(revcomp(r.ref_seq))
            ar, _ = loglik(revcomp(r.alt_seq))
            f.write(f"{r.ID}\t{rf:.6f}\t{af:.6f}\t{rr:.6f}\t{ar:.6f}\t{nr}\t{na}\t{args.model}\n")
            f.flush()
            if i % 25 == 0:
                print(f"{i + 1}/{len(d)} {r.ID} delta_fwd {af - rf:+.4f}", flush=True)
    print("done:", args.output)


if __name__ == "__main__":
    main()
