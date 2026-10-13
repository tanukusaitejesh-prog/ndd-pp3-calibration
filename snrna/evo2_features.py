#!/usr/bin/env python3
"""Evo 2 features for the selection model (amendment A2), run on the A100. Resumable; safe to re-run.

For every analysis-set gene: one 8,192-nt GRCh38 window centred on the gene (model_inputs/windows.tsv.gz).
  * reference pass, forward strand and reverse complement: summed log-likelihood, plus the layer-28 embedding
    (blocks.28.mlp.l3) at each transcript position; the forward and reverse-complement embeddings are averaged
    (the model is causal, so the two strands give left and right context)
  * every possible SNV in the gene: the same window with the variant applied, forward and reverse complement,
    summed log-likelihood (batched; same length as the reference)

Outputs (evo2_run/features/):
  deltas.tsv               ID, ll_ref_fwd, ll_alt_fwd, ll_ref_rc, ll_alt_rc     (delta = mean over strands of alt - ref)
  emb/<label>.npy          float16 array (L, 4096): embedding of transcript positions n = 1..L
  run_info.json            versions, GPU, timing

    python evo2_features.py                       # all genes, all variants
    python evo2_features.py --limit-genes 1       # smoke test
    SNRNA_EVO2_MOCK=1 python evo2_features.py     # CPU test with a fake Evo 2 (development only)
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
LAYER = "blocks.28.mlp.l3"


def revcomp(s):
    return s.translate(COMP)[::-1]


def load_model(name):
    import torch
    if os.environ.get("SNRNA_EVO2_MOCK"):
        sys.path.insert(0, os.path.join(HERE, "tests"))
        from mock_evo2 import Evo2  # noqa: deterministic fake with the Evo 2 v0.6.0 call signatures
        return Evo2(name), "cpu", "mock"
    if not torch.cuda.is_available():
        sys.exit("No CUDA device visible (set SNRNA_EVO2_MOCK=1 only for development tests)")
    from evo2 import Evo2
    try:
        from importlib.metadata import version
        v = version("evo2")
    except Exception:
        from evo2.version import __version__ as v
    return Evo2(name), "cuda:0", v


class Scorer:
    def __init__(self, model, device, batch):
        import torch
        self.torch, self.model, self.device, self.batch = torch, model, device, batch
        self.tok = model.tokenizer

    def ids(self, seqs):
        t = self.torch
        return t.tensor([self.tok.tokenize(s) for s in seqs], dtype=t.long).to(self.device)

    def _ll(self, logits, ids):
        t = self.torch
        lp = t.log_softmax(logits[:, :-1].float(), dim=-1)
        return lp.gather(2, ids[:, 1:].unsqueeze(-1)).squeeze(-1).sum(1)

    def ll_and_emb(self, seq):
        """Summed log-likelihood and layer embedding (L, D) for one sequence."""
        t = self.torch
        with t.inference_mode():
            ids = self.ids([seq])
            out, emb = self.model(ids, return_embeddings=True, layer_names=[LAYER])
            logits = out[0] if isinstance(out, (tuple, list)) else out
            ll = float(self._ll(logits, ids)[0])
            e = emb[LAYER][0].float().cpu().numpy()
        return ll, e

    def ll_batch(self, seqs):
        """Summed log-likelihoods for equal-length sequences. Every forward pass has the same shape (the last chunk is
        padded with copies) and the whole list is scored at one batch size: on OOM the batch is halved and the list
        restarted, so reference and alternative scores of a gene always come from identically shaped passes."""
        t = self.torch
        while True:
            out, i, b = [], 0, self.batch
            try:
                while i < len(seqs):
                    chunk = seqs[i:i + b]
                    padded = chunk + [chunk[-1]] * (b - len(chunk))
                    with t.inference_mode():
                        ids = self.ids(padded)
                        res = self.model.model(ids)
                        logits = res[0] if isinstance(res, (tuple, list)) else res
                        if isinstance(logits, (tuple, list)):
                            logits = logits[0]
                        out.extend(self._ll(logits, ids).tolist()[:len(chunk)])
                    i += len(chunk)
                return out
            except t.cuda.OutOfMemoryError:
                if b == 1:
                    raise
                t.cuda.empty_cache()
                self.batch = max(1, b // 2)
                print(f"  CUDA out of memory: batch size now {self.batch}; rescoring this gene", flush=True)


def preflight(sc, seq):
    """Check the embedding layer exists and has the expected shape, and measure numerical noise of batched scoring
    (identical sequences in one batch, and batch of copies vs single pass). Returns a dict for run_info.json."""
    short = seq[:512]
    ll, e = sc.ll_and_emb(short)
    if e.ndim != 2 or e.shape[0] != len(short) or not np.isfinite(e).all():
        sys.exit(f"preflight: layer {LAYER} returned shape {e.shape} (expected ({len(short)}, D)) or non-finite values")
    copies = sc.ll_batch([seq] * sc.batch)
    single, _ = sc.ll_and_emb(seq)
    res = {"embedding_dim": int(e.shape[1]), "batch": sc.batch, "spread_within_batch": float(np.ptp(copies)),
           "batched_minus_single": float(copies[0] - single), "ll_ref_window": float(copies[0])}
    print("preflight:", json.dumps(res), flush=True)
    if res["spread_within_batch"] > 0.05:
        sys.exit("preflight: identical sequences score differently within a batch (> 0.05 nats); use --batch 1")
    return res


def transcript_window_index(g, n, W):
    """Index in the forward window of transcript position n (1-based) of gene row g."""
    pos = g.start + n - 1 if g.strand == "+" else g.end - n + 1
    return pos - g.win_start


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--out", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--model", default="evo2_7b")
    ap.add_argument("--batch", type=int, default=4, help="variant windows per forward pass (halves automatically on OOM)")
    ap.add_argument("--limit-genes", type=int, default=0)
    args = ap.parse_args()

    genes = pd.read_csv(os.path.join(args.inputs, "genes.tsv.gz"), sep="\t")
    wins = pd.read_csv(os.path.join(args.inputs, "windows.tsv.gz"), sep="\t").set_index("label")
    var = pd.read_csv(os.path.join(args.inputs, "variants.tsv.gz"), sep="\t")
    os.makedirs(os.path.join(args.out, "emb"), exist_ok=True)
    dpath = os.path.join(args.out, "deltas.tsv")
    done = set(pd.read_csv(dpath, sep="\t").ID) if os.path.exists(dpath) else set()
    if args.limit_genes:
        genes = genes.head(args.limit_genes)

    model, device, ver = load_model(args.model)
    sc = Scorer(model, device, args.batch)
    W = len(wins.iloc[0]["seq"])
    pre = preflight(sc, wins.loc[genes.label.iloc[0], "seq"].upper())
    t0, n_scored = time.time(), 0
    todo_total = int((~var.ID.isin(done) & var.label.isin(genes.label)).sum())
    print(f"Evo 2 features: {len(genes)} genes, {todo_total} variants to score (batch {args.batch}), device {device}", flush=True)
    new = not os.path.exists(dpath)
    with open(dpath, "a") as fout:
        if new:
            fout.write("ID\tll_ref_fwd\tll_alt_fwd\tll_ref_rc\tll_alt_rc\n")
        for g in genes.itertuples():
            w = wins.loc[g.label]
            seq = w["seq"].upper()
            assert w["win_start"] == g.win_start and len(seq) == W
            embp = os.path.join(args.out, "emb", f"{g.label.replace('|', '_')}.npy")
            vg = var[(var.label == g.label) & ~var.ID.isin(done)]
            if os.path.exists(embp) and not len(vg):
                continue
            if not os.path.exists(embp):
                _, emb_f = sc.ll_and_emb(seq)
                _, emb_r = sc.ll_and_emb(revcomp(seq))
                idx = np.array([transcript_window_index(g, n, W) for n in range(1, g.length + 1)])
                assert idx.min() >= 0 and idx.max() < W
                e = (emb_f[idx] + emb_r[W - 1 - idx]) / 2
                np.save(embp, e.astype(np.float16))
            if not len(vg):
                continue
            alts = []
            for r in vg.itertuples():
                off = r.pos - g.win_start
                assert seq[off] == r.ref, (r.ID, seq[off])
                alts.append(seq[:off] + r.alt + seq[off + 1:])
            # reference scored in the same batched pass shape as its alternatives (first element of each list)
            lf = sc.ll_batch([seq] + alts)
            lr = sc.ll_batch([revcomp(seq)] + [revcomp(s) for s in alts])
            ll_f, ll_r = lf[0], lr[0]
            for r, a, b in zip(vg.itertuples(), lf[1:], lr[1:]):
                fout.write(f"{r.ID}\t{ll_f:.6f}\t{a:.6f}\t{ll_r:.6f}\t{b:.6f}\n")
            fout.flush()
            n_scored += len(vg)
            el = time.time() - t0
            rate = n_scored / el if el > 0 else 0
            print(f"  {g.label}: {len(vg)} variants, ref LL fwd {ll_f:.1f} rc {ll_r:.1f}; {n_scored}/{todo_total} done, "
                  f"{rate * 60:.0f}/min, ~{(todo_total - n_scored) / max(rate, 1e-9) / 60:.0f} min left", flush=True)
    info = {"model": args.model, "evo2": ver, "layer": LAYER, "window": W, "batch_final": sc.batch, "device": device,
            "variants_scored_this_run": n_scored, "preflight": pre, "seconds": round(time.time() - t0, 1),
            "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    try:
        import torch
        if torch.cuda.is_available():
            info.update({"gpu": torch.cuda.get_device_name(0), "torch": torch.__version__,
                         "peak_gpu_gb": round(torch.cuda.max_memory_allocated() / 2**30, 1)})
    except Exception:
        pass
    with open(os.path.join(args.out, "run_info.json"), "w") as f:
        json.dump(info, f, indent=1)
    print("Evo 2 features done:", json.dumps(info), flush=True)


if __name__ == "__main__":
    main()
