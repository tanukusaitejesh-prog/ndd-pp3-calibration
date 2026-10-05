"""Per-residue AlphaFold structure features for the human proteome (amendment A5).

Streams the AlphaFold DB human proteome tar from stdin (curl | python af_structure.py OUT) and never stores it.
For every single-fragment model (F1 only; proteins > 2,700 aa are split into fragments and are left missing):
  plddt     - model confidence (PDB B-factor column)
  n_ca10    - number of other C-alpha atoms within 10 A (burial / packing proxy), |i - j| >= 3
  plddt_w11 - mean pLDDT over the 11-residue window centred on the residue
Output: tab-separated, gzip: acc, pos, aa, plddt, n_ca10, plddt_w11 (one row per residue).
"""
import gzip
import io
import sys
import tarfile
import time

import numpy as np

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def parse_pdb(raw):
    pos, aa, pl, xyz = [], [], [], []
    for ln in raw.splitlines():
        if ln.startswith(b"ATOM") and ln[12:16] == b" CA ":
            pos.append(int(ln[22:26]))
            aa.append(AA3.get(ln[17:20].decode(), "X"))
            pl.append(float(ln[60:66]))
            xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.array(pos), aa, np.array(pl, np.float32), np.array(xyz, np.float32)


def features(pos, pl, xyz):
    n = len(pos)
    cnt = np.zeros(n, np.int32)
    for s in range(0, n, 512):                      # blockwise distances keep memory bounded for long chains
        d2 = ((xyz[s:s + 512, None, :] - xyz[None, :, :]) ** 2).sum(-1)
        near = d2 <= 100.0
        i = np.arange(s, min(s + 512, n))[:, None]
        near &= np.abs(i - np.arange(n)[None, :]) >= 3
        cnt[s:s + 512] = near.sum(1)
    k = np.ones(11, np.float32)
    num = np.convolve(pl, k, "same")
    den = np.convolve(np.ones_like(pl), k, "same")
    return cnt, num / den


def main(out_path):
    t0 = time.time()
    n_prot = n_skip = n_res = 0
    pending = None          # (acc, text, n_res): an F1 model is written only once we know no F2 follows it
    multi = set()
    with gzip.open(out_path, "wt") as out:
        out.write("acc\tpos\taa\tplddt\tn_ca10\tplddt_w11\n")

        def flush():
            nonlocal pending, n_prot, n_res
            if pending is not None and pending[0] not in multi:
                out.write(pending[1])
                n_prot += 1
                n_res += pending[2]
                if n_prot % 2000 == 0:
                    print(f"[{time.time() - t0:7.0f}s] {n_prot:,} proteins, {n_res:,} residues", flush=True)
            pending = None

        src = open(sys.argv[2], "rb") if len(sys.argv) > 2 else sys.stdin.buffer
        tar = tarfile.open(fileobj=src, mode="r|")
        for m in tar:
            if not m.name.endswith(".pdb.gz"):
                continue
            acc, frag = m.name.split("-")[1], m.name.split("-")[2]
            if frag != "F1":
                multi.add(acc)
                n_skip += 1
                continue
            flush()
            raw = gzip.decompress(tar.extractfile(m).read())
            pos, aa, pl, xyz = parse_pdb(raw)
            if len(pos) == 0:
                continue
            cnt, w = features(pos, pl, xyz)
            buf = io.StringIO()
            for p, a, x, c, ww in zip(pos.tolist(), aa, pl.tolist(), cnt.tolist(), w.tolist()):
                buf.write(f"{acc}\t{p}\t{a}\t{x:.2f}\t{c}\t{ww:.2f}\n")
            pending = (acc, buf.getvalue(), len(pos))
        flush()
    print(f"done: {n_prot:,} single-fragment models, {n_res:,} residues; {len(multi):,} multi-fragment proteins "
          f"({n_skip:,} non-F1 files) left missing; "
          f"{time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
