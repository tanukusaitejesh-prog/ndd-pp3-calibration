#!/usr/bin/env python3
"""Stand-in for `bcftools query -R` when bcftools is not installed (self-test only):
keeps the synthetic Roulette records that fall inside work/regions/chrN.regions.tsv."""
import os
import sys

import numpy as np
import pandas as pd

out = sys.argv[1]
work = os.environ.get("WORK", os.path.join(out, "work"))
os.makedirs(os.path.join(work, "roulette"), exist_ok=True)
for c in ("21", "22"):
    reg = pd.read_csv(os.path.join(work, "regions", f"chr{c}.regions.tsv"), sep="\t", header=None, names=["c", "s", "e"])
    r = pd.read_csv(os.path.join(out, "roulette_full", f"chr{c}.tsv.gz"), sep="\t", header=None, dtype=str,
                    keep_default_na=False)
    pos = r[0].astype(int).to_numpy()
    k = np.searchsorted(reg["s"].to_numpy(), pos, side="right") - 1
    inside = (k >= 0) & (pos <= reg["e"].to_numpy()[np.maximum(k, 0)])
    dest = os.path.join(work, "roulette", f"chr{c}.tsv.gz")
    r[inside].to_csv(dest, sep="\t", header=False, index=False)
    open(dest + ".ok", "w").close()
    print(f"chr{c}: kept {int(inside.sum())} of {len(r)} records (python stand-in for bcftools)")
