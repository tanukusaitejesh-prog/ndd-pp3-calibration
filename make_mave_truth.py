"""Task 1 (exploratory): PP3 x PM1 overlap and conditional PP3 points on MAVE truth (ProteinGym v1.3 DMS substitutions)
in DDG2P genes. Truth = DMS_score_bin == 0 (functionally damaging); single substitutions only. Each DMS mutant is mapped
to SNV-reachable column-store sites of the gene's analysed transcript by (prot_pos + offset, ref_aa, alt_aa); the offset
aligns the DMS target sequence to the transcript protein. Position-clustered bootstrap, 1,000, seed 20261009.
Run in the WSL repo: python make_mave_truth.py  (needs work/pg/dms.zip and work_pg.csv)"""
import io, math, zipfile
import numpy as np, pandas as pd

DATASETS = ["DLG4_HUMAN_Faure_2021", "PITX2_HUMAN_Tsuboyama_2023_2L7M", "PPM1D_HUMAN_Miller_2022",
            "PTEN_HUMAN_Matreyek_2021", "PTEN_HUMAN_Mighell_2018", "RAF1_HUMAN_Zinkus-Boltz_2019",
            "SHOC2_HUMAN_Kwon_2022", "YAP1_HUMAN_Araya_2012"]
tx = pd.read_pickle("work/gtf/transcripts.pkl")
idx = pd.read_csv("work_pg.csv").set_index("DMS_id")
zf = zipfile.ZipFile("work/pg/dms.zip")
names = {n.rsplit("/", 1)[-1]: n for n in zf.namelist()}


def store(gene):
    t = tx[(tx.gene_name == gene) & (tx.status == "ok")].iloc[0]
    od = f"work/fullstudy/col/chr{t.chrom}"
    ix = np.flatnonzero(np.load(f"{od}/tx.npy", mmap_mode="r") == t.tx)
    s = pd.DataFrame({c: np.load(f"{od}/{c}.npy", mmap_mode="r")[ix]
                      for c in ("prot_pos", "ref_aa", "alt_aa", "am", "rmc_oe", "ccr_pct")})
    s["ref_aa"] = s["ref_aa"].map(chr); s["alt_aa"] = s["alt_aa"].map(chr)
    return s.dropna(subset=["am"]).drop_duplicates(["prot_pos", "alt_aa"])


rows, out = [], []
for ds in DATASETS:
    gene = ds.split("_")[0]
    d = pd.read_csv(io.BytesIO(zf.read(names[f"{ds}.csv"])))
    d = d[~d["mutant"].str.contains(":")]
    m = d["mutant"].str.extract(r"^([A-Z])(\d+)([A-Z])$")
    d = d.assign(wt=m[0], pos=m[1].astype(int), mt=m[2]).dropna(subset=["wt"])
    s = store(gene)
    prot = dict(zip(s.prot_pos, s.ref_aa))
    best = max(range(-2000, 2001), key=lambda o: sum(prot.get(p + o) == w for p, w in zip(d.pos, d.wt)))
    wt_ok = np.mean([prot.get(p + best) == w for p, w in zip(d.pos, d.wt)])
    d["prot_pos"] = d["pos"] + best
    j = d.merge(s, left_on=["prot_pos", "wt", "mt"], right_on=["prot_pos", "ref_aa", "alt_aa"], how="inner")
    j = j.assign(ds=ds, gene=gene, truth=1 - j["DMS_score_bin"].astype(int))
    out.append(f"{ds}: {len(d)} single substitutions, offset {best:+d}, WT concordance {wt_ok:.3f}, "
               f"SNV-reachable & mapped {len(j)} (damaging {int(j.truth.sum())}); binarisation "
               f"{idx.loc[ds, 'DMS_binarization_method']}, assay {idx.loc[ds, 'coarse_selection_type']}")
    rows.append(j[["ds", "gene", "prot_pos", "wt", "mt", "DMS_score", "truth", "am", "rmc_oe", "ccr_pct"]])
x = pd.concat(rows, ignore_index=True)
x["F"] = (x["rmc_oe"] < 0.4) | (x["ccr_pct"] >= 95)
x.to_csv("results_fullstudy/mave_truth_sites.tsv", sep="\t", index=False)
rng = np.random.default_rng(20261009)


def lrs(z, cut):
    def lr(mk):
        return mk[z.truth == 1].mean() / max(mk[z.truth == 0].mean(), 1e-9)
    pp, f = z.am >= cut, z.F
    return lr(pp), lr(f), lr(pp & f), lr(pp & ~f)


def report(t, lab, boot=True):
    for cut, tier in ((0.792, "Supporting"), (0.972, "Strong")):
        lp, lf, lj, lnf = lrs(t, cut)
        R = lj / (lp * lf) if lp * lf > 0 else np.nan
        ci = ""
        if boot:
            cl = (t.ds + ":" + t.prot_pos.astype(str)).to_numpy()
            u = np.unique(cl); gi = pd.Series(np.arange(len(t))).groupby(cl).indices
            bs = []
            for _ in range(1000):
                b = lrs(t.iloc[np.concatenate([gi[g] for g in rng.choice(u, len(u))])], cut)
                bs.append(b[2] / (b[0] * b[1]) if b[0] * b[1] > 0 else np.nan)
            lo, hi = np.nanpercentile(bs, [2.5, 97.5]); ci = f" [{lo:.2f}, {hi:.2f}]"
        pts = math.log(max(lj / lf, 1e-9)) / math.log(2.406) if lf > 0 else np.nan
        out.append(f"{lab:34s} AM>={tier:10s} n {len(t)} LR(PP3) {lp:.2f} LR(F) {lf:.2f} joint {lj:.2f} R {R:.2f}{ci}; "
                   f"PP3 given F {lj / lf if lf else np.nan:.2f} ({pts:.2f} points) vs not F {lnf:.2f}; "
                   f"n(PP3&F) {int(((t.am >= cut) & t.F).sum())}, F frac {t.F.mean():.2f}")


report(x, "POOLED (all datasets)")
report(x[x.ds != "PTEN_HUMAN_Matreyek_2021"], "POOLED (PTEN activity only)")
for ds, g in x.groupby("ds"):
    if g.F.any() and (~g.F).any():
        report(g, ds, boot=False)
    else:
        out.append(f"{ds:34s} F is {'all' if g.F.all() else 'never'} true; R undefined")
open("results_fullstudy/mave_truth_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
