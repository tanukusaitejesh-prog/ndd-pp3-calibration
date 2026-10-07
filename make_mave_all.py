"""PP3 x PM1 dependence on functional truth across all human ProteinGym v1.3 single-substitution assays (review 15).
Truth: DMS_score_bin == 0 (damaging). PP3: AlphaMissense >= 0.792 (Supporting) or >= 0.972 (Strong). PM1-like: RMC o/e < 0.4
or CCR >= 95. Mutants are mapped to SNV-reachable column-store sites of the gene's analysed transcript after aligning the assay
sequence to the transcript protein. Within each assay (one gene), R = LR(PP3 & F) / [LR(PP3) LR(F)] with a
position-clustered bootstrap; assays with a PM1-like fraction of 5-95% are informative. Pooled log R: inverse-variance fixed
and DerSimonian-Laird random effects, for genes outside DDG2P monoallelic (primary) and all genes.
Run in the WSL repo: python make_mave_all.py"""
import io, math, os, sys, zipfile
import numpy as np, pandas as pd

sys.path.insert(0, ".")
import gate

idx = pd.read_csv("work_pg.csv")
idx = idx[(idx.taxon == "Human") & (~idx.includes_multiple_mutants.astype(str).str.upper().eq("TRUE") | (idx.DMS_number_single_mutants > 0))]
tx = pd.read_pickle("work/gtf/transcripts.pkl")
dd = pd.read_csv(os.path.join("data", gate.DDG2P_FILE), dtype=str)
q = dd["DDD category"].str.lower().isin(["confirmed", "probable"]) & dd["allelic requirement"].str.lower().str.contains("monoallelic", na=False)
ddg2p = set(dd.loc[q, "gene symbol"].str.upper())
zf = zipfile.ZipFile("work/pg/dms.zip")
names = {n.rsplit("/", 1)[-1]: n for n in zf.namelist()}
rng = np.random.default_rng(20261016)
um = pd.read_csv("work_pg_uniprot_genes.tsv", sep="	").dropna().drop_duplicates("UniProt_ID")     # UniProt REST, 2026-10-07
UNI = dict(zip(um.UniProt_ID, um.gene))


def store(gene):
    t = tx[(tx.gene_name == gene) & (tx.status == "ok")]
    if t.empty:
        return None
    t = t.iloc[0]
    od = f"work/fullstudy/col/chr{t.chrom}"
    if not os.path.exists(od):
        return None
    ix = np.flatnonzero(np.load(f"{od}/tx.npy", mmap_mode="r") == t.tx)
    if len(ix) == 0:
        return None
    s = pd.DataFrame({c: np.load(f"{od}/{c}.npy", mmap_mode="r")[ix] for c in ("prot_pos", "ref_aa", "alt_aa", "am", "rmc_oe", "ccr_pct")})
    s["ref_aa"] = s["ref_aa"].map(chr); s["alt_aa"] = s["alt_aa"].map(chr)
    return s.dropna(subset=["am"]).drop_duplicates(["prot_pos", "alt_aa"])


def lr(m, y):
    return m[y == 1].mean() / max(m[y == 0].mean(), 1e-9)


def Rf(z, cut):
    y, A, B = z.truth.to_numpy(), (z.am >= cut).to_numpy(), z.F.to_numpy()
    a, b, j = lr(A, y), lr(B, y), lr(A & B, y)
    return j / (a * b) if a * b > 0 else np.nan


rows, out = [], []
cache = {}
for r in idx.itertuples():
    gene = UNI.get(r.UniProt_ID, r.DMS_id.split("_HUMAN")[0])
    fn = names.get(f"{r.DMS_id}.csv")
    if fn is None:
        continue
    if gene not in cache:
        cache[gene] = store(gene)
    s = cache[gene]
    if s is None:
        continue
    d = pd.read_csv(io.BytesIO(zf.read(fn)))
    d = d[~d["mutant"].astype(str).str.contains(":")]
    m = d["mutant"].astype(str).str.extract(r"^([A-Z])(\d+)([A-Z])$")
    d = d.assign(wt=m[0], pos=pd.to_numeric(m[1]), mt=m[2]).dropna(subset=["wt", "pos"])
    if len(d) < 50 or "DMS_score_bin" not in d:
        continue
    prot = dict(zip(s.prot_pos, s.ref_aa))
    samp = d.drop_duplicates("pos").sample(min(300, d.pos.nunique()), random_state=1)
    cand = range(-3000, 3001)
    best = max(cand, key=lambda o: sum(prot.get(int(p) + o) == w for p, w in zip(samp.pos, samp.wt)))
    conc = np.mean([prot.get(int(p) + best) == w for p, w in zip(samp.pos, samp.wt) if prot.get(int(p) + best) is not None] or [0])
    if conc < 0.9:
        continue
    d["prot_pos"] = d["pos"].astype(int) + best
    j = d.merge(s, left_on=["prot_pos", "wt", "mt"], right_on=["prot_pos", "ref_aa", "alt_aa"])
    if len(j) < 50:
        continue
    j["truth"] = 1 - j["DMS_score_bin"].astype(int)
    j["F"] = (j.rmc_oe < 0.4) | (j.ccr_pct >= 95)
    fF = j.F.mean()
    rec = dict(ds=r.DMS_id, gene=gene, ddg2p=gene.upper() in ddg2p, n=len(j), damaging=int(j.truth.sum()), F_frac=fF,
               assay=r.coarse_selection_type)
    if 0.05 <= fF <= 0.95 and j.truth.sum() >= 10 and (1 - j.truth).sum() >= 10:
        pos = j.prot_pos.to_numpy(); up = np.unique(pos); gi = pd.Series(np.arange(len(j))).groupby(pos).indices
        for cut, tier in ((0.792, "sup"), (0.972, "str")):
            r0 = Rf(j, cut)
            bs = np.array([Rf(j.iloc[np.concatenate([gi[p] for p in rng.choice(up, len(up))])], cut) for _ in range(500)])
            lb = np.log(bs[np.isfinite(bs) & (bs > 0)])
            rec[f"R_{tier}"] = r0
            rec[f"lo_{tier}"], rec[f"hi_{tier}"] = (np.exp(np.percentile(lb, [2.5, 97.5])) if len(lb) > 50 else (np.nan, np.nan))
            rec[f"se_{tier}"] = lb.std() if len(lb) > 50 else np.nan
    rows.append(rec)
A = pd.DataFrame(rows)
A.to_csv("results_fullstudy/mave_all_sites_summary.tsv", sep="\t", index=False, float_format="%.4g")
inf = A.dropna(subset=["R_sup"])
out.append(f"Human ProteinGym single-substitution assays mapped: {len(A)} ({A.gene.nunique()} genes); informative (PM1-like 5-95%, "
           f">= 10 damaging and >= 10 not): {len(inf)} assays, {inf.gene.nunique()} genes; outside DDG2P: {(~inf.ddg2p).sum()} assays, "
           f"{inf[~inf.ddg2p].gene.nunique()} genes")


def pool(z, tier):
    z = z[np.isfinite(z[f"se_{tier}"]) & (z[f"R_{tier}"] > 0) & (z[f"se_{tier}"] > 0)]
    if len(z) < 2:
        return None
    y, v = np.log(z[f"R_{tier}"].to_numpy()), z[f"se_{tier}"].to_numpy() ** 2
    w = 1 / v
    fe = (w * y).sum() / w.sum(); sfe = math.sqrt(1 / w.sum())
    Q = (w * (y - fe) ** 2).sum(); k = len(y)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
    wr = 1 / (v + tau2)
    re = (wr * y).sum() / wr.sum(); sre = math.sqrt(1 / wr.sum())
    return k, math.exp(fe), math.exp(fe - 1.96 * sfe), math.exp(fe + 1.96 * sfe), math.exp(re), math.exp(re - 1.96 * sre), \
        math.exp(re + 1.96 * sre), (z[f"R_{tier}"] < 1).mean(), tau2


for lab, z in (("outside DDG2P (primary)", inf[~inf.ddg2p]), ("all genes", inf)):
    for tier, nm in (("sup", "Supporting"), ("str", "Strong")):
        p = pool(z, tier)
        if p:
            out.append(f"{lab:24s} AlphaMissense {nm:10s} x PM1-like: assays {p[0]}; fixed {p[1]:.2f} [{p[2]:.2f}, {p[3]:.2f}]; "
                       f"random {p[4]:.2f} [{p[5]:.2f}, {p[6]:.2f}]; share of assays with R < 1 {p[7]:.2f}; tau2 {p[8]:.3f}")
for rr in inf.sort_values("gene").itertuples():
    out.append(f"  {rr.ds:45s} {'DDG2P' if rr.ddg2p else 'other':5s} n {rr.n:5d} F {rr.F_frac:.2f} R_sup {rr.R_sup:.2f} "
               f"[{rr.lo_sup:.2f}, {rr.hi_sup:.2f}] R_str {rr.R_str:.2f}")
open("results_fullstudy/mave_all_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
