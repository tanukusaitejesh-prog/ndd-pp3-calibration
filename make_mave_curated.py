"""Review 20: PP3 x curated PM1 on functional truth (descriptive; not pre-specified).
Curated PM1 residues: protein positions at which any ClinGen expert panel applied PM1 (any strength) to a missense
variant in the Evidence Repository (export 2026-10-07); a gene is used only if >= 90% of these residues match the reference
amino acid of the analysed transcript. Assays: human ProteinGym v1.3 single-substitution assays, mapped as in
make_mave_all.py. Truth: DMS_score_bin == 0. PP3: AlphaMissense >= 0.792 (Supporting) or >= 0.972 (Strong).
F = curated PM1 residue. Informative assays: >= 10 assayed variants at curated residues, >= 10 damaging and >= 10 not.
R within assay with a position-clustered bootstrap; pooled by DerSimonian-Laird random effects.
Run in the WSL repo: python make_mave_curated.py"""
import io, math, os, re, sys, zipfile
import numpy as np, pandas as pd

sys.path.insert(0, ".")
AA3 = dict(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G", His="H", Ile="I", Leu="L", Lys="K",
           Met="M", Phe="F", Pro="P", Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V")
OUT = "results_fullstudy/mave_curated_summary.txt"
er = pd.read_csv("work/erepo/erepo_2026-10-07.tsv", sep="\t", dtype=str)
met = er["Applied Evidence Codes (Met)"].fillna("")
er = er[met.str.contains(r"\bPM1(?:_|\b)", regex=True)]
pm = er["#Variation"].str.extract(r"\(([A-Za-z0-9-]+)\):c\.[^ ]+ \(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
pm.columns = ["gene", "ref3", "pos", "alt3"]
pm = pm.dropna()
pm["ref"] = pm.ref3.map(AA3); pm["pos"] = pm.pos.astype(int)
pm = pm.dropna(subset=["ref"]).drop_duplicates(["gene", "pos"])
out = [f"Curated PM1: {len(er)} Evidence Repository classifications with PM1 met; {len(pm)} distinct missense residues in "
       f"{pm.gene.nunique()} genes"]
for g_ in ("TP53", "PTEN", "BRCA1"):
    out.append(f"  {g_} curated PM1 residues: {sorted(pm.loc[pm.gene == g_, 'pos'].tolist())}")

idx = pd.read_csv("work_pg.csv")
idx = idx[(idx.taxon == "Human") & (~idx.includes_multiple_mutants.astype(str).str.upper().eq("TRUE") | (idx.DMS_number_single_mutants > 0))]
tx = pd.read_pickle("work/gtf/transcripts.pkl")
zf = zipfile.ZipFile("work/pg/dms.zip")
names = {n.rsplit("/", 1)[-1]: n for n in zf.namelist()}
um = pd.read_csv("work_pg_uniprot_genes.tsv", sep="\t").dropna().drop_duplicates("UniProt_ID")
UNI = dict(zip(um.UniProt_ID, um.gene))
rng = np.random.default_rng(20261021)


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


cur_genes = set(pm.gene)
rows, cache, okgene = [], {}, {}
for r in idx.itertuples():
    gene = UNI.get(r.UniProt_ID, r.DMS_id.split("_HUMAN")[0])
    if gene not in cur_genes:
        continue
    fn = names.get(f"{r.DMS_id}.csv")
    if fn is None:
        continue
    if gene not in cache:
        cache[gene] = store(gene)
    s = cache[gene]
    if s is None:
        continue
    prot = dict(zip(s.prot_pos, s.ref_aa))
    cg = pm[pm.gene == gene]
    if gene not in okgene:
        conc_c = np.mean([prot.get(p) == a for p, a in zip(cg.pos, cg.ref)])
        okgene[gene] = conc_c
        out.append(f"  {gene}: {len(cg)} curated PM1 residues; reference concordance {conc_c:.2f}")
    if okgene[gene] < 0.9:
        continue
    d = pd.read_csv(io.BytesIO(zf.read(fn)))
    d = d[~d["mutant"].astype(str).str.contains(":")]
    m = d["mutant"].astype(str).str.extract(r"^([A-Z])(\d+)([A-Z])$")
    d = d.assign(wt=m[0], pos=pd.to_numeric(m[1]), mt=m[2]).dropna(subset=["wt", "pos"])
    if len(d) < 50 or "DMS_score_bin" not in d:
        continue
    samp = d.drop_duplicates("pos").sample(min(300, d.pos.nunique()), random_state=1)
    best = max(range(-3000, 3001), key=lambda o: sum(prot.get(int(p) + o) == w for p, w in zip(samp.pos, samp.wt)))
    conc = np.mean([prot.get(int(p) + best) == w for p, w in zip(samp.pos, samp.wt) if prot.get(int(p) + best) is not None] or [0])
    if conc < 0.9:
        continue
    d["prot_pos"] = d["pos"].astype(int) + best
    j = d.merge(s, left_on=["prot_pos", "wt", "mt"], right_on=["prot_pos", "ref_aa", "alt_aa"])
    if len(j) < 50:
        continue
    j["truth"] = 1 - j["DMS_score_bin"].astype(int)
    j["F"] = j.prot_pos.isin(set(cg.pos))
    j["Fp"] = (j.rmc_oe < 0.4) | (j.ccr_pct >= 95)
    y = j.truth.to_numpy()
    rec = dict(ds=r.DMS_id, gene=gene, n=len(j), nF=int(j.F.sum()), damaging=int(y.sum()), F_frac=j.F.mean(),
               Fp_frac=j.Fp.mean(), LR_F=lr(j.F.to_numpy(), y), LR_Fp=lr(j.Fp.to_numpy(), y),
               LR_PP3=lr((j.am >= 0.792).to_numpy(), y))
    if j.F.sum() >= 10 and y.sum() >= 10 and (1 - y).sum() >= 10:
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
A.to_csv("results_fullstudy/mave_curated_sites.tsv", sep="\t", index=False, float_format="%.4g")
inf = A.dropna(subset=["R_sup"]) if "R_sup" in A else A.iloc[0:0]
out.append(f"Assays in genes with curated PM1: {len(A)} ({A.gene.nunique()} genes); informative: {len(inf)} assays, {inf.gene.nunique()} genes")


def pool(z, tier):
    z = z[np.isfinite(z[f"se_{tier}"]) & (z[f"R_{tier}"] > 0) & (z[f"se_{tier}"] > 0)]
    if len(z) < 2:
        return None
    yv, v = np.log(z[f"R_{tier}"].to_numpy()), z[f"se_{tier}"].to_numpy() ** 2
    w = 1 / v
    fe = (w * yv).sum() / w.sum()
    Q = (w * (yv - fe) ** 2).sum(); k = len(yv)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
    wr = 1 / (v + tau2)
    re_ = (wr * yv).sum() / wr.sum(); sre = math.sqrt(1 / wr.sum())
    return k, z.gene.nunique(), math.exp(re_), math.exp(re_ - 1.96 * sre), math.exp(re_ + 1.96 * sre), (z[f"R_{tier}"] < 1).mean()


for tier, nm in (("sup", "Supporting"), ("str", "Strong")):
    p = pool(inf, tier) if len(inf) else None
    if p:
        out.append(f"AlphaMissense {nm:10s} x curated PM1: assays {p[0]} ({p[1]} genes); random effects {p[2]:.2f} [{p[3]:.2f}, {p[4]:.2f}]; "
                   f"share of assays with R < 1 {p[5]:.2f}")
if len(inf):
    out.append(f"Informative assays: median LR(curated PM1) {inf.LR_F.median():.2f}; median LR(PM1-like proxy) {inf.LR_Fp.median():.2f}; "
               f"median LR(PP3 Supporting) {inf.LR_PP3.median():.2f}")
for rr in A.sort_values("gene").itertuples():
    out.append(f"  {rr.ds:45s} n {rr.n:5d} curated F {rr.nF:4d} ({rr.F_frac:.2f}) proxy F {rr.Fp_frac:.2f} LR_F {rr.LR_F:.2f} "
               f"LR_proxy {rr.LR_Fp:.2f}" + (f" R_sup {rr.R_sup:.2f} [{rr.lo_sup:.2f}, {rr.hi_sup:.2f}] R_str {rr.R_str:.2f}"
                                             if hasattr(rr, "R_sup") and pd.notna(rr.R_sup) else ""))
open(OUT, "w").write("\n".join(out) + "\n")
print("\n".join(out))
