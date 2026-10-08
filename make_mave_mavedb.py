"""A21 (pre-registered, prereg_fullstudy/A21_MAVEDB_CURATED_PM1_PREREG.md): PP3 x curated PM1 on MaveDB functional truth.
Gene as unit; random-residue negative control. Inputs: work/erepo/curated_pm1_residues.tsv, work/mavedb/candidates.tsv,
work/am/AlphaMissense_aa_substitutions.tsv.gz; downloads score sets and UniProt sequences."""
import gzip
import io
import json
import math
import re
import urllib.request

import numpy as np
import pandas as pd

rng = np.random.default_rng(20261008)
AA3 = dict(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G", His="H", Ile="I", Leu="L", Lys="K",
           Met="M", Phe="F", Pro="P", Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V")
SUP, STR = 0.792, 0.972
pm = pd.read_csv("work/erepo/curated_pm1_residues.tsv", sep="\t")
cand = pd.read_csv("work/mavedb/candidates.tsv", sep="\t")


def get(url, js=True):
    with urllib.request.urlopen(url, timeout=120) as r:
        b = r.read()
    return json.loads(b) if js else b.decode()


uni_seq, meta = {}, {}
for urn in cand.urn:
    try:
        m = get(f"https://api.mavedb.org/api/v1/score-sets/{urn}")
    except Exception as e:
        print("meta fail", urn, e); continue
    tg = m["targetGenes"]
    if len(tg) != 1:
        continue
    uid = (tg[0].get("uniprotIdFromMappedMetadata") or "")
    org = json.dumps(tg[0])
    meta[urn] = dict(uniprot=uid, human=("Homo sapiens" in org or "9606" in org or bool(uid)))
uids = sorted({v["uniprot"] for v in meta.values() if v["uniprot"]})
for u in uids:
    try:
        uni_seq[u] = "".join(get(f"https://rest.uniprot.org/uniprotkb/{u}.fasta", js=False).split("\n")[1:])
    except Exception as e:
        print("uniprot fail", u, e)
# AlphaMissense for these UniProt IDs
am = {}
with gzip.open("work/am/AlphaMissense_aa_substitutions.tsv.gz", "rt") as fh:
    for line in fh:
        if line[0] == "#" or line.startswith("uniprot_id"):
            continue
        u, v, sc, _ = line.rstrip("\n").split("\t")
        if u in uni_seq:
            am[(u, v)] = float(sc)
print(f"score sets with metadata {len(meta)}; UniProt sequences {len(uni_seq)}; AlphaMissense entries {len(am)}")

PRO = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=)$")


def lr(m, y):
    return max(m[y == 1].mean(), 1e-9) / max(m[y == 0].mean(), 1e-9)


def Rv(d, cut, F):
    y, A = d.truth.to_numpy(), (d.am >= cut).to_numpy()
    return lr(A & F, y) / (lr(A, y) * lr(F, y))


rows, ctrl_rows, setinfo = [], [], []
for g, sets in cand.groupby("gene"):
    res = set(pm.loc[pm.gene == g, "pos"])
    refs = dict(zip(pm.loc[pm.gene == g, "pos"], pm.loc[pm.gene == g, "ref"]))
    for urn in sets.urn:
        if urn not in meta or not meta[urn]["uniprot"] or meta[urn]["uniprot"] not in uni_seq:
            continue
        u = meta[urn]["uniprot"]; seq = uni_seq[u]
        try:
            sc = pd.read_csv(io.StringIO(get(f"https://api.mavedb.org/api/v1/score-sets/{urn}/scores", js=False)))
        except Exception as e:
            print("scores fail", urn, e); continue
        if "hgvs_pro" not in sc or "score" not in sc:
            continue
        p = sc.hgvs_pro.astype(str).str.extract(PRO)
        p.columns = ["ref3", "pos", "alt3"]
        d = pd.concat([sc[["score"]], p], axis=1).dropna(subset=["ref3", "pos", "score"])
        d["pos"] = d.pos.astype(int)
        d["ref"] = d.ref3.map(AA3)
        d = d.dropna(subset=["ref"])
        best, bo = 0, 0
        for off in range(-50, 2001):
            pp = d.pos + off
            ok = (pp >= 1) & (pp <= len(seq))
            if ok.sum() < len(d) * 0.5:
                continue
            match = (d.ref[ok].to_numpy() == np.array([seq[i - 1] for i in pp[ok]])).mean()
            if match > best:
                best, bo = match, off
            if best > 0.999:
                break
        cur_ok = np.mean([seq[r - 1] == refs[r] for r in res if r <= len(seq)]) if res else 0
        info = dict(gene=g, urn=urn, uniprot=u, offset=bo, ref_match=best, curated_ref_match=cur_ok)
        if best < 0.9 or cur_ok < 0.9:
            setinfo.append({**info, "status": "mapping failed"}); continue
        d["upos"] = d.pos + bo
        ter = d[d.alt3 == "Ter"].score; syn = d[d.alt3 == "="].score
        if len(ter) < 10 or len(syn) < 10:
            setinfo.append({**info, "status": "too few nonsense/synonymous"}); continue
        mt, ms = ter.median(), syn.median()
        mis = d[d.alt3.isin(AA3)].copy()
        mis["alt"] = mis.alt3.map(AA3)
        mis["truth"] = (abs(mis.score - mt) < abs(mis.score - ms)).astype(int)
        mis["am"] = [am.get((u, f"{r}{p_}{a}")) for r, p_, a in zip(mis.ref, mis.upos, mis.alt)]
        mis = mis.dropna(subset=["am"])
        F = mis.upos.isin(res).to_numpy()
        ok = (F.sum() >= 10) and (mis.truth.sum() >= 10) and ((1 - mis.truth).sum() >= 10)
        setinfo.append({**info, "status": "informative" if ok else "not informative", "n": len(mis), "n_F": int(F.sum()),
                        "damaging": int(mis.truth.sum())})
        if not ok:
            continue
        pos_u = mis.upos.unique()
        gi = mis.groupby("upos").indices
        est = {}
        for nm, cut in [("sup", SUP), ("str", STR)]:
            r0 = Rv(mis, cut, F)
            bs = []
            for _ in range(500):
                pick = rng.choice(pos_u, len(pos_u))
                b = mis.iloc[np.concatenate([gi[x] for x in pick])]
                bs.append(math.log(max(Rv(b, cut, b.upos.isin(res).to_numpy()), 1e-6)))
            est[nm] = (r0, np.std(bs))
        A = (mis.am >= SUP).to_numpy(); y = mis.truth.to_numpy()
        share = math.log(max(lr(A & F, y) / lr(F, y), 1e-9)) / math.log(lr(A, y))
        rows.append(dict(gene=g, urn=urn, n=len(mis), n_F=int(F.sum()), LR_F=lr(F, y), LR_PP3=lr(A, y), R_sup=est["sup"][0],
                         se_sup=est["sup"][1], R_str=est["str"][0], se_str=est["str"][1], share_kept=share))
        # negative control: random residue sets of the same size among non-curated assayed residues
        noncur = np.setdiff1d(pos_u, list(res))
        k = len(set(mis.upos[F]))
        if len(noncur) > k:
            cr = []
            for _ in range(200):
                fake = set(rng.choice(noncur, k, replace=False))
                cr.append(Rv(mis, SUP, mis.upos.isin(fake).to_numpy()))
            ctrl_rows.append(dict(gene=g, urn=urn, R_ctrl_median=float(np.median(cr)), R_ctrl_lo=float(np.quantile(cr, .025)), R_ctrl_hi=float(np.quantile(cr, .975))))
S = pd.DataFrame(setinfo); A_ = pd.DataFrame(rows); C = pd.DataFrame(ctrl_rows)
S.to_csv("results_fullstudy/a21_mavedb_sets.tsv", sep="\t", index=False)
A_.to_csv("results_fullstudy/a21_mavedb_estimates.tsv", sep="\t", index=False)
C.to_csv("results_fullstudy/a21_mavedb_negative_control.tsv", sep="\t", index=False)


def dl(y, se):
    w = 1 / se ** 2
    mu = (w * y).sum() / w.sum()
    q = (w * (y - mu) ** 2).sum()
    k = len(y)
    tau2 = max(0.0, (q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0.0
    w2 = 1 / (se ** 2 + tau2)
    m = (w2 * y).sum() / w2.sum()
    s = math.sqrt(1 / w2.sum())
    return m, s


def by_gene(df, col, secol):
    out = []
    for g, d in df.groupby("gene"):
        y = np.log(d[col].clip(lower=1e-6)); se = d[secol].clip(lower=0.05)
        w = 1 / se ** 2
        out.append((g, float((w * y).sum() / w.sum()), float(math.sqrt(1 / w.sum()))))
    return pd.DataFrame(out, columns=["gene", "y", "se"])


out = [f"A21: score sets examined {len(S)}; status: " + ", ".join(f"{k} {v}" for k, v in S.status.value_counts().items())]
if len(A_):
    for band, col, secol in [("Supporting", "R_sup", "se_sup"), ("Strong", "R_str", "se_str")]:
        G = by_gene(A_, col, secol)
        m, s = dl(G.y.to_numpy(), G.se.to_numpy())
        out.append(f"AlphaMissense {band} x curated PM1: genes {len(G)}, score sets {len(A_)}; random effects R {math.exp(m):.2f} "
                   f"[{math.exp(m - 1.96 * s):.2f}, {math.exp(m + 1.96 * s):.2f}]; genes with R < 1: {(G.y < 0).sum()}/{len(G)}")
        if band == "Supporting":
            prim_hi = math.exp(m + 1.96 * s)
            loo = [math.exp(dl(G.y.drop(i).to_numpy(), G.se.drop(i).to_numpy())[0]) for i in G.index] if len(G) > 2 else []
            if loo:
                out.append(f"  leave-one-gene-out R range {min(loo):.2f}-{max(loo):.2f}")
            out.append(f"  PRIMARY decision: dependence confirmed if upper bound < 1 -> {'CONFIRMED' if prim_hi < 1 else 'NOT confirmed'}")
            for r in G.itertuples():
                out.append(f"    {r.gene:8s} R {math.exp(r.y):.2f} [{math.exp(r.y - 1.96 * r.se):.2f}, {math.exp(r.y + 1.96 * r.se):.2f}]")
    out.append(f"PP3 share of log-evidence kept at curated residues: median {A_.share_kept.median():.2f} (score sets {len(A_)})")
    if len(C):
        Gc = C.groupby("gene").R_ctrl_median.median()
        out.append(f"Negative control (random residue sets, same size): median R per gene " + ", ".join(f"{g} {v:.2f}" for g, v in Gc.items())
                   + f"; across genes median {Gc.median():.2f}")
open("results_fullstudy/a21_mavedb_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
print(S.to_string(index=False))

# ---- Sensitivity (NOT pre-registered; deviation 1 of A21): sparse-cell smoothing and inverted-truth exclusion
sens = ["Sensitivity (not pre-registered):"]
E = A_.copy()
inv = E.LR_PP3 < 1
sens.append(f"  score sets with inverted truth (LR(PP3) < 1, e.g. dominant-negative assays where nonsense behaves like wild type): "
            f"{', '.join(E.loc[inv, 'urn'])}")
degen = E.R_sup > 100
sens.append(f"  score sets with degenerate R (> 100; a joint cell without non-damaging variants): {', '.join(E.loc[degen, 'urn'])}")
for lab, d in [("excluding inverted-truth sets", E[~inv]), ("excluding inverted and degenerate sets", E[~inv & ~degen])]:
    Gs = by_gene(d, "R_sup", "se_sup")
    m, s = dl(Gs.y.to_numpy(), Gs.se.to_numpy())
    sens.append(f"  {lab}: genes {len(Gs)}, sets {len(d)}; random effects R {math.exp(m):.2f} [{math.exp(m - 1.96 * s):.2f}, {math.exp(m + 1.96 * s):.2f}]; "
                + ", ".join(f"{r.gene} {math.exp(r.y):.2f}" for r in Gs.itertuples()))
open("results_fullstudy/a21_mavedb_summary.txt", "a").write("\n".join(sens) + "\n")
print("\n".join(sens))
