"""Replication of the expert-practice stacking result in diagnostic-laboratory practice (review 14, descriptive).
ClinVar 2026-10 germline submissions (SCVs) whose comments list ACMG/AMP codes; missense variants (variant_summary name
carries a missense protein change); codes parsed with strengths, negations skipped. Truth: points from all codes except
PP3, BP4, PM1, PP2, BP1, BP3, PM4, PM5 (P/LP >= 6, B/LB <= -1). Credit of PP3 once PM1 applies, by PP3 strength; gene-clustered
bootstrap; ClinGen expert-panel SCVs excluded so the set is independent of the Evidence Repository analysis.
Run in the WSL repo: python make_lab_stacking.py"""
import gzip, math, re
import numpy as np, pandas as pd

CVT = "work/cvt/"
BASE = {"PVS": 8, "PS": 4, "PM": 2, "PP": 1, "BA": -8, "BS": -4, "BP": -1}
MOD = {"verystrong": 8, "strong": 4, "moderate": 2, "supporting": 1, "standalone": 8}
EXCL = {"PP3", "BP4", "PM1", "PP2", "BP1", "BP3", "PM4", "PM5"}
CODE = re.compile(r"\b(PVS|PS|PM|PP|BA|BS|BP)([1-9])(?:[_ ]?\(?(very[ _]?strong|strong|moderate|supporting|stand[ _-]?alone)\)?)?\b", re.I)
NEG_A = re.compile(r"^\W*(?:criteria\s+)?(?:was\s+|is\s+|were\s+)?(?:not\b|n/a|unmet|absent\b)", re.I)
NEG_B = re.compile(r"\b(?:not|no|without)\W+$", re.I)


def parse(txt):
    ev = {}
    for m in CODE.finditer(txt or ""):
        if NEG_A.match(txt[m.end():m.end() + 25]) or NEG_B.search(txt[max(0, m.start() - 12):m.start()]):
            continue
        fam, num, mod = m.group(1).upper(), m.group(2), m.group(3)
        v = BASE[fam]
        if mod:
            v = MOD[re.sub(r"[ _-]", "", mod.lower())] * (1 if v > 0 else -1)
        ev[fam + num] = v
    return ev


path = CVT + "submission_summary_2026-10.txt.gz"
with gzip.open(path, "rt") as fh:
    n = next(i for i, line in enumerate(fh) if line.startswith("#VariationID\t"))
s = pd.read_csv(path, sep="\t", skiprows=n, header=0, dtype=str, quoting=3, compression="gzip", on_bad_lines="skip",
                usecols=["#VariationID", "ClinicalSignificance", "Description", "OriginCounts", "ReviewStatus", "Submitter",
                         "SubmittedGeneSymbol"]).rename(columns={"#VariationID": "VariationID"})
s = s[s.OriginCounts.fillna("").str.contains("germline", case=False)]
s = s[~s.ReviewStatus.fillna("").str.contains("expert panel|practice guideline", case=False)]
s = s[s.Description.fillna("").str.contains(r"\b(?:PVS|PS|PM|PP|BA|BS|BP)[1-9]\b", regex=True)]   # any ACMG code: no selection on PP3/PM1
vs = pd.read_csv(CVT + "variant_summary_2026-10.txt.gz", sep="\t", dtype=str, compression="gzip",
                 usecols=["VariationID", "Name", "GeneSymbol", "Type", "Assembly"])
vs = vs[(vs.Assembly == "GRCh38") & (vs.Type == "single nucleotide variant")].drop_duplicates("VariationID")
MIS = re.compile(r"\(p\.(?:[A-Z][a-z]{2})\d+(?!Ter|=|del|ins|dup|fs|ext)[A-Z][a-z]{2}\)")
vs = vs[vs.Name.str.contains(MIS, na=False)]
s = s.merge(vs[["VariationID", "GeneSymbol"]], on="VariationID")
ev = s.Description.map(parse)
s["PP3"] = ev.map(lambda e: e.get("PP3", 0))
s["PM1"] = ev.map(lambda e: e.get("PM1", 0))
s["indep"] = ev.map(lambda e: sum(v for c, v in e.items() if c not in EXCL))
s["gene"] = s.GeneSymbol
out = [f"Lab submissions (non-expert-panel, germline, missense SNV) listing any ACMG code: {len(s)}; submitters {s.Submitter.nunique()}"]
both = (s.PP3 > 0) & (s.PM1 > 0)
out.append(f"PP3 applied {int((s.PP3 > 0).sum())}; PM1 applied {int((s.PM1 > 0).sum())}; both {int(both.sum())}; "
           f"PP3 strength when both: " + ", ".join(f"{k} {v}" for k, v in s.loc[both, "PP3"].value_counts().sort_index().items()))
T = s[(s.indep >= 6) | (s.indep <= -1)].copy()
T["truth"] = (T.indep >= 6).astype(int)
out.append(f"independent truth: {len(T)} (P/LP {T.truth.sum()}, B/LB {(1 - T.truth).sum()}, genes {T.gene.nunique()}, submitters {T.Submitter.nunique()})")


def lr(m, y):
    return m[y == 1].mean() / max(m[y == 0].mean(), 1e-9)


def credit(z, sel=None):
    y, B = z.truth.to_numpy(), (z.PM1 > 0).to_numpy()
    A = (z.PP3 > 0).to_numpy() if sel is None else (z.PP3 == sel).to_numpy()
    return math.log(max(lr(A & B, y) / max(lr(B, y), 1e-9), 1e-9)) / math.log(2.406)


def R(z):
    y, A, B = z.truth.to_numpy(), (z.PP3 > 0).to_numpy(), (z.PM1 > 0).to_numpy()
    return lr(A & B, y) / (lr(A, y) * lr(B, y))


rng = np.random.default_rng(20261015)
genes = T.gene.unique()
gi = T.groupby("gene").indices
boot = lambda f: np.nanpercentile([f(T.iloc[np.concatenate([gi[g] for g in rng.choice(genes, len(genes))])]) for _ in range(1000)], [2.5, 97.5])
r0 = R(T)
lo, hi = boot(R)
out.append(f"R(PP3, PM1) {r0:.2f} [{lo:.2f}, {hi:.2f}]")
for lab, sel in (("any strength", None), ("Supporting", 1), ("Moderate", 2), ("Strong", 4)):
    m = (T.PM1 > 0) & ((T.PP3 > 0) if sel is None else (T.PP3 == sel))
    if m.sum() < 5:
        out.append(f"PP3 {lab}: n {int(m.sum())} (too few)")
        continue
    c0 = credit(T, sel)
    lo, hi = boot(lambda z: credit(z, sel))
    out.append(f"PP3 {lab:12s} beside PM1: earns {c0:.2f} points [{lo:.2f}, {hi:.2f}]; applied {T.loc[m, 'PP3'].mean():.2f}; "
               f"n {int(m.sum())}; B/LB {int((m & (T.truth == 0)).sum())}")
open("results_fullstudy/lab_stacking_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
