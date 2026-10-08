"""Review 21, laboratory analyses (descriptive; not pre-specified).
(A) Held-out test of "PP2 not added when PM1 applies", tightened: variants used to estimate PP2 credit (2026 PM1-applied
    submissions with truth) are excluded; O2 (other submitters' 2026 majority VUS/LB/B) is primary, O1 (own submission no
    longer P/LP) secondary; one exact method for both (conditional ML odds ratio, exact central CI, central exact p);
    gene-stratified Mantel-Haenszel odds ratio with Robins-Breslow-Greenland CI; second control set restricted to
    Likely Pathogenic calls with PM1 but no PP2 at the same points.
(B) PP2 credit beside PM1 with frequency handled: (i) submissions without BA1/BS1/BS2; (ii) PM2-applied submissions;
    (iii) B/LB truth requiring at least one non-frequency benign criterion. Pooled (gene-clustered) and within-submitter
    (Mantel-Haenszel, submitter-clustered).
(C) Reach of the rule: share of classifications combining PM1 and PP2 that change class; genes; trend by year evaluated.
Run in the WSL repo: python make_review21_lab.py"""
import gzip, math, re
import numpy as np, pandas as pd
from scipy.stats import hypergeom
from scipy.stats.contingency import odds_ratio

CVT = "work/cvt/"
BASE = {"PVS": 8, "PS": 4, "PM": 2, "PP": 1, "BA": -8, "BS": -4, "BP": -1}
MOD = {"verystrong": 8, "strong": 4, "moderate": 2, "supporting": 1, "standalone": 8}
EXCL = {"PP3", "BP4", "PM1", "PP2", "BP1", "BP3", "PM4", "PM5"}
FREQ = {"BA1", "BS1", "BS2", "PM2"}
CODE = re.compile(r"\b(PVS|PS|PM|PP|BA|BS|BP)([1-9])(?:[_ ]?\(?(very[ _]?strong|strong|moderate|supporting|stand[ _-]?alone)\)?)?\b", re.I)
NEG_A = re.compile(r"^(?:/[A-Z]{2,3}\d)?\W*(?:criteria\s+|criterion\s+)?(?:was\s+|is\s+|were\s+|are\s+)?(?:not\b|n/a|unmet|absent\b)", re.I)
NEG_B = re.compile(r"(?:\b(?:not|no|without|neither|nor)\b|does\s+not\s+meet|do\s+not\s+meet|not\s+meet|insufficient\s+for|below\s+the\s+threshold\s+for)[^.;]{0,40}$", re.I)
MIS = re.compile(r"\(p\.(?:[A-Z][a-z]{2})\d+(?!Ter|=|del|ins|dup|fs|ext)[A-Z][a-z]{2}\)")
OUT = "results_fullstudy/review21_lab_summary.txt"
rng = np.random.default_rng(20261021)
LN = math.log(2.406)


def parse(txt):
    ev = {}
    for m in CODE.finditer(txt or ""):
        if NEG_A.match(txt[m.end():m.end() + 25]) or NEG_B.search(txt[max(0, m.start() - 60):m.start()]):
            continue
        fam, num, mod = m.group(1).upper(), m.group(2), m.group(3)
        v = BASE[fam]
        if mod:
            v = MOD[re.sub(r"[ _-]", "", mod.lower())] * (1 if v > 0 else -1)
        ev[fam + num] = v
    return ev


def load(path, cols):
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if line.startswith("#VariationID\t"):
                n, avail = i, line.rstrip("\n").split("\t")
                break
    return pd.read_csv(path, sep="\t", skiprows=n, header=0, dtype=str, quoting=3, compression="gzip", on_bad_lines="skip",
                       usecols=[c for c in cols if c in avail]).rename(columns={"#VariationID": "VariationID"})


COLS = ["#VariationID", "ClinicalSignificance", "Description", "OriginCounts", "ReviewStatus", "Submitter", "SCV", "DateLastEvaluated"]
vs = pd.read_csv(CVT + "variant_summary_2026-10.txt.gz", sep="\t", dtype=str, compression="gzip",
                 usecols=["VariationID", "Name", "GeneSymbol", "Type", "Assembly"])
vs = vs[(vs.Assembly == "GRCh38") & (vs.Type == "single nucleotide variant")].drop_duplicates("VariationID")
vs = vs[vs.Name.str.contains(MIS, na=False)]


def lab_set(path):
    s = load(path, COLS)
    s = s[s.OriginCounts.fillna("").str.contains("germline", case=False)]
    s = s[~s.ReviewStatus.fillna("").str.contains("expert panel|practice guideline", case=False)]
    return s.merge(vs[["VariationID", "GeneSymbol"]], on="VariationID")


def coded(s):
    s = s[s.Description.fillna("").str.contains(r"\b(?:PVS|PS|PM|PP|BA|BS|BP)[1-9]\b", regex=True)].copy()
    ev = s.Description.map(parse)
    for c in ("PP3", "PM1", "PP2", "PM2"):
        s[c] = ev.map(lambda e, c=c: e.get(c, 0))
    s["freqB"] = ev.map(lambda e: any(c in e for c in ("BA1", "BS1", "BS2")))
    s["nfB"] = ev.map(lambda e: any(v < 0 and c not in EXCL | FREQ for c, v in e.items()))
    s["indep"] = ev.map(lambda e: sum(v for c, v in e.items() if c not in EXCL))
    s["pts"] = ev.map(lambda e: sum(e.values())).astype(float)
    s["gene"] = s.GeneSymbol
    return s


def cls(p):
    return np.select([p >= 10, p >= 6, p >= 0, p >= -6], ["P", "LP", "VUS", "LB"], "B")


out = []
n_all = lab_set(CVT + "submission_summary_2026-10.txt.gz")
s = coded(n_all)
T = s[(s.indep >= 6) | (s.indep <= -1)].copy()
T["truth"] = (T.indep >= 6).astype(int)
P = T[T.PM1 > 0]

# ---------------- (B) PP2 credit beside PM1 with frequency handled
def lr(m, y):
    return m[y == 1].mean() / max(m[y == 0].mean(), 1e-9)


def credit_pooled(z):
    y = z.truth.to_numpy(); a = (z.PP2 > 0).to_numpy()
    return math.log(max(lr(a, y), 1e-9)) / LN


def credit_mh(z):
    num = den = 0.0
    for _, g in z.groupby("Submitter"):
        y = g.truth.to_numpy(); a = (g.PP2 > 0).to_numpy()
        n1, n0 = (y == 1).sum(), (y == 0).sum()
        if n1 == 0 or n0 == 0:
            continue
        N = n1 + n0
        num += (a & (y == 1)).sum() * n0 / N; den += (a & (y == 0)).sum() * n1 / N
    return math.log(num / den) / LN if num > 0 and den > 0 else np.nan


def boot(z, f, key, nb=1000):
    u = z[key].unique(); ix = z.groupby(key).indices
    return np.nanpercentile([f(z.iloc[np.concatenate([ix[c] for c in rng.choice(u, len(u))])]) for _ in range(nb)], [2.5, 97.5])


out.append("(B) PP2 credit beside PM1 (points; applied 1) with frequency evidence handled, laboratory truth 2026-10")
T_nf = s[(s.indep >= 6) | ((s.indep <= -1) & s.nfB)].copy()
T_nf["truth"] = (T_nf.indep >= 6).astype(int)
for lab, z in (("all PM1 submissions (reference)", P),
               ("(i) no BA1/BS1/BS2 applied", P[~P.freqB]),
               ("(ii) PM2 applied", P[P.PM2 > 0]),
               ("(iii) B/LB truth needs a non-frequency benign criterion", T_nf[T_nf.PM1 > 0])):
    nb_ = int(((z.PP2 > 0) & (z.truth == 0)).sum())
    if nb_ < 5:
        out.append(f"  {lab}: n {len(z)}; PP2 & B/LB {nb_} (too few)")
        continue
    c0, g_ = credit_pooled(z), boot(z, credit_pooled, "gene")
    m0, m_ = credit_mh(z), boot(z, credit_mh, "Submitter")
    out.append(f"  {lab}: n {len(z)} (B/LB {int((z.truth == 0).sum())}; PP2 & B/LB {nb_}); pooled {c0:.2f} [{g_[0]:.2f}, {g_[1]:.2f}]; "
               f"within-submitter {m0:.2f} [{m_[0]:.2f}, {m_[1]:.2f}]")
open(OUT, "w").write("\n".join(out) + "\n"); print("\n".join(out))

# ---------------- (C) reach of the rule
base = cls(s.pts.to_numpy())
new = cls((s.pts - np.where((s.PM1 > 0) & (s.PP2 > 0), s.PP2, 0)).to_numpy())
both = ((s.PM1 > 0) & (s.PP2 > 0)).to_numpy()
ch = base != new
lpb = both & (base == "LP")
out.append(f"(C) Reach: classifications combining PM1 and PP2 {both.sum()}; changed class {ch.sum()} ({ch.sum() / both.sum():.1%}); "
           f"of their Likely Pathogenic calls ({lpb.sum()}), {int((lpb & (new == 'VUS')).sum())} become VUS "
           f"({(lpb & (new == 'VUS')).sum() / max(lpb.sum(), 1):.1%}); of their Pathogenic calls ({int((both & (base == 'P')).sum())}), "
           f"{int((both & (base == 'P') & (new == 'LP')).sum())} become LP")
g = s.loc[ch, "gene"].value_counts()
out.append(f"  changed classifications span {len(g)} genes and {s.loc[ch, 'Submitter'].nunique()} laboratories; top genes: " +
           ", ".join(f"{k} {v}" for k, v in g.head(12).items()))
s["year"] = pd.to_datetime(s.DateLastEvaluated, errors="coerce").dt.year
pm = s[s.PM1 > 0]
tr = pm.groupby("year").agg(n=("PM1", "size"), pp2=("PP2", lambda x: (x > 0).mean()), pp3mod=("PP3", lambda x: (x >= 2).mean()),
                            pp3=("PP3", lambda x: (x > 0).mean()))
tr = tr[(tr.index >= 2015) & (tr.n >= 50)]
out.append("  by year last evaluated, PM1 submissions: " + "; ".join(
    f"{int(y)} n {int(r.n)} PP2 {r.pp2:.0%} PP3 {r.pp3:.0%} PP3>=Moderate {r.pp3mod:.1%}" for y, r in tr.iterrows()))
yy = pm[(pm.year >= 2015)]
from scipy.stats import spearmanr
for col, lab in (("PP2", "PP2 beside PM1"), ("PP3", "PP3>=Moderate beside PM1")):
    x = (yy[col] > 0) if col == "PP2" else (yy[col] >= 2)
    rho, p = spearmanr(yy.year, x)
    out.append(f"  trend {lab}: Spearman rho {rho:.3f}, p {p:.2g} (submission level)")
open(OUT, "w").write("\n".join(out) + "\n"); print("\n".join(out[-4:]))

# ---------------- (A) held-out test, tightened
o = coded(lab_set(CVT + "submission_summary_2023-01.txt.gz"))
o["LPcall"] = o.ClinicalSignificance.fillna("").str.strip().str.lower().eq("likely pathogenic")
stack = (o.PM1 > 0) & (o.PP2 > 0)
o["affected"] = o.LPcall & stack & (o.pts >= 6) & (o.pts < 10) & (o.pts - o.PP2 < 6)
o["ctrl_all"] = o.LPcall & ~stack & (o.pts >= 6) & (o.pts < 7)
o["ctrl_pm1"] = o["ctrl_all"] & (o.PM1 > 0)
used = set(P.VariationID)
n0 = int(o.affected.sum()), int(o.ctrl_all.sum())
o = o[~o.VariationID.isin(used)].copy()
out.append(f"(A) Held-out test, 2023-01 -> 2026-10. Excluding variants used to estimate PP2 credit ({len(used)} variants): "
           f"rule-affected {n0[0]} -> {int(o.affected.sum())}; controls {n0[1]} -> {int(o.ctrl_all.sum())} "
           f"(of which with PM1 {int(o.ctrl_pm1.sum())})")
n_all["acc"] = n_all.SCV.str.split(".").str[0]
o["acc"] = o.SCV.str.split(".").str[0]
plp = lambda x: x.fillna("").str.contains("pathogenic", case=False) & ~x.fillna("").str.contains("conflict|uncertain|benign", case=False)
nonplp = lambda x: x.fillna("").str.contains("uncertain|benign", case=False)
same = n_all.drop_duplicates("acc").set_index("acc").ClinicalSignificance
o["same_2026"] = o.acc.map(same)
o["O1"] = np.where(o.same_2026.isna(), np.nan, (~plp(o.same_2026.astype(str))).astype(float))
gg = n_all[["VariationID", "Submitter", "ClinicalSignificance"]].copy()
gg["plp"] = plp(gg.ClinicalSignificance).astype(int); gg["non"] = nonplp(gg.ClinicalSignificance).astype(int)
m = o[["VariationID", "Submitter"]].drop_duplicates().merge(gg, on="VariationID", suffixes=("", "_o"))
m = m[m.Submitter != m.Submitter_o]
agg = m.groupby(["VariationID", "Submitter"]).agg(n_o=("plp", "size"), plp_o=("plp", "sum"), non_o=("non", "sum")).reset_index()
o = o.merge(agg, on=["VariationID", "Submitter"], how="left")
o["O2"] = np.where(o.n_o.fillna(0) == 0, np.nan, (o.non_o > o.plp_o).astype(float))


def exact(a, b, c, d):
    r = odds_ratio([[a, b], [c, d]]); ci = r.confidence_interval(0.95)
    N, K, n = a + b + c + d, a + c, a + b
    p = min(1.0, 2 * min(hypergeom.cdf(a, N, K, n), hypergeom.sf(a - 1, N, K, n)))
    return r.statistic, ci.low, ci.high, p


def mh_rbg(z, ex, outc, key):
    R = S = PR = PSQR = QS = 0.0; k = 0
    for _, g_ in z.groupby(key):
        e = g_[ex].to_numpy(); y = g_[outc].to_numpy()
        a = (e & (y == 1)).sum(); b = (e & (y == 0)).sum(); c = (~e & (y == 1)).sum(); d = (~e & (y == 0)).sum()
        n = a + b + c + d
        if (a + b) == 0 or (c + d) == 0:
            continue
        k += 1
        r_, s_ = a * d / n, b * c / n; p_, q_ = (a + d) / n, (b + c) / n
        R += r_; S += s_; PR += p_ * r_; PSQR += p_ * s_ + q_ * r_; QS += q_ * s_
    if R == 0 or S == 0:
        return np.nan, np.nan, np.nan, k
    orr = R / S
    se = math.sqrt(PR / (2 * R * R) + PSQR / (2 * R * S) + QS / (2 * S * S))
    return orr, math.exp(math.log(orr) - 1.96 * se), math.exp(math.log(orr) + 1.96 * se), k


for outc, lab, role in (("O2", "other submitters in 2026 classify VUS/LB/B by majority", "PRIMARY"),
                        ("O1", "own submission no longer P/LP in 2026", "secondary")):
    for ctrl, cl in (("ctrl_all", "controls: LP at 6 points without PP2 beside PM1"), ("ctrl_pm1", "controls: LP at 6 points with PM1, without PP2")):
        z = o[(o.affected | o[ctrl])].dropna(subset=[outc]).copy()
        z["exp"] = z.affected.astype(bool)
        a = int((z.exp & (z[outc] == 1)).sum()); b = int((z.exp & (z[outc] == 0)).sum())
        c = int((~z.exp & (z[outc] == 1)).sum()); d = int((~z.exp & (z[outc] == 0)).sum())
        if a + b == 0 or c + d == 0:
            out.append(f"  {role} {outc}, {cl}: not estimable"); continue
        orr, lo, hi, p = exact(a, b, c, d)
        mg = mh_rbg(z, "exp", outc, "gene"); ms = mh_rbg(z, "exp", outc, "Submitter")
        out.append(f"  {role} {outc} ({lab}); {cl}: affected {a}/{a + b} ({a / (a + b):.1%}) vs controls {c}/{c + d} ({c / max(c + d, 1):.1%}); "
                   f"exact conditional OR {orr:.2f} [{lo:.2f}, {hi:.2f}], central exact p {p:.3f}; gene-stratified MH {mg[0]:.2f} "
                   f"[{mg[1]:.2f}, {mg[2]:.2f}] ({mg[3]} informative genes); submitter-stratified MH {ms[0]:.2f} [{ms[1]:.2f}, {ms[2]:.2f}] "
                   f"({ms[3]} informative submitters)")
open(OUT, "w").write("\n".join(out) + "\n"); print("\n".join(out[-5:]))
