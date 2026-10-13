#!/usr/bin/env python3
"""Evaluate the selection model: Aim 2, ablations, negative controls, gate G3 and the exploratory endpoints of
amendment A2 (PREREG s8, s10, s11). Reads model_runs/preds/*.tsv.gz written by train.py. Never reads sealed/.

    python evaluate.py                                                   # real run -> results/model_*
    python evaluate.py --runs model_runs_synthetic --results results_synthetic   # dry run (numbers meaningless)

Scores (higher = more damaging): the selection model's f, the ensemble mean over seeds (range across seeds reported);
d, the function head (ridge on frozen representations -> RNU4-2 SGE score), never used on RNU4-2 itself.
"""
import argparse
import glob
import hashlib
import json
import os
import re
import subprocess
import time

import numpy as np
import pandas as pd

import snrna as S

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = S.SEED
SNV_ID = re.compile(r"^chr[0-9XY]+-\d+-[ACGT]-[ACGT]$")
AIM2_BASELINES = ("CADD v1.7", "ViennaRNA U4/U6 ddG_bind", "Evo 2 7B (G1)")
CONTROL_RHO, CONTROL_AUC = 0.10, 0.5


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ------------------------------------------------------------------------------------------------ loading

def load_preds(runs):
    files = sorted(glob.glob(os.path.join(runs, "preds", "*.tsv.gz")))
    if not files:
        raise SystemExit(f"no predictions in {runs}/preds (run train.py first)")
    pred = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    info = {}
    for f in files:
        j = json.load(open(f.replace(".tsv.gz", ".json")))
        info[(j["scenario"], j["seed"])] = j
    return pred, info


def ensemble(pred, scenario):
    """Per-variant ensemble over seeds: mean f (the score), its range, mean gamma, mean d where fitted."""
    p = pred[pred.scenario == scenario]
    if not len(p):
        return None
    g = p.groupby("ID", sort=False)
    out = g.agg(label=("label", "first"), f=("f", "mean"), f_min=("f", "min"), f_max=("f", "max"),
                gamma=("gamma", "mean"), n_seeds=("seed", "nunique"))
    if "d" in p.columns and p.d.notna().any():
        out["d"] = g.d.mean()
    return out


def variant_counts(var, inputs, data_info):
    """The counts the model was fitted to: real AC_pass, or the planted synthetic counts of a dry run."""
    seed = data_info.get("synthetic_counts")
    if seed is None:
        return var.AC_pass.to_numpy(float)
    import model as M
    nodes = pd.read_csv(os.path.join(inputs, "nodes.tsv.gz"), sep="\t")
    genes = pd.read_csv(os.path.join(inputs, "genes.tsv.gz"), sep="\t")
    return M.planted_counts(var, nodes, genes, seed)[2].astype(float)


def load_l2(labels_dir):
    paths = [os.path.join(labels_dir, "L2_nonblind.tsv"), os.path.join(labels_dir, "L2_curated_additions.tsv")]
    L = pd.concat([pd.read_csv(p, sep="\t") for p in paths if os.path.exists(p)], ignore_index=True)
    L = L[L.ID.str.match(SNV_ID)]
    sealed = L.gene.isin(S.SEALED_GENES)
    if sealed.any():
        raise SystemExit("non-blind label file contains sealed-gene rows; stop and re-seal")
    return L


# ----------------------------------------------------------------------------------------------- metrics

def metric_row(name, d, t, B, renu_sets):
    """PREREG s6.3/s6.4 metrics on the RNU4-2 SGE SNVs for damage score d (aligned with t)."""
    fs = t.function_score_mean.to_numpy(float)
    dep = t.depleted.to_numpy(bool)
    ok = np.isfinite(d)
    rho = S.spearman(-d[ok], fs[ok])
    ci = S.boot_spearman(-d[ok], fs[ok], B)
    cpos = S.boot_spearman(-d[ok], fs[ok], B, clusters=t.n.to_numpy()[ok])
    auc = S.auroc(dep[ok], d[ok])
    aci = S.boot_auroc(dep[ok], d[ok], B)
    row = {"predictor": name, "n": int(ok.sum()), "rho": rho, "rho_lo": ci[0], "rho_hi": ci[1], "rho_pos_lo": cpos[0],
           "rho_pos_hi": cpos[1], "auroc_depleted": auc, "auroc_lo": aci[0], "auroc_hi": aci[1]}
    for key, (y, mask) in renu_sets.items():
        m = mask & ok
        row[f"renu_{key}"] = S.auroc(y[m], d[m])
        lo, hi = S.boot_auroc(y[m], d[m], B)
        row[f"renu_{key}_lo"], row[f"renu_{key}_hi"] = lo, hi
    return row


def renu_sets(t):
    renu = (t.category == "ReNU_syndrome").to_numpy()
    prim = renu | (~renu & (t.gnomad_AC.to_numpy(float) > 0))
    paper = renu | (t.category == "UKBB/AllofUs").to_numpy()
    return {"primary": (renu, prim), "paper": (renu, paper)}


def heldout_loglik(var, y, ens_by_seed, info, genes_mask_labels):
    """Held-out NB log-likelihood per site (alpha profiled per gene) for the model and for the mutation-only null
    (f = 0); returns mean per-site gain in millinats, averaged over seeds."""
    from scipy.special import gammaln
    sel_lab = set(genes_mask_labels)
    inc = var.include.to_numpy(bool) & var.log_mr.notna().to_numpy() & var.label.isin(sel_lab).to_numpy()
    out = []
    for (scen, seed), p in ens_by_seed.items():
        j = info[(scen, seed)]
        f = var.ID.map(dict(zip(p.ID, p.f))).to_numpy(float)
        gm = var.ID.map(dict(zip(p.ID, p.gamma))).to_numpy(float)
        if not np.isfinite(f[inc]).all():
            continue
        r = float(np.exp(j["log_r"]))
        res = []
        for use_f in (True, False):
            use_cov = use_f and j.get("cfg", {}).get("site_cov", True)
            eta = var.log_mr.to_numpy(float) + (j["beta"] * var.site_rate_cov.to_numpy(float) if use_cov else 0)
            if use_f:
                eta = eta - gm * f
            ll_all = []
            for lab in sel_lab:
                m = inc & (var.label == lab).to_numpy()
                yy = y[m]
                if yy.sum() <= 0:
                    continue
                e = eta[m]
                mu = np.exp(np.log(yy.sum()) - np.log(np.exp(e).sum()) + e)
                if j.get("cfg", {}).get("poisson"):
                    ll = yy * np.log(mu) - mu - gammaln(yy + 1)
                else:
                    ll = (gammaln(yy + r) - gammaln(r) - gammaln(yy + 1) + r * np.log(r / (r + mu))
                          + yy * np.log(mu / (r + mu) + 1e-30))
                ll_all.append(ll)
            res.append(np.concatenate(ll_all).mean() if ll_all else np.nan)
        out.append(1000 * (res[0] - res[1]))
    return float(np.mean(out)) if out else np.nan


# -------------------------------------------------------------------------------------------------- report

class Report:
    def __init__(self, path, banner):
        self.path, self.lines = path, ([banner, ""] if banner else [])

    def add(self, *lines):
        self.lines.extend(lines)

    def write(self):
        with open(self.path, "w") as f:
            f.write("\n".join(self.lines) + "\n")
        print("\n".join(self.lines), flush=True)


def fmt_row(r, renu=True):
    s = (f"{r['predictor']:44s} {r['rho']:6.3f} [{r['rho_lo']:6.3f},{r['rho_hi']:6.3f}] "
         f"AUROC dep {r['auroc_depleted']:.3f} [{r['auroc_lo']:.3f},{r['auroc_hi']:.3f}]")
    if renu:
        s += f"  ReNU {r['renu_primary']:.3f} [{r['renu_primary_lo']:.3f},{r['renu_primary_hi']:.3f}]"
    return s


def aim2(args, pred, t, rep, B, sge_all=None):
    """PREREG s8.2 as clarified by A2: the 435 RNU4-2 SGE SNVs, model vs the best baseline on the same SNVs."""
    rs = renu_sets(t)
    base = {"CADD v1.7": t.CADD_score.to_numpy(float), "ViennaRNA U4/U6 ddG_bind": t.ddG_bind.to_numpy(float)}
    ep = os.path.join(HERE, "evo2_run", "evo2_scores.tsv")
    if os.path.exists(ep):
        e = pd.read_csv(ep, sep="\t").drop_duplicates("ID")
        delta = ((e.ll_alt_fwd - e.ll_ref_fwd) + (e.ll_alt_rc - e.ll_ref_rc)) / 2
        base["Evo 2 7B (G1)"] = -t.ID.map(dict(zip(e.ID, delta))).to_numpy(float)
    dp = os.path.join(args.evo2, "deltas.tsv")
    if os.path.exists(dp):
        e = pd.read_csv(dp, sep="\t").drop_duplicates("ID")
        delta = ((e.ll_alt_fwd - e.ll_ref_fwd) + (e.ll_alt_rc - e.ll_ref_rc)) / 2
        x = -t.ID.map(dict(zip(e.ID, delta))).to_numpy(float)
        if np.isfinite(x).all():
            base["Evo 2 7B (gene window; exploratory)"] = x
    rows = [metric_row(k, v, t, B, rs) for k, v in base.items()]
    if sge_all is not None:  # supervised reference from G1 (fitted on SGE scores; not a competitor, not in the bar)
        floor = dict(zip(sge_all.ID, -S.position_floor(sge_all)))
        rows.append(metric_row("Position-only floor (LOPO; supervised ref.)", t.ID.map(floor).to_numpy(float), t, B, rs))
    scen = [("strict", "Model, strict (no U4 family) [primary]"), ("loose", "Model, loose (no RNU4-2) [secondary]"),
            ("lofo_C04", "Model, RNU4-1/RNU4-2 cluster held out [A2]")]
    model_rows = {}
    for sc, name in scen:
        en = ensemble(pred, sc)
        if en is None:
            continue
        d = t.ID.map(en.f).to_numpy(float)
        r = metric_row(name, d, t, B, rs)
        r["n_seeds"] = int(en.n_seeds.min())
        rows.append(r)
        model_rows[sc] = r
    m = pd.DataFrame(rows)
    m.to_csv(os.path.join(args.results, "model_aim2.tsv"), sep="\t", index=False, float_format="%.4f")
    avail = [k for k in AIM2_BASELINES if k in base]
    best = max(avail, key=lambda k: m.set_index("predictor").rho[k])
    best_rho = m.set_index("predictor").rho[best]
    rep.add("Aim 2 (PREREG s8.2; A2: SNVs only, 435 RNU4-2 SGE SNVs, baselines on the same SNVs)",
            f"  rho = Spearman(-score, function_score_mean); bootstrap {B}, seed {SEED}; ReNU = 12 ReNU SNVs vs gnomAD PASS AC>=1 SNVs",
            "")
    for r in rows:
        rep.add("  " + fmt_row(r))
    if "strict" in model_rows:
        r = model_rows["strict"]
        met_rho = r["rho"] >= best_rho + 0.10
        met_renu = r["renu_primary"] >= 0.85
        verdict = "AIM 2 MET" if (met_rho or met_renu) else "AIM 2 NOT MET"
        rep.add("", f"  Bar: rho >= best baseline ({best}, {best_rho:.3f}) + 0.10 = {best_rho + 0.10:.3f}: "
                    f"strict {r['rho']:.3f} -> {'yes' if met_rho else 'no'}; ReNU AUROC >= 0.85: {r['renu_primary']:.3f} -> "
                    f"{'yes' if met_renu else 'no'}", f"  Decision: {verdict} (strict model, point estimates)")
        if len(avail) < 3:
            rep.add(f"  Note: Evo 2 G1 scores missing; the comparison uses {', '.join(avail)} (PREREG s8.2)")
    # per-region rho, strict
    en = ensemble(pred, "strict")
    if en is not None:
        d = t.ID.map(en.f).to_numpy(float)
        reg = []
        for rg, sub in t.groupby("region"):
            idx = sub.index.to_numpy()
            reg.append({"region": rg, "n": len(sub), "rho_model_strict": S.spearman(-d[idx], sub.function_score_mean),
                        **{f"rho_{k}": S.spearman(-v[idx], sub.function_score_mean) for k, v in base.items() if np.isfinite(v[idx]).all()}})
        pd.DataFrame(reg).to_csv(os.path.join(args.results, "model_aim2_by_region.tsv"), sep="\t", index=False, float_format="%.4f")
        rep.add("", "  Per-region rho (strict model): " + "; ".join(f"{r['region']} {r['rho_model_strict']:.2f} (n={r['n']})" for r in reg))
        # seed spread
        rep.add(f"  Seed spread (strict): mean range of f across {int(en.n_seeds.min())} seeds = "
                f"{float((en.f_max - en.f_min).mean()):.3f} (mean f {float(en.f.mean()):.3f})")
    rep.add("")
    return m, model_rows


def ablations_and_controls(args, pred, info, t, var, y, rep, B):
    rs = renu_sets(t)
    g = pd.read_csv(os.path.join(args.inputs, "genes.tsv.gz"), sep="\t")
    u4 = list(g.label[g.u4_family])
    rows = []
    scen = ["strict"] + sorted(s for s in pred.scenario.unique() if s.startswith("abl_"))
    for sc in scen:
        en = ensemble(pred, sc)
        d = t.ID.map(en.f).to_numpy(float)
        r = metric_row(sc, d, t, B, rs)
        by_seed = {(sc, s): p for s, p in pred[pred.scenario == sc].groupby("seed")}
        r["heldout_u4_gain_mnats"] = heldout_loglik(var, y, by_seed, info, u4)
        rows.append(r)
    ab = pd.DataFrame(rows)
    if len(ab):
        full = ab.set_index("predictor").loc["strict"]
        ab["delta_rho_vs_full"] = ab.rho - full.rho
        ab["delta_renu_vs_full"] = ab.renu_primary - full.renu_primary
        ab.to_csv(os.path.join(args.results, "model_ablations.tsv"), sep="\t", index=False, float_format="%.4f")
        rep.add("Ablations (strict training set; A2). Held-out gain = NB log-lik per site on the U4-family genes, model minus "
                "mutation-only null, millinats", "")
        for r in ab.itertuples():
            rep.add(f"  {r.predictor:28s} rho {r.rho:6.3f} ({r.delta_rho_vs_full:+.3f})  AUROC dep {r.auroc_depleted:.3f}  "
                    f"ReNU {r.renu_primary:.3f} ({r.delta_renu_vs_full:+.3f})  held-out gain {r.heldout_u4_gain_mnats:7.2f}")
        rep.add("")
    rows = []
    for sc in sorted(s for s in pred.scenario.unique() if s.startswith("ctl_")):
        en = ensemble(pred, sc)
        d = t.ID.map(en.f).to_numpy(float)
        r = metric_row(sc, d, t, B, rs)
        r["collapsed"] = bool(abs(r["rho"]) < CONTROL_RHO and r["auroc_lo"] <= CONTROL_AUC <= r["auroc_hi"])
        rows.append(r)
    ct = pd.DataFrame(rows)
    if len(ct):
        ct.to_csv(os.path.join(args.results, "model_controls.tsv"), sep="\t", index=False, float_format="%.4f")
        rep.add("Negative controls (PREREG s11): each must give |rho| < 0.1 and an AUROC (depleted) interval containing 0.5", "")
        for r in ct.itertuples():
            rep.add(f"  {r.predictor:28s} rho {r.rho:6.3f}  AUROC dep {r.auroc_depleted:.3f} [{r.auroc_lo:.3f},{r.auroc_hi:.3f}] "
                    f"-> {'PASS (collapsed)' if r.collapsed else 'FAIL (did not collapse)'}")
        rep.add("")
    return ab, ct


def unit_scores(pred, var, genes, labels, scenario_of):
    """Concatenate cross-fitted scores for a set of gene labels; scenario_of(label) names the fit that held it out."""
    parts = []
    for lab in labels:
        en = ensemble(pred, scenario_of(lab))
        if en is None:
            continue
        v = var[var.label == lab][["ID", "label", "gene", "cadd_phred", "L3_benign_like"]].copy()
        v["f"] = v.ID.map(en.f)
        v["gamma_f"] = v.ID.map(en.f * en.gamma)
        v["d"] = v.ID.map(en.d) if "d" in en.columns else np.nan
        parts.append(v)
    return pd.concat(parts, ignore_index=True) if parts else None


def g3(args, pred, var, rep, B):
    """PREREG s8.3. Each gene is scored by the fit that held out its whole paralog cluster (A2: stricter than
    leave-one-gene-out). Units: RNU4ATAC, RNU6ATAC, RNU12, named RNU6 paralogs pooled; eligible with >= 10 L2 SNVs."""
    genes = pd.read_csv(os.path.join(args.inputs, "genes.tsv.gz"), sep="\t")
    cl = genes.set_index("label").cluster
    L = load_l2(args.labels)
    plp = L[L.label == "P/LP"].drop_duplicates("ID")
    blb = L[L.label == "B/LB"].drop_duplicates("ID")
    named = genes[genes.named]
    units = {"RNU4ATAC": ["RNU4ATAC"], "RNU6ATAC": ["RNU6ATAC"], "RNU12": ["RNU12"],
             "RNU6 paralogs (pooled)": sorted(x for x in named.gene if re.match(r"^RNU6-\d+$", x))}
    rows, boot_sets = [], []
    overlap = []
    for u, gl in units.items():
        labs = list(genes.label[genes.gene.isin(gl)])
        sc = unit_scores(pred, var, genes, labs, lambda lab: f"lofo_{cl[lab]}")
        if sc is None:
            rows.append({"unit": u, "status": "no predictions"})
            continue
        pos = sc.ID.isin(plp.ID)
        neg = sc.L3_benign_like.astype(bool) & ~pos
        overlap += list(sc.ID[pos & sc.L3_benign_like.astype(bool)])
        ok = (pos | neg) & sc.f.notna() & sc.cadd_phred.notna()
        y = pos[ok].to_numpy()
        n_pos = int(y.sum())
        r = {"unit": u, "genes": ",".join(gl), "n_L2_snv": n_pos, "n_L3": int((~y).sum()), "eligible": n_pos >= 10}
        if n_pos and (~y).sum():
            fm, cd = sc.f[ok].to_numpy(float), sc.cadd_phred[ok].to_numpy(float)
            r.update({"auroc_model": S.auroc(y, fm), "auroc_cadd": S.auroc(y, cd)})
            r["delta"] = r["auroc_model"] - r["auroc_cadd"]
            if sc.d[ok].notna().all():
                dd = sc.d[ok].to_numpy(float)
                r["auroc_function_head_d"] = S.auroc(y, dd)
                rk = (pd.Series(fm).rank(pct=True) + pd.Series(dd).rank(pct=True)).to_numpy() / 2
                r["auroc_combined_f_d"] = S.auroc(y, rk)
            if r["eligible"]:
                boot_sets.append((y, fm, cd))
        rows.append(r)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(args.results, "model_g3.tsv"), sep="\t", index=False, float_format="%.4f")
    rep.add("Gate G3 (PREREG s8.3): per-unit AUROC of L2 P/LP SNVs vs L3 (gnomAD PASS AC>=9) SNVs, cross-fitted scores", "")
    for r in rows:
        if "auroc_model" in r:
            extra = (f"; function head d {r['auroc_function_head_d']:.3f}, f+d {r['auroc_combined_f_d']:.3f} (exploratory)"
                     if "auroc_function_head_d" in r else "")
            rep.add(f"  {r['unit']:24s} L2 {r['n_L2_snv']:3d} vs L3 {r['n_L3']:3d}: model {r['auroc_model']:.3f}, CADD "
                    f"{r['auroc_cadd']:.3f}, delta {r['delta']:+.3f} {'(eligible)' if r['eligible'] else '(< 10 L2 SNVs: not in G3)'}{extra}")
        else:
            rep.add(f"  {r['unit']:24s} L2 {r.get('n_L2_snv', 0)} SNVs: {r.get('status', 'not evaluable')}")
    if overlap:
        rep.add(f"  Removed from L3 because also in L2: {', '.join(sorted(set(overlap)))}")
    total = int(sum(r["n_L2_snv"] for r in rows if r.get("eligible")))
    if not boot_sets:
        verdict = "G3 NOT EVALUABLE (fewer than 10 L2 SNVs in every week-6 unit): decided at the final evaluation"
        rep.add("", f"  Decision: {verdict}", "")
        return res, verdict
    w = np.array([y.sum() for y, _, _ in boot_sets], float)
    dl = np.array([S.auroc(y, a) - S.auroc(y, c) for y, a, c in boot_sets])
    point = float((w * dl).sum() / w.sum())
    rng = np.random.default_rng(SEED)
    bs = np.empty(B)
    for b in range(B):
        vals = []
        for y, a, c in boot_sets:
            i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
            idx = np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))])
            vals.append(S.auroc(y[idx], a[idx]) - S.auroc(y[idx], c[idx]))
        bs[b] = (w * np.array(vals)).sum() / w.sum()
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    passed = point >= 0.15 and lo > 0
    verdict = "G3 PASS" if passed else "G3 FAIL"
    rep.add("", f"  Weighted delta AUROC (model - CADD) over eligible units ({total} L2 SNVs): {point:+.3f} "
                f"[{lo:+.3f}, {hi:+.3f}] (bootstrap {B} within unit x class)",
            f"  Bar: >= 0.15 and interval excludes 0 -> Decision: {verdict}", "")
    # leave-family-out RNU4ATAC (reported, not gating)
    sc = unit_scores(pred, var, genes, list(genes.label[genes.gene == "RNU4ATAC"]), lambda lab: "strict")
    if sc is not None:
        y = sc.ID.isin(plp.ID)
        n = sc.ID.isin(blb.ID)
        ok = (y | n) & sc.f.notna()
        if y[ok].sum() and n[ok].sum():
            a = S.auroc(y[ok].to_numpy(), sc.f[ok].to_numpy(float))
            lo2, hi2 = S.boot_auroc(y[ok].to_numpy(), sc.f[ok].to_numpy(float), B)
            ac = S.auroc(y[ok].to_numpy(), sc.cadd_phred[ok].to_numpy(float))
            rep.add(f"  Leave-family-out (strict model), RNU4ATAC ClinVar P/LP ({int(y[ok].sum())}) vs ClinVar B/LB "
                    f"({int(n[ok].sum())}): AUROC {a:.3f} [{lo2:.3f}, {hi2:.3f}]; CADD {ac:.3f} (reported, not gating)", "")
    return res, verdict


def exploratory(args, pred, var, rep, B):
    """A2 upgrades: U1 inheritance (dominant vs recessive on RNU4-2), preliminary calibration on non-blind genes."""
    L = load_l2(args.labels)
    plp = L[L.label == "P/LP"].drop_duplicates("ID")
    rep.add("Exploratory (A2; not gating)", "")
    en = ensemble(pred, "strict")
    if en is not None:
        r42 = plp[(plp.gene == "RNU4-2") & plp.inheritance.isin(["dominant", "recessive"])].copy()
        r42["f"] = r42.ID.map(en.f)
        r42 = r42[r42.f.notna()]
        y = (r42.inheritance == "dominant").to_numpy()
        if y.sum() and (~y).sum():
            a = S.auroc(y, r42.f.to_numpy(float))
            lo, hi = S.boot_auroc(y, r42.f.to_numpy(float), B)
            rep.add(f"  U1 inheritance, RNU4-2 (strict model, never trained on U4 family): f separates dominant "
                    f"({int(y.sum())}) from recessive ({int((~y).sum())}) P/LP SNVs with AUROC {a:.3f} [{lo:.3f}, {hi:.3f}] "
                    f"(hypothesis: > 0.5)")
    # preliminary PP3/BP4 on pooled non-blind units with cross-fitted gamma*f
    genes = pd.read_csv(os.path.join(args.inputs, "genes.tsv.gz"), sep="\t")
    cl = genes.set_index("label").cluster
    labs = list(genes.label[genes.gene.isin(["RNU4ATAC", "RNU6ATAC", "RNU12"]) | genes.gene.str.match(r"^RNU6-\d+$")])
    sc = unit_scores(pred, var, genes, labs, lambda lab: f"lofo_{cl[lab]}")
    if sc is not None:
        pos = sc.ID.isin(plp.ID)
        neg = sc.L3_benign_like.astype(bool) & ~pos
        ok = (pos | neg) & sc.gamma_f.notna()
        y, s = pos[ok].to_numpy(), sc.gamma_f[ok].to_numpy(float)
        if y.sum() >= 10 and (~y).sum() >= 10:
            thr = S.pp3_threshold(s, y, 0.10, B=min(B, 1000))
            bp4 = S.pp3_threshold(-s, ~y, 0.10, B=min(B, 1000))
            c = sc.cadd_phred[ok].to_numpy(float)
            thr_c = S.pp3_threshold(c, y, 0.10, B=min(B, 1000))
            rep.add(f"  Preliminary calibration (pooled non-blind minor/U6 units, {int(y.sum())} P vs {int((~y).sum())} L3, "
                    f"prior 0.10): model PP3 Supporting threshold {thr if thr is not None else 'none'}, BP4 "
                    f"{-bp4 if bp4 is not None else 'none'}; CADD PP3 {thr_c if thr_c is not None else 'none'}. "
                    "Final Aim 4 is computed after the single unsealing (PREREG s9-s10).")
        else:
            rep.add(f"  Preliminary calibration: too few labelled SNVs ({int(y.sum())} P, {int((~y).sum())} L3)")
    rep.add("")


def write_scores(args, pred, var):
    """Cross-fitted score table for every variant: lofo (cluster held out) where available, strict otherwise."""
    genes = pd.read_csv(os.path.join(args.inputs, "genes.tsv.gz"), sep="\t")
    cl = genes.set_index("label").cluster
    out = var[["ID", "label", "gene", "n", "ref_t", "alt_t", "cadd_phred"]].copy()
    strict = ensemble(pred, "strict")
    out["f_strict"] = out.ID.map(strict.f) if strict is not None else np.nan
    f, lo, hi, gm, d = (np.full(len(out), np.nan) for _ in range(5))
    src = np.full(len(out), "", dtype=object)
    for lab in out.label.unique():
        en = ensemble(pred, f"lofo_{cl[lab]}")
        if en is None:
            continue
        m = (out.label == lab).to_numpy()
        ids = out.ID[m]
        f[m], lo[m], hi[m], gm[m] = ids.map(en.f), ids.map(en.f_min), ids.map(en.f_max), ids.map(en.gamma)
        if "d" in en.columns:
            d[m] = ids.map(en.d)
        src[m] = f"lofo_{cl[lab]}"
    out["f_crossfit"], out["f_min_seed"], out["f_max_seed"], out["gamma"], out["d_function_head"], out["crossfit_fit"] = f, lo, hi, gm, d, src
    path = os.path.join(args.results, "model_scores.tsv.gz")
    out.to_csv(path, sep="\t", index=False, float_format="%.5g")
    return path


def git_commit():
    try:
        return subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:
        return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=os.path.join(HERE, "model_runs"))
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--labels", default=os.path.join(HERE, "labels"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--boot", type=int, default=2000)
    args = ap.parse_args()
    os.makedirs(args.results, exist_ok=True)
    data_info = json.load(open(os.path.join(args.runs, "data_info.json")))
    synthetic = data_info.get("synthetic_counts") is not None
    banner = ("DRY RUN: model fitted to SYNTHETIC planted counts. Every number below is meaningless and must not be "
              "reported." if synthetic else "")
    if data_info.get("quick"):
        banner = (banner + " " if banner else "") + "QUICK settings (not the pre-registered model)."
    pred, info = load_preds(args.runs)
    var = pd.read_csv(os.path.join(args.inputs, "variants.tsv.gz"), sep="\t")
    y = variant_counts(var, args.inputs, data_info)
    sge = pd.read_csv(os.path.join(args.inputs, "sge_eval.tsv.gz"), sep="\t")
    t = sge[sge.is_snv].reset_index(drop=True)
    rep = Report(os.path.join(args.results, "model_report.txt"), banner)
    rep.add(f"Selection model evaluation; runs {os.path.relpath(args.runs, HERE)}; Evo 2 embeddings "
            f"{'yes' if data_info['has_evo2_emb'] else 'NO'}, deltas {'yes' if data_info['has_evo2_delta'] else 'NO'}; "
            f"device {data_info['device']}; scenarios {pred.scenario.nunique()}, fits {len(info)}", "")
    aim2(args, pred, t, rep, args.boot, sge.reset_index(drop=True))
    ablations_and_controls(args, pred, info, t, var, y, rep, args.boot)
    _, g3_verdict = g3(args, pred, var, rep, args.boot)
    exploratory(args, pred, var, rep, args.boot)
    sp = write_scores(args, pred, var)
    ens = []
    for sc in sorted(pred.scenario.unique()):  # compact per-scenario ensembles (the per-seed files stay on the GPU machine)
        e = ensemble(pred, sc).reset_index()
        e.insert(0, "scenario", sc)
        ens.append(e)
    pd.concat(ens).to_csv(os.path.join(args.results, "model_ensembles.tsv.gz"), sep="\t", index=False, float_format="%.5g")
    rep.add(f"Cross-fitted scores for all {len(var)} SNVs: {os.path.relpath(sp, HERE)}")
    rep.write()
    run = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "git_commit": git_commit(), "synthetic": synthetic,
           "quick": bool(data_info.get("quick")), "g3": g3_verdict, "boot": args.boot,
           "inputs_sha256": {os.path.basename(p): hashlib.sha256(open(p, "rb").read()).hexdigest()
                             for p in sorted(glob.glob(os.path.join(args.inputs, "*.gz")))}}
    with open(os.path.join(args.results, "model_run_info.json"), "w") as f:
        json.dump(run, f, indent=1)


if __name__ == "__main__":
    main()
