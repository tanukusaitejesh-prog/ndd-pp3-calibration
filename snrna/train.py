#!/usr/bin/env python3
"""Train the selection model for every pre-registered scenario (amendment A2). Resumable: each fit writes its own file.

    python train.py --stage all                 # hp -> main -> ablations -> controls (A100: a few hours)
    python train.py --stage all --quick         # tiny settings for a CPU dry run
    python train.py --stage hp|main|ablations|controls

Scenarios (sealed genes RNU2-2 and RNU5B-1 are never in any training set):
  strict      train without the whole U4 family (RNU4-1, RNU4-2, RNU4ATAC, unnamed U4 loci)   -> Aim 2 on RNU4-2
  loose       train without RNU4-2 only                                                    -> Aim 2, secondary
  lofo        one fit per paralog cluster, trained without that cluster                     -> cross-fitted scores
              for every gene (G3, transfer, inheritance; the sealed genes' final scores)
  ablations   strict scenario with one component removed
  controls    strict scenario with shuffled counts / a randomly rewired structure graph (must collapse)
Each fit also fits the function head (ridge on frozen representations -> RNU4-2 SGE score) unless its held-out set
contains RNU4-2; the head is used only to score genes other than RNU4-2 (major -> minor transfer).
"""
import argparse
import itertools
import json
import os
import time

import numpy as np
import pandas as pd
import torch

import model as M

HERE = os.path.dirname(os.path.abspath(__file__))
SEALED = ("RNU2-2", "RNU5B-1")
U4_FAMILY = ("RNU4-1", "RNU4-2", "RNU4ATAC", "U4")

HP_GRID = {"d": [64, 128], "layers": [2, 3], "dropout": [0.1, 0.3], "wd": [1e-4, 1e-2]}
BASE_CFG = {"lr": 3e-3, "epochs": 600, "patience": 12, "eval_every": 5, "heads": 4}
SEEDS = [11, 22, 33, 44, 55]
ABLATIONS = {
    "no_inter": ({"use_inter": False}, {}),
    "no_pairs": ({"use_pairs": False}, {}),
    "no_structure_graph": ({"use_inter": False, "use_pairs": False}, {}),
    "no_covariation": ({"use_covariation": False}, {}),
    "no_evo2_emb": ({"use_evo2_emb": False}, {}),
    "no_evo2_delta": ({"use_evo2_delta": False}, {}),
    "no_evo2": ({"use_evo2_emb": False, "use_evo2_delta": False}, {}),
    "no_gate": ({}, {"gate": False}),
    "no_site_cov": ({}, {"site_cov": False}),
    "poisson": ({}, {"poisson": True}),
    "no_attention": ({}, {"attention": False}),
    "mlp_only": ({}, {"graph": False}),
}
CONTROLS = {
    "shuffled_counts": ({"shuffle_counts_seed": 7}, {}),
    "random_graph": ({"random_graph_seed": 7}, {}),
}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def gene_sets(genes):
    g = genes.set_index("label")
    sealed = set(g.index[g.gene.isin(SEALED)])
    u4 = set(g.index[g.gene.isin(U4_FAMILY)])
    rnu42 = set(g.index[g.gene == "RNU4-2"])
    return g, sealed, u4, rnu42


def clusters_of(genes, labels):
    m = genes.set_index("label").cluster
    return sorted({m[x] for x in labels})


def fit_and_predict(data, train_labels, cfg, seed, device, genes, sge=None, need_function=True):
    """Fit with a cluster-level validation split; return predictions for all variants and fit info."""
    cl = clusters_of(genes, train_labels)
    fit_cl, val_cl = M.split_validation(cl, 0.15, seed)
    cmap = genes.set_index("label").cluster
    fit = [x for x in train_labels if cmap[x] in fit_cl]
    val = [x for x in train_labels if cmap[x] in val_cl]
    t0 = time.time()
    model, info = M.train_one(data, fit, val, cfg, seed, device)
    f, gamma, Z = M.predict(model, data)
    out = pd.DataFrame({"ID": data.var.ID, "label": data.var.label, "f": f, "gamma": gamma[data.v_gene.cpu().numpy()]})
    info.update({"seconds": round(time.time() - t0, 1), "n_fit_genes": len(fit), "n_val_genes": len(val)})
    info.pop("history", None)
    if need_function and sge is not None:
        Zfull = np.hstack([Z, data.v_feat.cpu().numpy()])
        s = sge[sge.is_snv].merge(pd.DataFrame({"ID": data.var.ID, "row": np.arange(len(data.var))}), on="ID")
        d_pred, finfo = M.function_head(Zfull[s.row.to_numpy()], s.function_score_mean.to_numpy(float),
                                        s.region.astype(str).to_numpy(), Zfull)
        out["d"] = -d_pred  # higher = more damaging (lower predicted function score)
        info["function_head"] = finfo
    return out, info


def save(df, path, info, info_path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path + ".tmp", sep="\t", index=False, float_format="%.6g")
    os.replace(path + ".tmp", path)
    with open(info_path, "w") as f:
        json.dump(info, f, indent=1)


def stage_hp(args, data, genes, device):
    path = os.path.join(args.out, "hp.json")
    if os.path.exists(path):
        log("hp: already done", path)
        return json.load(open(path))
    g, sealed, u4, _ = gene_sets(genes)
    train = [x for x in g.index if x not in sealed and x not in u4]
    cl = clusters_of(genes, train)
    rng = np.random.default_rng(20261013)
    cl = list(cl)
    rng.shuffle(cl)
    folds = [cl[i::3] for i in range(3)]
    grid = [dict(zip(HP_GRID, v)) for v in itertools.product(*HP_GRID.values())]
    if args.quick:
        grid = grid[:2]
    cmap = genes.set_index("label").cluster
    rows = []
    for gi, hp in enumerate(grid):
        cfg = {**BASE_CFG, **hp, **(QUICK if args.quick else {})}
        vals = []
        for k, vf in enumerate(folds):
            fit = [x for x in train if cmap[x] not in vf]
            val = [x for x in train if cmap[x] in vf]
            _, info = M.train_one(data, fit, val, cfg, SEEDS[0], device)
            vals.append(info["best_val_ll"])
        rows.append({**hp, "mean_val_ll": float(np.mean(vals)), "fold_val_ll": vals})
        log(f"hp {gi + 1}/{len(grid)} {hp}: mean held-out log-lik {np.mean(vals):.5f}")
    best = max(rows, key=lambda r: r["mean_val_ll"])
    res = {"grid": rows, "best": {k: best[k] for k in HP_GRID}, "folds": folds}
    with open(path, "w") as f:
        json.dump(res, f, indent=1)
    log("hp: best", res["best"])
    return res


QUICK = {"epochs": 25, "patience": 3, "eval_every": 5, "d": 16, "layers": 1}


def run_scenarios(args, data_cache, genes, sge, device, cfg, scenarios):
    g, sealed, u4, rnu42 = gene_sets(genes)
    seeds = SEEDS[:1] if args.quick else SEEDS
    for name, (data_kw, cfg_kw, train_labels, pred_labels, need_fn) in scenarios.items():
        data = data_cache(**data_kw)
        for seed in seeds:
            path = os.path.join(args.out, "preds", f"{name}__seed{seed}.tsv.gz")
            if os.path.exists(path):
                continue
            c = {**cfg, **cfg_kw}
            out, info = fit_and_predict(data, train_labels, c, seed, device, genes, sge, need_fn)
            out = out[out.label.isin(pred_labels)]
            out.insert(0, "scenario", name)
            out.insert(1, "seed", seed)
            save(out, path, {**info, "scenario": name, "seed": seed, "cfg": c, "n_train_genes": len(train_labels)},
                 path.replace(".tsv.gz", ".json"))
            log(f"{name} seed {seed}: best epoch {info['best_epoch']}, val ll {info['best_val_ll']:.5f}, {info['seconds']} s")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", default=os.path.join(HERE, "model_inputs"))
    ap.add_argument("--evo2", default=os.path.join(HERE, "evo2_run", "features"))
    ap.add_argument("--out", default=os.path.join(HERE, "model_runs"))
    ap.add_argument("--stage", default="all", choices=["hp", "main", "ablations", "controls", "all"])
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--quick", action="store_true", help="tiny settings for a dry run")
    ap.add_argument("--max-lofo", type=int, default=0, help="limit the number of leave-one-cluster-out folds (dry runs)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    genes = pd.read_csv(os.path.join(args.inputs, "genes.tsv.gz"), sep="\t")
    sge = pd.read_csv(os.path.join(args.inputs, "sge_eval.tsv.gz"), sep="\t")
    cache = {}

    def data_cache(**kw):
        key = json.dumps(kw, sort_keys=True)
        if key not in cache:
            log("building data", kw or "(full features)")
            cache[key] = M.Data(args.inputs, args.evo2, **kw).to(args.device)
            d = cache[key]
            log(f"  {len(d.labels)} genes, node features {d.d_node}, variant features {d.d_var}, "
                f"Evo 2 embeddings {'yes' if d.has_emb else 'NO'}, Evo 2 delta {'yes' if d.has_delta else 'NO'}")
        return cache[key]

    full = data_cache()
    with open(os.path.join(args.out, "data_info.json"), "w") as f:
        json.dump({"has_evo2_emb": full.has_emb, "has_evo2_delta": full.has_delta, "d_node": full.d_node, "d_var": full.d_var,
                   "node_features": full.node_cont, "variant_features": full.var_cont, "device": args.device,
                   "torch": torch.__version__, "quick": args.quick}, f, indent=1)
    if not full.has_emb or not full.has_delta:
        log("WARNING: Evo 2 features missing - run evo2_features.py first for the pre-registered model")

    hp = stage_hp(args, full, genes, args.device) if args.stage in ("hp", "all") else json.load(open(os.path.join(args.out, "hp.json")))
    cfg = {**BASE_CFG, **hp["best"], **(QUICK if args.quick else {})}
    g, sealed, u4, rnu42 = gene_sets(genes)
    allg = list(g.index)
    strict_train = [x for x in allg if x not in sealed and x not in u4]
    loose_train = [x for x in allg if x not in sealed and x not in rnu42]

    if args.stage in ("main", "all"):
        sc = {"strict": ({}, {}, strict_train, allg, True), "loose": ({}, {}, loose_train, allg, True)}
        clusters = sorted(g.cluster.unique())
        if args.max_lofo:
            clusters = clusters[:args.max_lofo]
        for c in clusters:
            held = list(g.index[g.cluster == c])
            train = [x for x in allg if x not in sealed and g.cluster[x] != c]
            sc[f"lofo_{c}"] = ({}, {}, train, held, not any(x in rnu42 for x in held))
        run_scenarios(args, data_cache, genes, sge, args.device, cfg, sc)
    if args.stage in ("ablations", "all"):
        sc = {f"abl_{k}": (dk, ck, strict_train, allg, False) for k, (dk, ck) in ABLATIONS.items()}
        if args.quick:
            sc = dict(list(sc.items())[:3])
        run_scenarios(args, data_cache, genes, sge, args.device, cfg, sc)
    if args.stage in ("controls", "all"):
        sc = {f"ctl_{k}": (dk, ck, strict_train, allg, False) for k, (dk, ck) in CONTROLS.items()}
        run_scenarios(args, data_cache, genes, sge, args.device, cfg, sc)
    log("training stages done:", args.stage)


if __name__ == "__main__":
    main()
