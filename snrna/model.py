#!/usr/bin/env python3
"""Selection-supervised snRNA variant-effect model (amendment A2). Needs numpy, pandas, torch (no scipy/sklearn).

Model
  * one graph per gene: target snRNA nucleotides + partner RNAs (U6 for U4, ...); edges = backbone, intramolecular
    consensus pairs, intermolecular pairs (RNAcofold)
  * encoder: L layers of [GINE-style message passing over typed edges] + [global self-attention with a learned
    per-head bias for the edge type between i and j and for their shortest-path distance] + feed-forward
  * variant head: MLP([h_node, alt one-hot, variant features]) -> softplus -> f(v) >= 0 (selection, log scale)
  * gene gate: gamma_g = sigmoid(linear(gene features)) in [0, 1] scales f (unexpressed copies -> no selection)
  * likelihood: AC_v ~ NegBin(mu_v, r), log mu_v = alpha_g + log m_v + beta * s_v - gamma_g * f(v), where m_v is the
    Roulette rate, s_v the site-rate covariate, and alpha_g is profiled out per gene (Poisson moment estimator), so
    only within-gene departures from the mutation model carry signal. f never sees allele counts.
"""
import math
import os

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

MOLECULES = ["target", "U6", "U4", "U2", "BPS", "5SS", "U6atac", "U4atac", "U12", "BPS12", "5SS12"]
EDGE_TYPES = ["backbone", "pair", "inter"]
FAMILIES = ["U1", "U2", "U4", "U4atac", "U5", "U6", "U6atac", "U11", "U12"]
BASES = "ACGU"
NODE_CONT = ["prof_A", "prof_C", "prof_G", "prof_U", "prof_gap", "entropy", "insert", "stack_len", "helix_end_dist",
             "paired_intra", "paired_inter", "rel_pos", "sm_motif"] + \
            [f"el_{e}" for e in ["stem", "hairpin", "bulge", "internal", "multiloop", "external"]]
COVARIATION_NODE = ["prof_A", "prof_C", "prof_G", "prof_U", "prof_gap", "entropy"]
VAR_CONT = ["transition", "cpg", "alt_col_freq", "ref_col_freq", "alt_pair_freq", "ref_pair_freq", "alt_inter_broken",
            "alt_inter_wobble", "n_inter_partners", "ddG_fold", "ddG_bind_max", "ddG_bind_sum"] + \
           [f"{s}_intra_{c}" for s in ("ref", "alt") for c in ("WC", "GU", "mismatch", "unpaired")]
COVARIATION_VAR = ["alt_col_freq", "ref_col_freq", "alt_pair_freq", "ref_pair_freq"]
GENE_CONT = ["log_kappa", "bits_per_col", "insert_frac", "len_ratio", "named"]
MAX_SPD = 16


# --------------------------------------------------------------------------------------------------- data

class Data:
    """All genes as tensors. Feature switches implement the pre-registered ablations."""

    def __init__(self, inputs, evo2_dir=None, pca_dim=64, use_inter=True, use_pairs=True, use_covariation=True,
                 use_evo2_emb=True, use_evo2_delta=True, random_graph_seed=None, shuffle_counts_seed=None,
                 synthetic_counts_seed=None):
        self.genes = pd.read_csv(os.path.join(inputs, "genes.tsv.gz"), sep="\t")
        nodes = pd.read_csv(os.path.join(inputs, "nodes.tsv.gz"), sep="\t")
        edges = pd.read_csv(os.path.join(inputs, "edges.tsv.gz"), sep="\t")
        var = pd.read_csv(os.path.join(inputs, "variants.tsv.gz"), sep="\t")
        if not use_inter:
            # drop partner molecules entirely
            nodes = nodes[nodes.molecule == "target"]
            edges = edges[edges.type != "inter"]
            tmax = nodes.groupby("label").node.max()
            edges = edges[(edges.i <= edges.label.map(tmax)) & (edges.j <= edges.label.map(tmax))]
        if not use_pairs:
            edges = edges[edges.type != "pair"]
        if random_graph_seed is not None:
            edges = randomise_pairs(edges, nodes, random_graph_seed)
        node_cont = [c for c in NODE_CONT if use_covariation or c not in COVARIATION_NODE]
        var_cont = [c for c in VAR_CONT if use_covariation or c not in COVARIATION_VAR]

        # Evo 2 features (optional)
        emb = {}
        if use_evo2_emb and evo2_dir and os.path.isdir(os.path.join(evo2_dir, "emb")):
            emb = load_pca_embeddings(evo2_dir, self.genes, pca_dim)
        self.has_emb = bool(emb)
        delta = None
        if use_evo2_delta and evo2_dir and os.path.exists(os.path.join(evo2_dir, "deltas.tsv")):
            d = pd.read_csv(os.path.join(evo2_dir, "deltas.tsv"), sep="\t").drop_duplicates("ID")
            d["evo2_delta"] = ((d.ll_alt_fwd - d.ll_ref_fwd) + (d.ll_alt_rc - d.ll_ref_rc)) / 2
            delta = var.ID.map(dict(zip(d.ID, d.evo2_delta)))
        if delta is not None and delta.notna().mean() > 0.99:
            var = var.assign(evo2_delta=delta.fillna(delta.median()))
            var_cont = var_cont + ["evo2_delta"]
        self.has_delta = "evo2_delta" in var_cont

        # standardise continuous features over all positions / variants (unsupervised, label-free)
        tn = nodes[nodes.molecule == "target"]
        self.node_mu = tn[node_cont].mean()
        self.node_sd = tn[node_cont].std().replace(0, 1)
        self.var_mu = var[var_cont].mean()
        self.var_sd = var[var_cont].std().replace(0, 1)
        g = self.genes.copy()
        g["log_kappa"] = g.log_kappa.fillna(g.log_kappa.median())
        g["named"] = g.named.astype(float)
        gc = g[GENE_CONT].astype(float)
        gfe = ((gc - gc.mean()) / gc.std().replace(0, 1)).to_numpy()
        fam = np.array([[float(f == x) for x in FAMILIES] for f in g.family])
        self.gene_feat = torch.tensor(np.hstack([gfe, fam]), dtype=torch.float32)

        self.labels = list(g.label)
        self.gidx = {lab: i for i, lab in enumerate(self.labels)}
        self.node_cont, self.var_cont = node_cont, var_cont
        self.graphs = []
        for lab in self.labels:
            nd = nodes[nodes.label == lab].sort_values("node").reset_index(drop=True)
            remap = {old: new for new, old in enumerate(nd.node)}
            ed = edges[edges.label == lab]
            ed = ed[ed.i.isin(remap) & ed.j.isin(remap)]
            x = self._node_features(nd, emb.get(lab), pca_dim)
            ei = np.array([[remap[i] for i in ed.i], [remap[j] for j in ed.j]]) if len(ed) else np.zeros((2, 0), int)
            et = np.array([EDGE_TYPES.index(t) for t in ed.type]) if len(ed) else np.zeros(0, int)
            self.graphs.append({"x": torch.tensor(x, dtype=torch.float32), "ei": torch.tensor(ei, dtype=torch.long),
                                "et": torch.tensor(et, dtype=torch.long), "n_target": int((nd.molecule == "target").sum()),
                                "remap": remap})
        # variants
        var = var.reset_index(drop=True)
        if synthetic_counts_seed is not None:
            # dry runs and tests: counts drawn from a planted model, so no real allele count is ever fitted
            var["AC_pass"] = planted_counts(var, nodes, self.genes, synthetic_counts_seed)[2]
        if shuffle_counts_seed is not None:  # after the synthetic draw, so the control also works in dry runs
            rng = np.random.default_rng(shuffle_counts_seed)
            var["AC_pass"] = var.groupby("label").AC_pass.transform(lambda s: rng.permutation(s.to_numpy()))
        self.var = var
        vf = ((var[var_cont] - self.var_mu) / self.var_sd).fillna(0).to_numpy()
        alt = np.zeros((len(var), 4))
        alt[np.arange(len(var)), [BASES.index(a) for a in var.alt_t]] = 1
        self.v_feat = torch.tensor(np.hstack([alt, vf]), dtype=torch.float32)
        self.v_gene = torch.tensor([self.gidx[x] for x in var.label], dtype=torch.long)
        self.v_node = torch.tensor([self.graphs[self.gidx[lab]]["remap"][n] for lab, n in zip(var.label, var.node)], dtype=torch.long)
        self.v_y = torch.tensor(var.AC_pass.to_numpy(), dtype=torch.float32)
        self.v_inc = torch.tensor(var.include.to_numpy().astype(bool) & var.log_mr.notna().to_numpy())
        self.v_logm = torch.tensor(var.log_mr.fillna(0).to_numpy(), dtype=torch.float32)
        self.v_s = torch.tensor(var.site_rate_cov.to_numpy(), dtype=torch.float32)
        self.d_node = self.graphs[0]["x"].shape[1]
        self.d_var = self.v_feat.shape[1]
        self.d_gene = self.gene_feat.shape[1]
        self._pad()

    def _node_features(self, nd, emb, pca_dim):
        n = len(nd)
        base = np.zeros((n, 4))
        for i, b in enumerate(nd.base):
            if b in BASES:
                base[i, BASES.index(b)] = 1
        mol = np.zeros((n, len(MOLECULES)))
        mol[np.arange(n), [MOLECULES.index(m) for m in nd.molecule]] = 1
        cont = ((nd[self.node_cont] - self.node_mu) / self.node_sd).to_numpy()
        cont[nd.molecule.to_numpy() != "target"] = 0  # partner nodes carry base + molecule identity only
        parts = [base, mol, cont]
        if self.has_emb:
            e = np.zeros((n, pca_dim))
            t = nd.molecule.to_numpy() == "target"
            e[t] = emb[nd.n.to_numpy()[t] - 1]
            parts.append(e)
        return np.hstack(parts).astype(np.float32)

    def _pad(self):
        """Padded batch tensors for attention: features, masks, edge-type matrix, shortest-path matrix."""
        B = len(self.graphs)
        N = max(g["x"].shape[0] for g in self.graphs)
        self.N = N
        X = torch.zeros(B, N, self.d_node)
        mask = torch.zeros(B, N, dtype=torch.bool)
        etm = torch.zeros(B, N, N, dtype=torch.long)       # 0 = no edge, 1 + type otherwise
        spd = torch.full((B, N, N), MAX_SPD + 1, dtype=torch.long)
        src, dst, typ = [], [], []
        for b, g in enumerate(self.graphs):
            n = g["x"].shape[0]
            X[b, :n] = g["x"]
            mask[b, :n] = True
            ei, et = g["ei"], g["et"]
            for (i, j), t in zip(ei.T.tolist(), et.tolist()):
                etm[b, i, j] = etm[b, j, i] = t + 1
                src += [b * N + i, b * N + j]
                dst += [b * N + j, b * N + i]
                typ += [t, t]
            spd[b, :n, :n] = torch.tensor(shortest_paths(n, ei.T.tolist()), dtype=torch.long)
        self.X, self.mask, self.etm, self.spd = X, mask, etm, spd
        self.src = torch.tensor(src, dtype=torch.long)
        self.dst = torch.tensor(dst, dtype=torch.long)
        self.etype = torch.tensor(typ, dtype=torch.long)
        self.v_flat = self.v_gene * N + self.v_node     # flat index of each variant's node

    def to(self, device):
        for k in ("X", "mask", "etm", "spd", "src", "dst", "etype", "v_flat", "v_feat", "v_gene", "v_y", "v_inc", "v_logm",
                  "v_s", "gene_feat"):
            setattr(self, k, getattr(self, k).to(device))
        return self


def planted_counts(var, nodes, genes, seed):
    """Synthetic truth: breaking a base pair (intra or inter) costs 2.5 log units, the Sm motif 1.5; genes without an
    HGNC name are under no selection (gamma = 0). Counts ~ NegBin(mu, r = 3). Returns (f_true, gamma, counts)."""
    sm = nodes[nodes.molecule == "target"].set_index(["label", "n"]).sm_motif
    smv = np.array([sm[(lab, n)] for lab, n in zip(var.label, var.n)])
    broken = ((var.alt_intra_mismatch == 1) | (var.alt_inter_broken == 1)).to_numpy()
    f_true = 2.5 * broken + 1.5 * smv
    named = genes.set_index("label").named.astype(bool)
    gamma = var.label.map(named).to_numpy().astype(float)
    rng = np.random.default_rng(seed)
    log_m = var.log_mr.fillna(var.log_mr.median()).to_numpy()
    mu = np.exp(np.log(30.0) + log_m - gamma * f_true)
    r = 3.0
    return f_true, gamma, rng.negative_binomial(r, r / (r + mu))


def shortest_paths(n, edge_list):
    adj = [[] for _ in range(n)]
    for i, j in edge_list:
        adj[i].append(j)
        adj[j].append(i)
    out = np.full((n, n), MAX_SPD + 1, dtype=np.int64)
    for s in range(n):
        out[s, s] = 0
        frontier, d = [s], 0
        while frontier and d < MAX_SPD:
            d += 1
            nxt = []
            for u in frontier:
                for w in adj[u]:
                    if out[s, w] > d:
                        out[s, w] = d
                        nxt.append(w)
            frontier = nxt
    return out


def randomise_pairs(edges, nodes, seed):
    """Negative control: rewire every pair/inter edge to random endpoints within the same gene."""
    rng = np.random.default_rng(seed)
    out = []
    nmax = nodes.groupby("label").node.max()
    for lab, e in edges.groupby("label"):
        e = e.copy()
        m = e.type != "backbone"
        k = int(m.sum())
        if k:
            e.loc[m, "i"] = rng.integers(0, nmax[lab] + 1, k)
            e.loc[m, "j"] = rng.integers(0, nmax[lab] + 1, k)
            e = e[e.i != e.j]
        out.append(e)
    return pd.concat(out, ignore_index=True)


_PCA_CACHE = {}


def load_pca_embeddings(evo2_dir, genes, k):
    """Per-gene (L, k) PCA projections of the Evo 2 layer embeddings, fitted on all target positions (label-free)."""
    key = (os.path.abspath(evo2_dir), tuple(genes.label), k)
    if key not in _PCA_CACHE:
        _PCA_CACHE[key] = _pca_embeddings(evo2_dir, genes, k)
    return _PCA_CACHE[key]


def _pca_embeddings(evo2_dir, genes, k):
    mats = {}
    for lab in genes.label:
        p = os.path.join(evo2_dir, "emb", f"{lab.replace('|', '_')}.npy")
        if not os.path.exists(p):
            return {}
        mats[lab] = np.load(p).astype(np.float32)
    allm = np.vstack(list(mats.values()))
    mu = allm.mean(0)
    u, s, vt = np.linalg.svd(allm - mu, full_matrices=False)
    comp = vt[:k].T
    scale = (allm - mu) @ comp
    sd = scale.std(0) + 1e-6
    return {lab: ((m - mu) @ comp) / sd for lab, m in mats.items()}


# -------------------------------------------------------------------------------------------------- model

class GraphLayer(nn.Module):
    def __init__(self, d, heads, dropout):
        super().__init__()
        self.h, self.dh = heads, d // heads
        self.mp = nn.Linear(d, d)
        self.edge_emb = nn.Embedding(len(EDGE_TYPES), d)
        self.mp_out = nn.Linear(d, d)
        self.qkv = nn.Linear(d, 3 * d)
        self.att_out = nn.Linear(d, d)
        self.b_edge = nn.Embedding(len(EDGE_TYPES) + 1, heads)
        self.b_spd = nn.Embedding(MAX_SPD + 2, heads)
        self.ff = nn.Sequential(nn.Linear(d, 2 * d), nn.GELU(), nn.Dropout(dropout), nn.Linear(2 * d, d))
        self.n1, self.n2, self.n3 = nn.LayerNorm(d), nn.LayerNorm(d), nn.LayerNorm(d)
        self.drop = nn.Dropout(dropout)

    def forward(self, h, data, use_attention=True):
        B, N, d = h.shape
        flat = h.reshape(B * N, d)
        msg = F.relu(self.mp(flat)[data.src] + self.edge_emb(data.etype))
        agg = torch.zeros_like(flat).index_add_(0, data.dst, msg)
        h = self.n1(h + self.drop(self.mp_out(agg).reshape(B, N, d)))
        if use_attention:
            q, k, v = self.qkv(h).reshape(B, N, 3, self.h, self.dh).permute(2, 0, 3, 1, 4)
            s = q @ k.transpose(-1, -2) / math.sqrt(self.dh)
            s = s + self.b_edge(data.etm).permute(0, 3, 1, 2) + self.b_spd(data.spd).permute(0, 3, 1, 2)
            s = s.masked_fill(~data.mask[:, None, None, :], float("-inf"))
            a = torch.softmax(s, dim=-1)
            o = (a @ v).transpose(1, 2).reshape(B, N, d)
            h = self.n2(h + self.drop(self.att_out(o)))
        h = self.n3(h + self.drop(self.ff(h)))
        return h * data.mask[..., None]


class SelectionModel(nn.Module):
    def __init__(self, d_node, d_var, d_gene, d=96, layers=3, heads=4, dropout=0.2, use_gate=True, use_graph=True,
                 use_attention=True):
        super().__init__()
        self.inp = nn.Linear(d_node, d)
        self.layers = nn.ModuleList([GraphLayer(d, heads, dropout) for _ in range(layers if use_graph else 0)])
        self.var_head = nn.Sequential(nn.Linear(d + d_var, d), nn.GELU(), nn.Dropout(dropout), nn.Linear(d, d // 2), nn.GELU(),
                                      nn.Linear(d // 2, 1))
        # zero-initialised output: f starts constant (absorbed by alpha_g), so any variation in f is learned from the
        # counts, and a fit to uninformative counts (negative controls) stays near-constant instead of echoing features
        nn.init.zeros_(self.var_head[-1].weight)
        nn.init.zeros_(self.var_head[-1].bias)
        self.gate = nn.Linear(d_gene, 1) if use_gate else None
        self.beta = nn.Parameter(torch.zeros(1))     # site-rate covariate (mutation side, not part of f)
        self.log_r = nn.Parameter(torch.tensor(1.0))  # NB size (log); large r -> Poisson
        self.use_attention = use_attention

    def node_embeddings(self, data):
        h = self.inp(data.X) * data.mask[..., None]
        for layer in self.layers:
            h = layer(h, data, self.use_attention)
        return h

    def forward(self, data):
        h = self.node_embeddings(data)
        hv = h.reshape(-1, h.shape[-1])[data.v_flat]
        f = F.softplus(self.var_head(torch.cat([hv, data.v_feat], dim=1)).squeeze(1))
        gamma = torch.sigmoid(self.gate(data.gene_feat).squeeze(1)) if self.gate is not None else torch.ones(data.gene_feat.shape[0], device=f.device)
        return f, gamma


def nb_loglik(y, mu, r):
    return (torch.lgamma(y + r) - torch.lgamma(r) - torch.lgamma(y + 1) + r * torch.log(r / (r + mu)) + y * torch.log(mu / (r + mu) + 1e-30))


def gene_loglik(model, data, f, gamma, gene_mask, poisson=False, use_site_cov=True):
    """Per-variant log-likelihood for included variants of genes in gene_mask (alpha_g profiled out)."""
    sel = data.v_inc & gene_mask[data.v_gene]
    g = data.v_gene[sel]
    eta = data.v_logm[sel] - gamma[g] * f[sel]
    if use_site_cov:
        eta = eta + model.beta * data.v_s[sel]
    y = data.v_y[sel]
    G = gene_mask.shape[0]
    sum_y = torch.zeros(G, device=y.device).index_add_(0, g, y)
    sum_e = torch.zeros(G, device=y.device).index_add_(0, g, torch.exp(eta))
    ok = sum_y[g] > 0                      # genes with no observed variant carry no within-gene information
    alpha = torch.log(sum_y.clamp(min=1e-12)) - torch.log(sum_e.clamp(min=1e-30))
    mu = torch.exp(alpha[g] + eta).clamp(min=1e-12)
    if poisson:
        ll = y * torch.log(mu) - mu - torch.lgamma(y + 1)
    else:
        ll = nb_loglik(y, mu, torch.exp(model.log_r))
    return ll[ok], g[ok]


# ----------------------------------------------------------------------------------------------- training

def train_one(data, train_genes, val_genes, cfg, seed, device="cpu", log=None):
    """Fit one model on train_genes; early stopping on val_genes' per-variant log-likelihood. Returns model, history."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = SelectionModel(data.d_node, data.d_var, data.d_gene, d=cfg["d"], layers=cfg["layers"], heads=cfg.get("heads", 4),
                           dropout=cfg["dropout"], use_gate=cfg.get("gate", True), use_graph=cfg.get("graph", True),
                           use_attention=cfg.get("attention", True)).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg["lr"], weight_decay=cfg["wd"])
    G = len(data.labels)
    tm = torch.zeros(G, dtype=torch.bool, device=device)
    tm[[data.gidx[x] for x in train_genes]] = True
    vm = torch.zeros(G, dtype=torch.bool, device=device)
    vm[[data.gidx[x] for x in val_genes]] = True
    best, best_state, best_ep, bad = -float("inf"), None, 0, 0
    hist = []
    for ep in range(cfg["epochs"]):
        model.train()
        opt.zero_grad()
        f, gamma = model(data)
        ll, _ = gene_loglik(model, data, f, gamma, tm, cfg.get("poisson", False), cfg.get("site_cov", True))
        loss = -ll.mean()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        opt.step()
        if ep % cfg.get("eval_every", 5) == 0 or ep == cfg["epochs"] - 1:
            model.eval()
            with torch.no_grad():
                f, gamma = model(data)
                vll, _ = gene_loglik(model, data, f, gamma, vm, cfg.get("poisson", False), cfg.get("site_cov", True))
                v = float(vll.mean()) if len(vll) else float(-loss)
            hist.append((ep, float(-loss.detach()), v))
            if not math.isfinite(v):
                break  # diverged: keep the best finite state (or the initial one)
            if v > best + 1e-5:
                best, best_ep, bad = v, ep, 0
                best_state = {k: t.detach().clone() for k, t in model.state_dict().items()}
            else:
                bad += 1
                if bad >= cfg.get("patience", 10):
                    break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return model, {"best_val_ll": best, "best_epoch": best_ep, "history": hist, "diverged": best_state is None}


def predict(model, data):
    with torch.no_grad():
        f, gamma = model(data)
        h = model.node_embeddings(data)
    return f.cpu().numpy(), gamma.cpu().numpy(), h.reshape(-1, h.shape[-1])[data.v_flat].cpu().numpy()


def split_validation(train_clusters, frac, seed):
    """Deterministic split of training clusters into fit / early-stopping validation clusters."""
    rng = np.random.default_rng(seed)
    cl = sorted(train_clusters)
    rng.shuffle(cl)
    k = max(1, int(round(frac * len(cl))))
    return cl[k:], cl[:k]


# ------------------------------------------------------------------------------------------ function head

def ridge_fit(X, y, lam):
    Xc = np.hstack([X, np.ones((len(X), 1))])
    A = Xc.T @ Xc + lam * np.diag(np.r_[np.ones(X.shape[1]), 0.0])
    return np.linalg.solve(A, Xc.T @ y)


def ridge_predict(w, X):
    return np.hstack([X, np.ones((len(X), 1))]) @ w


def function_head(Z_train, y_train, groups, Z_all, lams=(0.1, 1, 10, 100, 1000)):
    """Ridge regression of SGE function score on frozen model representations; lambda by leave-one-group-out R^2."""
    mu, sd = Z_train.mean(0), Z_train.std(0) + 1e-6
    Zt, Za = (Z_train - mu) / sd, (Z_all - mu) / sd
    best, best_lam = -np.inf, lams[0]
    for lam in lams:
        pred = np.zeros_like(y_train)
        for gname in np.unique(groups):
            te = groups == gname
            w = ridge_fit(Zt[~te], y_train[~te], lam)
            pred[te] = ridge_predict(w, Zt[te])
        r2 = 1 - ((pred - y_train) ** 2).sum() / ((y_train - y_train.mean()) ** 2).sum()
        if r2 > best:
            best, best_lam = r2, lam
    w = ridge_fit(Zt, y_train, best_lam)
    return ridge_predict(w, Za), {"lambda": best_lam, "cv_r2": best}
