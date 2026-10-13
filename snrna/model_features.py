#!/usr/bin/env python3
"""Model inputs for the selection model (amendment A2). CPU only; needs Infernal (cmalign, cmemit) and ViennaRNA 2.7.0.

    INFERNAL_BIN=/path/to/infernal/bin python model_features.py      # writes model_inputs/

Outputs (all keyed by gene `label`; transcript positions n are 1-based, 5'->3'):
  genes.tsv.gz      one row per analysis-set gene: family, paralog cluster, flags, gene-level features, Evo 2 window
  nodes.tsv.gz      one row per graph node: target positions (molecule "target") and partner-RNA positions
  edges.tsv.gz      undirected edges (i, j, type) between node indices within a gene graph
  variants.tsv.gz   one row per possible SNV at target positions: counts, inclusion flags, offset, variant features
  sge_eval.tsv.gz   RNU4-2 SGE labels for evaluation only (never read by training)
  windows.tsv.gz    8,192-nt GRCh38 reference window per gene, for Evo 2 (gene-centred)
No patient labels and no SGE data enter any training input.
"""
import json
import os
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import snrna

OUT = os.environ.get("SNRNA_MODEL_INPUTS", os.path.join(snrna.HERE, "model_inputs"))
EVO2_WINDOW = 8192
CMEMIT_N = 20000
BASES = "ACGU"
COMP_RNA = str.maketrans("ACGU", "UGCA")
ELEMENTS = ["stem", "hairpin", "bulge", "internal", "multiloop", "external"]

# partner RNAs per family (amendment A2): another snRNA gene (canonical copy) or a consensus motif it pairs with
PARTNERS = {
    "U4": [("U6", "gene:RNU6-1")],
    "U6": [("U4", "gene:RNU4-2"), ("U2", "gene:RNU2-1")],
    "U2": [("U6", "gene:RNU6-1"), ("BPS", "seq:UACUAAC")],
    "U1": [("5SS", "seq:CAGGUAAGU")],
    "U4atac": [("U6atac", "gene:RNU6ATAC")],
    "U6atac": [("U4atac", "gene:RNU4ATAC"), ("U12", "gene:RNU12")],
    "U12": [("U6atac", "gene:RNU6ATAC"), ("BPS12", "seq:UCCUUAAC")],
    "U11": [("5SS12", "seq:AUAUCCUUU")],
    "U5": [],
}
MOLECULES = ["target", "U6", "U4", "U2", "BPS", "5SS", "U6atac", "U4atac", "U12", "BPS12", "5SS12"]
EDGE_TYPES = ["backbone", "pair", "inter"]


def log(*a):
    print(*a, flush=True)


def infernal(exe):
    return os.path.join(snrna.INFERNAL_BIN, exe) if snrna.INFERNAL_BIN else exe


# ------------------------------------------------------------------------------------------- structure helpers

def pair_class(a, b):
    """Base-pair class of RNA bases a, b: 'WC', 'GU' or 'mismatch'."""
    p = a + b
    if p in ("AU", "UA", "GC", "CG"):
        return "WC"
    if p in ("GU", "UG"):
        return "GU"
    return "mismatch"


def nested_pairs(pairs):
    """Drop crossing (pseudoknotted) pairs, keeping the earlier-opened helix; pairs as (i, j), i < j."""
    kept = []
    for i, j in sorted(pairs):
        if all(not (a < i < b < j) and not (i < a < j < b) for a, b in kept):
            kept.append((i, j))
    return kept


def loop_types(L, pairs):
    """Element type and helix descriptors for positions 1..L from a pair list.

    Returns dict n -> (element, stack_len, dist_to_helix_end). Unpaired positions are typed by the loop that
    encloses them: hairpin (no inner pair), bulge or internal (one inner pair), multiloop (>= 2), external (none).
    """
    pr = nested_pairs([(min(i, j), max(i, j)) for i, j in pairs])
    partner = {}
    for i, j in pr:
        partner[i], partner[j] = j, i
    out = {}
    # helices: maximal runs of stacked pairs (i, j), (i+1, j-1), ...
    helix_of = {}
    for i, j in sorted(pr):
        if (i - 1, j + 1) in set(pr):
            continue
        run = [(i, j)]
        while (run[-1][0] + 1, run[-1][1] - 1) in set(pr):
            run.append((run[-1][0] + 1, run[-1][1] - 1))
        for k, (a, b) in enumerate(run):
            d = min(k, len(run) - 1 - k)
            helix_of[a] = helix_of[b] = (len(run), d)
    for n in range(1, L + 1):
        if n in partner:
            out[n] = ("stem",) + helix_of[n]
            continue
        enclosing = [(i, j) for i, j in pr if i < n < j]
        if not enclosing:
            out[n] = ("external", 0, 0)
            continue
        i0, j0 = max(enclosing, key=lambda p: p[0])  # innermost enclosing pair
        inner = [(i, j) for i, j in pr if i0 < i and j < j0 and not any(a < i and j < b for a, b in pr if (a, b) != (i0, j0) and i0 < a and b < j0)]
        if not inner:
            out[n] = ("hairpin", 0, 0)
        elif len(inner) == 1:
            i1, j1 = inner[0]
            out[n] = ("bulge" if (i1 == i0 + 1 or j1 == j0 - 1) else "internal", 0, 0)
        else:
            out[n] = ("multiloop", 0, 0)
    return out


def sm_site_positions(seq):
    """Positions in the canonical Sm-binding motif (AU3-6G) in the 3' half of the RNA."""
    pos = set()
    for m in re.finditer(r"(?=(A[U]{3,6}G))", seq):
        if m.start() >= len(seq) // 2:
            pos.update(range(m.start() + 1, m.start() + 1 + len(m.group(1))))
    return pos


# ------------------------------------------------------------------------------------- Rfam consensus profiles

def parse_pfam(path):
    """Pfam-format Stockholm: (seqs dict name->aligned, RF line, SS_cons line)."""
    seqs, rf, ss = {}, None, None
    for line in open(path):
        if line.startswith("#=GC RF"):
            rf = line.split()[-1]
        elif line.startswith("#=GC SS_cons"):
            ss = line.split()[-1]
        elif line.strip() and not line.startswith(("#", "//")):
            k, a = line.split()
            seqs[k] = a
    return seqs, rf, ss


def consensus_profile(cm, n=CMEMIT_N, seed=20261013):
    """Column (A, C, G, U, gap) frequencies and pair frequencies over consensus columns from cmemit samples."""
    with tempfile.TemporaryDirectory() as td:
        sto = os.path.join(td, "e.sto")
        subprocess.run([infernal("cmemit"), "-N", str(n), "-a", "--seed", str(seed), "-o", sto, cm], check=True, capture_output=True)
        # cmemit -a writes interleaved Stockholm; rewrite as Pfam with esl-free parsing
        seqs, rf, ss = {}, "", ""
        for line in open(sto):
            if line.startswith("#=GC RF"):
                rf += line.split()[-1]
            elif line.startswith("#=GC SS_cons"):
                ss += line.split()[-1]
            elif line.strip() and not line.startswith(("#", "//")):
                k, a = line.split()
                seqs[k] = seqs.get(k, "") + a
    cols = [c for c, ch in enumerate(rf) if ch not in ".~"]
    arr = np.array([[s[c].upper() for c in cols] for s in seqs.values()])
    freq = np.zeros((len(cols), 5))
    for b, ch in enumerate(BASES):
        freq[:, b] = (arr == ch).mean(0)
    freq[:, 4] = 1 - freq[:, :4].sum(1)
    partner_col = snrna.wuss_partners(ss)
    cons_partner = {}
    for k, c in enumerate(cols):
        pc = partner_col[c]
        if pc >= 0 and pc in cols:
            cons_partner[k] = cols.index(pc)
    pairfreq = {}
    for k, k2 in cons_partner.items():
        if k < k2:
            a, b = arr[:, k], arr[:, k2]
            ok = np.isin(a, list(BASES)) & np.isin(b, list(BASES))
            tab = {}
            for x in BASES:
                for y in BASES:
                    tab[x + y] = float(((a == x) & (b == y) & ok).sum() / max(ok.sum(), 1))
            pairfreq[(k, k2)] = tab
            pairfreq[(k2, k)] = {y + x: v for x in BASES for y in BASES for v in [tab[x + y]]}
    return freq, cons_partner, pairfreq


def cmalign_columns(cm, seqs):
    """Align sequences to a CM; return per-label list of consensus column (or -1 for insert) for each residue,
    the aligned consensus-column strings (for paralog identity), and bit scores."""
    with tempfile.TemporaryDirectory() as td:
        fa, sto, sc = (os.path.join(td, x) for x in ("in.fa", "out.sto", "scores.txt"))
        names = {f"s{i}": lab for i, lab in enumerate(seqs)}
        with open(fa, "w") as f:
            for k, lab in names.items():
                f.write(f">{k}\n{seqs[lab]}\n")
        subprocess.run([infernal("cmalign"), "--outformat", "Pfam", "--sfile", sc, "-o", sto, cm, fa], check=True, capture_output=True)
        aln, rf, _ = parse_pfam(sto)
        bits = {}
        for line in open(sc):
            if line.strip() and not line.startswith("#"):
                c = line.split()
                bits[names[c[1]]] = float(c[6])
    out, cons_str = {}, {}
    for k, a in aln.items():
        lab = names[k]
        colmap, kcons, cs = [], -1, []
        for col, ch in enumerate(a):
            is_cons = rf[col] not in ".~"  # '.' insert and '~' local-end columns are not consensus
            if is_cons:
                kcons += 1
                cs.append(ch.upper() if ch.isalpha() else "-")
            if ch.isalpha():
                colmap.append(kcons if is_cons else -1)
        assert len(colmap) == len(seqs[lab]), lab
        out[lab], cons_str[lab] = colmap, "".join(cs)
    return out, cons_str, bits


# ----------------------------------------------------------------------------------------- ViennaRNA energies

_MD = None


def _md():
    global _MD
    if _MD is None:
        import RNA
        _MD = RNA.md()
        _MD.temperature = 37.0
    return _MD


def mfe(seq):
    import RNA
    return RNA.fold_compound(seq, _md()).mfe()[1]


def dimer(a, b):
    import RNA
    return RNA.fold_compound(a + "&" + b, _md()).mfe_dimer()


def inter_pairs(target, partner):
    """Intermolecular base pairs (target n, partner m), 1-based, from the RNAcofold MFE dimer."""
    struct, _ = dimer(target, partner)
    s = struct.replace("&", "")
    L = len(target)
    stack, pairs = [], []
    for i, ch in enumerate(s):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            j = stack.pop()
            if j < L <= i:
                pairs.append((j + 1, i - L + 1))
    return pairs


def _variant_energies(args):
    """ddG_fold and ddG_bind per partner for every SNV of one gene (run in a worker process)."""
    label, seq, partner_seqs, snvs = args
    g0 = mfe(seq)
    bind0 = {p: dimer(seq, ps)[1] - g0 - mfe(ps) for p, ps in partner_seqs.items()}
    gp = {p: mfe(ps) for p, ps in partner_seqs.items()}
    rows = []
    for n, alt in snvs:
        s2 = seq[:n - 1] + alt + seq[n:]
        g = mfe(s2)
        r = {"label": label, "n": n, "alt_t": alt, "ddG_fold": g - g0}
        for p, ps in partner_seqs.items():
            r[f"ddG_bind_{p}"] = (dimer(s2, ps)[1] - g - gp[p]) - bind0[p]
        rows.append(r)
    return rows


# ------------------------------------------------------------------------------------------------- main build

def paralog_clusters(cons_str, family_of, threshold=0.90):
    """Single-linkage clusters of genes whose consensus-aligned sequences share >= threshold identity."""
    labels = sorted(cons_str)
    parent = {x: x for x in labels}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            if family_of[a] != family_of[b]:
                continue
            sa, sb = cons_str[a], cons_str[b]
            both = [(x, y) for x, y in zip(sa, sb) if x != "-" and y != "-"]
            if len(both) >= 0.5 * min(len(sa.replace("-", "")), len(sb.replace("-", ""))):
                ident = sum(x == y for x, y in both) / len(both)
                if ident >= threshold:
                    parent[find(a)] = find(b)
    roots = {}
    return {x: roots.setdefault(find(x), f"C{len(roots) + 1:02d}") for x in labels}


def build(workers=4):
    os.makedirs(OUT, exist_ok=True)
    sites = pd.read_csv(os.path.join(snrna.WORK, "sites.tsv.gz"), sep="\t")
    sites = sites[sites.in_analysis_set].reset_index(drop=True)
    gq = pd.read_csv(os.path.join(snrna.RESULTS, "week0_gene_qc.tsv"), sep="\t")
    allg = pd.read_csv(os.path.join(snrna.WORK, "genes_v49.tsv"), sep="\t")
    st = pd.read_csv(os.path.join(snrna.WORK, "structure.tsv"), sep="\t")
    seqs = snrna.transcript_seqs(sites)
    seqs = {k: v.replace("T", "U") for k, v in seqs.items()}
    meta = gq[gq.analysis_set].set_index("label")
    family_of = {lab: snrna.snrna_family(meta.loc[lab, "gene"]) for lab in seqs}

    # partner sequences (canonical gene copies, transcript orientation)
    def gene_seq(symbol):
        r = allg[allg.gene == symbol].iloc[0]
        s = snrna.ucsc_seq(r.chrom, int(r.start) - 1, int(r.end))
        s = s if r.strand == "+" else snrna.revcomp(s)
        return s.replace("T", "U")
    partner_seq = {}
    for fam, plist in PARTNERS.items():
        for name, src in plist:
            kind, val = src.split(":")
            partner_seq[name] = gene_seq(val) if kind == "gene" else val
    log("partners:", {k: len(v) for k, v in partner_seq.items()})

    # consensus profiles and column maps per family
    profiles, colmaps, cons_str, bits = {}, {}, {}, {}
    for fam, acc in snrna.RFAM.items():
        labs = [lab for lab in seqs if family_of[lab] == fam]
        if not labs:
            continue
        cm = os.path.join(snrna.DATA, "rfam", f"{acc}.cm")
        profiles[fam] = consensus_profile(cm)
        cmap, cstr, b = cmalign_columns(cm, {lab: seqs[lab] for lab in labs})
        colmaps.update(cmap)
        cons_str.update(cstr)
        bits.update(b)
        log(f"{fam}: {len(labs)} genes, consensus length {len(profiles[fam][0])}")
    clusters = paralog_clusters(cons_str, family_of)

    # gene table
    s_g2 = snrna.g2_sites()  # PREREG s7.1 inclusion rules, used for the per-gene kappa (same as frozen M0)
    kappa = {lab: snrna.fit_kappa(g.O.to_numpy(float), g.MR.to_numpy(float)) for lab, g in s_g2.groupby("label")}
    genes = []
    for lab in seqs:
        r = meta.loc[lab]
        fam = family_of[lab]
        clen = len(profiles[fam][0])
        L = len(seqs[lab])
        mid = (int(r.start) + int(r.end)) // 2
        genes.append({"label": lab, "gene": r.gene, "family": fam, "cluster": clusters[lab], "chrom": r.chrom,
                      "start": int(r.start), "end": int(r.end), "strand": r.strand, "length": L,
                      "sealed": bool(r.gene in snrna.SEALED_GENES), "u4_family": bool(r.gene in snrna.U4_FAMILY_NAMES),
                      "named": not bool(re.match(r"^U\d", r.gene) or r.gene.startswith("U4atac")),
                      "log_kappa": float(np.log(kappa.get(lab, np.nan))) if lab in kappa else np.nan,
                      "bits_per_col": bits[lab] / clen, "insert_frac": float(np.mean([c < 0 for c in colmaps[lab]])),
                      "len_ratio": L / clen, "win_start": mid - EVO2_WINDOW // 2 + 1})
    genes = pd.DataFrame(genes)

    # graphs
    node_rows, edge_rows, var_feat = [], [], {}
    inter_cache = {}
    for g in genes.itertuples():
        lab, fam = g.label, g.family
        seq = seqs[lab]
        L = len(seq)
        freq, cons_partner, pairfreq = profiles[fam]
        cmap = colmaps[lab]
        sg = st[st.label == lab].set_index("n")
        intra = [(n, int(sg.loc[n, "partner_n"])) for n in sg.index if sg.loc[n, "paired"] and int(sg.loc[n, "partner_n"]) > n]
        lt = loop_types(L, intra)
        sm = sm_site_positions(seq)
        partners = PARTNERS[fam]
        inter = {}
        for name, _ in partners:
            key = (lab, name)
            inter_cache[key] = inter_pairs(seq, partner_seq[name])
            for n, m in inter_cache[key]:
                inter.setdefault(n, []).append((name, m))
        # target nodes 0..L-1
        idx = 0
        for n in range(1, L + 1):
            k = cmap[n - 1]
            prof = freq[k] if k >= 0 else np.array([0, 0, 0, 0, 1.0])
            ent = float(-(prof[:4][prof[:4] > 0] * np.log2(prof[:4][prof[:4] > 0])).sum()) if k >= 0 else 2.0
            el, stack_len, dend = lt[n]
            row = {"label": lab, "node": idx, "molecule": "target", "n": n, "base": seq[n - 1], "cons_col": k,
                   "prof_A": prof[0], "prof_C": prof[1], "prof_G": prof[2], "prof_U": prof[3], "prof_gap": prof[4],
                   "entropy": ent, "insert": int(k < 0), "stack_len": stack_len, "helix_end_dist": dend,
                   "paired_intra": int(bool(sg.loc[n, "paired"])), "paired_inter": int(n in inter),
                   "rel_pos": (n - 1) / max(L - 1, 1), "sm_motif": int(n in sm)}
            for e in ELEMENTS:
                row[f"el_{e}"] = int(el == e)
            node_rows.append(row)
            idx += 1
        for n in range(1, L):
            edge_rows.append((lab, n - 1, n, "backbone"))
        for a, b in intra:
            edge_rows.append((lab, a - 1, b - 1, "pair"))
        # partner nodes
        offset = {}
        for name, _ in partners:
            ps = partner_seq[name]
            offset[name] = idx
            for m in range(1, len(ps) + 1):
                node_rows.append({"label": lab, "node": idx, "molecule": name, "n": m, "base": ps[m - 1], "cons_col": -1})
                idx += 1
            for m in range(1, len(ps)):
                edge_rows.append((lab, offset[name] + m - 1, offset[name] + m, "backbone"))
            for n, m in inter_cache[(lab, name)]:
                edge_rows.append((lab, n - 1, offset[name] + m - 1, "inter"))
        # per-position facts needed by variant features
        for n in range(1, L + 1):
            k = cmap[n - 1]
            pn = int(sg.loc[n, "partner_n"]) if sg.loc[n, "paired"] else -1
            var_feat[(lab, n)] = {"k": k, "partner_base": seq[pn - 1] if pn > 0 else None,
                                  "k_partner": cons_partner.get(k, -1) if k >= 0 else -1,
                                  "inter": [(nm, partner_seq[nm][m - 1]) for nm, m in inter.get(n, [])]}
    nodes = pd.DataFrame(node_rows).fillna(0)
    edges = pd.DataFrame(edge_rows, columns=["label", "i", "j", "type"])

    # variants
    v = sites.copy()
    v["ref_t"] = np.where(v.strand == "+", v.ref, v.ref.map(lambda b: b.translate(snrna.COMP))).astype(str)
    v["alt_t"] = np.where(v.strand == "+", v.alt, v.alt.map(lambda b: b.translate(snrna.COMP))).astype(str)
    v["ref_t"] = v.ref_t.str.replace("T", "U")
    v["alt_t"] = v.alt_t.str.replace("T", "U")
    v["ID"] = v.chrom + "-" + v.pos.astype(str) + "-" + v.ref + "-" + v.alt
    rec = v["filter"].notna()
    passrec = v["filter"] == "PASS"
    af = v.AF.fillna(0.0)
    v["AC_pass"] = np.where(passrec, v.AC.fillna(0), 0).astype(int)
    v["include"] = (~(passrec & (af >= 0.001)) & ~(rec & ~passrec) & v.MR.notna() & (v.roulette_filter != "low")).astype(int)
    v["L3_benign_like"] = (passrec & (v.AC.fillna(0) >= 9)).astype(int)
    v["log_mr"] = np.log(v.MR)
    obs_pos = v[passrec & (v.AC.fillna(0) >= 1)].groupby(["label", "n"]).size()
    v["site_rate_cov"] = [int(obs_pos.get((lab, n), 0)) - int(bool(p and ac >= 1)) for lab, n, p, ac in
                          zip(v.label, v.n, passrec, v.AC.fillna(0))]
    v["transition"] = [int({r, a} in ({"A", "G"}, {"C", "T"})) for r, a in zip(v.ref, v.alt)]
    v["cpg"] = [int((r == "C" and pn[3] == "G") or (r == "G" and pn[1] == "C")) if isinstance(pn, str) and len(pn) == 5 else 0
                for r, pn in zip(v.ref, v.PN)]
    fams = genes.set_index("label").family
    feats = []
    for lab, n, ref_t, alt_t in zip(v.label, v.n, v.ref_t, v.alt_t):
        f = var_feat[(lab, n)]
        freq, _, pairfreq = profiles[fams[lab]]
        k = f["k"]
        r = {"alt_col_freq": freq[k][BASES.index(alt_t)] if k >= 0 else 0.0,
             "ref_col_freq": freq[k][BASES.index(ref_t)] if k >= 0 else 0.0}
        pb = f["partner_base"]
        for side, base in (("ref", ref_t), ("alt", alt_t)):
            cls = pair_class(base, pb) if pb else "unpaired"
            for c in ("WC", "GU", "mismatch", "unpaired"):
                r[f"{side}_intra_{c}"] = int(cls == c)
        kp = f["k_partner"]
        r["alt_pair_freq"] = pairfreq[(k, kp)][alt_t + pb] if (pb and k >= 0 and kp >= 0 and (k, kp) in pairfreq) else 0.0
        r["ref_pair_freq"] = pairfreq[(k, kp)][ref_t + pb] if (pb and k >= 0 and kp >= 0 and (k, kp) in pairfreq) else 0.0
        inter_cls = [pair_class(alt_t, b) for _, b in f["inter"]]
        r["alt_inter_broken"] = int(any(c == "mismatch" for c in inter_cls))
        r["alt_inter_wobble"] = int(any(c == "GU" for c in inter_cls))
        r["n_inter_partners"] = len(f["inter"])
        feats.append(r)
    v = pd.concat([v, pd.DataFrame(feats)], axis=1)

    # ViennaRNA energies (parallel by gene)
    jobs = []
    for lab, g in v.groupby("label"):
        fam = fams[lab]
        jobs.append((lab, seqs[lab], {nm: partner_seq[nm] for nm, _ in PARTNERS[fam]}, list(zip(g.n, g.alt_t))))
    log(f"ViennaRNA energies for {len(v)} SNVs in {len(jobs)} genes ...")
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for out in ex.map(_variant_energies, jobs):
            rows.extend(out)
    en = pd.DataFrame(rows)
    bind_cols = [c for c in en.columns if c.startswith("ddG_bind_")]
    en["ddG_bind_max"] = en[bind_cols].max(axis=1, skipna=True).fillna(0.0) if bind_cols else 0.0
    en["ddG_bind_sum"] = en[bind_cols].sum(axis=1, skipna=True) if bind_cols else 0.0
    v = v.merge(en[["label", "n", "alt_t", "ddG_fold", "ddG_bind_max", "ddG_bind_sum"]], on=["label", "n", "alt_t"], how="left")

    # graph node index of each variant's position
    tnodes = nodes[nodes.molecule == "target"].set_index(["label", "n"]).node
    v["node"] = [int(tnodes[(lab, n)]) for lab, n in zip(v.label, v.n)]
    keep = ["ID", "label", "gene", "chrom", "pos", "ref", "alt", "n", "node", "strand", "ref_t", "alt_t", "AC_pass", "AN", "AF",
            "filter", "include", "L3_benign_like", "log_mr", "site_rate_cov", "transition", "cpg", "cadd_phred",
            "alt_col_freq", "ref_col_freq", "alt_pair_freq", "ref_pair_freq", "alt_inter_broken", "alt_inter_wobble",
            "n_inter_partners", "ddG_fold", "ddG_bind_max", "ddG_bind_sum"] + \
           [f"{s}_intra_{c}" for s in ("ref", "alt") for c in ("WC", "GU", "mismatch", "unpaired")]
    v = v[keep]

    # Evo 2 windows (gene-centred, GRCh38 + strand)
    wins = []
    for g in genes.itertuples():
        wseq = snrna.ucsc_seq(g.chrom, g.win_start - 1, g.win_start - 1 + EVO2_WINDOW)
        assert len(wseq) == EVO2_WINDOW and set(wseq) <= set("ACGTN"), g.label
        wins.append({"label": g.label, "chrom": g.chrom, "win_start": g.win_start, "strand": g.strand, "seq": wseq})
    wins = pd.DataFrame(wins)
    # REF check of every variant against its window
    wmap = wins.set_index("label")
    bad = [i for i, (lab, pos, ref) in enumerate(zip(v.label, v.pos, v.ref)) if wmap.loc[lab, "seq"][pos - wmap.loc[lab, "win_start"]] != ref]
    assert not bad, f"{len(bad)} variants disagree with their Evo 2 window"

    # SGE evaluation labels (evaluation only)
    sge = pd.read_csv(os.path.join(snrna.WORK, "sge.tsv"), sep="\t")
    sge = sge[sge.in_transcript][["ID", "n", "is_snv", "Type", "category", "function_score_mean", "function_score_dip", "depleted",
                                  "strong", "region", "critical", "gnomad_AC", "CADD_score", "ddG_bind", "HGVS_transcript"]]

    for name, df in (("genes", genes), ("nodes", nodes), ("edges", edges), ("variants", v), ("windows", wins), ("sge_eval", sge)):
        df.to_csv(os.path.join(OUT, f"{name}.tsv.gz"), sep="\t", index=False)
    summary = {"genes": len(genes), "clusters": int(genes.cluster.nunique()), "nodes": len(nodes), "edges": len(edges),
               "variants": len(v), "included": int(v.include.sum()), "observed_included": int(((v.AC_pass > 0) & (v.include == 1)).sum()),
               "partners": {k: len(s) for k, s in partner_seq.items()}, "cmemit_samples": CMEMIT_N,
               "edge_types": edges.type.value_counts().to_dict()}
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    log(json.dumps(summary, indent=1))
    return summary


if __name__ == "__main__":
    build()
