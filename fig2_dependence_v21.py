"""Figure 2 (draft 21): dependence of PP3 on PM1 on truth that does not come from the codes.
(a) De novo scale: R for AlphaMissense PP3 (Supporting) x PM1-like, by cohort, marginal / within DDG2P genes / within genes
    (values from Table 1; a12_gclust2.tsv). (b) Functional truth: R for AlphaMissense PP3 (Supporting) x curated PM1 per gene
    (MaveDB, pre-registered A21; a21_mavedb_summary.txt) with the random-residue negative control, and pooled estimates
    (MaveDB A21; ProteinGym). F9 is omitted from panel b: its estimate was degenerate (see Supplementary Note 13)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE, GREEN, GREY, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK})
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1, 1]})

# (a) de novo
rows = [("DD", [(0.44, 0.41, 0.47), (0.69, 0.65, 0.73), (0.75, 0.71, 0.78)]),
        ("Autism", [(0.56, 0.47, 0.66), (0.74, 0.67, 0.82), (0.96, 0.84, 1.10)]),
        ("Public", [(0.67, 0.48, 0.96), (0.75, 0.62, 0.89), (0.96, 0.74, 1.28)])]
kinds = [("all sites", BLUE, "o"), ("within DDG2P genes", ORANGE, "s"), ("within genes, all genes", GREEN, "D")]
y = 0
yt, yl = [], []
for coh, vals in rows:
    for k, ((r, lo, hi), (lab, col, mk)) in enumerate(zip(vals, kinds)):
        yy = y - k * 0.28
        a.plot([lo, hi], [yy, yy], color=col, lw=2, solid_capstyle="round")
        a.plot(r, yy, mk, color=col, ms=7, mec="white", mew=1.2, label=lab if coh == "DD" else None)
    yt.append(y - 0.28); yl.append(coh)
    y -= 1.3
a.axvline(1, color=MUTED, lw=1, ls="--")
a.set_xscale("log"); a.set_xticks([0.4, 0.6, 0.8, 1.0, 1.3]); a.set_xticklabels(["0.4", "0.6", "0.8", "1.0", "1.3"])
a.set_yticks(yt); a.set_yticklabels(yl)
a.set_xlabel("R = LR(PP3 ∧ PM1-like) / [LR(PP3) · LR(PM1-like)]")
a.set_title("a  De novo excess (label-free)", loc="left", fontweight="bold", fontsize=10)
a.legend(frameon=False, loc="upper right", fontsize=8)
a.spines[["top", "right"]].set_visible(False)

# (b) functional
genes = [("GCK", 0.49, 0.43, 0.54, 1.00), ("LDLR", 0.52, 0.44, 0.61, 0.99), ("PTEN", 0.57, 0.51, 0.62, 0.98),
         ("TP53", 0.61, 0.44, 0.84, 1.02)]
pooled = [("Pooled, MaveDB\n(5 genes; pre-registered)", 0.59, 0.44, 0.80, 1.00), ("Pooled, ProteinGym\n(3 genes)", 0.49, 0.40, 0.60, None)]
yy = 0
yt, yl = [], []
for g, r, lo, hi, ctrl in genes:
    b.plot([lo, hi], [yy, yy], color=BLUE, lw=2, solid_capstyle="round")
    b.plot(r, yy, "o", color=BLUE, ms=7, mec="white", mew=1.2, label="curated PM1 residues" if g == "GCK" else None)
    b.plot(ctrl, yy, "o", color=GREY, ms=7, mfc="white", mew=1.6, label="random residues, same number" if g == "GCK" else None)
    yt.append(yy); yl.append(g); yy -= 1
yy -= 0.4
for g, r, lo, hi, ctrl in pooled:
    b.plot([lo, hi], [yy, yy], color=INK, lw=2, solid_capstyle="round")
    b.plot(r, yy, "D", color=INK, ms=8, mec="white", mew=1.2)
    if ctrl:
        b.plot(ctrl, yy, "o", color=GREY, ms=7, mfc="white", mew=1.6)
    yt.append(yy); yl.append(g); yy -= 1
b.axvline(1, color=MUTED, lw=1, ls="--")
b.set_xscale("log"); b.set_xticks([0.4, 0.6, 0.8, 1.0, 1.3]); b.set_xticklabels(["0.4", "0.6", "0.8", "1.0", "1.3"])
b.set_xlim(0.35, 1.4)
b.set_yticks(yt); b.set_yticklabels(yl)
b.set_xlabel("R = LR(PP3 ∧ curated PM1) / [LR(PP3) · LR(curated PM1)]")
b.set_title("b  Functional assays (molecular damage)", loc="left", fontweight="bold", fontsize=10)
b.legend(frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.36), fontsize=8)
b.spines[["top", "right"]].set_visible(False)
fig.text(0.01, 0.01, "AlphaMissense PP3 at Supporting. Bars: 95% intervals (a: gene-clustered; b: position-clustered, pooled by random effects). "
         "R = 1: independent. F9 omitted from b (degenerate estimate).", fontsize=7.5, color=MUTED)
fig.tight_layout(rect=(0, 0.04, 1, 1))
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures", "fig2_dependence_v21.png")
fig.savefig(out, dpi=300)
print("written", out)
