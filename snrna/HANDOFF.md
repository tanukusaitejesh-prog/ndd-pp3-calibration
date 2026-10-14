# Handoff (14 Oct 2026, evening): state, decisions, remaining tasks

Branch `claude/snrna-vep-v2-framework-b0abng`, folder `snrna/`. **Every number is logged in `PREREG.md` s14** (dated
deviation log, newest rows at the bottom) and every amendment's text is in s15 (A1-A11). Read those first.
Plain-language explainer for non-experts: https://claude.ai/code/artifact/d74d8ab6-94f8-4562-9912-ee585007f77d
Manuscript draft: `paper/MANUSCRIPT.md`. Novelty audit: `NOVELTY_AUDIT_2026-10-14.md`.

## The paper in one paragraph
Benchmarks that use population-common variants as "benign" inflate variant-effect predictors by a score-specific
**rarity-only null** (A9 non-coding: 40-58% of above-chance AUROC for CADD, GPN-MSA, phyloP447, ReMM; A10 missense:
AlphaMissense 0.920 raw vs 0.694 rarity-matched, inflation +0.26). Measured beyond that null, **RNA-family (Rfam)
profiles** predict pathogenic variants in spliceosomal snRNA genes where genome-wide scores do not: two blind tests on
RNU2-2 (sealed 22 P/LP; literature-only 16 P/LP) and RNU4-2 SGE function. A frozen model (v4) gives ClinGen-calibrated
PP3 Supporting thresholds in RNU2-2, RNU4-2, RNU4ATAC (CADD: none). Pre-registered failures reported: v2 primary model;
A7 (no advantage on RMRP/TERC/SNORD118/RNU7-1).

## Key numbers (all in PREREG s14)
| Result | Number |
|---|---|
| Blind RNU2-2 (sealed, 22 P vs 127 L3) | v4 0.892, null 0.783, excess +0.108 (p .014); Rfam alt rarity excess +0.175 (p .001); CADD 0.501, excess -0.03 |
| A8 blind (16 new literature P/LP) | v4 0.910, excess +0.133 (p .012), v4-CADD +0.298 [.153,.426] - MET; broader label rules: v4>CADD holds, excess only a trend (p ~.07) |
| RNU4-2 SGE (435 SNVs) | Rfam ref col freq rho 0.459, alt rarity 0.389, CADD 0.088, phyloP -0.024 |
| Position-only baselines | 0.73-0.75 vs v4 0.892 (all deltas exclude 0) |
| Circularity | no PP3/BP4 on any blind positive; leave-position-out constraint 0.893 |
| A6 inheritance (literature annotation) | RNU2-2 AUROC 0.839, both genes 0.892 (small dominant n) |
| v4 calibration | pooled PP3 Sup >= 0.362, BP4 Sup <= -1.35; per gene RNU2-2 >= 0.216, RNU4-2 >= 0.189, RNU4ATAC >= 0.88; 37/279 VUS PP3 (own gene) |
| A7 (other RNA genes) | NOT MET; RMRP alt rarity 0.488 vs CADD 0.863 |
| A9 (1,787 ClinVar non-coding P/LP) | null CADD .632, GPN-MSA .676, phyloP .647, ReMM .618; H1 MET; rank changes are ties |
| A10 missense (1,484) | AlphaMissense raw .920 null .762 rarity-matched .694; GPN-MSA .870 -> .555; ClinVar B/LB negatives do NOT fix it; H1 MET |
| Mechanism | hypermutability does NOT explain where CADD fails (refuted); cause open |

## Claims withdrawn or qualified (do not reintroduce)
- "0.89 vs CADD 0.50" as the headline: report excess over the rarity null instead.
- "Family-profile advantage generalises to RNA genes" (A7 failed).
- "Failure is due to multicopy/hypermutable loci" (refuted by the mechanism table).
- v4 adds nothing over the parameter-free Rfam feature (+0.03 ns): lead with the family profile, v4 = clinical product.
- Gene discovery is largely prior art (Chen 2024, Messaoud medRxiv 2026).

## Still running / remaining tasks
1. **A10 SpliceAI** (DGX): `conda run -n spliceai pip install -q "setuptools<70"` then the spliceai, spliceai_parse, eval
   steps of `run_next.sh part5`. **A11 NCBoost re-analysis** (DGX, `part6`). Log both in PREREG s14.
2. **Get DGX outputs into git** (the DGX cannot push): bundle `results/ a9/fetch_info.json a10/*.json a11/fetch_info.json
   paper/figures release/` and copy it into the repo from another machine; or set up a GitHub token on the DGX.
3. Figures: `pip install matplotlib && python figures.py` on the DGX (8 figures); add an A10 missense panel.
4. `python make_lookup.py` -> `release/lookup.html`. Zenodo upload of `release/snrna-vep-v2-final-*.zip` + lookup page
   (user; set author in zenodo.json). Record the DOI in PREREG s14.
5. Optional: RarityBench (packaged benchmark sets + per-tool nulls + `rarity_null.py`) after A10/A11.
6. Manuscript: expand `paper/MANUSCRIPT.md` to journal length; framing decision pending (clinical-first for Human
   Mutation/GIM vs methods-first for Genome Medicine/Genome Biology). medRxiv preprint after.

## DGX notes
`~/snrna-run/snrna`, conda env `evo2`; `INFERNAL_BIN=/home/jovyan/.conda/envs/evo2/bin`; ViennaRNA 2.7.0 pip-installed;
separate env `spliceai`; data/hg38.fa and data/AlphaMissense_hg38.tsv.gz downloaded. Runner `run_next.sh part1..part6`
logs each step to logs/next_<step>.log. Git identity set; `git config pull.rebase false`.

## Rules that still apply
- Never re-run v4.py or edit results_v4/ (hash-frozen, A5). Post hoc analyses are labelled exploratory in PREREG s14.
- Commit messages: no co-author or model-name lines. Do not publish the seal passphrase.
