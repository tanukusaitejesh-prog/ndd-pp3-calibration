# Verification report, 4 November 2026

Scope: every claim in `paper/MANUSCRIPT.md` and `paper/SUPPLEMENT.md`, the integrity of the pre-registration record, the analysis code, the references, and novelty against work posted up to 4 November 2026. This report lists what was checked, how, what was wrong, and what remains open. Git commit times (UTC) are the authoritative record of order.

## 1. Integrity of the pre-registration record

| Check | Method | Result |
|---|---|---|
| Pre-registration, code, A1, G2 freeze, A2/A2b, A4, A5, A6 timestamps | `ots upgrade` + `ots verify`; Bitcoin block times (Supplementary Table S1, `paper/TIMELINE.md`) | All attested in Bitcoin blocks 971237–971363, 13–14 Oct 2026 |
| v4 frozen before blind labels were decrypted | A5 manifest in block 971362 (14 Oct 05:58 UTC) vs the decryption event log | Order holds |
| Frozen code unchanged since freeze | SHA-256 of `v4.py`, `v4_blind_eval.py` and related files vs the A5/A6 manifests | Match |
| Later amendments (A7–A16) specified before their results | git commit order of each specification vs each result | Order holds for every amendment |
| Post hoc checks of 4 Nov "specified before running" | `sge_classes.py` e2085a0 before 87fc8d8; `alignment_depth.py` b768536 before 0b2117f; `v4_rare_passrate.py` c183f40 (not yet run) | Order holds |
| Commit trailers | `git log` since 4 Nov | No Co-Authored-By lines |

Caveats recorded in the paper: the sealed-label event log (`sealed/SEALED_MANIFEST.jsonl`) is on the DGX and will be released with the scores; the A5 manifest stores output hashes truncated to a 62-character suffix.

## 2. Numbers

- **Recomputed here and matching:** RNU4-2 SGE Spearman correlations (all point estimates identical; intervals now from one saved run, `results/sge_spearman.tsv`); the A8 pattern (Rfam excess +0.16/+0.18, p = 0.001; CADD excess ≈ 0); the A10 and A11 report values; the alignment-depth Spearman (+0.90, recomputed by hand from ranks).
- **Checked against committed result files and the deviation log:** all 268 decimal numbers in the manuscript main text were traced to `paper/SUPPLEMENT.md`, `PREREG.md` §14, `results/*`, or the audited draft (automated check: 0 untraced).
- **Not recomputable here (DGX outputs):** the sealed-set AUROCs and nulls, v4 calibration thresholds, and the raw A9/A13 values. These come from DGX logs recorded in the deviation log at the time they were produced.

## 3. Code audit

An independent audit found no leakage (test suite 49/49 passing on 4 Nov, including the new cluster-bootstrap test) into the blind RNU2-2 result:
- RNU2-2 is alone in its paralog cluster;
- no training ID is in the blind genes;
- within-gene ranks use no labels;
- local-constraint self-counting changes the feature AUROC only from 0.745 to 0.748;
- the phyloP447 = GPN-MSA = 0.724 tie on the sealed set is a coincidence: the same Mann-Whitney U arises from different rank vectors.

## 4. Flaws found and fixed

| # | Flaw | Fix |
|---|---|---|
| 1 | Title and abstract said "calibrated" although the v4 calibration is post hoc with population-observed negatives | "Calibrated" dropped from the title; thresholds called provisional; the pre-specified rare-variant check is named as pending |
| 2 | L3 negatives were called "common" | Corrected: population-observed, median AF 1.8 × 10⁻⁴, 5.4% ≥ 1% |
| 3 | Claim that "40–58% of above-chance AUROC reflects rarity" | Withdrawn: the null pool is not cleaned or consequence-matched, so the share is an upper bound (A15 grid overstates it by +0.32 on average) |
| 4 | A6 inheritance rule presented as a predictor | Circular (the variant's own allele count alone gives AUROC 0.887); moved to Supplement S7 as non-informative |
| 5 | A15 window bootstrap collapsed repeated negatives and did not resample published negatives | Code fixed (`--cluster`) with a test; intervals removed from the main text; all A9–A13 rankings now descriptive; corrected intervals pending on the DGX |
| 6 | Discussion claimed the NCBoost reordering had "an interval excluding a tie" | Removed (it relied on the buggy bootstrap) |
| 7 | The A8 label-rule sensitivity table used 122 L3, not 125; it was described as differing by "≤ 0.004" | Explained in Supplement S3, with the exact differences: up to 0.007 in excess, 0.002 in Δ vs CADD |
| 8 | Sealed-set rarity tests were post hoc across 12 scores, with no multiplicity statement | Bonferroni values given: the Rfam features stay significant (0.012, 0.024); v4 does not (about 0.17); the A8 v4 test was pre-registered |
| 9 | The two blind tests share negatives (125 of 127) without saying so | Stated in Results and Limitations |
| 10 | Overclaims in the earlier draft | Fixed: "mostly dominant"; "CADD reached none" (true only in RNU2-2 and RNU4-2); "genome-wide scores at chance" (GPN-MSA reaches 0.60 on SGE); the tie inside n.1–67 (5′ proximity 0.860 vs v4 0.780) is now reported; the v4 SGE ρ of 0.382 is disclosed as computed without Evo 2 |
| 11 | Reference errors | Fixed: Greene 2026 is recessive RNU2-2; Johnson 2026 is RNU6ATAC/RNU4ATAC autoimmune diabetes; phyloP447 is cited as Kuderna 2024, not Zoonomia; Rfam 15 (Ontiveros-Palacios 2025); Evo 2 is Nature 2026; Jackson 2025 is RNU2-2 and RNU5B-1; Nava 2025 is RNU4-2, RNU5B-1 and the RNU5A-1 candidate; Lu et al. 2026 has 4 authors and is now in Nat Commun |
| 12 | Introduction cited the Roulette model (Pol III mutagenesis) for snRNA hotspots and claimed "hundreds of pseudogenes" | Cited instead: Jackson 2025 and Xue 2026 (Pol II snRNA hypermutability); "many pseudogenes", citing the guidance and Messaoud 2026 |
| 13 | Deviation-log rows and docstrings written on 4 Nov were dated "5 Nov" | Corrected text-only; erratum logged (PREREG §14); git times authoritative |
| 14 | SGE Spearman intervals differed slightly between the manuscript and the Supplement (separate bootstrap runs) | One seeded run saved to `results/sge_spearman.tsv`, cited by both |
| 15 | Rarity-matched A9 intervals in Table 3 came from the buggy bootstrap | Interval column removed; columns labelled upper bound (raw) and lower bound (rarity-matched) |

## 5. Novelty and scoop audit (to 4 Nov 2026)

**Search.**
- full bioRxiv listing 1 Sep–4 Nov 2026 (15,185 records), filtered for snRNA genes and variant-effect or calibration terms;
- medRxiv;
- Europe PMC, including full-text body searches;
- web search;
- full texts of the closest papers.

**Findings.**
- **No scoop.** No study proposes, blind-tests or calibrates an RNA-family (Rfam or covariance-model) profile as variant evidence for snRNA genes or other nuclear non-coding RNA genes. A Europe PMC body search for Rfam + snRNA + pathogenic returns 0 hits.
- **The gap is explicit.** The snRNA classification guidance (D'Souza et al., medRxiv, Aug 2026) states that "no tools have yet been developed that encompass features relevant to snRNA variant deleteriousness, and current genome-wide tools perform poorly". It calls for tools calibrated on snRNA genes and advises against PP3/BP4 until then.
- **Closest related work** (cited and positioned in the manuscript):
  - Miranda 2026: missense-calibrated CADD and AlphaGenome thresholds fail for non-coding variants;
  - Babazade 2026: AlphaGenome and SpliceAI2 calibrated in six non-coding classes, including population controls;
  - Tao & Marsh 2026: predictor evidence feeds back into ClinVar;
  - Lu et al. 2026: genomic heterogeneity inflates predictor performance;
  - Rastogi 2025: frequency-aware missense benchmark;
  - Xue 2026: germline hypermutability of snRNA genes;
  - Messaoud 2026: spliceosome gene and pseudogene prioritisation;
  - nAPOGEE: mitochondrial tRNA/rRNA, no Rfam, no nuclear genes.
- **The rarity-only null** quantifies and corrects a known caveat; it is not a new phenomenon. The manuscript says so.

## 6. Open items (DGX, morning)

1. `python v4_rare_passrate.py --blind sealed/blind_labels.tsv`. If any v4 threshold has LR_pool < 2.08, the paper reports that threshold as "not shown robust to rarity".
2. `python a15_null_validity.py ci a9 --cluster`, then the same for `a10`, `a11` and `a13:traitgym_mendelian`. These give corrected intervals for the descriptive rankings.
3. `sync_dgx_outputs.sh`. Include `results_v4/`, the A9 window files, `sealed/SEALED_MANIFEST.jsonl`, and the v4 thresholds at prior 0.0441.
4. `python a12_new_variants.py build --blind sealed/blind_labels.tsv`, to resolve the three pending RNU2-2/RNU5B-1 candidates.
5. Zenodo deposit of the scores; figures 1–5 from the synced outputs.

None of these can strengthen a claim beyond what the manuscript states. Items 1–2 can only weaken the calibration or ranking statements, which are already worded as provisional or descriptive.

## 6b. External review (5 Nov 2026) and what changed

An external AI review of the 4 Nov packet listed fatal, major and minor issues. Its factual claims were checked against the log and hold. Changes made without the analysis server (each new analysis committed before it was run):

| Review point | Response |
|---|---|
| Title and abstract oversell "pre-registered, blind-validated evidence" | Retitled around what survives: whole-genome conservation fails at key snRNA genes and family alignments recover position-level signal; v4 described as a post-failure model; Rfam features stated as identified after unsealing |
| Rarity null confounded by region and position ascertainment | New `rnu2_position_control.py` on 35 public RNU2-2 positives. The region baseline alone matches the Rfam N0 excess (+0.16 vs +0.18), which confirms the point. After region × pairing matching, Rfam keeps +0.07 (p 0.009) and genome-wide scores ≈ 0; under ±5-nt matching it is a trend (+0.04, p 0.07). Reported in Table 1 |
| SGE bootstrap treats SNVs as independent | New `sge_position_control.py`: position-cluster intervals (Rfam 0.80 [0.72, 0.88]; every difference still excludes 0, also outside the critical region). Within position, Rfam ties GPN-MSA (0.12 vs 0.14), so the claim is stated as position-level |
| Supervised position floor (G1) omitted | Reported: AUROC 0.859 on the same SNVs, as a ceiling that uses SGE data |
| G1, G2, A1 outcomes absent | Reported in Results |
| Pre-specified RNU2-2 score S not shown | S = 0.784 reported; its calibration flagged as pending on the server |
| A14 validity failure not in paper | Reported in Methods and Supplement S8 |
| BP4 unsafe with population negatives in recessive genes | BP4 removed from the main text; PP3 thresholds labelled exploratory |
| Rfam never shown for RNU12, RNU6ATAC | New `per_gene_rfam.py`: Rfam ≥ CADD in all 5 genes; gain where CADD fails (RNU2-2, RNU4-2, RNU12); ties in RNU4ATAC, RNU6ATAC; alignment-agreement prediction ρ −0.80, p 0.067 |
| Benchmark-inflation section is weak and secondary | Moved to Supplement S6; one descriptive paragraph in the main text |
| OpenTimestamps claimed for all amendments; seal log absent | Methods now state A1–A6 Bitcoin-anchored, A7–A16 git only, seal log pending release with a proof |
| Header claimed every number traces to committed files | Server-transcribed numbers now marked ‡ |
| Missing comparators (AlphaGenome, covariation, structure tools, ortholog-only baseline) | Listed as limitations; not run |

Still requiring the analysis server: v4 under the matched nulls on the sealed set; the rare-variant threshold check; the calibration of S and prior-0.0441 results; leave-own-position likelihood ratios; the seal log with an OpenTimestamps proof; corrected benchmark intervals.

## 7. Rating (revised 5 Nov 2026, after the external review)

The external review rated the 4 Nov version 3.0 / 3.0. That was fairer than my 4 / 4: the earlier framing claimed more than the evidence carried.

| | 4 Nov (own) | External (4 Nov version) | Now (own estimate, reframed version) |
|---|---|---|---|
| Novelty | 4 | 3.0 | 3.0–3.5 |
| Q1 potential | 4 | 3.0 | 3.5 |

**Why higher than 3.0 now.** The paper now claims only what survives the review's main test:
- position-level conservation recovered by family alignments, with position-cluster intervals on functional data;
- a residual Rfam excess after region and structure matching that no genome-wide score shows;
- a per-gene pattern consistent with an alignment failure.

The undisclosed results the review found are now reported.

**Why not higher.** It is still one gene for blind patient evidence and one gene for function. The Rfam features were found post hoc. The finding is a conservation repair, not a new principle. There is no new data and the paper has a single author.

**What would change it.** Functional data in a second snRNA gene, or an A12 U5-gene unit passing, would add about +0.5 to Q1. A failed rare-variant threshold check would remove the PP3 statement but not the main result.
