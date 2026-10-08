# Handoff update: draft 19 (8 October 2026, cloud session)

**Current manuscript: `MANUSCRIPT_v19.md` + `SUPPLEMENT_v19.md`** (A21 MaveDB assays added) (A20 de novo rule tests added; pre-registration A20_DENOVO_RULES_PREREG.md) (branch `claude/a18-temporal-reclassification`, with main's drafts 14-16 merged in).
Built on draft 16 from the review-22 roadmap. A cloud draft mistakenly built on v13 is in `archive/`. Pre-registration `prereg_fullstudy/A19_V14_PREREG.md` was committed (ceb1a11) before
H1–H4 were run; deviation log `A19_DEVIATIONS.md`. Reproduce with `bash reproduce_v14.sh` (public ClinVar only);
`REPRODUCE.md` maps every number to a script.

Results: H1 no rule beats the cap (cap retained); H2 PP2 credit beside PM1 0.35 [−0.12, 0.81] (confirmed < 1; PP2 alone
1.12); H3 gene-calibrated PP3 still dependent (R 0.41); H4 R < 1 in every era, BUT negative controls PP1/PP4 show
R ≈ 0.03–0.04 on laboratory truth -> classification-derived truth does not measure dependence; laboratory/panel R and
credit are now presented as descriptions of practice. Cost: uncapped Pathogenic 2,438 DD de novos at 84% vs ≥ 99%;
PP3 × PM1 × PP2 over-counts 2.3 points (DD).

**Not doable in the cloud (needs the WSL column store on X:):**
1. Share kept within genes on the de novo scale, and candidate rules on the de novo scale.
2. Functional-assay redo with gene as unit and VCEP-specification PM1 regions (PTEN numbering, BRCA1 domains).
3. Negative-control criteria on the de novo scale are not available (no de novo analogue of PP1/PP4); consider gene-level
   features unrelated to constraint as controls.

**Author actions:** two annotators fill `results_fullstudy/parser_annotation_sheet_v14.tsv` (then compute kappa);
OSF deposit of A19 (state that it postdates the git commit); obtain SVC v4 draft/pilot materials and verify the
"[verify against the v4 specification]" sentences; names, affiliations, funding, competing interests, contributions;
make the repository public + Zenodo DOI; journal formatting (Genetics in Medicine first target).

---

# Handoff: NDD PP3 study (6 October 2026, evening IST)

Read this first in a new chat, then continue the "Remaining" list in order. The user wants one-line replies, minimal tokens, full autonomy, and a commit and push after every step.

## Where things are

- **Manuscript:** `MANUSCRIPT_v2.md` (draft 3). Review responses: `REVIEW_RESPONSE.md`. Plan and amendments A1–A15: `FULL_STUDY_PLAN.md`.
- **Repository:** private GitHub `tanukusaitejesh-prog/ndd-pp3-calibration`. The WSL repo `/root/ndd_pp3_gate` (distro `ndd-gate`, disk on the X: USB drive) is the git working copy.
- **Commit:** `wsl -d ndd-gate -- bash /mnt/c/Users/saite/OneDrive/Desktop/Neuro/G-Test-LOCKIN/ndd_pp3_gate/ndd_pp3_gate/sync_commit.sh "msg" file...`
  - Run it from PowerShell, because Git Bash mangles `/mnt` paths and `$vars`.
  - It copies the named files from C: into WSL, then commits and pushes.
  - Result files under `results_fullstudy/` are gitignored: `git add -f` them first.
- **Python:** `/root/miniforge3/envs/ndd-gate/bin/python fullstudy.py <cmd>`, run in `/root/ndd_pp3_gate`.
  - Copy the edited `fullstudy.py` from C: first; `sync_commit.sh` does this.
  - Launch long jobs from PowerShell: `wsl -d ndd-gate -- bash -c "cd /root/ndd_pp3_gate && (nohup sh -c '<python ...> > work/x.log 2>&1' >/dev/null 2>&1 &)"`.
  - Run one job at a time; X: is flaky.

## Done this session (all in the manuscript and pushed)

- Simulation (Supplementary Note 3) and A12.4 classification filled in. Draft 3 restructure for reviews 5–6.
- **A12.1** decomposition: overlap is null-side driven.
- **A12.2** all predictors: 40/40 show R < 1 in DD (median 0.51). Population-informed vs other: p = 0.87.
- **A12.3** leave-top-50-genes-out: R 0.48 [0.42, 0.53].
- **A12.5** eRepo: 1,018 classifications with PP3 + PM1; the cap changes 24, the conditional rule 238.
- **A14/A15** models: v3, v4 and v5 all fall below AlphaMissense on ClinVar within gene (DDG2P 0.92/0.91/0.88 vs 0.97) and tie it on autism de novos. Presented as an independent test of the interpretation.
- **Estimand bridge:** ClinVar-scale R in DDG2P genes 0.32 (Supporting), 0.23 (Strong).
- **Independent expert truth** (eRepo, computational and position codes removed; `make_erepo_truth.py`): Strong R 0.17 [0.08, 0.61]. The cap is best calibrated: PP3 earns 1.83 points given PM1, against 2 allowed by the cap and 4 uncapped.
- Supplementary Table 4 (`supp_table4_conditional_pp3.tsv`); MANE concordance 74.4% (limitation 7); autism count reconciliation.

## 7 October 2026, morning: drafts 5 and 6 (pushed)

- **Current manuscript: `MANUSCRIPT_v6.md` + `SUPPLEMENT_v6.md`.** The paper is reframed: de novo benchmarks overstate predictor skill.
  - 39-predictor gene- and position-adjusted benchmark (`make_benchmark.py`, `supp_table7_benchmark.tsv`): marginal and within-window rankings unrelated (ρ −0.18); autism replication 0.85; inflation 3.6-fold.
  - PP3 + PM1 double counting is the clinical consequence.
- **Draft 5 (review 8, Phase 0):**
  - gene-clustered CIs for every R and f (`a12 --part gclust2`);
  - expert-truth PP3 points given PM1: 1.83 [0.84, 3.54], excluding 4;
  - population-frequency sensitivity analysis;
  - wording fixes.
- **A17 pre-registration** (`prereg_fullstudy/A17_PHASE2-5_PREREG.md`) is written. Task 13, the within-window model, was run: **all three tests failed** (supplement Note 7).
- **Still open:** the A17 tasks T5–T12 and T14 (needs an OSF deposit first); verify the Chao 2024 claim; author actions.
- **The X: drive disconnected once.** `wsl --shutdown` after re-plugging fixed it, and git fsck was clean.

## Overnight 7 October 2026 (00:15–02:30 IST): review-7 plan executed (pushed, last commit 86c4f93)

- **Current manuscript: `MANUSCRIPT_v4.md` (draft 4, ~5,000 words incl. tables) + `SUPPLEMENT_v4.md`.** Draft 3 (`MANUSCRIPT_v2.md`) is kept for history.
- **A16 (pre-registered, `prereg_fullstudy/A16_HIERARCHICAL_PREREG.md`, commit 34bac58):** label-free gene → window → AlphaMissense × region model fitted on DD (`fullstudy.py hier`). Validated on temporal ClinVar V1 (VUS 2023-01 → P/LP/B/LB 2026-10, same SCV, truth from the submitter's codes minus PP3/BP4/PM1/PP2/BP1/BP3/PM4/PM5; `fullstudy.py a16val`).
  - V1: 87 P, 887 B, 546 genes. The DDG2P stratum was too small, so the primary analysis used all genes.
  - **H1 refuted:** hierarchical − cap Δ log-loss +0.147 [+0.106, +0.190]. Also worse on eRepo. Label-free vs ClinVar gene-specific calibration: Spearman 0.01 (142 genes).
  - Reported as a pre-registered negative; the paper stays a measurement paper (the cap is the best rule).
- **Not on OSF before the analysis** (you were asleep). Logged in the deviation log and Supplementary Note 1. If you deposit on OSF now, it is a later copy; don't claim it predates the analysis.
- **Other additions:**
  - burial strata (DD R 0.55 buried / 0.36 exposed);
  - gene-clustered bootstrap (DD 0.44 [0.41, 0.47]);
  - MaveDB coverage (64 genes, mostly domain fragments);
  - review-7 response table in `REVIEW_RESPONSE.md`.
- **WSL note:** background `nohup` launches die when the wsl call exits. Run long jobs as tracked foreground `wsl` calls instead. The instance hung once on X: I/O; `wsl --terminate ndd-gate` fixed it and git fsck was clean.
- **Still yours:**
  - verify the Chao 2024 claim marked "[to verify]" in the Discussion;
  - names and affiliations; preprint; Zenodo; journal formatting.

## Done 7 October 2026 (pushed)

- Task 10 pLDDT: DD AlphaMissense R 0.53 (pLDDT>=70) / 0.39 (<70); Results paragraph added.
- Task 1 MAVE (`make_mave_truth.py`): only PPM1D and SHOC2 informative (flag constant elsewhere); PPM1D R 0.72 [0.65, 0.78], SHOC2 0.86 [0.71, 1.04], pooled 0.91 [0.82, 1.02]. Results paragraph + limitation 9.
- Final pass: amendments fifteen (Methods said eleven), abstract retention pair now 0.31 vs 0.51 (same sites), Supp Table 4 listed.

## Remaining, highest impact first

1. **DONE. Task 10, pLDDT strata.** `fullstudy.py a12 --part plddt` was relaunched; log at `work/a12_plddt.log`, output `results_fullstudy/a12_plddt_summary.txt`.
   - When it finishes: add one Results paragraph (R × PM1-like at pLDDT ≥ 70 vs < 70, DD and autism).
   - Then `git add -f` the outputs and commit.
2. **DONE. Task 1, MAVE validation (ProteinGym, 7 DDG2P genes:** DLG4, PITX2, PPM1D, PTEN, RAF1, SHOC2, YAP1).
   - Data: `https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip` (43 MB). The reference list is `work_pg.csv`.
   - Map `mutant` (e.g. A12V) to column-store sites using `tx`, `prot_pos`, `ref_aa` and `alt_aa` (stored as ord codes) for each gene's tx.
   - Truth: `DMS_score_bin` (0 = loss of function). Compute R(AlphaMissense PP3 × PM1-like) and conditional PP3 points, as in `make_erepo_truth.py`.
   - Report as exploratory: few genes.
3. **Task 7, label-free vs ClinVar gene-specific calibration** in label-poor genes. Optional; exploratory.
4. **DONE (7 Oct). Final pass.** Re-read `MANUSCRIPT_v2.md` for consistency:
   - abstract numbers vs Results;
   - amendment count is 15;
   - figure and table list includes Supplementary Table 4.
   Then push.

## Author actions (only the user can do these)

Names and affiliations; OSF pre-registration for any further validation; preprint; Zenodo/container; journal formatting. SPARK has been dropped.
