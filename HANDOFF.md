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

## Remaining, highest impact first

1. **Task 10, pLDDT strata.** `fullstudy.py a12 --part plddt` was relaunched; log at `work/a12_plddt.log`, output `results_fullstudy/a12_plddt_summary.txt`.
   - When it finishes: add one Results paragraph (R × PM1-like at pLDDT ≥ 70 vs < 70, DD and autism).
   - Then `git add -f` the outputs and commit.
2. **Task 1, MAVE validation (ProteinGym, 7 DDG2P genes:** DLG4, PITX2, PPM1D, PTEN, RAF1, SHOC2, YAP1).
   - Data: `https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip` (43 MB). The reference list is `work_pg.csv`.
   - Map `mutant` (e.g. A12V) to column-store sites using `tx`, `prot_pos`, `ref_aa` and `alt_aa` (stored as ord codes) for each gene's tx.
   - Truth: `DMS_score_bin` (0 = loss of function). Compute R(AlphaMissense PP3 × PM1-like) and conditional PP3 points, as in `make_erepo_truth.py`.
   - Report as exploratory: few genes.
3. **Task 7, label-free vs ClinVar gene-specific calibration** in label-poor genes. Optional; exploratory.
4. **Final pass.** Re-read `MANUSCRIPT_v2.md` for consistency:
   - abstract numbers vs Results;
   - amendment count is 15;
   - figure and table list includes Supplementary Table 4.
   Then push.

## Author actions (only the user can do these)

Names and affiliations; OSF pre-registration for any further validation; preprint; Zenodo/container; journal formatting. SPARK has been dropped.
