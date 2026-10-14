# Handoff (14 Oct 2026): where everything is and what comes next

All project state is in this folder (branch `claude/snrna-vep-v2-framework-b0abng`). Read in this order:
`PLAN.md` (idea) → `PREREG.md` (frozen plan, deviation log s14, amendments A1-A4 in s15) → `README.md` (status, runbooks) → `NOVELTY.md` → `paper/methods.md`.

## Done
- G1 PASSES (Evo 2 not enough), G2 PARTIAL (selection signal real). Model built, frozen (A2/A2b, timestamped), trained on the A100.
- Blind labels (RNU2-2, RNU5B-1) sealed on the A100 (`sealed/`); passphrase held by the author only.
- Comparators, calibration, final evaluation, release and v3 code written and tested (A3).
- UK Biobank exports for v3 in `cohorts/ukb/` (in progress; see `ukb_regions.txt`).

## Next steps, in order
1. A100 run finishes -> `mv results results_autorun_unopened; sha256sum results_autorun_unopened/* > autorun_unopened.sha256` (do NOT open).
2. Literature tally (prompt in chat history; output TSV) -> `python curate_literature.py tally.tsv` -> commit `labels/L2_curated_*`.
3. On the A100: `git pull; git checkout -- results; python evaluate.py` = official v2 result. Then `comparators_report.py`, `calibrate.py`, `release.py --stage scores` (Zenodo).
4. A4 analyses (PREREG s15 A4).
5. Final blind test: `SNRNA_SEAL_PASS=... python final_eval.py --reason "final evaluation"`, then `release.py --stage final`.
6. v3: finish UK Biobank exports -> combine -> `python v3_counts.py ukb.csv` -> `bash run_v3.sh`.
7. Write the paper (target GIM; AJHG if Aim 2 + blind test strong).
