# Handoff (14 Oct 2026, end of session): state and remaining tasks

Branch `claude/snrna-vep-v2-framework-b0abng`, folder `snrna/`. Every decision and result is logged in
`PREREG.md` s14 (deviation log, newest rows at the bottom of the table) and s15 (amendments). Read that first.
DGX: `~/snrna-run/snrna`, conda env `evo2`. The blind labels are decrypted at `sealed/blind_labels.tsv` (single
unsealing done; logged in `sealed/SEALED_MANIFEST.jsonl`).

## Results so far (all logged in PREREG s14)
| Item | Result |
|---|---|
| v2 Aim 2 (pre-registered, rho vs RNU4-2 SGE) | 0.144, FAILED (bar 0.346) |
| v2 G3 (week-6, RNU4ATAC) | delta +0.083, FAILED (bar 0.15) |
| **v2 blind (RNU2-2, 22 P/LP vs 127 L3)** | S 0.784, F 0.851, CADD 0.501; weighted delta +0.267 [+0.134, +0.399] |
| **v4 (A5, frozen before unsealing) blind RNU2-2** | **0.892 [0.817, 0.951] vs CADD 0.501; delta +0.374 [+0.262, +0.481]** |
| v4 internal LOCO (non-blind) | delta vs CADD +0.186 [+0.108, +0.267] |
| v4 vs SGE, RNU4 cluster held out (post hoc) | rho 0.382; Rfam ref_col_freq alone 0.459 |
| v4 ablation on blind RNU2-2 (post hoc) | drop-one 0.875-0.908; Rfam features alone 0.85-0.86 |
| Aim 4 final calibration | PP3 Supporting in RNU2-2 (F >= 0.402), RNU4-2, RNU4ATAC (+BP4), dominant pool; no Moderate at prior 0.1 |
| Gene discovery (exploratory) | within-gene window constraint ranks the 5 NDD genes top (AUROC 1.0, p 5e-5) |
| A6 inheritance prediction | non-blind AUROC 0.98; blind not evaluable (no inheritance annotation) |
| Geuvadis carriers (pre-specified) | null (H1 p 0.29) |
| v3 (+UK Biobank) | Aim 2 rho 0.142 (no gain); RNU4ATAC f 0.646 -> 0.687 |
| RNU5B-1 blind | only 1 P/LP: uninformative |

## Remaining tasks, in order
1. **Release:** on the DGX `python release.py --stage final 2>&1 | tail -10` (if it errors, fix release.py; likely it
   expects results/final_run_info.json, which final_eval did not write because calibrate stopped there; it can be
   written by hand with the date, reason and git commit). Upload `release/*final*.zip` to Zenodo; record the DOI in PREREG s14.
2. **Circularity check (high priority):** re-score blind RNU2-2 with v4, S and CADD using only labels whose
   classification did not rest on in-silico/conservation evidence (use `source`/`review_status` columns of
   sealed/blind_labels.tsv). Write `circularity_check.py` (copy v4_blind_eval.py, add a label filter); report AUROCs.
3. **v4 calibration (post hoc):** run calibrate.py's local-LR method on results_v4/v4_scores.tsv.gz (all labels
   incl. blind) -> PP3/BP4 thresholds for v4. Label "post hoc" in PREREG s14.
4. **ClinVar VUS table:** score all ClinVar snRNA VUS (labels.py fetches ClinVar) with v4 and the calibrated
   thresholds; count reclassifiable ones.
5. **Novelty check:** paste NOVELTY_CHECK_PROMPT.md into a web-search AI; add the line "Also list any published
   SGE/MPRA/functional screens of RNU4ATAC, RNU2-2, RNU5B-1, RNU12, RNU6ATAC, RNU6-x or other snRNAs (data
   availability)". If an untouched screen exists, test frozen v4 on it (third validation).
6. **Timestamps/registry:** `ots upgrade *.ots` (needs the `ots` client), deposit PREREG.md + amendment_*_manifest.txt
   + .ots files on OSF.
7. **Lookup webpage** (static, from the release scores + thresholds).
8. **Figures:** gene-discovery ranking; blind ROC (v4, S, CADD on RNU2-2); calibration (local LR curves); ablation
   bar chart; v2/v3/v4 Aim 2 comparison.
9. **Paper** (draft in paper/methods.md): target GIM; fallback Genome Medicine / npj Genomic Medicine.
   Framing: pre-registered study; the primary selection model failed per variant (v2, v3, informative negative);
   the frozen v4 beat CADD blind on RNU2-2 (0.89 vs 0.50); RNA-family profiles carry the signal; calibrated
   PP3/BP4 for labs; gene discovery and inheritance as secondary findings. State limitations: one blind gene (22
   variants), v4 designed after v2 results (frozen before unsealing), no functional/clinical validation, solo author.
10. Preprint on medRxiv after the draft.

## Rules that still apply
- Never re-run v4.py or edit results_v4/ (hash-frozen, A5). Post hoc analyses are labelled exploratory in PREREG s14.
- Commits: end messages with the Claude-Session line only (no co-author line). No model names in commits.
- Do not publish the seal passphrase.
