# Reproducing draft 17 (pre-registered laboratory analyses and cost)

Environment: `environment.yml` / `environment.lock.yml`. De novo analyses need the column store built by `run_all.sh`
(Roulette, AlphaMissense, dbNSFP; large). The draft-14 laboratory analyses need only public ClinVar files.

## One command for the draft-14 laboratory analyses
```
bash reproduce_v14.sh
```
Downloads ClinVar submission_summary (2020-01, 2022-01, 2024-01, 2026-10) and variant_summary (2026-10), then runs:
`make_lab_dataset.py` → `make_rules_cv.py` (H1, pre-registered and `--acmg`) → `make_pp2_confirm.py` (H2) →
`make_genecal_pp3.py` (H3) → `make_vintage_replication.py` (H4) → `make_validity_v14.py` → `make_cost_headline.py`.

## Number / table / figure → script → output

| Item in MANUSCRIPT_v17 | Script | Output (results_fullstudy/) |
|---|---|---|
| Table 1 rows 1–5 (R, de novo) | `fullstudy.py a10` / `a12 --part gclust2` / `a13` | `a10_overlap.tsv`, `a12_gclust2_summary.txt`, `a13_summary.txt`, `pm1pp2_summary.txt` |
| Table 1 row 6 (R3) | `fullstudy.py` (triple stacking) | `triple_summary.txt` |
| Table 2, Figure 3 (cost) | `make_cost_headline.py` | `cost_headline_summary.txt`, `figures/fig_cost_calibration.png` |
| 7,101 VUS | `fullstudy.py a10_post` | `a10_post_summary.txt` |
| Section 3 (H1) | `make_rules_cv.py [--acmg]` | `a19_h1_rules_2026-10[_acmg]_summary.txt` |
| Section 4 (H2) | `make_pp2_confirm.py` | `a19_h2_pp2_2026-10_summary.txt` |
| Section 5 (H4, negative controls) | `make_vintage_replication.py` | `a19_h4_summary.txt`, `a19_h4_eras.tsv` |
| Section 6 (H3) | `make_genecal_pp3.py` | `a19_h3_summary.txt` |
| Section 7, Table 3 | `make_lab_stacking.py`, `make_validity_v14.py`, `make_rules_cv.py` | `lab_stacking_summary.txt`, `validity_v14_summary.txt` |
| Figure 2 | `fig_manuscript.py` | `fig2_dependence.png` |
| Parser annotation sheet | `make_validity_v14.py` | `parser_annotation_sheet_v14.tsv` |

dbNSFP licence: dbNSFP scores are used under its academic licence and are not redistributed; scripts re-derive them.
