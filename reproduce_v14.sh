#!/usr/bin/env bash
# Draft-14 laboratory analyses end to end (public ClinVar only). Safe to re-run: finished downloads are skipped.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p work/cvt
B=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive
get() { [ -s "work/cvt/$2" ] || curl -sSfo "work/cvt/$2" "$B/$1"; }
get submission_summary_2026-10.txt.gz submission_summary_2026-10.txt.gz
get variant_summary_2026-10.txt.gz variant_summary_2026-10.txt.gz
for v in 2024-01 2022-01 2020-01; do get "${v%%-*}/submission_summary_$v.txt.gz" "submission_summary_$v.txt.gz"; done
python make_lab_dataset.py 2026-10 2024-01 2022-01 2020-01
python make_rules_cv.py 2026-10
python make_rules_cv.py 2026-10 --acmg
python make_pp2_confirm.py 2026-10
python make_genecal_pp3.py
python make_vintage_replication.py
python make_validity_v14.py
python make_cost_headline.py
python make_denovo_rules.py
