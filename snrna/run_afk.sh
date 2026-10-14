#!/usr/bin/env bash
# Unattended wrap-up on the DGX (~5-6 h). Run in tmux AFTER final_eval.py has finished:
#   tmux new -s afk 'bash run_afk.sh'
# 1. release.py --stage final  (scores + calibrated PP3/BP4 thresholds -> release/*.zip for Zenodo)
# 2. v4_ablation.py            (post hoc: which features carry the blind RNU2-2 result)
# 3. run_v3.sh                 (exploratory: the v2 graph model retrained with gnomAD + UK Biobank counts, ~4-5 h)
# 4. packs everything into snrna_afk_<date>.tar.gz
set -uo pipefail
cd "$(dirname "$0")"
LOG="afk_$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee -a "$LOG") 2>&1
step() { echo; echo "=== $(date -u +%H:%M:%S) $*"; }
pgrep -f final_eval.py > /dev/null && { echo "final_eval.py is still running; start this after it finishes"; exit 1; }
[[ -f sealed/blind_labels.tsv ]] || echo "WARNING: sealed/blind_labels.tsv missing (ablation will be skipped)"

step "1/3 release.py --stage final"; python release.py --stage final || echo "!!! release failed"
step "2/3 v4_ablation.py"
[[ -f sealed/blind_labels.tsv ]] && { python v4_ablation.py --blind sealed/blind_labels.tsv || echo "!!! ablation failed"; }
step "3/3 run_v3.sh"
[[ -f model_inputs_v3/v3_counts.json ]] || python v3_counts.py cohorts/ukb/*.csv > /dev/null
bash run_v3.sh || echo "!!! run_v3 failed"
step "pack"
tar czf "snrna_afk_$(date -u +%Y%m%dT%H%M%SZ).tar.gz" release/*.zip results/v4_ablation.tsv results/final_* results/calibration* \
  results_v3 "$LOG" 2>/dev/null
echo "AFK DONE. Paste: cat results/v4_ablation.tsv; ls release/; tail -40 results_v3/model_report.txt"
