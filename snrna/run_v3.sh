#!/usr/bin/env bash
# v3 (exploratory follow-up, NOT the pre-registered v2): retrain with extra cohort counts, reusing the v2 Evo 2 features.
#   python v3_counts.py cohorts/ukb.tsv [...]    # first: build model_inputs_v3/
#   tmux new -s v3 'bash run_v3.sh; bash'        # ~3-4 h on the A100; never touches model_runs/ or results/
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
[[ -f model_inputs_v3/v3_counts.json ]] || { echo "run v3_counts.py first"; exit 1; }
[[ -f evo2_run/features/deltas.tsv ]] || { echo "Evo 2 features missing (they come from the v2 run)"; exit 1; }
mkdir -p logs; LOG="logs/run_v3_$(date -u +%Y%m%dT%H%M%SZ).log"; exec > >(tee -a "$LOG") 2>&1
python train.py --inputs model_inputs_v3 --out model_runs_v3 --stage all
python evaluate.py --runs model_runs_v3 --inputs model_inputs_v3 --results results_v3
tar czf "snrna_v3_results_$(date -u +%Y%m%dT%H%M%SZ).tar.gz" results_v3 model_runs_v3/*.json model_runs_v3/preds/*.json model_inputs_v3/v3_counts.json "$LOG"
echo "DONE (v3)"
