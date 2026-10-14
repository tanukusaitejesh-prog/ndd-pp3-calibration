#!/usr/bin/env bash
# A100/DGX: build v4 (amendment A5), print its internal validation and the freeze manifest. No sealed labels touched.
# In snrna/, conda env evo2:   bash run_v4.sh     (a few minutes, CPU)
# Then paste results_v4/v4_manifest.txt and results_v4/v4_report.txt into the chat; the manifest is committed and
# timestamped (A5) BEFORE final_eval.py is run. Do not re-run v4.py after that (the scores would no longer match).
set -euo pipefail
cd "$(dirname "$0")"
LOG="v4_$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee -a "$LOG") 2>&1
[[ -f results_v4/v4_manifest.txt ]] && { echo "results_v4/v4_manifest.txt exists: v4 already frozen, not re-running"; exit 1; }
python -c "import sklearn" 2>/dev/null || pip install -q scikit-learn
[[ -f model_inputs_v3/variants.tsv.gz ]] || python v3_counts.py cohorts/ukb/*.csv > /dev/null
[[ -f evo2_run/features/deltas.tsv ]] || echo "WARNING: evo2_run/features/deltas.tsv missing (the evo2 feature will be neutral)"
python v4.py
echo; echo "=== paste these two into the chat:"; cat results_v4/v4_manifest.txt; echo; cat results_v4/v4_report.txt
