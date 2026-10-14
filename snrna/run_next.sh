#!/usr/bin/env bash
# DGX, after `git pull` in ~/snrna-run/snrna (conda env evo2; needs internet for ClinVar/UCSC/gnomAD/CADD/Rfam).
#   bash run_next.sh part1     # post hoc checks + v4 calibration + VUS + final release + A7 label-free scoring (~1 h)
#   (then: git add a7/ results/; git commit; git push; ots stamp a7/a7_manifest.txt)
#   bash run_next.sh part2     # A7 evaluation (refuses unless a7/ is committed)
# Each step logs to logs/next_<step>.log and a failure does not stop the others.
set -u
cd "$(dirname "$0")"
mkdir -p logs results
B=sealed/blind_labels.tsv
step() { local name=$1; shift; echo "== $name"; if "$@" > "logs/next_$name.log" 2>&1; then tail -n 25 "logs/next_$name.log"; else echo "!! $name FAILED (see logs/next_$name.log)"; tail -n 15 "logs/next_$name.log"; fi; }

case "${1:-}" in
  part1)
    step circularity python circularity_check.py --blind $B --fetch
    step robustness  python blind_robustness.py --blind $B
    step v4_calib    python v4_calibrate.py --blind $B
    step release     python release.py --stage final
    step a7_score    python a7_ncrna.py score
    echo; echo "NEXT: git add a7/ results/ && git commit -m 'A7 scores frozen; post hoc checks' && git push"
    echo "      ots stamp a7/a7_manifest.txt && git add a7/a7_manifest.txt.ots && git commit -m 'A7 timestamp' && git push"
    echo "      then: bash run_next.sh part2"
    ;;
  part2)
    step a7_eval python a7_ncrna.py eval --acmg
    echo; echo "NEXT: git add results/a7_* && git commit -m 'A7 evaluation' && git push"
    ;;
  *) echo "usage: bash run_next.sh part1|part2"; exit 1 ;;
esac
