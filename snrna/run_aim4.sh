#!/usr/bin/env bash
# A100/DGX: everything after the official v2 evaluation that does NOT touch the sealed labels.
#   1. comparators_report.py   (A3, reported only)
#   2. calibrate.py            (Aim 4, non-blind genes; preliminary until final_eval.py)
#   3. release.py --stage scores (scores only, starts the pseudo-prospective clock; upload the zip to Zenodo)
#   4. gene_discovery.py       (gnomAD, then gnomAD + UK Biobank) + robustness checks
#   5. geuvadis_carriers.py    (downloads Geuvadis junctions + hg19.2bit, ~1.3 GB; 1000G genotypes streamed)
# Run inside tmux, in snrna/, conda env evo2:   tmux new -s aim4 'bash run_aim4.sh'
# Needs: python evaluate.py already run (results/model_scores.tsv.gz). Never decrypts sealed/.
set -uo pipefail
cd "$(dirname "$0")"
LOG="aim4_$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee -a "$LOG") 2>&1
step() { echo; echo "=== $(date -u +%H:%M:%S) $*"; }
ok() { echo "--- done: $1"; }
fail() { echo "!!! FAILED: $1 (continuing)"; }

[[ -f results/model_scores.tsv.gz ]] || { echo "results/model_scores.tsv.gz missing: run python evaluate.py first"; exit 1; }
python -c "import pysam, py2bit" 2>/dev/null || pip install -q pysam py2bit
python -c "import scipy" 2>/dev/null || pip install -q scipy

step "1/5 comparators_report.py"; python comparators_report.py && ok comparators || fail comparators
step "2/5 calibrate.py (non-blind)"; python calibrate.py && ok calibrate || fail calibrate
step "3/5 release.py --stage scores"; python release.py --stage scores && ok release || fail release

step "4/5 gene_discovery.py"
python gene_discovery.py && ok "gene_discovery (gnomAD)" || fail gene_discovery
if [[ ! -f model_inputs_v3/variants.tsv.gz ]]; then python v3_counts.py cohorts/ukb/*.csv > /dev/null || fail v3_counts; fi
if [[ -f model_inputs_v3/variants.tsv.gz ]]; then
  python gene_discovery.py --inputs model_inputs_v3 --out results/gene_discovery_v3 && ok "gene_discovery (+UKB)" || fail gene_discovery_v3
  python gene_discovery_robustness.py > results/gene_discovery_robustness.txt 2>/dev/null && ok robustness || fail robustness
fi

step "5/5 geuvadis_carriers.py"
G=geuv_data; mkdir -p "$G"
B="https://ftp.ebi.ac.uk/biostudies/fire/E-GEUV-/001/E-GEUV-1/Files/E-GEUV-1/analysis_results"
[[ -s $G/junc.txt.gz ]] || curl -fsS --retry 4 -o $G/junc.txt.gz "$B/GD462.JunctionQuantCount.45N.50FN.samplename.resk10.txt.gz" || fail "download junctions"
[[ -s $G/hg19.2bit ]] || curl -fsS --retry 4 -o $G/hg19.2bit https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit || fail "download hg19.2bit"
python geuvadis_carriers.py --junctions $G/junc.txt.gz --twobit $G/hg19.2bit && ok geuvadis || fail geuvadis

step "pack"
tar czf "snrna_aim4_$(date -u +%Y%m%dT%H%M%SZ).tar.gz" results/comparators_* results/calibration* results/gene_discovery* \
  results/geuvadis_carriers* release/*.zip "$LOG" 2>/dev/null
echo "AIM4 DONE. Paste: cat results/comparators_report.txt results/calibration*.txt results/geuvadis_carriers.json | head -120"
