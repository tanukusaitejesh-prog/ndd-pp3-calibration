#!/usr/bin/env bash
# Self-test on synthetic data with known answers (about 4 minutes). Run this before the big downloads.
# It runs the real pipeline (including bcftools if installed) on two synthetic chromosomes and checks
# that it recovers the true coordinates, counts and likelihood ratios, and that the null check
# fails when detection is deliberately made to depend on protein context.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
PY="${PYTHON:-python}"

run_case() {
  local T="$1" FLAG="$2"
  rm -rf "$T"
  "$PY" tests/make_synthetic.py "$T" $FLAG
  export DATA="$T/data" WORK="$T/work" RESULTS="$T/results"
  "$PY" gate.py prep-dnv > "$T/log.txt"
  "$PY" gate.py prep-am --chunksize 200000 >> "$T/log.txt"
  "$PY" gate.py prep-gtf >> "$T/log.txt"
  if command -v bcftools >/dev/null && command -v bgzip >/dev/null; then
    for N in 21 22; do
      f="$DATA/roulette/${N}_rate_v5.2_TFBS_correction_all.vcf.bgz"
      { cat "$T/roulette_vcf/header.txt"
        gzip -dc "$T/roulette_full/chr$N.tsv.gz" | awk -F'\t' -v OFS='\t' -v C="$N" \
          '{info="PN="$5; if ($6 != ".") info=info";MR="$6; print C,$1,".",$2,$3,".",$4,info";AR="$7";MG="$8";MC="$9}'
      } | bgzip -c > "$f"
      bcftools index -c -f "$f"
      touch "$f.ok" "$f.csi.ok"
    done
    CHROMS="21 22" KEEP_VCF=1 bash get_roulette.sh >> "$T/log.txt"
    echo "  (Roulette extraction ran through bcftools)"
  else
    "$PY" tests/fake_extract.py "$T" >> "$T/log.txt"
    echo "  (bcftools not found: used the python stand-in for the extraction step)"
  fi
  "$PY" gate.py check-inputs >> "$T/log.txt"
  "$PY" gate.py sites --chroms 21-22 >> "$T/log.txt"
  "$PY" gate.py run --allow-partial --n-boot 400 > "$T/report_stdout.txt"
  "$PY" tests/check_synthetic.py "$T" $FLAG $( [ "$FLAG" = "--cpg-factor 0.75" ] && echo --cpg )
}

echo "== clean scenario: the pipeline must recover the truth and say GO"
run_case "$HERE/selftest_out/clean" ""
echo
echo "== biased scenario: de novos are under-detected in high-AlphaMissense codons; the gate must say KILL"
run_case "$HERE/selftest_out/biased" "--biased"
echo
echo "== CpG scenario: true CpG transition rate is 0.75x the model; pooled scaling fails, class scaling must rescue it"
run_case "$HERE/selftest_out/cpg" "--cpg-factor 0.75"
echo
echo "ALL SELF-TESTS PASSED - the pipeline is working on this computer."
