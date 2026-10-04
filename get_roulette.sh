#!/usr/bin/env bash
# Roulette per-site mutation rates (Seplyarskiy et al. 2023, GRCh38), one chromosome at a time:
# download the VCF + index, keep only the coding sites listed in work/regions/chrN.regions.tsv
# (written by `python gate.py prep-gtf`), then delete the VCF unless KEEP_VCF=1.
#
#   bash get_roulette.sh                     all 22 autosomes (~73 GB downloaded in total, ~7 GB on disk at a time)
#   CHROMS="21 22" bash get_roulette.sh      just two small chromosomes (smoke test, ~1.9 GB)
#   PARALLEL=3 bash get_roulette.sh          three chromosomes at once (needs ~20 GB free disk)
#
# Safe to re-run: finished chromosomes are skipped and interrupted downloads resume.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATA="${DATA:-$HERE/data}"
export WORK="${WORK:-$HERE/work}"
export BASE="${ROULETTE_URL:-https://genetics.bwh.harvard.edu/downloads/Vova/Roulette}"
export KEEP_VCF="${KEEP_VCF:-0}"
CHROMS="${CHROMS:-$(seq 1 22)}"
PARALLEL="${PARALLEL:-1}"

export CURL_Q="-sS"; [ "$PARALLEL" = "1" ] && CURL_Q="--progress-bar"   # progress bar only when not parallel
command -v bcftools >/dev/null || { echo "bcftools not found - run: conda activate ndd-gate"; exit 1; }
mkdir -p "$DATA/roulette" "$WORK/roulette"

fetch() {  # fetch URL DEST  - resumable; DEST.ok marks a finished download
  local url="$1" dest="$2"
  if [ -f "$dest.ok" ]; then return 0; fi
  if ! curl -fL --retry 20 --retry-delay 10 --retry-all-errors -C - $CURL_Q -o "$dest" "$url"; then
    rm -f "$dest"
    curl -fL --retry 20 --retry-delay 10 --retry-all-errors $CURL_Q -o "$dest" "$url"
  fi
  touch "$dest.ok"
}

one_chrom() {
  set -euo pipefail
  local N="$1"
  local out="$WORK/roulette/chr$N.tsv.gz"
  if [ -f "$out.ok" ]; then echo "chr$N: already extracted"; return 0; fi
  local reg="$WORK/regions/chr$N.regions.tsv"
  if [ ! -s "$reg" ]; then echo "chr$N: $reg missing - run 'python gate.py prep-gtf' first"; return 1; fi
  local f="$DATA/roulette/${N}_rate_v5.2_TFBS_correction_all.vcf.bgz"
  echo "chr$N: downloading $(basename "$f")"
  fetch "$BASE/$(basename "$f")" "$f"
  fetch "$BASE/$(basename "$f").csi" "$f.csi"
  echo "chr$N: extracting coding sites"
  if ! bcftools query -R "$reg" \
    -f '%POS\t%REF\t%ALT\t%FILTER\t%INFO/PN\t%INFO/MR\t%INFO/AR\t%INFO/MG\t%INFO/MC\n' "$f" | gzip -1 > "$out.tmp"; then
    echo "chr$N: bcftools failed - the download may be corrupt; its .ok markers are removed so a re-run fetches it again"
    rm -f "$f.ok" "$f.csi.ok"; return 1
  fi
  local n last_pos last_reg
  n=$(gzip -dc "$out.tmp" | wc -l | tr -d ' ')
  if [ "$n" -lt 1000 ]; then echo "chr$N: extraction gave only $n lines - something is wrong"; return 1; fi
  last_pos=$(gzip -dc "$out.tmp" | tail -n 1 | cut -f1)
  last_reg=$(tail -n 1 "$reg" | cut -f2)
  if [ "$last_pos" -lt "$last_reg" ]; then   # truncated input: the last coding region was never reached
    echo "chr$N: extraction stops at $last_pos, before the last region ($last_reg) - removing the download's .ok markers"
    rm -f "$f.ok" "$f.csi.ok"; return 1
  fi
  mv "$out.tmp" "$out"
  touch "$out.ok"
  echo "chr$N: $n coding-site records kept"
  if [ "$KEEP_VCF" != "1" ]; then rm -f "$f" "$f.ok" "$f.csi" "$f.csi.ok"; fi
}
export -f fetch one_chrom

printf '%s\n' $CHROMS | xargs -P "$PARALLEL" -I{} bash -c 'one_chrom "$1"' _ {}
echo "Roulette extraction finished for: $(echo $CHROMS)"
