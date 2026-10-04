#!/usr/bin/env bash
# Steps 6-10 of README.md in one go (after the self-test and the pre-registration commit).
# Safe to re-run: finished downloads, chromosomes and preparation steps are skipped.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
bash get_data.sh
python gate.py check-inputs
[ -f work/dnv_grch38.tsv ] || python gate.py prep-dnv
python gate.py prep-am
[ -f work/gtf/transcripts.pkl ] || python gate.py prep-gtf
PARALLEL="${PARALLEL:-1}" bash get_roulette.sh
python gate.py sites
python gate.py run
