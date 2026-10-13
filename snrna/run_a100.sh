#!/usr/bin/env bash
# snRNA-VEP v2 on one GPU: seal -> Evo 2 features -> train -> evaluate -> pack. Resumable: re-run the same command
# after any interruption and finished work is skipped.
#
#   conda activate evo2
#   bash run_a100.sh --dry-run    # ~20-30 min. Real Evo 2 on one gene, GPU timing benchmark, quick training on
#                                 # SYNTHETIC counts, evaluation. Fits no real count and needs no passphrase.
#   tmux new -s snrna 'bash run_a100.sh; bash'   # full pre-registered run (amendment A2), unattended after the
#                                 # sealing passphrase (asked once at the start; or export SNRNA_SEAL_PASS first)
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
MODE=full
[[ "${1:-}" == "--dry-run" ]] && MODE=dry
# SNRNA_TEST_CPU=1: development check of this script without a GPU (mock Evo 2, CPU training; dry run only)
CPU_TEST=${SNRNA_TEST_CPU:-}
DEV=()
if [[ -n $CPU_TEST ]]; then
  [[ $MODE == dry ]] || { echo "SNRNA_TEST_CPU is for --dry-run only"; exit 1; }
  export SNRNA_EVO2_MOCK=1; DEV=(--device cpu)
fi
STAMP=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p logs
LOG="logs/run_${MODE}_${STAMP}.log"
exec > >(tee -a "$LOG") 2>&1
trap 'echo; echo "!!! stopped at line $LINENO (exit $?). Fix the cause and re-run the same command: finished work is kept."' ERR

stage() { echo; echo "=== $(date -u +%H:%M:%S) UTC  $*"; }

stage "0. environment ($MODE run, log $LOG)"
rc=0
python - <<'PY' || rc=$?
import importlib.util, os, sys
cpu_test = bool(os.environ.get("SNRNA_TEST_CPU"))
missing = [m for m in ["numpy", "pandas", "torch"] + ([] if cpu_test else ["evo2"]) if importlib.util.find_spec(m) is None]
if missing:
    sys.exit("MISSING python packages: " + " ".join(missing) + " (activate the evo2 environment)")
extra = [m for m in ["scipy", "pytest"] if importlib.util.find_spec(m) is None]
if extra:
    print("missing (will install):", " ".join(extra))
    sys.exit(3)
import torch
print("torch", torch.__version__, "| CUDA", torch.cuda.is_available(), "|", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "no GPU")
free, total = torch.cuda.mem_get_info() if torch.cuda.is_available() else (0, 0)
print(f"GPU memory free {free / 2**30:.1f} of {total / 2**30:.1f} GB")
if not torch.cuda.is_available() and not cpu_test:
    sys.exit("no CUDA device: activate the evo2 environment on the GPU machine")
if free < 30 * 2**30 and not cpu_test:
    print("WARNING: less than 30 GB free on the GPU; another process is using it (batches will shrink automatically)")
PY
if [[ $rc -eq 3 ]]; then
  echo "installing scipy and pytest (the only extra dependencies)"; pip install -q scipy pytest
elif [[ $rc -ne 0 ]]; then
  exit $rc
fi
[[ -n $CPU_TEST ]] || python -c "import torch; assert torch.cuda.is_available(), 'no CUDA device'"
command -v openssl >/dev/null || { echo "openssl not found"; exit 1; }
echo "git commit: $(git rev-parse HEAD 2>/dev/null || echo unknown)"
if [[ -n "$(git status --porcelain --untracked-files=no 2>/dev/null)" ]]; then
  echo "WARNING: tracked files differ from the commit:"; git status --short --untracked-files=no
fi

# the frozen code must match the latest A2 manifest (full run); A2b = logging-only patch of A2
MANIFEST=amendment_A2_manifest.txt
[[ -f amendment_A2b_manifest.txt ]] && MANIFEST=amendment_A2b_manifest.txt
if [[ -f $MANIFEST ]]; then
  # code and inputs only (.md files may gain later log entries, e.g. the sealing hashes)
  if grep -E '^[0-9a-f]{64}  ' "$MANIFEST" | grep -v '\.md$' | sha256sum --quiet -c -; then
    echo "code and inputs match $MANIFEST"
  elif [[ $MODE == full ]]; then
    echo "code or inputs differ from $MANIFEST: git pull / git checkout the frozen commit"; exit 1
  else
    echo "WARNING: code differs from $MANIFEST (allowed for a dry run)"
  fi
elif [[ $MODE == full ]]; then
  echo "amendment_A2_manifest.txt missing: git pull first"; exit 1
fi

stage "0b. fast tests (CPU)"
python -m pytest -q -x tests/test_labels.py tests/test_evo2_features.py tests/test_homology.py
if [[ $MODE == dry && -z $CPU_TEST ]]; then python -m pytest -q -x tests/test_model.py; fi

if [[ $MODE == full ]]; then
  stage "1. seal the blind labels (PREREG s5, A2.10)"
  if [[ -f sealed/SEALED_MANIFEST.jsonl ]] && grep -q '"event": "SEALED"' sealed/SEALED_MANIFEST.jsonl; then
    echo "already sealed:"; grep '"event": "SEALED"' sealed/SEALED_MANIFEST.jsonl | tail -1
  else
    if [[ -z "${SNRNA_SEAL_PASS:-}" ]]; then
      [[ -t 0 ]] || { echo "no terminal for the passphrase: export SNRNA_SEAL_PASS=... and re-run"; exit 1; }
      read -r -s -p "Sealing passphrase (keep it; the analyst must not know it): " p1; echo
      read -r -s -p "Repeat: " p2; echo
      [[ "$p1" == "$p2" && -n "$p1" ]] || { echo "passphrases differ or empty"; exit 1; }
      export SNRNA_SEAL_PASS="$p1"; unset p1 p2
    fi
    python labels.py seal
    unset SNRNA_SEAL_PASS
  fi
  if [[ ! -f sealed/SEALED_MANIFEST.jsonl.ots ]]; then
    command -v ots >/dev/null || timeout 600 pip install -q opentimestamps-client || true
    if command -v ots >/dev/null && timeout 300 ots stamp sealed/SEALED_MANIFEST.jsonl; then
      echo "manifest timestamped: sealed/SEALED_MANIFEST.jsonl.ots"
    else
      echo "WARNING: OpenTimestamps unavailable; the manifest's own UTC time and SHA-256 stand (stamped when results are committed)"
    fi
  fi
  sha256sum sealed/SEALED_MANIFEST.jsonl sealed/blind_labels.tsv.enc
fi

stage "2. Evo 2 features"
if [[ $MODE == dry ]]; then
  python evo2_features.py --out evo2_run/features_dryrun --limit-genes 1
else
  python evo2_features.py
  python - <<'PY'
import os, pandas as pd
v = pd.read_csv("model_inputs/variants.tsv.gz", sep="\t", usecols=["ID"])
d = pd.read_csv("evo2_run/features/deltas.tsv", sep="\t")
g = pd.read_csv("model_inputs/genes.tsv.gz", sep="\t")
emb = [f for f in os.listdir("evo2_run/features/emb") if f.endswith(".npy")]
assert set(v.ID) <= set(d.ID), f"deltas incomplete: {len(set(v.ID) - set(d.ID))} missing"
assert len(emb) == len(g), f"embeddings for {len(emb)} of {len(g)} genes"
print(f"Evo 2 features complete: {d.ID.nunique()} variants, {len(emb)} genes")
PY
fi

stage "3. train"
if [[ $MODE == dry ]]; then
  [[ -n $CPU_TEST ]] || python train.py --benchmark
  python train.py --synthetic-counts 1 --quick --max-lofo 3 --evo2 evo2_run/features_dryrun ${DEV[@]+"${DEV[@]}"}
else
  python train.py --stage all
fi

stage "4. evaluate"
if [[ $MODE == dry ]]; then
  python evaluate.py --runs model_runs_synthetic --results results_dryrun --evo2 evo2_run/features_dryrun --boot 200
else
  python evaluate.py
fi

stage "5. pack"
if [[ $MODE == dry ]]; then
  OUT="snrna_dryrun_${STAMP}.tar.gz"
  tar czf "$OUT" results_dryrun model_runs_synthetic/*.json model_runs_synthetic/preds/*.json \
      evo2_run/features_dryrun/run_info.json "$LOG"
else
  OUT="snrna_a100_results_${STAMP}.tar.gz"
  tar czf "$OUT" results/model_* model_runs/*.json model_runs/preds/*.json evo2_run/features/deltas.tsv \
      evo2_run/features/run_info.json sealed/SEALED_MANIFEST.jsonl* sealed/blind_labels.tsv.enc logs
fi
ls -lh "$OUT"
echo
echo "DONE ($MODE). Upload $OUT to the chat (or paste results*/model_report.txt). Per-seed predictions and Evo 2"
echo "embeddings stay in model_runs*/preds and evo2_run/features/emb on this machine."
