#!/usr/bin/env bash
# DGX: put every DGX-only output that the paper's numbers rest on into git, then push (or pack a tarball if push fails).
# Never adds decrypted labels (sealed/*.tsv). Files over 95 MB are listed for Zenodo instead (GitHub's limit is 100 MB).
#   cd ~/snrna-run/snrna && bash sync_dgx_outputs.sh            # dry list first: bash sync_dgx_outputs.sh --list
set -u
cd "$(dirname "$0")"
LIST=0; [ "${1:-}" = "--list" ] && LIST=1
CANDIDATES=(
  results_v4/v4_manifest.txt results_v4/v4_report.txt results_v4/v4_model.json results_v4/v4_scores.tsv.gz results_v4/v4_cv.tsv
  sealed/blind_labels.tsv.enc sealed/SEALED_MANIFEST.jsonl sealed/SEALED_MANIFEST.jsonl.ots
  a7/a7_manifest.txt a7/a7_manifest.txt.ots a7/a7_scores.tsv.gz a7/a7_score_info.json
  a9/window_snvs.tsv.gz a9/fetch_info.json a9_phylop/window_snvs.tsv.gz
  a10/missense_window_snvs.tsv.gz a10/missense_fetch_info.json a10/spliceai_scores.tsv
  a11/positives_GRCh38.tsv a11/negatives_GRCh38.tsv a11/window_snvs.tsv.gz a11/fetch_info.json
  a13/*/window_snvs.tsv.gz a13/*/fetch_info.json a14/scores.tsv.gz a14/fetch_info.json
  a16/unlabelled.tsv.gz a16/fetch_info.json a10/spliceai_out.vcf paper/tables.md results/*.tsv.gz
  evo2_run/features/deltas.tsv evo2_run/features/run_info.json
  model_runs/*.json model_runs/preds/*.json
  model_inputs_v3/v3_counts.json
  results/*.txt results/*.tsv results/*.tsv.gz results/*.json results/a9_phylop/* results/a8_a6/*
  paper/figures/*.png paper/figures/*.pdf
  release/*.zip release/lookup.html
)
ADD=(); BIG=()
for pat in "${CANDIDATES[@]}"; do
  for f in $pat; do
    [ -f "$f" ] || continue
    case "$f" in sealed/*.tsv) echo "SKIP decrypted labels: $f"; continue ;; esac
    sz=$(stat -c %s "$f")
    if [ "$sz" -gt 95000000 ]; then BIG+=("$f ($((sz / 1000000)) MB)"); else ADD+=("$f"); fi
  done
done
echo "== ${#ADD[@]} files to add; ${#BIG[@]} too large for GitHub (put these on Zenodo):"
for b in "${BIG[@]}"; do echo "   $b"; done
echo "== full SHA-256 of the frozen v4 outputs (log them in PREREG s14; the A5 manifest holds 62-character suffixes):"
sha256sum results_v4/v4_scores.tsv.gz results_v4/v4_model.json 2>/dev/null
grep -E "v4_scores|v4_model" amendment_A5_manifest.txt
if [ "$LIST" = 1 ]; then printf '   %s\n' "${ADD[@]}"; exit 0; fi
git add -f -- "${ADD[@]}"
git diff --cached --name-only | grep -E 'sealed/[^/]*\.tsv$' && { echo "decrypted labels staged: abort"; git reset -q; exit 1; }
git commit -m "DGX outputs: frozen v4, sealed-label record, A7-A14 tables, model runs, figures (hash-checked files unchanged)" || echo "nothing new to commit"
if git push; then
  echo "PUSHED"
else
  T="dgx_outputs_$(date -u +%Y%m%dT%H%M%SZ).tar.gz"
  tar czf "$T" "${ADD[@]}"
  echo "push failed: set a GitHub token (see below) or copy $T to a machine that can push and run there:"
  echo "   tar xzf $T -C <repo>/snrna && cd <repo> && git add -f snrna && git commit -m 'DGX outputs' && git push"
  echo "token: GitHub > Settings > Developer settings > Fine-grained token (this repo, Contents: read/write), then"
  echo "   git config credential.helper store && git push    # username = your GitHub name, password = the token"
fi
