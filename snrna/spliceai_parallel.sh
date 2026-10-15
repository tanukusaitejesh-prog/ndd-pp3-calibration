#!/usr/bin/env bash
# Run SpliceAI on a10/spliceai_in.vcf in N parallel chunks and merge to a10/spliceai_out.vcf (same variants, same
# settings as run_next.sh part5: GRCh38 annotation, distance 50). Operational only; the A10 analysis is unchanged.
#   bash spliceai_parallel.sh gpu      # TensorFlow with CUDA in env spliceai-gpu, 4 chunks sharing the GPU(s)
#   bash spliceai_parallel.sh cpu 16   # CPU env spliceai, 16 chunks, one thread each
set -u
cd "$(dirname "$0")"
MODE=${1:-cpu}; N=${2:-16}
IN=a10/spliceai_in.vcf; WORK=a10/spliceai_chunks; mkdir -p "$WORK" logs
if [ "$MODE" = gpu ]; then
  ENV=spliceai-gpu; N=${2:-4}
  if ! conda env list | grep -q "^spliceai-gpu "; then
    echo "== creating conda env spliceai-gpu (python 3.10, tensorflow 2.15 with CUDA)"
    conda create -y -q -n spliceai-gpu -c conda-forge python=3.10 && \
      conda run -n spliceai-gpu pip install -q "tensorflow[and-cuda]==2.15.1" "setuptools<70" spliceai
  fi
  echo "== checking that TensorFlow sees a GPU (first CUDA start can take 1-3 min; gives up after 5 min) $(date +%H:%M:%S)"
  timeout 300 conda run --no-capture-output -n spliceai-gpu python -c "import tensorflow as tf; g=tf.config.list_physical_devices('GPU'); print('GPUs visible to TensorFlow:', len(g)); assert g, 'no GPU'"     || { echo "!! TensorFlow GPU check failed or hung: run 'bash spliceai_parallel.sh cpu 16' instead"; exit 1; }
  THREADS=""
else
  ENV=spliceai; THREADS="OMP_NUM_THREADS=1 TF_NUM_INTRAOP_THREADS=1 TF_NUM_INTEROP_THREADS=1"
fi
rm -f "$WORK"/part_*   # stale chunks from an earlier run (another N) would otherwise be merged twice
grep "^#" "$IN" > "$WORK/header.vcf"
grep -v "^#" "$IN" > "$WORK/body.txt"
TOTAL=$(wc -l < "$WORK/body.txt")
split -n l/"$N" -d -a 2 "$WORK/body.txt" "$WORK/part_"
echo "== $TOTAL variants in $N chunks, mode $MODE ($(date +%H:%M:%S))"
for f in "$WORK"/part_??; do
  cat "$WORK/header.vcf" "$f" > "$f.vcf"
  ( env $THREADS conda run --no-capture-output -n "$ENV" spliceai -I "$f.vcf" -O "$f.out.vcf" -R data/hg38.fa -A grch38 -D 50 \
      > "logs/spliceai_$(basename "$f").log" 2>&1 && echo "   chunk $(basename "$f") done ($(date +%H:%M:%S))" ) &
done
while [ "$(jobs -r | wc -l)" -gt 0 ]; do
  sleep 60
  DONE=$(cat "$WORK"/part_??.out.vcf 2>/dev/null | grep -vc "^#")
  echo "   $(date +%H:%M:%S) scored about $DONE / $TOTAL (written in buffered blocks)"
done
wait
grep "^#" "$WORK/part_00.out.vcf" > a10/spliceai_out.vcf
for f in "$WORK"/part_??.out.vcf; do grep -v "^#" "$f"; done >> a10/spliceai_out.vcf
echo "== merged: $(grep -vc '^#' a10/spliceai_out.vcf) of $TOTAL variants in a10/spliceai_out.vcf ($(date +%H:%M:%S))"
echo "NEXT: python a10_extend.py spliceai_parse && python a10_extend.py eval"
