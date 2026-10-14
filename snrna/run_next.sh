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
step() { local name=$1; shift; echo; echo "== [$(date +%H:%M:%S)] START $name"; local t0=$SECONDS
  if "$@" 2>&1 | tee "logs/next_$name.log"; [ "${PIPESTATUS[0]}" -eq 0 ]; then echo "== [$(date +%H:%M:%S)] DONE $name ($((SECONDS-t0)) s)"
  else echo "!! [$(date +%H:%M:%S)] FAILED $name ($((SECONDS-t0)) s) - see logs/next_$name.log"; fi; }

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
  part3)  # A9: genome-wide rarity inflation (positives are committed in a9/; needs internet; ~1-2 h)
    pip install -q pyBigWig 2>/dev/null || echo "pyBigWig not installed: phyloP447 will be skipped"
    step a9_fetch python a9_inflation.py fetch --workers 8
    step a9_eval  python a9_inflation.py eval
    echo; echo "NEXT: git add a9/window_snvs.tsv.gz results/a9_* && git commit -m 'A9 evaluation' && git push"
    ;;
  part4)  # post hoc robustness after A8/A9 (~25 min, mostly the UCSC phyloP requests)
    step a8_sens   python extras.py a8_sensitivity --blind $B
    step mechanism python extras.py mechanism --blind $B
    step a9_phylop python extras.py a9_phylop
    step figures   python figures.py
    echo; echo "NEXT: git add results/ a9_phylop/ paper/figures && git commit -m 'Post hoc robustness and figures' && (git push || true)"
    ;;
  part5)  # A10: missense + SpliceAI + fixes (~2-4 h: AlphaMissense download, phyloP API, SpliceAI on CPU)
    pip install -q matplotlib 2>/dev/null
    step a10_missense python a10_extend.py missense_fetch
    step spliceai_prep python a10_extend.py spliceai_prep
    if ! conda env list | grep -q "^spliceai "; then
      echo "== creating conda env spliceai (python 3.10, tensorflow-cpu, spliceai)"
      conda create -y -q -n spliceai -c conda-forge python=3.10 && conda run -n spliceai pip install -q "tensorflow-cpu==2.15.1" "setuptools<70" spliceai
    fi
    if [ ! -s data/hg38.fa ]; then
      echo "== downloading hg38.fa (~950 MB gz)"
      curl -sSL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip > data/hg38.fa
    fi
    step spliceai conda run -n spliceai spliceai -I a10/spliceai_in.vcf -O a10/spliceai_out.vcf -R data/hg38.fa -A grch38 -D 50
    step spliceai_parse python a10_extend.py spliceai_parse
    step a10_eval python a10_extend.py eval
    echo; echo "NEXT: git add a10/ results/ && git commit -m 'A10 evaluation' && (git push || true)"
    ;;
  part6)  # A11: NCBoost published benchmark re-scored with the rarity-only null (~1-2 h)
    pip install -q pyliftover 2>/dev/null
    step a11_fetch python a11_ncboost.py fetch
    step a11_eval  python a11_ncboost.py eval
    echo; echo "NEXT: git add a11/ results/ && git commit -m 'A11 evaluation' && (git push || true)"
    ;;
  *) echo "usage: bash run_next.sh part1|part2|part3|part4|part5|part6"; exit 1 ;;
esac
