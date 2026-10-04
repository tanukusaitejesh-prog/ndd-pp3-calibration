#!/usr/bin/env bash
# Downloads every input except the Roulette mutation rates (about 0.8 GB in total).
# Safe to re-run: finished files are skipped and interrupted downloads resume.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${DATA:-$HERE/data}"
mkdir -p "$DATA"

# DeNovoWEST commit this gate was written and tested against (Kaplanis et al. 2020 de novo list)
DNW_COMMIT=49068922a840ff14351b0614f1e9e4be38b7df4f

fetch() {  # fetch URL DEST  - resumable; DEST.ok marks a finished download
  local url="$1" dest="$2"
  if [ -f "$dest.ok" ]; then echo "  already have $(basename "$dest")"; return 0; fi
  echo "  downloading $(basename "$dest")"
  if ! curl -fL --retry 20 --retry-delay 10 --retry-all-errors -C - -o "$dest" "$url"; then
    echo "  resume failed, restarting $(basename "$dest") from scratch"
    rm -f "$dest"
    curl -fL --retry 20 --retry-delay 10 --retry-all-errors -o "$dest" "$url"
  fi
  touch "$dest.ok"
}

echo "1/4 DeNovoWEST repository at the pinned commit (published de novo list, gene table, DDG2P; ~55 MB)"
if [ "$(git -C "$DATA/DeNovoWEST" rev-parse HEAD 2>/dev/null || true)" != "$DNW_COMMIT" ]; then
  rm -rf "$DATA/DeNovoWEST"
  git init -q "$DATA/DeNovoWEST"
  git -C "$DATA/DeNovoWEST" remote add origin https://github.com/HurlesGroupSanger/DeNovoWEST   # formerly queenjobo/DeNovoWEST
  if ! git -C "$DATA/DeNovoWEST" fetch -q --depth 1 origin "$DNW_COMMIT"; then
    git -C "$DATA/DeNovoWEST" fetch -q origin            # full history if the shallow fetch is refused
  fi
  git -C "$DATA/DeNovoWEST" -c advice.detachedHead=false checkout -q "$DNW_COMMIT"
fi
echo "  DeNovoWEST at $(git -C "$DATA/DeNovoWEST" rev-parse --short HEAD)"

echo "2/4 AlphaMissense hg38 (643 MB, CC BY 4.0, Zenodo 10.5281/zenodo.10813168)"
fetch "https://zenodo.org/records/10813168/files/AlphaMissense_hg38.tsv.gz?download=1" "$DATA/AlphaMissense_hg38.tsv.gz"
AM_MD5=9fd167735f16a1b87da6eb3e4c25fcb5          # published by Zenodo for this file
got=$(python -c "import hashlib,sys; h=hashlib.md5(); f=open(sys.argv[1],'rb'); [h.update(b) for b in iter(lambda: f.read(1<<20), b'')]; print(h.hexdigest())" "$DATA/AlphaMissense_hg38.tsv.gz")
if [ "$got" != "$AM_MD5" ]; then
  rm -f "$DATA/AlphaMissense_hg38.tsv.gz" "$DATA/AlphaMissense_hg38.tsv.gz.ok"
  echo "  AlphaMissense md5 mismatch ($got) - removed; re-run get_data.sh"; exit 1
fi
echo "  AlphaMissense md5 OK"

echo "3/4 GENCODE v32 annotation (the release AlphaMissense hg38 was mapped with)"
fetch "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz" \
      "$DATA/gencode.v32.annotation.gtf.gz"

echo "4/4 UCSC hg19 -> hg38 liftOver chain"
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" "$DATA/hg19ToHg38.over.chain.gz"

echo "Done. Inputs are in $DATA"
