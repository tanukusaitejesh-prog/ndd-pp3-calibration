#!/bin/bash
# Usage: sync_commit.sh "message" file...   (copies files from the C: working folder into the WSL repo, commits, pushes)
set -e
SRC=/mnt/c/Users/saite/OneDrive/Desktop/Neuro/G-Test-LOCKIN/ndd_pp3_gate/ndd_pp3_gate
cd /root/ndd_pp3_gate
msg="$1"; shift
for f in "$@"; do
  if [ -f "$SRC/$f" ]; then mkdir -p "$(dirname "$f")"; cp "$SRC/$f" "$f"; fi
  git add "$f"
done
git commit -q -m "$msg

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push -q 2>&1 | tail -1
git log --oneline -1
