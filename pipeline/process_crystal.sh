#!/bin/bash
# Process one crystal from data/manifest.csv with DIALS 3.30.1, following the DIALS
# multi-crystal tutorial (import, find_spots, index, refine, integrate; defaults throughout).
# Usage: process_crystal.sh OFFSET   (row = OFFSET + CF_SHARD, 0-based over data rows)
# Outputs: .venv-data/proc/<collection>/<crystal>/integrated.{expt,refl}
#          results/processing/<collection>/<crystal>.json (summary)
# Idempotent and safe to restart.
set -euo pipefail
ROOT=$(pwd)
export PATH=$ROOT/.venv-dials/env/bin:$PATH
IDX=$(( ${1:-0} + ${CF_SHARD:-0} ))
LINE=$(sed -n "$((IDX + 2))p" data/manifest.csv)
IFS=, read -r COLL CRYSTAL SPECIES RECORD FILES MD5S BYTES URLS <<< "$LINE"
OUT=$ROOT/.venv-data/proc/$COLL/$CRYSTAL
SUM=$ROOT/results/processing/$COLL/$CRYSTAL.json
if [ -f "$OUT/integrated.refl" ] && [ -f "$SUM" ]; then echo "already done: $CRYSTAL"; exit 0; fi
NPROC=${NPROC:-1}
mkdir -p "$ROOT/.venv-data/tmp"
SCR=$(mktemp -d "$ROOT/.venv-data/tmp/process_${CRYSTAL}_XXXX")
trap 'rm -rf "$SCR"' EXIT
mkdir -p "$SCR/raw" "$SCR/work"
IFS=';' read -r -a FA <<< "$FILES"; IFS=';' read -r -a MA <<< "$MD5S"; IFS=';' read -r -a UA <<< "$URLS"
for i in "${!FA[@]}"; do
  f="$SCR/raw/${FA[$i]}"
  for attempt in 1 2 3 4 5; do
    curl -sSL --retry 5 -o "$f" "${UA[$i]}" && \
      [ "md5:$(python -c 'import hashlib,sys; print(hashlib.file_digest(open(sys.argv[1], "rb"), "md5").hexdigest())' "$f")" = "${MA[$i]}" ] && break
    echo "download retry $attempt for ${FA[$i]}"; sleep $((attempt * 20))
  done
  [ "md5:$(python -c 'import hashlib,sys; print(hashlib.file_digest(open(sys.argv[1], "rb"), "md5").hexdigest())' "$f")" = "${MA[$i]}" ] || { echo "md5 mismatch ${FA[$i]}"; exit 2; }
done
cd "$SCR/work"
if [ "$COLL" = cryo ]; then
  (cd "$SCR/raw" && tar xf *.tar && rm -f *.tar)
  dials.import $(find "$SCR/raw" -name "*.cbf*" | sort) > /dev/null
else
  dials.import "$SCR/raw/$CRYSTAL.nxs" > /dev/null
fi
dials.find_spots imported.expt nproc=$NPROC > /dev/null
if ! dials.index imported.expt strong.refl > /dev/null 2>&1; then
  echo "default indexing failed for $CRYSTAL; retrying with the insulin primitive cell as prior"
  dials.index imported.expt strong.refl unit_cell=67.5,67.5,67.5,109.47,109.47,109.47 > /dev/null
  echo "index_retry_with_cell" > index_note.txt
fi
dials.refine indexed.expt indexed.refl > /dev/null
dials.integrate refined.expt refined.refl nproc=$NPROC > /dev/null
python "$ROOT/pipeline/summarize_crystal.py" "$COLL" "$CRYSTAL" "$SPECIES" > summary.json
mkdir -p "$OUT.tmp" "$(dirname "$SUM")"
cp integrated.expt integrated.refl *.log summary.json "$OUT.tmp"/
[ -f index_note.txt ] && cp index_note.txt "$OUT.tmp"/
rm -rf "$OUT"; mv "$OUT.tmp" "$OUT"
cp summary.json "$SUM"
echo "done: $CRYSTAL"
