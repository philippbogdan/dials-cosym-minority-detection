#!/bin/bash
# Run xia2.multiplex (DIALS 3.30.1 / xia2 3.30, every default) on integrated data of the
# crystals listed in a file (one crystal name per line, in input order).
# Usage: multiplex.sh COLLECTION LISTFILE TAG NPROC [extra xia2.multiplex options...]
# With IMAGES=N in the environment, every crystal is first cut to its first N images
# (dials.slice_sequence), as the paper's reduced-rotation runs (section 3.3.3).
set -euo pipefail
ROOT=$(pwd)
export PATH=$ROOT/.venv-dials/env/bin:$PATH
COLL=$1; LIST=$2; TAG=$3; NPROC=$4; shift 4
W=$ROOT/.venv-data/multiplex/$TAG
OUT=$ROOT/results/fidelity/$TAG
rm -rf "$W"; mkdir -p "$W" "$OUT"
cd "$W"
ARGS=()
while read -r c; do
  [ -z "$c" ] && continue
  if [ -n "${IMAGES:-}" ]; then
    mkdir -p "$W/cut/$c"
    (cd "$W/cut/$c" && dials.slice_sequence "$ROOT/.venv-data/proc/$COLL/$c/integrated.expt" \
       "$ROOT/.venv-data/proc/$COLL/$c/integrated.refl" "image_range=1 $IMAGES" \
       output.experiments=cut.expt output.reflections=cut.refl > slice.log 2>&1)
    ARGS+=("$W/cut/$c/cut.expt" "$W/cut/$c/cut.refl")
  else
    ARGS+=("$ROOT/.venv-data/proc/$COLL/$c/integrated.expt" "$ROOT/.venv-data/proc/$COLL/$c/integrated.refl")
  fi
done < "$ROOT/$LIST"
START=$(date +%s)
xia2.multiplex "${ARGS[@]}" nproc=$NPROC "$@" > multiplex.stdout 2>&1 || echo "xia2.multiplex failed rc=$?"
END=$(date +%s)
echo "wall_s=$((END-START))" > "$OUT/timing.txt"
cp "$ROOT/$LIST" "$OUT/input_order.txt"
cp xia2.multiplex.log multiplex.stdout "$OUT/" 2>/dev/null || true
cp LogFiles/*dials.correlation_matrix.log "$OUT/" 2>/dev/null || true
cp LogFiles/*dials.cosym.log LogFiles/*dials.scale.log "$OUT/" 2>/dev/null || true
grep -A400 "OPTICS identified" xia2.multiplex.log | head -200 > "$OUT/optics_summary.txt" || true
