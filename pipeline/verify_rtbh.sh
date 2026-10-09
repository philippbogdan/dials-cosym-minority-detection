#!/bin/bash
# Verify the rtbh labels (PREREGISTRATION.md amendment 19) before any trial: dials.cosym (defaults) on
# the 16 rtbh crystals with the labelled room-temperature bovine and human crystals, then
# cf.verify_rtbh (scaling, DIALS's correlation matrix, mean correlation with each labelled species).
# Exits non-zero if any crystal's label disagrees. Usage: verify_rtbh.sh NPROC
set -euo pipefail
ROOT=$(pwd)
export PATH=$ROOT/.venv-dials/env/bin:$PATH
NPROC=${1:-8}
J=$ROOT/.venv-data/joint/rtbh_verify
mkdir -p "$J"
cd "$J"
ARGS=()
while IFS=, read -r COLL CRYSTAL SPECIES REST; do
  if { [ "$COLL" = rt ] && { [ "$SPECIES" = bovine ] || [ "$SPECIES" = human ]; }; } || [ "$COLL" = rtbh ]; then
    d=$ROOT/.venv-data/proc/$COLL/$CRYSTAL
    [ -f "$d/integrated.refl" ] && ARGS+=("$d/integrated.expt" "$d/integrated.refl")
  fi
done < <(tail -n +2 "$ROOT/data/manifest.csv")
echo "$(( ${#ARGS[@]} / 2 )) crystals"
[ -f symmetrized.refl ] || dials.cosym "${ARGS[@]}" nproc=$NPROC > /dev/null
cd "$ROOT"
PYTHONPATH=$ROOT python -m cf.verify_rtbh
