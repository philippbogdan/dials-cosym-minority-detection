#!/bin/bash
# Resolve the indexing ambiguity once per collection: dials.cosym (DIALS 3.30.1 defaults) on
# the full integrated data of every crystal of the collection, then split back per crystal.
# Usage: joint.sh COLLECTION NPROC
set -euo pipefail
ROOT=$(pwd)
export PATH=$ROOT/.venv-dials/env/bin:$PATH
COLL=$1; NPROC=${2:-8}
J=$ROOT/.venv-data/joint/$COLL
mkdir -p "$J" "$ROOT/results/joint/$COLL"
cd "$J"
ARGS=()
for c in $(awk -F, -v c="$COLL" 'NR>1 && $1==c {print $2}' "$ROOT/data/manifest.csv"); do
  d=$ROOT/.venv-data/proc/$COLL/$c
  if [ -f "$d/integrated.refl" ]; then ARGS+=("$d/integrated.expt" "$d/integrated.refl"); else echo "missing $c"; fi
done
echo "${#ARGS[@]} files"
dials.cosym "${ARGS[@]}" nproc=$NPROC > /dev/null
cp dials.cosym.log dials.cosym.html "$ROOT/results/joint/$COLL/" 2>/dev/null || true
PYTHONPATH=$ROOT python -m cf.split_joint "$COLL"
