#!/bin/bash
# Analysis of one collection after today's DIALS has run on all its trials (arrays saved):
# compact matrices, analytic noise replicates, the frozen procedures (cf.newproc), noise
# calibration against repeat segments, refined limits, pre-registered limits, cause ablations.
# Runs on one machine with P parallel workers (each stage is sharded and idempotent).
# Usage: analysis_stage.sh DESIGN P [STAGE ...]   (stages: cmat analytic new noisecal limits3 limits cause)
set -euo pipefail
DESIGN=$1; P=$2; shift 2
STAGES=("$@"); [ ${#STAGES[@]} -eq 0 ] && STAGES=(cmat analytic new noisecal limits3 limits cause)
shards() {  # NSHARDS MODULE ARGS... : NSHARDS shards, P at a time
  local n=$1; shift
  seq 0 $((n - 1)) | xargs -P "$P" -I{} env CF_SHARD={} CF_NSHARDS="$n" bash pipeline/run_py.sh "$@"
}
for s in "${STAGES[@]}"; do
  echo "=== $s $(date)"
  case $s in
    cmat)     shards 16 cf.cmatrices "$DESIGN" ;;
    analytic) shards 16 cf.bootstrap "$DESIGN" analytic ;;
    new)      shards $((P * 6)) cf.newproc "$DESIGN" ;;
    noisecal) bash pipeline/run_py.sh cf.noisecal "$DESIGN" ;;
    limits3)  shards $((P * 4)) cf.limits3 "$DESIGN" ;;
    limits)   bash pipeline/run_py.sh cf.limits "$DESIGN" ;;
    cause)    shards $((P * 6)) cf.cause "$DESIGN" ;;
  esac
done
echo "=== done $(date)"
