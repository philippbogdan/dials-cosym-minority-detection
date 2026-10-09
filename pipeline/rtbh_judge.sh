#!/bin/bash
# The third collection's judgement (amendments 19-21): today's DIALS on every trial (arrays saved),
# the frozen procedures on every trial, and today's DIALS on all 16 crystals at each wedge.
# Usage: rtbh_judge.sh P PROCEDURES (comma-separated cf.newproc configurations)
set -euo pipefail
P=${1:-8}; PROCS=${2:?procedures}
echo "=== trials $(date)"
seq 0 $((P - 1)) | xargs -P $P -I{} env CF_SHARD={} CF_NSHARDS=$P bash pipeline/run_py.sh cf.run_trials results/design/rtbh.jsonl today --save-arrays
echo "=== procedures $(date)"
seq 0 $((P * 2 - 1)) | xargs -P $P -I{} env CF_SHARD={} CF_NSHARDS=$((P * 2)) bash pipeline/run_py.sh cf.newproc rtbh "$PROCS"
echo "=== all $(date)"; bash pipeline/run_py.sh cf.run_trials results/design/rtbh_all.jsonl today
echo "=== done $(date)"
