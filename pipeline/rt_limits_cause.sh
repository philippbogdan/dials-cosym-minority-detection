#!/bin/bash
# rt: noise calibration, refined limits and the cause replays, alongside the frozen procedures
# (rt_analysis_once.sh). Shard layouts are those of analysis_stage.sh with P=8 (limits3: 32,
# cause: 48), so work done by either is skipped by the other.
set -uo pipefail
P=${1:-4}
shards() { local n=$1; shift; seq 0 $((n - 1)) | xargs -P "$P" -I{} env CF_SHARD={} CF_NSHARDS="$n" bash pipeline/run_py.sh "$@"; }
echo "=== noisecal $(date)"; bash pipeline/run_py.sh cf.noisecal rt
echo "=== limits3 $(date)"; shards 32 cf.limits3 rt
echo "=== cause $(date)"; shards 48 cf.cause rt
echo "=== done $(date)"
