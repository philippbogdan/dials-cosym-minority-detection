#!/bin/bash
# Split the current shard into P sub-shards.
set -euo pipefail
P=$1; RUNNER=$2; shift 2
T=${CF_SHARD:-0}; N=${CF_NSHARDS:-1}
PIDS=(); FAIL=0
for ((j=0; j<P; j++)); do
  CF_SHARD=$((T * P + j)) CF_NSHARDS=$((N * P)) bash "$RUNNER" "$@" &
  PIDS+=("$!")
done
for pid in "${PIDS[@]}"; do wait "$pid" || FAIL=1; done
exit "$FAIL"
