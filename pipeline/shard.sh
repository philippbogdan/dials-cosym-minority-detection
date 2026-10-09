#!/bin/bash
# Run N shards, with at most PARALLEL concurrent processes.
set -euo pipefail
N=$1; shift
P=${PARALLEL:-1}
FAIL=0
PIDS=()
for ((i=0; i<N; i++)); do
  CF_SHARD=$i CF_NSHARDS=$N "$@" &
  PIDS+=("$!")
  if [ ${#PIDS[@]} -ge "$P" ] || [ "$i" -eq "$((N-1))" ]; then
    for pid in "${PIDS[@]}"; do wait "$pid" || FAIL=1; done
    PIDS=()
  fi
done
exit "$FAIL"
