#!/bin/bash
# Finish today's DIALS on the rt trials still missing (in parallel), merge them into the rt
# results, then run the rt analysis stage. Usage: rt_finish_and_analyse.sh P IDS
set -uo pipefail
P=$1; IDS=$2
N=$(echo "$IDS" | tr ',' '\n' | wc -l)
seq 0 $((N - 1)) | xargs -P "$N" -I{} bash pipeline/run_py.sh cf.run_trials results/design/rt.jsonl today \
  --only "$IDS" --shard {} --nshards "$N" --tag today_extra --save-arrays
cat results/trials/rt/today_extra/shard_*.jsonl >> results/trials/rt/today/shard_114.jsonl
rm -rf results/trials/rt/today_extra
bash pipeline/analysis_stage.sh rt "$P" cmat analytic new noisecal limits3 cause
