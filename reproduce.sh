#!/bin/bash
# Regenerate the experiment, or selected stages. Report needs only the committed records.
set -euo pipefail
cd "$(dirname "$0")"
export PARALLEL=${PARALLEL:-1}
case "$PARALLEL" in ''|*[!0-9]*|0) echo 'PARALLEL must be a positive integer' >&2; exit 2;; esac
REPORT_PY=${REPORT_PYTHON:-python3}
if ! "$REPORT_PY" -c 'import numpy, scipy, jinja2' >/dev/null 2>&1; then
  if [ ! -x .venv-report/bin/python ]; then python3 -m venv .venv-report; fi
  .venv-report/bin/python -m pip install -r requirements-report.txt
  REPORT_PY=.venv-report/bin/python
fi
PY=(bash pipeline/run_py.sh)
run_shards() { bash pipeline/shard.sh "$@"; }

stage_env() { bash infra/setup_env.sh; }
stage_process() {
  "$REPORT_PY" pipeline/manifest.py
  run_shards 130 bash pipeline/process_crystal.sh 0
}
stage_joint() { for c in cryo rt; do bash pipeline/joint.sh "$c" "$PARALLEL"; done; }
stage_fidelity() {
  bash pipeline/multiplex.sh cryo results/fidelity/cryo_all.txt cryo_all_default "$PARALLEL"
  IMAGES=10 bash pipeline/multiplex.sh cryo results/fidelity/cryo_all.txt cryo_all_10images "$PARALLEL"
  IMAGES=8 bash pipeline/multiplex.sh cryo results/fidelity/cryo_all.txt cryo_all_8images "$PARALLEL"
  "${PY[@]}" cf.fidelity_coords
}
stage_pilot() { for c in cryo rt; do "${PY[@]}" cf.pilot "$c"; done; }
stage_design() { for c in cryo rt rtbh; do "$REPORT_PY" -m cf.design "$c"; done; }
stage_trials() {
  run_shards 60 "${PY[@]}" cf.run_trials results/design/cryo.jsonl today --save-arrays
  run_shards 60 "${PY[@]}" cf.run_trials results/design/cryo.jsonl today --save-arrays --retry-errors
  run_shards 120 "${PY[@]}" cf.run_trials results/design/rt.jsonl today --save-arrays
  run_shards 6 "${PY[@]}" cf.run_trials results/design/cryo_all.jsonl today --save-arrays
  run_shards 12 "${PY[@]}" cf.run_trials results/design/rt_all.jsonl today --save-arrays
}
stage_analysis() {
  bash pipeline/analysis_stage.sh cryo "$PARALLEL"
  bash pipeline/analysis_stage.sh rt "$PARALLEL" cmat analytic new noisecal limits3 limits
  CF_CAUSE_VARIANTS=default,floor2 CF_CAUSE_WEDGES=60 CF_CAUSE_KMAX=6 bash pipeline/analysis_stage.sh rt "$PARALLEL" cause
  for c in cryo_amp rt_amp; do bash pipeline/analysis_stage.sh "$c" "$PARALLEL" cmat analytic noisecal limits3; done
  for c in cryo_xcal rt_xcal; do bash pipeline/analysis_stage.sh "$c" "$PARALLEL" limits3; done
  "${PY[@]}" cf.newproc cryo v2_f0.035,v2_f0.04,v2_f0.045,v2_f0.05
  run_shards 40 "${PY[@]}" cf.bootstrap cryo
}
stage_heldout() { bash pipeline/rtbh_all.sh "$PARALLEL" v6,v7,v8,v9,v10,v11; }
stage_limits() { for c in cryo rt; do "$REPORT_PY" -m cf.limits "$c"; done; }
stage_dev() {
  while read -r tag cfg; do run_shards 30 "${PY[@]}" cf.devrun cryo "$tag" "$cfg"; done < results/dev/configurations.txt
  bash pipeline/devalt_all.sh "$PARALLEL"
  "${PY[@]}" cf.devselect cryo
  "${PY[@]}" cf.devselect cryo_amp
  "${PY[@]}" cf.devjack
  for ((i=0; i<4; i++)); do
    "${PY[@]}" cf.devfast cryo "$i" 4
    "${PY[@]}" cf.devsparse "$i" 4 0.001
    "${PY[@]}" cf.devsparse_iter "$i" 4 0.001 10 3
  done
  for ((i=0; i<6; i++)); do "${PY[@]}" cf.devsparse_local "$i" 6 0.001 10 3; done
  for c in cryo rt; do "$REPORT_PY" -m cf.devfloor "$c" 3 4 5 6 8; done
  for ((i=0; i<5; i++)); do "${PY[@]}" cf.devsmall cryo "$i" 5; done
  for ((i=0; i<8; i++)); do "${PY[@]}" cf.devsmall rt "$i" 8; done
  "$REPORT_PY" -m cf.devrules
}
stage_e2e() {
  PATCH_MAX=15 FROZEN_PHIL="$PWD/dials-patch/frozen.phil" OVERLAY="$PWD/.venv-dials/overlay_v11" bash pipeline/par.sh "$PARALLEL" pipeline/run_py_patched.sh cf.run_trials results/design/e2e_sample.jsonl groups --tag groups_v11
  local max ver tag
  for spec in '14 10 groups_v10' '13 9 groups_v9' '11 8 groups_v8' '10 7 groups' '7 6 groups_v6_0007'; do
    read -r max ver tag <<< "$spec"
    PATCH_MAX=$max FROZEN_PHIL="$PWD/dials-patch/frozen-v$ver.phil" OVERLAY="$PWD/.venv-dials/overlay_v$ver" bash pipeline/par.sh "$PARALLEL" pipeline/run_py_patched.sh cf.run_trials results/design/e2e_sample.jsonl groups --tag "$tag"
  done
}
stage_tests() { bash pipeline/test_patch.sh; }
stage_report() { "$REPORT_PY" -m cf.report; "$REPORT_PY" -m cf.render; }
STAGES=(env process joint fidelity pilot design trials analysis heldout limits dev e2e tests report)
if [ $# -gt 0 ]; then STAGES=("$@"); fi
for s in "${STAGES[@]}"; do
  case "$s" in env|process|joint|fidelity|pilot|design|trials|analysis|heldout|limits|dev|e2e|tests|report) ;;
    *) echo "Unknown stage: $s" >&2; exit 2;;
  esac
  if [ "$s" != report ]; then git rev-parse HEAD > .cf_commit; fi
  echo "=== stage $s"
  "stage_$s"
done
