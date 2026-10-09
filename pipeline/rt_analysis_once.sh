#!/bin/bash
# The rt analysis stage, run once even when submitted to several places: the first to start
# takes the lock (mkdir is atomic on the shared filesystem), the others exit.
# Usage: rt_analysis_once.sh P
P=$1
mkdir .venv-data/rt_analysis.lock 2>/dev/null || { echo "rt analysis already running elsewhere"; exit 0; }
rm -rf results/cmat/rt .venv-data/boot_analytic/rt results/trials/rt/v[0-9]* results/limits/rt_refined results/noisecal/rt.json
bash pipeline/analysis_stage.sh rt "$P" cmat analytic new noisecal limits3 cause
