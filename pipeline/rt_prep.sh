#!/bin/bash
# rt analysis, first part (after today's DIALS): compact matrices aligned with DIALS's rows,
# analytic noise replicates, noise calibration against repeat segments. Restart-safe.
set -uo pipefail
P=${1:-2}
rm -rf results/cmat/rt .venv-data/boot_analytic/rt results/noisecal/rt.json
seq 0 7 | xargs -P "$P" -I{} env CF_SHARD={} CF_NSHARDS=8 bash pipeline/run_py.sh cf.cmatrices rt
seq 0 7 | xargs -P "$P" -I{} env CF_SHARD={} CF_NSHARDS=8 bash pipeline/run_py.sh cf.bootstrap rt analytic
bash pipeline/run_py.sh cf.noisecal rt
echo prep-done > .venv-data/rt_prep.done
