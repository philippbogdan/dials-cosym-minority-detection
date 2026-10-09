#!/bin/bash
set -euo pipefail
ROOT=$(pwd)
PY=$ROOT/.venv-dials/env/bin/python
[ -x "$PY" ] || { echo 'DIALS environment missing: run bash reproduce.sh env first' >&2; exit 2; }
export PATH=$ROOT/.venv-dials/env/bin:$PATH
export PYTHONPATH=$ROOT${PYTHONPATH:+:$PYTHONPATH}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
exec "$PY" -m "$@"
