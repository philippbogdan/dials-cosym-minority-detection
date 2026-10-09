#!/bin/bash
set -euo pipefail
ROOT=$(pwd)
PY=${REPORT_PYTHON:-python3}
if [ -x "$ROOT/.venv-report/bin/python" ]; then PY=$ROOT/.venv-report/bin/python; fi
export PYTHONPATH=$ROOT${PYTHONPATH:+:$PYTHONPATH}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
exec "$PY" -m "$@"
