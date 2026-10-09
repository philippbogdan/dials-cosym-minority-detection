#!/bin/bash
# Fresh verification output is separate from the original run's test records.
set -euo pipefail
ROOT=$(pwd)
OV=${OVERLAY:-$ROOT/.venv-dials/overlay}
export OVERLAY=$OV
[ -d external/dials ] || { echo 'DIALS source missing: run bash reproduce.sh env first' >&2; exit 2; }
bash pipeline/make_overlay.sh
mkdir -p .venv-dials/test-source .venv-data/dials_data .venv-data/verification
SRC=$ROOT/.venv-dials/test-source/dials
rm -rf "$SRC"
cp -R external/dials "$SRC"
for p in "$ROOT"/dials-patch/${PATCHES:-*.patch}; do
  n=$(basename "$p" | cut -c1-4)
  if [ -n "${PATCH_MAX:-}" ] && [ "$((10#$n))" -gt "$PATCH_MAX" ]; then continue; fi
  (cd "$SRC" && patch --batch -p1 --forward < "$p" > /dev/null)
done
(
  export PATH=$ROOT/.venv-dials/env/bin:$PATH
  export PYTHONPATH=$OV
  export DIALS_DATA=$ROOT/.venv-data/dials_data
  cd "$SRC"
  "$ROOT/.venv-dials/env/bin/python" -m pytest -q -rs -p no:cacheprovider --regression tests/algorithms/correlation tests/command_line/test_correlation_matrix.py
) 2>&1 | tee .venv-data/verification/dials_patch.txt
bash pipeline/run_py.sh pytest -q -rs -p no:cacheprovider tests 2>&1 | tee .venv-data/verification/cf.txt
