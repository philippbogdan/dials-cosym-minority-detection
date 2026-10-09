#!/bin/bash
# Install pinned DIALS and xia2, and fetch the matching DIALS source for its tests.
set -euo pipefail
ROOT=$(pwd)
DIALS_ENV=$ROOT/.venv-dials
mkdir -p "$DIALS_ENV/bin" external
if [ ! -x "$DIALS_ENV/bin/micromamba" ]; then
  case "$(uname -s)-$(uname -m)" in
    Linux-x86_64) PLAT=linux-64;;
    Darwin-arm64) PLAT=osx-arm64;;
    Darwin-x86_64) PLAT=osx-64;;
    *) echo 'No configured binary platform for this host' >&2; exit 2;;
  esac
  curl -fLsS "https://micro.mamba.pm/api/micromamba/$PLAT/latest" | tar -xj -C "$DIALS_ENV" bin/micromamba
fi
export MAMBA_ROOT_PREFIX=$DIALS_ENV/root
if [ ! -x "$DIALS_ENV/env/bin/dials.version" ]; then
  "$DIALS_ENV/bin/micromamba" create -y -p "$DIALS_ENV/env" -c conda-forge python=3.12 dials=3.30.1 xia2=3.30.0 pytest dials-data jinja2
fi
"$DIALS_ENV/env/bin/dials.version"
if [ ! -d external/dials ]; then
  curl -fLsS https://github.com/dials/dials/archive/refs/tags/v3.30.1.tar.gz | tar -xz -C external
  mv external/dials-3.30.1 external/dials
fi
