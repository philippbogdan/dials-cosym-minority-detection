#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
exec python3 pipeline/build_overlay.py
