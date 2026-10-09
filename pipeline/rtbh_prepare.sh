#!/bin/bash
# The third collection (PREREGISTRATION.md amendment 19) up to its trials: processing of its 16
# crystals, the indexing ambiguity resolved jointly, the label verification (stops here if a label
# disagrees). No trial: pipeline/rtbh_judge.sh runs them once every judged procedure is frozen.
# Usage: rtbh_prepare.sh P
set -euo pipefail
P=${1:-8}
FIRST=$(awk -F, 'NR>1 && $1=="rtbh" {print NR-2; exit}' data/manifest.csv)
echo "=== process $(date)"; bash pipeline/process_range.sh $FIRST $((FIRST + 15)) $((P > 1 ? P / 2 : 1)) 2
echo "=== joint $(date)"; bash pipeline/joint.sh rtbh $P
echo "=== verify $(date)"; bash pipeline/verify_rtbh.sh $P
echo "=== prepared $(date)"
