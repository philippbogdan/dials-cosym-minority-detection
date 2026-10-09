#!/bin/bash
# The third collection (PREREGISTRATION.md amendments 19-21), start to judgement: processing, joint
# cosym and the label verification (rtbh_prepare.sh), then today's DIALS, the frozen procedures and
# all 16 crystals at each wedge (rtbh_judge.sh). Usage: rtbh_all.sh P PROCEDURES
set -euo pipefail
bash pipeline/rtbh_prepare.sh "${1:-8}"
bash pipeline/rtbh_judge.sh "${1:-8}" "${2:-v6,v7,v8}"
