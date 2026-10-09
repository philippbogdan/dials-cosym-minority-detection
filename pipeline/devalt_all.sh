#!/bin/bash
# Development runs on amplitude correlations (cf/devalt.py, amendment 14), every cryo trial: each
# variant as one shard (results already present are skipped), P variants at a time.
P=${1:-4}
printf '%s\n' "amp" "ampw" "nscore" "amp ratio" "amp c0.1" "amp c0.08" |
  xargs -P "$P" -I{} bash -c 'bash pipeline/run_py.sh cf.devalt $(echo {} | cut -d" " -f1) 0 1 $(echo {} | cut -s -d" " -f2)'
