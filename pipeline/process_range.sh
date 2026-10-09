#!/bin/bash
# Process manifest rows FIRST..LAST (0-based, inclusive), PAR crystals at a time, NPROC cores each.
# Usage: process_range.sh FIRST LAST PAR NPROC
FIRST=$1; LAST=$2; PAR=$3; export NPROC=$4
seq $FIRST $LAST | xargs -P $PAR -I{} bash -c 'CF_SHARD=0 bash pipeline/process_crystal.sh {} > .venv-data/log_row_{}.txt 2>&1; echo "row {} rc=$?"; tail -1 .venv-data/log_row_{}.txt'
