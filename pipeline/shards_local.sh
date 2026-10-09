#!/bin/bash
# Run shards FIRST..LAST of a design with PAR parallel workers on one machine.
# Usage: shards_local.sh DESIGN PROCEDURE NSHARDS FIRST LAST PAR [extra run_trials args]
DESIGN=$1; PROC=$2; NSH=$3; FIRST=$4; LAST=$5; PAR=$6; shift 6
seq $FIRST $LAST | xargs -P $PAR -I{} bash pipeline/run_py.sh cf.run_trials $DESIGN $PROC --shard {} --nshards $NSH "$@"
