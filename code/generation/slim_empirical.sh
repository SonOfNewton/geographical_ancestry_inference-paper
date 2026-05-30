#!/bin/bash

PARALLEL=$1

if [ -z "$PARALLEL" ]; then
    PARALLEL=1
fi
echo "Using $PARALLEL cores"

# specify your number of cores for parallel computation
seq 1 200 | xargs -I {} -P "$PARALLEL" bash -c '
    echo "Running the $2 th replication..."
    slim -d "WD='\''$1'\''" -d "REP=$2" code/generation/run_world1.slim
' _ "$CWD" {}

echo "All replications finished!"