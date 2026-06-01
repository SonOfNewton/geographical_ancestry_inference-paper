#!/bin/bash

PARALLEL=$1
SOURCE_POP=$2
END_GEN=$3

if [ -z "$PARALLEL" ]; then
    PARALLEL=1
fi
echo "Using $PARALLEL cores"

# specify your number of cores for parallel computation
# seq 1 200 | xargs -I {} -P "$PARALLEL" bash -c '
#     echo "Running the $2 th replication..."
#     slim -d "WD='\''$1'\''" -d "REP=$2" code/generation/run_world1.slim
# ' _ "$CWD" {}
seq 1 200 | xargs -I {} -P "$PARALLEL" bash -c '
    echo "Running replication $2 ..."
    slim \
        -d "WD='\''$1'\''" \
        -d "REP=$2" \
        -d "SOURCE_POP=$3" \
        -d "END_GEN=$4" \
        code/generation/run_world1.slim
' _ "$CWD" {} "$SOURCE_POP" "$END_GEN"
echo "All replications finished!"