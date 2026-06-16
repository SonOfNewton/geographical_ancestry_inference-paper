#!/bin/bash

PARALLEL=$1
SOURCE_POP=$2
END_GEN=$3
SLIM=$4

if [ -z "$PARALLEL" ]; then
    PARALLEL=1
fi
echo "Using SLiM at: $SLIM"
echo "Using $PARALLEL cores"
export SLIM

# specify your number of cores for parallel computation
# seq 1 200 | xargs -I {} -P "$PARALLEL" bash -c '
#     echo "Running the $2 th replication..."
#     slim -d "WD='\''$1'\''" -d "REP=$2" code/generation/run_world1.slim
# ' _ "$CWD" {}
seq 1 200 | xargs -I {} -P "$PARALLEL" bash -c '
WD=$1
REP=$2
SOURCE_POP=$3
END_GEN=$4
SLIM=$5

echo "Running replication $REP ..."

"$SLIM" \
    -d "WD='\''$WD'\''" \
    -d "REP=$REP" \
    -d "SOURCE_POP=$SOURCE_POP" \
    -d "END_GEN=$END_GEN" \
    code/generation/run_world1.slim
' _ "$CWD" {} "$SOURCE_POP" "$END_GEN" "$SLIM"