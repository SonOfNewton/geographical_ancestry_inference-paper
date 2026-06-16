#!/bin/bash

TOPOLOGY=$1
MODEL=$2
PARALLEL=$3
SLIM=$4

if [ -z "$TOPOLOGY" ] || [ -z "$MODEL" ]; then
    echo "Usage: bash slim_math.sh [TOPOLOGY] [MODEL] [PARALLEL] [SLIM_PATH]"
    exit 1
fi

if [ -z "$PARALLEL" ]; then
    PARALLEL=1
fi

if [ -z "$SLIM" ]; then
    echo "ERROR: No SLiM executable path supplied."
    exit 1
fi

if [ ! -x "$SLIM" ]; then
    echo "ERROR: SLiM executable not found or not executable:"
    echo "$SLIM"
    exit 1
fi

echo "Using SLiM at: $SLIM"
echo "Using $PARALLEL cores"

export SLIM

seq 1 50 | xargs -I {} -P "$PARALLEL" bash -c '
INDEX=$1
TOPOLOGY=$2
MODEL=$3
SLIM=$4

echo "Running replication $INDEX with topology=$TOPOLOGY model=$MODEL"

"$SLIM" \
    -d "INDEX=$INDEX" \
    -d "TOPOLOGY=\"$TOPOLOGY\"" \
    -d "MODEL=\"$MODEL\"" \
    -s "$INDEX" \
    code/generation/math.slim
' _ {} "$TOPOLOGY" "$MODEL" "$SLIM"