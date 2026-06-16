#!/bin/bash

export PATH=/usr/local/bin:$PATH
TOPOLOGY=$1
MODEL=$2
PARALLEL=$3

if [ -z "$TOPOLOGY" ] || [ -z "$MODEL" ]; then
    echo "Usage: bash slim_math.sh [line|square|cube|annulus|annulus2|annulus3] [MODEL]"
    exit 1
fi
if [ -z "$PARALLEL" ]; then
    PARALLEL=1
fi
echo "Using $PARALLEL cores"

seq 1 50 | xargs -I {} -P "$PARALLEL" bash -c '
INDEX=$1
TOPOLOGY=$2
MODEL=$3

echo "Running replication $INDEX with topology=$TOPOLOGY model=$MODEL"

slim \
    -d "INDEX=$INDEX" \
    -d "TOPOLOGY=\"$TOPOLOGY\"" \
    -d "MODEL=\"$MODEL\"" \
    -s $INDEX \
    code/generation/math.slim
' _ {} "$TOPOLOGY" "$MODEL"