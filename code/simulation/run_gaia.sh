#!/bin/bash

# 检查参数
if [ "$#" -ne 3 ]; then
    echo "Usage: bash ooa_selected.sh <world> <map> <parallel>"
    echo "Example: bash ooa_selected.sh afro-eurasia friction 15"
    exit 1
fi

WORLD=$1
MAP=$2
PARALLEL=$3

if [ -z "$PARALLEL" ]; then
    PARALLEL=1
fi
echo "Using $PARALLEL cores"

SELECTION_FILE="output/tables/selected_reps_${WORLD}.csv"

if [ ! -f "$SELECTION_FILE" ]; then
    echo "Error: cannot find file $SELECTION_FILE "
    exit 1
fi

echo "Starting parallel replications for sampled replications in $SELECTION_FILE ..."

tail -n +2 "$SELECTION_FILE" | tr -d '"' | tr -d '\r' | xargs -I {} -P "$PARALLEL" bash -c '
    echo "Running selected replication {} with world=$1, map=$2 ..."
    Rscript --vanilla code/simulation/gaia.R "{}" "$1" "$2"
' _ "$WORLD" "$MAP"

# merge flux_strait files (necessary for parallel)
if [ "$WORLD" = "afro-eurasia" ]; then
    FINAL_CSV="data/flux/flux_strait_${WORLD}_${MAP}.csv"
    
    # fetch one file for header
    FIRST_REP=$(head -n 2 "$SELECTION_FILE" | tail -n 1 | tr -d '"\r')
    head -n 1 "data/flux/flux_strait_${WORLD}_${MAP}_${FIRST_REP}.csv" > "$FINAL_CSV"
    
    # merge
    for rep in $(tail -n +2 "$SELECTION_FILE" | tr -d '"\r'); do
        tail -n +2 "data/flux/flux_strait_${WORLD}_${MAP}_${rep}.csv" >> "$FINAL_CSV"
    done
    rm data/flux/flux_strait_${WORLD}_${MAP}_[0-9]*.csv
fi

echo "All selected replications finished!"