seq 1 50 | xargs -I {} -P 14 bash -c '
    INDEX={}
    echo "Running Replication $INDEX..."
    slim -d "INDEX=$INDEX" -s $INDEX data/toy/toy.slim
'

