#!/bin/bash
n_time_splits=15
echo "Processing $n_time_splits time splits"
for chunk_idx in $(seq 5 10); do
    echo "Processing chunk_idx: $chunk_idx"
    python selected_events.py --chunk_idx $chunk_idx --n_time_splits $n_time_splits 
done
