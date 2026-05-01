#!/bin/bash
n_time_splits=15

python selected_events.py --chunk_idx 3 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 5 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 8 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 10 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 13 --n_time_splits $n_time_splits & \
