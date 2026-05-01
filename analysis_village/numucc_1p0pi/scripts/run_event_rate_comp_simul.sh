#!/bin/bash
n_time_splits=15

python selected_events.py --chunk_idx 1 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 2 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 4 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 6 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 7 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 9 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 11 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 12 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 14 --n_time_splits $n_time_splits & \
python selected_events.py --chunk_idx 15 --n_time_splits $n_time_splits & 
