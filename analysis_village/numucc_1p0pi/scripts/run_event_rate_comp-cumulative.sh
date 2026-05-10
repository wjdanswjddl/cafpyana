#!/bin/bash
# Cumulative event-rate sweep.
#
# This calls the Python script ONCE; the script internally loops over chunks
# 0..(n_time_splits-1) and reuses the already-loaded dataframes, so the heavy
# I/O (get_ana_dfs and the systematic file loads) is only paid once.

n_time_splits=15
echo "Processing $n_time_splits cumulative chunks in a single Python invocation"
python selected_events_cumulative.py --n_time_splits $n_time_splits "$@"
