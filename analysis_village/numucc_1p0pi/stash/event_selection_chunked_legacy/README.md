# Legacy per-file chunked event selection (stashed)

These files implemented a map/reduce workflow that processed **one `.df` file per job**
via `event_selection_chunk.py`. That path accumulated many edge-case bugs when merged
shards had sparse track/PID columns.

**Replaced by:** `event_selection_batched.py` and `scripts/run_event_selection_batched.sh`,
which group input files into **≤1 GB jobs**, load them with the same HDF helpers as the
notebook, run the notebook pipeline (`build_runner` / `ChunkRunner`), and aggregate with
`scripts/event_selection_aggregate.py`.

Stashed on request — kept for reference only.
