"""Staged data access vocabulary for νμ CC 1p0π workflows.

Terminology (avoid overloading *chunk*)
---------------------------------------

**Map shard / CAF shard**
    One production ``.df`` file processed in isolation (e.g. ``syst_multisim_chunk.py``,
    ``event_selection_chunk.py``). Outputs are intermediate pickles; this is *not* a time
    slice of the detector exposure.

**HDF split**
    Internal key ``evt_<i>`` inside a single ``.df`` file — memory-bounded sequential reads.

**Exposure batch**
    A discrete, time-ordered slice of *real data* used for unblinding studies. Scripts such as
    ``selected_events.py`` split the sorted data dataframe into ``n_time_splits`` batches;
    ``selected_events_cumulative.py`` integrates batches ``0 .. K`` for each index ``K``.

**Data access stages** (analysis policy — configure inputs under ``dataset_locations``):

1. **Stage 1 — Fixed Dev sample:** use the development / calibration dataset only (e.g. data
   glob pointing at ``_Fixed_all.df`` or equivalent policy path).

2. **Stage 2 — Gen 1 independent batches:** Spring Gen 1 data processed one exposure batch at a
   time (no cumulative mixing across batches).

3. **Stage 3 — Gen 1 cumulative:** same dataset as Stage 2, but each run includes all exposure
   batches from the start of the run through batch ``K`` (monotonically increasing integrated
   POT).

4. **Stage 4 — Full Gen 1:** entire Gen 1 statistics in one pass (single batch spanning all
   data), subject to analysis blind/release policy.

CLI flags ``--exposure-batch-index`` / ``--exposure-batch-indices`` mirror the legacy
``--chunk_idx`` / ``--chunk_idxs`` names in the selection drivers.
"""

from __future__ import annotations

from enum import IntEnum


class DataAccessStage(IntEnum):
    """Sequential unblinding / staged-access policy (see module docstring)."""

    FIXED_DEV_ONLY = 1
    GEN1_INDEPENDENT_BATCHES = 2
    GEN1_CUMULATIVE_BATCHES = 3
    GEN1_FULL = 4
