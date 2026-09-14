#!/usr/bin/env bash
# Product B — cosmics measurement systematics (final-selected offbeam/intime).
#
# Forces INPUT_STAGE=final. To attach contamination-scaled ``SelectedRate`` in
# the aggregate NPZ (required by ``systematics-summary.ipynb``), set::
#
#   export COSMICS_SELECTED_MC_DF=/path/to/sel_mup.df   # HDF key ``evt``
#
# or pass ``--selected-mc-df`` through to aggregate via the chunked driver once
# that flag is wired; otherwise run ``syst_cosmics_aggregate.py`` with the flag.
#
# Usage: run_syst_measurement_cosmics.sh [args passed to run_syst_cosmics_chunked.sh]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export INPUT_STAGE=final
exec bash "$ROOT/scripts/run_syst_cosmics_chunked.sh" "$@"
