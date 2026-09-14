#!/usr/bin/env bash
# Legacy Product A1: cosmics systematics from existing sel_all DFs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
export MC_DF_STAGE=sel_all
exec bash "$ROOT/analysis_village/numucc_1p0pi/scripts/run_syst_cosmics_chunked.sh" "$@"
