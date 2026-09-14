#!/usr/bin/env bash
# Legacy Product A1: GENIE systematics from existing sel_all weight DFs.
# Cut-stage vars = rate only; final measurement vars still get GENIE xsec.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
export MC_DF_STAGE=sel_all
exec bash "$ROOT/analysis_village/numucc_1p0pi/scripts/run_syst_genie_chunked.sh" "$@"
