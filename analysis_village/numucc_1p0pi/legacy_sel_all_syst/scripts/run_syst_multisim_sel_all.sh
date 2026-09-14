#!/usr/bin/env bash
# Legacy Product A1: Flux/G4(/MCstat) systematics from existing sel_all weight DFs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
export MC_DF_STAGE=sel_all
export VAR_SET="${VAR_SET:-sel_all}"
exec bash "$ROOT/analysis_village/numucc_1p0pi/scripts/run_syst_multisim_chunked.sh" "$@"
