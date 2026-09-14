#!/usr/bin/env bash
# Product B — GENIE measurement systematics (sel_mup + mcnu → rate + xsec).
#
# Forces MC_DF_STAGE=final so chunk-map reads already-selected evt+mcnu and
# builds the full GENIE xsec path (R_u @ N_gen^CV + Δbg).
#
# Usage: run_syst_measurement_genie.sh [args passed to run_syst_genie_chunked.sh]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MC_DF_STAGE=final
exec bash "$ROOT/scripts/run_syst_genie_chunked.sh" "$@"
