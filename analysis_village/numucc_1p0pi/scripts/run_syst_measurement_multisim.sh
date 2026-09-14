#!/usr/bin/env bash
# Product B — Flux / G4 / MCstat measurement systematics (final-selected dfs).
#
# Forces MC_DF_STAGE=final so chunks use ``get_univ_rates`` on sel_mup-style
# inputs (no selection re-walk). GENIE xsec is *not* produced here — use
# ``run_syst_measurement_genie.sh`` for that.
#
# Usage: run_syst_measurement_multisim.sh [args passed to run_syst_multisim_chunked.sh]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MC_DF_STAGE=final
exec bash "$ROOT/scripts/run_syst_multisim_chunked.sh" "$@"
