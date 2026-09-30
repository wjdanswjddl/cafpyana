#!/usr/bin/env bash
# Stacked inclusive joint CC for PRL Product B constraint kinematics.
# Outputs under ``prl_syst_disk_root("B")/JointCC``.
#
# GENIE map: FSI_compare → GENIE_slim_v3 only. Do not sum VecFF/Ar23p/CCQE/MEC
# (those knobs are already inside slim_v3). align_joint_cc_genie.py keeps that
# slim_v3 joint as the combined cell (full off-diagonals).
#
# Usage:
#   bash analysis_village/numucc_1p0pi/scripts/run_cc_systs.sh
#   MAX_FILES=2 bash .../run_cc_systs.sh     # smoke

set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
cd "$THIS_DIR"

VENV="${REPO_ROOT}/envs/venv_py310_cafpyana/bin/activate"
if [[ -f "$VENV" ]]; then
    # shellcheck disable=SC1090
    source "$VENV"
fi
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
unset NUMUCC_SYST_DISK_CC_ROOT

FLUX_WORKERS="${FLUX_WORKERS:-8}"
G4_WORKERS="${G4_WORKERS:-8}"
GENIE_WORKERS="${GENIE_WORKERS:-12}"

bash run_syst_cc_joint_multisim_chunked.sh --syst-types Flux --workers "$FLUX_WORKERS"
bash run_syst_cc_joint_multisim_chunked.sh --syst-types G4 --workers "$G4_WORKERS"
# bash run_syst_cc_joint_multisim_chunked.sh --syst-types MCstat --workers 5

GENIE_RUN_GROUPS="${GENIE_RUN_GROUPS:-FSI_compare}" \
  bash run_syst_cc_joint_genie_chunked.sh --workers "$GENIE_WORKERS"
