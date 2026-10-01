#!/bin/bash
# One GENIE dataframe submit. The nominal systematic is base × FSI v3
# (configs/numucc_1p0pi/sel_mup-genieslimwgts.py, GROUP=slim).
# Per-group breakdowns use sel_{mup,all}-geniewgts-knobgroups.py.
#
#   STAGE=sel_mup GROUP=slim   bash submit_mc_jobs_GENIE.sh
#   STAGE=sel_mup GROUP=Ar23p  bash submit_mc_jobs_GENIE.sh
#   STAGE=sel_all GROUP=slim   bash submit_mc_jobs_GENIE.sh
#   STAGE=sel_all GROUP=Ar23p  bash submit_mc_jobs_GENIE.sh
#   DRY_RUN=1 prints the run_df_maker line and does not submit.
#
# STAGE and GROUP are required so a bare invocation does not launch a grid job.
set -euo pipefail

STAGE="${STAGE:-}"
GROUP="${GROUP:-}"
DRY_RUN="${DRY_RUN:-0}"

if [[ -z "$STAGE" || -z "$GROUP" ]]; then
  echo "usage: STAGE=sel_mup|sel_all GROUP=slim|Ar23p|CCQE|MEC|RES|nonRES|DIS|Other bash $0" >&2
  exit 2
fi
if [[ "$STAGE" != "sel_mup" && "$STAGE" != "sel_all" ]]; then
  echo "STAGE must be sel_mup or sel_all (got $STAGE)" >&2
  exit 2
fi

# Slim sel_mup used the Ar23+ respin list. Knob-group and sel_all jobs use the
# 2026A Ar23+ xrootd list. Slim sel_all uses the Spring CV list.
RESPIN_LIST=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/Ar23+/ar23p_respin-xrootd.list
AR23_LIST=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_SBND2026A_AR23plus_knobs_BNBLight_CV_v1_00_01_flatcaf_sbnd_xrootd.list
SPRING_LIST=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list

if [[ "$GROUP" == "slim" && "$STAGE" == "sel_mup" ]]; then
  cfg=configs/numucc_1p0pi/sel_mup-genieslimwgts.py
  list="$RESPIN_LIST"
  out=sel_mup-wgts_genie_slim
  ngrid=4000
  extra_env=()
elif [[ "$GROUP" == "slim" && "$STAGE" == "sel_all" ]]; then
  cfg=configs/numucc_1p0pi/sel_all-geniewgts-knobgroups.py
  list="$SPRING_LIST"
  out=sel_all-wgts_genie_slim
  ngrid=3000
  extra_env=(GENIE_KNOB_GROUP=slim)
elif [[ "$STAGE" == "sel_mup" ]]; then
  cfg=configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py
  list="$AR23_LIST"
  out="sel_mup-wgts_genie_${GROUP}"
  ngrid=2000
  extra_env=(GENIE_KNOB_GROUP="$GROUP")
else
  cfg=configs/numucc_1p0pi/sel_all-geniewgts-knobgroups.py
  list="$AR23_LIST"
  out="sel_all-wgts_genie_${GROUP}"
  ngrid=2000
  extra_env=(GENIE_KNOB_GROUP="$GROUP")
fi

cmd=(python run_df_maker.py -c "$cfg" -l "$list" -o "$out" -ngrid "$ngrid")
echo "[submit-genie] STAGE=$STAGE GROUP=$GROUP"
echo "[submit-genie] ${extra_env[*]} ${cmd[*]}"
if [[ "$DRY_RUN" == "1" ]]; then
  exit 0
fi
env "${extra_env[@]}" "${cmd[@]}"
