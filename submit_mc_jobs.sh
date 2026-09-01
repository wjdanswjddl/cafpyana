#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-mc.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list -o sel_all-mc-BNB_cosmics -ngrid 1000
#
#inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/DENT
#list="${inputdir}/EField_R00_xrootd.txt"
#python run_df_maker.py \
#    -c configs/numucc_1p0pi/sel_all-mc.py \
#    -l $list \
#    -o sel_all-mc-BNB_cosmics-EField_R00 -ngrid 100
#
#list="${inputdir}/EField_R30_Short_xrootd.txt"
#python run_df_maker.py \
#    -c configs/numucc_1p0pi/sel_all-mc.py \
#    -l $list \
#    -o sel_all-mc-BNB_cosmics-EField_R30_Short -ngrid 100

#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-wgts-mc.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/Ar23+/ar23p_respin-xrootd.list -o sel_2prong-mc-BNB_cosmics -ngrid 1000
# python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-wgts-mc.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list -o sel_2prong-mc-BNB_cosmics -ngrid 2000
#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-mc.py -l /exp/sbnd/app/users/nrowe/cafpyana/new_joseph.list -o sel_all-mc-BNB_cosmics-josephsim -ngrid 1000


# WireMod
#inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_10/morestats
#
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong_updatecalo.py \
#       -l ${inputdir}/wiremod_bnb_20260211_yz_xrootd.list \
#       -o sel_2prong_wiremod_bnb_20260211_yz_updatecalo -ngrid 1000
#
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong_updatecalo.py \
#       -l ${inputdir}/wiremod_bnb_20260211_xtxw_xrootd.list \
#       -o sel_2prong_wiremod_bnb_20260211_xtxw_updatecalo -ngrid 1000


## GENIE vars (knob lists: makedf/geniesyst.GENIE_KNOB_GROUPS; unified config: sel_mup-geniewgts-knobgroups.py)
## Single group: prefix with GENIE_KNOB_GROUP=CCQE (valid keys: Ar23p CCQE ZExp MEC RES nonRES DIS Other)
## With -ngrid, run_df_maker forwards GENIE_KNOB_GROUP into each worker (else all evt_<Group> keys).
#inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
#
## CCQE (same physics via unified config + env)
#GENIE_KNOB_GROUP=CCQE python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py \
## or the thin alias config:
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-geniewgts_CCQE.py \
#    -l "${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list" \
#    -o sel_mup-wgts_genie_CCQE -ngrid 200
#
#
## MEC
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-geniewgts_MEC.py \
#    -l "${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list" \
#    -o sel_mup-wgts_genie_MEC -ngrid 200
#
#
## RES
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-geniewgts_RES.py \
#    -l "${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list" \
#    -o sel_mup-wgts_genie_RES -ngrid 200
#
#
## nonRES
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-geniewgts_nonRES.py \
#    -l "${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list" \
#    -o sel_mup-wgts_genie_nonRES -ngrid 200
#
#
## DIS
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-geniewgts_DIS.py \
#    -l "${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list" \
#    -o sel_mup-wgts_genie_DIS -ngrid 200
#
#
## Other
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-geniewgts_Other.py \
#    -l "${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list" \
#    -o sel_mup-wgts_genie_Other -ngrid 200
#
#
## AR23
#inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/Ar23+
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-geniewgts_Ar23p.py \
#    -l "${inputdir}/ar23p_respin-xrootd.list" \
#    -o sel_mup-wgts_genie_AR23p -ngrid 200


## BNB flux (all knobs in one pass: sel_mup-fluxwgts-knobgroups.py → HDF keys evt, hdr)
#inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
#mc="${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list"
#
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-fluxwgts-knobgroups.py \
#    -l "$mc" -o sel_mup-wgts_flux_all -ngrid 200
#
## Subset of knobs only: use thin configs configs/numucc_1p0pi/sel_mup-wgts_flux_<cat>.py



#list="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/old/MCP2025_GiBUU_flatcaf_xrootd.list"
#python run_df_maker.py \
#    -c configs/numucc_1p0pi/sel_all-mc.py \
#    -l $list \
#    -o sel_all-mc-GiBUU -ngrid 100
#
#list="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/old/MCP2025_GiBUU_flatcaf_xrootd.list"
#python run_df_maker.py \
#    -c configs/numucc_1p0pi/sel_mup.py \
#    -l $list \
#    -o sel_mup-mc-GiBUU -ngrid 100

#list="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list"
#list="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/DENT/aurora_MCP2026A_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_DENT_caf_flat_caf_sbnd_xrootd.list"
list="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_CV_caf_flat_caf_sbnd_xrootd.list"
#python run_df_maker.py \
#    -c configs/numucc_1p0pi/sel_all-mc.py \
#    -l $list \
#    -o sel_all-mc-CV -ngrid 2000

python run_df_maker.py \
    -c configs/numucc_1p0pi/sel_mup.py \
    -l $list \
    -o sel_mup-mc-CV-fvfix-chi2fix -ngrid 1000
