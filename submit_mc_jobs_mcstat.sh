# MC statistical uncertainty weights (Poisson universes; sel_mup-mcstatwgts.py).
# Tune multisim_nuniv in configs/numucc_1p0pi/sel_mup-mcstatwgts.py if needed.

#inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
#list="${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list"
#
#python run_df_maker.py \
#  -c configs/numucc_1p0pi/sel_mup-mcstatwgts.py \
#  -l "$list" \
#  -o sel_mup-wgts_mcstat \
#  -ngrid 1000


inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/DENT
list="${inputdir}/EField_R00_xrootd.txt"

python run_df_maker.py \
  -c configs/numucc_1p0pi/sel_mup-mcstatwgts.py \
  -l "$list" \
  -o sel_mup-wgts_mcstat-EField_R00 \
  -ngrid 100

list="${inputdir}/EField_R30_Short_xrootd.txt"

python run_df_maker.py \
  -c configs/numucc_1p0pi/sel_mup-mcstatwgts.py \
  -l "$list" \
  -o sel_mup-wgts_mcstat-EField_R30_Short \
  -ngrid 100
