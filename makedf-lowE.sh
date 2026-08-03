inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/lowE
outputdir=/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/lowE

. ~/get_token.sh
inputlist=${inputdir}/mc_MCP2025B_v10_06_00_09_prodgenie_corsika_proton_rockbox_lowenergydirt_sbnd_CV_caf_flat_caf_sbnd.list

python run_df_maker.py -c configs/numucc_1p0pi/sel_mup.py \
    -l "${inputlist}" \
    -o "${outputdir}/sel_mup"
