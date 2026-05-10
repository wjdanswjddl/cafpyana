inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/lowE
outputdir=/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/lowE

chunktags=("aa" "ab" "ac" "ad" "ae" "af" "ag" "aj" "ai" "ak" "al" "am")
for tag in "${chunktags[@]}"; do
    . ~/get_token.sh
    chunkinput=${inputdir}/mc_MCP2025B_v10_06_00_09_prodgenie_corsika_proton_rockbox_lowenergydirt_sbnd_CV_caf_flat_caf_sbnd_${tag}

    echo "Processing chunk: ${tag}"

    python run_df_maker.py -c configs/numucc_1p0pi/sel_all_comparison.py \
        -l "${chunkinput}" \
        -o "${outputdir}/all${tag}"
done
