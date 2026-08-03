inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/Ar23+
list="${inputdir}/ar23p_respin-xrootd.list"

python run_df_maker.py \
 -c configs/numucc_1p0pi/sel_mup-genieslimwgts.py \
 -l "$list" \
 -o "sel_mup-wgts_genie_slim" \
 -ngrid 4000
