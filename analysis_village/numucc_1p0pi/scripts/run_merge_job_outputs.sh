basedir="/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"

python merge_grid_job_dfs.py --df-dir $basedir/2026_05_23_235202__sel_mup-wgts_genie_slim --filename sel_mup --fv perTPC &
python merge_grid_job_dfs.py --df-dir $basedir/2026_05_23_215228__sel_mup-wgts_genie_slim --filename sel_mup --fv perTPC

#
#python merge_grid_job_dfs_trk12.py --df-dir $basedir/2026_05_12_113644__sel_2prong-mc-BNB_cosmics-0xSCE --filename sel_2prong --fv perTPC & \
#
#python merge_grid_job_dfs_only.py --df-dir $basedir/2026_05_11_035756__sel_all-data-OffBeamLight --filename-str sel_all & \
