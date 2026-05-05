df_dir=/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_01_102644__sel_all-mc-BNB_cosmics
save_name=$df_dir/ana

for file in $df_dir/*.df; do
    echo $file
    filename=$(basename $file)
    # echo "saving to $save_name/$filename.pkl"
    python event_selection_hist_chunk.py --df_dir $df_dir --filename $filename --file_type mc --save_content --save_name $save_name/$filename
done