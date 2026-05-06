"""
Helper functions for working with the HDF5 df files.  

Each dataset (identified by a *key*) is stored as a set of smaller dfs (one per split),
so that very large samples can be handled in chunks.  
"""

import pandas as pd
from os import path
from tqdm import tqdm
import numpy as np

# utils for the case where events are split into keys within a single file
def get_n_split(file):
    this_split_df = pd.read_hdf(file, key="split")
    this_n_split = this_split_df.n_split.iloc[0]
    return this_n_split

def print_keys(file):
    with pd.HDFStore(file, mode='r') as store:
        keys = store.keys()       # list of all keys in the file
        # print("Keys:", keys)

def load_dfs(file, keys2load, n_max_concat=100):
    out_df_dict = {}
    this_n_keys = get_n_split(file)
    n_concat = min(n_max_concat, this_n_keys)
    for key in keys2load:
        dfs = []  # collect all splits for this key
        for i in range(n_concat):
            this_df = pd.read_hdf(file, key=f"{key}_{i}")
            dfs.append(this_df)
        out_df_dict[key] = pd.concat(dfs, ignore_index=False)

    return out_df_dict


# utils for the case where events are split into different files (which also has split keys like above)
import glob

def dfs_from_dir(
    search_dir,
    filename_str=None,
    keys2load=['hdr', 'evt'],
    n_max_concat=3,
    ):
    """
    Loops over files in the directory that match each chunk_tag as a substring/glob ("*<chunk_tag>*").
    Loads the specified dfs from each file, concats them, and keeps the first level of multiindex value
    unique (__ntuple) by offsetting it by the previously summed lengths.
    """

    df_lists = {k: [] for k in keys2load}
    ntuple_offset = np.int64(0)

    # Compose the search directory
    matched_files = set()

    if filename_str is not None:
        pattern = path.join(search_dir, f"*{filename_str}*.df")
        matched_files = glob.glob(pattern)
    else:
        raise ValueError("No filename string provided, aborting...")

    # Sort files for reproducibility and then limit to n_max_concat files
    files_to_process = sorted(matched_files)[:n_max_concat]
    print(f"Found {len(files_to_process)} files to process")
    print(f"Files to process: {files_to_process}")

    for mc_file in tqdm(files_to_process):
        # print_keys(mc_file)
        mc_n_split = get_n_split(mc_file)
        # print(f"Reading file {mc_file}, mc_n_split: {mc_n_split}")
        try:
            mc_dfs = load_dfs(mc_file, keys2load, n_max_concat=n_max_concat)
        except Exception as e:
            print(f"Error loading file {mc_file}: {e}")
            continue

        # Build a dense remapping of this file's ntuple values to avoid overflow.
        # Using max()+1 as the bump causes the offset to grow with the magnitude of
        # the raw indices; remapping to a contiguous range keeps growth proportional
        # to the actual number of unique ntuples instead.
        ref_df = mc_dfs[keys2load[0]]
        if isinstance(ref_df.index, pd.MultiIndex):
            raw_ntuple_vals = ref_df.index.get_level_values(0)
        else:
            raw_ntuple_vals = ref_df.index
        unique_ntuples = np.array(sorted(raw_ntuple_vals.unique()))
        ntuple_remap = {old: np.int64(ntuple_offset + i) for i, old in enumerate(unique_ntuples)}
        n_unique = np.int64(len(unique_ntuples))

        # Make __ntuple unique by remapping first level index (if exists)
        for df_key in keys2load:
            df = mc_dfs[df_key]
            if isinstance(df.index, pd.MultiIndex):
                names = df.index.names
                # __ntuple should be at level 0
                if "__ntuple" in names:
                    idx_loc = names.index("__ntuple")
                else:
                    idx_loc = 0
                new_tuples = []
                for tup in df.index:
                    tup = list(tup)
                    tup[idx_loc] = ntuple_remap[tup[idx_loc]]
                    new_tuples.append(tuple(tup))
                df.index = pd.MultiIndex.from_tuples(new_tuples, names=names)
            else:
                if df.index.name == "__ntuple":
                    df.index = df.index.map(ntuple_remap)

            df_lists[df_key].append(df)

        # advance offset by the number of unique ntuples in this file
        ntuple_offset += n_unique

    concat_dfs = {k: pd.concat(df_lists[k], axis=0, sort=False) for k in df_lists.keys()}
    return concat_dfs
