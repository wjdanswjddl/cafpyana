import os
import sys
import pandas as pd
import uproot
import pyanalib.pandas_helpers as ph
import awkward as ak
import numpy as np

def get_keys(file):
    with pd.HDFStore(file, mode='r') as store:
        keys = store.keys()       # list of all keys in the file
        print("Keys:", keys)
    return keys

def get_n_split(file):
    this_split_df = pd.read_hdf(file, key="split")
    this_n_split = this_split_df.n_split.iloc[0]
    return this_n_split

def make_gump_ttree_mc(dfname):
    mc_split_df = pd.read_hdf(dfname, key="split")
    mc_n_split = get_n_split(dfname)
    print("mc_n_split: %d" %(mc_n_split))
    keys = get_keys(dfname)

    for k in keys:
        if('evt' in k):
            recodf = pd.read_hdf(dfname, key=k)
        elif('hdr' in k):
            hdrdf = pd.read_hdf(dfname, key=k)
        elif('mcnuwgtslim' in k):
            mcnuwgtdf = pd.read_hdf(dfname, key=k)

    ## Collect POT and scale factor to the target POT
    this_pot = sum(hdrdf.pot)
    target_POT = 4.6e18
    POT_scale = target_POT / this_pot

    wgt_columns = [c for c in list(set(mcnuwgtdf.columns.get_level_values(0))) if c.startswith("GENIEReWeight")]

    truedf_wgt_out = pd.DataFrame({}, index=mcnuwgtdf.index)
    for col in wgt_columns:
        truedf_wgt_out[col] = np.array([mcnuwgtdf[col][u].values for u in mcnuwgtdf[col].columns]).T.tolist()

    non_syst_columns = [col for col in mcnuwgtdf.columns if not col[0].startswith("GENIEReWeight")]
    truedf_out = mcnuwgtdf[non_syst_columns]
    truedf_out.columns = truedf_out.columns.get_level_values(0)
    truedf_out = pd.concat([truedf_out, truedf_wgt_out], axis = 1)
    
    return recodf, truedf_out

def make_gump_ttree_data(dfname):
    recodf = pd.read_hdf(dfname, key='evt')
    hdrdf = pd.read_hdf(dfname, key='hdr')

    ## Collect POT and scale factor to the target POT
    POT_scale = sum(hdrdf.pot)

    return recodf
