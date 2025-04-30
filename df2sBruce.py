#!/usr/bin/env python3 

import os
import sys
import pandas as pd
import uproot
import pyanalib.pandas_helpers as ph
import awkward as ak
import numpy as np

def convert(dfname):
    cohpidf = pd.read_hdf(dfname, key='cohpi')
    hdrdf = pd.read_hdf(dfname, key='hdr')
    mcnuwgtdf = pd.read_hdf(dfname, key='mcnuwgt')

    ## Collect POT and scale factor to the target POT
    this_pot = sum(hdrdf.pot)
    target_POT = 3.0e18
    POT_scale = target_POT / this_pot

    ## Make a df that matches mcnudf to slc
    matchdf = ph.multicol_merge(cohpidf.reset_index(), mcnuwgtdf.reset_index(),
                            left_on=[("entry", "",""), ("rec", "slc","tmatch", "idx")],
                            right_on=[("entry", "",""), ("rec.mc.nu..index", "","")], 
                            how="left") ## -- save all sllices
    matchdf.columns = ['.'.join(str(part) for part in col if part) for col in matchdf.columns]

    for col in matchdf.columns:
        if col.startswith("GENIEReWeight"):
            series = matchdf[col]
            first_valid = series.dropna().iloc[0]
            fill_value = 1.0
            matchdf[col] = matchdf[col].apply(lambda x: x if isinstance(x, list) else fill_value)

    # Convert the mcnuwgtdf format to an awkward array
    # This will make uproot save the output in the necessary format for sBruce trees
    mcnuwgtdf_out = pd.DataFrame({}, index=mcnuwgtdf.index)
    wgt_columns = [c for c in list(set(mcnuwgtdf.columns.get_level_values(0))) if c.startswith("GENIEReWeight")]
    for col in wgt_columns:
        mcnuwgtdf_out[col] = np.array([mcnuwgtdf[col][u].values for u in mcnuwgtdf[col].columns]).T.tolist()
        # mcnuwgtdf_out[col] = ak.to_dataframe(ak.Array(mcnuwgtdf_out[col]))
    return matchdf, mcnuwgtdf_out.reset_index()


def save(outf, matchdf, mcnuwgtdf):
    with uproot.recreate(outf) as f:
        f["mcnu"] = mcnuwgtdf
        f["matched"] = matchdf

def main(output, *files):
    # load
    matchdfs, mcnuwgtdfs = zip(*list(map(convert, files)))
    matchdf = pd.concat(matchdfs)
    mcnuwgtdf = pd.concat(mcnuwgtdfs)

    save(output, matchdf, mcnuwgtdf)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python df2sBruce.py <output.root> <input.df,>") 
        sys.exit(1)
    main(sys.argv[1], *sys.argv[2:])
