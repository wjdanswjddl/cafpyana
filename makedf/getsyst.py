import uproot
import numpy as np
import pandas as pd
import awkward as ak
from pyanalib.pandas_helpers import *

def getsyst(f, systematics, nuind):
    if "globalTree" not in f:
        return pd.DataFrame(index=nuind.index)

    nuidx = pd.MultiIndex.from_arrays([nuind.index.get_level_values(0), nuind])

    syst_branches = ["wgts.name", "wgts.type", "wgts.nuniv"]
    systdf = loadbranches(f["globalTree/global/wgts"], syst_branches)
    systdf = systdf.reset_index().drop(columns =['entry', 'subentry'])
    systdf = systdf.wgts

    wgtdf = loadbranches(f["recTree"], ["rec.mc.nu.wgt.univ"]).rec.mc.nu.wgt
    wgtdf = wgtdf.rename(columns={"univ": "wgt"})
    wgtdf = wgtdf.rename_axis(['entry', 'rec.mc.nu..index', 'isyst', 'iuniv'])
    wgtdf = wgtdf.reset_index().set_index(['entry', 'rec.mc.nu..index'])

    systs = []
    for s in systematics:
        try:
            isyst = systdf.index[systdf['name'] == s][0]
        except IndexError:
            continue
        this_systs = []

        # Get weight type
        # +/- 1,2,3 sigma
        if systdf.type[isyst] == 3 and systdf.nuniv[isyst] == 1: # morph unisim
            s_morph = wgtdf[wgtdf.isyst == isyst].wgt.groupby(level=[0,1]).first()
            s_morph.name = (s, "morph")

            this_systs.append(s_morph)
        elif systdf.type[isyst] == 3 and systdf.nuniv[isyst] > 1: # +/- sigma unisim
            nsigma = wgtdf[wgtdf.isyst == isyst].wgt.groupby(level=[0,1]).size().values[0] // 2
            for isigma in range(nsigma):
                s_ps = wgtdf[wgtdf.isyst == isyst].wgt.groupby(level=[0,1]).nth(2*isigma)
                s_ps.name = (s, "ps%i" % (isigma+1))
                s_ms = wgtdf[wgtdf.isyst == isyst].wgt.groupby(level=[0,1]).nth(2*isigma+1)
                s_ms.name = (s, "ms%i" % (isigma+1))
 
                this_systs.append(s_ps)
                this_systs.append(s_ms)

        elif systdf.type[isyst] == 0: # multisi
            this_wgts =  wgtdf[wgtdf.isyst == isyst].wgt.groupby(level=[0,1]).head(250) # limit to 250 universes
            this_wgts = this_wgts.reset_index()
            this_wgts = this_wgts.pivot_table(values="wgt", index=["entry", "rec.mc.nu..index"], columns="iuniv")
            this_wgts.columns = pd.MultiIndex.from_tuples([(s, "univ_%i"% i) for i in range(len(this_wgts.columns))])
            for c in this_wgts.columns:
                this_systs.append(this_wgts[c])

        else:
            raise Exception("Cannot decode systematic uncertainty: %s" % s)

        for syst in this_systs:
            systs.append(syst)

    systs = pd.DataFrame(systs).T

    s_idx = systs.index.get_indexer(nuidx)
    systs_match = systs.iloc[s_idx]
    systs_match.loc[s_idx < 0, :] = 1.
    systs_match.index = nuind.index

    return systs_match

def getsyst_gundam(f, systematics, nuind):
    if "globalTree" not in f:
        return pd.DataFrame(index=nuind.index)

    nuidx = pd.MultiIndex.from_arrays([nuind.index.get_level_values(0), nuind])

    syst_branches = ["wgts.name", "wgts.type", "wgts.nuniv"]
    systdf = loadbranches(f["globalTree/global/wgts"], syst_branches)
    systdf = systdf.reset_index().drop(columns =['entry', 'subentry'])
    systdf = systdf.wgts

    wgtdf = loadbranches(f["recTree"], ["rec.mc.nu.wgt.univ"]).rec.mc.nu.wgt
    wgtdf = wgtdf.rename(columns={"univ": "wgt"})
    wgtdf = wgtdf.rename_axis(['entry', 'rec.mc.nu..index', 'isyst', 'iuniv'])
    wgtdf = wgtdf.reset_index().set_index(['entry', 'rec.mc.nu..index'])

    systs = []
    for s in systematics:
        try:
            isyst = systdf.index[systdf['name'] == s][0]
        except IndexError:
            continue
        this_systs = []

        if systdf.type[isyst] == 3 and systdf.nuniv[isyst] == 1: # morph unisim
            s_morph = wgtdf[wgtdf.isyst == isyst].wgt.groupby(level=[0,1]).agg(lambda x: list(x))
            s_morph.name = (s)
            this_systs.append(s_morph)

        elif systdf.type[isyst] == 3 and systdf.nuniv[isyst] > 1: # +/- sigma unisim
            def insert_middle(lst, value):
                mid = len(lst) // 2
                return lst[:mid] + [value] + lst[mid:]
            s_nsigma = wgtdf[wgtdf.isyst == isyst].wgt.groupby(level=[0,1]).agg(
                lambda x: insert_middle(list(x), 1.)
            )
            s_nsigma.name = (s)
            this_systs.append(s_nsigma)

        elif systdf.type[isyst] == 0: # multisim
            s_multisim = wgtdf[wgtdf.isyst == isyst].wgt.groupby(level=[0,1]).agg(lambda x: list(x)).head(250) # limit to 250 universes
            s_multisim.name = (s)
            this_systs.append(s_multisim)

        else:
            raise Exception("Cannot decode systematic uncertainty: %s" % s)

        for syst in this_systs:
            systs.append(syst)

    systs = pd.DataFrame(systs).T
    systs.columns = pd.MultiIndex.from_tuples([(col, '', '') for col in systs.columns])
    return systs
