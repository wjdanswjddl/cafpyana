# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_metadf]

ARGS = [{}]

NAMES = ["meta"]


for ccal_var in ["p", "cv", "m"]:
    for alpha_var in ["p", "cv", "m"]:
        for beta_var in ["p", "cv", "m"]:
            for R_var in ["p", "cv", "m"]:
                DFS.append(make_pandora_evtdf_2prong_vtxdist_updatecalo)
                ARGS.append({"updatecalo": f"ccal_{ccal_var}-alpha_{alpha_var}-beta_{beta_var}-R_{R_var}"})
                NAMES.append(f"evt_ccal_{ccal_var}-alpha_{alpha_var}-beta_{beta_var}-R_{R_var}")
