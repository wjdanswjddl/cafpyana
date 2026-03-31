# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

# DFS = [make_pandora_evtdf_all_updatecalo,
#        make_metadf]

# ARGS = [{"updatecalo": "ccal_p-alpha_p-beta_p-R_p"}, 
# {}]

# NAMES = ["evt_ccal_p-alpha_p-beta_p-R_p", "meta"]


DFS = [make_metadf]

ARGS = [{}]

NAMES = ["meta"]


for ccal_var in ["p", "cv", "m"]:
    for alpha_var in ["p", "cv", "m"]:
        for beta_var in ["p", "cv", "m"]:
            for R_var in ["p", "cv", "m"]:
                DFS.append(make_pandora_evtdf_all_updatecalo)
                ARGS.append({"updatecalo": f"ccal_{ccal_var}-alpha_{alpha_var}-beta_{beta_var}-R_{R_var}"})
                NAMES.append(f"evt_ccal_{ccal_var}-alpha_{alpha_var}-beta_{beta_var}-R_{R_var}")
