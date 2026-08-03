# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_metadf, make_pandora_evtdf_mup_updatecalo]

ARGS = [{}, {"updatecalo": f"ccal_cv-alpha_cv-beta_cv-R_cv"}]

NAMES = ["meta", "evt_cv"]


for ccal_var in ["p", "m"]:
    DFS.append(make_pandora_evtdf_mup_updatecalo)
    ARGS.append({"updatecalo": f"ccal_{ccal_var}-alpha_cv-beta_cv-R_cv"})
    NAMES.append(f"evt_ccal_{ccal_var}")

for alpha_var in ["p", "m"]:
    DFS.append(make_pandora_evtdf_mup_updatecalo)
    ARGS.append({"updatecalo": f"ccal_cv-alpha_{alpha_var}-beta_cv-R_cv"})
    NAMES.append(f"evt_alpha_{alpha_var}")

for beta_var in ["p", "m"]:
    DFS.append(make_pandora_evtdf_mup_updatecalo)
    ARGS.append({"updatecalo": f"ccal_cv-alpha_cv-beta_{beta_var}-R_cv"})
    NAMES.append(f"evt_beta_{beta_var}")

for R_var in ["p", "m"]:
    DFS.append(make_pandora_evtdf_mup_updatecalo)
    ARGS.append({"updatecalo": f"ccal_cv-alpha_cv-beta_cv-R_{R_var}"})
    NAMES.append(f"evt_R_{R_var}")

