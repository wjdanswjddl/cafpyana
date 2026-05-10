# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_mup_updatecalo_alpha_m,
       make_pandora_evtdf_mup_updatecalo_alpha_p,
       make_pandora_evtdf_mup_updatecalo_beta_m,
       make_pandora_evtdf_mup_updatecalo_beta_p,
       make_pandora_evtdf_mup_updatecalo_R_m,
       make_pandora_evtdf_mup_updatecalo_R_p,
       make_pandora_evtdf_mup_updatecalo_ccal_m,
       make_pandora_evtdf_mup_updatecalo_ccal_p,
       make_metadf]

ARGS = [{},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {}
        ]

NAMES = ["evt_alpha_m", 
         "evt_alpha_p", 
         "evt_beta_m", 
         "evt_beta_p", 
         "evt_R_m", 
         "evt_R_p", 
         "evt_ccal_m", 
         "evt_ccal_p", 
         "meta"]
