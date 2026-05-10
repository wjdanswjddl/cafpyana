# selected event rates for MC files
# Full calorimetry uncertainty grid (3^4 = 81 evt dfs), keys aligned with makedf/chi2pid.CALO_VARIATIONS.
#
# updatecalo is only applied inside make_trkdf → sel_level must run past the "all" early return in
# make_pandora_evtdf. We use 2prong (same stage as sel_2prong_updatecalo.py). For variation at
# final mup selection, switch to make_pandora_evtdf_mup_updatecalo and adjust NAMES prefix if needed.
from itertools import product

from analysis_village.numucc_1p0pi.makedf.makedf import *

VARS = ("p", "cv", "m")

DFS = [make_metadf, make_hdrdf]
ARGS = [{}, {}]
NAMES = ["meta", "hdr"]

for ccal_var, alpha_var, beta_var, R_var in product(VARS, repeat=4):
    DFS.append(make_pandora_evtdf_2prong_updatecalo)
    ARGS.append(
        {"updatecalo": f"ccal_{ccal_var}-alpha_{alpha_var}-beta_{beta_var}-R_{R_var}"}
    )
    NAMES.append(
        f"evt_ccal_{ccal_var}-alpha_{alpha_var}-beta_{beta_var}-R_{R_var}"
    )

assert len(DFS) == len(ARGS) == len(NAMES)
