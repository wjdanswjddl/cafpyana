# CV updatecalo + dE/dx resolution smear χ² on the production CV CAF list.
#
# Per plane, writes on ``trk_cv``:
#   chi2_{muon,proton}_new      — CV redo (no smear; same as updatecalo-cvonly)
#   chi2_{muon,proton}_smear13  — dE/dx × N(1, 0.13)  (GUMP nominal)
#   chi2_{muon,proton}_smear26  — dE/dx × N(1, 0.26)  (2× nominal)
#
# Selection / PID still uses score_tag="_new" (CV). Smear columns are stored for
# offline unisim histcounts. Writes: hdr, evt_cv, trk_cv.
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_hdrdf]
ARGS = [{}]
NAMES = ["hdr"]

_CV = "ccal_cv-alpha_cv-beta_cv-R_cv"
# Relative Gaussian widths passed to chi2pid.dedx(..., smear=).
_SMEAR = (0.13, 0.26)

_KW = {"updatecalo": _CV, "updatesmear": _SMEAR}

DFS += [make_pandora_evtdf_all_updatecalo, make_trkdf_updatecalo]
ARGS += [dict(_KW), dict(_KW)]
NAMES += ["evt_cv", "trk_cv"]
