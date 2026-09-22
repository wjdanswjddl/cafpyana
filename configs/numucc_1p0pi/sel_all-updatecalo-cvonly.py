# CV-only updatecalo at sel_all: rewrite chi2_{muon,proton}_new (and ndof_*)
# with nominal calo params. No ±1σ calo universes, no E-field redo.
#
# Writes: hdr, evt_cv, trk_cv
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_hdrdf]
ARGS = [{}]
NAMES = ["hdr"]

_CV = "ccal_cv-alpha_cv-beta_cv-R_cv"

DFS += [make_pandora_evtdf_all_updatecalo, make_trkdf_updatecalo]
ARGS += [{"updatecalo": _CV}, {"updatecalo": _CV}]
NAMES += ["evt_cv", "trk_cv"]
