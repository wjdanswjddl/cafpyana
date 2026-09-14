# Detector-variation dataframes at sel_all (before any event selection).
# Writes:
#   hdr
#   evt_cv / trk_cv              — nominal calo (CV params, CAF E-field)
#   evt_<calo>_{p,m} / trk_…    — eight ±1σ calo universes
#   evt_efield / trk_efield     — CV calo params + double-anode E-field redo
#
# E-field tables use the distinct ``efield`` stem so they never overwrite
# nominal ``cv`` keys.
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_hdrdf]
ARGS = [{}]
NAMES = ["hdr"]

_CV = "ccal_cv-alpha_cv-beta_cv-R_cv"


def _add_univ(stem: str, updatecalo: str, *, updateefield: bool = False) -> None:
    kwargs = {"updatecalo": updatecalo}
    if updateefield:
        kwargs["updateefield"] = True
    DFS.append(make_pandora_evtdf_all_updatecalo)
    ARGS.append(dict(kwargs))
    NAMES.append(f"evt_{stem}")
    DFS.append(make_trkdf_updatecalo)
    ARGS.append(dict(kwargs))
    NAMES.append(f"trk_{stem}")


_add_univ("cv", _CV)

for ccal_var in ("p", "m"):
    _add_univ(f"ccal_{ccal_var}", f"ccal_{ccal_var}-alpha_cv-beta_cv-R_cv")

for alpha_var in ("p", "m"):
    _add_univ(f"alpha_{alpha_var}", f"ccal_cv-alpha_{alpha_var}-beta_cv-R_cv")

for beta_var in ("p", "m"):
    _add_univ(f"beta_{beta_var}", f"ccal_cv-alpha_cv-beta_{beta_var}-R_cv")

for R_var in ("p", "m"):
    _add_univ(f"R_{R_var}", f"ccal_cv-alpha_cv-beta_cv-R_{R_var}")

# Distinct stem — do not reuse ``cv`` / bare ``evt`` / ``trk``.
_add_univ("efield", _CV, updateefield=True)
