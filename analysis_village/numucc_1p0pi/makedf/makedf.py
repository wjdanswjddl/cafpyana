from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *

from analysis_village.numucc_1p0pi.makedf.selections import *
from analysis_village.numucc_1p0pi.categories import DETECTOR, TRK_CALO_DET, in_fv
from analysis_village.numucc_1p0pi.event_selection_pipeline_def import (
    CAF_SEL_LEVEL_TO_STAGE,
    CAF_SEL_LEVELS_NEED_TRACKS,
    apply_selection_pipeline,
)
from makedf.geniesyst import *
from makedf.bnbsyst import *

def make_spine_evtdf(f):
    # load slices and particles
    partdf = make_epartdf(f)

    df = make_eslcdf(f)

    # load the proton and muon candidates
    primary = partdf.is_primary
    mudf = partdf[primary & (partdf.pid == 2)].sort_values(partdf.index.names[:2] + [("length", "", "")]).groupby(level=[0,1]).last()
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])

    pdf = partdf[primary & (partdf.pid == 4)].sort_values(partdf.index.names[:2] + [("length", "", "")]).groupby(level=[0,1]).last()
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])

    df = multicol_merge(df, mudf, left_index=True, right_index=True, how="left", validate="one_to_one")
    df = multicol_merge(df, pdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    # in case we want to cut out other objects -- save the highest energy of each other particle
    lead_gamma_energy = partdf.ke[primary & (partdf.pid == 0)].groupby(level=[0,1]).max().rename("lead_gamma_energy")
    df = multicol_add(df, lead_gamma_energy)

    lead_elec_energy = partdf.ke[primary & (partdf.pid == 1)].groupby(level=[0,1]).max().rename("lead_elec_energy")
    df = multicol_add(df, lead_elec_energy)

    lead_pion_length = partdf.length[primary & (partdf.pid == 3)].groupby(level=[0,1]).max().rename("lead_pion_length")
    df = multicol_add(df, lead_pion_length)

    subl_muon_length = partdf[primary & (partdf.pid == 2)].sort_values(partdf.index.names[:2] + [("length", "", "")]).length.groupby(level=[0,1]).nth(-2).rename("subl_muon_length")
    df = multicol_add(df, subl_muon_length)

    subl_proton_length = partdf[primary & (partdf.pid == 4)].sort_values(partdf.index.names[:2] + [("length", "", "")]).length.groupby(level=[0,1]).nth(-2).rename("subl_proton_length")
    df = multicol_add(df, subl_proton_length)

    # Apply pre-selection: Require fiducial vertex, at least one muon, at least one proton

    # require both muon and proton to be present
    df = df[~np.isnan(df.mu.pid) & ~np.isnan(df.p.pid)]

    # require fiducial vertex (Gen-1 FV when DETECTOR is perTPC)
    df = df[in_fv(df.vertex, detector=DETECTOR)]

    return df

# ===== selection stages =====
def make_pandora_evtdf_all(f, sel_level="all", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=False, 
                       trkScoreCut=False, trkDistCut=1000., cutClearCosmic=False, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    #print("LEN OF DF:", len(df))
    #print("CHECKPOINT LAST")
    return df

def make_pandora_evtdf_2prong(f, sel_level="2prong", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_pandora_evtdf_2prong_wcandidates(f, sel_level="2prong_wcandidates", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_pandora_evtdf_mup(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

# ===== syst weights =====

def make_pandora_evtdf_wgts(f, include_weights=True, multisim_nuniv=1000, wgt_types=["bnb","g4"], slim=True, 
                       trkScoreCut=False, trkDistCut=10., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_pandora_evtdf_2prong_wgts(f, sel_level="2prong", include_weights=True, multisim_nuniv=100, wgt_types=["bnb","genie","g4"], slim=True, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, genie_systematics=None,
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_pandora_evtdf_mup_wgts(f, sel_level="mup", include_weights=True, multisim_nuniv=200, wgt_types=["bnb","genie","g4"], slim=True, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_pandora_evtdf_mup_genieslimwgts(f, sel_level="mup", include_weights=True, genie_multisim_nuniv=100, wgt_types=["genie"], slim=True, 
                        genie_systematics=gen1_systematics,
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, genie_multisim_nuniv=genie_multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            genie_systematics=genie_systematics, trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_mcnudf_genieslimwgts(f, multisim_nuniv=100, genie_multisim_nuniv=100, genie_systematics=gen1_systematics, slim=True):
    return make_mcnudf(f, include_weights=True, multisim_nuniv=multisim_nuniv, genie_multisim_nuniv=genie_multisim_nuniv, wgt_types=["genie"], slim=slim, genie_systematics=genie_systematics)


def make_pandora_evtdf_mup_mc_multisim(
    f,
    sel_level="mup",
    include_weights=True,
    multisim_nuniv=1000,
    wgt_types=None,
    slim=False,
    genie_systematics=None,
    flux_systematics=None,
    trkScoreCut=False,
    trkDistCut=100.0,
    cutClearCosmic=True,
    **trkArgs,
):
    """mup selection + configurable multisim weights (wgt_types / slim / systematics lists)."""
    if wgt_types is None:
        raise ValueError("make_pandora_evtdf_mup_mc_multisim requires wgt_types, e.g. ['genie'], ['bnb'], ['g4']")
    return make_pandora_evtdf(
        f,
        sel_level=sel_level,
        include_weights=include_weights,
        multisim_nuniv=multisim_nuniv,
        wgt_types=wgt_types,
        slim=slim,
        genie_systematics=genie_systematics,
        flux_systematics=flux_systematics,
        trkScoreCut=trkScoreCut,
        trkDistCut=trkDistCut,
        cutClearCosmic=cutClearCosmic,
        **trkArgs,
    )


def make_pandora_evtdf_all_mc_multisim(
    f,
    sel_level="all",
    include_weights=True,
    multisim_nuniv=1000,
    wgt_types=None,
    slim=False,
    genie_systematics=None,
    flux_systematics=None,
    trkScoreCut=False,
    trkDistCut=1000.0,
    cutClearCosmic=False,
    **trkArgs,
):
    """Loose ``sel_all`` selection + configurable multisim weights (matches ``sel_all-mc`` cuts)."""
    if wgt_types is None:
        raise ValueError(
            "make_pandora_evtdf_all_mc_multisim requires wgt_types, e.g. ['genie']"
        )
    return make_pandora_evtdf(
        f,
        sel_level=sel_level,
        include_weights=include_weights,
        multisim_nuniv=multisim_nuniv,
        wgt_types=wgt_types,
        slim=slim,
        genie_systematics=genie_systematics,
        flux_systematics=flux_systematics,
        trkScoreCut=trkScoreCut,
        trkDistCut=trkDistCut,
        cutClearCosmic=cutClearCosmic,
        **trkArgs,
    )


def make_pandora_evtdf_all_genieslimwgts(
    f,
    sel_level="all",
    include_weights=True,
    genie_multisim_nuniv=100,
    wgt_types=None,
    slim=True,
    genie_systematics=None,
    trkScoreCut=False,
    trkDistCut=1000.0,
    cutClearCosmic=False,
    **trkArgs,
):
    """Loose ``sel_all`` + slim GENIE weights (``mc.GENIE.univ_*`` product; see ``getsyst.slim``)."""
    if wgt_types is None:
        wgt_types = ["genie"]
    if genie_systematics is None:
        from makedf.geniesyst import regen_systematics

        genie_systematics = regen_systematics
    return make_pandora_evtdf(
        f,
        sel_level=sel_level,
        include_weights=include_weights,
        genie_multisim_nuniv=genie_multisim_nuniv,
        wgt_types=wgt_types,
        slim=slim,
        genie_systematics=genie_systematics,
        trkScoreCut=trkScoreCut,
        trkDistCut=trkDistCut,
        cutClearCosmic=cutClearCosmic,
        **trkArgs,
    )


def make_mcnudf_mc_multisim(
    f,
    include_weights=True,
    multisim_nuniv=100,
    genie_multisim_nuniv=None,
    wgt_types=None,
    slim=False,
    genie_systematics=None,
    flux_systematics=None,
):
    if wgt_types is None:
        raise ValueError("make_mcnudf_mc_multisim requires wgt_types")
    if genie_multisim_nuniv is None:
        genie_multisim_nuniv = multisim_nuniv
    return make_mcnudf(
        f,
        include_weights=include_weights,
        multisim_nuniv=multisim_nuniv,
        genie_multisim_nuniv=genie_multisim_nuniv,
        wgt_types=wgt_types,
        slim=slim,
        genie_systematics=genie_systematics,
        flux_systematics=flux_systematics,
    )


def build_flux_knobgroup_config(group_filter=None):
    """
    Build DFS / ARGS / NAMES for BNB flux multisim (all knobs in one pass).

    HDF keys evt, hdr — same layout as sel_mup-g4wgts.py.

    Per-knob category bundles for *other* tooling live in makedf.bnbsyst.BNB_FLUX_GROUPS;
    dataframe production no longer splits outputs by those groups.

    group_filter: deprecated. If set, raises ValueError (use thin flux configs under
    configs/ if you need a subset of knobs in isolation).
    """
    if group_filter is not None:
        raise ValueError(
            "group_filter / FLUX_GROUP is no longer supported for flux df configs: "
            "all BNB flux multisim knobs are written to a single evt table. "
            "Unset FLUX_GROUP, or use a dedicated sel_mup-wgts_flux_*.py config."
        )

    from makedf.bnbsyst import regen_systematics

    evt_kw = dict(
        include_weights=True,
        wgt_types=["bnb"],
        slim=False,
        multisim_nuniv=1000,
        flux_systematics=regen_systematics,
        trkScoreCut=False,
        trkDistCut=100.0,
        cutClearCosmic=True,
    )

    DFS = [make_pandora_evtdf_mup_mc_multisim, make_hdrdf]
    ARGS = [evt_kw, {}]
    NAMES = ["evt", "hdr"]
    assert len(DFS) == len(ARGS) == len(NAMES)
    return DFS, ARGS, NAMES


def build_flux_knobgroup_config_sel_all(group_filter=None):
    """
    Build DFS / ARGS / NAMES for loose ``sel_all`` + BNB flux multisim.

    HDF keys evt, trk, hdr — required by ``syst_multisim_chunk.py`` with
    ``--input-stage sel_all``. Same weight knobs as ``build_flux_knobgroup_config``.

    group_filter: deprecated; if set, raises ValueError (see mup flux builder).
    """
    if group_filter is not None:
        raise ValueError(
            "group_filter / FLUX_GROUP is no longer supported for flux df configs: "
            "all BNB flux multisim knobs are written to a single evt table. "
            "Unset FLUX_GROUP, or use a dedicated sel_all-wgts_flux_*.py config."
        )

    from makedf.bnbsyst import regen_systematics

    evt_kw = dict(
        include_weights=True,
        wgt_types=["bnb"],
        slim=False,
        multisim_nuniv=1000,
        flux_systematics=regen_systematics,
    )

    DFS = [make_pandora_evtdf_all_mc_multisim, make_trkdf, make_hdrdf]
    ARGS = [evt_kw, {}, {}]
    NAMES = ["evt", "trk", "hdr"]
    assert len(DFS) == len(ARGS) == len(NAMES)
    return DFS, ARGS, NAMES


# --- GENIE knob groups: registry-driven configs (see build_genie_knobgroup_config) ---


def build_genie_knobgroup_config(group_filter=None):
    """
    Build DFS / ARGS / NAMES for GENIE multisim grouped by physics knob set.

    group_filter: None  -> all groups in GENIE_KNOB_GROUPS (HDF keys evt_<Name>, mcnu_<Name>, ...)
                  str   -> single group (HDF keys evt, mcnu, hdr for backward compatibility)

    Knob lists live in makedf.geniesyst (GENIE_KNOB_GROUPS).
    """
    from makedf.geniesyst import GENIE_KNOB_GROUPS

    if group_filter is not None:
        if group_filter not in GENIE_KNOB_GROUPS:
            raise ValueError(
                "Unknown GENIE knob group %r; valid keys: %s"
                % (group_filter, tuple(GENIE_KNOB_GROUPS))
            )
        groups = {group_filter: GENIE_KNOB_GROUPS[group_filter]}
    else:
        groups = GENIE_KNOB_GROUPS

    single = len(groups) == 1
    evt_kw = dict(
        include_weights=True,
        multisim_nuniv=200,
        wgt_types=["genie"],
        slim=False,
    )
    mcn_kw = dict(
        include_weights=True,
        multisim_nuniv=100,
        genie_multisim_nuniv=100,
        wgt_types=["genie"],
        slim=False,
    )

    DFS, ARGS, NAMES = [], [], []
    for name, syst in groups.items():
        DFS.append(make_pandora_evtdf_mup_mc_multisim)
        ARGS.append({**evt_kw, "genie_systematics": syst})
        NAMES.append("evt" if single else "evt_%s" % name)

        DFS.append(make_mcnudf_mc_multisim)
        ARGS.append({**mcn_kw, "genie_systematics": syst})
        NAMES.append("mcnu" if single else "mcnu_%s" % name)

    DFS.append(make_hdrdf)
    ARGS.append({})
    NAMES.append("hdr")

    assert len(DFS) == len(ARGS) == len(NAMES)
    return DFS, ARGS, NAMES


def build_genie_slim_config_sel_all(
    genie_systematics=None,
    genie_multisim_nuniv=100,
):
    """
    Build DFS / ARGS / NAMES for loose ``sel_all`` + **slim** GENIE weights.

    Multisim knobs (CAF type 0) are multiplied universe-by-universe into
    ``mc.GENIE.univ_*``; multisigma / morph knobs stay as separate per-knob columns
    (see ``makedf.getsyst.getsyst(..., slim=True)`` and ``geniesyst._slim_genie_weight_columns``).

    Default knob list: ``regen_systematics`` (Spring regen CAFs). Pass e.g.
    ``gen1_systematics`` for legacy gen1 samples, or extend with Ar23p knobs on
    respin CAFs.
    """
    if genie_systematics is None:
        from makedf.geniesyst import regen_systematics

        genie_systematics = regen_systematics

    evt_kw = dict(
        include_weights=True,
        genie_multisim_nuniv=genie_multisim_nuniv,
        wgt_types=["genie"],
        slim=True,
        genie_systematics=genie_systematics,
    )
    mcnu_kw = dict(
        genie_systematics=genie_systematics,
        genie_multisim_nuniv=genie_multisim_nuniv,
        slim=True,
    )

    DFS = [
        make_pandora_evtdf_all_genieslimwgts,
        make_trkdf,
        make_mcnudf_genieslimwgts,
        make_hdrdf,
    ]
    ARGS = [evt_kw, {}, mcnu_kw, {}]
    NAMES = ["evt", "trk", "mcnu", "hdr"]
    assert len(DFS) == len(ARGS) == len(NAMES)
    return DFS, ARGS, NAMES


def build_genie_knobgroup_config_sel_all(group_filter=None):
    """
    Build DFS / ARGS / NAMES for loose ``sel_all`` + GENIE multisim by knob group.

    Output layout (single group via ``GENIE_KNOB_GROUP``):
      ``evt``, ``trk``, ``mcnu``, ``hdr`` — required by ``get_systematics_genie.py``
      with ``--input-stage sel_all``.

    Multi-group (``group_filter=None``): ``evt_<Name>``, ``mcnu_<Name>`` per group,
    plus one shared ``trk`` and ``hdr``.
    """
    from makedf.geniesyst import GENIE_KNOB_GROUPS

    if group_filter is not None:
        if group_filter not in GENIE_KNOB_GROUPS:
            raise ValueError(
                "Unknown GENIE knob group %r; valid keys: %s"
                % (group_filter, tuple(GENIE_KNOB_GROUPS))
            )
        groups = {group_filter: GENIE_KNOB_GROUPS[group_filter]}
    else:
        groups = GENIE_KNOB_GROUPS

    single = len(groups) == 1
    evt_kw = dict(
        include_weights=True,
        multisim_nuniv=200,
        wgt_types=["genie"],
        slim=False,
    )
    mcn_kw = dict(
        include_weights=True,
        multisim_nuniv=100,
        genie_multisim_nuniv=100,
        wgt_types=["genie"],
        slim=False,
    )

    DFS, ARGS, NAMES = [], [], []
    for name, syst in groups.items():
        DFS.append(make_pandora_evtdf_all_mc_multisim)
        ARGS.append({**evt_kw, "genie_systematics": syst})
        NAMES.append("evt" if single else "evt_%s" % name)

        if single:
            DFS.append(make_trkdf)
            ARGS.append({})
            NAMES.append("trk")

        DFS.append(make_mcnudf_mc_multisim)
        ARGS.append({**mcn_kw, "genie_systematics": syst})
        NAMES.append("mcnu" if single else "mcnu_%s" % name)

    if not single:
        DFS.append(make_trkdf)
        ARGS.append({})
        NAMES.append("trk")

    DFS.append(make_hdrdf)
    ARGS.append({})
    NAMES.append("hdr")

    assert len(DFS) == len(ARGS) == len(NAMES)
    return DFS, ARGS, NAMES

# ================================================


# ===== Calo / E-field variations =====
# Configs should pass ``updatecalo=...`` / ``updateefield=...`` kwargs (see
# ``configs/numucc_1p0pi/sel_mup-updatecalo.py``). Named wrappers below are thin
# aliases kept for older job scripts.

def _pandora_updatecalo_maker(sel_level, updatecalo="CV", updateefield=False):
    def _maker(f, sel_level=sel_level, include_weights=False, multisim_nuniv=0, wgt_types=[],
               slim=True, trkScoreCut=False, trkDistCut=100., cutClearCosmic=True,
               updatecalo=updatecalo, updateefield=updateefield, **trkArgs):
        return make_pandora_evtdf(
            f, sel_level=sel_level, include_weights=include_weights,
            multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim,
            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic,
            updatecalo=updatecalo, updateefield=updateefield, **trkArgs,
        )
    return _maker


def _trkdf_updatecalo_maker(updatecalo=True):
    def _maker(f, trkScoreCut=False, trkDistCut=100., updatecalo=updatecalo, **trkArgs):
        return make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, **trkArgs)
    return _maker


make_pandora_evtdf_mup_updateefield = _pandora_updatecalo_maker("mup", updatecalo="CV", updateefield=True)
make_pandora_evtdf_2prong_updateefield = _pandora_updatecalo_maker("2prong", updatecalo="CV", updateefield=True)
make_pandora_evtdf_mup_updatecalo = _pandora_updatecalo_maker("mup", updatecalo="CV")
make_pandora_evtdf_2prong_updatecalo = _pandora_updatecalo_maker("2prong", updatecalo="CV")
make_pandora_evtdf_2prong_vtxdist_updatecalo = _pandora_updatecalo_maker("2prong_vtxdist", updatecalo="CV")
make_pandora_evtdf_2prong_wcandidates_updatecalo = _pandora_updatecalo_maker("2prong_wcandidates", updatecalo="CV")
make_pandora_evtdf_all_updatecalo = _pandora_updatecalo_maker("all", updatecalo="CV")
make_trkdf_updatecalo = _trkdf_updatecalo_maker(True)

# Legacy per-knob aliases (prefer ARGS={"updatecalo": ...} on the generic makers).
for _calo_tag in ("ccal_p", "ccal_m", "alpha_p", "alpha_m", "beta_p", "beta_m", "R_p", "R_m"):
    globals()[f"make_pandora_evtdf_mup_updatecalo_{_calo_tag}"] = _pandora_updatecalo_maker(
        "mup", updatecalo=_calo_tag
    )
    globals()[f"make_trkdf_updatecalo_{_calo_tag}"] = _trkdf_updatecalo_maker(_calo_tag)


# for SystVar samples
def make_metadf(f):
    mcdf = make_mcnudf(f, include_weights=False)
    metabranches = ["rec.hdr.pot",
                     "rec.hdr.nbnbinfo",
                     "rec.hdr.first_in_subrun",
                     "rec.hdr.ismc",
                     "rec.hdr.run",
                     "rec.hdr.subrun",
                     "rec.hdr.ngenevt",
                     "rec.hdr.evt"]
    hdrdf = loadbranches(f["recTree"], metabranches).rec.hdr
    df = multicol_merge(mcdf.reset_index(), hdrdf.reset_index(), 
                            left_on=["entry"], right_on=["entry"], 
                            how="left")
    return df

# ================================================

def make_pandora_evtdf(f, sel_level="all", 
                       include_weights=True, multisim_nuniv=1000, genie_multisim_nuniv=100, wgt_types=[], slim=True, genie_systematics=None, flux_systematics=None,
                       trkScoreCut=False, trkDistCut=100., updatecalo=None, updateefield=False,
                       cutClearCosmic=True, **trkArgs):

    """Build a pandora event dataframe at the requested selection depth.

    Cuts come from ``event_selection_pipeline_def.build_pipeline`` (same chain as
    the batched notebook and syst walker). Thresholds live in ``makedf/selections.py``.

    sel_level (CAF product names → pipeline stages via ``CAF_SEL_LEVEL_TO_STAGE``):
        "all": all slices, no cuts
        "clearcosmic": cosmic rejection
        "fv": vertex in Gen-1 FV
        "nu": nu-score cut
        "2prong" / "2prong_contained" / "2prong_trackscore" / "2prong_vtxdist":
            mid-selection products
        "2prong_wcandidates": after vtxdist + μ/p candidate columns (no has_μ/p cuts)
        "muX": muon candidate + muon kinematics
        "mup": final μ+p selection (+ reco TKI)
    """

    if sel_level not in CAF_SEL_LEVEL_TO_STAGE:
        raise ValueError(
            f"Invalid sel_level: {sel_level!r}; "
            f"expected one of {sorted(CAF_SEL_LEVEL_TO_STAGE)}"
        )

    def truth_match(this_evtdf, this_mcdf):
        # ---- truth match ----
        bad_tmatch = np.invert(this_evtdf.slc.tmatch.eff > 0.5) & (this_evtdf.slc.tmatch.idx >= 0)
        this_evtdf.loc[bad_tmatch, pad_column_name(("slc","tmatch","idx"), this_evtdf)] = np.nan

        nlevels = this_evtdf.columns.nlevels

        this_mcdf.columns = pd.MultiIndex.from_tuples([tuple(["mc"] + list(c) +[""] * (nlevels-len(c)-1)) for c in this_mcdf.columns])     # match # of column levels
        df = multicol_merge(this_evtdf.reset_index(), 
                    this_mcdf.reset_index(),
                    left_on=[("entry", "", "",), 
                            ("slc", "tmatch", "idx")], 
                    right_on=[("entry", "", ""), 
                                ("rec.mc.nu..index", "", "")], 
                    how="left"
                    ) 

        df = df.set_index(this_evtdf.index.names, verify_integrity=True) 
        return df

    mcdf = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, genie_systematics=genie_systematics, flux_systematics=flux_systematics)

    # calculate TKI for MC (attached before truth-match)
    tki_var_names = ["del_alpha", "del_phi", "del_Tp", "del_p", "del_Tp_x", "del_Tp_y"]
    mc_mudf = mcdf.mu
    mc_pdf = mcdf.p
    mc_P_mu_col = pad_column_name(("totp",), mc_mudf)
    mc_P_p_col = pad_column_name(("totp",), mc_pdf)
    tki_mc = get_cc1p0pi_tki(mc_mudf, mc_pdf, mc_P_mu_col, mc_P_p_col)
    for var_name in tki_var_names:
        mcdf = multicol_add(mcdf, tki_mc[var_name].rename("{}".format(var_name)))

    slcdf = make_slcdf(f)
    stop_at = CAF_SEL_LEVEL_TO_STAGE[sel_level]
    if stop_at is None:
        return truth_match(slcdf, mcdf)

    trkdf = None
    if sel_level in CAF_SEL_LEVELS_NEED_TRACKS:
        trkdf = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, updateefield=updateefield, **trkArgs)
        trkdf = get_valid_trks(trkdf)

    # Calorimetry variations write chi2 columns with a "_new" suffix; thresholds
    # still come from selections.get_mu_p_candidate defaults.
    mu_p_candidate_kwargs = (
        {"score_tag": "_new"} if updatecalo is not None else None
    )

    state = apply_selection_pipeline(
        {"evt": slcdf, "trk": trkdf, "hdr": None},
        stop_at=stop_at,
        sample="mc",
        mu_p_candidate_kwargs=mu_p_candidate_kwargs,
    )
    evtdf = state["evt"]

    # CAF-only mid product: μ/p candidate columns without has_μ / has_p / kinematics.
    if sel_level == "2prong_wcandidates":
        evtdf = get_mu_p_candidate(evtdf, **(mu_p_candidate_kwargs or {}))

    return truth_match(evtdf, mcdf)


# ===========================================================================
# Systematic histogram counts (grid): selection + xsec variables per knob/univ
# ===========================================================================
# Instead of writing heavy weight tables and re-histogramming offline, these
# makers walk the event-selection pipeline on the CAF and store long-format
# bin counts under HDF key ``syst_hists`` (see ``syst_histcounts.py``).
#
# GENIE: rate histcounts + xsec response tensors (rate ≠ xsec — see module doc).
# Flux / G4: rate histcounts only (multisim).
# WireMod / DENT / intime / offbeam: CV counts only (unisim; pair samples in notebook).


def make_var_config_snapshot(f, **_kwargs):
    """Freeze VariableConfigs used by histcounts into HDF key ``var_configs``.

    Independent of the CAF contents so every job output carries the binning
    that was live at fill / submit time (notebook must not import the live module).
    """
    from analysis_village.numucc_1p0pi.syst_histcounts import histcounts_var_configs_df

    return histcounts_var_configs_df()


def make_syst_histcounts(
    f,
    wgt_types=None,
    family=None,
    multisim_nuniv=200,
    genie_multisim_nuniv=100,
    slim=False,
    include_slim=True,
    genie_systematics=None,
    flux_systematics=None,
    knob_names=None,
    sample="mc",
    do_xsec=None,
    trkScoreCut=False,
    trkDistCut=1000.0,
    cutClearCosmic=False,
    **trkArgs,
):
    """Load CAF with weights, run event selection, return ``syst_hists`` DataFrame.

    Parameters
    ----------
    wgt_types : list
        e.g. ``['genie']``, ``['bnb']``, ``['g4']``. Required for weighted modes.
    family : str or None
        ``GENIE`` / ``Flux`` / ``G4``. Inferred from ``wgt_types`` if omitted.
    slim : bool
        If True, load via ``getsyst(..., slim=True)`` (product + leftover ±σ/morph).
    include_slim : bool
        If True (default), also histogram slim products:
        * ``slim_multisim`` (GENIE) / ``Flux_slim_multisim`` / ``G4_slim_multisim`` —
          product of **true multisim** knobs only.
        * ``slim`` (GENIE) / ``Flux_slim`` / ``G4_slim`` — slim_multisim × Gaussian
          throws of multisigma (``ps1``) and morph (see notebook recipe).
        Do not sum both with per-knob covs (double-counting).
    knob_names : sequence or None
        Restrict to these knobs under ``mc``; default = auto-detect on the frame.
    do_xsec : bool or None
        Store GENIE xsec tensors. Default True only for ``family=='GENIE'``.
    """
    from analysis_village.numucc_1p0pi.syst_histcounts import (
        attach_family_slim_products,
        discover_syst_names_on_df,
        empty_histcounts_df,
        fill_syst_histcounts,
    )

    if wgt_types is None:
        raise ValueError("make_syst_histcounts requires wgt_types, e.g. ['genie']")
    wgt_types = list(wgt_types)
    if family is None:
        if "genie" in wgt_types:
            family = "GENIE"
        elif "bnb" in wgt_types:
            family = "Flux"
        elif "g4" in wgt_types:
            family = "G4"
        else:
            family = "SYST"
    if do_xsec is None:
        do_xsec = family.upper() == "GENIE"

    evt_kw = dict(
        include_weights=True,
        multisim_nuniv=multisim_nuniv,
        wgt_types=wgt_types,
        slim=slim,
        genie_systematics=genie_systematics,
        flux_systematics=flux_systematics,
        trkScoreCut=trkScoreCut,
        trkDistCut=trkDistCut,
        cutClearCosmic=cutClearCosmic,
    )
    if "genie" in wgt_types:
        evt_kw["genie_multisim_nuniv"] = genie_multisim_nuniv

    try:
        if slim and "genie" in wgt_types:
            evt = make_pandora_evtdf_all_genieslimwgts(
                f,
                genie_multisim_nuniv=genie_multisim_nuniv,
                genie_systematics=genie_systematics,
                trkScoreCut=trkScoreCut,
                trkDistCut=trkDistCut,
                cutClearCosmic=cutClearCosmic,
                **trkArgs,
            )
        else:
            evt = make_pandora_evtdf_all_mc_multisim(f, **evt_kw, **trkArgs)

        trk = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, **trkArgs)

        mcnu = None
        if do_xsec or "genie" in wgt_types:
            if slim and "genie" in wgt_types:
                mcnu = make_mcnudf_genieslimwgts(
                    f,
                    genie_multisim_nuniv=genie_multisim_nuniv,
                    genie_systematics=genie_systematics,
                    slim=True,
                )
            else:
                mcnu = make_mcnudf_mc_multisim(
                    f,
                    include_weights=True,
                    multisim_nuniv=multisim_nuniv,
                    genie_multisim_nuniv=genie_multisim_nuniv,
                    wgt_types=wgt_types,
                    slim=slim,
                    genie_systematics=genie_systematics,
                    flux_systematics=flux_systematics,
                )
    except Exception as ex:
        # Common on CAFs missing the requested weight branches (empty wgtdf → MultiIndex error).
        print(
            "[make_syst_histcounts] failed building evt/trk/mcnu with weights (%s); "
            "returning empty syst_hists" % ex
        )
        return empty_histcounts_df()


    if evt is None or len(evt) == 0:
        return empty_histcounts_df()

    n_slim = int(genie_multisim_nuniv) if str(family).upper() == "GENIE" else int(multisim_nuniv)
    want_knobs = list(knob_names) if knob_names is not None else None
    if include_slim:
        evt, slim_names = attach_family_slim_products(
            evt, family=family, n_univ=n_slim, knob_names=want_knobs
        )
        if mcnu is not None and len(mcnu) > 0:
            mcnu, _ = attach_family_slim_products(
                mcnu, family=family, n_univ=n_slim, knob_names=want_knobs
            )
        if want_knobs is not None:
            for sn in slim_names:
                if sn not in want_knobs:
                    want_knobs.append(sn)
        # Auto-detect path: slim products are already on the frame.

    syst_names = discover_syst_names_on_df(evt, family=family, knob_names=want_knobs)
    if not syst_names:
        print(
            "[make_syst_histcounts] no %s weight knobs found on evt (requested %r); "
            "returning empty syst_hists" % (family, want_knobs)
        )
        return empty_histcounts_df()

    try:
        return fill_syst_histcounts(
            evt,
            trk,
            mcnu,
            family=family,
            syst_names=syst_names,
            sample=sample,
            do_xsec=bool(do_xsec),
        )
    except Exception as ex:
        print("[make_syst_histcounts] fill failed (%s); returning empty syst_hists" % ex)
        return empty_histcounts_df()



def make_syst_histcounts_nowgt(
    f,
    sample="mc",
    trkScoreCut=False,
    trkDistCut=1000.0,
    cutClearCosmic=False,
    **trkArgs,
):
    """Unisim path: no weights — CV histogram counts only (WireMod / DENT / intime / offbeam).

    Pair CV vs variation sample counts in the notebook via
    ``syst_histcounts.unisim_cov_from_cv_and_var``.
    """
    from analysis_village.numucc_1p0pi.syst_histcounts import (
        empty_histcounts_df,
        fill_syst_histcounts,
    )

    evt = make_pandora_evtdf_all(
        f,
        include_weights=False,
        multisim_nuniv=0,
        wgt_types=[],
        slim=True,
        trkScoreCut=trkScoreCut,
        trkDistCut=trkDistCut,
        cutClearCosmic=cutClearCosmic,
        **trkArgs,
    )
    trk = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, **trkArgs)
    # Truth categories for signal/background bookkeeping when MC/dirt/intime.
    mcnu = None
    try:
        mcnu = make_mcnudf(f, include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True)
    except Exception:
        mcnu = None

    if evt is None or len(evt) == 0:
        return empty_histcounts_df()

    return fill_syst_histcounts(
        evt,
        trk,
        mcnu,
        family="UNISIM",
        syst_names=[],
        sample=sample,
        do_xsec=False,
    )


def make_syst_histcounts_all(
    f,
    sample="mc",
    include_ar23p=True,
    include_slim=True,
    genie_multisim_nuniv=100,
    flux_multisim_nuniv=1000,
    g4_multisim_nuniv=1000,
    trkScoreCut=False,
    trkDistCut=1000.0,
    cutClearCosmic=False,
    **trkArgs,
):
    """One CAF pass: all GENIE knobs + slim products + Flux + G4 histcounts.

    Histcounts are stored **per individual knob** (e.g.
    ``GENIEReWeight_SBN_v1_multisim_CoulombCCQE``,
    ``GENIEReWeight_SBN_v1_multisim_NormCCMEC``), never under mode labels
    like ``CCQE`` / ``MEC``.

    With ``include_slim=True`` (default), also attaches and histograms:

    * ``slim_multisim`` / ``Flux_slim_multisim`` / ``G4_slim_multisim`` —
      product of **true multisim** (CAF type 0) knobs only.
    * ``slim`` / ``Flux_slim`` / ``G4_slim`` — slim_multisim × Gaussian throws of
      multisigma (``ps1``) and morph (notebook / historical getsyst recipe;
      weights clipped ``≥ 0``).

    Per-knob ±σ/morph columns remain available separately. Do **not** sum slim
    (or slim_multisim) fractional cov on top of the sum of per-knob multisim covs
    in the notebook (double-counting). ``slim`` already includes ``slim_multisim``.

    ``include_ar23p``: include Ar23p template knobs (default True).
    """
    from makedf.bnbsyst import regen_systematics
    from makedf.g4syst import g4_systematics
    from makedf.geniesyst import GENIE_KNOB_GROUPS

    from analysis_village.numucc_1p0pi.syst_histcounts import (
        attach_family_slim_products,
        discover_syst_names_on_df,
        empty_histcounts_df,
        fill_syst_histcounts,
        sum_histcounts_dfs,
    )

    # Flatten group registry → per-knob CAF weight names (not group labels).
    genie_syst = []
    genie_knobs = []
    for name, lst in GENIE_KNOB_GROUPS.items():
        if name == "Ar23p" and not include_ar23p:
            continue
        genie_syst.extend(lst)
        genie_knobs.extend(lst)
    flux_knobs = list(regen_systematics)
    g4_knobs = [k for k in g4_systematics if "neutron" not in k]

    try:
        # slim=False keeps every per-knob column; we build slim products ourselves so
        # multisigma/morph stay separate and never enter the multisim product.
        evt = make_pandora_evtdf_all_mc_multisim(
            f,
            include_weights=True,
            multisim_nuniv=max(int(flux_multisim_nuniv), int(g4_multisim_nuniv), 200),
            wgt_types=["genie", "bnb", "g4"],
            slim=False,
            genie_systematics=genie_syst,
            flux_systematics=regen_systematics,
            trkScoreCut=trkScoreCut,
            trkDistCut=trkDistCut,
            cutClearCosmic=cutClearCosmic,
            genie_multisim_nuniv=genie_multisim_nuniv,
            **trkArgs,
        )
        trk = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, **trkArgs)
        mcnu = make_mcnudf_mc_multisim(
            f,
            include_weights=True,
            multisim_nuniv=200,
            genie_multisim_nuniv=genie_multisim_nuniv,
            wgt_types=["genie", "bnb", "g4"],
            slim=False,
            genie_systematics=genie_syst,
            flux_systematics=regen_systematics,
        )
    except Exception as ex:
        print(
            "[make_syst_histcounts_all] failed building weighted tables (%s); "
            "returning empty syst_hists" % ex
        )
        return empty_histcounts_df()

    if evt is None or len(evt) == 0:
        return empty_histcounts_df()

    if include_slim:
        evt, genie_slim = attach_family_slim_products(
            evt, family="GENIE", n_univ=int(genie_multisim_nuniv), knob_names=genie_knobs
        )
        evt, flux_slim = attach_family_slim_products(
            evt, family="Flux", n_univ=int(flux_multisim_nuniv), knob_names=flux_knobs
        )
        evt, g4_slim = attach_family_slim_products(
            evt, family="G4", n_univ=int(g4_multisim_nuniv), knob_names=g4_knobs
        )
        if mcnu is not None and len(mcnu) > 0:
            mcnu, _ = attach_family_slim_products(
                mcnu, family="GENIE", n_univ=int(genie_multisim_nuniv), knob_names=genie_knobs
            )
            mcnu, _ = attach_family_slim_products(
                mcnu, family="Flux", n_univ=int(flux_multisim_nuniv), knob_names=flux_knobs
            )
            mcnu, _ = attach_family_slim_products(
                mcnu, family="G4", n_univ=int(g4_multisim_nuniv), knob_names=g4_knobs
            )
    else:
        genie_slim, flux_slim, g4_slim = [], [], []

    parts = []
    # GENIE: per-knob (+ optional slim products)
    genie_want = list(genie_knobs) + list(genie_slim)
    genie_names = discover_syst_names_on_df(evt, family="GENIE", knob_names=genie_want)
    if genie_names:
        parts.append(
            fill_syst_histcounts(
                evt, trk, mcnu, family="GENIE", syst_names=genie_names, sample=sample, do_xsec=True
            )
        )
    else:
        print("[make_syst_histcounts_all] no GENIE knobs found on evt")

    # Flux / G4: rate only
    flux_want = list(flux_knobs) + list(flux_slim)
    flux_names = discover_syst_names_on_df(evt, family="Flux", knob_names=flux_want)
    if flux_names:
        parts.append(
            fill_syst_histcounts(
                evt, trk, None, family="Flux", syst_names=flux_names, sample=sample, do_xsec=False
            )
        )
    else:
        print("[make_syst_histcounts_all] no Flux knobs found on evt")

    g4_want = list(g4_knobs) + list(g4_slim)
    g4_names = discover_syst_names_on_df(evt, family="G4", knob_names=g4_want)
    if g4_names:
        parts.append(
            fill_syst_histcounts(
                evt, trk, None, family="G4", syst_names=g4_names, sample=sample, do_xsec=False
            )
        )
    else:
        print("[make_syst_histcounts_all] no G4 knobs found on evt")

    return sum_histcounts_dfs(parts) if parts else empty_histcounts_df()


def build_syst_histcounts_config(
    mode=None,
    group_filter=None,
    sample="mc",
):
    """Build ``DFS, ARGS, NAMES`` for grid histcount jobs.

    ``mode`` (or env ``SYST_HIST_MODE``):
      ``all`` | ``genie`` | ``genie_slim`` | ``flux`` | ``g4`` | ``nowgt``.

    ``all`` — every GENIE knob in ``GENIE_KNOB_GROUPS`` (including Ar23p) + Flux + G4,
    plus slim products ``slim_multisim`` / ``slim`` (and Flux/G4 analogues) by default.
    Histcounts are keyed by **individual knob name**, not by group (CCQE/MEC/…).
    Set ``SYST_HIST_EXCLUDE_AR23P=1`` to drop Ar23p knobs only.
    Set ``SYST_HIST_EXCLUDE_SLIM=1`` to skip slim product histcounts.

    ``group_filter`` (or env ``GENIE_KNOB_GROUP``): optional GENIE *group* used only
    to select which knobs to load when ``mode=genie``; output is still one entry
    per knob inside that group (+ ``slim_multisim`` / ``slim`` unless excluded).
    """
    import os

    if mode is None:
        mode = os.environ.get("SYST_HIST_MODE", "genie").strip().lower()
    if group_filter is None:
        group_filter = os.environ.get("GENIE_KNOB_GROUP", "").strip() or None
    # Ar23p off by default for Spring CV histcounts; set SYST_HIST_EXCLUDE_AR23P=0
    # or GENIE_KNOB_GROUP=Ar23p on the AR23plus CAF sample to load those knobs.
    _raw_excl = os.environ.get("SYST_HIST_EXCLUDE_AR23P", "1").strip().lower()
    exclude_ar23p = _raw_excl not in ("0", "false", "no")
    include_ar23p = not exclude_ar23p
    exclude_slim = os.environ.get("SYST_HIST_EXCLUDE_SLIM", "").strip() in (
        "1",
        "true",
        "True",
        "yes",
        "YES",
    )
    include_slim = not exclude_slim

    def _env_int(name, default):
        raw = os.environ.get(name, "").strip()
        if not raw:
            return int(default)
        return int(raw)

    # Smoke / override universe counts (defaults match production).
    genie_nuniv = _env_int("SYST_HIST_GENIE_NUNIV", 100)
    flux_nuniv = _env_int("SYST_HIST_FLUX_NUNIV", 1000)
    g4_nuniv = _env_int("SYST_HIST_G4_NUNIV", 1000)

    def _with_var_config_snapshot(DFS, ARGS, NAMES):
        return (
            list(DFS) + [make_var_config_snapshot],
            list(ARGS) + [{}],
            list(NAMES) + ["var_configs"],
        )

    if mode in ("nowgt", "unisim", "none"):
        DFS = [make_syst_histcounts_nowgt, make_hdrdf]
        ARGS = [dict(sample=sample), {}]
        NAMES = ["syst_hists", "hdr"]
        return _with_var_config_snapshot(DFS, ARGS, NAMES)

    if mode in ("all", "full"):
        DFS = [make_syst_histcounts_all, make_hdrdf]
        ARGS = [
            dict(
                sample=sample,
                include_ar23p=include_ar23p,
                include_slim=include_slim,
                genie_multisim_nuniv=genie_nuniv,
                flux_multisim_nuniv=flux_nuniv,
                g4_multisim_nuniv=g4_nuniv,
            ),
            {},
        ]
        NAMES = ["syst_hists", "hdr"]
        return _with_var_config_snapshot(DFS, ARGS, NAMES)

    if mode in ("genie_slim", "slim"):
        from makedf.geniesyst import regen_systematics

        DFS = [make_syst_histcounts, make_hdrdf]
        ARGS = [
            dict(
                wgt_types=["genie"],
                family="GENIE",
                slim=True,
                include_slim=True,
                genie_systematics=regen_systematics,
                genie_multisim_nuniv=genie_nuniv,
                do_xsec=True,
                sample=sample,
            ),
            {},
        ]
        NAMES = ["syst_hists", "hdr"]
        return _with_var_config_snapshot(DFS, ARGS, NAMES)

    if mode == "genie":
        from makedf.geniesyst import GENIE_KNOB_GROUPS

        if group_filter is not None:
            if group_filter not in GENIE_KNOB_GROUPS:
                raise ValueError(
                    "Unknown GENIE knob group %r; valid: %s"
                    % (group_filter, tuple(GENIE_KNOB_GROUPS))
                )
            syst = list(GENIE_KNOB_GROUPS[group_filter])
            knobs = list(syst)
        else:
            syst = []
            knobs = []
            for name, lst in GENIE_KNOB_GROUPS.items():
                if name == "Ar23p" and not include_ar23p:
                    continue
                syst.extend(lst)
                knobs.extend(lst)
        DFS = [make_syst_histcounts, make_hdrdf]
        ARGS = [
            dict(
                wgt_types=["genie"],
                family="GENIE",
                slim=False,
                include_slim=include_slim,
                genie_systematics=syst,
                multisim_nuniv=200,
                genie_multisim_nuniv=genie_nuniv,
                knob_names=knobs,
                do_xsec=True,
                sample=sample,
            ),
            {},
        ]
        NAMES = ["syst_hists", "hdr"]
        return _with_var_config_snapshot(DFS, ARGS, NAMES)

    if mode == "flux":
        from makedf.bnbsyst import regen_systematics

        DFS = [make_syst_histcounts, make_hdrdf]
        ARGS = [
            dict(
                wgt_types=["bnb"],
                family="Flux",
                slim=False,
                include_slim=include_slim,
                multisim_nuniv=flux_nuniv,
                flux_systematics=regen_systematics,
                knob_names=list(regen_systematics),
                do_xsec=False,
                sample=sample,
            ),
            {},
        ]
        NAMES = ["syst_hists", "hdr"]
        return _with_var_config_snapshot(DFS, ARGS, NAMES)

    if mode == "g4":
        from makedf.g4syst import g4_systematics

        DFS = [make_syst_histcounts, make_hdrdf]
        ARGS = [
            dict(
                wgt_types=["g4"],
                family="G4",
                slim=False,
                include_slim=include_slim,
                multisim_nuniv=g4_nuniv,
                knob_names=[k for k in g4_systematics if "neutron" not in k],
                do_xsec=False,
                sample=sample,
            ),
            {},
        ]
        NAMES = ["syst_hists", "hdr"]
        return _with_var_config_snapshot(DFS, ARGS, NAMES)

    raise ValueError(
        "Unknown SYST_HIST_MODE %r; use all, genie, genie_slim, flux, g4, or nowgt" % (mode,)
    )
