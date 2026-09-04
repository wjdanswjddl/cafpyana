from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *

from analysis_village.numucc_1p0pi.makedf.selections import *
from analysis_village.numucc_1p0pi.categories import DETECTOR, TRK_CALO_DET, in_fv
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


# ===== Calo variations =====

def make_pandora_evtdf_mup_updateefield(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="CV", updateefield=True, cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, updatecalo=updatecalo, updateefield=updateefield, **trkArgs)
    return df

def make_pandora_evtdf_2prong_updateefield(f, sel_level="2prong", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="CV", updateefield=True, cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, updatecalo=updatecalo, updateefield=updateefield, **trkArgs)
    return df

def make_pandora_evtdf_mup_updatecalo(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="CV", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_pandora_evtdf_2prong_updatecalo(f, sel_level="2prong", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="CV", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df


def make_pandora_evtdf_2prong_vtxdist_updatecalo(f, sel_level="2prong_vtxdist", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="CV", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_pandora_evtdf_2prong_wcandidates_updatecalo(f, sel_level="2prong_wcandidates", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="CV", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_pandora_evtdf_all_updatecalo(f, sel_level="all", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="CV", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_trkdf_updatecalo(f, trkScoreCut=False, trkDistCut=100., updatecalo=True, **trkArgs):
    df = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, **trkArgs)
    return df

def make_pandora_evtdf_mup_updatecalo_ccal_p(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="ccal_p", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_trkdf_updatecalo_ccal_p(f, trkScoreCut=False, trkDistCut=100., updatecalo="ccal_p", **trkArgs):
    df = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, **trkArgs)
    return df

def make_pandora_evtdf_mup_updatecalo_ccal_m(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="ccal_m", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_trkdf_updatecalo_ccal_m(f, trkScoreCut=False, trkDistCut=100., updatecalo="ccal_m", **trkArgs):
    df = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, **trkArgs)
    return df

def make_pandora_evtdf_mup_updatecalo_alpha_p(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="alpha_p", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_trkdf_updatecalo_alpha_p(f, trkScoreCut=False, trkDistCut=100., updatecalo="alpha_p", **trkArgs):
    df = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, **trkArgs)
    return df

def make_pandora_evtdf_mup_updatecalo_alpha_m(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="alpha_m", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_trkdf_updatecalo_alpha_m(f, trkScoreCut=False, trkDistCut=100., updatecalo="alpha_m", **trkArgs):
    df = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, **trkArgs)
    return df

def make_pandora_evtdf_mup_updatecalo_beta_p(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="beta_p", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_trkdf_updatecalo_beta_p(f, trkScoreCut=False, trkDistCut=100., updatecalo="beta_p", **trkArgs):
    df = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, **trkArgs)
    return df

def make_pandora_evtdf_mup_updatecalo_beta_m(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="beta_m", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_trkdf_updatecalo_beta_m(f, trkScoreCut=False, trkDistCut=100., updatecalo="beta_m", **trkArgs):
    df = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, **trkArgs)
    return df

def make_pandora_evtdf_mup_updatecalo_R_p(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="R_p", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_trkdf_updatecalo_R_p(f, trkScoreCut=False, trkDistCut=100., updatecalo="R_p", **trkArgs):
    df = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, **trkArgs)
    return df

def make_pandora_evtdf_mup_updatecalo_R_m(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., updatecalo="R_m", cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, updatecalo=updatecalo, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_trkdf_updatecalo_R_m(f, trkScoreCut=False, trkDistCut=100., updatecalo="R_m", **trkArgs):
    df = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, **trkArgs)
    return df

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

    """
    sel_level:
        "all": all slices, no cuts
        "clearcosmic": cosmic rejection
        "fv": vertex in FV
        "nu": n-score cut
        "2prong": 2-prong slices
        "2prong_qual": 2-prong slices with quality cuts on the 2 prongs
        "muX": muon-X cuts
        "mup": final selection
    """

    if sel_level not in ["all", "clearcosmic", "fv", "nu", "2prong", "2prong_contained", "2prong_trackscore", "2prong_vtxdist",  "2prong_wcandidates", "muX", "mup"]:
        raise ValueError("Invalid sel_level: {}".format(sel_level))


    def truth_match(this_evtdf, this_mcdf):
        # ---- truth match ----
        bad_tmatch = np.invert(this_evtdf.slc.tmatch.eff > 0.5) & (this_evtdf.slc.tmatch.idx >= 0)
        # this_evtdf.loc[bad_tmatch, ("slc","tmatch","idx", "", "", "", "")] = np.nan
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

    # TODO: read from caf
    # event selection cuts (from ``makedf.selections``)
    nu_score_th = NU_SCORE_TH
    save_ntrks = SAVE_NTRKS
    trackscore_th = TRACKSCORE_TH
    vtxdist_th = VTXDIST_TH
    mu_chi2mu_th = MU_CHI2MU_TH
    mu_chi2p_th = MU_CHI2P_TH
    mu_len_th = MU_LEN_TH
    qual_th = QUAL_TH
    p_chi2p_th = P_CHI2P_TH
    p_len_th = P_LEN_TH
    mu_Plo_th = MU_PLO_TH
    mu_Phi_th = MU_PHI_TH
    p_Plo_th = P_PLO_TH
    p_Phi_th = P_PHI_TH

    mcdf = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, genie_systematics=genie_systematics, flux_systematics=flux_systematics)


    # calculate TKI for MC 
    tki_var_names = ["del_alpha", "del_phi", "del_Tp", "del_p", "del_Tp_x", "del_Tp_y"]
    mc_mudf = mcdf.mu
    mc_pdf = mcdf.p
    mc_P_mu_col = pad_column_name(("totp",), mc_mudf)
    mc_P_p_col = pad_column_name(("totp",), mc_pdf)
    tki_mc = get_cc1p0pi_tki(mc_mudf, mc_pdf, mc_P_mu_col, mc_P_p_col)
    for var_name in tki_var_names:
        mcdf = multicol_add(mcdf, tki_mc[var_name].rename("{}".format(var_name)))

    
    slcdf = make_slcdf(f)


    if sel_level == "all":
        return truth_match(slcdf, mcdf)

    slcdf = cut_clear_cosmic(slcdf)
    if sel_level == "clearcosmic":
        return truth_match(slcdf, mcdf)

    slcdf = cut_vertex_in_fv(slcdf, det=DETECTOR)
    if sel_level == "fv":
        return truth_match(slcdf, mcdf)

    slcdf = cut_nu_score(slcdf, nu_score_th)
    if sel_level == "nu":
        return truth_match(slcdf, mcdf)

    trkdf = make_trkdf(f, det=TRK_CALO_DET, scoreCut=trkScoreCut, updatecalo=updatecalo, updateefield=updateefield, **trkArgs)
    trkdf = get_valid_trks(trkdf)
    trkdf = match_trkdf_to_slcdf(trkdf, slcdf)
    evtdf = get_trk_info(slcdf, trkdf, save_ntrks)

    evtdf = cut_2prong(evtdf)
    if sel_level == "2prong":
        ret = truth_match(evtdf, mcdf)
        return ret

    evtdf = cut_2prong_contained(evtdf, det=DETECTOR)
    if sel_level == "2prong_contained":
        return truth_match(evtdf, mcdf)

    evtdf = cut_2prong_trackscore(evtdf, trackscore_th)
    if sel_level == "2prong_trackscore":
        return truth_match(evtdf, mcdf)

    evtdf = cut_2prong_vtxdist(evtdf, vtxdist_th)
    if sel_level == "2prong_vtxdist":
        return truth_match(evtdf, mcdf)

    if updatecalo is not None:
        evtdf = get_mu_p_candidate(evtdf, 
                                    mu_chi2mu_th=mu_chi2mu_th, mu_chi2p_th=mu_chi2p_th, mu_len_th=mu_len_th, qual_th=qual_th, 
                                    p_chi2mu_th=-1, p_chi2p_th=p_chi2p_th, p_len_th=p_len_th, score_tag="_new")

    else:
        evtdf = get_mu_p_candidate(evtdf, 
                                    mu_chi2mu_th=mu_chi2mu_th, mu_chi2p_th=mu_chi2p_th, mu_len_th=mu_len_th, qual_th=qual_th, 
                                    p_chi2mu_th=-1, p_chi2p_th=p_chi2p_th, p_len_th=p_len_th, score_tag="")

    evtdf = cut_2prong_vtxdist(evtdf, vtxdist_th)

    if sel_level == "2prong_wcandidates":
        return truth_match(evtdf, mcdf)

    evtdf = cut_has_mu(evtdf)
    evtdf = cut_mu_kinematics(evtdf, mu_Plo_th=mu_Plo_th, mu_Phi_th=mu_Phi_th)
    if sel_level == "muX":
        return truth_match(evtdf, mcdf)

    evtdf = cut_has_p(evtdf)
    evtdf = cut_p_kinematics(evtdf, p_Plo_th=p_Plo_th, p_Phi_th=p_Phi_th)

    # calculate TKI for reco slices
    slc_mudf = evtdf.mu.pfp.trk
    slc_pdf = evtdf.p.pfp.trk
    slc_P_mu_col = pad_column_name(("P", "p_muon"), slc_mudf)
    slc_P_p_col = pad_column_name(("P", "p_proton"), slc_pdf)
    tki_reco = get_cc1p0pi_tki(slc_mudf, slc_pdf, slc_P_mu_col, slc_P_p_col)
    for var_name in tki_var_names:
        evtdf = multicol_add(evtdf, tki_reco[var_name].rename(var_name))

    if sel_level == "mup":
        return truth_match(evtdf, mcdf)
