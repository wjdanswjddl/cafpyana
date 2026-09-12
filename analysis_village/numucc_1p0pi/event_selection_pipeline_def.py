"""Pipeline definition for the numuCC 1p0pi event selection.

This file is the SINGLE place to edit when you want to change **cut order**,
add or remove a plot, or follow a new variable through the efficiency curve.
Thresholds and cut formulas live in ``makedf/selections.py``.

Consumers (all import ``build_pipeline`` / ``apply_selection_pipeline``):

* CAF maker ``make_pandora_evtdf`` (maps ``sel_level`` → stage key)
* Batched / live event selection
* Interactive notebook (via ``SampleBundle.run_pipeline``)
* Systematics ``syst_pipeline_walker`` (cut-stage histograms)

Adding new things
-----------------
* New cut       : add a Stage(...) at the right point in ``build_pipeline``.
* New plot      : append a PlotSpec(...) to a stage's ``plots`` list.
* New eff. var  : extend ``CORE_SELECTED_EVT_VARIABLE_CONFIGS`` or the extras passed
                  into ``with_final_selected_evt_variables`` for ``EFFICIENCY_VARS`` below,
                  or add to ``analysis_village/numucc_1p0pi/final_selected_evt_vars.py``.
"""
from __future__ import annotations

import sys
from os import path
from typing import Any, Dict, Iterator, List, Optional, Callable, Tuple

import numpy as np
import pandas as pd

sys.path.append(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))

from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.makedf.selections import (
    cut_clear_cosmic, cut_vertex_in_fv, cut_nu_score, cut_2prong, cut_2prong_contained,
    cut_2prong_trackscore, cut_2prong_vtxdist, cut_has_mu, cut_has_p,
    cut_mu_kinematics, cut_p_kinematics, get_mu_p_candidate, get_valid_trks,
    get_trk_info, evt_has_trk1_trk2,
    NU_SCORE_TH, TRACKSCORE_TH, VTXDIST_TH, SAVE_NTRKS,
    MU_CHI2MU_TH, MU_CHI2P_TH, MU_LEN_TH, QUAL_TH, P_CHI2P_TH, P_LEN_TH,
    MU_PLO_TH, MU_PHI_TH, P_PLO_TH, P_PHI_TH,
)
from pyanalib.variable_calculator import (
    add_mc_cc1p0pi_tki_mcnu,
    add_reco_cc1p0pi_tki_evtdf,
    add_truth_cc1p0pi_tki_evtdf,
)
from makedf.util import avg_chi2, match_trkdf_to_slcdf

from analysis_village.numucc_1p0pi.selection_framework import (
    Stage, PlotSpec, ChunkRunner,
)


SAMPLES = ("mc", "data", "intime", "offbeam", "dirt")

# CAF ``sel_level`` (underscore / short names) → pipeline stage key (hyphenated).
# ``None`` means return before any cut (raw slices).
# ``2prong_wcandidates`` is CAF-only: stop after vtxdist PID columns, then attach
# μ/p candidates without applying has_μ / has_p / kinematics (handled in the maker).
CAF_SEL_LEVEL_TO_STAGE: Dict[str, Optional[str]] = {
    "all": None,
    "clearcosmic": "is_clear_cosmic",
    "fv": "vertex_in_fv",
    "nu": "nu_score",
    "2prong": "2prong",
    "2prong_contained": "2prong-contained",
    "2prong_trackscore": "2prong-trackscore",
    "2prong_vtxdist": "2prong-vtxdist",
    "2prong_wcandidates": "2prong-vtxdist",  # + get_mu_p_candidate in the maker
    "muX": "2prong-muX",
    "mup": "2prong-mup",
}

# Stages at / after which the CAF maker must load the track table.
# ``nu`` stops after the nu-score cut without attaching tracks (legacy CAF product).
CAF_SEL_LEVELS_NEED_TRACKS = frozenset({
    "2prong", "2prong_contained", "2prong_trackscore", "2prong_vtxdist",
    "2prong_wcandidates", "muX", "mup",
})


def apply_selection_pipeline(
    state: Dict[str, Any],
    *,
    stop_at: Optional[str] = None,
    sample: str = "mc",
    mu_p_candidate_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Apply ``build_pipeline()`` cuts to ``state`` until (and including) ``stop_at``.

    ``state`` is ``{"evt", "trk", "hdr", ...}``. Optional ``mu_p_candidate_kwargs``
    (e.g. ``score_tag="_new"`` for calorimetry variations) is forwarded to
    ``get_mu_p_candidate`` at the ``2prong-muX`` stage.

    If ``stop_at`` is ``None``, no cuts are applied (caller wants raw slices).
    """
    if stop_at is None:
        return state
    if mu_p_candidate_kwargs:
        state = dict(state)
        state["_mu_p_candidate_kwargs"] = dict(mu_p_candidate_kwargs)
    found = False
    for stage in build_pipeline():
        if stage.cut is not None:
            state = stage.cut(state, sample=sample)
        if stage.key == stop_at:
            found = True
            break
    if not found:
        raise ValueError(
            f"Unknown pipeline stop stage {stop_at!r}; "
            f"valid keys: {[s.key for s in build_pipeline()]}"
        )
    return state


def iter_pipeline_stages(
    state: Dict[str, Any],
    sample: str = "mc",
) -> Iterator[Tuple[str, Dict[str, Any]]]:
    """Yield ``(stage_key, state)`` after each stage cut (same as syst walker)."""
    cur = dict(state)
    for stage in build_pipeline():
        if stage.cut is not None:
            cur = stage.cut(cur, sample=sample)
        yield stage.key, cur


# ===========================================================================
# Sample-aware state ops
# ---------------------------------------------------------------------------
# A "state" is a flat dict carrying both per-event and per-track dataframes
# for the sample currently being processed:
#   { "evt": df, "trk": df, "hdr": df }
# We keep the SAME state structure regardless of sample so cuts are uniform.
# Cuts are written as small functions ``(state, sample) -> state`` that mutate
# (and return) the state dict.
# ===========================================================================
def _apply_to_evt(fn: Callable[[pd.DataFrame], pd.DataFrame]) -> Callable:
    """Return a stage-cut that applies ``fn`` to state['evt'] only."""
    def _cut(state, sample):
        if state.get("evt") is not None:
            state["evt"] = fn(state["evt"])
        return state
    return _cut


# ---- track-related helpers -------------------------------------------------
def _refresh_tracks_and_attach_ntrks(state, sample):
    """Re-extract valid tracks, match them to the current evt_df, and attach
    n_trks / trk1 / trk2 columns (mirrors the notebook between cuts).
    """
    if state.get("evt") is None or state.get("trk") is None:
        return state
    trk = get_valid_trks(state["trk"])
    trk = match_trkdf_to_slcdf(trk, state["evt"])
    state["trk"] = trk
    state["evt"] = get_trk_info(state["evt"], trk, SAVE_NTRKS)
    return state


def _attach_chi2_avgs(trk_df: pd.DataFrame) -> pd.DataFrame:
    """Add the per-plane-averaged chi2 columns expected by the PID cuts."""
    chimu = avg_chi2(trk_df, "chi2_muon")
    trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_muon", "")] = chimu
    chip = avg_chi2(trk_df, "chi2_proton")
    trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_proton", "")] = chip
    return trk_df


def _attach_mcs_range_diff(trk_df: pd.DataFrame) -> pd.DataFrame:
    diff = (trk_df.pfp.trk.rangeP.p_muon - trk_df.pfp.trk.mcsP.fwdP_muon) / trk_df.pfp.trk.rangeP.p_muon
    trk_df[("pfp", "trk", "mcs_range_diff", "", "", "")] = diff
    return trk_df


def _attach_evt_prim_trk_cols(state):
    """Compute ``prim_trk_*`` columns on evt_df from the longest valid track per slice.

    This is the consolidated equivalent of the multi-block notebook code that
    repeatedly sorted ``trk_df`` and copied 5+ different columns onto evt_df.
    Doing it once is much faster than the notebook style.
    """
    if state.get("evt") is None or state.get("trk") is None:
        return state
    trk = state["trk"].copy()
    if len(trk) == 0:
        return state
    trk[("pfp", "trk", "phi", "", "", "")] = np.degrees(np.arctan2(
        trk["pfp", "trk", "dir", "x", ""],
        trk["pfp", "trk", "dir", "y", ""],
    ))

    nlevels_trk = len(trk.index.names)
    # take longest track per slice
    prim = (
        trk.sort_values(("pfp", "trk", "len", "", "", ""), ascending=False)
           .groupby(level=list(range(nlevels_trk - 1)))
           .head(1)
           .reset_index(level=[nlevels_trk - 1], drop=True)
    )
    evt = state["evt"]
    cols = {
        "prim_trk_phi":               prim.pfp.trk.phi,
        "prim_trk_start_x":           prim.pfp.trk.start.x,
        "prim_trk_end_x":             prim.pfp.trk.end.x,
        "prim_trk_dir_y":             prim.pfp.trk.dir.y,
        "prim_trk_dir_z":             prim.pfp.trk.dir.z,
    }
    for name, series in cols.items():
        evt.loc[:, name] = np.nan
        evt.loc[:, name] = series

    # also frac diff in P (used by some cosmic-rejection cuts)
    try:
        frac_diff = (prim.pfp.trk.rangeP.p_muon - prim.pfp.trk.mcsP.fwdP_muon) / prim.pfp.trk.rangeP.p_muon
        evt.loc[:, "prim_trk_P_frac_diff"] = np.nan
        evt.loc[:, "prim_trk_P_frac_diff"] = frac_diff
    except Exception:
        # not all df flavours have these columns -- skip silently
        pass

    state["evt"] = evt
    return state


def _add_reco_cc1p0pi_tki_evt(evtdf: pd.DataFrame) -> pd.DataFrame:
    """Pipeline hook; implementation in :func:`pyanalib.variable_calculator.add_reco_cc1p0pi_tki_evtdf`."""
    return add_reco_cc1p0pi_tki_evtdf(evtdf)


def _add_mc_cc1p0pi_tki_mcnu(mc_nu_df: pd.DataFrame) -> pd.DataFrame:
    """GENIE / HDF hook; implementation in :func:`pyanalib.variable_calculator.add_mc_cc1p0pi_tki_mcnu`."""
    return add_mc_cc1p0pi_tki_mcnu(mc_nu_df)


def _add_truth_cc1p0pi_tki_evt(evtdf: pd.DataFrame) -> pd.DataFrame:
    """Truth TKI hook; implementation in :func:`pyanalib.variable_calculator.add_truth_cc1p0pi_tki_evtdf`."""
    return add_truth_cc1p0pi_tki_evtdf(evtdf)


# ===========================================================================
# Selectors used by PlotSpec to pluck the right rows for a given plot.
# Each takes (state, sample) and returns a DataFrame (or None to skip).
# ===========================================================================
def sel_evt(state, sample):
    """Whole-event dataframe (current state)."""
    return state.get("evt")


def sel_trks_concat(state, sample):
    """Concatenation of the two surviving tracks (per-track plots)."""
    evt = state.get("evt")
    if evt is None or len(evt) == 0:
        return None
    if "trk1" not in evt or "trk2" not in evt:
        # not all stages have trk1/trk2 yet
        return None
    return pd.concat([evt.trk1, evt.trk2])


def sel_trks_concat_not_mu(state, sample):
    """All tracks that are NOT the muon candidate.

    Mirrors notebook ``is_not_mu_candidate``: takes both tracks per slice
    and removes the one tagged as the muon candidate.
    """
    evt = state.get("evt")
    if evt is None or len(evt) == 0:
        return None
    if "trk1" not in evt or "trk2" not in evt:
        return None
    trks = pd.concat([evt.trk1, evt.trk2])
    # bug-safe: derive nlevels from trks itself, not from a global mc_df
    nlevels = len(trks.index.names) - 1  # number of *event* levels
    mcs_range_diff = np.abs(
        (trks.pfp.trk.rangeP.p_muon - trks.pfp.trk.mcsP.fwdP_muon) / trks.pfp.trk.rangeP.p_muon
    )
    chimu_avg = trks.pfp.trk.chi2pid.avg.chi2_muon
    chip_avg = trks.pfp.trk.chi2pid.avg.chi2_proton
    mu_cut = (
        (chimu_avg > 0) & (chimu_avg < MU_CHI2MU_TH) &
        (chip_avg > MU_CHI2P_TH) &
        (trks.pfp.trk.len > MU_LEN_TH) &
        (mcs_range_diff < QUAL_TH)
    )
    return pd.concat([trks[~mu_cut], trks[mu_cut].groupby(level=list(range(nlevels))).nth(1)])


# ===========================================================================
# Variables for the efficiency curve (cell 82-86 in the notebook)
# ===========================================================================
EFFICIENCY_VARS = with_final_selected_evt_variables(
    list(CORE_SELECTED_EVT_VARIABLE_CONFIGS) + [VariableConfig.neutrino_energy()]
)


# ===========================================================================
# Pipeline definition
# ===========================================================================
def build_pipeline() -> List[Stage]:
    """Return the list of stages run by the chunk processor.

    Edit this function to add/remove cuts and plots. Both the per-chunk
    processor and the aggregator import this same list, so changes here
    propagate end-to-end.
    """
    stages: List[Stage] = []

    # ------------------------------------------------------------------
    # Stage 0: all reconstructed slices (no cut)
    # ------------------------------------------------------------------
    stages.append(Stage(
        key="allreco",
        label="All reconstructed slices",
        cut=None,
        plots=[],
        save_for_efficiency=True,
        save_for_breakdown=True,  # stage counts (summary bar plot still skips this stage)
    ))

    # ------------------------------------------------------------------
    # Stage 1: not clear cosmic
    # ------------------------------------------------------------------
    stages.append(Stage(
        key="is_clear_cosmic",
        label="Not clear cosmic",
        cut=_apply_to_evt(cut_clear_cosmic),
        plots=[],
        save_for_efficiency=True,
        save_for_breakdown=True,
    ))

    # ------------------------------------------------------------------
    # Stage 2: vertex in Gen-1 fiducial volume (per-TPC analysis)
    # ------------------------------------------------------------------
    stages.append(Stage(
        key="vertex_in_fv",
        label="Vertex in Gen-1 fiducial volume",
        cut=_apply_to_evt(cut_vertex_in_fv),
        plots=[
            PlotSpec(
                var_config=VariableConfig.nu_score(),
                breakdown_type="topology",
                selector=sel_evt,
                plot_label_template=("Neutrino Score", "Events (POT={pot})", ""),
                save_kwargs={"ratio": True, "ax_ylim_ratio": 1.8, "vline": [[NU_SCORE_TH, 1]]},
            ),
        ],
        save_for_efficiency=True,
        save_for_breakdown=True,
    ))

    # ------------------------------------------------------------------
    # Stage 3: nu-score cut (then re-attach n_trks for plotting before 2prong)
    # ------------------------------------------------------------------
    def _nu_score_cut_then_refresh(state, sample):
        if state.get("evt") is not None:
            state["evt"] = cut_nu_score(state["evt"])
        # CAF early sel_level="nu" may omit the track table; skip refresh then.
        if state.get("trk") is not None:
            state = _refresh_tracks_and_attach_ntrks(state, sample)
        return state

    stages.append(Stage(
        key="nu_score",
        label=f"Nu-score > {NU_SCORE_TH}",
        cut=_nu_score_cut_then_refresh,
        plots=[
            PlotSpec(
                var_config=VariableConfig.n_trks(),
                breakdown_type="topology",
                selector=sel_evt,
                plot_label_template=("Number of tracks", "Events (POT={pot})", ""),
                save_kwargs={"ratio": True, "ax_ylim_ratio": 1.8},
            ),
            PlotSpec(
                var_config=VariableConfig.n_trks(),
                breakdown_type="genie",
                selector=sel_evt,
                plot_label_template=("Number of tracks", "Events (POT={pot})", ""),
                save_kwargs={"ratio": True, "ax_ylim_ratio": 1.8},
            ),
        ],
        save_for_efficiency=True,
        save_for_breakdown=True,
    ))

    # ------------------------------------------------------------------
    # Stage 4: exactly 2 PFPs (slice has two tracks)
    # ------------------------------------------------------------------
    stages.append(Stage(
        key="2prong",
        label="Has exactly 2 PFPs",
        cut=_apply_to_evt(cut_2prong),
        plots=[],
        save_for_efficiency=True,
        save_for_breakdown=True,
    ))

    # ------------------------------------------------------------------
    # Stage 5: both PFPs per-TPC contained
    # ------------------------------------------------------------------
    stages.append(Stage(
        key="2prong-contained",
        label="Both PFPs per-TPC contained",
        cut=_apply_to_evt(cut_2prong_contained),
        plots=[
            PlotSpec(
                var_config=VariableConfig.track_score(),
                breakdown_type="pdg",
                selector=sel_trks_concat,
                plot_label_template=("Track Score", "Tracks / Bin (POT={pot})", ""),
                save_kwargs={"ratio": True, "ax_ylim_ratio": 1.8, "vline": [[TRACKSCORE_TH, 1]]},
            ),
        ],
        save_for_efficiency=True,
        save_for_breakdown=True,
    ))

    # ------------------------------------------------------------------
    # Stage 6: track-score cut on both tracks
    # ------------------------------------------------------------------
    stages.append(Stage(
        key="2prong-trackscore",
        label=f"Both PFPs track score > {TRACKSCORE_TH}",
        cut=_apply_to_evt(cut_2prong_trackscore),
        plots=[
            PlotSpec(
                var_config=VariableConfig.vtx_dist(),
                breakdown_type="pdg",
                selector=sel_trks_concat,
                plot_label_template=("Vertex Distance (cm)", "Tracks / Bin (POT={pot})", ""),
                save_kwargs={"ratio": True, "vline": [[VTXDIST_TH, 0]]},
            ),
        ],
        save_for_efficiency=True,
        save_for_breakdown=True,
    ))

    # ------------------------------------------------------------------
    # Stage 7: vtxdist cut on both tracks (then attach chi2 averages + mcs diff)
    # ------------------------------------------------------------------
    def _vtxdist_cut_and_attach_pid_cols(state, sample):
        if state.get("evt") is not None:
            state["evt"] = cut_2prong_vtxdist(state["evt"])
        # Re-match tracks to surviving slices, attach chi2/MCS on trk, merge once onto evt.
        # Do not call _refresh_tracks_and_attach_ntrks first — a second get_trk_info on evt
        # that already has trk1/trk2 suffixes columns to trk1_x and breaks .trk1 access.
        if state.get("evt") is None or state.get("trk") is None:
            return state
        trk = get_valid_trks(state["trk"])
        trk = match_trkdf_to_slcdf(trk, state["evt"])
        state["trk"] = trk
        if len(trk) > 0:
            state["trk"] = _attach_chi2_avgs(state["trk"])
            state["trk"] = _attach_mcs_range_diff(state["trk"])
            state["evt"] = get_trk_info(state["evt"], state["trk"], SAVE_NTRKS)
        return state

    stages.append(Stage(
        key="2prong-vtxdist",
        label=f"Both tracks have (start - vertex) < {VTXDIST_TH} cm",
        cut=_vtxdist_cut_and_attach_pid_cols,
        plots=[
            PlotSpec(
                var_config=VariableConfig.trk_len(),
                breakdown_type="pdg",
                selector=sel_trks_concat,
                plot_label_template=(VariableConfig.trk_len().var_labels[0], "Tracks / Bin (POT={pot})", ""),
                save_kwargs={"ratio": True, "vline": [[50, 1]]},
            ),
            PlotSpec(
                var_config=VariableConfig.mcs_range_diff(),
                breakdown_type="pdg",
                selector=sel_trks_concat,
                plot_label_template=(VariableConfig.mcs_range_diff().var_labels[0], "Tracks / Bin (POT={pot})", ""),
                save_kwargs={"ratio": True, "vline": [[-QUAL_TH, 0], [QUAL_TH, 1]]},
            ),
            PlotSpec(
                var_config=VariableConfig.chi2_mu(),
                breakdown_type="pdg",
                selector=sel_trks_concat,
                plot_label_template=(VariableConfig.chi2_mu().var_labels[0], "Tracks / Bin (POT={pot})", ""),
                save_kwargs={"ratio": True, "vline": [[MU_CHI2MU_TH, 0]]},
            ),
            PlotSpec(
                var_config=VariableConfig.chi2_proton(),
                breakdown_type="pdg",
                selector=sel_trks_concat,
                plot_label_template=(VariableConfig.chi2_proton().var_labels[0], "Tracks / Bin (POT={pot})", ""),
                save_kwargs={"ratio": True, "vline": [[MU_CHI2P_TH, 1]], "ax_ylim_ratio": 1.8},
            ),
            # tracks that are NOT muon candidates
            PlotSpec(
                var_config=VariableConfig.chi2_mu(),
                breakdown_type="pdg",
                selector=sel_trks_concat_not_mu,
                name_suffix="not_mu",
                plot_label_template=(VariableConfig.chi2_mu().var_labels[0], "Events (POT={pot})", ""),
                save_kwargs={"ratio": True},
            ),
            PlotSpec(
                var_config=VariableConfig.chi2_proton(),
                breakdown_type="pdg",
                selector=sel_trks_concat_not_mu,
                name_suffix="not_mu",
                plot_label_template=(VariableConfig.chi2_proton().var_labels[0], "Events (POT={pot})", ""),
                save_kwargs={"ratio": True, "vline": [[MU_CHI2P_TH, 1]]},
            ),
        ],
        save_for_efficiency=True,
        save_for_breakdown=True,
    ))

    # ------------------------------------------------------------------
    # Stage 8: get mu/p candidates and apply muX cut + mu kinematics
    # ------------------------------------------------------------------
    def _muX_cut(state, sample):
        evt = state.get("evt")
        if evt is None or len(evt) == 0:
            return state
        if not evt_has_trk1_trk2(evt) and state.get("trk") is not None:
            trk = get_valid_trks(state["trk"])
            trk = match_trkdf_to_slcdf(trk, evt)
            if len(trk) > 0:
                trk = _attach_chi2_avgs(trk)
                trk = _attach_mcs_range_diff(trk)
                state["evt"] = get_trk_info(evt, trk, SAVE_NTRKS)
                evt = state["evt"]
        if not evt_has_trk1_trk2(evt):
            raise KeyError(
                "evt is missing trk1/trk2 at 2prong-muX — "
                "get_trk_info did not attach track blocks (check trk–evt matching on this shard)"
            )
        # CAF calorimetry variations pass score_tag="_new" via apply_selection_pipeline.
        pid_kw = dict(state.get("_mu_p_candidate_kwargs") or {})
        df = get_mu_p_candidate(evt, **pid_kw)
        df = cut_has_mu(df)
        df = cut_mu_kinematics(df)
        state["evt"] = df
        return state

    stages.append(Stage(
        key="2prong-muX",
        label="One track is muon-like",
        cut=_muX_cut,
        plots=[],
        save_for_efficiency=True,
        save_for_breakdown=True,
    ))

    # ------------------------------------------------------------------
    # Stage 9: mup cut + p kinematics  (= final selection)
    # ------------------------------------------------------------------
    def _mup_cut(state, sample):
        evt = state.get("evt")
        if evt is None or len(evt) == 0:
            return state
        df = cut_has_p(evt)
        df = cut_p_kinematics(df)
        df = _add_reco_cc1p0pi_tki_evt(df)
        state["evt"] = df
        return state

    # final-stage summary plots (cell 80 in the notebook)
    final_summary_plots: List[PlotSpec] = []
    _final_stage_evt_vcs = list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)
    for vc in with_final_selected_evt_variables(_final_stage_evt_vcs):
        final_summary_plots.append(PlotSpec(
            var_config=vc,
            breakdown_type="topology",
            selector=sel_evt,
            plot_label_template=(vc.var_labels[1], "Events (POT={pot})", ""),
            save_kwargs={"ratio": True, "ax_ylim_ratio": 1.6},
            name_suffix="final",
        ))

    stages.append(Stage(
        key="2prong-mup",
        label="The other is proton-like",
        cut=_mup_cut,
        plots=final_summary_plots,
        save_for_efficiency=True,
        save_for_breakdown=True,
    ))

    return stages


def build_runner(
    sample: str,
    mc_univ_syst_tags: tuple[str, ...] | None = None,
) -> ChunkRunner:
    """Build a ChunkRunner instance for a given sample.

    Always uses the same pipeline definition, but knows which sample-slot to
    fill in the histogram accumulators.

    ``mc_univ_syst_tags``: optional tuple of MC multi-universe syst names (e.g.
    ``("Flux", "G4", "GENIE")``) whose columns ``mc[s]['univ_i']`` are summed into
    chunked histograms for later fractional covariance (see aggregate flag).
    """
    return ChunkRunner(
        sample=sample,
        stages=build_pipeline(),
        efficiency_vars=EFFICIENCY_VARS,
        mc_univ_syst_tags=mc_univ_syst_tags,
    )
