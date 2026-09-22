"""Walk the numuCC 1p0pi event-selection pipeline on a single sample.

Used by ``syst_cosmics_chunk`` / ``syst_multisim_chunk`` when the input is a
**sel_all**-style ``.df`` (raw ``evt`` / ``trk`` / ``hdr``) so they can record
systematic histograms

* at every cut stage (variables from ``build_pipeline()`` PlotSpecs), and
* at the final stage (all final-selected event-level variables).

Cuts are imported from :func:`event_selection_pipeline_def.build_pipeline`,
which is the same definition used by CAF / batched selection, so the
selection here is bit-for-bit identical.

For *final*-style inputs the chunk scripts skip this module and just
histogram the columns already present in the saved ``evt`` table.
"""
from __future__ import annotations

import os
import sys
from os import path
from typing import Any, Callable, Dict, Iterator, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# Make the repo importable when this module is loaded by a script using a
# relative sys.path entry.
sys.path.append(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))

from analysis_village.numucc_1p0pi.event_selection_pipeline_def import (
    build_pipeline,
    sel_evt,
    sel_trks_concat,
)
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.selection_framework import multicol_get_series
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.makedf.selections import (
    MU_CHI2MU_TH,
    MU_CHI2P_TH,
    MU_LEN_TH,
    QUAL_TH,
)


# ---------------------------------------------------------------------------
# Stage → variable mapping derived from ``build_pipeline()`` PlotSpecs.
# ``target`` is either "evt" (column on the event df) or "trk" (column on the
# concat of evt.trk1 / evt.trk2), or track subsets ``trk_len50`` / ``trk_not_mu``.
# ---------------------------------------------------------------------------
class CutStageVarSpec(NamedTuple):
    stage_key: str
    var_config: Any
    target: str


# Chi2 plane×stage grid is the only intentional extra beyond pipeline PlotSpecs
# (syst needs I0/I1/I2/avg at several stages; pipeline plots only show avg once).
CHI2_CUT_STAGES: Tuple[str, ...] = ("2prong-vtxdist", "2prong-muX", "2prong-mup")
CHI2_PLANES: Tuple[str, ...] = ("I0", "I1", "I2", "avg")

# Skip PlotSpecs that are not cut-driving diagnostics for syst rate histos.
# ``not_mu`` PlotSpecs are skipped here; dedicated subset specs below fill them.
_SKIP_NAME_SUFFIXES = frozenset({"final", "not_mu"})

# Stage for the dedicated track-subset χ² packs (avg plane only).
CHI2_SUBSET_STAGE_KEY: str = "2prong-vtxdist"

# Slugs for len>50 / not_mu avg-χ² at vtxdist (added to Product A without replacing
# existing all-track ``chi2_avg_*__at_2prong-vtxdist`` packs).
CHI2_TRACK_SUBSET_SLUGS: Tuple[str, ...] = (
    "chi2_avg_mu_len50__at_2prong-vtxdist",
    "chi2_avg_p_len50__at_2prong-vtxdist",
    "chi2_avg_mu_not_mu__at_2prong-vtxdist",
    "chi2_avg_p_not_mu__at_2prong-vtxdist",
)


def clone_var_config(vc: VariableConfig, *, var_save_name: str) -> VariableConfig:
    """Copy a ``VariableConfig`` with a distinct ``var_save_name`` (stage / plane tags)."""
    return VariableConfig(
        var_save_name=var_save_name,
        var_plot_name=vc.var_plot_name,
        var_labels=list(vc.var_labels),
        bins=np.asarray(vc.bins, dtype=float).copy(),
        var_evt_reco_col=vc.var_evt_reco_col,
        var_evt_truth_col=vc.var_evt_truth_col,
        var_nu_col=vc.var_nu_col,
        xsec_label=vc.xsec_label,
        category_syst_var_save_name=getattr(vc, "category_syst_var_save_name", None),
    )


def _stage_tagged_vc(vc: VariableConfig, stage_key: str) -> VariableConfig:
    """Unique slug so the same observable can be saved at multiple selection steps."""
    return clone_var_config(vc, var_save_name=f"{vc.var_save_name}__at_{stage_key}")


def _target_for_selector(selector: Callable) -> Optional[str]:
    if selector is sel_evt:
        return "evt"
    if selector is sel_trks_concat:
        return "trk"
    return None


def _build_cut_stage_specs() -> List[CutStageVarSpec]:
    """Cut-stage histograms from ``build_pipeline()`` plots + chi2 plane×stage grid.

    Primary specs are taken from each stage's ``PlotSpec`` list (evt / trk1+trk2
    selectors only). Chi2 is stored once per ``(plane, species, stage)`` with
    save names like ``chi2_mu_I0__at_2prong-vtxdist`` so covariances from
    different selection steps never collide.
    """
    specs: List[CutStageVarSpec] = []
    seen: set[Tuple[str, str]] = set()

    for stage in build_pipeline():
        for plot in stage.plots:
            if (plot.name_suffix or "") in _SKIP_NAME_SUFFIXES:
                continue
            target = _target_for_selector(plot.selector)
            if target is None:
                continue
            vc = plot.var_config
            slug = getattr(vc, "var_save_name", None)
            if not slug:
                continue
            # Avg chi2 PlotSpecs are superseded by the plane×stage grid below.
            if "chi2" in str(slug):
                continue
            key = (stage.key, str(slug))
            if key in seen:
                continue
            seen.add(key)
            specs.append(CutStageVarSpec(stage.key, vc, target))

    for stage_key in CHI2_CUT_STAGES:
        for plane in CHI2_PLANES:
            specs.append(
                CutStageVarSpec(
                    stage_key,
                    _stage_tagged_vc(VariableConfig.chi2_plane_mu(plane), stage_key),
                    "trk",
                )
            )
            specs.append(
                CutStageVarSpec(
                    stage_key,
                    _stage_tagged_vc(VariableConfig.chi2_plane_proton(plane), stage_key),
                    "trk",
                )
            )

    # Track-subset avg-χ² at vtxdist (len>50 and not_mu). Distinct slugs so they
    # never collide with all-track ``chi2_avg_*__at_2prong-vtxdist``.
    stage = CHI2_SUBSET_STAGE_KEY
    specs.extend(
        [
            CutStageVarSpec(
                stage,
                clone_var_config(
                    VariableConfig.chi2_avg_mu(),
                    var_save_name="chi2_avg_mu_len50__at_2prong-vtxdist",
                ),
                "trk_len50",
            ),
            CutStageVarSpec(
                stage,
                clone_var_config(
                    VariableConfig.chi2_avg_proton(),
                    var_save_name="chi2_avg_p_len50__at_2prong-vtxdist",
                ),
                "trk_len50",
            ),
            CutStageVarSpec(
                stage,
                clone_var_config(
                    VariableConfig.chi2_avg_mu(),
                    var_save_name="chi2_avg_mu_not_mu__at_2prong-vtxdist",
                ),
                "trk_not_mu",
            ),
            CutStageVarSpec(
                stage,
                clone_var_config(
                    VariableConfig.chi2_avg_proton(),
                    var_save_name="chi2_avg_p_not_mu__at_2prong-vtxdist",
                ),
                "trk_not_mu",
            ),
        ]
    )
    return specs


CUT_STAGE_VAR_SPECS: Tuple[CutStageVarSpec, ...] = tuple(_build_cut_stage_specs())

# ``var_save_name`` values for cut-stage observables (GENIE / multisim sel_all: rate-only).
CUT_STAGE_RATE_ONLY_SLUGS: frozenset = frozenset(
    spec.var_config.var_save_name for spec in CUT_STAGE_VAR_SPECS
)


def active_cut_stage_specs() -> Tuple[CutStageVarSpec, ...]:
    """``CUT_STAGE_VAR_SPECS``, optionally filtered by ``NUMUCC_CUT_STAGE_SLUGS``.

    Set ``NUMUCC_CUT_STAGE_SLUGS=slug1,slug2,...`` to histogram only those cut-stage
    variables (used for focused chi2 subset campaigns without redoing all packs).
    """
    raw = os.environ.get("NUMUCC_CUT_STAGE_SLUGS", "").strip()
    if not raw:
        return CUT_STAGE_VAR_SPECS
    allow = {p.strip() for p in raw.split(",") if p.strip()}
    return tuple(s for s in CUT_STAGE_VAR_SPECS if s.var_config.var_save_name in allow)

# Stage at which final-selected variables are histogrammed (after the full chain).
FINAL_STAGE_KEY: str = "2prong-mup"


def final_stage_var_configs(extra: Optional[Sequence[VariableConfig]] = None) -> List[VariableConfig]:
    """Final-selection variables (CORE + with_final_selected_evt_variables(...))."""
    base: List[VariableConfig] = list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)
    if extra:
        for vc in extra:
            base.append(vc)
    return with_final_selected_evt_variables(base)


# ---------------------------------------------------------------------------
# Pipeline walking
# ---------------------------------------------------------------------------
def walk_pipeline(
    state: Dict[str, Any],
    sample: str,
    trace: Optional[Callable[[str], None]] = None,
) -> Iterator[Tuple[str, Dict[str, Any]]]:
    """Yield ``(stage_key, state)`` for each pipeline stage after its cut is applied.

    ``state`` must contain at least ``evt`` / ``trk`` / ``hdr`` (``mcnu`` may be
    ``None``). Delegates to ``build_pipeline()`` so cuts match CAF / batched selection.
    """
    from analysis_village.numucc_1p0pi.event_selection_pipeline_def import (
        iter_pipeline_stages,
    )

    for stage_key, cur in iter_pipeline_stages(state, sample=sample):
        if trace is not None:
            trace(f"[walker] stage={stage_key!r}")
        yield stage_key, cur


# ---------------------------------------------------------------------------
# Column extraction (evt or per-trk concat) with MultiIndex-depth handling.
# ---------------------------------------------------------------------------
def _evt_has_trk1_trk2(evt: pd.DataFrame) -> bool:
    if evt is None or len(evt) == 0:
        return False
    try:
        top = evt.columns.get_level_values(0).unique()
    except Exception:
        return False
    return ("trk1" in top) and ("trk2" in top)


def get_var_series(
    state: Dict[str, Any],
    var_config: VariableConfig,
    target: str,
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Return ``(values, evt_index_for_each_value)`` or ``None`` if not available.

    ``evt_index_for_each_value`` is an integer position array mapping each
    output row to its source event index (so per-universe weights can be
    broadcast). For ``target == "trk"`` each event contributes 2 rows
    (trk1 + trk2) so the returned indices repeat.

    Track subsets:
    * ``trk_len50`` — concat trk1+trk2 with ``pfp.trk.len > MU_LEN_TH`` (50 cm)
    * ``trk_not_mu`` — same population as :func:`sel_trks_concat_not_mu`
    """
    evt = state.get("evt")
    if evt is None or len(evt) == 0:
        return None
    if target == "evt":
        try:
            s = multicol_get_series(evt, var_config.var_evt_reco_col)
        except KeyError:
            return None
        v = np.asarray(s, dtype=float)
        idx = np.arange(len(v), dtype=np.int64)
        return v, idx
    if target == "trk":
        return _trk_concat_series(evt, var_config)
    if target == "trk_len50":
        return _trk_len50_series(evt, var_config)
    if target == "trk_not_mu":
        return _trk_not_mu_series(state, var_config)
    raise ValueError(
        f"target must be 'evt', 'trk', 'trk_len50', or 'trk_not_mu', got {target!r}"
    )


def _trk_concat_series(
    evt: pd.DataFrame, var_config: VariableConfig
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    if not _evt_has_trk1_trk2(evt):
        return None
    try:
        s1 = multicol_get_series(evt.trk1, var_config.var_evt_reco_col)
        s2 = multicol_get_series(evt.trk2, var_config.var_evt_reco_col)
    except KeyError:
        return None
    v1 = np.asarray(s1, dtype=float)
    v2 = np.asarray(s2, dtype=float)
    n = len(evt)
    idx = np.arange(n, dtype=np.int64)
    v = np.concatenate([v1, v2])
    idx_full = np.concatenate([idx, idx])
    return v, idx_full


def _trk_len_series(evt: pd.DataFrame) -> Optional[np.ndarray]:
    """Concatenated ``pfp.trk.len`` for trk1+trk2, or None."""
    if not _evt_has_trk1_trk2(evt):
        return None
    try:
        l1 = multicol_get_series(evt.trk1, ("pfp", "trk", "len"))
        l2 = multicol_get_series(evt.trk2, ("pfp", "trk", "len"))
    except KeyError:
        try:
            l1 = evt.trk1.pfp.trk.len
            l2 = evt.trk2.pfp.trk.len
        except Exception:
            return None
    return np.concatenate([np.asarray(l1, dtype=float), np.asarray(l2, dtype=float)])


def _trk_len50_series(
    evt: pd.DataFrame, var_config: VariableConfig
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    got = _trk_concat_series(evt, var_config)
    if got is None:
        return None
    v, idx = got
    lengths = _trk_len_series(evt)
    if lengths is None or lengths.shape[0] != v.shape[0]:
        return None
    mask = np.asarray(lengths, dtype=float) > float(MU_LEN_TH)
    if not mask.any():
        return (
            np.zeros(0, dtype=np.float64),
            np.zeros(0, dtype=np.int64),
        )
    return v[mask], idx[mask]


def _trk_not_mu_series(
    state: Dict[str, Any], var_config: VariableConfig
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Avg/plane χ² on the not-muon-candidate track population (same as plot selector)."""
    evt = state.get("evt")
    if evt is None or len(evt) == 0 or not _evt_has_trk1_trk2(evt):
        return None
    got = _trk_concat_series(evt, var_config)
    if got is None:
        return None
    v, idx = got
    n = len(evt)
    if v.shape[0] != 2 * n:
        return None

    try:
        trks = pd.concat([evt.trk1, evt.trk2])
        mcs_range_diff = np.abs(
            (trks.pfp.trk.rangeP.p_muon - trks.pfp.trk.mcsP.fwdP_muon)
            / trks.pfp.trk.rangeP.p_muon
        )
        chimu_avg = trks.pfp.trk.chi2pid.avg.chi2_muon
        chip_avg = trks.pfp.trk.chi2pid.avg.chi2_proton
        lengths = trks.pfp.trk.len
    except Exception:
        return None

    pid_kw = dict(state.get("_mu_p_candidate_kwargs") or {})
    mu_chi2mu_th = float(pid_kw.get("mu_chi2mu_th", MU_CHI2MU_TH))
    mu_chi2p_th = float(pid_kw.get("mu_chi2p_th", MU_CHI2P_TH))
    mu_len_th = float(pid_kw.get("mu_len_th", MU_LEN_TH))
    qual_th = float(pid_kw.get("qual_th", QUAL_TH))
    mu_cut = np.asarray(
        (chimu_avg > 0)
        & (chimu_avg < mu_chi2mu_th)
        & (chip_avg > mu_chi2p_th)
        & (lengths > mu_len_th)
        & (mcs_range_diff < qual_th),
        dtype=bool,
    )
    if mu_cut.shape[0] != 2 * n:
        return None
    mu1 = mu_cut[:n]
    mu2 = mu_cut[n:]
    # Match sel_trks_concat_not_mu: keep non-mu tracks; if both tracks are mu,
    # keep the second (trk2) only.
    keep1 = ~mu1
    keep2 = (~mu2) | (mu1 & mu2)
    keep = np.concatenate([keep1, keep2])
    if not keep.any():
        return (
            np.zeros(0, dtype=np.float64),
            np.zeros(0, dtype=np.int64),
        )
    return v[keep], idx[keep]


# ---------------------------------------------------------------------------
# Helpers used by both cosmics and multisim chunk drivers
# ---------------------------------------------------------------------------
def histogram_var(
    values: np.ndarray,
    bins: np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    """np.histogram wrapper that clips values, sanitises non-finite weights."""
    if values is None or len(values) == 0:
        return np.zeros(len(bins) - 1, dtype=np.float64)
    v = np.asarray(values, dtype=float)
    if weights is None:
        w = None
    else:
        w = np.asarray(weights, dtype=float)
        w = np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)
    finite = np.isfinite(v)
    if not finite.all():
        v = v[finite]
        if w is not None:
            w = w[finite]
    eps = (float(bins[-1]) - float(bins[0])) * 1e-9
    v = np.clip(v, bins[0], bins[-1] - eps)
    h, _ = np.histogram(v, bins=bins, weights=w)
    return h.astype(np.float64, copy=False)
