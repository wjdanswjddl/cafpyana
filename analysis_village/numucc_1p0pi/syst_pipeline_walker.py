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


# ---------------------------------------------------------------------------
# Stage → variable mapping derived from ``build_pipeline()`` PlotSpecs.
# ``target`` is either "evt" (column on the event df) or "trk" (column on the
# concat of evt.trk1 / evt.trk2).
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
_SKIP_NAME_SUFFIXES = frozenset({"final", "not_mu"})


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
    return specs


CUT_STAGE_VAR_SPECS: Tuple[CutStageVarSpec, ...] = tuple(_build_cut_stage_specs())

# ``var_save_name`` values for cut-stage observables (GENIE / multisim sel_all: rate-only).
CUT_STAGE_RATE_ONLY_SLUGS: frozenset = frozenset(
    spec.var_config.var_save_name for spec in CUT_STAGE_VAR_SPECS
)

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
    raise ValueError(f"target must be 'evt' or 'trk', got {target!r}")


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
