"""Walk the numuCC 1p0pi event-selection pipeline on a single sample.

Used by ``syst_cosmics_chunk`` / ``syst_multisim_chunk`` when the input is a
**sel_all**-style ``.df`` (raw ``evt`` / ``trk`` / ``hdr``) so they can record
systematic histograms

* at every cut stage (the variable that drives the *next* cut), and
* at the final stage (all final-selected event-level variables).

Cuts are imported from :func:`event_selection_pipeline_def.build_pipeline`,
which is the same definition used by ``event_selection_chunk.py``, so the
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

from analysis_village.numucc_1p0pi.event_selection_pipeline_def import build_pipeline
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.selection_framework import multicol_get_series
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig


# ---------------------------------------------------------------------------
# Stage → variable mapping (mirrors the ``plots`` lists in build_pipeline()).
# ``target`` is either "evt" (column on the event df) or "trk" (column on the
# concat of evt.trk1 / evt.trk2).
# ---------------------------------------------------------------------------
class CutStageVarSpec(NamedTuple):
    stage_key: str
    var_config: Any
    target: str


def _build_cut_stage_specs() -> List[CutStageVarSpec]:
    """Cut-stage histograms mirror the PlotSpec list in build_pipeline()."""
    return [
        CutStageVarSpec("vertex_in_fv", VariableConfig.nu_score(), "evt"),
        CutStageVarSpec("nu_score", VariableConfig.n_trks(), "evt"),
        CutStageVarSpec("2prong-contained", VariableConfig.track_score(), "trk"),
        CutStageVarSpec("2prong-trackscore", VariableConfig.vtx_dist(), "trk"),
        CutStageVarSpec("2prong-vtxdist", VariableConfig.trk_len(), "trk"),
        CutStageVarSpec("2prong-vtxdist", VariableConfig.mcs_range_diff(), "trk"),
        CutStageVarSpec("2prong-vtxdist", VariableConfig.chi2_mu(), "trk"),
        CutStageVarSpec("2prong-vtxdist", VariableConfig.chi2_proton(), "trk"),
    ]


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
    ``None``). ``sample`` is one of ``"mc" / "data" / "intime" / "offbeam" / "dirt"``
    -- only relevant for sample-aware cuts (none of the cuts in the current
    pipeline branch on it, but we forward it for forward-compat with the
    pipeline definition).
    """
    cur = dict(state)
    for stage in build_pipeline():
        if trace is not None:
            trace(f"[walker] stage={stage.key!r} (cut={'yes' if stage.cut else 'no'})")
        if stage.cut is not None:
            cur = stage.cut(cur, sample=sample)
        yield stage.key, cur


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
