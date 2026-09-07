"""Memory-efficient event-selection framework for numuCC 1p0pi.

This module turns the all-in-memory event_selection-CV.ipynb notebook into a
chunked pipeline that

  1. loads ONE input file at a time (per sample),
  2. runs the full selection pipeline on that single file,
  3. at every "plot point" in the pipeline, bins the relevant variables into
     per-category histograms (sum-of-weights AND sum-of-weights^2) instead of
     producing a plot,
  4. saves all the histogram contents (one pickle per (sample, chunk)).

A second pass (see scripts/event_selection_aggregate.py) reads back the per-chunk
pickles, sums intrinsic-weight histograms across chunks, applies **global** POT /
gates factors from summed chunk metadata, and then renders
the final plots through ``overlay_hists_from_histdata`` -- which produces a plot
that is bit-for-bit identical to ``overlay_hists`` called with the raw
dataframes.

Memory note
-----------
After this refactor the scripts only ever hold one file's events in RAM at a
time, plus the (negligible) accumulated histograms. This is what makes it run
on the full sample.

Systematics on overlay plots
----------------------------
Consumer paths load **pre-saved** fractional covariances from the syst-disk tree
via ``utils.get_syst_unc`` (``NUMUCC_SYST_DISK_ROOT`` / ``--syst-disk-root``), or
the pre-summed ``CategorySummary/category_syst_summary.npz`` when overlay helpers
default to that. On-the-fly covariances from MC multi-universe weights at plot
time are retired (see ``cafpyana_trash/…/legacy_get_frac_unc.py`` and
``legacy_frac_cov_from_mc_univ_histdata.py``). Producer pipelines that *write*
those pre-saved files still use ``utils.get_univ_rates`` / DETVAR scripts.

Adding new things
-----------------
* A new selection cut    -> add a new ``Stage`` to ``build_pipeline()``.
* A new variable to plot -> add a new ``PlotSpec`` to a stage's ``plots`` list.
* A new variable to follow through every stage -> add a ``VariableConfig`` to
  the ``efficiency_vars`` list passed to the pipeline.

Compat
------
The module supports either MultiIndex per-event dataframes (``mc_df``) or
per-track dataframes built via ``pd.concat([df.trk1, df.trk2])`` -- a
``PlotSpec`` simply specifies which one to use through its ``selector`` field.
"""

from __future__ import annotations

import sys
import os
from os import path, makedirs
from dataclasses import dataclass, field
from typing import Callable, List, Dict, Optional, Any, Tuple
import pickle
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

# Local imports (kept minimal here; the pipeline definition imports the rest)
sys.path.append(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))

from pyanalib.pandas_helpers import pad_column_name

from analysis_village.numucc_1p0pi.categories import (
    get_topo_category, get_genie_category, get_genie_sb_category, get_pdg_category,
    topology_labels, genie_mode_labels, genie_sb_mode_labels,
    pdg_labels,
    IsNuInFV_NumuCC_1p0pi, DETECTOR,
)
# Keep in sync with ``analysis_village.numucc_1p0pi.constants.EPSILON`` — duplicated here so we do not
# import that module (it pulls uproot/matplotlib at import time).
_EPS_CLIP = 1e-6

# ---------------------------------------------------------------------------
# MC multi-universe weights (Flux / G4 / GENIE …) — map-phase optional bookkeeping.
# Consumer-side systematic bands load pre-saved files via utils.get_syst_unc;
# on-the-fly cov from these histograms was retired (see cafpyana_trash
# legacy_frac_cov_from_mc_univ_histdata.py).
# ---------------------------------------------------------------------------
def mc_univ_weight_matrix(df: pd.DataFrame, syst_tag: str, max_univ: int = 512) -> Optional[np.ndarray]:
    """Return ``(n_univ, n_evt)`` multiplicative weights for ``mc[syst_tag]['univ_i']`` columns.

    Rows are consecutive ``univ_0 … univ_{n-1}`` until the first missing column.
    Missing / non-finite weights are sanitised (NaN and infinities treated as 1).
    """
    if df is None or len(df) == 0:
        return None
    rows = []
    for i in range(max_univ):
        key = multicol_resolve_column_key(df, ("mc", syst_tag, f"univ_{i}"))
        if key is None:
            break
        w = np.asarray(df.loc[:, key], dtype=float).reshape(-1)
        w = np.nan_to_num(w, nan=1.0, posinf=10.0, neginf=0.0)
        w = np.clip(w, 0.0, 10.0)
        rows.append(w)
    if not rows:
        return None
    return np.stack(rows, axis=0)


def _squeeze_cfg_key_to_depth(parts: List[Any], n: int) -> Tuple[Any, ...]:
    """Match VariableConfig tuples to on-disk MultiIndex depth ``n``.

    HDF/event frames often omit the intermediate ``p`` particle level used in
    accessor paths (e.g. ``… truth.p.totp`` is stored as ``… truth.totp``).
    Strip trailing ``''``, collapse ``(…, 'p', leaf)`` and ``(…, 'p', mid, leaf)``
    when still too long, then truncate if needed.
    """
    q = list(parts)
    for _ in range(max(len(q), n)):
        if len(q) <= n:
            break
        if len(q) >= 2 and q[-2] == "p":
            q = q[:-2] + [q[-1]]
            continue
        if len(q) >= 3 and q[-3] == "p":
            q = q[:-3] + [q[-1]]
            continue
        break
    if len(q) > n:
        q = q[:n]
    return tuple(q)


def multicol_resolve_column_key(df: pd.DataFrame, var_col):
    """Return a column key present on ``df``, or ``None`` if no matching column.

    Handles MultiIndex depth / ``p``-level mismatches the same way as
    ``multicol_get_series``.
    """
    if df is None or len(df.columns) == 0:
        return None
    if not isinstance(df.columns, pd.MultiIndex) or not isinstance(var_col, tuple):
        return var_col if var_col in df.columns else None
    n = df.columns.nlevels
    parts = list(var_col)
    while parts and parts[-1] == "":
        parts.pop()
    if not parts:
        return None
    squeezed = _squeeze_cfg_key_to_depth(parts, n)
    key = pad_column_name(squeezed, n)
    if key in df.columns:
        return key
    key2 = pad_column_name(tuple(parts[:n]), n)
    if key2 != key and key2 in df.columns:
        return key2
    return None


def multicol_get_series(df: pd.DataFrame, var_col):
    """Resolve a column label for MultiIndex frames when tuple depth does not match.

    Variable configs often use trailing ``''`` padding or deeper tuples than the
    on-disk HDF5 column index. See ``_squeeze_cfg_key_to_depth``.
    """
    key = multicol_resolve_column_key(df, var_col)
    if key is None:
        raise KeyError(var_col)
    return df.loc[:, key]


def get_clipped_evts(df, var_col, bins, verbose=False):
    """Clip variable to bin range and return weights (same logic as ``utils.get_clipped_evts``).

    Kept here so the chunk pipeline does not import ``utils.py`` (which pulls matplotlib +
    statsmodels at module import time and can destabilize headless batch jobs).

    Non-finite weights are zeroed: ``numpy.histogram(..., weights=...)`` propagates a single
    NaN weight to **every bin**, which wipes whole MC topology slices in aggregates.
    """
    var = multicol_get_series(df, var_col)
    var = np.asarray(var, dtype=float)
    var = np.clip(var, bins[0], bins[-1] - _EPS_CLIP)
    if "pot_weight" in df.columns:
        weights = df.loc[:, "pot_weight"]
    else:
        if verbose:
            print("No pot_weight column found, return 1 as pot scale (expected for data)")
        weights = np.ones_like(var)
    weights = np.asarray(weights, dtype=float)
    weights = np.nan_to_num(weights, nan=0.0, posinf=0.0, neginf=0.0)
    return var, weights


# ---------------------------------------------------------------------------
# Map breakdown_type -> (n_categories, get_cuts_fn)
#   - n_categories must match the LENGTH of the list returned by
#     ``get_*_category(df, ret_cuts=True)``.
#   - The histograms in OverlayHistData.mc_hist are stored in the SAME order as
#     the cuts returned by ``get_*_category(...ret_cuts=True)`` -- which is the
#     order overlay_hists / overlay_hists_from_histdata expect for stacking.
# ---------------------------------------------------------------------------
BREAKDOWN_REGISTRY: Dict[str, Tuple[int, Callable]] = {
    "topology": (len(topology_labels), get_topo_category),
    "genie":    (len(genie_mode_labels), get_genie_category),
    "genie_sb": (len(genie_sb_mode_labels), get_genie_sb_category),
    "pdg":      (len(pdg_labels), get_pdg_category),
}


# ===========================================================================
# Pre-binned histogram container -- the unit of data exchanged between the
# per-chunk processor and the aggregator.
# ===========================================================================
@dataclass
class OverlayHistData:
    """Per-(stage, variable, breakdown_type) accumulator.

    Histograms in ``mc_hist`` and ``mc_err2`` are indexed [category, bin] in
    the SAME order as the cuts returned by ``get_*_category(ret_cuts=True)``
    so the aggregator can pass them straight to ``overlay_hists_from_histdata``.

    All histograms hold POT-weighted sums (i.e. each event contributes its
    ``pot_weight`` to ``*_hist`` and its ``pot_weight**2`` to ``*_err2``).
    Histograms are additive across chunks AND across samples handled at
    different times.
    """
    var_save_name: str
    breakdown_type: str
    bins: np.ndarray

    mc_hist: np.ndarray = field(default=None)        # (n_cat, n_bin)
    mc_err2: np.ndarray = field(default=None)        # (n_cat, n_bin)
    intime_hist: np.ndarray = field(default=None)    # (n_bin,)
    intime_err2: np.ndarray = field(default=None)
    offbeam_hist: np.ndarray = field(default=None) # (n_bin,) scaled cosmics from offbeam (optional alt.)
    offbeam_err2: np.ndarray = field(default=None)
    dirt_hist: np.ndarray = field(default=None)
    dirt_err2: np.ndarray = field(default=None)
    data_hist: np.ndarray = field(default=None)
    data_err2: np.ndarray = field(default=None)
    # MC-only: optional per-syst universe histograms for chunked syst covariance:
    #   mc_univ_hist[syst_tag] -> (n_univ, n_cat, n_bin), same category order as mc_hist.
    mc_univ_hist: Optional[Dict[str, np.ndarray]] = None
    has_mc: bool = False
    has_intime: bool = False
    has_offbeam: bool = False
    has_dirt: bool = False
    has_data: bool = False

    def __post_init__(self):
        n_cat = BREAKDOWN_REGISTRY[self.breakdown_type][0]
        n_bin = len(self.bins) - 1
        if self.mc_hist is None:
            self.mc_hist = np.zeros((n_cat, n_bin))
            self.mc_err2 = np.zeros((n_cat, n_bin))
        if self.intime_hist is None:
            self.intime_hist = np.zeros(n_bin)
            self.intime_err2 = np.zeros(n_bin)
        if self.offbeam_hist is None:
            self.offbeam_hist = np.zeros(n_bin)
            self.offbeam_err2 = np.zeros(n_bin)
        if self.dirt_hist is None:
            self.dirt_hist = np.zeros(n_bin)
            self.dirt_err2 = np.zeros(n_bin)
        if self.data_hist is None:
            self.data_hist = np.zeros(n_bin)
            self.data_err2 = np.zeros(n_bin)

    # ----- accumulation ---------------------------------------------------
    def fill_from_df(
        self,
        df: pd.DataFrame,
        var_col,
        sample: str,
        mc_univ_syst_tags: Optional[Tuple[str, ...]] = None,
    ):
        """Fill from one chunk of one sample. ``sample`` is mc/intime/dirt/data.

        Uses ``df.pot_weight`` (already attached upstream) unless omitted — see
        ``fill_from_df_intrinsic`` used by the deferred-global-normalisation workflow.

        ``mc_univ_syst_tags``: when ``sample=='mc'``, also accumulate weighted histograms for
        each multi-universe syst (columns ``mc[s]['univ_i']``), used at aggregation time to
        build a fractional covariance across universes at aggregation time.
        """
        if df is None:
            return
        # mark sample as "used" even if empty so that the renderer knows what
        # to draw vs. what to skip.
        if sample == "mc":
            self.has_mc = True
        elif sample == "intime":
            self.has_intime = True
        elif sample == "offbeam":
            self.has_offbeam = True
        elif sample == "dirt":
            self.has_dirt = True
        elif sample == "data":
            self.has_data = True
        else:
            raise ValueError(f"unknown sample: {sample!r}")

        if len(df) == 0:
            return

        var, weights = get_clipped_evts(df, var_col, self.bins)
        weights = np.asarray(weights, dtype=float)

        if sample == "mc":
            _, get_cuts_fn = BREAKDOWN_REGISTRY[self.breakdown_type]
            cuts = get_cuts_fn(df, ret_cuts=True)
            for ic, cut in enumerate(cuts):
                cut = np.asarray(cut, dtype=bool)
                if not cut.any():
                    continue
                v = var[cut]
                w = weights[cut]
                h, _  = np.histogram(v, bins=self.bins, weights=w)
                e2, _ = np.histogram(v, bins=self.bins, weights=np.square(w))
                self.mc_hist[ic] += h
                self.mc_err2[ic] += e2

            if mc_univ_syst_tags:
                idx_evt = np.arange(len(df))
                if getattr(self, "mc_univ_hist", None) is None:
                    self.mc_univ_hist = {}
                for syst in mc_univ_syst_tags:
                    uw_full = mc_univ_weight_matrix(df, syst)
                    if uw_full is None:
                        continue
                    n_univ = uw_full.shape[0]
                    if syst not in self.mc_univ_hist:
                        self.mc_univ_hist[syst] = np.zeros(
                            (n_univ, len(cuts), len(self.bins) - 1), dtype=float
                        )
                    elif self.mc_univ_hist[syst].shape != (
                        n_univ,
                        len(cuts),
                        len(self.bins) - 1,
                    ):
                        raise ValueError(
                            f"mc_univ_hist shape mismatch for syst={syst!r} "
                            f"(got {self.mc_univ_hist[syst].shape}, expected "
                            f"{(n_univ, len(cuts), len(self.bins) - 1)})"
                        )
                    for ic, cut in enumerate(cuts):
                        cut = np.asarray(cut, dtype=bool)
                        if not cut.any():
                            continue
                        v = var[cut]
                        w = weights[cut]
                        idx = idx_evt[cut]
                        uw_sub = uw_full[:, idx]
                        for iu in range(n_univ):
                            h_u, _ = np.histogram(
                                v,
                                bins=self.bins,
                                weights=w * uw_sub[iu],
                            )
                            self.mc_univ_hist[syst][iu, ic] += h_u
        else:
            h, _  = np.histogram(var, bins=self.bins, weights=weights)
            e2, _ = np.histogram(var, bins=self.bins, weights=np.square(weights))
            target_h, target_e2 = {
                "intime": (self.intime_hist, self.intime_err2),
                "offbeam": (self.offbeam_hist, self.offbeam_err2),
                "dirt":   (self.dirt_hist, self.dirt_err2),
                "data":   (self.data_hist, self.data_err2),
            }[sample]
            target_h += h
            target_e2 += e2

    def __iadd__(self, other: "OverlayHistData"):
        assert self.var_save_name == other.var_save_name
        assert self.breakdown_type == other.breakdown_type
        assert np.array_equal(self.bins, other.bins)
        self.mc_hist += other.mc_hist
        self.mc_err2 += other.mc_err2
        self.intime_hist += other.intime_hist
        self.intime_err2 += other.intime_err2
        self.offbeam_hist += other.offbeam_hist
        self.offbeam_err2 += other.offbeam_err2
        self.dirt_hist += other.dirt_hist
        self.dirt_err2 += other.dirt_err2
        self.data_hist += other.data_hist
        self.data_err2 += other.data_err2
        ou = getattr(other, "mc_univ_hist", None)
        if ou:
            if self.mc_univ_hist is None:
                self.mc_univ_hist = {
                    k: np.asarray(v, dtype=float).copy() for k, v in ou.items()
                }
            else:
                for k, v in ou.items():
                    v = np.asarray(v, dtype=float)
                    if k not in self.mc_univ_hist:
                        self.mc_univ_hist[k] = v.copy()
                    else:
                        if self.mc_univ_hist[k].shape != v.shape:
                            raise ValueError(
                                f"mc_univ_hist[{k!r}] shape mismatch in merge: "
                                f"{self.mc_univ_hist[k].shape} vs {v.shape}"
                            )
                        self.mc_univ_hist[k] += v
        self.has_mc     = self.has_mc     or other.has_mc
        self.has_intime = self.has_intime or other.has_intime
        self.has_offbeam = self.has_offbeam or other.has_offbeam
        self.has_dirt   = self.has_dirt   or other.has_dirt
        self.has_data   = self.has_data   or other.has_data
        return self


# ===========================================================================
# Auxiliary accumulators
# ===========================================================================
@dataclass
class BarBreakdown:
    """Per-stage event-count breakdown, for the summary bar plot.

    Stores POT-weighted event counts per topology/genie category for MC and
    plain POT-weighted counts for intime/dirt. Per-stage and per-breakdown.
    """
    breakdown_type: str
    mc_counts: np.ndarray   # length n_cat, in cuts order
    intime_count: float = 0.0
    offbeam_count: float = 0.0
    dirt_count: float = 0.0
    data_count: float = 0.0

    @classmethod
    def empty(cls, breakdown_type: str) -> "BarBreakdown":
        n_cat = BREAKDOWN_REGISTRY[breakdown_type][0]
        return cls(breakdown_type=breakdown_type, mc_counts=np.zeros(n_cat))

    def fill_mc(self, df: pd.DataFrame):
        if df is None or len(df) == 0:
            return
        weights = df["pot_weight"] if "pot_weight" in df.columns else np.ones(len(df))
        weights = np.asarray(weights, dtype=float)
        cuts = BREAKDOWN_REGISTRY[self.breakdown_type][1](df, ret_cuts=True)
        for ic, cut in enumerate(cuts):
            cut = np.asarray(cut, dtype=bool)
            self.mc_counts[ic] += float(np.sum(weights[cut]))

    def fill_intime(self, df: pd.DataFrame):
        if df is None or len(df) == 0:
            return
        weights = df["pot_weight"] if "pot_weight" in df.columns else np.ones(len(df))
        self.intime_count += float(np.asarray(weights, dtype=float).sum())

    def fill_offbeam(self, df: pd.DataFrame):
        if df is None or len(df) == 0:
            return
        weights = df["pot_weight"] if "pot_weight" in df.columns else np.ones(len(df))
        self.offbeam_count += float(np.asarray(weights, dtype=float).sum())

    def fill_dirt(self, df: pd.DataFrame):
        if df is None or len(df) == 0:
            return
        weights = df["pot_weight"] if "pot_weight" in df.columns else np.ones(len(df))
        self.dirt_count += float(np.asarray(weights, dtype=float).sum())

    def fill_data(self, df: pd.DataFrame):
        if df is None or len(df) == 0:
            return
        weights = df["pot_weight"] if "pot_weight" in df.columns else np.ones(len(df))
        self.data_count += float(np.asarray(weights, dtype=float).sum())

    def total_count(self, sample: str) -> float:
        """POT-weighted slice count at this stage (sum of ``pot_weight``)."""
        if sample == "mc":
            return float(np.sum(self.mc_counts))
        if sample == "intime":
            return self.intime_count
        if sample == "offbeam":
            return self.offbeam_count
        if sample == "dirt":
            return self.dirt_count
        if sample == "data":
            return self.data_count
        raise ValueError(f"unknown sample: {sample!r}")

    def __iadd__(self, other: "BarBreakdown"):
        assert self.breakdown_type == other.breakdown_type
        self.mc_counts += other.mc_counts
        self.intime_count += other.intime_count
        self.offbeam_count += other.offbeam_count
        self.dirt_count += other.dirt_count
        self.data_count += other.data_count
        return self


@dataclass
class EfficiencyAccumulator:
    """Per-(stage, variable) data needed to compute the efficiency curve.

    * **Denominator (total MC truth)** — histogram of generated signal on the
      ``mcnu`` table using ``VariableConfig.var_nu_col``. Filled only at the **first**
      efficiency stage (typically ``allreco``) so summed neutrino denominators match
      one POT pass per chunk.

    * **Numerator (selected truth)** — histogram of reconstructed signal slices on
      ``evt`` using truth carried on the slice. Prefer ``VariableConfig.var_nu_col``
      (available from the first stage for μ/p kinematics), otherwise fall back to
      ``VariableConfig.var_evt_truth_col`` (matched-track truth, available after mu/p ID).

    Raw Wilson intervals compare per-bin ``n_signal_raw`` to ``n_truth_nu_raw`` from the
    first stage's denominator block.
    """
    var_save_name: str
    bins: np.ndarray
    # Selected νμ CC 1p0π in FV: evt × var_evt_truth_col (per stage)
    n_signal_pot: np.ndarray
    n_signal_raw: np.ndarray
    # Same signal definition on mcnu × var_nu_col (first efficiency stage only)
    n_truth_nu_pot: np.ndarray
    n_truth_nu_raw: np.ndarray
    n_total_signal_int: float = 0.0  # total integral of signal POT-weighted (evt)
    n_total_signal_int_raw: float = 0.0
    n_at_stage_int: float = 0.0     # total events (any topology) POT-weighted (evt)
    n_at_stage_int_raw: float = 0.0  # raw event count at stage (notebook parity)

    @classmethod
    def empty(cls, var_config) -> "EfficiencyAccumulator":
        n_bin = len(var_config.bins) - 1
        return cls(
            var_save_name=var_config.var_save_name,
            bins=np.asarray(var_config.bins).copy(),
            n_signal_pot=np.zeros(n_bin),
            n_signal_raw=np.zeros(n_bin),
            n_truth_nu_pot=np.zeros(n_bin),
            n_truth_nu_raw=np.zeros(n_bin),
        )

    def fill_denominator_from_mcnu(self, mcnu_df, var_config):
        """Incremental truth histogram on generated neutrinos (call once/chunk/stage-loop)."""
        if mcnu_df is None or len(mcnu_df) == 0:
            return
        nu_col = getattr(var_config, "var_nu_col", None)
        if nu_col is None:
            return
        if multicol_resolve_column_key(mcnu_df, nu_col) is None:
            return
        weights = (
            mcnu_df["pot_weight"]
            if "pot_weight" in mcnu_df.columns
            else np.ones(len(mcnu_df))
        )
        weights = np.asarray(weights, dtype=float)
        sig_mask = IsNuInFV_NumuCC_1p0pi(mcnu_df, detector=DETECTOR)
        sig_df = mcnu_df[sig_mask]
        if len(sig_df) == 0:
            return
        sm = np.asarray(sig_mask.values, dtype=bool) if hasattr(sig_mask, "values") \
            else np.asarray(sig_mask, dtype=bool)
        sig_w = weights[sm]
        var_sig, _ = get_clipped_evts(sig_df, nu_col, self.bins)
        h_pot, _ = np.histogram(var_sig, bins=self.bins, weights=sig_w)
        h_raw, _ = np.histogram(var_sig, bins=self.bins)
        self.n_truth_nu_pot += h_pot
        self.n_truth_nu_raw += h_raw

    def fill_numerator_from_evt(self, evt_df, var_config, signal_mask_evt):
        """Signal spectrum on reco slices × truth-on-slice column."""
        if evt_df is None or len(evt_df) == 0:
            return
        # For "full selection" efficiency curves we want stage-by-stage spectra even
        # before mu/p candidate assignment. Those stages do not have matched track-truth
        # columns yet (var_evt_truth_col), but they DO carry generator truth in the
        # evt.mc block (var_nu_col). Use that when present; otherwise fall back.
        truth_col = getattr(var_config, "var_nu_col", None)
        if truth_col is None or multicol_resolve_column_key(evt_df, truth_col) is None:
            truth_col = var_config.var_evt_truth_col
        if multicol_resolve_column_key(evt_df, truth_col) is None:
            return
        weights = evt_df["pot_weight"] if "pot_weight" in evt_df.columns else np.ones(len(evt_df))
        weights = np.asarray(weights, dtype=float)
        sig_df = evt_df[signal_mask_evt]
        sm = np.asarray(signal_mask_evt.values, dtype=bool) if hasattr(signal_mask_evt, "values") \
            else np.asarray(signal_mask_evt, dtype=bool)
        sig_w = weights[sm]
        if len(sig_df) == 0:
            self.n_at_stage_int += float(weights.sum())
            self.n_at_stage_int_raw += float(len(evt_df))
            return
        var_sig, _ = get_clipped_evts(sig_df, truth_col, self.bins)
        h_pot, _ = np.histogram(var_sig, bins=self.bins, weights=sig_w)
        h_raw, _ = np.histogram(var_sig, bins=self.bins)
        self.n_signal_pot += h_pot
        self.n_signal_raw += h_raw
        self.n_total_signal_int += float(sig_w.sum())
        self.n_total_signal_int_raw += float(len(sig_df))
        self.n_at_stage_int += float(weights.sum())
        self.n_at_stage_int_raw += float(len(evt_df))

    def scale_pot_components(self, factor: float):
        """Multiply POT-weighted fields after deferred MC POT normalization."""
        if factor == 1.0:
            return self
        self.n_signal_pot *= factor
        self.n_truth_nu_pot *= factor
        self.n_total_signal_int *= factor
        self.n_at_stage_int *= factor
        return self

    def __iadd__(self, other: "EfficiencyAccumulator"):
        assert self.var_save_name == other.var_save_name
        assert np.array_equal(self.bins, other.bins)
        self.n_signal_pot += other.n_signal_pot
        self.n_signal_raw += other.n_signal_raw
        self.n_truth_nu_pot += other.n_truth_nu_pot
        self.n_truth_nu_raw += other.n_truth_nu_raw
        self.n_total_signal_int += other.n_total_signal_int
        self.n_total_signal_int_raw += other.n_total_signal_int_raw
        self.n_at_stage_int += other.n_at_stage_int
        self.n_at_stage_int_raw += other.n_at_stage_int_raw
        return self


# ===========================================================================
# Pipeline / Stage / PlotSpec
# ===========================================================================
@dataclass
class PlotSpec:
    """A single plot to make at the END of a stage.

    A ``selector`` is an arbitrary function that receives the current state
    dict ``{"mc":..., "data":..., "intime":..., "dirt":..., "offbeam":...,
    "mc_trk":..., ...}`` and returns the dataframe to histogram for THIS plot
    (per-event or per-track). It can also apply additional inline cuts that
    should NOT propagate to subsequent stages (e.g. the ``cut_nu_like`` blocks
    in the notebook).
    """
    var_config: Any                  # VariableConfig
    breakdown_type: str              # 'topology' | 'genie' | 'genie_sb' | 'pdg'
    selector: Callable[[Dict[str, pd.DataFrame], str], Optional[pd.DataFrame]]
    name_suffix: str = ""            # for save_name disambiguation
    save_kwargs: Dict[str, Any] = field(default_factory=dict)  # kwargs for overlay_hists_from_histdata
    plot_label_template: Optional[Tuple[str, str, str]] = None
    # ``plot_label_template`` overrides the default labels at plot time. Use {pot}
    # placeholder to fill in the POT string at plot time.


@dataclass
class Stage:
    """One stage of the pipeline.

    ``cut`` is the cut function applied to evt-level dataframes, the form
    ``df -> df``. After applying the cut to all samples, the ``plots``
    associated with this stage are evaluated and stored.

    ``before_cut`` and ``after_cut`` hooks let you do bookkeeping that doesn't
    belong in the cut itself (e.g. computing ``prim_trk_*`` derived columns
    from track dfs, or rebuilding mu/p candidates).
    """
    key: str
    label: str                       # human-readable label for summary plots
    cut: Optional[Callable] = None   # state-dict -> state-dict (or None for no cut)
    plots: List[PlotSpec] = field(default_factory=list)
    save_for_efficiency: bool = False  # if True, run the efficiency accumulator at this stage
    save_for_breakdown: bool = False   # if True, run the bar-plot accumulator at this stage


# ===========================================================================
# The chunk-runner: orchestrates one file (one sample) through the pipeline.
# ===========================================================================
class ChunkRunner:
    """Drives one sample's chunk through the pipeline.

    It maintains per-(stage, plot) ``OverlayHistData`` accumulators; for the
    sample being processed it fills only the matching slot
    (e.g. for ``sample='mc'`` it only fills ``hd.mc_hist``).

    After all chunks of all samples have been processed, the per-chunk pickles
    are aggregated by ``aggregate_chunks_to_histdata`` -- which simply calls
    ``+=``.
    """

    def __init__(
        self,
        sample: str,
        stages: List[Stage],
        efficiency_vars: List[Any],
        mc_univ_syst_tags: Optional[Tuple[str, ...]] = None,
    ):
        if sample not in {"mc", "data", "intime", "dirt", "offbeam"}:
            raise ValueError(f"unknown sample: {sample!r}")
        self.sample = sample
        self.stages = stages
        self.efficiency_vars = efficiency_vars
        self.mc_univ_syst_tags = mc_univ_syst_tags
        self._first_eff_stage_key = next(
            (s.key for s in stages if s.save_for_efficiency), None
        )

        # accumulators -- keyed by ((stage_key, plot_name) -> OverlayHistData)
        self.histdata: Dict[Tuple[str, str], OverlayHistData] = {}
        # bar plots (per-stage breakdown counts) for both topology and genie
        self.bar: Dict[str, Dict[str, BarBreakdown]] = {}  # bar[stage_key][breakdown_type]
        self.eff: Dict[str, Dict[str, EfficiencyAccumulator]] = {}  # eff[stage_key][var_save_name]

    @staticmethod
    def plot_key(stage_key: str, ps: PlotSpec) -> str:
        """Stable string identifying a plot across runs."""
        suffix = ("_" + ps.name_suffix) if ps.name_suffix else ""
        return f"{stage_key}__{ps.breakdown_type}__{ps.var_config.var_save_name}{suffix}"

    def _fill_plot(self, stage_key: str, ps: PlotSpec, state: Dict[str, pd.DataFrame]):
        df = ps.selector(state, self.sample)
        if df is None:
            return
        key = (stage_key, self.plot_key(stage_key, ps))
        if key not in self.histdata:
            self.histdata[key] = OverlayHistData(
                var_save_name=ps.var_config.var_save_name,
                breakdown_type=ps.breakdown_type,
                bins=np.asarray(ps.var_config.bins).copy(),
            )
        tags = self.mc_univ_syst_tags if self.sample == "mc" else None
        self.histdata[key].fill_from_df(
            df,
            ps.var_config.var_evt_reco_col,
            self.sample,
            mc_univ_syst_tags=tags,
        )

    def _fill_breakdown(self, stage_key: str, state: Dict[str, pd.DataFrame]):
        if stage_key not in self.bar:
            self.bar[stage_key] = {bt: BarBreakdown.empty(bt) for bt in ("topology", "genie")}
        evt = state.get("evt")
        if self.sample == "mc":
            for bt in ("topology", "genie"):
                self.bar[stage_key][bt].fill_mc(evt)
        elif self.sample == "intime":
            for bt in ("topology", "genie"):
                self.bar[stage_key][bt].fill_intime(evt)
        elif self.sample == "offbeam":
            for bt in ("topology", "genie"):
                self.bar[stage_key][bt].fill_offbeam(evt)
        elif self.sample == "dirt":
            for bt in ("topology", "genie"):
                self.bar[stage_key][bt].fill_dirt(evt)
        elif self.sample == "data":
            for bt in ("topology", "genie"):
                self.bar[stage_key][bt].fill_data(evt)

    def _fill_efficiency(self, stage_key: str, state: Dict[str, pd.DataFrame]):
        if self.sample != "mc":
            return
        # Numerator: reco ``evt`` signal at this stage.
        # Denominator: generated signal on ``mcnu`` when available (preferred); otherwise
        # leave ``n_truth_nu_*`` empty and let render fall back to first-stage evt signal
        # (legacy ``plot_efficiency`` / ``sel_all`` without an ``mcnu`` table).
        mc_df = state.get("evt")
        if mc_df is None:
            return
        mcnu_df = state.get("mcnu")
        if stage_key not in self.eff:
            self.eff[stage_key] = {}
        fill_mcnu = (
            stage_key == self._first_eff_stage_key
            and mcnu_df is not None
            and len(mcnu_df) > 0
        )
        for var_config in self.efficiency_vars:
            if var_config.var_save_name not in self.eff[stage_key]:
                self.eff[stage_key][var_config.var_save_name] = EfficiencyAccumulator.empty(var_config)
            ea = self.eff[stage_key][var_config.var_save_name]
            if fill_mcnu:
                ea.fill_denominator_from_mcnu(mcnu_df, var_config)
            if len(mc_df) > 0:
                sig_mask_evt = IsNuInFV_NumuCC_1p0pi(mc_df, detector=DETECTOR)
                ea.fill_numerator_from_evt(mc_df, var_config, sig_mask_evt)

    def run(
        self,
        initial_state: Dict[str, pd.DataFrame],
        pipeline_trace: Optional[Callable[[str], None]] = None,
    ):
        """Run the pipeline on a single sample's loaded dfs.

        ``initial_state`` is a dict carrying per-event and per-track dfs, e.g.
            ``{"evt": ..., "trk": ..., "hdr": ..., "mcnu": ...}``.
        For MC efficiency, pass ``mcnu`` when the HDF has an ``mcnu_*`` table so the
        denominator uses generated neutrinos. Without ``mcnu`` (e.g. plain
        ``sel_all``), numerators are still filled from ``evt`` and render falls
        back to the first-stage evt signal as the denominator.

        Only entries whose first key matches ``self.sample`` are touched (the
        rest can be set to None when running per-sample).

        ``pipeline_trace`` if set is called with human-readable lines (flush in the
        callback) to localize stalls/crashes when debugging batch jobs.
        """
        def _tr(msg: str) -> None:
            if pipeline_trace is not None:
                pipeline_trace(msg)

        state = dict(initial_state)  # shallow copy; cut funcs do their own copies

        for stage in self.stages:
            _tr(f"[pipeline] >>> stage={stage.key!r} plots={len(stage.plots)} "
                f"breakdown={stage.save_for_breakdown} eff={stage.save_for_efficiency}")
            if stage.cut is not None:
                cname = getattr(stage.cut, "__name__", type(stage.cut).__name__)
                _tr(f"[pipeline]     executing cut {cname!r} …")
                state = stage.cut(state, sample=self.sample)
                ne = len(state["evt"]) if state.get("evt") is not None else 0
                nt = len(state["trk"]) if state.get("trk") is not None else 0
                _tr(f"[pipeline]     after cut: len(evt)={ne} len(trk)={nt}")
            for j, ps in enumerate(stage.plots):
                pk = self.plot_key(stage.key, ps)
                _tr(f"[pipeline]     histogram plot {j + 1}/{len(stage.plots)} → {pk}")
                self._fill_plot(stage.key, ps, state)
            if stage.save_for_breakdown:
                _tr(f"[pipeline]     bar breakdown …")
                self._fill_breakdown(stage.key, state)
            if stage.save_for_efficiency:
                if self.sample == "mc" and state.get("mcnu") is None:
                    _tr("[pipeline]     efficiency accumulators (evt-only; no mcnu) …")
                else:
                    _tr("[pipeline]     efficiency accumulators …")
                self._fill_efficiency(stage.key, state)
            _tr(f"[pipeline] <<< stage={stage.key!r} finished")

    # -------- save/load -------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample": self.sample,
            "stage_keys": [s.key for s in self.stages],
            "stage_labels": [s.label for s in self.stages],
            "histdata": self.histdata,
            "bar": self.bar,
            "eff": self.eff,
        }

    def save(self, out_path: str, extra_meta: Optional[Dict[str, Any]] = None):
        out = self.to_dict()
        if extra_meta:
            out["meta"] = extra_meta
        makedirs(path.dirname(out_path), exist_ok=True)
        with open(out_path, "wb") as f:
            pickle.dump(out, f)


# ===========================================================================
# Aggregation helpers
# ===========================================================================
def _bins_fingerprint(bins) -> Tuple[float, ...]:
    return tuple(np.asarray(bins, dtype=float).ravel().tolist())


def _majority_bins_fingerprint(fingerprints: List[Tuple[float, ...]]) -> Tuple[float, ...]:
    """Return the most common bins fingerprint (stable tie-break: first seen)."""
    if not fingerprints:
        raise ValueError("empty fingerprints")
    counts = Counter(fingerprints)
    best_n = max(counts.values())
    for fp in fingerprints:
        if counts[fp] == best_n:
            return fp
    return fingerprints[0]


def _warn_bin_skip(context: str, key, n_keep: int, n_skip: int, bins_keep) -> None:
    if n_skip <= 0:
        return
    nb = max(len(bins_keep) - 1, 0)
    print(
        f"[aggregate] WARN: {context} {key!r}: skipped {n_skip} contribution(s) "
        f"with mismatched bins (kept {n_keep} @ nbins={nb}). "
        f"Re-run map jobs so all samples share the same VariableConfig bins.",
        flush=True,
    )


def aggregate_chunk_files(chunk_files: List[str]) -> Dict[str, Any]:
    """Sum the contents of multiple per-chunk pickles.

    Pickles must be from the SAME sample. Returns a single dict-of-accumulators
    with histograms summed bin-by-bin.

    If VariableConfig bin edges changed mid-campaign, chunks with a minority
    bin scheme for a given plot are skipped (majority wins) so live aggregation
    does not assert-fail.
    """
    if not chunk_files:
        raise ValueError("no chunk files passed to aggregate_chunk_files")

    loaded: List[Dict[str, Any]] = []
    for cf in chunk_files:
        with open(cf, "rb") as f:
            loaded.append(pickle.load(f))

    sample0 = loaded[0]["sample"]
    for d in loaded[1:]:
        assert d["sample"] == sample0

    # Majority bin scheme per histdata / eff key across chunks.
    hist_fps: Dict[Any, List[Tuple[float, ...]]] = defaultdict(list)
    eff_fps: Dict[Tuple[str, str], List[Tuple[float, ...]]] = defaultdict(list)
    for d in loaded:
        for key, hd in d["histdata"].items():
            hist_fps[key].append(_bins_fingerprint(hd.bins))
        for stage_key, by_v in (d.get("eff") or {}).items():
            for v, ea in by_v.items():
                eff_fps[(stage_key, v)].append(_bins_fingerprint(ea.bins))

    hist_pref = {k: _majority_bins_fingerprint(fps) for k, fps in hist_fps.items()}
    eff_pref = {k: _majority_bins_fingerprint(fps) for k, fps in eff_fps.items()}

    out: Optional[Dict[str, Any]] = None
    hist_keep: Dict[Any, int] = defaultdict(int)
    hist_skip: Dict[Any, int] = defaultdict(int)
    eff_keep: Dict[Tuple[str, str], int] = defaultdict(int)
    eff_skip: Dict[Tuple[str, str], int] = defaultdict(int)

    for d in loaded:
        if out is None:
            out = {
                "sample": d["sample"],
                "stage_keys": d["stage_keys"],
                "stage_labels": d["stage_labels"],
                "histdata": {},
                "bar": {},
                "eff": {},
                "meta": d.get("meta"),
            }
        # histdata
        for key, hd in d["histdata"].items():
            fp = _bins_fingerprint(hd.bins)
            if fp != hist_pref[key]:
                hist_skip[key] += 1
                continue
            hist_keep[key] += 1
            if key in out["histdata"]:
                out["histdata"][key] += hd
            else:
                out["histdata"][key] = hd
        # bar (no bins)
        for stage_key, by_bt in d["bar"].items():
            if stage_key not in out["bar"]:
                out["bar"][stage_key] = dict(by_bt)
            else:
                for bt, bb in by_bt.items():
                    if bt in out["bar"][stage_key]:
                        out["bar"][stage_key][bt] += bb
                    else:
                        out["bar"][stage_key][bt] = bb
        # eff
        for stage_key, by_v in (d.get("eff") or {}).items():
            if stage_key not in out["eff"]:
                out["eff"][stage_key] = {}
            for v, ea in by_v.items():
                ek = (stage_key, v)
                fp = _bins_fingerprint(ea.bins)
                if fp != eff_pref[ek]:
                    eff_skip[ek] += 1
                    continue
                eff_keep[ek] += 1
                if v in out["eff"][stage_key]:
                    out["eff"][stage_key][v] += ea
                else:
                    out["eff"][stage_key][v] = ea

    for key, n_skip in hist_skip.items():
        seed = out["histdata"].get(key)
        _warn_bin_skip(
            f"sample={sample0!r} histdata",
            key,
            hist_keep.get(key, 0),
            n_skip,
            seed.bins if seed is not None else [],
        )
    for ek, n_skip in eff_skip.items():
        stage_key, v = ek
        seed = (out["eff"].get(stage_key) or {}).get(v)
        _warn_bin_skip(
            f"sample={sample0!r} eff",
            ek,
            eff_keep.get(ek, 0),
            n_skip,
            seed.bins if seed is not None else [],
        )
    return out


def merge_samples(samples: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Combine per-sample aggregated results into one set of OverlayHistData.

    ``samples`` maps sample-name -> aggregated-dict (output of
    aggregate_chunk_files for that sample). The output's ``histdata`` is a
    dict (stage_key, plot_key) -> OverlayHistData, with all sample slots
    populated.

    Intime and offbeam cosmics are kept in separate arrays until plotting so the
    central prediction can be chosen at aggregation time.

    Samples whose bin edges disagree with the majority scheme for a plot are
    skipped for that plot (with a warning) so mixed VariableConfig campaigns
    still render.
    """
    # union of all (stage_key, plot_key) entries across samples
    all_keys = set()
    for s in samples.values():
        all_keys.update(s["histdata"].keys())

    merged_histdata: Dict[Tuple[str, str], OverlayHistData] = {}
    for key in all_keys:
        fps = []
        holders = []
        for sample_name, s in samples.items():
            if key not in s["histdata"]:
                continue
            hd = s["histdata"][key]
            fps.append(_bins_fingerprint(hd.bins))
            holders.append((sample_name, hd))
        if not holders:
            continue
        pref = _majority_bins_fingerprint(fps)
        keep = [(n, hd) for (n, hd), fp in zip(holders, fps) if fp == pref]
        skip_n = len(holders) - len(keep)
        seed = keep[0][1]
        _warn_bin_skip("merge_samples histdata", key, len(keep), skip_n, seed.bins)
        merged = OverlayHistData(
            var_save_name=seed.var_save_name,
            breakdown_type=seed.breakdown_type,
            bins=seed.bins.copy(),
        )
        for _sample_name, hd in keep:
            merged += hd
        merged_histdata[key] = merged

    # bar / eff: only mc samples have bar mc-counts, intime/dirt have intime/dirt
    # counts; merge by adding stage-by-stage.
    merged_bar: Dict[str, Dict[str, BarBreakdown]] = {}
    for s in samples.values():
        for stage_key, by_bt in s["bar"].items():
            merged_bar.setdefault(stage_key, {})
            for bt, bb in by_bt.items():
                if bt in merged_bar[stage_key]:
                    merged_bar[stage_key][bt] += bb
                else:
                    # take a copy to avoid mutating someone else's data
                    merged_bar[stage_key][bt] = BarBreakdown(
                        breakdown_type=bb.breakdown_type,
                        mc_counts=bb.mc_counts.copy(),
                        intime_count=bb.intime_count,
                        offbeam_count=bb.offbeam_count,
                        dirt_count=bb.dirt_count,
                        data_count=bb.data_count,
                    )

    # eff: majority bins per (stage, var)
    eff_votes: Dict[Tuple[str, str], List[Tuple[float, ...]]] = defaultdict(list)
    for s in samples.values():
        for stage_key, by_v in (s.get("eff") or {}).items():
            for v, ea in by_v.items():
                eff_votes[(stage_key, v)].append(_bins_fingerprint(ea.bins))
    eff_pref = {k: _majority_bins_fingerprint(fps) for k, fps in eff_votes.items()}

    merged_eff: Dict[str, Dict[str, EfficiencyAccumulator]] = {}
    for s in samples.values():
        for stage_key, by_v in (s.get("eff") or {}).items():
            merged_eff.setdefault(stage_key, {})
            for v, ea in by_v.items():
                ek = (stage_key, v)
                if _bins_fingerprint(ea.bins) != eff_pref[ek]:
                    _warn_bin_skip(
                        "merge_samples eff",
                        ek,
                        0,
                        1,
                        np.asarray(eff_pref[ek]),
                    )
                    continue
                if v in merged_eff[stage_key]:
                    merged_eff[stage_key][v] += ea
                else:
                    zpot = np.zeros(len(ea.bins) - 1)
                    merged_eff[stage_key][v] = EfficiencyAccumulator(
                        var_save_name=ea.var_save_name,
                        bins=ea.bins.copy(),
                        n_signal_pot=ea.n_signal_pot.copy(),
                        n_signal_raw=ea.n_signal_raw.copy(),
                        n_truth_nu_pot=np.asarray(
                            getattr(ea, "n_truth_nu_pot", zpot), dtype=float
                        ).copy(),
                        n_truth_nu_raw=np.asarray(
                            getattr(ea, "n_truth_nu_raw", zpot), dtype=float
                        ).copy(),
                        n_total_signal_int=ea.n_total_signal_int,
                        n_total_signal_int_raw=ea.n_total_signal_int_raw,
                        n_at_stage_int=ea.n_at_stage_int,
                        n_at_stage_int_raw=getattr(ea, "n_at_stage_int_raw", 0.0),
                    )

    # Assume stage_keys/labels are the same across samples (they should be: same pipeline)
    any_sample = next(iter(samples.values()))
    return {
        "stage_keys": any_sample["stage_keys"],
        "stage_labels": any_sample["stage_labels"],
        "histdata": merged_histdata,
        "bar": merged_bar,
        "eff": merged_eff,
    }


def sanitize_merged_histdata_finite(merged: Dict[str, Any]) -> Tuple[int, int]:
    """Zero out NaN/inf bins in merged histogram arrays.

    ``numpy.histogram(..., weights=...)`` turns **every** bin into NaN if any event
    weight is NaN; stale pickles from before genweight sanitization need this pass so
    plots and exposure scaling stay finite.

    Returns (n_overlay_histdata_patched, n_bar_breakdowns_patched).
    """
    n_hd = 0
    for hd in merged["histdata"].values():
        patched_here = False
        for attr in (
            "mc_hist",
            "mc_err2",
            "intime_hist",
            "intime_err2",
            "offbeam_hist",
            "offbeam_err2",
            "dirt_hist",
            "dirt_err2",
            "data_hist",
            "data_err2",
        ):
            arr = getattr(hd, attr)
            if arr is None:
                continue
            if np.issubdtype(arr.dtype, np.number) and not np.all(np.isfinite(arr)):
                setattr(
                    hd,
                    attr,
                    np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0),
                )
                patched_here = True
        if patched_here:
            n_hd += 1

    n_bar = 0
    for by_bt in merged.get("bar", {}).values():
        for bb in by_bt.values():
            bar_dirty = False
            if not np.all(np.isfinite(bb.mc_counts)):
                bb.mc_counts = np.nan_to_num(
                    bb.mc_counts, nan=0.0, posinf=0.0, neginf=0.0
                )
                bar_dirty = True
            for fld in ("intime_count", "offbeam_count", "dirt_count", "data_count"):
                v = getattr(bb, fld)
                if not np.isfinite(v):
                    setattr(bb, fld, 0.0)
                    bar_dirty = True
            if bar_dirty:
                n_bar += 1

    return n_hd, n_bar


@dataclass
class ExposureTotals:
    """Sums of POT / gates across every chunk file (filled before global scaling)."""

    data_pot: float = 0.0
    data_gates_bnb: float = 0.0
    mc_pot: float = 0.0
    dirt_pot: float = 0.0
    intime_gates: float = 0.0
    offbeam_gates: float = 0.0


def apply_global_exposure_scales(
    merged: Dict[str, Any],
    totals: ExposureTotals,
    f_offbeam_coincident: float = 0.075,
) -> Dict[str, float]:
    """Scale histograms / breakdown / MC efficiency after summing intrinsic-weight chunks.

    Mutates ``merged`` in place (histdata, bar, eff). Data histograms stay raw counts.

    Returns the dict of scale factors applied for logging.
    """
    def _safe_ratio(num: float, den: float, default: float) -> float:
        # No target exposure yet (e.g. live run before any data chunk) → keep
        # intrinsic weights instead of scaling MC/dirt to zero.
        if num <= 0 or den <= 0 or not np.isfinite(num) or not np.isfinite(den):
            return default
        r = num / den
        return float(r) if np.isfinite(r) else default

    sm = {
        "scale_mc": _safe_ratio(totals.data_pot, totals.mc_pot, 1.0),
        "scale_dirt": _safe_ratio(totals.data_pot, totals.dirt_pot, 1.0),
    }
    gate_fac = (1.0 - f_offbeam_coincident)
    # Cosmic scales default to 0 when on-beam gates are missing (nothing to normalize to).
    sm["scale_intime"] = _safe_ratio(
        gate_fac * totals.data_gates_bnb, totals.intime_gates, 0.0
    )
    sm["scale_offbeam"] = _safe_ratio(
        gate_fac * totals.data_gates_bnb, totals.offbeam_gates, 0.0
    )

    for hd in merged["histdata"].values():
        if hd.has_mc:
            hd.mc_hist *= sm["scale_mc"]
            hd.mc_err2 *= sm["scale_mc"] ** 2
        if hd.has_intime:
            hd.intime_hist *= sm["scale_intime"]
            hd.intime_err2 *= sm["scale_intime"] ** 2
        if hd.has_offbeam:
            hd.offbeam_hist *= sm["scale_offbeam"]
            hd.offbeam_err2 *= sm["scale_offbeam"] ** 2
        if hd.has_dirt:
            hd.dirt_hist *= sm["scale_dirt"]
            hd.dirt_err2 *= sm["scale_dirt"] ** 2
        if hd.has_mc and getattr(hd, "mc_univ_hist", None):
            sf = sm["scale_mc"]
            for _k in hd.mc_univ_hist:
                hd.mc_univ_hist[_k] = np.asarray(hd.mc_univ_hist[_k], dtype=float) * sf

    for by_bt in merged["bar"].values():
        for bb in by_bt.values():
            bb.mc_counts *= sm["scale_mc"]
            bb.intime_count *= sm["scale_intime"]
            bb.offbeam_count *= sm["scale_offbeam"]
            bb.dirt_count *= sm["scale_dirt"]

    sf = sm["scale_mc"]
    if sf != 1.0:
        for by_v in merged["eff"].values():
            for ea in by_v.values():
                ea.scale_pot_components(sf)

    return sm
