"""In-memory multi-sample helper for the legacy event-selection notebook.

Holds event / track / header dataframes for ``mc``, ``data``, ``intime``,
``offbeam``, and ``dirt``, and applies the same cut or transform to every
sample in one call so notebook cells can focus on cut thresholds.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from os import path
from typing import Any, Callable, Dict, Iterable, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from analysis_village.numucc_1p0pi.makedf.selections import SAVE_NTRKS, get_trk_info, get_valid_trks
from makedf.util import avg_chi2, match_trkdf_to_slcdf

SAMPLES: Tuple[str, ...] = ("mc", "data", "intime", "offbeam", "dirt")
EvtFn = Callable[..., pd.DataFrame]
TrkFn = Callable[..., pd.DataFrame]


def _ensure_samples(samples: Optional[Iterable[str]] = None) -> Tuple[str, ...]:
    if samples is None:
        return SAMPLES
    return tuple(samples)


def _own(df: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
    """Return an owned DataFrame copy.

    Selection cuts typically return views (``df[mask]``). Storing those views and
    later assigning columns triggers ``SettingWithCopyWarning`` and can silently
    drop writes. Always keep owned frames inside ``SampleBundle``.
    """
    if df is None:
        return None
    return df.copy()


@dataclass
class SampleBundle:
    """Container for the five analysis samples used by the legacy notebook.

    Typical stage cell::

        samples.apply_evt(cut_nu_score)  # defaults from selections.py
        samples.refresh_tracks(attach_ntrks=True)
        samples.record_stage("nu_score", bars=True, breakdown_dict=breakdown_dict,
                             plot_bar_plots=plot_bar_plots, show_plot=show_plot)

    All cut / transform entry points store **owned** copies so column writes after
    filtering never hit a pandas view.
    """

    evt: Dict[str, pd.DataFrame] = field(default_factory=dict)
    trk: Dict[str, pd.DataFrame] = field(default_factory=dict)
    hdr: Dict[str, pd.DataFrame] = field(default_factory=dict)
    stages: Dict[str, Dict[str, pd.DataFrame]] = field(default_factory=dict)
    pot_str: str = "dummy"
    data_tot_pot: float = 0.0

    def __post_init__(self) -> None:
        if not self.stages:
            self.stages = {s: {} for s in SAMPLES}

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_dfs_dicts(
        cls,
        *,
        mc: Mapping[str, pd.DataFrame],
        data: Mapping[str, pd.DataFrame],
        intime: Mapping[str, pd.DataFrame],
        offbeam: Mapping[str, pd.DataFrame],
        dirt: Mapping[str, pd.DataFrame],
    ) -> "SampleBundle":
        loaded = {"mc": mc, "data": data, "intime": intime, "offbeam": offbeam, "dirt": dirt}
        return cls(
            evt={s: _own(loaded[s]["evt"]) for s in SAMPLES},
            trk={s: _own(loaded[s]["trk"]) for s in SAMPLES},
            hdr={s: _own(loaded[s]["hdr"]) for s in SAMPLES},
        )

    @classmethod
    def load_from_dirs(
        cls,
        base_dir: str,
        sample_dirs: Mapping[str, Tuple[str, Sequence[str]]],
        *,
        filename_str: str = "sel_all",
        n_max_concat: Optional[int] = None,
        dfs_from_dir: Optional[Callable[..., Mapping[str, pd.DataFrame]]] = None,
    ) -> "SampleBundle":
        """Load all samples via ``dfs_from_dir``.

        ``sample_dirs`` maps sample name → ``(subdir, keys2load)``.
        """
        if dfs_from_dir is None:
            from pyanalib.split_df_helpers_new import dfs_from_dir as _dfs_from_dir
            dfs_from_dir = _dfs_from_dir

        loaded: Dict[str, Mapping[str, pd.DataFrame]] = {}
        for sample, (subdir, keys) in sample_dirs.items():
            loaded[sample] = dfs_from_dir(
                search_dir=path.join(base_dir, subdir),
                filename_str=filename_str,
                keys2load=list(keys),
                n_max_concat=n_max_concat,
            )
        missing = [s for s in SAMPLES if s not in loaded]
        if missing:
            raise KeyError(f"load_from_dirs missing samples: {missing}")
        return cls.from_dfs_dicts(
            mc=loaded["mc"],
            data=loaded["data"],
            intime=loaded["intime"],
            offbeam=loaded["offbeam"],
            dirt=loaded["dirt"],
        )

    # ------------------------------------------------------------------
    # Exposure / weights
    # ------------------------------------------------------------------
    def _set_evt_columns(self, sample: str, **columns: Any) -> None:
        """Assign columns on an owned copy of ``evt[sample]``."""
        evt = _own(self.evt[sample])
        for col, value in columns.items():
            evt[col] = value
        self.evt[sample] = evt

    def _set_trk_columns(self, sample: str, **columns: Any) -> None:
        """Assign columns on an owned copy of ``trk[sample]``."""
        trk = self.trk.get(sample)
        if trk is None:
            return
        trk = _own(trk)
        for col, value in columns.items():
            trk[col] = value
        self.trk[sample] = trk

    def _sync_trk_exposure_weights(
        self, samples: Optional[Iterable[str]] = None
    ) -> "SampleBundle":
        """Copy ``pot_weight`` / ``gates_weight`` from evt onto trk (event→track).

        Track-level plotters (``overlay_hists`` on concat'd trk1/trk2) need these
        columns. They are lost when tracks are rebuilt from ``evt.trk1``/``trk2``.
        """
        for s in _ensure_samples(samples):
            evt = self.evt.get(s)
            trk = self.trk.get(s)
            if evt is None or trk is None or len(trk) == 0:
                continue
            trk = _own(trk)
            n_evt = len(evt.index.names)
            n_trk = len(trk.index.names)
            if n_trk > n_evt:
                evt_key = trk.index.droplevel(list(range(n_evt, n_trk)))
            else:
                evt_key = trk.index
            for col in ("pot_weight", "gates_weight"):
                if col not in evt.columns:
                    continue
                trk[col] = evt[col].reindex(evt_key).to_numpy()
            self.trk[s] = trk
        return self

    def assign_exposure(self, f: float = 0.08) -> "SampleBundle":
        """Attach POT / gates weights so MC and cosmics match on-beam data.

        Weights are written on both ``evt`` and ``trk`` (plotters read
        ``df.pot_weight`` for track-level overlays too).
        """
        from analysis_village.numucc_1p0pi.utils import get_pot_str

        data_hdr = self.hdr["data"]
        data_tot_pot = float(data_hdr["pot"].sum())
        self.data_tot_pot = data_tot_pot
        self.pot_str = get_pot_str(data_tot_pot)
        self._set_evt_columns("data", pot_weight=1.0)
        self._set_trk_columns("data", pot_weight=1.0)
        data_gates = float(data_hdr.nbnbinfo.sum())
        print(f"data_tot_pot: {data_tot_pot:.3e}  data tot gates: {data_gates:.3e}")

        for sample, pot_key in (("mc", "pot"), ("dirt", "pot")):
            tot = float(self.hdr[sample][pot_key].sum())
            scale = data_tot_pot / tot if tot > 0 else 0.0
            self._set_evt_columns(sample, pot_weight=scale)
            self._set_trk_columns(sample, pot_weight=scale)
            print(f"{sample}_tot_pot: {tot:.3e}  {sample}_pot_scale: {scale:.3e}")

        offbeam_hdr = self.hdr["offbeam"]
        offbeam_gates = float(
            offbeam_hdr.loc[offbeam_hdr["first_in_subrun"] == 1, "noffbeambnb"].sum()
        )
        scale_offbeam = (1 - f) * data_gates / offbeam_gates if offbeam_gates > 0 else 0.0
        self._set_evt_columns("offbeam", gates_weight=scale_offbeam, pot_weight=scale_offbeam)
        self._set_trk_columns("offbeam", gates_weight=scale_offbeam, pot_weight=scale_offbeam)
        print(f"offbeam cosmics data gates: {offbeam_gates:.2e}  goal scale: {scale_offbeam:.2f}")

        intime_hdr = self.hdr["intime"]
        intime_gates = float(
            intime_hdr.loc[intime_hdr["first_in_subrun"] == 1, "ngenevt"].sum()
        )
        scale_intime = (1 - f) * data_gates / intime_gates if intime_gates > 0 else 0.0
        self._set_evt_columns("intime", gates_weight=scale_intime, pot_weight=scale_intime)
        self._set_trk_columns("intime", gates_weight=scale_intime, pot_weight=scale_intime)
        print(f"intime cosmics data gates: {intime_gates:.2e}  goal scale: {scale_intime:.2f}")
        return self

    # ------------------------------------------------------------------
    # Uniform transforms
    # ------------------------------------------------------------------
    def apply_evt(
        self,
        fn: EvtFn,
        *args: Any,
        samples: Optional[Iterable[str]] = None,
        **kwargs: Any,
    ) -> "SampleBundle":
        for s in _ensure_samples(samples):
            self.evt[s] = _own(fn(self.evt[s], *args, **kwargs))
        return self

    def apply_trk(
        self,
        fn: TrkFn,
        *args: Any,
        samples: Optional[Iterable[str]] = None,
        **kwargs: Any,
    ) -> "SampleBundle":
        for s in _ensure_samples(samples):
            self.trk[s] = _own(fn(self.trk[s], *args, **kwargs))
        return self

    def apply_paired(
        self,
        fn: Callable[[pd.DataFrame, pd.DataFrame], Tuple[pd.DataFrame, pd.DataFrame]],
        samples: Optional[Iterable[str]] = None,
    ) -> "SampleBundle":
        """Apply ``fn(evt, trk) -> (evt, trk)`` to every sample."""
        for s in _ensure_samples(samples):
            evt, trk = fn(self.evt[s], self.trk[s])
            self.evt[s], self.trk[s] = _own(evt), _own(trk)
        return self

    def map_evt_trk(
        self,
        fn: Callable[[str, pd.DataFrame, pd.DataFrame], None],
        samples: Optional[Iterable[str]] = None,
    ) -> "SampleBundle":
        """Call ``fn(sample, evt, trk)`` for side effects (e.g. column writes).

        ``fn`` must treat the passed frames as owned (or copy before writing).
        Prefer ``apply_evt`` / ``apply_trk`` / the ``attach_*`` helpers instead.
        """
        for s in _ensure_samples(samples):
            fn(s, self.evt[s], self.trk[s])
        return self

    # ------------------------------------------------------------------
    # Stage bookkeeping
    # ------------------------------------------------------------------
    def record_stage(
        self,
        key: str,
        *,
        bars: bool = False,
        breakdown_dict: Optional[Dict[str, dict]] = None,
        plot_bar_plots: Optional[Callable[..., dict]] = None,
        show_plot: bool = True,
        samples: Optional[Iterable[str]] = None,
    ) -> "SampleBundle":
        for s in _ensure_samples(samples):
            # Snapshot an owned copy so later in-place updates to self.evt cannot
            # retroactively change efficiency / breakdown history.
            self.stages[s][key] = _own(self.evt[s])
        if bars:
            if plot_bar_plots is None:
                raise ValueError("bars=True requires plot_bar_plots")
            ret = plot_bar_plots(
                key, self.evt["mc"], self.evt["intime"], self.evt["dirt"], show_plot=show_plot
            )
            if breakdown_dict is not None:
                for bt in ("topology", "genie"):
                    breakdown_dict[bt][key] = ret[bt]["perc_list"]
        return self

    def cut_stage(
        self,
        key: str,
        fn: EvtFn,
        *args: Any,
        bars: bool = True,
        breakdown_dict: Optional[Dict[str, dict]] = None,
        plot_bar_plots: Optional[Callable[..., dict]] = None,
        show_plot: bool = True,
        **kwargs: Any,
    ) -> "SampleBundle":
        """Apply an event cut and optionally record + bar-plot the stage."""
        self.apply_evt(fn, *args, **kwargs)
        return self.record_stage(
            key,
            bars=bars,
            breakdown_dict=breakdown_dict,
            plot_bar_plots=plot_bar_plots,
            show_plot=show_plot,
        )

    def run_pipeline(
        self,
        *,
        bars: bool = True,
        breakdown_dict: Optional[Dict[str, dict]] = None,
        plot_bar_plots: Optional[Callable[..., dict]] = None,
        show_plot: bool = True,
        samples: Optional[Iterable[str]] = None,
        stop_at: Optional[str] = None,
    ) -> "SampleBundle":
        """Apply :func:`build_pipeline` cuts to every sample (same chain as CAF / syst).

        Records a stage snapshot after each pipeline stage. Prefer this over
        hand-written ``cut_stage`` / ``apply_evt`` sequences so thresholds and
        cut order stay in sync with ``event_selection_pipeline_def``.
        """
        from analysis_village.numucc_1p0pi.event_selection_pipeline_def import (
            build_pipeline,
        )

        sample_names = _ensure_samples(samples)
        for stage in build_pipeline():
            if stage.cut is not None:
                for s in sample_names:
                    state = {
                        "evt": self.evt[s],
                        "trk": self.trk.get(s),
                        "hdr": self.hdr.get(s),
                    }
                    state = stage.cut(state, sample=s)
                    self.evt[s] = _own(state.get("evt"))
                    if state.get("trk") is not None:
                        self.trk[s] = _own(state["trk"])
                self._sync_trk_exposure_weights(samples=sample_names)
            self.record_stage(
                stage.key,
                bars=bars and stage.save_for_breakdown,
                breakdown_dict=breakdown_dict,
                plot_bar_plots=plot_bar_plots,
                show_plot=show_plot,
                samples=sample_names,
            )
            if stop_at is not None and stage.key == stop_at:
                break
        return self

    # ------------------------------------------------------------------
    # Track / column helpers (shared across stages)
    # ------------------------------------------------------------------
    def refresh_tracks(
        self,
        *,
        attach_ntrks: bool = False,
        save_ntrks: int = SAVE_NTRKS,
        samples: Optional[Iterable[str]] = None,
    ) -> "SampleBundle":
        """``get_valid_trks`` → match to current evt; optionally attach ``trk1``/``trk2``."""
        for s in _ensure_samples(samples):
            trk = _own(get_valid_trks(self.trk[s]))
            trk = _own(match_trkdf_to_slcdf(trk, self.evt[s]))
            self.trk[s] = trk
            if attach_ntrks:
                self.evt[s] = _own(get_trk_info(self.evt[s], trk, save_ntrks))
        return self._sync_trk_exposure_weights(samples=samples)

    def attach_chi2_avgs(self, samples: Optional[Iterable[str]] = None) -> "SampleBundle":
        for s in _ensure_samples(samples):
            trk = _own(self.trk[s])
            trk[("pfp", "trk", "chi2pid", "avg", "chi2_muon", "")] = avg_chi2(trk, "chi2_muon")
            trk[("pfp", "trk", "chi2pid", "avg", "chi2_proton", "")] = avg_chi2(trk, "chi2_proton")
            self.trk[s] = trk
        return self

    def attach_mcs_range_diff(self, samples: Optional[Iterable[str]] = None) -> "SampleBundle":
        for s in _ensure_samples(samples):
            trk = _own(self.trk[s])
            diff = (trk.pfp.trk.rangeP.p_muon - trk.pfp.trk.mcsP.fwdP_muon) / trk.pfp.trk.rangeP.p_muon
            trk[("pfp", "trk", "mcs_range_diff", "", "", "")] = diff
            self.trk[s] = trk
        return self

    def attach_prim_trk_cols(
        self,
        cols: Optional[Sequence[str]] = None,
        samples: Optional[Iterable[str]] = None,
    ) -> "SampleBundle":
        """Copy quantities from the longest track onto the event dataframe.

        Default columns: ``phi``, ``dir_y``, ``start_x``, ``end_x``, ``P_frac_diff``.
        """
        wanted = set(cols) if cols is not None else {
            "phi", "dir_y", "start_x", "end_x", "P_frac_diff",
        }
        for s in _ensure_samples(samples):
            trk = self.trk[s]
            if trk is None or len(trk) == 0:
                continue
            trk = _own(trk)
            evt = _own(self.evt[s])
            trk[("pfp", "trk", "phi", "", "", "")] = np.degrees(np.arctan2(
                trk["pfp", "trk", "dir", "x", ""],
                trk["pfp", "trk", "dir", "y", ""],
            ))
            nlevels = len(trk.index.names)
            prim = (
                trk.sort_values(("pfp", "trk", "len", "", "", ""), ascending=False)
                   .groupby(level=list(range(nlevels - 1)))
                   .head(1)
                   .reset_index(level=[nlevels - 1], drop=True)
            )
            mapping = {
                "phi": ("prim_trk_phi", prim.pfp.trk.phi),
                "dir_y": ("prim_trk_dir_y", prim.pfp.trk.dir.y),
                "start_x": ("prim_trk_start_x", prim.pfp.trk.start.x),
                "end_x": ("prim_trk_end_x", prim.pfp.trk.end.x),
            }
            for key, (col, series) in mapping.items():
                if key in wanted:
                    evt[col] = series.reindex(evt.index)
            if "P_frac_diff" in wanted:
                try:
                    frac = (
                        (prim.pfp.trk.rangeP.p_muon - prim.pfp.trk.mcsP.fwdP_muon)
                        / prim.pfp.trk.rangeP.p_muon
                    )
                    evt["prim_trk_P_frac_diff"] = frac.reindex(evt.index)
                except Exception:
                    pass
            self.trk[s] = trk
            self.evt[s] = evt
        return self

    def concat_trk12(self, samples: Optional[Iterable[str]] = None) -> Dict[str, pd.DataFrame]:
        """``pd.concat([evt.trk1, evt.trk2])`` per sample (for track-level plots)."""
        out: Dict[str, pd.DataFrame] = {}
        for s in _ensure_samples(samples):
            evt = self.evt[s]
            out[s] = pd.concat([evt.trk1, evt.trk2])
        return out

    def set_trk_from_concat(self, samples: Optional[Iterable[str]] = None) -> "SampleBundle":
        for s, trk in self.concat_trk12(samples=samples).items():
            self.trk[s] = _own(trk)
        return self._sync_trk_exposure_weights(samples=samples)

    # ------------------------------------------------------------------
    # Plotting adapters
    # ------------------------------------------------------------------
    def as_overlay_kwargs(
        self,
        source: Union[str, Mapping[str, pd.DataFrame]] = "evt",
        *,
        include_offbeam: bool = False,
        selector: Optional[Callable[[pd.DataFrame], pd.DataFrame]] = None,
    ) -> Dict[str, pd.DataFrame]:
        """Keyword args for ``overlay_hists`` / ``bar_plot``.

        ``source``: ``"evt"``, ``"trk"``, ``"trk12"``, or an explicit sample→df map.
        ``selector``: optional per-df transform (e.g. ``lambda df: df.trk1``).
        When a selector drops ``pot_weight``, it is re-attached from the parent evt.
        """
        if isinstance(source, str):
            if source == "evt":
                src = self.evt
            elif source == "trk":
                src = self.trk
            elif source == "trk12":
                src = self.concat_trk12()
            else:
                raise ValueError(f"unknown source {source!r}")
        else:
            src = source

        def _pick(name: str) -> pd.DataFrame:
            parent = src[name]
            df = selector(parent) if selector is not None else parent
            # Track blocks (trk1/trk2) and fresh concat'd tracks may lack pot_weight.
            if df is not None and "pot_weight" not in getattr(df, "columns", []):
                evt = self.evt.get(name)
                if evt is not None and "pot_weight" in evt.columns and len(df) > 0:
                    df = df.copy()
                    n_evt = len(evt.index.names)
                    n_df = len(df.index.names)
                    evt_key = (
                        df.index.droplevel(list(range(n_evt, n_df)))
                        if n_df > n_evt
                        else df.index
                    )
                    df["pot_weight"] = evt["pot_weight"].reindex(evt_key).to_numpy()
            return df

        kw = {
            "mc_df": _pick("mc"),
            "data_df": _pick("data"),
            "intime_df": _pick("intime"),
            "dirt_df": _pick("dirt"),
        }
        if include_offbeam:
            kw["offbeam_df"] = _pick("offbeam")
        return kw

    # ------------------------------------------------------------------
    # Convenience aliases (notebook / efficiency)
    # ------------------------------------------------------------------
    @property
    def df_dict(self) -> Dict[str, pd.DataFrame]:
        """MC stage snapshots (efficiency curves)."""
        return self.stages["mc"]

    @property
    def df_dict_data(self) -> Dict[str, pd.DataFrame]:
        return self.stages["data"]

    @property
    def df_dict_intime(self) -> Dict[str, pd.DataFrame]:
        return self.stages["intime"]

    @property
    def df_dict_offbeam(self) -> Dict[str, pd.DataFrame]:
        return self.stages["offbeam"]

    @property
    def df_dict_dirt(self) -> Dict[str, pd.DataFrame]:
        return self.stages["dirt"]

    def evt_tuple(self) -> Tuple[pd.DataFrame, ...]:
        return tuple(self.evt[s] for s in SAMPLES)

    def trk_tuple(self) -> Tuple[pd.DataFrame, ...]:
        return tuple(self.trk[s] for s in SAMPLES)

    def unpack_evt(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        return self.evt["mc"], self.evt["data"], self.evt["intime"], self.evt["offbeam"], self.evt["dirt"]

    def unpack_trk(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        return self.trk["mc"], self.trk["data"], self.trk["intime"], self.trk["offbeam"], self.trk["dirt"]
