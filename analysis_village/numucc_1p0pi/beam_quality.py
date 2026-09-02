"""Beam-quality cuts for SBND 1e20 data (from ``notebooks/beam_quality.ipynb``).

Identifies bad triggered events from:
- bad beam spills (FOM = 0 or 4, or 0 < FOM < ``fom_cut``)
- short runs (duration < ``min_run_duration_min``)

Returns ``evt_good`` and a hdr table restricted to good triggered events.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

BEAM_COLS = ("TOR875", "TOR860", "FOM", "THCURR", "spill_time")


@dataclass
class BeamQualitySummary:
    n_evt_prefilter: int
    n_evt_good: int
    n_hdr_total: int
    n_hdr_bad: int
    n_short_runs: int
    pot_prefilter: float
    pot_good: float


def match_beam_spills(data_hdr_df: pd.DataFrame, pot_df: pd.DataFrame) -> pd.DataFrame:
    """Attach nearest preceding spill quantities to each triggered hdr row."""
    hdr = data_hdr_df.copy()
    pot_reset = pot_df.reset_index()
    pot_reset["spill_time"] = (
        pot_reset["spill_time_sec"] + pot_reset["spill_time_nsec"] * 1e-9
    )

    hdr_reset = hdr.reset_index()[["__ntuple", "entry", "evt", "global_trigger_time"]]
    hdr_reset["trigger_time_s"] = hdr_reset["global_trigger_time"] * 1e-9
    hdr_reset["_orig_idx"] = np.arange(len(hdr_reset))

    pot_by_ntuple = {
        ntuple: grp.sort_values("spill_time").reset_index(drop=True)
        for ntuple, grp in pot_reset.groupby("__ntuple")
    }
    all_spills = pot_reset.sort_values("spill_time").reset_index(drop=True)

    results = []
    for ntuple, hdr_grp in hdr_reset.groupby("__ntuple"):
        spills = pot_by_ntuple.get(ntuple, all_spills)
        times = spills["spill_time"].values
        trigger_times = hdr_grp["trigger_time_s"].values
        idx = np.searchsorted(times, trigger_times, side="right") - 1
        idx = np.clip(idx, 0, len(spills) - 1)
        matched_spills = spills.iloc[idx][list(BEAM_COLS)].reset_index(drop=True)
        result = hdr_grp[["_orig_idx"]].reset_index(drop=True)
        result = pd.concat([result, matched_spills], axis=1)
        results.append(result)

    matched_ordered = pd.concat(results, ignore_index=True).sort_values("_orig_idx")
    for col in BEAM_COLS:
        hdr[col] = matched_ordered[col].values
    return hdr


def short_run_ids(
    data_hdr_df: pd.DataFrame,
    *,
    min_run_duration_min: float = 20.0,
) -> set:
    """Return run IDs whose triggered-event span is shorter than the threshold."""
    short_runs: set = set()
    for run, run_df in data_hdr_df.groupby("run", sort=False):
        t_min = run_df["global_trigger_time"].min()
        t_max = run_df["global_trigger_time"].max()
        duration_min = (t_max - t_min) / 60e9
        if duration_min < min_run_duration_min:
            short_runs.add(run)
    return short_runs


def bad_hdr_mask(
    data_hdr_df: pd.DataFrame,
    *,
    fom_cut: float = 0.98,
    short_runs: set | None = None,
) -> pd.Series:
    fom = data_hdr_df["FOM"]
    fom_failure = (fom == 0) | (fom == 4)
    fom_bad_value = (fom < fom_cut) & (fom > 0)
    bad_spill = fom_failure | fom_bad_value
    if short_runs is None:
        short_runs = short_run_ids(data_hdr_df)
    in_short_run = data_hdr_df["run"].isin(short_runs)
    return bad_spill | in_short_run


def apply_beam_quality_cuts(
    data_evt_df: pd.DataFrame,
    data_hdr_df: pd.DataFrame,
    pot_df: pd.DataFrame,
    *,
    fom_cut: float = 0.98,
    min_run_duration_min: float = 20.0,
) -> tuple[pd.DataFrame, pd.DataFrame, BeamQualitySummary]:
    """Return ``(evt_good, hdr_good, summary)`` after beam-quality selection."""
    hdr_matched = match_beam_spills(data_hdr_df, pot_df)
    short_runs = short_run_ids(
        hdr_matched, min_run_duration_min=min_run_duration_min
    )
    bad_mask = bad_hdr_mask(
        hdr_matched, fom_cut=fom_cut, short_runs=short_runs
    )

    bad_keys = set(hdr_matched.index[bad_mask])
    evt_prefilter = data_evt_df.copy()
    if not isinstance(evt_prefilter.index, pd.MultiIndex):
        raise ValueError("Expected evt MultiIndex with slice level")
    evt_hdr_index = evt_prefilter.index.droplevel(-1)
    evt_good = evt_prefilter.loc[~evt_hdr_index.isin(bad_keys)].copy()
    hdr_good = hdr_matched.loc[~bad_mask].copy()

    summary = BeamQualitySummary(
        n_evt_prefilter=len(evt_prefilter),
        n_evt_good=len(evt_good),
        n_hdr_total=len(hdr_matched),
        n_hdr_bad=int(bad_mask.sum()),
        n_short_runs=len(short_runs),
        pot_prefilter=float(hdr_matched["pot"].sum()),
        pot_good=float(hdr_good["pot"].sum()),
    )
    return evt_good, hdr_good, summary
