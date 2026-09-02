#!/usr/bin/env python3
"""MC/data overlay plots for final-selected cross-section variables.

Configurable input directories (``PLOT_SETS`` below). Used by
``notebooks/selected_xsec_overlay.ipynb`` and runnable headless for batch jobs.
"""

from __future__ import annotations

import gc
import glob
import os
import sys
import warnings
from datetime import datetime
from functools import partial
from os import makedirs, path

import numpy as np
import pandas as pd

sys.path.append("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")

from pyanalib.split_df_helpers_new import (  # noqa: E402
    _concat_hdf_frames,
    _remap_ntuple_index,
    _unique_ntuple_values_across_keys,
    dfs_from_dir,
    get_n_split,
    load_dfs,
)
from analysis_village.numucc_1p0pi.beam_quality import apply_beam_quality_cuts  # noqa: E402
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
)
from analysis_village.numucc_1p0pi import utils as numucc_utils  # noqa: E402
from analysis_village.numucc_1p0pi.utils import (  # noqa: E402
    get_pot_str,
    overlay_hists,
)
from analysis_village.numucc_1p0pi.utils import get_syst_unc as load_syst_disk_unc  # noqa: E402
from analysis_village.numucc_1p0pi.utils import _DEFAULT_SYST_DISK_ROOT  # noqa: E402
from analysis_village.numucc_1p0pi.syst_disk_layout import SYST_DISK_ENV  # noqa: E402

import matplotlib.pyplot as plt  # noqa: F401, E402

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
# ``utils`` import above already applies ``notebooks/presentation.mplstyle``.
_ = numucc_utils

# ---------------------------------------------------------------------------
# Configurable inputs
# ---------------------------------------------------------------------------
DFS_ROOT = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"
OUTPUT_BASE = "/exp/sbnd/data/users/munjung/PRL_data"

KEYS2LOAD = ("hdr", "evt")
DATA_KEYS2LOAD = ("hdr", "evt", "bnbpot", "trigger")
N_MAX_CONCAT = 999
MEMORY_LIMIT_FRAC = 0.60  # stop MC load and retry with half the files if exceeded

# Beam-quality selection (``notebooks/beam_quality.ipynb``)
APPLY_BEAM_QUALITY = True
FOM_CUT = 0.98
MIN_RUN_DURATION_MIN = 20.0
# When beam quality is off, optional scalar on full hdr POT (legacy approximation)
DATA_POT_SCALE = 1.0

SYST_DISK_ROOT = os.environ.get(SYST_DISK_ENV) or _DEFAULT_SYST_DISK_ROOT
LOAD_SYST = True

PLOT_SETS = [
    {
        "tag": "selected_chi2fix_qual_cut",
        "output_dir": path.join(OUTPUT_BASE, "selected_chi2fix_qual_cut"),
        "mc_dir": path.join(
            DFS_ROOT, "2026_09_01_063545__sel_mup-mc-fvfix-chi2fix-real"
        ),
        "mc_filename_str": "sel_mup-mc-fvfix-chi2fix-real",
        "data_dir": path.join(
            DFS_ROOT, "2026_09_01_063308__sel_mup-data-1e20-fvfix-chi2fix-real"
        ),
        "data_filename_str": "sel_mup-data-1e20-fvfix-chi2fix-real",
    },
    {
        "tag": "selected_fvfix_qual_cut",
        "output_dir": path.join(OUTPUT_BASE, "selected_fvfix_qual_cut"),
        "mc_dir": path.join(DFS_ROOT, "2026_09_01_063924__sel_mup-mc-fvfix"),
        "mc_filename_str": "sel_mup-mc-fvfix",
        "data_dir": path.join(
            DFS_ROOT, "2026_09_01_064250__sel_mup-data-1e20-fvfix"
        ),
        "data_filename_str": "sel_mup-data-1e20-fvfix",
    },
]

# Cross-section measurement variables (same list as unfolding.ipynb)
VAR_CONFIGS = [
    vc
    for vc in CORE_SELECTED_EVT_VARIABLE_CONFIGS
    if vc.var_save_name != "integrated"
]

BREAKDOWN_TYPES = ("topology", "genie_sb")

# Plot style (event_selection.ipynb / selected_events.py)
RATIO = True
AX_YLIM_RATIO = 1.9
TEXTLOC = [0.03, 0.55]
APPROVAL = "internal"
SAVE_FIG = True
PLOT = False
TEXTCHI2 = True


class MemoryLimitExceeded(RuntimeError):
    """Raised when system memory use exceeds ``MEMORY_LIMIT_FRAC``."""


def get_memory_used_frac() -> float:
    """Fraction of total RAM in use (1 - MemAvailable/MemTotal)."""
    with open("/proc/meminfo") as f:
        info = {}
        for line in f:
            key, rest = line.split(":", 1)
            info[key.strip()] = int(rest.split()[0])  # kB
    total = info["MemTotal"]
    avail = info.get("MemAvailable", info["MemFree"])
    return 1.0 - (avail / total)


def count_matching_files(sample_dir: str, filename_str: str) -> int:
    pattern = path.join(sample_dir, f"*{filename_str}*.df")
    return len(glob.glob(pattern))


def check_memory(label: str) -> None:
    frac = get_memory_used_frac()
    print(f"  [{label}] memory used: {frac * 100:.1f}%", flush=True)
    if frac >= MEMORY_LIMIT_FRAC:
        raise MemoryLimitExceeded(
            f"{label}: memory {frac * 100:.1f}% >= {MEMORY_LIMIT_FRAC * 100:.0f}% limit"
        )


def dfs_from_dir_monitored(
    search_dir: str,
    filename_str: str,
    *,
    n_max_concat: int,
    label: str = "load",
):
    """Like ``dfs_from_dir`` but checks memory after each file."""
    from tqdm import tqdm

    df_lists = {k: [] for k in KEYS2LOAD}
    ntuple_offset = np.int64(0)
    pattern = path.join(search_dir, f"*{filename_str}*.df")
    files_to_process = sorted(glob.glob(pattern))[: int(n_max_concat)]
    print(
        f"  [{label}] loading {len(files_to_process)} file(s) from {search_dir}",
        flush=True,
    )
    if not files_to_process:
        raise FileNotFoundError(f"No files matching {pattern!r}")

    check_memory(f"{label} before files")

    for mc_file in tqdm(files_to_process, desc=label):
        mc_n_split = get_n_split(mc_file)
        mc_dfs = load_dfs(mc_file, list(KEYS2LOAD), n_max_concat=int(mc_n_split))
        unique_ntuples = _unique_ntuple_values_across_keys(mc_dfs, KEYS2LOAD)
        ntuple_remap = {
            old: np.int64(ntuple_offset + i) for i, old in enumerate(unique_ntuples)
        }
        for df_key in KEYS2LOAD:
            df = mc_dfs[df_key]
            _remap_ntuple_index(df, ntuple_remap)
            df_lists[df_key].append(df)
        ntuple_offset += np.int64(len(ntuple_remap))
        check_memory(f"{label} after {path.basename(mc_file)}")

    concat_dfs = {
        k: _concat_hdf_frames(df_lists[k], label=k) for k in KEYS2LOAD if df_lists[k]
    }
    print("REMEMBER TO RECALCULATE TKI AND CHECK FV!!", flush=True)
    return concat_dfs


def _finalize_evt_hdr(dfs: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    evt = dfs["evt"]
    hdr = dfs["hdr"]
    if "mc" in evt.columns.get_level_values(0):
        evt.loc[evt.mc.iscc.isna(), ("mc", "iscc")] = 999
    return evt, hdr


def load_sample(
    sample_dir: str,
    filename_str: str,
    *,
    n_max_concat: int | None = None,
    label: str = "",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cap = N_MAX_CONCAT if n_max_concat is None else n_max_concat
    dfs = dfs_from_dir(
        sample_dir,
        filename_str=filename_str,
        keys2load=list(KEYS2LOAD),
        n_max_concat=cap,
    )
    return _finalize_evt_hdr(dfs)


def load_mc_sample(
    mc_dir: str,
    filename_str: str,
) -> tuple[pd.DataFrame, pd.DataFrame, int, int]:
    """Load MC with memory monitoring; retry with half the files if limit hit."""
    n_total = count_matching_files(mc_dir, filename_str)
    n_try = n_total
    min_files = max(1, n_total // 2)

    while True:
        gc.collect()
        check_memory("MC before load")
        try:
            dfs = dfs_from_dir_monitored(
                mc_dir,
                filename_str,
                n_max_concat=n_try,
                label=f"MC ({n_try}/{n_total} files)",
            )
            check_memory(f"MC after concat ({n_try}/{n_total} files)")
            evt, hdr = _finalize_evt_hdr(dfs)
            return evt, hdr, n_try, n_total
        except MemoryLimitExceeded as exc:
            next_try = max(min_files, n_try // 2)
            if n_try <= min_files or next_try >= n_try:
                raise RuntimeError(
                    f"MC load exceeded memory limit even with {n_try}/{n_total} files"
                ) from exc
            n_try = next_try
            print(
                f"  Memory limit hit — retrying MC with {n_try}/{n_total} files",
                flush=True,
            )
            gc.collect()


def _get_syst_cov(var_config):
    if not LOAD_SYST:
        return None
    try:
        frac_unc, cov = load_syst_disk_unc(
            var_config,
            syst_disk_root=SYST_DISK_ROOT,
            skip_missing_vars=True,
        )
        if cov is None:
            return None
        return cov
    except Exception as ex:
        print(f"  syst skip {var_config.var_save_name}: {ex}")
        return None


def load_data_sample(
    data_dir: str,
    filename_str: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load data evt/hdr, optionally applying beam-quality cuts."""
    keys = list(DATA_KEYS2LOAD if APPLY_BEAM_QUALITY else KEYS2LOAD)
    dfs = dfs_from_dir(
        data_dir,
        filename_str=filename_str,
        keys2load=keys,
        n_max_concat=N_MAX_CONCAT,
    )
    data_evt = dfs["evt"]
    data_hdr = dfs["hdr"]
    if APPLY_BEAM_QUALITY:
        data_hdr = data_hdr.join(dfs["trigger"])
        data_evt, data_hdr, summary = apply_beam_quality_cuts(
            data_evt,
            data_hdr,
            dfs["bnbpot"],
            fom_cut=FOM_CUT,
            min_run_duration_min=MIN_RUN_DURATION_MIN,
        )
        print(
            f"  beam quality: evt {summary.n_evt_good:,}/{summary.n_evt_prefilter:,} "
            f"({100 * summary.n_evt_good / summary.n_evt_prefilter:.3f}%)  "
            f"POT {summary.pot_good:.3e} / {summary.pot_prefilter:.3e}",
            flush=True,
        )
    data_evt[("mc", "iscc")] = 999
    return data_evt, data_hdr


def setup_pot_weights(
    mc_evt: pd.DataFrame,
    mc_hdr: pd.DataFrame,
    data_evt: pd.DataFrame,
    data_hdr: pd.DataFrame,
) -> str:
    if APPLY_BEAM_QUALITY:
        data_tot_pot = data_hdr["pot"].sum()
    else:
        data_tot_pot = data_hdr["pot"].sum() * DATA_POT_SCALE
    pot_str = get_pot_str(data_tot_pot)
    data_evt["pot_weight"] = np.ones(len(data_evt))

    mc_tot_pot = mc_hdr["pot"].sum()
    mc_scale = data_tot_pot / mc_tot_pot
    mc_evt["pot_weight"] = mc_scale * np.ones(len(mc_evt))
    print(
        f"  data POT={data_tot_pot:.3e}  MC POT={mc_tot_pot:.3e}  scale={mc_scale:.3e}"
    )
    print(f"  evt rows: data={len(data_evt):,}  mc={len(mc_evt):,}")
    return f"Events / Bin (POT={pot_str})"


def run_plot_set(plot_set: dict) -> None:
    tag = plot_set["tag"]
    out_dir = plot_set["output_dir"]
    print(f"\n{'=' * 72}\nPlot set: {tag}\n  MC:   {plot_set['mc_dir']}\n  data: {plot_set['data_dir']}")
    if SAVE_FIG:
        makedirs(out_dir, exist_ok=True)
        print(f"  saving -> {out_dir}")

    mc_evt, mc_hdr, n_mc_loaded, n_mc_total = load_mc_sample(
        plot_set["mc_dir"], plot_set["mc_filename_str"]
    )
    if n_mc_loaded < n_mc_total:
        print(
            f"  MC subsample: {n_mc_loaded}/{n_mc_total} files "
            f"(pot_weight scaled in setup_pot_weights)",
            flush=True,
        )
    data_evt, data_hdr = load_data_sample(
        plot_set["data_dir"], plot_set["data_filename_str"]
    )
    pot_label = setup_pot_weights(mc_evt, mc_hdr, data_evt, data_hdr)

    plotter = partial(
        overlay_hists,
        mc_df=mc_evt,
        data_df=data_evt,
        intime_df=None,
        dirt_df=None,
        ax_ylim_ratio=AX_YLIM_RATIO,
        ratio=RATIO,
        textloc=TEXTLOC,
        approval=APPROVAL,
        save_fig=SAVE_FIG,
        plot=PLOT,
    )

    for var_config in VAR_CONFIGS:
        syst = _get_syst_cov(var_config)
        for breakdown_type in BREAKDOWN_TYPES:
            save_name = path.join(
                out_dir, f"{var_config.var_save_name}_{breakdown_type}"
            )
            plot_labels = [var_config.var_labels[1], pot_label, ""]
            print(f"  plot {var_config.var_save_name} ({breakdown_type})")
            plotter(
                breakdown_type=breakdown_type,
                var_config=var_config,
                plot_labels=plot_labels,
                syst=syst,
                textchi2=TEXTCHI2,
                save_name=save_name,
            )


def main():
    t0 = datetime.now()
    print(f"selected_xsec_overlay start {t0.isoformat()}", flush=True)
    print(
        f"SYST_DISK_ROOT={SYST_DISK_ROOT}  LOAD_SYST={LOAD_SYST}  "
        f"APPLY_BEAM_QUALITY={APPLY_BEAM_QUALITY}  "
        f"MEMORY_LIMIT={MEMORY_LIMIT_FRAC * 100:.0f}%  "
        f"mem_now={get_memory_used_frac() * 100:.1f}%",
        flush=True,
    )
    for plot_set in PLOT_SETS:
        run_plot_set(plot_set)
    print(f"\nDone in {datetime.now() - t0}")


if __name__ == "__main__":
    main()
