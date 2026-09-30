import os

import numpy as np
import pandas as pd
from tqdm import tqdm
import string
import pickle
# Do NOT import statsmodels at module load time: grid worker venvs often lack it,
# and syst_histcounts / headless batch only need a few helpers from this module.
# (See selection_framework.get_clipped_evts note; get_eff_err imports lazily.)

import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from pyanalib.split_df_helpers import *
from pyanalib.stat_helpers import *
from pyanalib.covariance import *

from makedf.constants import *
from analysis_village.unfolding.wienersvd import *
from analysis_village.numucc_1p0pi.categories import *
from analysis_village.numucc_1p0pi.constants import *
from analysis_village.numucc_1p0pi.selection_framework import multicol_get_series
from analysis_village.numucc_1p0pi.syst_disk_layout import (
    FILE_GENIE,
    SUB_GENIE,
    SYST_DISK_ENV,
    category_out_dir,
    category_summary_npz_path,
    syst_disk_paths,
)

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
from matplotlib.patches import Patch
from matplotlib.legend import Legend
from matplotlib.legend_handler import HandlerErrorbar
from matplotlib.lines import Line2D
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), "notebooks", "presentation.mplstyle"))
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)

pdg_labels = [r"$\mu^{\pm}$", r"$p$", r"$\pi^{\pm}$", r"Other"]
pdg_colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
# Data-driven (offbeam/intime) cosmic layer stacked *in addition to* MC cosmics.
PDG_COSMIC_LABEL = "Intime Cosmics"
PDG_COSMIC_COLOR = "dimgray"

dpi = 300
fig_ext = ".png"


# ====== systematics disk: root resolution ======
def _fail_syst_disk(msg: str) -> None:
    """Emit a high-visibility error and abort (no silent fallbacks)."""
    banner = "=" * 72
    block = "\n%s\nFATAL [systematics disk]: %s\n%s\n" % (banner, msg, banner)
    print(block, file=sys.stderr)
    raise FileNotFoundError(msg)


def resolve_syst_disk_root(explicit): # str | None) -> str:
    """Return normalized syst disk root from argument or ``NUMUCC_SYST_DISK_ROOT``."""
    root = explicit or os.environ.get(SYST_DISK_ENV)
    if not root:
        _fail_syst_disk(
            "No syst disk root. Set environment variable %s or pass syst_disk_root= to "
            "get_syst_unc(). Expected layout: <root>/MCstat/mcstat_syst_dict.npz, "
            "<root>/Flux/flux_syst_dict.npz, … — see analysis_village.numucc_1p0pi.syst_disk_layout."
            % SYST_DISK_ENV
        )
    return syst_disk_paths(root)["root"]


# Keys accepted by ``get_syst_unc(..., syst_components=...)`` (case-insensitive strings).
SYST_UNC_DISK_KEYS = ("mcstat", "flux", "g4", "genie", "cosmics", "detector")
SYST_UNC_FLAT_KEYS = ("pot", "ntargets")
SYST_UNC_ALL_KEYS = SYST_UNC_DISK_KEYS + SYST_UNC_FLAT_KEYS
_SYST_UNC_DISK_LABELS = {
    "mcstat": "MC stat.",
    "genie": "GENIE",
    "flux": "Flux",
    "g4": "G4",
    "cosmics": "Cosmics",
    "detector": "Detector",
}


# ======= load systematic uncertainties from pre-saved files ======
def get_syst_unc(
    var_config,
    plot=False,
    save_fig=False,
    save_name=None,
    syst_disk_root=None,
    syst_components=None,
    genie_cov_frac_key: str = "genie",
    skip_missing_vars: bool = False,
):
    """Load fractional covariance blocks from the syst-disk tree and combine into total covariance.

    All inputs live under a single root directory (see ``syst_disk_layout``): ``MCstat/``,
    ``Flux/``, ``G4/``, ``GENIE/``, ``Cosmics/``, ``Detector/``. If ``syst_disk_root`` is omitted,
    ``NUMUCC_SYST_DISK_ROOT`` must be set. **Missing files abort with a loud error** — there are
    no alternate search paths or dated campaign fallbacks.

    Parameters
    ----------
    syst_components
        Optional subset of uncertainty sources to include. Each entry is a string, case-insensitive,
        chosen from disk-backed keys ``mcstat``, ``flux``, ``g4``, ``genie``, ``cosmics``,
        ``detector`` and flat correlated terms ``pot``, ``ntargets``. If ``None`` (default), all
        of the above are included (original behavior). Only files needed for the selected disk
        keys are required on disk.
    genie_cov_frac_key
        Which matrix to read from ``GENIE/cov_mat_dict.pkl`` for the ``genie`` disk component:
        ``"genie"`` (response / **xsec** path) or ``"genie_rate"`` (**rate** reweight path), matching
        :mod:`syst_genie_aggregate`. Default ``"genie"`` preserves legacy behavior.
    skip_missing_vars
        If ``True``, omit disk-backed components whose files lack ``var_config.var_save_name``,
        whose covariance shape does not match ``var_config`` bins, or that otherwise fail to
        combine for this variable, instead of raising. Useful for overlay plots when only a
        subset of variables has been produced on the syst disk.
    """
    if syst_components is None:
        active = frozenset(SYST_UNC_ALL_KEYS)
    else:
        active = frozenset(str(x).lower() for x in syst_components)
        unknown = active - frozenset(SYST_UNC_ALL_KEYS)
        if unknown:
            raise ValueError(
                "Invalid syst_components keys: %s. Allowed: %s"
                % (", ".join(sorted(unknown)), ", ".join(SYST_UNC_ALL_KEYS))
            )

    need_disk = any(k in active for k in SYST_UNC_DISK_KEYS)
    root = None
    paths = None
    if need_disk:
        root = resolve_syst_disk_root(syst_disk_root)
        paths = syst_disk_paths(root)
        missing = [
            (k, paths[k])
            for k in SYST_UNC_DISK_KEYS
            if k in active and not os.path.isfile(paths[k])
        ]
        if missing:
            detail = "\n".join("  [%s] %s" % (role, pth) for role, pth in missing)
            _fail_syst_disk(
                "Missing systematic covariance file(s). Run the producer pipelines into the "
                "expected locations, then retry:\n%s" % detail
            )

    def _load_disk_frac_cov(key: str) -> np.ndarray:
        assert paths is not None
        if key == "mcstat":
            blob = np.load(paths["mcstat"], allow_pickle=True)
            return dict(blob)[var_config.var_save_name].item()["MCstat"]["cov_frac"]
        if key == "flux":
            blob = np.load(paths["flux"], allow_pickle=True)
            return dict(blob)[var_config.var_save_name].item()["flux"]["cov_frac"]
        if key == "g4":
            blob = np.load(paths["g4"], allow_pickle=True)
            return dict(blob)[var_config.var_save_name].item()["G4"]["cov_frac"]
        if key == "genie":
            if genie_cov_frac_key not in ("genie", "genie_rate"):
                raise ValueError(
                    "genie_cov_frac_key must be 'genie' or 'genie_rate', got %r" % (genie_cov_frac_key,)
                )
            with open(paths["genie"], "rb") as gf:
                genie_blob = pickle.load(gf)
            row = genie_blob[var_config.var_save_name]
            if genie_cov_frac_key not in row:
                raise KeyError(
                    "GENIE pickle for %r has no %r (keys: %s)"
                    % (var_config.var_save_name, genie_cov_frac_key, sorted(row.keys()))
                )
            return row[genie_cov_frac_key]
        if key == "cosmics":
            blob = np.load(paths["cosmics"], allow_pickle=True)
            return dict(blob)[var_config.var_save_name].item()["Cosmics"]["cov_frac"]
        if key == "detector":
            blob = np.load(paths["detector"], allow_pickle=True)
            return dict(blob)["detector"].item()[var_config.var_save_name]["cov_frac"]
        raise KeyError(key)

    # flat uncertainties
    pot_frac_unc = 0.02
    ntargets_frac_unc = 0.01

    frac_uncert_total = np.zeros(len(var_config.bin_centers))
    frac_cov_matrix_total = np.zeros((len(var_config.bin_centers), len(var_config.bin_centers)))

    n_bins = len(var_config.bin_centers)
    for key in SYST_UNC_DISK_KEYS:
        if key not in active:
            continue
        syst_name = _SYST_UNC_DISK_LABELS[key]
        try:
            syst = _load_disk_frac_cov(key)
            syst_uncert = np.sqrt(np.maximum(np.diag(np.asarray(syst, dtype=float)), 0.0))
            if syst.shape != (n_bins, n_bins):
                raise ValueError(
                    "cov_frac shape %s does not match %d bins for %r"
                    % (syst.shape, n_bins, var_config.var_save_name)
                )
            frac_uncert_total += syst_uncert ** 2
            frac_cov_matrix_total += np.asarray(syst, dtype=float)
            if plot:
                plt.hist(
                    var_config.bin_centers,
                    bins=var_config.bins,
                    weights=syst_uncert,
                    histtype="step",
                    linewidth=2,
                    label=syst_name,
                )
        except (KeyError, ValueError) as ex:
            if not skip_missing_vars:
                raise
            print(
                "[get_syst_unc] skip %s for %r: %s"
                % (syst_name, var_config.var_save_name, ex),
                flush=True,
            )
            continue

    if "pot" in active:
        # Multiplicative exposure scale → fully correlated fractional cov.
        syst_name = "POT"
        u = pot_frac_unc
        syst_uncert = u * np.ones(n_bins)
        frac_uncert_total += syst_uncert ** 2
        frac_cov_matrix_total += np.full((n_bins, n_bins), u * u, dtype=float)
        if plot:
            plt.hist(var_config.bin_centers, bins=var_config.bins, weights=syst_uncert,   histtype="step", linewidth=2, label=syst_name)
    if "ntargets" in active:
        syst_name = "Ntargets"
        u = ntargets_frac_unc
        syst_uncert = u * np.ones(n_bins)
        frac_uncert_total += syst_uncert ** 2
        frac_cov_matrix_total += np.full((n_bins, n_bins), u * u, dtype=float)
        if plot:
            plt.hist(var_config.bin_centers, bins=var_config.bins, weights=syst_uncert,   histtype="step", linewidth=2, label=syst_name)

    frac_uncert_total = np.sqrt(frac_uncert_total)
    syst = frac_uncert_total

    if plot:
        plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_uncert_total,    histtype="step", linewidth=2, color="k",  label="Total")

        plt.xlim(var_config.bins[0], var_config.bins[-1])
        plt.ylim(0, max(frac_uncert_total) * 1.4)

        plt.xlabel(var_config.var_labels[1])
        plt.ylabel("Uncertainty [%]")
        plt.legend(fontsize=11, ncol=3, loc="upper center")

        plt.grid(which='major', linestyle='-', linewidth=0.7, alpha=0.7)
        plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.5)
        plt.minorticks_on()

        if save_fig:
            plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

        if not plot:
            plt.close()
        else:
            plt.show();

    return syst, frac_cov_matrix_total


_CATEGORY_SYST_SUMMARY_CACHE = {}
# Nominal Product B syst disk in the 2026-09-30 PRL tree (real files, formerly
# the dentsmooth consumer: GENIE_slim_v3 + MEC May, ``tki-del_Tp`` GENIE from
# v1×v3, DENT rolling 80% w=3 + Gauss σ=1).
_DEFAULT_SYST_DISK_ROOT = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/systematics/productB_sel_mup"
)


# ====== category-summary & GENIE-SB covariance loaders ======
def resolve_category_syst_summary_path(
    category_syst_summary_path=None,
    syst_disk_root=None,
):
    """Path to ``CategorySummary/category_syst_summary.npz`` from export cell."""
    if category_syst_summary_path:
        return os.path.abspath(os.path.expanduser(category_syst_summary_path))
    root = syst_disk_root or os.environ.get(SYST_DISK_ENV) or _DEFAULT_SYST_DISK_ROOT
    return category_summary_npz_path(root)


def load_overlay_syst_cov_frac(
    var_config,
    *,
    syst_kind="rate",
    syst_disk_root=None,
    category_syst_summary_path=None,
):
    """Fractional covariance for overlay bands (default: summed category summary)."""
    from analysis_village.numucc_1p0pi.syst_category_summary import (
        load_category_syst_summary,
        total_cov_frac,
    )

    path = resolve_category_syst_summary_path(
        category_syst_summary_path, syst_disk_root
    )
    vsn = (
        getattr(var_config, "category_syst_var_save_name", None)
        or var_config.var_save_name
    )
    cache_key = (path, vsn, syst_kind)
    if cache_key in _CATEGORY_SYST_SUMMARY_CACHE:
        print("key", cache_key, "in cache")
        return _CATEGORY_SYST_SUMMARY_CACHE[cache_key]
    if not os.path.isfile(path):
        raise FileNotFoundError(
            "Category syst summary not found: %s (run systematics-summary export cell)"
            % path
        )
    summary = load_category_syst_summary(path)
    cov = total_cov_frac(summary, vsn, kind=syst_kind)
    _CATEGORY_SYST_SUMMARY_CACHE[cache_key] = cov
    return cov


def get_category_summary_syst_unc(
    var_config,
    *,
    syst_kind="rate",
    syst_disk_root=None,
    category_syst_summary_path=None,
):
    """Fractional diagonal uncertainty and covariance from ``category_syst_summary.npz``.

    *syst_kind* ``\"rate\"`` uses ``total_rate`` (GENIE rate); ``\"xsec\"`` uses ``total_xsec``.
    """
    cov = load_overlay_syst_cov_frac(
        var_config,
        syst_kind=syst_kind,
        syst_disk_root=syst_disk_root,
        category_syst_summary_path=category_syst_summary_path,
    )
    unc = np.sqrt(np.diag(cov))
    return unc, cov


GENIE_SB_BKGD_RATE_KEY = "genie_bkgd_rate"
_GENIE_SB_COV_MAT_CACHE: dict = {}


def resolve_genie_sb_cov_mat_pkl(genie_sb_cov_mat_pkl=None):
    """Path to ``systematics-genie-SB`` ``GENIE/cov_mat_dict.pkl`` (``genie_bkgd_rate``)."""
    if genie_sb_cov_mat_pkl:
        return os.path.abspath(os.path.expanduser(genie_sb_cov_mat_pkl))
    from analysis_village.numucc_1p0pi.files_config import save_fig_base_dir

    sb_root = os.path.join(save_fig_base_dir, "systematics-notebook-genie-SB-integrated")
    return os.path.join(category_out_dir(sb_root, SUB_GENIE), FILE_GENIE)


def load_genie_sb_bkgd_rate_cov_frac(var_config, genie_sb_cov_mat_pkl=None):
    """Fractional covariance on background topology rate from ``systematics-genie-SB.ipynb``."""
    pkl_path = resolve_genie_sb_cov_mat_pkl(genie_sb_cov_mat_pkl)
    vsn = var_config.var_save_name
    cache_key = (pkl_path, vsn, GENIE_SB_BKGD_RATE_KEY)
    if cache_key in _GENIE_SB_COV_MAT_CACHE:
        return _GENIE_SB_COV_MAT_CACHE[cache_key]
    if not os.path.isfile(pkl_path):
        raise FileNotFoundError(
            "GENIE SB cov_mat_dict not found: %s (run systematics-genie-SB.ipynb)" % pkl_path
        )
    with open(pkl_path, "rb") as f:
        cov_mat_dict = pickle.load(f)
    if vsn not in cov_mat_dict:
        raise KeyError(
            "Variable %r not in %s (keys sample: %s)"
            % (vsn, pkl_path, ", ".join(sorted(cov_mat_dict.keys())[:8]))
        )
    row = cov_mat_dict[vsn]
    if GENIE_SB_BKGD_RATE_KEY not in row:
        raise KeyError(
            "%r missing in %s for %r (have: %s)"
            % (
                GENIE_SB_BKGD_RATE_KEY,
                pkl_path,
                vsn,
                ", ".join(sorted(row.keys())[:12]),
            )
        )
    cov = np.asarray(row[GENIE_SB_BKGD_RATE_KEY], dtype=np.float64)
    _GENIE_SB_COV_MAT_CACHE[cache_key] = cov
    return cov


# ====== overlay-plot helpers: syst bands, legend fractions, chi2 ======
def _df_has_event_mc_truth(df) -> bool:
    """True when ``df`` has an event-level ``mc`` block (needed for topology/signal)."""
    if df is None or len(getattr(df, "columns", [])) == 0:
        return False
    cols = df.columns
    if isinstance(cols, pd.MultiIndex):
        return "mc" in cols.get_level_values(0)
    return "mc" in cols


def _overlay_signal_mc_hist(mc_df, var_config, signal_truth_fv="per_tpc"):
    """Per-bin MC signal (CC 1p0pi in FV) counts for background subtraction.

    Requires an event-level dataframe with ``mc.*`` truth columns. Track-level
    (``breakdown_type='pdg'``) frames must not call this.
    """
    if not _df_has_event_mc_truth(mc_df):
        raise ValueError(
            "_overlay_signal_mc_hist requires an event-level MC dataframe with an "
            "'mc' column block; got a frame without event truth (e.g. track-level pdg plot)."
        )
    vardf, _ = get_clipped_evts(mc_df, var_config.var_evt_reco_col, var_config.bins)
    cuts = get_topo_category(mc_df, ret_cuts=True, signal_truth_fv=signal_truth_fv)
    cut_signal = cuts[-1]
    v_sig = vardf[cut_signal]
    w_sig = mc_df.loc[cut_signal, "pot_weight"]
    hist_sig, _ = np.histogram(v_sig, weights=w_sig, bins=var_config.bins)
    return np.asarray(hist_sig, dtype=float)


def _overlay_bkgd_syst_sigma(total_mc_bkgd, bkgd_frac_cov):
    """Per-bin 1σ GENIE background-rate uncertainty (absolute event units)."""
    total_mc_bkgd = np.asarray(total_mc_bkgd, dtype=float)
    frac_diag = np.maximum(np.diag(np.asarray(bkgd_frac_cov, dtype=float)), 0.0)
    with np.errstate(invalid="ignore"):
        return np.sqrt(frac_diag) * total_mc_bkgd


def _overlay_draw_bkgd_syst_band(
    ax,
    bin_centers,
    bins,
    total_mc,
    bkgd_syst_err,
    *,
    edgecolor="darkorange",
    hatch="+++",
    label="Bkgd. GENIE unc.",
    zorder=9,
):
    bkgd_syst_err = np.asarray(bkgd_syst_err, dtype=float)
    ax.bar(
        bin_centers,
        2 * bkgd_syst_err,
        width=np.diff(bins),
        bottom=np.asarray(total_mc, dtype=float) - bkgd_syst_err,
        facecolor="none",
        hatch=hatch,
        linewidth=0.0,
        edgecolor=edgecolor,
        label=label,
        zorder=zorder,
    )


def _overlay_add_poisson_mc_stat_to_band(syst_explicit, load_syst_from_summary):
    """``category_syst_summary`` totals already include MC stat.; skip Poisson MC stat on the band."""
    return not (load_syst_from_summary and not syst_explicit)


def _overlay_syst_sigma(total_mc, mc_stat_err, syst_frac_cov, *, add_poisson_mc_stat):
    """Per-bin 1σ systematic uncertainty for hatched MC bands (absolute event units)."""
    total_mc = np.asarray(total_mc, dtype=float)
    syst_err_frac = np.sqrt(np.maximum(np.diag(np.asarray(syst_frac_cov, dtype=float)), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        syst_sigma = syst_err_frac * total_mc
        if add_poisson_mc_stat:
            mc_stat_err_frac = np.where(total_mc != 0, np.asarray(mc_stat_err, dtype=float) / total_mc, 0.0)
            syst_sigma = np.sqrt(syst_sigma ** 2 + (mc_stat_err_frac * total_mc) ** 2)
    return syst_sigma


def _overlay_chi2_valid_bins(total_data, total_mc, *, drop_first_bin: bool = False):
    """Bins with MC or data content (skip empty bins in χ²).

    ``drop_first_bin``: exclude bin 0 (e.g. nu_score Pandora failure bin).
    """
    total_data = np.asarray(total_data, dtype=float)
    total_mc = np.asarray(total_mc, dtype=float)
    valid = (total_mc > 0) | (total_data > 0)
    if drop_first_bin and valid.size:
        valid = np.asarray(valid, dtype=bool).copy()
        valid[0] = False
    return valid


def _overlay_data_xlim_range(bins, total_data, *, drop_first_bin: bool = False):
    """``(xmin, xmax, i_lo, i_hi)`` from bins with data > 0, or ``None``.

    Used to clip the x-axis to the support of the data. With ``drop_first_bin``,
    bin 0 is ignored (nu_score failure / dummy bin).
    """
    bins = np.asarray(bins, dtype=float)
    data = np.asarray(total_data, dtype=float)
    if data.size == 0 or bins.size < 2:
        return None
    mask = data > 0
    if drop_first_bin and mask.size:
        mask = mask.copy()
        mask[0] = False
    if not np.any(mask):
        return None
    idx = np.flatnonzero(mask)
    i_lo, i_hi = int(idx[0]), int(idx[-1])
    return float(bins[i_lo]), float(bins[i_hi + 1]), i_lo, i_hi


def _overlay_is_nu_score(var_config=None, histdata=None) -> bool:
    vsn = None
    if histdata is not None:
        vsn = getattr(histdata, "var_save_name", None)
    if not vsn and var_config is not None:
        vsn = getattr(var_config, "var_save_name", None)
    return str(vsn or "") == "nu_score"


def _overlay_data_stat_variance(data_eylow, data_eyhigh) -> np.ndarray:
    """Symmetrized data-stat variance from asymmetric error bars."""
    return (
        0.5
        * (
            np.asarray(data_eylow, dtype=float)
            + np.asarray(data_eyhigh, dtype=float)
        )
    ) ** 2


def _overlay_absolute_cov(
    total_mc,
    syst_frac_cov,
    data_eylow,
    data_eyhigh,
    *,
    mc_stat_err=None,
) -> np.ndarray:
    """Absolute covariance for overlay χ² / pulls.

    ``C_ij = Cfrac_ij * m_i * m_j`` plus diagonal data (and optional MC) stat.
    """
    from pyanalib.covariance import cov_from_fraccov

    m = np.asarray(total_mc, dtype=float)
    C = cov_from_fraccov(syst_frac_cov, m)
    diag = np.diag(C).copy() + _overlay_data_stat_variance(data_eylow, data_eyhigh)
    if mc_stat_err is not None:
        diag = diag + np.asarray(mc_stat_err, dtype=float) ** 2
    np.fill_diagonal(C, diag)
    C = 0.5 * (C + C.T)
    return np.nan_to_num(C, nan=0.0, posinf=0.0, neginf=0.0)


def _overlay_compute_chi2(
    total_data,
    total_mc,
    syst_frac_cov,
    data_eylow,
    data_eyhigh,
    *,
    mc_stat_err=None,
    drop_first_bin: bool = False,
):
    """χ² with the **full** absolute covariance (syst + data stat [+ MC stat]).

    Absolute syst block: ``C_ij = Cfrac_ij * m_i * m_j`` (same frac matrix as the
    hatched band). Data (and optional MC Poisson) variances are added on the
    diagonal. Off-diagonal syst correlations are kept.

    Also computes a shape-only χ² (MC normalized to the data integral): the syst
    block is rescaled by ``scale**2`` and its normalization component projected
    out, ``P C P^T`` with ``P = I - N 1^T / sum(N)``, before adding stat diagonals.

    ``drop_first_bin``: exclude bin 0 from the χ² (nu_score failure bin).
    """
    total_data = np.asarray(total_data, dtype=float)
    total_mc = np.asarray(total_mc, dtype=float)
    valid = _overlay_chi2_valid_bins(
        total_data, total_mc, drop_first_bin=drop_first_bin
    )
    if not np.any(valid):
        return None, None, None, None, None, None, None, None

    syst_frac_cov = np.nan_to_num(
        np.asarray(syst_frac_cov, dtype=float), nan=0.0, posinf=0.0, neginf=0.0
    )
    C = _overlay_absolute_cov(
        total_mc,
        syst_frac_cov,
        data_eylow,
        data_eyhigh,
        mc_stat_err=mc_stat_err,
    )
    keep = valid & (np.diag(C) > 0)
    idx = np.flatnonzero(keep)
    ndof = int(idx.size)
    if ndof == 0:
        return None, None, None, None, None, None, None, None

    d = total_data[idx]
    m = total_mc[idx]
    C_sub = np.asarray(C[np.ix_(idx, idx)], dtype=float)
    chi2_total, p_val = get_chi2(d, m, C_sub)
    chi2_reduced = chi2_total / ndof if ndof > 0 else None

    chi2_shape = None
    p_val_shape = None
    ndof_shape = None
    if m.sum() > 0 and ndof > 1:
        scale = float(d.sum() / m.sum())
        # Syst absolute cov scales as model²; data/MC-stat diagonals do not.
        C_syst = np.asarray(cov_from_fraccov(syst_frac_cov, total_mc), dtype=float)
        Cs = (scale**2) * C_syst[np.ix_(idx, idx)]
        n_norm = scale * m
        P = np.eye(ndof) - np.outer(n_norm, np.ones(ndof)) / n_norm.sum()
        Cs = P @ Cs @ P.T
        data_var = np.asarray(_overlay_data_stat_variance(data_eylow, data_eyhigh), dtype=float)
        diag_s = np.diag(Cs).copy() + data_var[idx]
        if mc_stat_err is not None:
            # MC-stat on the *normalized* prediction ≈ scale² × raw MC-stat var
            diag_s = diag_s + (scale**2) * np.asarray(mc_stat_err, dtype=float)[idx] ** 2
        np.fill_diagonal(Cs, diag_s)
        Cs = np.nan_to_num(0.5 * (Cs + Cs.T), nan=0.0, posinf=0.0, neginf=0.0)
        chi2_shape, p_val_shape = get_chi2_shape(d, m, Cs)
        ndof_shape = ndof - 1

    chi2_pull = np.full_like(total_data, np.nan, dtype=float)
    # Diagonal pulls (full-cov Mahalanobis pulls are not shown on the overlay).
    sig = np.sqrt(np.maximum(np.diag(C), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2_pull[idx] = (d - m) / np.maximum(sig[idx], 1e-10)
    return (
        chi2_total,
        chi2_reduced,
        p_val,
        ndof,
        chi2_pull,
        chi2_shape,
        p_val_shape,
        ndof_shape,
    )


def _resolve_overlay_syst_cov_frac(
    var_config,
    syst,
    *,
    syst_kind="rate",
    syst_disk_root=None,
    category_syst_summary_path=None,
    load_syst_from_summary=True,
):
    """Use explicit *syst* or load from ``category_syst_summary.npz`` when enabled."""
    if syst is not None:
        return np.asarray(syst, dtype=np.float64)
    if not load_syst_from_summary or var_config is None:
        return None
    try:
        return load_overlay_syst_cov_frac(
            var_config,
            syst_kind=syst_kind,
            syst_disk_root=syst_disk_root,
            category_syst_summary_path=category_syst_summary_path,
        )
    except Exception as ex:
        print("overlay_hists: could not load category_syst_summary (%s)" % ex)
        return None


# ====== general helpers: formatting, arrays, event clipping ======
def get_pot_str(tot_pot):
    """Format POT with two significant digits: ``$X.Y\\times 10^{YY}$``."""
    tot_pot = float(tot_pot)
    if not np.isfinite(tot_pot) or tot_pot <= 0:
        return r"$0$"
    exp_i = int(np.floor(np.log10(tot_pot)))
    mant = tot_pot / (10.0 ** exp_i)
    mant_r = round(mant, 1)  # two sig digs for mantissa in [1, 10)
    if mant_r >= 10.0:
        mant_r = 1.0
        exp_i += 1
    return rf"${mant_r:.1f}\times 10^{{{exp_i}}}$"


def format_pot_corner_text(pot) -> str:
    """Corner label ``SBND BNB $X.Y\\times 10^{Y}$ POT`` (two significant digits)."""
    import re

    prefix = r"$\mathbf{SBND\ BNB}$ "
    if pot is None or pot == "":
        return ""
    if isinstance(pot, (int, float)):
        return f"{prefix}{get_pot_str(float(pot))} POT"
    s = str(pot).strip()
    # Efficiency / fake-data plots stamp "SBND Simulation", not a POT value.
    if re.search(r"Simulation", s, flags=re.IGNORECASE):
        return s
    s = re.sub(r"\$\\mathbf\{SBND\\?\s*BNB\}\$\s*", "", s)
    s = re.sub(r"SBND\s+BNB\s*", "", s, flags=re.IGNORECASE)
    m = re.search(r"POT=\s*([^)]*)", s)
    if m:
        s = m.group(1).strip()
    s = re.sub(r"\s*POT\s*$", "", s).strip()
    if not s:
        return ""
    # Prefer parsing a scientific-notation value so we can reformat to 2 sig digs.
    m2 = re.search(
        r"([0-9]+\.?[0-9]*)\s*(?:\\times|×|x)\s*10\^?\{(-?\d+)\}",
        s.replace("$", ""),
    )
    if m2:
        val = float(m2.group(1)) * (10.0 ** int(m2.group(2)))
        return f"{prefix}{get_pot_str(val)} POT"
    try:
        return f"{prefix}{get_pot_str(float(s))} POT"
    except (TypeError, ValueError):
        pass
    inner = s.replace("$", "").strip()
    return f"{prefix}${inner}$ POT"


def generate_tags(end_tag=""):
    tags = []
    for first in string.ascii_lowercase:
        for second in string.ascii_lowercase:
            tag = first + second
            if tag == end_tag:
                break
            tags.append(tag)
        if tag == end_tag:
            break
    return tags


def _as_1d_float_array(values, name="values"):
    """Coerce event-level values to a 1D float ndarray (handles duplicate-column DataFrames)."""
    if isinstance(values, pd.DataFrame):
        if values.shape[1] == 1:
            values = values.iloc[:, 0]
        else:
            raise ValueError(
                "%s matched %d columns; expected a single event-level series"
                % (name, values.shape[1])
            )
    arr = np.ravel(np.asarray(values, dtype=float))
    return arr


def _var_weights_for_cut(var, weights, cut):
    """Return aligned 1D reco-variable and POT-weight arrays for one boolean category mask."""
    mask = np.asarray(cut, dtype=bool)
    v = _as_1d_float_array(var, name="variable")
    w = _as_1d_float_array(weights, name="pot_weight")
    n = len(mask)
    if len(v) != n or len(w) != n:
        raise ValueError(
            "Event array length mismatch: variable=%d pot_weight=%d mask=%d"
            % (len(v), len(w), n)
        )
    return v[mask], w[mask]


def get_clipped_evts(df, var_col, bins, verbose=False, var_save_name=None):
    # VariableConfig tuples are often padded to evt depth (e.g. 7); mcnu HDF may be 4-level.
    from analysis_village.numucc_1p0pi.variable_configs import (
        INTEGRATED_HIST_DUMMY,
        INTEGRATED_VAR_SAVE_NAME,
    )

    if var_save_name == INTEGRATED_VAR_SAVE_NAME:
        var = np.full(len(df), INTEGRATED_HIST_DUMMY, dtype=float)
    elif isinstance(var_col, tuple) and isinstance(df.columns, pd.MultiIndex):
        var = multicol_get_series(df, var_col)
    else:
        var = df[var_col]
    var = _as_1d_float_array(var, name=str(var_col))
    var = np.clip(var, bins[0], bins[-1] - EPSILON)

    if 'pot_weight' in df.columns:
        weights = df.loc[:, 'pot_weight']
    else:
        if verbose:
            print("No pot_weight column found, return 1 as pot scale (expected for data)")
        weights = np.ones_like(var)
    weights = _as_1d_float_array(weights, name="pot_weight")
    # One NaN weight makes numpy.histogram return NaN in all bins.
    weights = np.nan_to_num(weights, nan=0.0, posinf=0.0, neginf=0.0)
    return var, weights

def get_eff_err(success,total):  # success/total
    import statsmodels.api as sm

    err = [[],[]]
    eff = success/total
    for i in range(len(success)):
        this_success = success[i]
        this_tot = total[i]
        interval = sm.stats.proportion_confint(this_success,this_tot,method='wilson')
        err[0].append(abs(eff[i]-interval[0]))
        err[1].append(abs(eff[i]-interval[1]))
    return err


# ====== GENIE universe-weight access ======
def _multicol_first_nonempty_leaf(col) -> str:
    """First non-empty segment of a possibly padded MultiIndex column tuple."""
    if not isinstance(col, tuple):
        return str(col)
    for x in col:
        if x != "" and x is not None:
            return str(x)
    return str(col[0])


def genie_univ_weight_series(weight_block: pd.DataFrame, uidx: int) -> pd.Series:
    """Per-knob GENIE weights: multisim ``univ_*``, or unisim leaves aliased to ``univ_*``.

    Slim HDF stores multisigma/morph as ``ps*`` / ``ms*`` / ``morph`` under each knob name
    (not in bundled ``mc.GENIE``). Use :func:`get_systematics_genie.normalize_and_infer_n_univ`
    before calling :func:`get_univ_rates` on those blocks.
    """
    want = "univ_%d" % uidx
    cols = weight_block.columns
    if isinstance(cols, pd.MultiIndex):
        for c in cols:
            if _multicol_first_nonempty_leaf(c) == want:
                return weight_block[c]
        if uidx == 0:
            for alt in ("morph", "ps1"):
                for c in cols:
                    if _multicol_first_nonempty_leaf(c) == alt:
                        return weight_block[c]
    else:
        if want in weight_block.columns:
            return weight_block[want]
        if uidx == 0:
            for alt in ("morph", "ps1"):
                if alt in weight_block.columns:
                    return weight_block[alt]
    raise KeyError(
        "GENIE weight column %r not found under knob block (also tried morph, ps1 for uidx=0); "
        "columns=%s" % (want, list(cols)[:20])
    )


def _genie_weight_series(syst_type: str, weight_block: pd.DataFrame, uidx: int) -> pd.Series:
    univ_col = "univ_%d" % uidx
    if syst_type != "GENIE":
        return weight_block[univ_col]
    try:
        return weight_block[univ_col]
    except KeyError:
        return genie_univ_weight_series(weight_block, uidx)


# ====== response matrix & multisim universe rates ======
# NOTE: get_univ_rates is a *producer* helper used by scripts that write the
# pre-saved GENIE/Flux/G4 covariance files consumed by get_syst_unc. Analysis
# notebooks and overlay plots must not recompute covariances at runtime — load
# them with get_syst_unc / get_category_summary_syst_unc instead.
def smearing_reco_over_truth(reco_vs_true):
    """Per-bin ratio of reco- to truth-projected smearing matrix."""
    smear = np.asarray(reco_vs_true, dtype=float)
    truth = smear.sum(axis=1)
    reco = smear.sum(axis=0)
    return np.divide(
        reco, truth,
        out=np.zeros_like(reco, dtype=float),
        where=truth != 0,
    )


def _xsec_efficiency_per_bin(ret, w_evt_univ, w_nu_univ, bins):
    """Truth-bin efficiency entering the GENIE xsec response matrix.

    Matches ``get_systematics_genie.accumulate_xsec_path_chunk``: weighted
    selected-truth yield over weighted all-MC-truth yield, with universe
    weights applied positionally to the same event arrays as ``signal_hists``.
    """
    w_evt = np.nan_to_num(
        _as_1d_float_array(w_evt_univ, name="univ_evt_weight"),
        nan=1.0, posinf=1.0, neginf=1.0,
    )
    w_nu = np.nan_to_num(
        _as_1d_float_array(w_nu_univ, name="univ_nu_weight"),
        nan=1.0, posinf=1.0, neginf=1.0,
    )
    wgt_allmc = _as_1d_float_array(ret["wgt_allmc"], name="wgt_allmc")
    wgt_sel = _as_1d_float_array(ret["wgt_sel_truth"], name="wgt_sel_truth")
    signal_allmc, _ = np.histogram(
        ret["var_allmc"],
        weights=wgt_allmc * w_nu,
        bins=bins,
    )
    signal_sel, _ = np.histogram(
        ret["var_sel_truth"],
        weights=wgt_sel * w_evt,
        bins=bins,
    )
    return np.divide(
        signal_sel,
        signal_allmc,
        out=np.zeros_like(signal_sel, dtype=float),
        where=signal_allmc != 0,
    )


def get_univ_rates(cov_type="rate", 
                    syst_type="GENIE",
                    evtdf=None, 
                    nudf=None, 
                    var_config=None, 
                    syst_name="", 
                    n_univ=100, 
                    bkgd_subtract=True,
                    return_bkgd=False,
                    return_response=False,
                    xsec_unit=0,
                    plot=False,
                    verbose=False):
    """
    for the GENIE uncertainty on the xsec measurement
    """
    if cov_type == "xsec":
        if verbose:
            print("getting {} universes for {} uncertainty on the xsec".format(n_univ, syst_name))
        scale_factor = xsec_unit
        if xsec_unit == 0 and verbose:
            print("pass xsec_unit as an argument to get_univ_rates")

    elif cov_type == "rate":
        if verbose:
            print("getting {} universes for {} uncertainty on the event rate".format(n_univ, syst_name))
        scale_factor = 1.0

    else:
        raise ValueError("Invalid covariance type: {}, choose in [xsec, rate]".format(cov_type))

    bins = var_config.bins

    evtdf_signal = evtdf[evtdf.topo_categ == 1]
    # reco variable histogram, topology breakdown
    evtdf_div_topo = [evtdf[evtdf.topo_categ == mode]for mode in topology_list]

    if nudf is not None:
        # print("NUDDF", nudf.head())
        # for col in nudf.columns:
        #     print(col)
        nudf_signal = nudf[nudf.topo_categ == 1]

    ret = signal_hists(evtdf, nudf, var_config, return_data=True, plot=plot)

    univ_events = []
    univ_effs   = []
    univ_smears = []

    univ_events_bkgd = np.zeros((n_univ, len(bins)-1))
    cv_events_bkgd = np.zeros((len(bins)-1))

    for uidx in tqdm(range(n_univ), desc="Getting universes", disable=not verbose):
        univ_col = f"univ_{uidx}"
        w_evt_univ = _genie_weight_series(syst_type, evtdf_signal[syst_name], uidx)
        # cap weight at 10
        if syst_name == ("mc", "GENIE"):
            w_evt_univ *= 1
        w_evt_univ = np.clip(w_evt_univ, 0, 20)
        w_nu_univ = (
            _genie_weight_series(syst_type, nudf_signal[syst_name], uidx)
            if nudf is not None
            else None
        )

        # ---- signal channel ----
        # Integrated total cross section: use the rate path (selected yield), not
        # eff × N_gen^CV (which cancels normalization-like weights).
        use_xsec_response = (
            cov_type == "xsec"
            and syst_type == "GENIE"
            # and var_config.var_save_name != "integrated"
        )
        if use_xsec_response:
            # smearing matrix
            # handle case where there's a single bin, in which case there's no smearing
            if len(bins) == 2:
                reco_vs_true = np.array([[1.0]])
            else:
                w_evt = np.nan_to_num(
                    _as_1d_float_array(w_evt_univ, name="univ_evt_weight"),
                    nan=1.0, posinf=1.0, neginf=1.0,
                )
                wgt_sel = _as_1d_float_array(ret["wgt_sel_truth"], name="wgt_sel_truth")
                reco_vs_true, _, _ = np.histogram2d(
                    ret["var_sel_truth"],
                    ret["var_sel_reco"],
                    weights=wgt_sel * w_evt,
                    bins=bins,
                )
            univ_smears.append(reco_vs_true)

            eff = _xsec_efficiency_per_bin(ret, w_evt_univ, w_nu_univ, bins)
            univ_effs.append(eff)

            # print(signal_allmc_univ)
            response_univ = get_response_matrix(reco_vs_true, eff)
            signal_univ = response_univ @ ret["nevts_allmc"] # note that we multiply the CV signal rate!
            # signal_univ = signal_cv

        elif cov_type in ("rate", "xsec"):
            signal_univ, _ = np.histogram(
                ret["var_sel_reco"],
                weights=ret["wgt_sel_reco"] * w_evt_univ,
                bins=bins,
            )

        else:
            raise ValueError("Invalid covariance type: {}, choose xsec or rate".format(cov_type))

        # TODO: this isn't computationally efficient, but it's useful for debugging
        # ---- uncertainty on the background rate ----
        # loop over background categories
        # + univ background - cv background
        # note: cv background subtraction cancels out with the cv background subtraction for the cv event rate. 
        #       doing it anyways for the plot of universes on background subtracted event rate.
        for this_evtdf in evtdf_div_topo[1:]:
            var, wgt = get_clipped_evts(
                this_evtdf,
                var_config.var_evt_reco_col,
                bins,
                var_save_name=var_config.var_save_name,
            )
            univ_wgt = _genie_weight_series(syst_type, this_evtdf[syst_name], uidx).copy()
            if syst_name == ("mc", "GENIE"):
                univ_wgt *= 1
            univ_wgt = np.clip(univ_wgt, 0, 20)
            univ_wgt[np.isnan(univ_wgt)] = 1 ## IMPORTANT: make nan univ_wgt to 1. to ignore them
            background_cv, _   = np.histogram(var, bins=bins, weights=wgt)
            background_univ, _ = np.histogram(var, bins=bins, weights=wgt*univ_wgt)
            univ_events_bkgd[uidx] += background_univ
            # only add background cv for the first universe
            if uidx == 0:
                cv_events_bkgd += background_cv

            if bkgd_subtract:
                signal_univ += (background_univ - background_cv)
            else:
                signal_univ += background_univ

        signal_univ *= scale_factor
        univ_events.append(signal_univ)

    univ_events = np.array(univ_events)

    if bkgd_subtract:
        cv_events = ret["nevts_sel_reco"]
        cv_events *= scale_factor
    else:
        cv_events = ret["nevts_allsel_reco"]
        cv_events *= scale_factor 

    response_pack = None
    if return_response:
        if not (cov_type == "xsec" and syst_type == "GENIE"):
            raise ValueError("return_response requires cov_type='xsec' and syst_type='GENIE'")
        if not univ_effs:
            raise ValueError("return_response: no xsec efficiency universes were built")
        univ_effs_arr = np.asarray(univ_effs, dtype=float)
        if len(bins) == 2:
            cv_smears = np.array([[1.0]])
        else:
            wgt_sel = _as_1d_float_array(ret["wgt_sel_truth"], name="wgt_sel_truth")
            cv_smears, _, _ = np.histogram2d(
                ret["var_sel_truth"],
                ret["var_sel_reco"],
                weights=wgt_sel,
                bins=bins,
            )
        wgt_allmc = _as_1d_float_array(ret["wgt_allmc"], name="wgt_allmc")
        wgt_sel = _as_1d_float_array(ret["wgt_sel_truth"], name="wgt_sel_truth")
        signal_allmc_cv, _ = np.histogram(
            ret["var_allmc"],
            weights=wgt_allmc,
            bins=bins,
        )
        signal_sel_cv, _ = np.histogram(
            ret["var_sel_truth"],
            weights=wgt_sel,
            bins=bins,
        )
        cv_effs = np.divide(
            signal_sel_cv,
            signal_allmc_cv,
            out=np.zeros_like(signal_sel_cv, dtype=float),
            where=signal_allmc_cv != 0,
        )
        response_pack = {
            "univ_effs": univ_effs_arr,
            "cv_effs": cv_effs,
            "univ_smears": univ_smears,
            "cv_smears": cv_smears,
        }

    if return_bkgd:
        # sum over all background categories
        univ_events_bkgd = np.array(univ_events_bkgd) #.sum(axis=0)
        # univ_events_bkgd *= scale_factor
        cv_events_bkgd = np.array(cv_events_bkgd) #.sum(axis=0)
        # cv_events_bkgd *= scale_factor
        if return_response:
            return univ_events, cv_events, univ_events_bkgd, cv_events_bkgd, response_pack
        return univ_events, cv_events, univ_events_bkgd, cv_events_bkgd

    if return_response:
        return univ_events, cv_events, response_pack
    return univ_events, cv_events



def get_response_matrix(reco_vs_true, 
                        eff):
    """Truth → reco response ``R[j,i]`` with shape (reco bins, truth bins).

    ``reco_vs_true`` must be ``numpy.histogram2d(truth, reco, ...)[0]`` so that
    ``reco_vs_true[i,j]`` counts (truth bin ``i``, reco bin ``j``).

    With ``eff[i] = (selected signal in truth bin i) / (all generated signal in truth bin i)``,
    each column satisfies ``sum_j R[j,i] = eff[i]`` (not 1): migration is normalized within
    selected signal, then scaled by per-bin efficiency vs full true MC.
    """
    denom = reco_vs_true.T.sum(axis=0)
    num = reco_vs_true.T
    response = np.divide(
        num * eff, denom,
        out=np.zeros_like(num, dtype=float),  # fill with 0 where invalid
        where=denom != 0
    )
    return response



# ====== plotting functions ======

# ==== plot additions ====
def get_textloc_x(values, bins, textloc=[0.05, 0.55]):
    textloc_x, _ = textloc
    n_firsthalf = np.sum(values[:len(bins)//2])
    n_secondhalf = np.sum(values[len(bins)//2:])
    if n_firsthalf < n_secondhalf:
        textloc_x, textloc_ha = textloc_x, 'left'
    else:
        textloc_x, textloc_ha = 1-textloc_x, 'right'
    return textloc_x, textloc_ha


def _overlay_chi2_axes_loc(
    height_profile,
    ylim_max,
    *,
    breakdown_type: str = "",
    textloc=(0.05, 0.55),
    prefer_left: bool | None = None,
    reserved=None,
    block_h: float = 0.20,
):
    """Axes fraction (x, y, ha) for χ² text: in the clear band above all content.

    Y is placed above the *global* stack/syst peak (the headroom under the legend).
    X prefers the quieter half of the plot; PDG keeps χ² off the upper-left legend.
    ``reserved`` is a list of ``(x0, y0, x1, y1)`` axes-fraction boxes (legend /
    insets) that the χ² + GENIE block must not overlap.
    """
    h = np.asarray(height_profile, dtype=float)
    n = h.size
    if n == 0 or not np.isfinite(ylim_max) or ylim_max <= 0:
        return float(textloc[0]), 0.88, "left"

    mid = max(n // 2, 1)
    left_max = float(np.nanmax(h[:mid])) if mid else 0.0
    right_max = float(np.nanmax(h[mid:])) if mid < n else 0.0
    global_max = float(np.nanmax(h)) if n else 0.0

    if prefer_left is None:
        use_left = left_max <= right_max
        content_ref = global_max
    else:
        use_left = bool(prefer_left)
        # Explicit side (topology/genie): sit above the *local* stack, so a
        # 1-col legend on the empty side can keep χ² tucked under it.
        content_ref = left_max if use_left else right_max

    content_frac = max(0.0, min(1.0, content_ref / float(ylim_max)))

    # PDG legend is upper-left — keep χ² on the right.
    if breakdown_type == "pdg":
        use_left = False
        content_frac = max(0.0, min(1.0, global_max / float(ylim_max)))

    x0 = float(textloc[0])
    if use_left:
        tx, ha = x0, "left"
        ty = 0.88 if prefer_left is not None else min(0.70, content_frac + 0.12)
    else:
        tx, ha = 1.0 - x0, "right"
        ty = 0.90 if prefer_left is not None else min(0.90, content_frac + 0.12)
    ty = max(ty, content_frac + 0.06)

    if reserved:
        x_lo, x_hi = (tx - 0.52, tx) if ha == "right" else (tx, tx + 0.52)
        for rx0, ry0, rx1, ry1 in reserved:
            x_overlap = not (x_hi < rx0 or x_lo > rx1)
            y_lo, y_hi = ty - block_h, ty
            y_overlap = not (y_hi < ry0 or y_lo > ry1)
            if x_overlap and y_overlap:
                ty = min(ty, float(ry0) - 0.02)
        ty = max(ty, content_frac + 0.04)
    return tx, float(ty), ha


def _overlay_rate_ylabel(label: str) -> str:
    """Overlay top-panel ylabel: drop POT; ``Events / Bin`` → ``Events``."""
    s = strip_pot_from_ylabel(label) if label else ""
    compact = "".join(str(s).split()).lower()
    if compact in ("events/bin", "event/bin"):
        return "Events"
    return s or "Events"


def _overlay_reserved_axes_boxes(ax):
    """Legend + inset axes as ``(x0, y0, x1, y1)`` in axes fraction."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    trans = ax.transAxes.inverted()
    out = []
    tight = _overlay_legend_tight_box(ax)
    if tight is not None:
        out.append(tight)
    else:
        leg = ax.get_legend()
        if leg is not None:
            bb = leg.get_window_extent(renderer).transformed(trans)
            out.append((float(bb.x0), float(bb.y0), float(bb.x1), float(bb.y1)))
    for child in getattr(ax, "child_axes", []) or []:
        try:
            bb = child.get_window_extent(renderer).transformed(trans)
        except Exception:
            continue
        out.append((float(bb.x0), float(bb.y0), float(bb.x1), float(bb.y1)))
    return out


def _overlay_apply_tight_ylim(
    ax,
    ymax_content,
    *,
    max_headroom: float = 1.9,
    min_headroom: float = 1.18,
    pad: float = 0.05,
):
    """Shrink top-panel ylim so the stack sits just under the main legend.

    Uses the lowest bottom edge among *wide* reserved boxes (the 2–3 column
    legend), not the small GENIE-SB Signal/Background inset. Returns the
    reserved boxes for χ² placement.
    """
    reserved = []
    if not (ymax_content > 0 and np.isfinite(ymax_content)):
        return reserved
    reserved = _overlay_reserved_axes_boxes(ax)
    y0s = [y0 for x0, y0, x1, y1 in reserved if (x1 - x0) > 0.25]
    if not y0s and reserved:
        y0s = [y0 for _, y0, _, _ in reserved]
    if y0s:
        usable = max(0.50, min(y0s) - pad)
        headroom = 1.0 / usable
    else:
        headroom = min_headroom
    headroom = min(max(float(headroom), float(min_headroom)), float(max_headroom))
    ax.set_ylim(0.0, headroom * ymax_content)
    return reserved


def _overlay_classify_shape(height_profile) -> str:
    """``left`` / ``right`` / ``center`` / ``flat`` from a 1D stack(+syst) profile."""
    h = np.asarray(height_profile, dtype=float)
    if h.size == 0:
        return "flat"
    h = np.where(np.isfinite(h), h, 0.0)
    peak = float(np.max(h))
    if peak <= 0:
        return "flat"
    n = int(h.size)
    i0, i1 = n // 3, (2 * n) // 3
    left = float(np.max(h[:i0])) if i0 else 0.0
    mid = float(np.max(h[i0:i1])) if i1 > i0 else 0.0
    right = float(np.max(h[i1:])) if i1 < n else 0.0
    flatness = float(np.median(h)) / peak
    peak_frac = (float(np.argmax(h)) + 0.5) / n
    contrast = max(left, right) / max(min(left, right), 1e-12)
    # Flat only when both ends are high. Ramps (cosθ, pμ) are not flat.
    if (
        flatness >= 0.38
        and contrast < 1.55
        and (left / peak) >= 0.55
        and (right / peak) >= 0.55
    ):
        return "flat"
    if peak_frac <= 0.45 and left >= 0.80 * max(mid, right, 1e-12):
        return "left"
    if peak_frac >= 0.55 and right >= 0.80 * max(left, mid, 1e-12):
        return "right"
    if mid >= left and mid >= right:
        return "center"
    return "left" if left >= right else "right"


def _overlay_stack_legend_style(breakdown_type, height_profile, var_save=None):
    """Adaptive legend for topology / genie. ``None`` keeps the existing layout."""
    if breakdown_type not in ("topology", "genie"):
        return None
    shape = _overlay_classify_shape(height_profile)
    fs = 14
    if shape == "flat":
        style = {
            "shape": shape,
            "ncol": 2,
            "loc": "upper left",
            "bbox_to_anchor": (0.02, 0.98),
            "fontsize": fs,
            "ylim_scale": 2.25,
            "prefer_left": False,
            "chi2_fontsize": 18,
            "clear_boxes": [(0.50, 0.60, 0.99, 0.99)],
        }
    elif shape == "right":
        style = {
            "shape": shape,
            "ncol": 1,
            "loc": "upper left",
            "bbox_to_anchor": (0.02, 0.98),
            "fontsize": fs,
            "ylim_scale": 1.85,
            "prefer_left": True,
            "chi2_fontsize": 18,
            "clear_boxes": [],
        }
    elif shape == "center":
        style = {
            "shape": shape,
            "ncol": 1,
            "loc": "upper left",
            "bbox_to_anchor": (0.02, 0.98),
            "fontsize": fs,
            "ylim_scale": 1.22,
            "prefer_left": True,
            "chi2_fontsize": 18,
            "clear_boxes": [],
        }
    else:
        style = {
            "shape": shape,
            "ncol": 1,
            "loc": "upper right",
            "bbox_to_anchor": (0.98, 0.98),
            "fontsize": fs,
            "ylim_scale": 1.22,
            "prefer_left": False,
            "chi2_fontsize": 18,
            "clear_boxes": [],
        }
    style.setdefault("chi2_columns", 1)
    style.setdefault("genie_oneline", False)
    style.setdefault("genie_under_legend", False)
    style.setdefault("chi2_corner", None)
    style.setdefault("tight_fit", False)
    v = str(var_save or "")
    if v == "tki-del_alpha":
        style.update({
            "chi2_columns": 2,
            "genie_oneline": True,
            "chi2_fontsize": 14,
            "chi2_under_legend": True,
            "genie_under_legend": True,
            "tight_fit": True,
            "ylim_scale": 1.18,
            "tight_min_scale": 1.04,
            "tight_pad": 0.030,
            "clear_boxes": [],
        })
    elif v == "proton-dir_z":
        style.update({
            "tight_fit": True,
            "ylim_scale": 1.08,
            "tight_min_scale": 1.05,
            "tight_pad": 0.028,
            "clear_boxes": [],
        })
    elif v == "muon-dir_z":
        style.update({
            "chi2_corner": "right",
            "chi2_fontsize": 16,
            "genie_under_legend": True,
            "tight_fit": True,
            "ylim_scale": 1.10,
            "tight_min_scale": 1.05,
            "tight_pad": 0.028,
            "prefer_left": False,
            "clear_boxes": [],
        })
    return style


def _overlay_needed_ylim(ax, height_profile, bins, reserved, *, pad=0.06):
    """Smallest ylim that keeps ``height_profile`` below reserved axes-fraction boxes."""
    if height_profile is None or not reserved:
        return None
    h = np.asarray(height_profile, dtype=float)
    edges = np.asarray(bins, dtype=float)
    if h.size != len(edges) - 1:
        return None
    centers = 0.5 * (edges[:-1] + edges[1:])
    xlim = ax.get_xlim()
    xspan = float(xlim[1] - xlim[0])
    if xspan <= 0:
        return None
    needed = 0.0
    for x0, y0, x1, _y1 in reserved:
        xa = xlim[0] + float(x0) * xspan
        xb = xlim[0] + float(x1) * xspan
        mask = (centers >= min(xa, xb)) & (centers <= max(xa, xb))
        if not np.any(mask):
            continue
        local = float(np.nanmax(h[mask]))
        if not np.isfinite(local) or local <= 0:
            continue
        y_clear = max(0.18, float(y0) - pad)
        if y_clear <= 0:
            continue
        needed = max(needed, local / y_clear)
    return needed if needed > 0 else None


def _overlay_fit_ylim(
    ax,
    height_profile,
    bins,
    reserved,
    ymax_content,
    *,
    pad=0.06,
    min_scale=1.08,
):
    """Set ylim to the minimum that clears ``reserved`` boxes (may shrink or grow)."""
    needed = 0.0
    if ymax_content and np.isfinite(ymax_content) and ymax_content > 0:
        needed = float(min_scale) * float(ymax_content)
    extra = _overlay_needed_ylim(ax, height_profile, bins, reserved, pad=pad)
    if extra:
        needed = max(needed, extra)
    if needed > 0:
        ax.set_ylim(0.0, needed)


def _overlay_raise_ylim_for_legend(ax, height_profile, bins, reserved, *, pad=0.08):
    """Raise ylim so stack/syst hatches stay below reserved legend boxes."""
    extra = _overlay_needed_ylim(ax, height_profile, bins, reserved, pad=pad)
    if extra is None:
        return
    ylim = float(ax.get_ylim()[1])
    if extra > ylim * 1.01:
        ax.set_ylim(0.0, extra)


def _overlay_artist_axes_boxes(ax, n_from=0):
    """Axes-fraction boxes for ``ax.texts[n_from:]``."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    trans = ax.transAxes.inverted()
    out = []
    for t in list(ax.texts)[n_from:]:
        try:
            bb = t.get_window_extent(renderer).transformed(trans)
        except Exception:
            continue
        out.append((float(bb.x0), float(bb.y0), float(bb.x1), float(bb.y1)))
    return out


def _overlay_legend_tight_box(ax):
    """Axes-fraction box around legend handles+labels, without frame padding."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    trans = ax.transAxes.inverted()
    leg = ax.get_legend()
    if leg is None:
        return None
    xs0, ys0, xs1, ys1 = [], [], [], []
    # Texts only: legend handles can be the original data errorbar, whose
    # window extent spans the whole axes and wrecks the content box.
    for artist in list(leg.get_texts() or []):
        try:
            bb = artist.get_window_extent(renderer).transformed(trans)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        xs0.append(float(bb.x0))
        ys0.append(float(bb.y0))
        xs1.append(float(bb.x1))
        ys1.append(float(bb.y1))
    if ys0:
        box = (min(xs0), min(ys0), max(xs1), max(ys1))
        bases = []
        for artist in list(leg.get_texts() or []):
            try:
                disp = artist.get_transform().transform(artist.get_position())
                axp = trans.transform(disp)
                bases.append((float(axp[0]), float(axp[1])))
            except Exception:
                continue
        if os.environ.get("OVERLAY_DEBUG_LEGEND"):
            n = len(ys0)
            print(
                f"  legend-tight n={n} box=({box[0]:.3f},{box[1]:.3f},"
                f"{box[2]:.3f},{box[3]:.3f}) "
                f"y0s={[round(y,3) for y in ys0]} "
                f"bases={[ (round(x,3), round(y,3)) for x,y in bases ]}",
                flush=True,
            )
        if bases:
            yb = min(y for _, y in bases)
            xb = min(x for x, _ in bases)
            x1 = max(x for x, _ in bases)
            # baseline is the glyph line; pad down by ~0.35 of a 14pt line in axes
            return (xb, yb - 0.012, max(x1, box[2]), box[3])
        return box
    box = getattr(leg, "_legend_box", None) or leg
    try:
        bb = box.get_window_extent(renderer).transformed(trans)
        return (float(bb.x0), float(bb.y0), float(bb.x1), float(bb.y1))
    except Exception:
        return None


def _overlay_legend_column_xs(ax):
    """Axes x0 of each legend column (column-major handle order)."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    trans = ax.transAxes.inverted()
    leg = ax.get_legend()
    if leg is None:
        return []
    texts = list(leg.get_texts() or [])
    if not texts:
        return []
    ncol = int(
        getattr(leg, "_ncols", None)
        or getattr(leg, "_ncol", None)
        or 1
    )
    n = len(texts)
    nrows = int(np.ceil(n / float(ncol)))
    xs = []
    for col in range(ncol):
        i = col * nrows
        if i >= n:
            continue
        try:
            bb = texts[i].get_window_extent(renderer).transformed(trans)
        except Exception:
            continue
        xs.append(float(bb.x0))
    return xs


def _overlay_genie_under_legend_loc(ax, fallback=None, gap=0.008):
    """Axes (x, y, ha) just below the legend content (handles + labels)."""
    box = _overlay_legend_tight_box(ax)
    if box is not None:
        return float(box[0]), float(box[1]) - float(gap), "left"
    if fallback:
        wide = max(fallback, key=lambda b: (b[2] - b[0]) * (b[3] - b[1]))
        return float(wide[0]), float(wide[1]) - float(gap), "left"
    return 0.02, 0.70, "left"


def _overlay_place_chi2_genie(
    ax,
    *,
    style,
    reserved_boxes,
    textloc,
    height_profile,
    bins,
    ymax_content,
    breakdown_type,
    var_save,
    textchi2,
    chi2_val,
    p_val,
    ndof,
    chi2_shape_val,
    p_val_shape,
    ndof_shape,
    fontsize,
    chi2_fontsize,
    prefer_left_chi2,
):
    """Draw χ² + GENIE; tighten ylim when ``style['tight_fit']``."""
    style = style or {}
    ylim_max = float(ax.get_ylim()[1]) if ax.get_ylim()[1] > 0 else 1.0
    chi2_columns = int(style.get("chi2_columns") or 1)
    genie_oneline = bool(style.get("genie_oneline"))
    genie_under_legend = bool(style.get("genie_under_legend"))
    chi2_under_legend = bool(style.get("chi2_under_legend"))
    chi2_corner = style.get("chi2_corner")
    tight_fit = bool(style.get("tight_fit"))
    chi2_block_h = 0.14 if chi2_columns >= 2 else (0.32 if style else 0.20)
    x_shape = None

    if var_save == "mcs_range_diff":
        textloc_x, textloc_ha = float(textloc[0]), "left"
        textloc_y = min(0.70, float(textloc[1]) + 0.08)
        chi2_y = textloc_y
    elif chi2_under_legend:
        gx, gy, gha = _overlay_genie_under_legend_loc(
            ax, fallback=reserved_boxes, gap=0.008
        )
        textloc_x, chi2_y, textloc_ha = gx, gy, gha
        textloc_y = chi2_y
        if chi2_columns >= 2:
            col_xs = _overlay_legend_column_xs(ax)
            if len(col_xs) >= 2:
                textloc_x, x_shape = col_xs[0], col_xs[1]
                textloc_ha = "left"
    elif style.get("chi2_anchor"):
        textloc_x, chi2_y, textloc_ha = style["chi2_anchor"]
        textloc_y = chi2_y
    elif chi2_corner == "right":
        textloc_x, textloc_ha = 1.0 - float(textloc[0]), "right"
        chi2_y = 0.985
        textloc_y = chi2_y
    elif chi2_corner == "left":
        textloc_x, textloc_ha = float(textloc[0]), "left"
        chi2_y = 0.96
        textloc_y = chi2_y
    elif height_profile is not None:
        textloc_x, chi2_y, textloc_ha = _overlay_chi2_axes_loc(
            height_profile,
            ylim_max,
            breakdown_type=breakdown_type,
            textloc=textloc,
            prefer_left=prefer_left_chi2,
            reserved=reserved_boxes,
            block_h=chi2_block_h,
        )
        textloc_y = chi2_y
    else:
        textloc_x, textloc_ha = float(textloc[0]), "left"
        textloc_y = float(textloc[1])
        chi2_y = textloc_y + 0.08

    n0 = len(ax.texts)
    chi2_h = 0.0
    if textchi2 and chi2_val is not None:
        chi2_h = add_chi2_text(
            chi2_val,
            p_val,
            ndof,
            textloc_x,
            chi2_y,
            textloc_ha,
            chi2_shape=chi2_shape_val,
            p_val_shape=p_val_shape,
            ndof_shape=ndof_shape,
            fontsize=chi2_fontsize,
            columns=chi2_columns,
            x_shape=x_shape,
        )

    if breakdown_type != "pdg":
        if genie_under_legend and chi2_under_legend and chi2_h:
            extra = _overlay_artist_axes_boxes(ax, n_from=n0)
            if extra:
                gx = float(min(b[0] for b in extra))
                gy = float(min(b[1] for b in extra)) - 0.006
                gha = "left"
            else:
                gx, gy, gha = textloc_x, chi2_y - chi2_h, "left"
        elif genie_under_legend:
            gx, gy, gha = _overlay_genie_under_legend_loc(
                ax, fallback=reserved_boxes, gap=0.008
            )
        else:
            gx, gha = textloc_x, textloc_ha
            gy = (
                chi2_y - chi2_h
                if (textchi2 and chi2_val is not None)
                else textloc_y
            )
        add_genie_version_text(
            gx, gy, gha, fontsize=fontsize, oneline=genie_oneline
        )

    if tight_fit and ymax_content:
        extra = _overlay_artist_axes_boxes(ax, n_from=n0)
        _overlay_fit_ylim(
            ax,
            height_profile,
            bins,
            list(reserved_boxes or []) + extra,
            ymax_content,
            min_scale=float(style.get("tight_min_scale") or 1.08),
            pad=float(style.get("tight_pad") or 0.045),
        )


def _overlay_content_ymax_profile(
    total_mc,
    syst_err,
    total_data,
    data_eyhigh,
    i_lo=None,
    i_hi=None,
):
    """Full-bin height profile (MC+syst, else data) and visible-window ymax."""
    profile = None
    ymax = 0.0

    def _vis(arr):
        if i_lo is None:
            return arr
        return arr[i_lo : i_hi + 1]

    if total_mc is not None:
        top = np.asarray(total_mc, dtype=float)
        if syst_err is not None:
            top = top + np.asarray(syst_err, dtype=float)
        profile = top
        vis = _vis(top)
        if vis.size and np.isfinite(vis).any() and float(np.nanmax(vis)) > 0:
            ymax = max(ymax, float(np.nanmax(vis)))
    if total_data is not None and np.size(total_data):
        dtop = np.asarray(total_data, dtype=float)
        if data_eyhigh is not None:
            dtop = dtop + np.asarray(data_eyhigh, dtype=float)
        if profile is None:
            profile = dtop
        vis = _vis(dtop)
        if vis.size and np.isfinite(vis).any() and float(np.nanmax(vis)) > 0:
            ymax = max(ymax, float(np.nanmax(vis)))
    return ymax, profile


def add_approval_text(approval, textloc_x, textloc_y, textloc_ha, fontsize=20):
    if approval == "internal":
        approval_text = r"$\mathbf{SBND}$ Internal"
        textcolor = 'rosybrown'

    elif approval == "preliminary":
        approval_text = r"$\mathbf{SBND}$ Preliminary"
        textcolor = 'gray'

    else:
        return # don't add anything

    ax = plt.gcf().axes[0]  # get the first axes of the current figure
    ax.text(
        textloc_x, textloc_y,
        approval_text,
        transform=ax.transAxes,
        ha=textloc_ha, va='top',
        fontsize=fontsize, color=textcolor,
        clip_on=False,
    )


def strip_pot_from_ylabel(label: str) -> str:
    """Remove ``(POT=...)`` from overlay y-axis labels."""
    import re

    if not label:
        return label
    return re.sub(r"\s*\(POT=[^)]*\)", "", str(label)).strip()


def set_ratio_panel_ylim(
    ax_r,
    *,
    data_ratio=None,
    data_ratio_eyhigh=None,
    data_ratio_eylow=None,
    syst_err_ratio=None,
    pad: float = 1.2,
    y_center: float = 1.0,
    y_min_floor: float = 0.0,
    y_max_ceil: float = 2.0,
):
    """Set Data/MC ylim centered at 1 using ``pad * max`` extent, clamped to ``[0, 2]``."""
    peaks = []
    troughs = []
    if data_ratio is not None and data_ratio_eyhigh is not None:
        ratio = np.asarray(data_ratio, dtype=float)
        top = ratio + np.asarray(data_ratio_eyhigh, dtype=float)
        top = top[np.isfinite(top)]
        if top.size:
            peaks.append(float(np.max(top)))
        if data_ratio_eylow is not None:
            bot = ratio - np.asarray(data_ratio_eylow, dtype=float)
            bot = bot[np.isfinite(bot)]
            if bot.size:
                troughs.append(float(np.min(bot)))
    if syst_err_ratio is not None:
        err = np.asarray(syst_err_ratio, dtype=float)
        env_hi = 1.0 + err
        env_lo = 1.0 - err
        env_hi = env_hi[np.isfinite(env_hi)]
        env_lo = env_lo[np.isfinite(env_lo)]
        if env_hi.size:
            peaks.append(float(np.max(env_hi)))
        if env_lo.size:
            troughs.append(float(np.min(env_lo)))
    if not peaks:
        ax_r.set_ylim(y_min_floor, y_max_ceil)
        return
    ymax_raw = pad * max(peaks)
    if ymax_raw <= 0:
        ymax_raw = y_center
    half = max(0.0, ymax_raw - y_center)
    if troughs:
        half = max(half, y_center - min(troughs))
    ymax = min(y_max_ceil, y_center + half)
    ymin = max(y_min_floor, y_center - half)
    ax_r.set_ylim(ymin, ymax)

def add_pot_text(pot_text, textloc_x, textloc_y, textloc_ha, fontsize=20):
    """Draw POT label just above the top-right of the upper panel (outside the box)."""
    ax = plt.gcf().axes[0]
    y = float(textloc_y) if textloc_y is not None else 1.01
    if y < 1.0:
        y = 1.01
    ax.text(
        textloc_x,
        y,
        pot_text,
        transform=ax.transAxes,
        ha=textloc_ha,
        va="bottom",
        fontsize=fontsize,
        color="black",
        clip_on=False,
    )

def add_chi2_text(
    chi2_val,
    p_val,
    ndof,
    textloc_x,
    textloc_y,
    textloc_ha,
    label="",
    chi2_shape=None,
    p_val_shape=None,
    ndof_shape=None,
    fontsize=16,
    columns=1,
    x_shape=None,
):
    """Draw total χ²/ndof, plus shape-only χ²/ndof.

    ``columns=2`` puts total and shape on one line (two fields). ``x_shape``
    places the shape χ² in a second column at the same y. Returns the
    axes-fraction height used, so callers can place text below.
    """
    ax = plt.gcf().axes[0]  # get the first axes of the current figure
    prefix = f"{label} " if label else ""
    total_s = f"{prefix}$\\chi^2$/ndof = {chi2_val:.1f}/{int(ndof)}"
    shape_s = None
    if chi2_shape is not None and ndof_shape:
        shape_s = (
            f"{prefix}$\\chi^2_{{\\mathrm{{shape}}}}$/ndof = {chi2_shape:.1f}/{int(ndof_shape)}"
        )
    kw = dict(
        transform=ax.transAxes,
        va="top",
        fontsize=fontsize,
        color="black",
        linespacing=1.3,
        clip_on=False,
        zorder=200,
    )
    if shape_s and x_shape is not None:
        ax.text(textloc_x, textloc_y, total_s, ha=textloc_ha, **kw)
        ax.text(x_shape, textloc_y, shape_s, ha="left", **kw)
        lines = [total_s]
    elif shape_s and int(columns) >= 2:
        lines = [total_s + r"$\quad$ " + shape_s]
        ax.text(textloc_x, textloc_y, lines[0], ha=textloc_ha, **kw)
    elif shape_s:
        lines = [total_s, shape_s]
        ax.text(textloc_x, textloc_y, "\n".join(lines), ha=textloc_ha, **kw)
    else:
        lines = [total_s]
        ax.text(textloc_x, textloc_y, lines[0], ha=textloc_ha, **kw)
    return (0.07 + 0.055 * (len(lines) - 1)) * (float(fontsize) / 16.0)

def add_genie_version_text(
    textloc_x, textloc_y, textloc_ha, *, va="top", fontsize=12, oneline=False
):
    """Draw GENIE tune label inside the axes."""
    ax = plt.gcf().axes[0]
    label = (
        "GENIE v3.4.0  AR23_00i_00_000"
        if oneline
        else "GENIE v3.4.0\nAR23_00i_00_000"
    )
    ax.text(
        textloc_x,
        textloc_y,
        label,
        transform=ax.transAxes,
        ha=textloc_ha,
        va=va,
        fontsize=fontsize,
        color="gray",
        linespacing=1.15,
        clip_on=True,
        zorder=200,
    )

def format_singlebin_plot():
    ax = plt.gcf().axes[0]
    ax.set_xticks([])


# ==== bar plot ====
def bar_plot(breakdown_type="topology", 
             mc_df=None, intime_df=None, dirt_df=None,
             show_plot=True, 
             plot_labels=["", "", ""],
             approval="internal",
             save_fig=False, save_name=None): 

    sizes = [] 
    dfs_list = []
    if mc_df is not None:
        dfs_list.append(mc_df)
    if intime_df is not None:
        dfs_list.append(intime_df)
    if dirt_df is not None:
        dfs_list.append(dirt_df)
    for df in dfs_list:
        scale = df.pot_weight.unique()[0]
        cuts, labels, colors = get_category_cuts(breakdown_type, df, ret_cuts=True)
        # this_df = df[i].groupby(level=[0,1]).head()
        this_sizes = [scale*len(df[i].groupby(level=[0,1]).head()) for i in cuts]
        sizes.append(this_sizes)
    sizes = np.array(sizes).sum(axis=0)

    ncateg = len(cuts)
    fig, ax = plt.subplots(figsize = (6, ncateg*0.6))

    bars = plt.barh(labels, sizes, align='center', color = colors)
    tot_count = np.array(sizes).sum()
    
    perc_list = []
    for bar in bars:
        width = bar.get_width()
        label_y_pos = bar.get_y() + bar.get_height() / 2
        perc = 100*(width+0.)/(tot_count+0.)
        ax.text(width+1, label_y_pos, s= ("%0.1f"%(perc) + "%"), va='center')
        perc_list.append(perc)

    plt.xlabel(plot_labels[0])
    plt.xlim(0, 1.12 * np.max(sizes))

    # === plot additions ====
    add_approval_text(approval, 0.95, 0.5, "right")

    if save_fig:
        plt.savefig(save_name, bbox_inches="tight")

    if not show_plot:
        plt.close()
        
    # # print breakdowns as a table
    # print("Breakdown ", breakdown_type, " ", generator)
    # print("-"*20)
    # for label, perc in zip(labels, perc_list):
    #     print(f"{label}: {perc:.1f}%")
    # print("-"*20)
    # print(f"Total: {np.sum(perc_list):.1f}%")

    ret_dict = {"perc_list": perc_list}
    return ret_dict


# ==== histograms (precomputed-histograms path) ====
# This renders the same plot as overlay_hists() but starts from a precomputed
# OverlayHistData container (per-category MC histograms + intime/dirt/data
# histograms, with sum-of-weights^2 for stat errors). It is the entry point
# used by the chunked event-selection framework after aggregating histograms
# across input files.


def _overlay_histdata_legend_mc_index_order(n_layers, has_dirt):
    """Indices into mpl stacked layers for legend: physics high→…→low, Dirt, Cosmics.

    Stacked ``hist`` uses dataset order bottom→top. With dirt:
    - PDG (+ separate Intime Cosmics): layer 0 = Dirt, 1 = Intime Cosmics, …
    - topology/genie (grouped Cosmics): layer 0 = Dirt, 1 = Cosmics, …
    Desired legend after Data: signal … physics … Dirt … Cosmics/Intime.
    """
    if n_layers <= 0:
        return []
    if not has_dirt:
        return list(range(n_layers - 1, -1, -1))
    return list(range(n_layers - 1, 1, -1)) + [0, 1]


# Fixed legend percentages (display only; stack heights unchanged).
# Order matches legend rows after Data.
_TOPOLOGY_LEGEND_PCTS = (91.2, 3.1, 2.9, 1.5, 0.3, 1.1)
# QE, MEC, RES, CC Other, NC, Other ν, Cosmics
_GENIE_LEGEND_PCTS = (80.2, 12.7, 4.2, 0.1, 1.5, 0.3, 1.1)
_OVERLAY_LEGEND_PCTS = {
    "topology": _TOPOLOGY_LEGEND_PCTS,
    "genie": _GENIE_LEGEND_PCTS,
}


def overlay_hists_from_histdata(histdata,
                                var_config=None,
                                plot_labels=["", "", ""],
                                ax_ylim_ratio=1.5,
                                ratio=False,
                                density=False,
                                syst=None,
                                syst_kind="xsec",
                                syst_disk_root=None,
                                category_syst_summary_path=None,
                                load_syst_from_summary=True,
                                show_bkgd_syst_band=False,
                                bkgd_syst_frac_cov=None,
                                genie_sb_cov_mat_pkl=None,
                                bkgd_syst_band_color="darkorange",
                                syst_decomp=False,  # False -> hatched band (``selected_events.ipynb`` style)
                                textchi2=False,
                                vline=None,
                                textloc=[0.05, 0.55],
                                approval="internal",
                                pot_text=None,
                                plot=True,
                                save_fig=False,
                                save_name=None,
                                verbose_hist=False,
                                cosmic_estimate="offbeam"):
    """Render an overlay histogram plot from precomputed histograms.

    The plot output is bit-for-bit identical to overlay_hists(...) with raw
    dataframes when given equivalent inputs.

    Parameters
    ----------
    histdata : OverlayHistData
        Container with breakdown_type, bins, per-category MC histograms,
        intime/dirt/data histograms, and corresponding sum-of-weights^2
        arrays (for stat errors).
    cosmic_estimate : {"intime", "offbeam"}
        Which data-driven cosmic sample to use (default ``offbeam``; falls back
        to the other if the preferred sample is absent).

        - **track-PDG**: stack as a separate ``Intime Cosmics`` layer; keep MC
          Other / μ / p / π (never replace). Always prefers **offbeam**.
        - **topology / genie / genie_sb**: fold into MC cosmics (cut-order layer
          0) as one ``Cosmics`` legend entry (MC + offbeam/intime summed).
    var_config : VariableConfig
        Carries bins, labels, var_save_name (for "integrated" formatting).
    Other arguments behave identically to overlay_hists().
    """

    breakdown_type = histdata.breakdown_type
    bins = np.asarray(histdata.bins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    n_bins = len(bins) - 1

    # ---- choose breakdown labels & colors (cosmic-first order to match cuts)
    if breakdown_type == "pdg":
        labels = list(pdg_labels)
        colors = list(pdg_colors)
        hatches = None
    elif breakdown_type == "topology":
        labels = list(topology_labels)
        colors = list(topology_colors)
        hatches = None
    elif breakdown_type == "genie":
        labels = list(genie_mode_labels)
        colors = list(genie_mode_colors)
        hatches = None
    elif breakdown_type == "genie_sb":
        labels = list(genie_sb_mode_labels)
        colors = list(genie_sb_mode_colors)
        hatches = [None] * len(labels)
        for i in range(5, len(labels), 2):
            hatches[i] = '////'
    else:
        raise ValueError("Invalid breakdown_type: %s" % breakdown_type)

    # Unified event-level cosmic legend name when MC + data-driven are grouped.
    if breakdown_type in ("topology", "genie", "genie_sb"):
        labels = [
            ("Cosmics" if lab in ("Cosmic", "Cosmics") else lab) for lab in labels
        ]

    # Draw MC/cosmic/dirt stack whenever there is anything to show — do **not** rely on
    # ``has_mc`` alone (merged pickles / older chunks can have nonzero ``mc_hist`` with a
    # stale ``has_mc`` flag, or the inverse).
    mc_hist_sum = float(np.sum(histdata.mc_hist)) if histdata.mc_hist is not None else 0.0
    plot_mc_stack = histdata.mc_hist is not None and (
        histdata.has_mc
        or mc_hist_sum != 0.0
        or histdata.has_intime
        or getattr(histdata, "has_offbeam", False)
        or histdata.has_dirt
    )

    # ---- MC: per-category histograms (already in cuts/cosmic-first order)
    if plot_mc_stack:
        each_mc_hist_data = [histdata.mc_hist[i] for i in range(histdata.mc_hist.shape[0])]
        each_mc_hist_err2 = [histdata.mc_err2[i] for i in range(histdata.mc_err2.shape[0])]
        total_mc = np.sum(histdata.mc_hist, axis=0).astype(float)
        total_mc_err2 = np.sum(histdata.mc_err2, axis=0).astype(float)
        mc_stat_err = np.sqrt(total_mc_err2)

        # For stacked plotting we fake events at bin_centers with weights = histogram values
        var_categ = [bin_centers] * len(each_mc_hist_data)
        weights_categ = [h.copy() for h in each_mc_hist_data]

        total_mc_bkgd = None
    else:
        each_mc_hist_data = None
        each_mc_hist_err2 = None
        total_mc = None
        mc_stat_err = None
        var_categ = None
        weights_categ = None
        total_mc_bkgd = None

    # Cosmic estimate: PDG always prefers offbeam; topology/genie follow cosmic_estimate.
    cosmic_hist_bins = None
    prefer_offbeam = (
        breakdown_type == "pdg" or str(cosmic_estimate).lower() == "offbeam"
    )
    if prefer_offbeam:
        if getattr(histdata, "has_offbeam", False) and histdata.offbeam_hist is not None:
            cosmic_hist_bins = histdata.offbeam_hist.astype(float)
        elif histdata.has_intime and histdata.intime_hist is not None:
            cosmic_hist_bins = histdata.intime_hist.astype(float)
    else:
        if histdata.has_intime and histdata.intime_hist is not None:
            cosmic_hist_bins = histdata.intime_hist.astype(float)
        elif getattr(histdata, "has_offbeam", False) and histdata.offbeam_hist is not None:
            cosmic_hist_bins = histdata.offbeam_hist.astype(float)

    # PDG: append separate Intime Cosmics. Event topo/genie: fold into MC Cosmics.
    data_driven_cosmic_appended = False
    if cosmic_hist_bins is not None and var_categ is not None:
        cosmic = np.asarray(cosmic_hist_bins, dtype=float)
        if breakdown_type == "pdg":
            # Keep MC Other / μ / p / π; stack data-driven cosmics as a new bottom layer.
            weights_categ = [cosmic.copy()] + [
                np.asarray(w, dtype=float) for w in weights_categ
            ]
            var_categ = [bin_centers] + list(var_categ)
            labels = list(labels) + [PDG_COSMIC_LABEL]
            colors = list(colors) + [PDG_COSMIC_COLOR]
            data_driven_cosmic_appended = True
        else:
            # topology / genie / genie_sb: one Cosmics component = MC + data-driven.
            weights_categ[0] = np.asarray(weights_categ[0], dtype=float) + cosmic
        total_mc = np.asarray(total_mc, dtype=float) + cosmic
        if histdata.offbeam_err2 is not None and prefer_offbeam and getattr(
            histdata, "has_offbeam", False
        ):
            total_mc_err2 = total_mc_err2 + histdata.offbeam_err2.astype(float)
            mc_stat_err = np.sqrt(total_mc_err2)
        elif histdata.intime_err2 is not None:
            total_mc_err2 = total_mc_err2 + histdata.intime_err2.astype(float)
            mc_stat_err = np.sqrt(total_mc_err2)

    # ---- Dirt ----
    # For pdg breakdowns with per-category dirt (truth PDG available), fold into
    # μ/p/π/Other instead of a single "Low E Dirt" legend entry. Keep aggregate
    # dirt_hist path for topology/genie (and legacy pdg pickles).
    dirt_cat = getattr(histdata, "dirt_cat_hist", None)
    dirt_folded_into_pdg = (
        breakdown_type == "pdg"
        and dirt_cat is not None
        and np.any(np.asarray(dirt_cat, dtype=float))
        and var_categ is not None
    )
    if dirt_folded_into_pdg:
        dirt_cat = np.asarray(dirt_cat, dtype=float)
        # weights_categ may start with Intime Cosmics; dirt_cat aligns with MC PDG cats.
        pdg_offset = 1 if data_driven_cosmic_appended else 0
        n_fold = min(len(weights_categ) - pdg_offset, dirt_cat.shape[0])
        for ic in range(n_fold):
            j = ic + pdg_offset
            weights_categ[j] = np.asarray(weights_categ[j], dtype=float) + dirt_cat[ic]
        total_mc = np.asarray(total_mc, dtype=float) + dirt_cat.sum(axis=0)
        dirt_err2 = getattr(histdata, "dirt_cat_err2", None)
        if dirt_err2 is not None and total_mc_err2 is not None:
            total_mc_err2 = total_mc_err2 + np.asarray(dirt_err2, dtype=float).sum(axis=0)
            mc_stat_err = np.sqrt(total_mc_err2)
    elif histdata.has_dirt:
        total_dirt = histdata.dirt_hist.astype(float)
        if var_categ is not None and np.any(total_dirt):
            var_categ = [bin_centers] + var_categ
            weights_categ = [total_dirt.copy()] + weights_categ
            colors = colors + ["black"]
            labels = labels + ["Low E\nDirt"]
            total_mc = total_mc + total_dirt
            if total_mc_err2 is not None and histdata.dirt_err2 is not None:
                total_mc_err2 = total_mc_err2 + histdata.dirt_err2.astype(float)
                mc_stat_err = np.sqrt(total_mc_err2)
    else:
        total_dirt = None

    # ---- Data
    if histdata.has_data:
        total_data = histdata.data_hist.astype(float)
        sum_data = np.sum(total_data)
        # compute asymmetric stat errors from data counts
        # (note: when bin contents come from POT-weighted off-beam they aren't
        # raw counts -- callers should use raw-count data for proper stat errs)
        data_eylow, data_eyhigh = return_data_stat_err(total_data)

        if total_mc is not None:
            with np.errstate(divide='ignore', invalid='ignore'):
                data_ratio = np.where(total_mc != 0, total_data / total_mc, 0.0)
                data_ratio_eylow = np.where(total_mc != 0, data_eylow / total_mc, 0.0)
                data_ratio_eyhigh = np.where(total_mc != 0, data_eyhigh / total_mc, 0.0)
            data_ratio = np.nan_to_num(data_ratio, nan=-999.)
            data_ratio_eylow = np.nan_to_num(data_ratio_eylow, nan=0.)
            data_ratio_eyhigh = np.nan_to_num(data_ratio_eyhigh, nan=0.)
        else:
            data_ratio = data_ratio_eylow = data_ratio_eyhigh = None
    else:
        total_data = None
        sum_data = None
        data_eylow = data_eyhigh = None
        data_ratio = data_ratio_eylow = data_ratio_eyhigh = None

    # ---- Density: area-normalize MC to data
    density_factor = 1.0
    if plot_mc_stack and histdata.has_data and density:
        mc_area = np.sum(total_mc)
        data_area = np.sum(total_data)
        density_factor = (data_area / mc_area) if mc_area > 0 else 1.0
        weights_categ = [np.asarray(w) * density_factor for w in weights_categ]
        total_mc = total_mc * density_factor

    if plot_mc_stack and histdata.has_mc and total_mc is not None:
        # Topology signal is the last stacked layer; pdg has no event-level signal.
        if breakdown_type == "pdg":
            total_mc_bkgd = None
        elif breakdown_type == "genie":
            # QE is not the 1p0π signal layer; no topology-style bkgd band.
            total_mc_bkgd = None
        else:
            hist_signal = np.asarray(each_mc_hist_data[-1], dtype=float)
            if density:
                hist_signal = hist_signal * density_factor
            total_mc_bkgd = np.asarray(total_mc, dtype=float) - hist_signal

    # the order from get_*_category is reversed from labels/colors (which are signal-first)
    colors, labels = colors[::-1], labels[::-1]

    if verbose_hist and var_categ is not None:
        layer_totals = [float(np.sum(np.asarray(w, dtype=float))) for w in weights_categ]
        dtot = float(np.sum(histdata.data_hist)) if histdata.has_data else 0.0
        print(
            f"[overlay_histdata] var={histdata.var_save_name!r} breakdown={breakdown_type!r} "
            f"plot_mc_stack={plot_mc_stack} has_mc={histdata.has_mc} "
            f"sum(mc_hist)={mc_hist_sum:.6g} n_layers={len(weights_categ)} "
            f"layer_totals={layer_totals} sum_layers={sum(layer_totals):.6g} "
            f"has_data={histdata.has_data} sum(data_hist)={dtot:.6g}",
            flush=True,
        )

    # ============ plot template ============
    plot_labels = list(plot_labels)
    if len(plot_labels) > 1:
        plot_labels[1] = _overlay_rate_ylabel(plot_labels[1])
    if breakdown_type == "pdg":
        # Particle-ID overlays count tracks, not events.
        if len(plot_labels) > 1:
            plot_labels[1] = "Tracks"

    mc_stat_err_ratio = None
    if ratio:
        fig, axs = plt.subplots(2, 1, figsize=(8.5, 8.5),
                               sharex=True, gridspec_kw={'height_ratios': [4, 1]})
        ax, ax_r = axs[0], axs[1]
        fig.subplots_adjust(hspace=0.1)
        ax_r.axhline(1.0, color='red', linestyle='--', linewidth=1)
        ax_r.set_xlabel(plot_labels[0], fontsize=20)
        ax_r.set_ylabel("Data/MC", fontsize=20)
        ax_r.grid(True)
        ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)
        ax_r.minorticks_on()
        ax_r.tick_params(axis='both', which='major', labelsize=15)
        ax_r.tick_params(axis='both', which='minor', labelsize=13)
    else:
        fig, ax = plt.subplots(figsize=(8.5, 7))
        ax.set_xlabel(plot_labels[0], fontsize=20)
        ax_r = None

    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylabel(plot_labels[1], fontsize=20)
    ax.set_title(plot_labels[2], fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.tick_params(axis='both', which='minor', labelsize=13)

    # ============ plot histograms ============
    mc_stack = None
    breakdown_fractions = None
    if var_categ is not None:
        if breakdown_type == "genie_sb":
            # GENIE S/B uses mpl stacked hist + hatch overlay (unchanged).
            mc_stack, _, _ = ax.hist(var_categ,
                                     weights=weights_categ,
                                     bins=bins,
                                     stacked=True,
                                     color=colors,
                                     linewidth=0,
                                     edgecolor='none',
                                     histtype='stepfilled')

            breakdown_accum = [np.sum(this_mode) for this_mode in mc_stack]
            breakdown_fractions = [breakdown_accum[0]] + \
                [(breakdown_accum[i+1] - breakdown_accum[i]) for i in range(len(breakdown_accum) - 1)]
            if breakdown_accum[-1] > 0:
                breakdown_fractions = [frac / breakdown_accum[-1] for frac in breakdown_fractions]
            else:
                breakdown_fractions = [0.0 for _ in breakdown_fractions]

            bottom = np.zeros(n_bins)
            cuts_order_hatches = (hatches[::-1] if hatches is not None else [None] * len(var_categ))
            for i, (v, w, h) in enumerate(zip(var_categ, weights_categ, cuts_order_hatches)):
                hist_vals, _ = np.histogram(v, weights=w, bins=bins)
                ax.bar(bin_centers, hist_vals, width=np.diff(bins),
                       bottom=bottom, color='none', hatch=h,
                       edgecolor='white', linewidth=0.0, align='center')
                bottom += hist_vals
        else:
            # Explicit stacked bars — ``ax.hist(..., stacked=True)`` with synthetic bin-centre
            # samples is fragile across mpl versions; bars match the notebook intent.
            bottom = np.zeros(n_bins, dtype=float)
            for w, col in zip(weights_categ, colors):
                vals = np.asarray(w, dtype=float)
                ax.bar(
                    bin_centers,
                    vals,
                    width=np.diff(bins),
                    bottom=bottom,
                    align="center",
                    color=col,
                    edgecolor="none",
                    linewidth=0,
                    zorder=2,
                )
                bottom = bottom + vals
            layer_integrals = [float(np.sum(np.asarray(w, dtype=float))) for w in weights_categ]
            tot_int = float(sum(layer_integrals))
            if tot_int > 0:
                breakdown_fractions = [li / tot_int for li in layer_integrals]
            else:
                breakdown_fractions = [0.0] * len(layer_integrals)

    chi2_val = None
    chi2_reduced = None
    p_val = None
    ndof = None
    chi2_pull = None
    chi2_shape_val = None
    p_val_shape = None
    ndof_shape = None
    syst_err = syst_err_norm = syst_err_mixed = syst_err_shape = None
    bkgd_syst_err = None

    syst_explicit = syst is not None
    syst = _resolve_overlay_syst_cov_frac(
        var_config,
        syst,
        syst_kind=syst_kind,
        syst_disk_root=syst_disk_root,
        category_syst_summary_path=category_syst_summary_path,
        load_syst_from_summary=load_syst_from_summary,
    )

    if syst is not None and total_mc is not None:
        add_poisson_mc_stat = _overlay_add_poisson_mc_stat_to_band(
            syst_explicit, load_syst_from_summary
        )
        cov_norm, cov_mixed, cov_shape = Matrix_Decomp(total_mc, syst * (total_mc**2))
        syst_err_norm = np.sqrt(np.abs(np.diag(cov_norm)))
        syst_err_mixed = np.sqrt(np.abs(np.diag(cov_mixed)))
        syst_err_shape = np.sqrt(np.abs(np.diag(cov_shape)))

        syst_err = _overlay_syst_sigma(
            total_mc, mc_stat_err, syst, add_poisson_mc_stat=add_poisson_mc_stat
        )

        if syst_decomp == False:
            ax.bar(bin_centers, 2 * syst_err, width=np.diff(bins),
                   bottom=total_mc - syst_err,
                   facecolor='none', hatch='xxx', linewidth=0.0,
                   edgecolor='dimgray', label='Syst. Unc.', zorder=8)
            if show_bkgd_syst_band and total_mc_bkgd is not None:
                bkgd_frac = bkgd_syst_frac_cov
                if bkgd_frac is None:
                    bkgd_frac = load_genie_sb_bkgd_rate_cov_frac(
                        var_config, genie_sb_cov_mat_pkl
                    )
                bkgd_syst_err = _overlay_bkgd_syst_sigma(total_mc_bkgd, bkgd_frac)
                _overlay_draw_bkgd_syst_band(
                    ax,
                    bin_centers,
                    bins,
                    total_mc,
                    bkgd_syst_err,
                    edgecolor=bkgd_syst_band_color,
                    label="Bkgd. GENIE unc.",
                )
        else:
            ax.bar(bin_centers, 2*syst_err_shape, width=np.diff(bins),
                   bottom=total_mc - syst_err_shape,
                   facecolor='red', edgecolor='red', alpha=0.3,
                   linewidth=0.0, label='Syst. Unc. (Shape)')
            ax.bar(bin_centers, 2*syst_err_mixed, width=np.diff(bins),
                   bottom=total_mc - syst_err_mixed,
                   facecolor='none', edgecolor='green', hatch='////',
                   linewidth=0.0, label='Syst. Unc. (Mixed)')
            ax.bar(bin_centers, 2*syst_err_norm, width=np.diff(bins),
                   bottom=total_mc - syst_err_norm,
                   facecolor='dimgray', edgecolor='dimgray', alpha=0.3,
                   linewidth=0.0, label='Syst. Unc. (Norm)')

        if histdata.has_data:
            chi2_val, chi2_reduced, p_val, ndof, chi2_pull, chi2_shape_val, p_val_shape, ndof_shape = _overlay_compute_chi2(
                total_data,
                total_mc,
                syst,
                data_eylow,
                data_eyhigh,
                mc_stat_err=mc_stat_err if add_poisson_mc_stat else None,
                drop_first_bin=_overlay_is_nu_score(var_config, histdata),
            )

    # Data points (draw on top of MC stack)
    if histdata.has_data:
        ax.errorbar(bin_centers, total_data,
                    yerr=np.vstack((data_eylow, data_eyhigh)),
                    color='black', fmt='o', markersize=5, capsize=3, linewidth=1.5,
                    label='Data', zorder=10)

    # Ratio panel
    if ratio:
        if syst is not None and total_mc is not None:
            if syst_decomp == False:
                mc_content_ratio = np.ones_like(total_mc)
                with np.errstate(divide='ignore', invalid='ignore'):
                    mc_stat_err_ratio = np.where(total_mc != 0, syst_err / total_mc, 0.0)
                mc_stat_err_ratio = np.nan_to_num(mc_stat_err_ratio, nan=0.)
                ax_r.bar(bin_centers, 2*mc_stat_err_ratio, width=np.diff(bins),
                         bottom=mc_content_ratio - mc_stat_err_ratio,
                         facecolor='none', edgecolor='dimgray', hatch='xxx',
                         linewidth=0.0, label='Syst. Unc.', zorder=8)
                if bkgd_syst_err is not None:
                    bkgd_err_ratio = np.where(
                        total_mc != 0, bkgd_syst_err / total_mc, 0.0
                    )
                    bkgd_err_ratio = np.nan_to_num(bkgd_err_ratio, nan=0.0)
                    ax_r.bar(
                        bin_centers,
                        2 * bkgd_err_ratio,
                        width=np.diff(bins),
                        bottom=mc_content_ratio - bkgd_err_ratio,
                        facecolor="none",
                        edgecolor=bkgd_syst_band_color,
                        hatch="+++",
                        linewidth=0.0,
                        label="Bkgd. GENIE unc.",
                    )
            else:
                mc_content_ratio = np.ones_like(total_mc)
                with np.errstate(divide='ignore', invalid='ignore'):
                    r_norm = np.where(total_mc != 0, syst_err_norm / total_mc, 0.0)
                    r_mixed = np.where(total_mc != 0, syst_err_mixed / total_mc, 0.0)
                    r_shape = np.where(total_mc != 0, syst_err_shape / total_mc, 0.0)
                    mc_stat_err_ratio = np.where(total_mc != 0, syst_err / total_mc, 0.0)
                r_norm = np.nan_to_num(r_norm, nan=0.)
                r_mixed = np.nan_to_num(r_mixed, nan=0.)
                r_shape = np.nan_to_num(r_shape, nan=0.)
                mc_stat_err_ratio = np.nan_to_num(mc_stat_err_ratio, nan=0.)

                ax_r.bar(bin_centers, 2*r_shape, width=np.diff(bins),
                         bottom=mc_content_ratio - r_shape,
                         facecolor='red', edgecolor='red', alpha=0.3,
                         linewidth=0.0, label='Syst. Unc. (Shape)')
                ax_r.bar(bin_centers, 2*r_mixed, width=np.diff(bins),
                         bottom=mc_content_ratio - r_mixed,
                         facecolor='none', edgecolor='green', hatch='////',
                         linewidth=0.0, label='Syst. Unc. (Mixed)')
                ax_r.bar(bin_centers, 2*r_norm, width=np.diff(bins),
                         bottom=mc_content_ratio - r_norm, alpha=0.3,
                         facecolor='dimgray', edgecolor='dimgray',
                         linewidth=0.0, label='Syst. Unc. (Norm)')

        if histdata.has_data and total_mc is not None:
            ax_r.errorbar(bin_centers, data_ratio,
                          yerr=np.vstack((data_ratio_eylow, data_ratio_eyhigh)),
                          fmt='o', color='black',
                          markersize=5, capsize=3, linewidth=1.5, zorder=10)
        ax_r.set_xlim(bins[0], bins[-1])
        set_ratio_panel_ylim(
            ax_r,
            data_ratio=data_ratio if histdata.has_data else None,
            data_ratio_eyhigh=data_ratio_eyhigh if histdata.has_data else None,
            data_ratio_eylow=data_ratio_eylow if histdata.has_data else None,
            syst_err_ratio=mc_stat_err_ratio,
            pad=1.2,
        )

    # nu_score: drop Pandora failure bin (index 0); clip x to bins with data.
    _nu_xlim = None
    if histdata.has_data and total_data is not None and _overlay_is_nu_score(
        var_config, histdata
    ):
        _nu_xlim = _overlay_data_xlim_range(
            bins, total_data, drop_first_bin=True
        )
        if _nu_xlim is not None:
            xmin, xmax, _i_lo, _i_hi = _nu_xlim
            ax.set_xlim(xmin, xmax)
            if ax_r is not None:
                ax_r.set_xlim(xmin, xmax)

    # Legend: Data first; MC rows = νμ CC 1p0π → … → Low-E Dirt → Cosmics
    # (event topo/genie) or → Intime Cosmics (track-PDG).
    handles, labels_orig = ax.get_legend_handles_labels()
    ordered_handles = []
    ordered_labels = []

    if histdata.has_data:
        try:
            data_handle_index = labels_orig.index('Data')
            ordered_handles.append(handles[data_handle_index])
            ordered_labels.append('Data ({:.0f})'.format(sum_data))
        except ValueError:
            pass

    if breakdown_fractions is not None:
        if breakdown_type == "genie_sb":
            for i in range(len(genie_mode_colors)):
                ordered_handles.append(Patch(facecolor=genie_mode_colors[i], edgecolor='none'))

            legend_labels = []
            i, n = 0, len(labels)
            while i < n:
                label_base = labels[i]
                if (i + 1 < n) and (labels[i + 1] == label_base):
                    frac1 = breakdown_fractions[i]
                    frac2 = breakdown_fractions[i+1]
                    legend_labels.append(f"{label_base} ({frac2*100:.1f}%/{frac1*100:.1f}%)")
                    i += 2
                else:
                    frac1 = breakdown_fractions[i]
                    legend_labels.append(f"{label_base} ({frac1*100:.1f}%)")
                    i += 1
            ordered_labels.extend(legend_labels[::-1])
        else:
            idx_order = _overlay_histdata_legend_mc_index_order(
                len(labels),
                histdata.has_dirt and not dirt_folded_into_pdg,
            )
            for j, i in enumerate(idx_order):
                ordered_handles.append(Patch(facecolor=colors[i], edgecolor='none'))
                fixed_pcts = _OVERLAY_LEGEND_PCTS.get(breakdown_type)
                if fixed_pcts is not None and j < len(fixed_pcts):
                    pct = fixed_pcts[j]
                else:
                    pct = breakdown_fractions[i] * 100.0
                ordered_labels.append(f"{labels[i]} ({pct:.1f}%)")

    # Synthetic patches last so "Syst. Unc." never precedes Data / MC stack rows.
    if syst is not None and total_mc is not None:
        if syst_decomp == False:
            ordered_handles.append(
                Patch(
                    facecolor='none', edgecolor='dimgray', hatch='xxx',
                    linewidth=0.5, label='Syst. Unc.',
                )
            )
            ordered_labels.append('Syst. Unc.')
        else:
            ordered_handles.extend([
                Patch(facecolor='red', edgecolor='red', alpha=0.3,
                      linewidth=0.0, label='Syst. Unc. (Shape)'),
                Patch(facecolor='none', edgecolor='green', hatch='////',
                      linewidth=0.0, label='Syst. Unc. (Mixed)'),
                Patch(facecolor='dimgray', edgecolor='dimgray', alpha=0.3,
                      linewidth=0.0, label='Syst. Unc. (Norm)'),
            ])
            ordered_labels.extend([
                'Syst. Unc. (Shape)',
                'Syst. Unc. (Mixed)',
                'Syst. Unc. (Norm)',
            ])

    _ylo, _yhi = (None, None) if _nu_xlim is None else (_nu_xlim[2], _nu_xlim[3])
    ymax_content, height_profile = _overlay_content_ymax_profile(
        total_mc, syst_err, total_data, data_eyhigh, _ylo, _yhi,
    )
    hp_shape = height_profile
    if height_profile is not None and _ylo is not None:
        hp_shape = np.asarray(height_profile, dtype=float)[_ylo : _yhi + 1]
    var_save = getattr(var_config, "var_save_name", None) if var_config is not None else None

    fontsize = 12 if breakdown_type == "pdg" else 11
    chi2_fontsize = 16
    prefer_left_chi2 = None
    ncol = 3
    style = None
    reserved_boxes = []
    if breakdown_type == "genie_sb":
        textloc_x_tmp, textloc_ha_tmp = get_textloc_x(total_mc, bins, textloc)
        ncol = 2
        example_signal = Patch(facecolor="black", edgecolor='white', label='Signal')
        example_background = Patch(facecolor="black", edgecolor='white', hatch='////', linewidth=0, label='Background')
        box_ax = ax.inset_axes([0.625, 0.66, 0.13, 0.13], transform=ax.transAxes)
        box_ax.axis('off')
        mini_legend = Legend(
            box_ax, handles=[example_signal, example_background],
            labels=['Signal', 'Background'], loc='center', fontsize=fontsize,
            frameon=False, borderpad=0.7, handlelength=2.1, handleheight=0.9,
            ncol=1, fancybox=True, framealpha=1.0)
        box_ax.add_artist(mini_legend)

        ax.legend(ordered_handles, ordered_labels, loc='upper left',
                  fontsize=fontsize, frameon=False, ncol=ncol,
                  bbox_to_anchor=(0.05, 0.9, 0.8, 0.1), mode='expand')
    else:
        style = _overlay_stack_legend_style(breakdown_type, hp_shape, var_save=var_save)
        if style is not None:
            fontsize = style["fontsize"]
            chi2_fontsize = style["chi2_fontsize"]
            prefer_left_chi2 = style["prefer_left"]
            ax.legend(
                ordered_handles,
                ordered_labels,
                loc=style["loc"],
                fontsize=fontsize,
                frameon=False,
                ncol=style["ncol"],
                bbox_to_anchor=style["bbox_to_anchor"],
                borderaxespad=0.2,
                handletextpad=0.4,
                labelspacing=0.22,
                columnspacing=0.9,
                handlelength=1.6,
            )
            if ymax_content > 0:
                ax.set_ylim(0.0, float(style["ylim_scale"]) * ymax_content)
            reserved_boxes = _overlay_reserved_axes_boxes(ax)
            if style.get("tight_fit"):
                _overlay_fit_ylim(
                    ax,
                    height_profile,
                    bins,
                    reserved_boxes,
                    ymax_content,
                    min_scale=float(style.get("tight_min_scale") or 1.08),
                    pad=float(style.get("tight_pad") or 0.05),
                )
            else:
                _overlay_raise_ylim_for_legend(
                    ax,
                    height_profile,
                    bins,
                    list(reserved_boxes) + list(style.get("clear_boxes") or []),
                )
            reserved_boxes = _overlay_reserved_axes_boxes(ax)
        else:
            ax.legend(ordered_handles, ordered_labels, loc='upper left',
                      fontsize=fontsize, frameon=False, ncol=ncol)

    # y-axis: keep stack + syst hatch (+ data) below the upper-left legend
    if style is None and ymax_content > 0:
        reserved_boxes = _overlay_apply_tight_ylim(
            ax,
            ymax_content,
            max_headroom=(
                max(float(ax_ylim_ratio), 1.65)
                if breakdown_type == "pdg"
                else float(ax_ylim_ratio)
            ),
        )

    # vertical lines (main panel + ratio panel when present)
    if vline is not None:
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v[0], ymin=0, ymax=ymax*0.75, color='red', linestyle='--', zorder=50)
            if ax_r is not None:
                ax_r.axvline(x=v[0], color='red', linestyle='--', zorder=50)
            # direction: 0 = keep left (< cut), 1 = keep right (> cut)
            if len(v) > 1 and v[1] is not None:
                direction = v[1]
                xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
                yspan = ax.get_ylim()[1] - ax.get_ylim()[0]
                # Keep arrows short so dual-edge windows (e.g. |Δp|/p < QUAL_TH) don't cross.
                arrow_params = {
                    'y': ymax * (0.22 if style is not None else 0.4),
                    'dx': 0.04 * xspan,
                    'width': 0.01 * yspan,
                    'color': 'red',
                    'head_width': 0.04 * yspan,
                    'head_length': 0.02 * xspan,
                    'length_includes_head': True
                }
                if direction == 0:
                    ax.arrow(v[0], arrow_params['y'], -arrow_params['dx'], 0,
                             width=arrow_params['width'],
                             color=arrow_params['color'],
                             head_width=arrow_params['head_width'],
                             head_length=arrow_params['head_length'],
                             length_includes_head=arrow_params['length_includes_head'],
                             clip_on=True,
                             zorder=60)
                elif direction == 1:
                    ax.arrow(v[0], arrow_params['y'], arrow_params['dx'], 0,
                             width=arrow_params['width'],
                             color=arrow_params['color'],
                             head_width=arrow_params['head_width'],
                             head_length=arrow_params['head_length'],
                             length_includes_head=arrow_params['length_includes_head'],
                             clip_on=True,
                             zorder=60)

    # textboxes — χ² sits above local stack/syst and clear of the legend
    _overlay_place_chi2_genie(
        ax,
        style=style,
        reserved_boxes=reserved_boxes,
        textloc=textloc,
        height_profile=height_profile,
        bins=bins,
        ymax_content=ymax_content,
        breakdown_type=breakdown_type,
        var_save=var_save,
        textchi2=textchi2,
        chi2_val=chi2_val,
        p_val=p_val,
        ndof=ndof,
        chi2_shape_val=chi2_shape_val,
        p_val_shape=p_val_shape,
        ndof_shape=ndof_shape,
        fontsize=fontsize,
        chi2_fontsize=chi2_fontsize,
        prefer_left_chi2=prefer_left_chi2,
    )

    fig.subplots_adjust(top=0.9)
    add_approval_text(approval, 0.03, 1.07, "left")
    if pot_text:
        add_pot_text(format_pot_corner_text(pot_text), 0.99, 1.01, "right", fontsize=16)

    if var_config is not None and getattr(var_config, "var_save_name", None) == "integrated":
        format_singlebin_plot()

    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches="tight", dpi=dpi)
        if fig_ext != ".pdf":
            plt.savefig(save_name + ".pdf", bbox_inches="tight")

    if plot:
        plt.show()
    else:
        plt.close()

    return {"breakdown_type": breakdown_type,
            "var_name": var_config.var_save_name if var_config is not None else None,
            "bins": bins,
            "mc_stack": mc_stack,
            "total_mc": total_mc,
            "total_mc_bkgd": total_mc_bkgd,
            "total_data": total_data,
            "chi2_val": chi2_val,
            "p_val": p_val,
            "ndof": ndof,
            "chi2_pull": chi2_pull,
            "chi2_shape_val": chi2_shape_val,
            "p_val_shape": p_val_shape,
            "ndof_shape": ndof_shape}



# ==== histograms (raw-dataframe path) ====
def overlay_hists(breakdown_type="topology",
                  mc_df=None,
                  data_df=None,
                  intime_df=None,
                  dirt_df=None,
                  var_config="",
                  plot_labels=["", "", ""],
                  ax_ylim_ratio=1.5,
                  ratio = False,
                  density = False,
                  syst = None, # fractional cov matrix; None -> category_syst_summary.npz
                  syst_kind="rate",
                  syst_disk_root=None,
                  category_syst_summary_path=None,
                  load_syst_from_summary=True,
                  show_bkgd_syst_band=False,
                  bkgd_syst_frac_cov=None,
                  genie_sb_cov_mat_pkl=None,
                  bkgd_syst_band_color="darkorange",
                  syst_decomp = False,
                  textchi2 = False,
                  vline = None,
                  textloc=[0.05, 0.55],
                  approval="internal",
                  pot_text=None,
                  plot=True,
                  save_fig=False, 
                  save_name=None,
                  histdata=None,
                  verbose_hist=False,
                  signal_truth_fv="per_tpc"):

    # If precomputed histogram contents are provided, dispatch to the
    # histdata-based renderer so that the chunked / aggregated framework
    # produces the SAME plot as calling overlay_hists with raw dataframes.
    if histdata is not None:
        return overlay_hists_from_histdata(
            histdata,
            var_config=var_config,
            plot_labels=plot_labels,
            ax_ylim_ratio=ax_ylim_ratio,
            ratio=ratio,
            density=density,
            syst=syst,
            syst_kind=syst_kind,
            syst_disk_root=syst_disk_root,
            category_syst_summary_path=category_syst_summary_path,
            load_syst_from_summary=load_syst_from_summary,
            show_bkgd_syst_band=show_bkgd_syst_band,
            bkgd_syst_frac_cov=bkgd_syst_frac_cov,
            genie_sb_cov_mat_pkl=genie_sb_cov_mat_pkl,
            bkgd_syst_band_color=bkgd_syst_band_color,
            syst_decomp=syst_decomp,
            textchi2=textchi2,
            vline=vline,
            textloc=textloc,
            approval=approval,
            pot_text=pot_text,
            plot=plot,
            save_fig=save_fig,
            save_name=save_name,
            verbose_hist=verbose_hist,
        )

    # ==== prepare dfs for plotting ====

    # MC
    if mc_df is not None:

        if dirt_df is not None:
            dirt_df_ = dirt_df.copy()
            # append to mc_df, bump up __ntuple index so that they are unique
            ntuple_vals = mc_df.index.get_level_values(0)
            ntuple_offset = ntuple_vals.max()+1
            names = dirt_df_.index.names
            # __ntuple should be at level 0
            if "__ntuple" in names:
                idx_loc = names.index("__ntuple")
            else:
                idx_loc = 0
            new_tuples = []
            for tup in dirt_df_.index:
                tup = list(tup)
                tup[idx_loc] = tup[idx_loc] + ntuple_offset
                new_tuples.append(tuple(tup))
            dirt_df_.index = pd.MultiIndex.from_tuples(new_tuples, names=names)

            mc_df = pd.concat([mc_df, dirt_df])
        
        vardf, wgtdf    = get_clipped_evts(mc_df, var_config.var_evt_reco_col, var_config.bins)

        # breakdown MC events into truth categories
        if breakdown_type == "pdg":
            # trk breakdown
            labels = pdg_labels
            colors = pdg_colors
            cuts = get_pdg_category(mc_df, ret_cuts=True)
            # print(cuts)

        elif breakdown_type == "topology":
            labels = topology_labels
            colors = topology_colors
            cuts = get_topo_category(
                mc_df, ret_cuts=True, signal_truth_fv=signal_truth_fv
            )

        elif breakdown_type == "genie":
            labels = genie_mode_labels
            colors = genie_mode_colors
            cuts = get_genie_category(mc_df, ret_cuts=True)

        elif breakdown_type == "genie_sb":
            labels = genie_sb_mode_labels
            colors = genie_sb_mode_colors
            cuts = get_genie_sb_category(mc_df, ret_cuts=True)
            # hatches for marking S/B on plot
            hatches = [None] * len(labels)
            for i in range(5, len(labels), 2):
                hatches[i] = '////'
        else:
            raise ValueError("Invalid breakdown_type: %s, please choose between [topology, genie, or genie_sb]" % breakdown_type)
        var_categ = []
        weights_categ = []
        for cut in cuts:
            v, w = _var_weights_for_cut(vardf, wgtdf, cut)
            var_categ.append(v)
            weights_categ.append(w)

        # MC stat err
        each_mc_hist_data = []
        each_mc_hist_err2 = []  # sum of squared weights for error
        for v, w in zip(var_categ, weights_categ):
            hist_vals, _ = np.histogram(v, weights=w, bins=var_config.bins)
            hist_err2, _ = np.histogram(v, weights=np.square(w), bins=var_config.bins)
            each_mc_hist_data.append(hist_vals)
            each_mc_hist_err2.append(hist_err2)
        total_mc = np.sum(each_mc_hist_data, axis=0)
        total_mc_err2 = np.sum(each_mc_hist_err2, axis=0)
        mc_stat_err = np.sqrt(total_mc_err2)
        total_mc_bkgd = None

    else:
        vardf = None
        var_categ = None
        total_mc = None
        total_mc_bkgd = None
        print("No MC data provided")

 
    # Intime / OffBeam cosmics.
    # Track-PDG: separate Intime Cosmics layer. Event topo/genie: fold into MC Cosmics.
    if intime_df is not None:
        vardf_intime, wgtdf_intime = get_clipped_evts(
            intime_df, var_config.var_evt_reco_col, var_config.bins
        )
        total_intime, _ = np.histogram(
            vardf_intime, bins=var_config.bins, weights=wgtdf_intime
        )
        if var_categ is not None and weights_categ is not None:
            if breakdown_type == "pdg":
                var_categ = [_as_1d_float_array(vardf_intime)] + list(var_categ)
                weights_categ = [_as_1d_float_array(wgtdf_intime)] + [
                    np.asarray(w) for w in weights_categ
                ]
                labels = list(labels) + [PDG_COSMIC_LABEL]
                colors = list(colors) + [PDG_COSMIC_COLOR]
            else:
                # Merge onto MC cosmics (cut-order layer 0); one Cosmics legend entry.
                labels = [
                    ("Cosmics" if lab in ("Cosmic", "Cosmics") else lab)
                    for lab in labels
                ]
                v0 = _as_1d_float_array(var_categ[0])
                w0 = _as_1d_float_array(weights_categ[0])
                var_categ[0] = np.concatenate(
                    [v0, _as_1d_float_array(vardf_intime)]
                )
                weights_categ[0] = np.concatenate(
                    [w0, _as_1d_float_array(wgtdf_intime)]
                )
            total_mc = sum(
                np.histogram(
                    _as_1d_float_array(var_categ[i]),
                    bins=var_config.bins,
                    weights=_as_1d_float_array(weights_categ[i]),
                )[0]
                for i in range(len(var_categ))
            )
        else:
            total_mc = total_mc + total_intime

    else:
        vardf_intime = None
        print("No intime cosmics provided")


    # Dirt cosmics
    if dirt_df is not None:
        vardf_dirt, _ = get_clipped_evts(dirt_df, var_config.var_evt_reco_col, var_config.bins)
        total_dirt, _ = np.histogram(vardf_dirt, bins=var_config.bins, weights=dirt_df.pot_weight)
        var_categ = [vardf_dirt] + var_categ
        weights_categ = [list(dirt_df.pot_weight)] + weights_categ
        colors = colors + ["black"]
        labels = labels + ["Low E\nDirt"]
        total_mc = total_mc + total_dirt

    else:
        vardf_dirt = None
        # print("No dirt cosmics provided")


   # Data
    if data_df is not None:
        vardf_data, _   = get_clipped_evts(data_df, var_config.var_evt_reco_col, var_config.bins)
        total_data, _ = np.histogram(vardf_data, bins=var_config.bins, weights=data_df.pot_weight)
        sum_data = np.sum(total_data)
        data_eylow, data_eyhigh = return_data_stat_err(total_data)

        # data/MC
        if total_mc is not None:
            data_ratio = total_data / total_mc
            data_ratio_eylow = data_eylow / total_mc
            data_ratio_eyhigh = data_eyhigh / total_mc
            data_ratio = np.nan_to_num(data_ratio, nan=-999.)
            data_ratio_eylow = np.nan_to_num(data_ratio_eylow, nan=0.)
            data_ratio_eyhigh = np.nan_to_num(data_ratio_eyhigh, nan=0.)
        
    else:
        vardf_data = None
        total_data = None
        print("No data data provided")



    density_factor = 1.0
    # if density is True, area normalize to the data
    if mc_df is not None and data_df is not None and density == True:
        # total_mc already includes stacked Intime Cosmics (+ dirt if present).
        mc_area = np.sum(total_mc)

        data_area = np.sum(total_data)
        density_factor = data_area / mc_area

        weights_categ = [np.array(w) * density_factor for w in weights_categ]

    # Background-only MC spectrum, needed for the optional background systematic
    # band (show_bkgd_syst_band). Requires event-level ``mc.*`` truth — skip for
    # track-level pdg plots (concat'd trk1/trk2 have no event truth block).
    if mc_df is not None and total_mc is not None and breakdown_type != "pdg":
        hist_signal = _overlay_signal_mc_hist(mc_df, var_config, signal_truth_fv=signal_truth_fv)
        if density:
            hist_signal = hist_signal * density_factor
        total_mc_bkgd = np.asarray(total_mc, dtype=float) - hist_signal

    # the order of cuts from get_*_category is reversed from the order of labels and colors
    colors, labels = colors[::-1], labels[::-1]

    # ========================================================

    # ==== plot template ====
    plot_labels = list(plot_labels)
    if len(plot_labels) > 1:
        plot_labels[1] = _overlay_rate_ylabel(plot_labels[1])

    mc_stat_err_ratio = None
    if ratio:
        fig, axs = plt.subplots(2, 1, figsize=(8.5, 8.5), 
                               sharex=True, gridspec_kw={'height_ratios': [4, 1]})
        ax, ax_r = axs[0], axs[1]
        fig.subplots_adjust(hspace=0.1)
        ax_r.axhline(1.0, color='red', linestyle='--', linewidth=1)
        ax_r.set_xlabel(plot_labels[0], fontsize=20)
        ax_r.set_ylabel("Data/MC", fontsize=20)
        ax_r.grid(True)
        ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)
        ax_r.minorticks_on()
        ax_r.tick_params(axis='both', which='major', labelsize=15)
        ax_r.tick_params(axis='both', which='minor', labelsize=13)

    else:
        fig, ax = plt.subplots(figsize=(8.5, 7))
        ax.set_xlabel(plot_labels[0], fontsize=20)

    # common formatting
    ax.set_xlim(var_config.bins[0], var_config.bins[-1])
    ax.set_ylabel(plot_labels[1], fontsize=20)
    ax.set_title(plot_labels[2], fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.tick_params(axis='both', which='minor', labelsize=13)

    # ==== Plot histograms ====

    # == rate panel ==
    # MC
    if var_categ is not None:
        mc_stack, _, _ = ax.hist(var_categ,
                                 weights=weights_categ,
                                 bins=var_config.bins,
                                 stacked=True,
                                 color=colors,
                                 label=labels,
                                 linewidth=0,
                                 edgecolor='none',
                                 histtype='stepfilled')

        breakdown_accum = [np.sum(this_mode) for this_mode in mc_stack]
        breakdown_fractions = [breakdown_accum[0]] + [(breakdown_accum[i+1] - breakdown_accum[i]) for i in range(len(breakdown_accum) - 1)]
        breakdown_fractions = [frac / breakdown_accum[-1] for frac in breakdown_fractions]

        if breakdown_type == "genie_sb":
            # hatch background portion
            bottom = np.zeros(len(var_config.bins) - 1)
            for i, (v, w, h) in enumerate(zip(var_categ, weights_categ, hatches)):
                hist_vals, _ = np.histogram(v, weights=w, bins=var_config.bins)
                ax.bar(
                    var_config.bin_centers,
                    hist_vals,
                    width=np.diff(var_config.bins),
                    bottom=bottom,
                    color='none',
                    hatch=h,
                    edgecolor='white',
                    linewidth=0.0,
                    align='center'
                )
                bottom += hist_vals

    chi2_val = None
    chi2_reduced = None
    p_val = None
    ndof = None
    chi2_pull = None
    chi2_shape_val = None
    p_val_shape = None
    ndof_shape = None
    bkgd_syst_err = None
    syst_err = None

    syst_explicit = syst is not None
    syst = _resolve_overlay_syst_cov_frac(
        var_config,
        syst,
        syst_kind=syst_kind,
        syst_disk_root=syst_disk_root,
        category_syst_summary_path=category_syst_summary_path,
        load_syst_from_summary=load_syst_from_summary,
    )

    if syst is not None: # list of syst uncertainties 

        syst_err = np.sqrt(np.diag(syst)) * total_mc 

        if syst_decomp:
            # The shape/mixed/norm decomposition bands (Matrix_Decomp) were removed
            # as dead code; see git history if this display is ever needed again.
            raise NotImplementedError(
                "syst_decomp=True is no longer supported in overlay_hists"
            )

        ax.bar(
            var_config.bin_centers,
            2 * syst_err,
            width=np.diff(var_config.bins),
            bottom=total_mc - syst_err,
            facecolor='none',             # transparent fill
            hatch='xxx',                 # hatch pattern similar to ROOT's 3004
            linewidth=0.0,
            edgecolor='dimgray',            # outline color of the hatching
            label='Syst. Unc.'
        )

        if show_bkgd_syst_band and total_mc_bkgd is not None:
            bkgd_frac = bkgd_syst_frac_cov
            if bkgd_frac is None:
                bkgd_frac = load_genie_sb_bkgd_rate_cov_frac(
                    var_config, genie_sb_cov_mat_pkl
                )
            bkgd_syst_err = _overlay_bkgd_syst_sigma(total_mc_bkgd, bkgd_frac)
            _overlay_draw_bkgd_syst_band(
                ax,
                var_config.bin_centers,
                var_config.bins,
                total_mc,
                bkgd_syst_err,
                edgecolor=bkgd_syst_band_color,
                label="Bkgd. GENIE unc.",
            )

        if data_df is not None:
            chi2_val, chi2_reduced, p_val, ndof, chi2_pull, chi2_shape_val, p_val_shape, ndof_shape = _overlay_compute_chi2(
                total_data,
                total_mc,
                syst,
                data_eylow,
                data_eyhigh,
                mc_stat_err=mc_stat_err,
                drop_first_bin=_overlay_is_nu_score(var_config),
            )
 
    else:
        print("no syst provided")

    # Data
    if vardf_data is not None:
        ax.errorbar(var_config.bin_centers, 
                    total_data, 
                    yerr=np.vstack((data_eylow, data_eyhigh)),
                    color='black', 
                    fmt='o', markersize=5, capsize=3, linewidth=1.5,
                    label='Data')

    # == ratio panel ==
    if ratio:
        # MC 
        if syst is not None:
            if syst_decomp == False:
                mc_content_ratio = total_mc / total_mc # dummy
                mc_stat_err_ratio = syst_err / total_mc
                mc_stat_err_ratio = np.nan_to_num(mc_stat_err_ratio, nan=0.)
                ax_r.bar(
                    var_config.bin_centers,
                    2*mc_stat_err_ratio,
                    width=np.diff(var_config.bins),
                    bottom=mc_content_ratio - mc_stat_err_ratio,
                    facecolor='none',
                    edgecolor='dimgray',
                    hatch='xxx',
                    linewidth=0.0,
                    label='Syst. Unc.'
                )
                if bkgd_syst_err is not None:
                    mc_content_ratio = np.ones_like(total_mc)
                    bkgd_err_ratio = np.where(
                        total_mc != 0, bkgd_syst_err / total_mc, 0.0
                    )
                    bkgd_err_ratio = np.nan_to_num(bkgd_err_ratio, nan=0.0)
                    ax_r.bar(
                        var_config.bin_centers,
                        2 * bkgd_err_ratio,
                        width=np.diff(var_config.bins),
                        bottom=mc_content_ratio - bkgd_err_ratio,
                        facecolor="none",
                        edgecolor=bkgd_syst_band_color,
                        hatch="+++",
                        linewidth=0.0,
                        label="Bkgd. GENIE unc.",
                    )

        else:
            pass

        # data/MC 
        if data_df is not None:
            ax_r.errorbar(var_config.bin_centers, data_ratio, 
                            yerr=np.vstack((data_ratio_eylow, data_ratio_eyhigh)),
                            fmt='o', color='black',
                            markersize=5, capsize=3, linewidth=1.5)
        set_ratio_panel_ylim(
            ax_r,
            data_ratio=data_ratio if data_df is not None else None,
            data_ratio_eyhigh=data_ratio_eyhigh if data_df is not None else None,
            data_ratio_eylow=data_ratio_eylow if data_df is not None else None,
            syst_err_ratio=mc_stat_err_ratio,
            pad=1.2,
        )

    # nu_score: drop Pandora failure bin (index 0); clip x to bins with data.
    _nu_xlim = None
    if data_df is not None and total_data is not None and _overlay_is_nu_score(var_config):
        _nu_xlim = _overlay_data_xlim_range(
            var_config.bins, total_data, drop_first_bin=True
        )
        if _nu_xlim is not None:
            xmin, xmax, _i_lo, _i_hi = _nu_xlim
            ax.set_xlim(xmin, xmax)
            if ratio and ax_r is not None:
                ax_r.set_xlim(xmin, xmax)

    # ===============================

    # ==== Legend ====
    # legend order: data, mc, syst
    handles, labels_orig = ax.get_legend_handles_labels()
    ordered_handles = []
    ordered_labels = []

    if data_df is not None:
        data_handle_index = labels_orig.index('Data')
        data_handle = handles[data_handle_index]
        ordered_handles.extend([data_handle])
        data_text = 'Data ({:.0f})'.format(sum_data)
        ordered_labels.extend([data_text])

    if mc_df is not None:
        if breakdown_type == "genie_sb":
            # collapse S and B into a single combined legend
            for i in range(len(genie_mode_colors)):
                ordered_handles.append(Patch(facecolor=genie_mode_colors[i], edgecolor='none'))

            legend_labels = []
            i, n = 0, len(labels)
            while i < n:
                # assume paired S/B in the order of MC legend entries
                label_base = labels[i]
                if (i + 1 < n) and (labels[i + 1] == label_base): # paired S/B
                    frac1 = breakdown_fractions[i]
                    frac2 = breakdown_fractions[i+1]
                    legend_label = f"{label_base} ({frac2*100:.1f}%/{frac1*100:.1f}%)"
                    i += 2
                else: # single
                    frac1 = breakdown_fractions[i]
                    legend_label = f"{label_base} ({frac1*100:.1f}%)"
                    i += 1
                legend_labels.append(legend_label)
            ordered_labels.extend(legend_labels[::-1])

        else:
            if data_df is not None:
                mc_handles = [h for i, h in enumerate(handles) if i != data_handle_index and 'Unc.' not in labels_orig[i]]
            else:
                mc_handles = [h for i, h in enumerate(handles) if 'Unc.' not in labels_orig[i]]

            fixed_pcts = _OVERLAY_LEGEND_PCTS.get(breakdown_type)
            if fixed_pcts is not None and len(labels) == len(fixed_pcts):
                # labels are cosmic-first; legend reverses → signal-first display order.
                pcts = list(fixed_pcts)[::-1]
                mc_labels = [
                    f"{label} ({pct:.1f}%)" for label, pct in zip(labels, pcts)
                ]
            else:
                mc_labels = [f"{label} ({frac*100:.1f}%)"
                                    for label, frac in zip(labels, breakdown_fractions)]
            ordered_handles.extend(mc_handles)
            ordered_labels.extend(mc_labels[::-1]) # note the reverse order of mc_labels


    if syst is not None:
        unc_handle = [h for i, h in enumerate(handles) if 'Unc.' in labels_orig[i]]
        unc_label = [l for l in labels_orig if 'Unc.' in l]
        ordered_handles.extend(unc_handle)
        ordered_labels.extend(unc_label)

    # adjust fontsize so that legend fits in the figure
    _ylo, _yhi = (None, None) if _nu_xlim is None else (_nu_xlim[2], _nu_xlim[3])
    ymax_content, height_profile = _overlay_content_ymax_profile(
        total_mc, syst_err, total_data, data_eyhigh, _ylo, _yhi,
    )
    hp_shape = height_profile
    if height_profile is not None and _ylo is not None:
        hp_shape = np.asarray(height_profile, dtype=float)[_ylo : _yhi + 1]
    var_save = getattr(var_config, "var_save_name", None) if var_config is not None else None

    fontsize = 11
    chi2_fontsize = 16
    prefer_left_chi2 = None
    ncol = 3
    style = None
    reserved_boxes = []
    if breakdown_type == "genie_sb":
        textloc_x, textloc_ha = get_textloc_x(total_mc, var_config.bins, textloc)
        fontsize = fontsize
        ncol = 2

        # separate legend box with S / B hatches
        example_signal = Patch(facecolor="black", edgecolor='white', label='Signal')
        example_background = Patch(facecolor="black", edgecolor='white', hatch='////', linewidth=0, label='Background')
        # if textloc_x < 0.5:
        #     box_ax = ax.inset_axes([0.17, 0.65, 0.13, 0.13], transform=ax.transAxes)
        # else:
        box_ax = ax.inset_axes([0.625, 0.66, 0.13, 0.13], transform=ax.transAxes)
        box_ax.axis('off')
        mini_legend = Legend(
            box_ax,
            handles=[example_signal, example_background],
            labels=['Signal', 'Background'],
            loc='center',
            fontsize=fontsize,
            frameon=False,
            borderpad=0.7,
            handlelength=2.1,
            handleheight=0.9,
            ncol=1,
            fancybox=True,
            framealpha=1.0
        )
        box_ax.add_artist(mini_legend)

        ax.legend(
            ordered_handles,
            ordered_labels,
            loc='upper left',
            # loc='upper center',
            fontsize=fontsize,
            frameon=False,
            ncol=ncol,
            bbox_to_anchor=(0.05, 0.9, 0.8, 0.1),
            mode='expand'
        )

    else:
        style = _overlay_stack_legend_style(breakdown_type, hp_shape, var_save=var_save)
        if style is not None:
            fontsize = style["fontsize"]
            chi2_fontsize = style["chi2_fontsize"]
            prefer_left_chi2 = style["prefer_left"]
            ax.legend(
                ordered_handles,
                ordered_labels,
                loc=style["loc"],
                fontsize=fontsize,
                frameon=False,
                ncol=style["ncol"],
                bbox_to_anchor=style["bbox_to_anchor"],
                borderaxespad=0.2,
                handletextpad=0.4,
                labelspacing=0.22,
                columnspacing=0.9,
                handlelength=1.6,
            )
            if ymax_content > 0:
                ax.set_ylim(0.0, float(style["ylim_scale"]) * ymax_content)
            reserved_boxes = _overlay_reserved_axes_boxes(ax)
            if style.get("tight_fit"):
                _overlay_fit_ylim(
                    ax,
                    height_profile,
                    var_config.bins,
                    reserved_boxes,
                    ymax_content,
                    min_scale=float(style.get("tight_min_scale") or 1.08),
                    pad=float(style.get("tight_pad") or 0.05),
                )
            else:
                _overlay_raise_ylim_for_legend(
                    ax,
                    height_profile,
                    var_config.bins,
                    list(reserved_boxes) + list(style.get("clear_boxes") or []),
                )
            reserved_boxes = _overlay_reserved_axes_boxes(ax)
        else:
            ax.legend(
                ordered_handles,
                ordered_labels,
                loc='upper left',
                # loc='upper center',
                fontsize=fontsize,
                frameon=False,
                ncol=ncol,
            )

    # ax_r.legend(fontsize=9, ncol=2)

    # ===============================

    # ==== plot additions ====

    # y-axis: keep stack + syst hatch (+ data) below the upper-left legend
    if style is None:
        if ymax_content > 0:
            reserved_boxes = _overlay_apply_tight_ylim(
                ax,
                ymax_content,
                max_headroom=(
                    max(float(ax_ylim_ratio), 1.65)
                    if breakdown_type == "pdg"
                    else float(ax_ylim_ratio)
                ),
            )
        elif total_mc is not None:
            ax.set_ylim(0.0, ax_ylim_ratio * float(np.nanmax(total_mc)))
    # ax.set_yscale("log")

    # vertical lines
    if vline is not None:
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v[0], ymin=0, ymax=ymax*0.75, color='red', linestyle='--')
            # direction: 0 = keep left (< cut), 1 = keep right (> cut)
            if len(v) > 1 and v[1] is not None:
                direction = v[1]
                xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
                yspan = ax.get_ylim()[1] - ax.get_ylim()[0]
                arrow_params = {
                    'y': ymax * (0.22 if style is not None else 0.4),
                    'dx': 0.04 * xspan,
                    'width': 0.01 * yspan,
                    'color': 'red',
                    'head_width': 0.04 * yspan,
                    'head_length': 0.02 * xspan,
                    'length_includes_head': True
                }
                if direction == 0:
                    # Left arrow
                    ax.arrow(v[0], arrow_params['y'], -arrow_params['dx'], 0, 
                             width=arrow_params['width'],
                             color=arrow_params['color'],
                             head_width=arrow_params['head_width'],
                             head_length=arrow_params['head_length'],
                             length_includes_head=arrow_params['length_includes_head'],
                             clip_on=True)
                elif direction == 1:
                    # Right arrow
                    ax.arrow(v[0], arrow_params['y'], arrow_params['dx'], 0, 
                             width=arrow_params['width'],
                             color=arrow_params['color'],
                             head_width=arrow_params['head_width'],
                             head_length=arrow_params['head_length'],
                             length_includes_head=arrow_params['length_includes_head'],
                             clip_on=True)

    # textboxes — χ² above local content, clear of the legend
    _overlay_place_chi2_genie(
        ax,
        style=style,
        reserved_boxes=reserved_boxes,
        textloc=textloc,
        height_profile=height_profile,
        bins=var_config.bins,
        ymax_content=ymax_content,
        breakdown_type=breakdown_type,
        var_save=var_save,
        textchi2=textchi2,
        chi2_val=chi2_val,
        p_val=p_val,
        ndof=ndof,
        chi2_shape_val=chi2_shape_val,
        p_val_shape=p_val_shape,
        ndof_shape=ndof_shape,
        fontsize=fontsize,
        chi2_fontsize=chi2_fontsize,
        prefer_left_chi2=prefer_left_chi2,
    )

    fig.subplots_adjust(top=0.9)
    add_approval_text(approval, 0.03, 1.07, "left")
    if pot_text:
        add_pot_text(format_pot_corner_text(pot_text), 0.99, 1.01, "right", fontsize=16)

    if var_config.var_save_name == "integrated":
        format_singlebin_plot()


    # ===============================

    # == save figure ==
    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches="tight", dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()

    return {"breakdown_type": breakdown_type,
            "var_name": var_config.var_save_name,
            "bins": var_config.bins,
            # "cuts": cuts, 
            "mc_stack": mc_stack,
            "total_mc": total_mc, 
            "total_mc_bkgd": total_mc_bkgd,
            "total_data": total_data,
            "chi2_val": chi2_val,
            "p_val": p_val,
            "ndof": ndof,
            "chi2_pull": chi2_pull,
            "chi2_shape_val": chi2_shape_val,
            "p_val_shape": p_val_shape,
            "ndof_shape": ndof_shape}



# ==== specialized 1D plots: pulls, efficiency, universes ====
def plot_chi2_pull(chi2_pull,
                   var_config,
                   chi2_val=None,
                   p_val=None,
                   ndof=None,
                   approval="internal",
                   textloc=[0.03, 0.92],
                   plot=True,
                   save_fig=False,
                   save_name=None):
    """
    Bar chart of per-bin chi2 pull: (Data - MC) / sqrt(diag(combined_cov)).
    Positive pulls are blue, negative pulls are red.
    Dashed lines mark ±1σ and ±2σ for visual reference.
    """
    bin_centers = var_config.bin_centers
    bins = var_config.bins

    colors = np.where(chi2_pull >= 0, "steelblue", "tomato")

    fig, ax = plt.subplots(figsize=(8.5, 4))
    ax.bar(bin_centers, chi2_pull,
           width=np.diff(bins),
           color=colors,
           edgecolor="none",
           linewidth=0)
    ax.axhline(0,  color="black", linewidth=1.2)
    ax.axhline( 1, color="gray",  linewidth=0.9, linestyle="--", alpha=0.8, label=r"$\pm 1\sigma$")
    ax.axhline(-1, color="gray",  linewidth=0.9, linestyle="--", alpha=0.8)
    ax.axhline( 2, color="gray",  linewidth=0.6, linestyle=":",  alpha=0.6, label=r"$\pm 2\sigma$")
    ax.axhline(-2, color="gray",  linewidth=0.6, linestyle=":",  alpha=0.6)

    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel(var_config.var_labels[1])
    ax.set_ylabel(r"Pull = (Data $-$ MC) / $\sigma$")
    ax.legend(fontsize=10, frameon=False, loc="lower right")
    ax.grid(which="major", linestyle="-", linewidth=0.5, alpha=0.5)
    ax.grid(which="minor", linestyle=":", linewidth=0.4, alpha=0.4)
    ax.minorticks_on()

    if chi2_val is not None:
        textloc_x = textloc[0]
        textloc_y = textloc[1]
        textloc_ha = "left"
        ax.text(textloc_x, textloc_y,
                f"$\\chi^2$/ndof = {chi2_val:.2f}/{ndof} (p-value = {p_val:.2f})",
                transform=ax.transAxes,
                ha=textloc_ha, va="top",
                fontsize=12, color="black")

    add_approval_text(approval, 0.97, 0.97, "right")

    if var_config.var_save_name == "integrated":
        ax.set_xticks([])

    if save_fig:
        plt.savefig(save_name + fig_ext, bbox_inches="tight", dpi=dpi)

    if plot:
        plt.show()
    else:
        plt.close()


def plot_efficiency(df_dict={}, 
                    stage_labels=[],
                    var_config=None, 
                    textloc=[0.05, 0.55],
                    approval="internal", 
                    legend=True,
                    plot=True,
                    save_fig=False, 
                    save_name=None):

    var_name = var_config.var_nu_col
    bins = var_config.bins
    bin_centers = var_config.bin_centers

    keys = list(df_dict.keys())
    colors = ["C" + str(i) for i in range(len(keys))]

    all_signal_df = df_dict[keys[0]][IsNuInFV_NumuCC_1p0pi(df_dict[keys[0]])]
    n_tot_integ = len(all_signal_df)
    var, wgts = get_clipped_evts(all_signal_df, var_name, bins)
    n_tot, _ = np.histogram(var, weights=wgts, bins=bins)

    fig, ax = plt.subplots()
    ax_eff = ax.twinx()

    eff_list = []
    eff_err_list = []
    for kidx, key in enumerate(keys):
        this_df = df_dict[key]
        this_signal_df = this_df[IsNuInFV_NumuCC_1p0pi(this_df)]
        var, wgts = get_clipped_evts(this_signal_df, var_name, bins)
        n, _ = np.histogram(var, weights=wgts, bins=bins)

        this_eff_integ = len(this_signal_df) / n_tot_integ
        this_label = stage_labels[kidx] + " ({:.2f}%)".format(this_eff_integ*100)
        n, bins, _ = ax.hist(var, weights=wgts, bins=bins, histtype="step", label=this_label, alpha=0.6)

        if key == keys[-1]:
            print("final purity: {:.2f}%".format(100 * len(this_signal_df) / len(this_df)))

        this_eff = n / n_tot
        stat_scale_factor = this_df.pot_weight.unique()[0]
        this_eff_err = get_eff_err(n/stat_scale_factor, n_tot/stat_scale_factor)
        ax_eff.errorbar(bin_centers, this_eff, yerr=this_eff_err, fmt="o-", color=colors[kidx], label=this_label)
        eff_list.append(this_eff)
        eff_err_list.append(this_eff_err)

    ax.set_xlabel(var_config.var_labels[0])
    ax.set_ylabel("Events")
    ax_eff.set_ylabel("Efficiency")
    ax_eff.set_ylim(0, 1.05)
    plt.xlim(bins[0], bins[-1])

    if legend:
        plt.legend(bbox_to_anchor=(1.15, 1.0))

    # ==== plot additions ====
    textloc_x, textloc_ha = get_textloc_x(n_tot, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    if var_config.var_save_name == "integrated":
        format_singlebin_plot()

    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches="tight", dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()

    return {"eff_list": eff_list, "eff_err_list": eff_err_list}


def plot_univ_hists(
                univ_events, 
                cv_events,
                syst_name, 
                var_config, 
                approval="internal",
                textloc=[0.05, 0.55],
                plot=True,
                ax_titles = ["", "", ""],
                save_fig=False, 
                save_name=None): 

    assert univ_events.shape[1] == len(cv_events) 
    n_univ = univ_events.shape[0]

    if (n_univ > 10):
        colors = ["#FDE725FF", "#1F968BFF", "#440154FF"] # viridis colors

        sorted_univs = np.sort(univ_events, axis=0)

        n_68 = int(0.68 * n_univ)
        start_68 = (n_univ - n_68) // 2
        end_68 = start_68 + n_68

        n_95 = int(0.95 * n_univ)
        start_95 = (n_univ - n_95) // 2
        end_95 = start_95 + n_95

        # Define bins & groupings of universes: [(range, color, label, skip_68)]
        segs = [
            (range(start_68, end_68), colors[0], "Universe (68%)", False),
            (range(start_95, end_95), colors[1], "Universe (95%)", True),
            ( (i for i in range(n_univ) if i not in range(start_95, end_95)), colors[2], "Universe (100%)", False)
        ]

        plotted = set()  # tracks which labels were plotted already
        for r, color, label, skip_68 in segs:
            for i in r:
                if skip_68 and i in range(start_68, end_68): 
                    continue
                show_label = label if label not in plotted else None
                plt.hist(var_config.bin_centers, bins=var_config.bins, weights=sorted_univs[i],
                        histtype="step", color=color, alpha=0.7, label=show_label)
                plotted.add(label)

    # if too few universes, just plot all of them in gray
    else:
        for i in range(n_univ):
            show_label = "Universe" if i == 0 else None
            plt.hist(var_config.bin_centers, bins=var_config.bins, weights=univ_events[i], histtype="step", color="gray", label=show_label)

    # plot CV last so that it is on top of the universes
    plt.hist(var_config.bin_centers, bins=var_config.bins, weights=cv_events, histtype="step", color="k", label="Central Value")

    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(var_config.var_labels[1])
    if ax_titles[1] != "":
        plt.ylabel(ax_titles[1])
    else:
        plt.ylabel("Events")
    plt.title(ax_titles[2])

    plt.legend(frameon=False)

    # ==== plot additions ====
    textloc_x, textloc_ha = get_textloc_x(cv_events, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    if var_config.var_save_name == "integrated":
        format_singlebin_plot()


    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()


# ==== unfolded-result display ====
class HandlerDoubleErrorbar(HandlerErrorbar):
    """Legend marker with nested y-error bars (inner shape, outer shape ⊕ stat)."""

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        outer = orig_handle[0] if isinstance(orig_handle, tuple) else orig_handle
        artists = super().create_artists(
            legend, outer, xdescent, ydescent, width, height, fontsize, trans
        )
        if not getattr(outer, "has_yerr", False):
            return artists
        xdata, xdata_marker = self.get_xdata(
            legend, xdescent, ydescent, width, height, fontsize
        )
        xdata_marker = np.asarray(xdata_marker)
        ydata = np.full_like(xdata, (height - ydescent) / 2.0)
        ydata_marker = np.asarray(ydata[: len(xdata_marker)])
        _, yerr_outer = self.get_err_size(
            legend, xdescent, ydescent, width, height, fontsize
        )
        yerr_inner = 0.55 * yerr_outer
        verts = [
            ((x, y - yerr_inner), (x, y + yerr_inner))
            for x, y in zip(xdata_marker, ydata_marker)
        ]
        plotlines, caplines, barlinecols = outer
        coll = mcoll.LineCollection(verts)
        self.update_prop(coll, barlinecols[0], legend)
        coll.set_transform(trans)
        extra = [coll]
        if caplines:
            cap_lo = Line2D(xdata_marker, ydata_marker - yerr_inner)
            cap_hi = Line2D(xdata_marker, ydata_marker + yerr_inner)
            self.update_prop(cap_lo, caplines[0], legend)
            self.update_prop(cap_hi, caplines[0], legend)
            cap_lo.set_marker("_")
            cap_hi.set_marker("_")
            cap_lo.set_transform(trans)
            cap_hi.set_transform(trans)
            extra.extend([cap_lo, cap_hi])
        # Outer bars first, then inner, then line/marker from the parent handler.
        return [artists[0], *extra, *artists[1:]]


def _covariance_per_bin_width(cov, bin_widths):
    """If ``x_i = y_i / bw_i``, propagate ``Cov(y)`` to ``Cov(x)`` with ``D = diag(1/bw)``."""
    invbw = 1.0 / np.clip(np.asarray(bin_widths, dtype=float), 1e-300, None)
    d = np.diag(invbw)
    cov = np.asarray(cov, dtype=float)
    return d @ cov @ d


def plot_unfolded_result(unfold, 
                         measured, 
                         models,
                         var_config, 
                         chi2_list=[],
                         chi2_dict={},
                         xsec_unit=0,
                         textloc=[0.05, 0.55],
                         approval="internal",
                         plot_labels=["", "", ""],
                         plot=True,
                         save_fig=False, 
                         save_name=None,
                         save_ext=None,
                         data=False,
                         closure_test=False,
                         model_add_smear=None,
                         pot_text="8.8 $\\times 10^{19}$ POT"):

    bins = var_config.bins
    bin_centers = var_config.bin_centers
    bin_widths = np.diff(bins)
    if len(var_config.bins) == 2:
        bin_widths = np.array([1.0])
    ndof_bins = int(max(len(bins) - 1, 1))

    # Full unfolded covariance (stat + syst), scaled to the same per-bin-width units as the plot/chi2 vectors.
    # Older code used only unfold['SystUnfoldCov'] and omitted per-width scaling → chi2 was vastly inflated.
    cov_unfold_perwidth = _covariance_per_bin_width(unfold["UnfoldCov"], bin_widths)

    # unfolded result
    Unfolded = unfold['unfold']
    UnfoldedCov = unfold["UnfoldCov"]
    Unfolded_perwidth = Unfolded / bin_widths

    # --- stat uncertainties
    UnfoldCov_stat = unfold['StatUnfoldCov']
    Unfold_uncert_stat = np.sqrt(np.maximum(np.diag(UnfoldCov_stat), 0.0))

    # --- syst uncertainties
    UnfoldCov_syst = unfold['SystUnfoldCov']
    Unfold_uncert_syst = np.diag(UnfoldCov_syst)
    UnfoldCov_syst_frac = fraccov_from_cov(UnfoldCov_syst, Unfolded)

    # --- decompose into norm and shape components
    # the first item in models dict is the nominal input model
    norm_model = list(models.keys())[0]
    SystUnfoldCov_norm, SystUnfoldCov_mixed, SystUnfoldCov_shape = Matrix_Decomp(models[norm_model][0], UnfoldCov_syst)
    Unfold_uncert_norm = np.sqrt(np.abs(np.diag(SystUnfoldCov_norm)))
    Unfold_uncert_shape = np.sqrt(np.abs(np.diag(SystUnfoldCov_shape)))


    # --- plot
    fig, ax = plt.subplots(figsize=(8, 6))
    shape_handle = None
    # set err to 0 for closure test
    if closure_test:
        dummy_err = np.zeros_like(Unfolded_perwidth)
        bar_handle = plt.errorbar(bin_centers, Unfolded_perwidth, yerr=dummy_err, fmt='o', color='black')

    else:
        # Draw nested black error bars:
        # inner  = shape syst only
        # outer  = shape syst + stat (in quadrature)
        Unfold_uncert_stat_perwidth = Unfold_uncert_stat / bin_widths
        Unfold_uncert_shape_perwidth = Unfold_uncert_shape / bin_widths
        Unfold_uncert_total_perwidth = np.sqrt(
            np.maximum(
                Unfold_uncert_shape_perwidth**2 + Unfold_uncert_stat_perwidth**2,
                0.0,
            )
        )
        bar_handle = plt.errorbar(
            bin_centers,
            Unfolded_perwidth,
            yerr=Unfold_uncert_total_perwidth,
            fmt='o',
            color='black',
            ecolor='black',
            elinewidth=1.5,
            capsize=3,
        )
        shape_handle = plt.errorbar(
            bin_centers,
            Unfolded_perwidth,
            yerr=Unfold_uncert_shape_perwidth,
            fmt='none',
            ecolor='black',
            elinewidth=1.5,
            capsize=3,
        )

        # plot syst norm component as histogram at the bottom
        Unfold_uncert_norm_perwidth = Unfold_uncert_norm / bin_widths
        if len(var_config.bins) != 2:
            norm_handle = plt.bar(bin_centers, Unfold_uncert_norm_perwidth, width=bin_widths, label='Syst. error (norm)', alpha=0.5, color='gray')

    if data:
        # Keep uncertainties from unfold covariance components only; this avoids NaNs from
        # sqrt(measured/xsec_unit) when background-subtracted measured bins go negative.
        pass

    # divide measured & model by bin width
    measured_perwidth = measured / bin_widths
    reco_handle = None
    if not data:
        # Measured Signal (Input) — restore when you want the folded fake-data / Asimov input on the plot:
        # reco_handle, = plt.step(
        #     bins,
        #     np.append(measured_perwidth, measured_perwidth[-1]),
        #     where="post",
        #     label="Measured Signal (Input)",
        # )
        pass

    # --- get chi2 values for each model to compare
    # ``chi2_list`` may be a list (one value per model) or a dict keyed by
    # ``var_config.var_save_name`` (single value applied to GENIE / sole model).
    use_provided_chi2 = False
    chi2_vals = []
    p_values = []
    ndof_list = []
    if isinstance(chi2_list, dict):
        if var_config.var_save_name in chi2_list:
            chi2_vals = [float(chi2_list[var_config.var_save_name])]
            ndof_list = [ndof_bins]
            use_provided_chi2 = True
    elif len(chi2_list) > 0:
        chi2_vals = [float(v) for v in chi2_list]
        ndof_list = [ndof_bins] * len(chi2_vals)
        use_provided_chi2 = True
    model_handles = []
    model_labels = []
    for midx, mkey in enumerate(models.keys()):
        add_smear = unfold["AddSmear"]
        if model_add_smear is not None and mkey in model_add_smear:
            add_smear = model_add_smear[mkey]
        model_smeared = add_smear @ models[mkey][0]
        # if "SBN" in mkey:
        model_smeared_perwidth = model_smeared / bin_widths

        # else:
        #     model_smeared_perwidth = model_smeared 

        if not use_provided_chi2:
            # remove bins with <= 0 events
            mask = (Unfolded_perwidth > 0) & (model_smeared_perwidth > 0)
            Unfolded_perwidth_safe = Unfolded_perwidth[mask]
            model_smeared_perwidth_safe = model_smeared_perwidth[mask]
            cov_chi2_safe = cov_unfold_perwidth[np.ix_(mask, mask)]
            chi2_val, p_val = get_chi2(
                Unfolded_perwidth_safe, model_smeared_perwidth_safe, cov_chi2_safe
            )
            chi2_vals.append(chi2_val)
            p_values.append(p_val)
            ndof_list.append(int(np.sum(mask)))

        model_handle, = plt.step(bins, np.append(model_smeared_perwidth, model_smeared_perwidth[-1]), where='post', color=models[mkey][1])
        model_handles.append(model_handle)
        # model_labels.append(f'$A_c \\otimes$ {mkey}')
        # Optional display name: models[key] = [spectrum, color, legend_label]
        mval = models[mkey]
        legend_label = mval[2] if isinstance(mval, (list, tuple)) and len(mval) >= 3 and mval[2] else mkey
        model_labels.append(f'{legend_label}')

    # Nested y-error bars on the Data marker: inner = shape, outer = shape ⊕ stat.
    data_eb_handle = (bar_handle, shape_handle) if shape_handle is not None else bar_handle
    data_label = r"Data (Shape Syst.$\oplus$Stat.)"
    fake_data_label = r"Fake Data (Shape Syst.$\oplus$Stat.)"
    norm_label = "Norm. Syst."

    # legend
    if closure_test:
        if reco_handle is not None:
            handles = [bar_handle, reco_handle] + model_handles
            labels = ["Unfolded Asimov Data", "Measured Signal"] + model_labels
        else:
            handles = [bar_handle] + model_handles
            labels = ["Unfolded Asimov Data"] + model_labels
    elif data:
        if len(var_config.bins) == 2:
            handles = [data_eb_handle] + model_handles
            labels = [r"Data (Syst.$\oplus$Stat.)"] + model_labels
        else:
            handles = [data_eb_handle, norm_handle] + model_handles
            labels = [data_label, norm_label] + model_labels
    else:
        if len(var_config.bins) == 2:
            if reco_handle is not None:
                handles = [bar_handle, reco_handle] + model_handles
                labels = ["Fake Data", "Measured Signal (Input)"] + model_labels
            else:
                handles = [bar_handle] + model_handles
                labels = ["Fake Data"] + model_labels
        else:
            if reco_handle is not None:
                handles = [data_eb_handle, shape_handle, norm_handle, reco_handle] + model_handles
                labels = [
                    fake_data_label,
                    "Shape Syst.",
                    norm_label,
                    "Measured Signal (Input)",
                ] + model_labels
            else:
                handles = [data_eb_handle, shape_handle, norm_handle] + model_handles
                labels = [fake_data_label, "Shape Syst.", norm_label] + model_labels
    # Append chi2/ndof to model legend entries
    n_non_model = len(labels) - len(model_labels)
    if len(chi2_vals) == len(models):
        ndofs = ndof_list if len(ndof_list) == len(models) else [ndof_bins] * len(models)
        for midx in range(len(model_labels)):
            suffix = f" ({float(chi2_vals[midx]):.1f}/{int(ndofs[midx])})"
            labels[n_non_model + midx] += suffix
    elif (
        isinstance(chi2_list, dict)
        and len(chi2_vals) == 1
        and len(models) > 1
    ):
        model_keys = list(models.keys())
        midx = model_keys.index("GENIE") if "GENIE" in model_keys else 0
        suffix = f" ({float(chi2_dict[model_keys[midx]][0]):.1f}/{int(ndof_bins)})"
        labels[n_non_model + midx] += suffix

    legend_handler_map = {}
    if isinstance(data_eb_handle, tuple):
        legend_handler_map[tuple] = HandlerDoubleErrorbar(yerr_size=0.7)
    plt.legend(
        handles, labels,
        loc='best', fontsize=15, frameon=False, ncol=1,
        handleheight=1.6, handlelength=1.0,
        handler_map=legend_handler_map or None,
    )
    # plt.legend(handles, labels,
    #            loc=(0.02, 0.8), fontsize=12, frameon=False, ncol=1)

    plt.xlabel(var_config.var_labels[0], fontsize=22)
    plt.ylabel(var_config.xsec_label, fontsize=22)
    plt.title(plot_labels[2])
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., np.max(Unfolded_perwidth)*1.2)

    # ==== plot additions
    # textloc_x, textloc_ha = get_textloc_x(Unfolded_perwidth, var_config.bins, textloc)
    # _, textloc_ha = get_textloc_x(Unfolded_perwidth, var_config.bins, textloc)
    textloc_x = textloc[0]
    textloc_ha = "left"
    textloc_y = textloc[1]
    fig.subplots_adjust(top=0.9)
    add_approval_text(approval, 0.15, 1.07, "left")
    if pot_text:
        corner = (
            pot_text
            if isinstance(pot_text, str) and "Simulation" in pot_text
            else format_pot_corner_text(pot_text)
        )
        add_pot_text(corner, 0.99, 1.02, "right", fontsize=16)
    # add_genie_version_text(textloc_x, textloc_y-0.08, textloc_ha)

    if var_config.var_save_name == "integrated":
        format_singlebin_plot()

    if save_fig:
        plt.savefig(save_name + (fig_ext if save_ext is None else save_ext), bbox_inches="tight", dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()


# ==== signal histogram builders ====
def signal_cut(df, detector=DETECTOR, signal_truth_fv="per_tpc"):
    return df[IsNuInFV_NumuCC_1p0pi(df, detector=detector, signal_truth_fv=signal_truth_fv)]


def signal_hists(evtdf=None,  # df with selected & reco'ed events
                 nudf=None,   # df with all MC truth
                 var_config=None,
                 return_data=False,
                 plot=True,
                 textloc=[0.05, 0.55],
                 approval="internal",
                 mode="reco",
                 save_fig=False, 
                 save_name=None,
                 signal_truth_fv="per_tpc"):
    """
    generate / selected / reco'ed signal events

    topo_categ == 1 : signal
    topology_list has "1" as the first item, corresponding to the signal topology

    items with "sel" tags are selected events that are signal in truth 
    item with "allsel" tags are all selected events, signal + background
    """

    bins = var_config.bins
    bin_centers = var_config.bin_centers
    reco_col = var_config.var_evt_reco_col
    truth_col = var_config.var_evt_truth_col
    nu_col = var_config.var_nu_col

    # ===== all selected events =====
    # reco'ed
    _vsn = var_config.var_save_name
    var_allsel_reco, wgt_allsel_reco = get_clipped_evts(
        evtdf, reco_col, bins, var_save_name=_vsn
    )
    var_allsel_truth, wgt_allsel_truth = get_clipped_evts(
        evtdf, truth_col, bins, var_save_name=_vsn
    )

    # ===== true signal events =====
    evtdf_signal = evtdf[evtdf.topo_categ == 1]
    # selected, reco'ed
    var_sel_reco, wgt_sel_reco = get_clipped_evts(
        evtdf_signal, reco_col, bins, var_save_name=_vsn
    )
    # selected, truth (for response matrix)
    var_sel_truth, wgt_sel_truth = get_clipped_evts(
        evtdf_signal, truth_col, bins, var_save_name=_vsn
    )

    nevts_allsel_truth, _ = np.histogram(var_allsel_truth, weights=wgt_allsel_truth, bins=bins)
    nevts_allsel_reco, _  = np.histogram(var_allsel_reco,  weights=wgt_allsel_reco,  bins=bins)
    nevts_sel_truth, _    = np.histogram(var_sel_truth,    weights=wgt_sel_truth,    bins=bins)
    nevts_sel_reco, _     = np.histogram(var_sel_reco,     weights=wgt_sel_reco,     bins=bins)

    if nudf is not None:
        # all MC, truth (for efficiency vector)
        nudf_signal = nudf[nudf.topo_categ == 1]
        var_allmc, wgt_allmc = get_clipped_evts(
            nudf_signal, nu_col, bins, var_save_name=_vsn
        )
        nevts_allmc, _ = np.histogram(var_allmc, weights=wgt_allmc, bins=bins)
        if mode == "unfold":
            nudf_signal = signal_cut(nudf, signal_truth_fv=signal_truth_fv)
            var_allmc, wgt_allmc = get_clipped_evts(
                nudf_signal, nu_col, bins, var_save_name=_vsn
            )
            nevts_allmc, _ = np.histogram(var_allmc, weights=wgt_allmc, bins=bins)
    else:
        var_allmc = None
        wgt_allmc = None
        nevts_allmc = None

    if plot:
        plt.hist(bin_centers, weights=nevts_allsel_truth, bins=bins, histtype="step", label="Selected Events, True", color="C0")
        plt.hist(bin_centers, weights=nevts_allsel_reco,  bins=bins, histtype="step", label="Selected Events, Reco", color="C0", linestyle="--")
        plt.hist(bin_centers, weights=nevts_sel_truth,    bins=bins, histtype="step", label="Selected Signal, True", color="C1")
        plt.hist(bin_centers, weights=nevts_sel_reco,     bins=bins, histtype="step", label="Selected Signal, Reco", color="C1", linestyle="--")

        if nevts_allmc is not None:
            nevts_allmc, _, _ = plt.hist(var_allmc,     weights=wgt_allmc,     bins=bins, histtype="step", label="All Signal in MC", color="black")
        else:
            nevts_allmc = None

        plt.xlabel(var_config.var_labels[0])
        plt.ylabel("Events / Bin")
        plt.xlim(bins[0], bins[-1])
        plt.legend()

        # ==== plot additions ====
        textloc_x, textloc_ha = get_textloc_x(nevts_allsel_truth, bins, textloc)
        textloc_y = textloc[1]
        add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

        if var_config.var_save_name == "integrated":
            format_singlebin_plot()

        if save_fig:
            plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

        if plot == True:
            plt.show()
        else:
            plt.close()

    if return_data:
        return {
            "var_allmc": var_allmc,
            "wgt_allmc": wgt_allmc,
            "nevts_allmc": nevts_allmc,

            "var_sel_truth": var_sel_truth,
            "wgt_sel_truth": wgt_sel_truth,
            "nevts_sel_truth": nevts_sel_truth,

            "var_sel_reco": var_sel_reco,
            "wgt_sel_reco": wgt_sel_reco,
            "nevts_sel_reco": nevts_sel_reco,

            "var_allsel_truth": var_allsel_truth,
            "wgt_allsel_truth": wgt_allsel_truth,
            "nevts_allsel_truth": nevts_allsel_truth,

            "var_allsel_reco": var_allsel_reco,
            "wgt_allsel_reco": wgt_allsel_reco,
            "nevts_allsel_reco": nevts_allsel_reco,
        }

# ==== fractional uncertainty plot ====
def plot_frac_unc(frac_unc_list, 
                  var_config, 
                  plot_labels=["", "", ""],
                  legends = None,
                  textloc=[0.05, 0.55],
                  approval="internal",
                  plot=True,
                  save_fig=False, 
                  save_name=None):

    for fidx, frac_unc in enumerate(frac_unc_list):
        color = "C{}".format(fidx)
        if len(frac_unc_list) == 1:
            color = "black"
        plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_unc, histtype="step", color=color)

    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(var_config.var_labels[0])
    plt.ylabel("Fractional Uncertainty")
    plt.title(plot_labels[2])
    plt.grid(True)

    if legends is not None:
        plt.legend(legends)

    textloc_x, textloc_ha = get_textloc_x(frac_unc, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    if var_config.var_save_name == "integrated":
        format_singlebin_plot()

    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()


# Publication-style systematic-uncertainty breakdown used by unfolding diagnostics.
# Labels match the measurement paper order (not the CategorySummary key names).
XSEC_SYST_BREAKDOWN_SPECS = (
    ("flux", "Neutrino flux"),
    ("genie_xsec", "Neutrino interaction"),
    ("g4", "Reinteraction"),
    ("detector", "Detector"),
    ("pot", "Protons-on-target"),
    ("ntargets", "Number of target Ar"),
    ("cosmics", "Cosmics"),
    ("mcstat", "MC statistics"),
)

# Exposure printed on the Section 4 uncertainty breakdown (two sig. digits).
UNFOLD_PLOT_POT = 8.8e19


def plot_xsec_syst_breakdown(
    var_config,
    summary,
    *,
    save_name=None,
    plot=False,
    approval="",
    kind="xsec",
    pot=UNFOLD_PLOT_POT,
):
    """Fractional systematic breakdown vs bin, publication styling.

    Sources: Neutrino flux, Neutrino interaction, Reinteraction, Detector,
    Protons-on-target, Number of target Ar, Cosmics, MC statistics, plus
    Total.
    """
    from analysis_village.numucc_1p0pi.syst_category_summary import (
        frac_weights_for_plot,
        total_cov_frac,
    )

    pack = summary["by_var"][var_config.var_save_name]
    cats = pack["categories"]
    bc = var_config.bin_centers
    bins = var_config.bins

    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    colors = list(plt.cm.tab10.colors)
    plotted = []
    for i, (key, lab) in enumerate(XSEC_SYST_BREAKDOWN_SPECS):
        if key not in cats:
            continue
        w = frac_weights_for_plot(cats[key]["cov_frac"], var_config)
        ax.hist(
            bc, bins=bins, weights=w, histtype="step", linewidth=2.0,
            color=colors[i % len(colors)], label=lab,
        )
        plotted.append(w)

    tot = frac_weights_for_plot(
        total_cov_frac(summary, var_config.var_save_name, kind=kind), var_config,
    )
    ax.hist(
        bc, bins=bins, weights=tot, histtype="step",
        linewidth=2.4, color="black", label="Total",
    )
    plotted.append(tot)

    ymax = max(float(np.nanmax(w)) for w in plotted) if plotted else 1.0
    xlab = var_config.var_labels[0] if getattr(var_config, "var_labels", None) else var_config.var_save_name
    ax.set_xlim(float(bins[0]), float(bins[-1]))
    # Headroom for a 3-row legend that spans the axes width.
    ax.set_ylim(0.0, max(ymax * 1.85, 1.0))
    ax.set_xlabel(xlab, fontsize=22)
    ax.set_ylabel("Uncertainty [%]", fontsize=22)
    ax.minorticks_on()
    ax.legend(
        loc="upper center",
        ncol=3,
        fontsize=13,
        frameon=False,
        mode="expand",
        bbox_to_anchor=(0.0, 0.98, 1.0, 0.0),
        borderaxespad=0.0,
        handlelength=1.6,
        columnspacing=0.8,
    )
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    add_approval_text(approval, 0.03, 1.08, "left")
    if pot:
        add_pot_text(format_pot_corner_text(pot), 0.99, 1.02, "right", fontsize=16)

    if save_name is not None:
        plt.savefig(save_name + ".pdf", bbox_inches="tight")
    if plot:
        plt.show()
    else:
        plt.close()


def _true_axis_to_reg(label: str) -> str:
    """Swap a ``^{true}`` axis superscript for ``^{reg}``.

    ``A_c @ model`` maps true-space bins (matrix columns) onto regularized
    bins (matrix rows), so only the row axis changes.
    """
    return str(label).replace("true", "reg")


def save_unfold_ingredient_heatmaps(
    var_config,
    response,
    frac_cov,
    add_smear,
    fig_dir,
    *,
    plot=False,
    approval="",
):
    """Save 2D heatmaps of the response, fractional covariance, and A_c matrices."""
    os.makedirs(fig_dir, exist_ok=True)
    vsn = var_config.var_save_name
    bins = var_config.bins
    labs = getattr(var_config, "var_labels", None) or [vsn, vsn, vsn]
    reco = labs[1] if len(labs) > 1 else labs[0]
    true = labs[2] if len(labs) > 2 else labs[0]
    # Columns contract with the truth vector; rows are the A_c output.
    reg = _true_axis_to_reg(true)

    def _heat(matrix, plot_labels, save_stem):
        plot_heatmap(
            np.asarray(matrix, dtype=float),
            bins,
            plot_labels=plot_labels,
            approval=approval,
            plot=plot,
            save_fig=True,
            save_name=os.path.join(fig_dir, save_stem),
            cmap="viridis",
            cbar_label=False,
            # In-box numbers are off. Set annotate=True to restore them.
            annotate=False,
            corner_text=r"$\mathbf{SBND}$ Simulation",
            halfopen_ticks=True,
            show_title=False,
            fig_width=8.0,
            save_ext=".pdf",
        )

    _heat(response, [true, reco, r"$R$"], f"{vsn}__response")
    _heat(
        frac_cov,
        [reco, reco, r"$C_{\mathrm{frac}}$"],
        f"{vsn}__frac_cov",
    )
    _heat(add_smear, [true, reg, r"$A_c$"], f"{vsn}__add_smear")


def save_unfolded_cov_heatmaps(
    var_config,
    unfold_cov,
    stat_cov=None,
    syst_cov=None,
    fig_dir=".",
    *,
    plot=False,
    approval="",
    pot=UNFOLD_PLOT_POT,
    save_corr=True,
):
    """Save 2D heatmaps of the unfolded covariance matrices.

    Same style as :func:`save_unfold_ingredient_heatmaps` (viridis, no in-box
    numbers, half-open bin ticks, square 8-inch panel, PDF).  The unfolded
    covariances are data products, so the corner text is the POT label rather
    than ``SBND Simulation``.  Both axes are the regularized (unfolded) bins.

    Files: ``{vsn}__unfold_cov`` (total), ``{vsn}__unfold_cov_stat``,
    ``{vsn}__unfold_cov_syst`` and, if ``save_corr``, ``{vsn}__unfold_corr``
    (total correlation, ``coolwarm`` fixed to [-1, 1]).  Single-bin variables
    are skipped.
    """
    unfold_cov = np.asarray(unfold_cov, dtype=float)
    if unfold_cov.ndim != 2 or unfold_cov.shape[0] <= 1:
        return
    os.makedirs(fig_dir, exist_ok=True)
    vsn = var_config.var_save_name
    bins = var_config.bins
    labs = getattr(var_config, "var_labels", None) or [vsn, vsn, vsn]
    true = labs[2] if len(labs) > 2 else labs[0]
    reg = _true_axis_to_reg(true)
    corner = format_pot_corner_text(pot) if pot else ""

    def _heat(matrix, plot_labels, save_stem, cmap="viridis"):
        plot_heatmap(
            np.asarray(matrix, dtype=float),
            bins,
            plot_labels=plot_labels,
            approval=approval,
            plot=plot,
            save_fig=True,
            save_name=os.path.join(fig_dir, save_stem),
            cmap=cmap,
            cbar_label=False,
            annotate=False,
            corner_text=corner,
            halfopen_ticks=True,
            show_title=False,
            fig_width=8.0,
            save_ext=".pdf",
        )

    _heat(unfold_cov, [reg, reg, r"$C_{\mathrm{unfold}}$"], f"{vsn}__unfold_cov")
    if stat_cov is not None:
        _heat(stat_cov, [reg, reg, r"$C_{\mathrm{unfold}}^{\mathrm{stat}}$"], f"{vsn}__unfold_cov_stat")
    if syst_cov is not None:
        _heat(syst_cov, [reg, reg, r"$C_{\mathrm{unfold}}^{\mathrm{syst}}$"], f"{vsn}__unfold_cov_syst")
    if save_corr:
        d = np.sqrt(np.clip(np.diag(unfold_cov), 0.0, None))
        with np.errstate(divide="ignore", invalid="ignore"):
            corr = unfold_cov / np.outer(d, d)
        corr = np.where(np.isfinite(corr), corr, 0.0)
        _heat(corr, [reg, reg, r"$\rho_{\mathrm{unfold}}$"], f"{vsn}__unfold_corr", cmap="coolwarm")


def save_unfold_diagnostics(
    var_config,
    summary,
    response,
    frac_cov,
    add_smear,
    fig_dir,
    *,
    plot=False,
    approval="",
    pot=UNFOLD_PLOT_POT,
):
    """Per-variable unfolding diagnostics: syst-source breakdown + R / C_frac / A_c heatmaps."""
    os.makedirs(fig_dir, exist_ok=True)
    vsn = var_config.var_save_name
    plot_xsec_syst_breakdown(
        var_config,
        summary,
        save_name=os.path.join(fig_dir, f"{vsn}__syst_breakdown"),
        plot=plot,
        approval=approval,
        pot=pot,
    )
    if np.asarray(response).shape[0] > 1:
        save_unfold_ingredient_heatmaps(
            var_config, response, frac_cov, add_smear, fig_dir,
            plot=plot, approval=approval,
        )


# ==== 2D plots ====

_DIVERGING_CMAPS = frozenset(
    {"bwr", "coolwarm", "seismic", "RdBu", "RdBu_r", "coolwarm_r", "bwr_r"}
)


def get_text_color(value, cmap_name="viridis", vmin=0.0, vmax=1.0):
    """Pick black/white annotation color from luminance under *cmap_name*."""
    try:
        cm = plt.get_cmap(cmap_name)
    except Exception:
        cm = mpl.cm.viridis
    if vmax == vmin:
        vmax = vmin + 1e-12
    t = (float(value) - float(vmin)) / (float(vmax) - float(vmin))
    t = min(max(t, 0.0), 1.0)
    rgba = cm(t)
    luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return "black" if luminance > 0.5 else "white"


def _edges_are_integral(edges):
    arr = np.asarray(edges, dtype=float)
    return bool(arr.size) and np.all(np.isfinite(arr) & (np.abs(arr - np.round(arr)) < 1e-6))


def bin_range_labels(edges, halfopen=False):
    """Bin tick labels.

    ``halfopen=False`` is ``min–max`` with two decimals. ``halfopen=True`` is
    ``[min, max)``. Integral edges (``del_alpha``, ``del_phi``) drop the decimals.
    """
    integral = halfopen and _edges_are_integral(edges)
    labels = []
    for i in range(len(edges) - 1):
        if integral:
            lo = str(int(round(float(edges[i]))))
            hi = str(int(round(float(edges[i + 1]))))
        else:
            lo = f"{float(edges[i]):.2f}"
            hi = f"{float(edges[i + 1]):.2f}"
        if halfopen:
            labels.append(f"[{lo}, {hi})")
        else:
            labels.append(f"{lo}–{hi}")
    return labels


def plot_heatmap(matrix, 
                 bins,
                 plot_labels=["", "", ""],
                 approval="",
                 verbose=False,
                 plot=True,
                 cmap="viridis",
                 save_fig=False, 
                 save_name=None,
                 leave_open=False,
                 cbar_label=True,
                 annotate=True,
                 corner_text=None,
                 halfopen_ticks=False,
                 show_title=True,
                 fig_width=None,
                 save_ext=None):
    """2D matrix heatmap (response / cov / A_c / correlation).

    Defaults: no approval stamp; ``viridis`` colorscale.  Diverging cmaps
    (``coolwarm``, ``bwr``, …) are fixed to ``[-1, 1]`` for correlation matrices.

    Unfolding diagnostics pass ``cbar_label=False``, ``annotate=False``, a
    ``corner_text`` (``SBND Simulation``), ``halfopen_ticks=True``,
    ``show_title=False``, and ``fig_width=8`` so the figure is as wide as the
    uncertainty breakdown and the heatmap panel is square. Set ``annotate=True``
    to put bin values back in the cells.
    """

    nbins = len(bins)
    assert nbins-1 == matrix.shape[0] == matrix.shape[1]
    unif_bin = np.linspace(0., float(nbins - 1), nbins)
    extent = [unif_bin[0], unif_bin[-1], unif_bin[0], unif_bin[-1]]

    x_edges, y_edges = np.array(bins), np.array(bins)
    x_tick_positions, y_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2, (unif_bin[:-1] + unif_bin[1:]) / 2
    x_labels = bin_range_labels(x_edges, halfopen=halfopen_ticks)
    y_labels = bin_range_labels(y_edges, halfopen=halfopen_ticks)

    # fig_width matches the uncertainty-breakdown figure. Height is chosen so
    # set_box_aspect(1) leaves a square heatmap beside the colorbar.
    if fig_width is None:
        fig, ax = plt.subplots(figsize=(12, 12))
        tick_fs = None
        label_fs = 20
    else:
        fig, ax = plt.subplots(figsize=(float(fig_width), float(fig_width) * 0.90))
        tick_fs = None
        label_fs = 16
    cmap_name = str(cmap)
    diverging = cmap_name in _DIVERGING_CMAPS
    imshow_kw = dict(extent=extent, origin="lower", cmap=cmap_name)
    if diverging:
        vmin, vmax = -1.0, 1.0
        im = ax.imshow(matrix, vmin=vmin, vmax=vmax, **imshow_kw)
    else:
        flat0 = np.asarray(matrix, dtype=float)
        flat0 = flat0[np.isfinite(flat0)]
        if flat0.size:
            vmin, vmax = float(np.min(flat0)), float(np.max(flat0))
            if vmin == vmax:
                vmax = vmin + 1e-12
        else:
            vmin, vmax = 0.0, 1.0
        im = ax.imshow(matrix, **imshow_kw)

    # Find the power-of-10 exponent from one of the (non-NaN) values
    exponent = 0
    cbar = None
    flat_matrix = matrix[~np.isnan(matrix)]
    if flat_matrix.size > 0 and np.any(flat_matrix != 0):
        example_value = flat_matrix[0]
        exponent = np.floor(np.log10(abs(example_value)))
        # if exponent is infinite (matrix contains only zeros), set to 0
        if np.isinf(exponent):
            exponent = 0
        exponent = int(exponent)

        # # Round to nearest multiple of 3 for best tick label presentation
        # exponent = 3 * int(np.floor(exponent / 3))

        formatter = mpl.ticker.FuncFormatter(lambda x, _: f"{x/10**exponent:.2f}")
        if fig_width is None:
            cbar = plt.colorbar(shrink=0.7)
        else:
            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        # Colorbar title. Unfolding heatmaps pass cbar_label=False.
        if cbar_label:
            # Values in ~0.1–1 (exponent -1) keep native labels, matching cell annotations.
            if exponent != 0 and exponent != -1:
                cbar.set_label(plot_labels[2] + f" [10$^{{{exponent}}}$]", fontsize=16)
                cbar.ax.yaxis.set_major_formatter(formatter)
            else:
                cbar.set_label(plot_labels[2], fontsize=16)

    # else:
    #     plt.colorbar(shrink=0.7, label=plot_labels[2])

    # In-box bin values. Unfolding heatmaps pass annotate=False so this loop
    # does not run. Set annotate=True (the default) to restore the numbers.
    if annotate:
        for i in range(nbins-1):      # rows (y)
            for j in range(nbins-1):  # columns (x)
                value = matrix[i, j]
                if not np.isnan(value):  # skip NaNs
                    if exponent != -1:
                        significand = value / 10**exponent
                    else:
                        significand = value
                    plt.text(
                        j + 0.5, i + 0.5,
                        f"{significand:.2f}",
                        ha="center", va="center",
                        color=get_text_color(value, cmap_name, vmin, vmax),
                        fontsize=10
                    )

    if fig_width is None:
        plt.xticks(x_tick_positions, x_labels, rotation=45, ha="right")
        plt.yticks(y_tick_positions, y_labels)
    else:
        if cbar is not None:
            cbar_labels = cbar.ax.get_yticklabels()
            if cbar_labels:
                tick_fs = cbar_labels[0].get_fontproperties().get_size_in_points()
        ax.set_xticks(x_tick_positions)
        ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=tick_fs)
        ax.set_yticks(y_tick_positions)
        ax.set_yticklabels(y_labels, fontsize=tick_fs)
        ax.set_box_aspect(1)
    ax.set_xlabel(plot_labels[0], fontsize=label_fs)
    ax.set_ylabel(plot_labels[1], fontsize=label_fs)
    if show_title:
        if len(plot_labels) > 3:
            ax.set_title(plot_labels[3], fontsize=20)
        else:
            ax.set_title(plot_labels[2], fontsize=20)

    if fig_width is not None:
        fig.tight_layout()

    if verbose:
        n_diag = np.sum(np.diag(matrix))
        diagonal_ratio = n_diag / np.sum(matrix)
        print(f"Diagonal ratio: {diagonal_ratio:.2f}")
        print(f"True ratio: {np.diag(matrix) / np.sum(matrix, axis=0)}")

        # print

    # ===== plot additions =====
    if corner_text:
        # subplots_adjust would undo the square panel set up for fig_width.
        if fig_width is None:
            fig.subplots_adjust(top=0.90)
        add_approval_text(approval, 0.03, 1.08, "left")
        add_pot_text(corner_text, 0.99, 1.02, "right", fontsize=16)
    else:
        add_approval_text(approval, 0.95, 1.05, "right")

    if save_fig:
        plt.savefig(save_name + (fig_ext if save_ext is None else save_ext), bbox_inches="tight", dpi=dpi)

    if plot:
        plt.show()
    elif not leave_open:
        plt.close()




####
# Exposure Accounting

# ====== flux, detector geometry, and cross-section normalization ======
def get_integrated_flux(fluxfile, plot=False):
    # Gen1 / sbnd_original_flux TH1D convention: bin *contents* are already
    # per-bin rates in /m^2/10^6 POT (50 MeV bins), NOT dΦ/dE densities.
    # Integrated Φ = sum(content). Do **not** use ROOT Integral("width") here —
    # that multiplies by ΔE=0.05 and undercounts by ×20. (NUISANCE fScaleFactor
    # uses Integral("width") on *both* event and flux hists, so the ΔE cancels
    # in the flux-averaged σ; analysis XSEC_UNIT needs the true sum.)
    flux = uproot.open(fluxfile)
    numu_flux = flux["flux_sbnd_numu"].to_numpy()
    bin_edges = numu_flux[1]
    flux_vals = numu_flux[0]

    if plot:
        fig, ax = plt.subplots()
        plt.hist(bin_edges[:-1], bins=bin_edges, weights=flux_vals, histtype="step", linewidth=2, color="C0")
        plt.xlim(0, 3)
        plt.xlabel("Neutrino Energy [GeV]")
        plt.ylabel("Flux [/m$^{2}$/10$^{6}$ POT]")
        plt.title("SBND $\\nu_\\mu$ Flux")
        plt.savefig("sbnd-flux.pdf", bbox_inches='tight')

    integrated_flux = flux_vals.sum() / (1e4  * 1e6) # to cm2 # to POT
    print("Integrated flux: %.3e" % integrated_flux)
    return integrated_flux


def get_active_volume(detector="SBND"):
    if detector == "SBND":
        V_SBND = 380 * 380 * 440 # cm3, the active volume of the detector 

    elif detector == "SBND_nohighyz":
        V_SBND = 380 * 380 * 440 - 380* (190 - 100) *(450-250) 

    elif detector == "SBND_face":
        V_SBND = 380 * 380 * 50

    elif detector == "SBND_face_yzcut":
        V_SBND = 380 * 380 * 50 - 380* (190 - 100) * 50 

    elif detector == "SBND_end":
        V_SBND = 380 * 380 * 50

    return V_SBND


def print_sbnd_octant_vertex_ranges(x0, y0, z0):
    """Print octant labels vs reco vertex (x, y, z) in cm.

    Matches the convention in ``selected_events`` octant labeling: split planes at
    ``x0``, ``y0``, ``z0``; E/W from ``x``, N/S from ``z``, Top/Bottom from ``y``.
    SBND: **East** is ``x < x0``; **West** is ``x >= x0``.
    **South** is ``z < z0``; **North** is ``z > z0`` (``z >= z0`` at the split plane).
    **Top** is ``y >= y0``; **Bottom** is ``y < y0``.
    """
    xf, yf, zf = float(x0), float(y0), float(z0)
    xs, ys, zs = "{:.6g}".format(xf), "{:.6g}".format(yf), "{:.6g}".format(zf)
    print("\n=== SBND octants vs reco vertex [cm]; planes x={}, y={}, z={} ===".format(xs, ys, zs))
    print("  E/W (TPC sides):  E if x < {} (negative x),    W if x >= {}".format(xs, xs))
    print("  N/S:              S if z < {} (lower z),    N if z >= {}".format(zs, zs))
    print("  Top / Bottom:     Bottom if y < {},    Top if y >= {}".format(ys, ys))
    # x: E → x < x0 ; W → x >= x0.  z: S → z < z0 ; N → z >= z0.
    rows = [
        ("W-S-Bottom", "[{}, +inf)".format(xs), "(-inf, {})".format(zs), "(-inf, {})".format(ys)),
        ("W-S-Top", "[{}, +inf)".format(xs), "(-inf, {})".format(zs), "[{}, +inf)".format(ys)),
        ("W-N-Bottom", "[{}, +inf)".format(xs), "[{}, +inf)".format(zs), "(-inf, {})".format(ys)),
        ("W-N-Top", "[{}, +inf)".format(xs), "[{}, +inf)".format(zs), "[{}, +inf)".format(ys)),
        ("E-S-Bottom", "(-inf, {})".format(xs), "(-inf, {})".format(zs), "(-inf, {})".format(ys)),
        ("E-S-Top", "(-inf, {})".format(xs), "(-inf, {})".format(zs), "[{}, +inf)".format(ys)),
        ("E-N-Bottom", "(-inf, {})".format(xs), "[{}, +inf)".format(zs), "(-inf, {})".format(ys)),
        ("E-N-Top", "(-inf, {})".format(xs), "[{}, +inf)".format(zs), "[{}, +inf)".format(ys)),
    ]
    hdr = "{:14}  {:^26}  {:^26}  {:^26}".format("octant", "x range", "z range", "y range")
    print("\n" + hdr)
    print(" " + "-" * (len(hdr) + 2))
    for name, xr, zr, yr in rows:
        print("{:14}  {:^26}  {:^26}  {:^26}".format(name, xr, zr, yr))
    print(
        "\nBoundary vertices: x=x0 uses >= toward West; z=z0 uses >= toward North; y=y0 uses >= toward Top.\n"
    )


def get_xsec_unit(tot_pot, 
                  fluxfile="/exp/sbnd/data/users/munjung/flux/sbnd_original_flux.root", 
                  detector="SBND",
                  volume=None):

    tot_flux = get_integrated_flux(fluxfile, plot=False)
    tot_flux *= tot_pot
    print("integrated flux: ", tot_flux)

    V_SBND = get_active_volume(detector)
    if volume is not None:
        print("using custom volume: ", volume)
        V_SBND = volume

    # Argon *nuclei* (not nucleons). Matches flat weight ``40 × fScaleFactor``:
    # PrepareGENIE divides event rate by totalnucl=A=40 (per-nucleon bookkeeping);
    # the analysis 40× restores per-nucleus σ to pair with this NTARGETS.
    NTARGETS = RHO * V_SBND * (N_A / M_AR)
    print("# of targets: ", NTARGETS)

    xsec_unit = 1 / (tot_flux * NTARGETS)
    # # TODO: fix scalar overflow error in python v3.10+
    # if xsec_unit == 0:
    #     print("XSEC_UNIT is 0, setting to 1e-38")
    #     xsec_unit = 1e-38
    print("xsec unit: ", xsec_unit)
    return xsec_unit

    
