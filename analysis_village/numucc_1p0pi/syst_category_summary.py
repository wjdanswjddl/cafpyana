"""Export / load per-category fractional systematics for plotting and unfolding.

Written by ``systematics-summary.ipynb`` into ``<SYST_DISK_ROOT>/CategorySummary/``.
Each variable entry stores fractional covariance matrices and per-bin uncertainties in percent,
using the same conventions as the summary breakdown plots (cosmics = contamination-scaled
``SelectedRate``; GENIE rate and xsec totals available separately).

``categories`` keep unrebased source fracs (Flux/G4/MCstat/GENIE rate = frac vs signal CV).
``total_rate`` optionally rebases those signal-CV sources onto total-selected CV before
summing (see :func:`rebase_fraccov_signal_to_total`); ``total_xsec`` is never rebased.

Example
-------
>>> from analysis_village.numucc_1p0pi.syst_category_summary import (
...     load_category_syst_summary,
...     category_cov_frac,
...     category_frac_unc_pct,
...     total_frac_unc_pct,
... )
>>> summary = load_category_syst_summary()
>>> vsn = "muon-p"
>>> c_cosmics = category_cov_frac(summary, vsn, "cosmics")
>>> u_total = total_frac_unc_pct(summary, vsn, kind="xsec")
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, Mapping, MutableMapping, Optional, Sequence

import numpy as np

from pyanalib.covariance import corr_from_fraccov, cov_from_fraccov, fraccov_from_cov

from analysis_village.numucc_1p0pi.syst_disk_layout import SYST_DISK_ENV, normalized_root

# Local path helpers (mirror ``syst_disk_layout``; kept here so a stale cached
# ``syst_disk_layout`` in a long-lived notebook kernel cannot break import).
_SUB_CATEGORY_SUMMARY = "CategorySummary"
_FILE_CATEGORY_SUMMARY = "category_syst_summary.npz"
_FILE_CATEGORY_SUMMARY_MANIFEST = "category_syst_summary_manifest.json"


def category_summary_npz_path(root: str) -> str:
    return os.path.join(
        normalized_root(root), _SUB_CATEGORY_SUMMARY, _FILE_CATEGORY_SUMMARY
    )


def category_summary_manifest_path(npz_path: str) -> str:
    d = os.path.dirname(os.path.abspath(npz_path))
    return os.path.join(d, _FILE_CATEGORY_SUMMARY_MANIFEST)


SCHEMA = "numucc_category_syst_summary_v1"

# Stable keys for ``pack["categories"]`` and ``load_category_syst_summary`` lookups.
CAT_FLUX = "flux"
CAT_G4 = "g4"
CAT_MCSTAT = "mcstat"
CAT_DETECTOR = "detector"
CAT_COSMICS = "cosmics"  # SelectedRate (contamination-scaled), not raw template
CAT_GENIE_RATE = "genie_rate"
CAT_GENIE_XSEC = "genie_xsec"
CAT_POT = "pot"
CAT_NTARGETS = "ntargets"

CATEGORY_KEYS = (
    CAT_FLUX,
    CAT_G4,
    CAT_MCSTAT,
    CAT_DETECTOR,
    CAT_COSMICS,
    CAT_GENIE_RATE,
    CAT_GENIE_XSEC,
    CAT_POT,
    CAT_NTARGETS,
)

TOTAL_RATE = "total_rate"
TOTAL_XSEC = "total_xsec"
TOTAL_KEYS = (TOTAL_RATE, TOTAL_XSEC)

# Flat normalization scales (percent): one global factor → 100% correlated across bins.
# (Exposure and Targets are independent of each other, added as separate categories.)
POT_FRAC_UNC_PCT = 2.0   # Exposure / POT
NTARGETS_FRAC_UNC_PCT = 1.0  # Number of targets

_RATE_TOTAL_CATEGORIES = (
    CAT_FLUX,
    CAT_G4,
    CAT_MCSTAT,
    CAT_DETECTOR,
    CAT_COSMICS,
    CAT_GENIE_RATE,
    CAT_POT,
    CAT_NTARGETS,
)
_XSEC_TOTAL_CATEGORIES = (
    CAT_FLUX,
    CAT_G4,
    CAT_MCSTAT,
    CAT_DETECTOR,
    CAT_COSMICS,
    CAT_GENIE_XSEC,
    CAT_POT,
    CAT_NTARGETS,
)

# Multisim / GENIE **rate** fracs are built with ``bkgd_subtract=True`` (fractional
# vs signal CV). Overlay bands apply ``frac × total_mc``, so these sources must be
# rebased onto total-selected CV before entering ``total_rate``.
# Do **not** rebase detector, cosmics ``SelectedRate``, pot, or ntargets.
RATE_REBASE_TO_TOTAL_CATEGORIES = (
    CAT_FLUX,
    CAT_G4,
    CAT_MCSTAT,
    CAT_GENIE_RATE,
)


def rebase_fraccov_signal_to_total(
    frac_sig: np.ndarray,
    n_signal: np.ndarray,
    n_total: np.ndarray,
) -> np.ndarray:
    """Rebase a fractional covariance from signal CV onto total-selected CV.

    Flux / G4 / MCstat / GENIE-rate matrices are produced with
    ``bkgd_subtract=True``, so ``cov_frac[i,j]`` is
    ``(δ_i / n_signal_i)(δ_j / n_signal_j)``. Overlay hatch / χ² paths convert
    frac → absolute with ``cov_from_fraccov(frac, total_mc)``, which under-scales
    the absolute uncertainty unless the frac is first rewritten vs ``n_total``::

        abs = cov_from_fraccov(frac_sig, n_signal)
        frac_tot = fraccov_from_cov(abs, n_total)

    Parameters
    ----------
    frac_sig
        Fractional covariance relative to the signal CV.
    n_signal, n_total
        Per-bin signal and total-selected (stack) counts; same length as the
        matrix dimension.

    Returns
    -------
    np.ndarray
        Fractional covariance relative to ``n_total``.

    Notes
    -----
    Do **not** apply this to detector, cosmics ``SelectedRate``, pot, or
    ntargets — those are already fractional vs the total-selected (or are
    pure multiplicative scales on the full stack).
    """
    frac_sig = np.asarray(frac_sig, dtype=np.float64)
    n_signal = np.asarray(n_signal, dtype=np.float64).reshape(-1)
    n_total = np.asarray(n_total, dtype=np.float64).reshape(-1)
    if frac_sig.ndim != 2 or frac_sig.shape[0] != frac_sig.shape[1]:
        raise ValueError(
            "frac_sig must be square; got shape %s" % (frac_sig.shape,)
        )
    n = frac_sig.shape[0]
    if n_signal.shape[0] != n or n_total.shape[0] != n:
        raise ValueError(
            "n_signal/n_total length (%d, %d) != cov dim %d"
            % (n_signal.shape[0], n_total.shape[0], n)
        )
    abs_cov = cov_from_fraccov(frac_sig, n_signal)
    return np.asarray(fraccov_from_cov(abs_cov, n_total), dtype=np.float64)


def assemble_total_rate_cov_frac(
    category_covs: Mapping[str, np.ndarray],
    *,
    n_signal: Optional[np.ndarray] = None,
    n_total: Optional[np.ndarray] = None,
    rebase: bool = True,
) -> Optional[np.ndarray]:
    """Sum rate-total category fracs, optionally rebasing signal-CV sources.

    ``category_covs`` keys match :data:`CATEGORY_KEYS` (``flux``, ``g4``, …).
    When ``rebase`` is True and both ``n_signal`` / ``n_total`` are supplied,
    entries in :data:`RATE_REBASE_TO_TOTAL_CATEGORIES` are passed through
    :func:`rebase_fraccov_signal_to_total` before summing; other categories
    (detector, cosmics, pot, ntargets) are left unchanged.

    Returns ``None`` if no category matrices are present.
    """
    parts: list[np.ndarray] = []
    do_rebase = bool(rebase) and n_signal is not None and n_total is not None
    for key in _RATE_TOTAL_CATEGORIES:
        if key not in category_covs or category_covs[key] is None:
            continue
        cf = np.asarray(category_covs[key], dtype=np.float64)
        if do_rebase and key in RATE_REBASE_TO_TOTAL_CATEGORIES:
            cf = rebase_fraccov_signal_to_total(cf, n_signal, n_total)
        parts.append(cf)
    return sum_cov_frac_matrices(parts)


def frac_unc_diag(cov_frac: np.ndarray) -> np.ndarray:
    c = np.asarray(cov_frac, dtype=np.float64)
    return np.nan_to_num(np.sqrt(np.diag(c)), nan=0.0, posinf=0.0, neginf=0.0)


def frac_unc_pct(cov_frac: np.ndarray) -> np.ndarray:
    """Per-bin fractional uncertainty in percent (100 × sqrt(diag(cov_frac)))."""
    return 100.0 * frac_unc_diag(cov_frac)


def frac_weights_for_plot(cov_frac: np.ndarray, var_config: Any) -> np.ndarray:
    """Per-bin uncertainty [%] with bin-count alignment and integrated flat maximum."""
    w = frac_unc_pct(cov_frac)
    n = len(var_config.bin_centers)
    if len(w) != n:
        if len(w) == 1 and n > 1:
            w = np.full(n, float(w[0]))
        elif len(w) > n:
            w = np.asarray(w[:n], dtype=np.float64)
        elif len(w) < n:
            out = np.zeros(n, dtype=np.float64)
            out[: len(w)] = w
            w = out
    if getattr(var_config, "var_save_name", None) == "integrated" and len(w):
        w = np.full(n, float(np.nanmax(w)))
    return w


def integrated_rate_frac_variance(cov_frac: np.ndarray) -> float:
    c = np.asarray(cov_frac, dtype=np.float64)
    ones = np.ones(c.shape[0], dtype=np.float64)
    return float(ones @ c @ ones)


def sum_cov_frac_matrices(matrices: Iterable[np.ndarray]) -> Optional[np.ndarray]:
    total = None
    for m in matrices:
        arr = np.asarray(m, dtype=np.float64)
        total = arr if total is None else total + arr
    return total


def cosmics_selected_rate_cov_frac(cosmics_npz: Any, var_name: str) -> Optional[np.ndarray]:
    if cosmics_npz is None or var_name not in cosmics_npz:
        return None
    cell = cosmics_npz[var_name].item()
    if isinstance(cell, dict) and "SelectedRate" in cell:
        rate = cell["SelectedRate"].get("rate")
        if isinstance(rate, dict) and "cov_frac" in rate:
            return np.asarray(rate["cov_frac"], dtype=np.float64)
    return None


def genie_var_dict(genie_blob: Optional[Mapping], var_name: str) -> Optional[Mapping]:
    if genie_blob is None or var_name not in genie_blob:
        return None
    cell = genie_blob[var_name]
    if hasattr(cell, "item"):
        cell = cell.item()
    return cell if isinstance(cell, dict) else None


def _genie_multisim_bundle_pack(gd: Mapping) -> Optional[Dict[str, Any]]:
    """Bundled GENIE block from ``save_neutrino_multisim_npzs`` (``{GENIE: {cov_frac, ...}}``).

    ``systematics-genie.ipynb`` writes one combined matrix per variable; ``cov_type`` in the
    pack (or default ``rate``) selects ``rate_total`` vs ``xsec_total``.
    """
    pay = gd.get("GENIE")
    if not isinstance(pay, dict) or pay.get("cov_frac") is None:
        return None
    cf = np.asarray(pay["cov_frac"], dtype=np.float64)
    cov_type = str(pay.get("cov_type", "rate")).lower()
    empty = {"rate_parts": {}, "xsec_parts": {}}
    if cov_type == "xsec":
        return {**empty, "rate_total": None, "xsec_total": cf}
    return {**empty, "rate_total": cf, "xsec_total": None}


def genie_knob_covs(gd: Optional[Mapping]) -> Optional[Dict[str, Any]]:
    if gd is None:
        return None
    if hasattr(gd, "item"):
        gd = gd.item()
    if not isinstance(gd, dict):
        return None
    bundled = _genie_multisim_bundle_pack(gd)
    if bundled is not None:
        return bundled
    rate_parts: Dict[str, np.ndarray] = {}
    xsec_parts: Dict[str, np.ndarray] = {}
    rate_total = xsec_total = None
    for k, v in gd.items():
        if not isinstance(v, np.ndarray) or v.ndim != 2:
            continue
        if k == "genie":
            xsec_total = v
        elif k == "genie_rate":
            rate_total = v
        elif k.endswith("_rate"):
            rate_parts[k[: -len("_rate")]] = v
        else:
            xsec_parts[k] = v
    return {
        "rate_parts": rate_parts,
        "xsec_parts": xsec_parts,
        "rate_total": rate_total,
        "xsec_total": xsec_total,
    }


def _detector_subsystem_total_cov_frac(
    det_npz: Any, var_name: str
) -> Optional[np.ndarray]:
    """Combined detector fractional covariance from one ``detector_syst_dict.npz``."""
    if det_npz is None:
        return None
    z = dict(det_npz)
    total = None
    combined = z.get("detector")
    if combined is not None:
        comb_item = combined.item() if hasattr(combined, "item") else combined
        if isinstance(comb_item, dict) and var_name in comb_item:
            total = comb_item[var_name]["cov_frac"]
    dbw = z.get("detector_by_wiremod")
    if dbw is not None:
        cell = dbw.item() if hasattr(dbw, "item") else dbw
        if isinstance(cell, dict) and var_name in cell:
            per_var = cell[var_name]
            parts = [
                pack["cov_frac"]
                for pack in per_var.values()
                if isinstance(pack, dict) and "cov_frac" in pack
            ]
            if parts:
                subtotal = sum_cov_frac_matrices(parts)
                if total is None:
                    total = subtotal
            if total is not None:
                return np.asarray(total, dtype=np.float64)
    out = {}
    for k in sorted(z.keys()):
        if not k.startswith("detector-"):
            continue
        tag = k[len("detector-") :]
        item = z[k].item() if hasattr(z[k], "item") else z[k]
        if var_name not in item:
            continue
        out[tag] = item[var_name]["cov_frac"]
    if total is None and out:
        total = sum_cov_frac_matrices(out.values())
    if total is None:
        return None
    return np.asarray(total, dtype=np.float64)


def detector_total_cov_frac(
    var_name: str,
    *,
    detector_npz: Any = None,
    wiremod_npz: Any = None,
    dent_npz: Any = None,
    sce_npz: Any = None,
) -> Optional[np.ndarray]:
    """Detector category total from the combined Product **B** NPZ.

    Canonical input is ``detector_npz`` from ``systematics-detector.ipynb``
    (``Detector/detector_syst_dict.npz`` = WireMod YZ + XTXW + DENT).

    ``wiremod_npz`` / ``dent_npz`` / ``sce_npz`` remain as optional legacy
    fallbacks when the combined file is absent.
    """
    combined = _detector_subsystem_total_cov_frac(detector_npz, var_name)
    if combined is not None:
        return combined

    totals = []
    for det_npz in (wiremod_npz, dent_npz):
        cov = _detector_subsystem_total_cov_frac(det_npz, var_name)
        if cov is not None:
            totals.append(cov)
    if totals:
        out = sum_cov_frac_matrices(totals)
        return None if out is None else np.asarray(out, dtype=np.float64)

    if sce_npz is not None:
        totals = []
        for det_npz in (wiremod_npz, sce_npz):
            cov = _detector_subsystem_total_cov_frac(det_npz, var_name)
            if cov is not None:
                totals.append(cov)
        if totals:
            out = sum_cov_frac_matrices(totals)
            return None if out is None else np.asarray(out, dtype=np.float64)
    return None


def mcstat_cov_frac(mcstat_npz: Any, var_name: str) -> Optional[np.ndarray]:
    if mcstat_npz is None or var_name not in mcstat_npz:
        return None
    cell = mcstat_npz[var_name].item()
    if isinstance(cell, dict) and "MCstat" in cell:
        pay = cell["MCstat"]
        if isinstance(pay, dict) and pay.get("cov_frac") is not None:
            return np.asarray(pay["cov_frac"], dtype=np.float64)
    return None


def genie_category_cov_frac(genie_pack: Optional[Mapping], kind: str) -> Optional[np.ndarray]:
    if genie_pack is None:
        return None
    if kind == "rate":
        if genie_pack["rate_total"] is not None:
            return np.asarray(genie_pack["rate_total"], dtype=np.float64)
        parts = genie_pack.get("rate_parts") or {}
        if not parts:
            return None
        return sum_cov_frac_matrices(parts.values())
    if genie_pack["xsec_total"] is not None:
        return np.asarray(genie_pack["xsec_total"], dtype=np.float64)
    parts = genie_pack.get("xsec_parts") or {}
    if not parts:
        return None
    return sum_cov_frac_matrices(parts.values())


def _flat_cov_frac(nbins: int, frac_unc_pct_val: float) -> np.ndarray:
    """Fractional cov for a single multiplicative scale (POT, N_targets, …).

    Same fractional shift in every bin → ``cov_frac[i,j] = (pct/100)²`` for all i, j.
    Per-bin curves still show ``pct``%; off-diagonals are non-zero in heatmaps/totals.
    """
    u = float(frac_unc_pct_val) / 100.0
    v = u * u
    return np.full((int(nbins), int(nbins)), v, dtype=np.float64)


def _category_block(
    cov_frac: np.ndarray,
    var_config: Any,
    *,
    nominal_mc: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    cov_frac = np.asarray(cov_frac, dtype=np.float64)
    w = frac_weights_for_plot(cov_frac, var_config)
    block: Dict[str, np.ndarray] = {
        "cov_frac": cov_frac,
        "corr": np.asarray(corr_from_fraccov(cov_frac), dtype=np.float64),
        "frac_unc_pct": w,
        "integrated_frac_variance": np.array(integrated_rate_frac_variance(cov_frac)),
    }
    if nominal_mc is not None:
        mc = np.asarray(nominal_mc, dtype=np.float64)
        block["cov"] = np.asarray(cov_from_fraccov(cov_frac, mc), dtype=np.float64)
    return block


def build_category_cov_frac(
    vsn: str,
    *,
    flux_npz: Any,
    g4_npz: Any,
    cosmics_npz: Any,
    detector_npz: Any = None,
    wiremod_npz: Any = None,
    dent_npz: Any = None,
    sce_npz: Any = None,
    mcstat_npz: Any = None,
    genie_blob: Optional[Mapping] = None,
) -> Dict[str, np.ndarray]:
    """Return ``{category_key: cov_frac}`` for one variable slug."""
    covs: Dict[str, np.ndarray] = {
        CAT_FLUX: np.asarray(flux_npz[vsn].item()["flux"]["cov_frac"], dtype=np.float64),
        CAT_G4: np.asarray(g4_npz[vsn].item()["G4"]["cov_frac"], dtype=np.float64),
    }
    mc = mcstat_cov_frac(mcstat_npz, vsn)
    if mc is not None:
        covs[CAT_MCSTAT] = mc
    det = detector_total_cov_frac(
        vsn,
        detector_npz=detector_npz,
        wiremod_npz=wiremod_npz,
        dent_npz=dent_npz,
        sce_npz=sce_npz,
    )
    if det is not None:
        covs[CAT_DETECTOR] = det
    cc = cosmics_selected_rate_cov_frac(cosmics_npz, vsn)
    if cc is not None:
        covs[CAT_COSMICS] = cc
    gp = genie_knob_covs(genie_var_dict(genie_blob, vsn))
    if gp is not None:
        gr = genie_category_cov_frac(gp, "rate")
        if gr is not None:
            covs[CAT_GENIE_RATE] = gr
        gx = genie_category_cov_frac(gp, "xsec")
        if gx is not None:
            covs[CAT_GENIE_XSEC] = gx
    return covs


def build_variable_pack(
    var_config: Any,
    *,
    flux_npz: Any,
    g4_npz: Any,
    cosmics_npz: Any,
    detector_npz: Any = None,
    wiremod_npz: Any = None,
    dent_npz: Any = None,
    sce_npz: Any = None,
    mcstat_npz: Any = None,
    genie_blob: Optional[Mapping] = None,
    include_flat: bool = True,
    nominal_mc: Optional[np.ndarray] = None,
    nominal_mc_signal: Optional[np.ndarray] = None,
    nominal_mc_total: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """Full export payload for one ``VariableConfig``.

    ``categories`` store **unrebased** source fractional covariances (Flux / G4 /
    MCstat / GENIE rate remain frac-vs-signal CV as written by the multisim
    producers). ``total_rate`` rebases those sources onto total-selected CV when
    ``nominal_mc_signal`` and ``nominal_mc_total`` are both provided; detector /
    cosmics / pot / ntargets are summed as-is. ``total_xsec`` always uses the
    unrebased category fracs (identical to the pre-rebase convention).
    """
    vsn = var_config.var_save_name
    nbins = len(var_config.bin_centers)
    covs = build_category_cov_frac(
        vsn,
        flux_npz=flux_npz,
        g4_npz=g4_npz,
        cosmics_npz=cosmics_npz,
        detector_npz=detector_npz,
        wiremod_npz=wiremod_npz,
        dent_npz=dent_npz,
        sce_npz=sce_npz,
        mcstat_npz=mcstat_npz,
        genie_blob=genie_blob,
    )
    if include_flat:
        covs[CAT_POT] = _flat_cov_frac(nbins, POT_FRAC_UNC_PCT)
        covs[CAT_NTARGETS] = _flat_cov_frac(nbins, NTARGETS_FRAC_UNC_PCT)

    # Per-category blocks keep unrebased source fracs (documented above).
    categories: Dict[str, Dict[str, np.ndarray]] = {
        k: _category_block(c, var_config, nominal_mc=nominal_mc)
        for k, c in covs.items()
    }

    def _total_block(
        keys: Sequence[str],
        *,
        cov_frac: np.ndarray,
        abs_mc: Optional[np.ndarray],
    ) -> Dict[str, np.ndarray]:
        blk = _category_block(cov_frac, var_config, nominal_mc=abs_mc)
        blk["category_keys"] = np.array(list(keys), dtype=object)
        return blk

    pack: Dict[str, Any] = {
        "schema": SCHEMA,
        "var_save_name": vsn,
        "bins": np.asarray(var_config.bins, dtype=np.float64),
        "bin_centers": np.asarray(var_config.bin_centers, dtype=np.float64),
        "categories": categories,
        "cosmics_kind": "selected_rate",
    }

    rate_total = assemble_total_rate_cov_frac(
        covs,
        n_signal=nominal_mc_signal,
        n_total=nominal_mc_total,
        rebase=(nominal_mc_signal is not None and nominal_mc_total is not None),
    )
    if rate_total is not None:
        # Absolute cov for total_rate prefers the stack CV used after rebase.
        abs_mc_rate = (
            nominal_mc_total
            if nominal_mc_total is not None
            else nominal_mc
        )
        pack[TOTAL_RATE] = _total_block(
            _RATE_TOTAL_CATEGORIES, cov_frac=rate_total, abs_mc=abs_mc_rate
        )

    xsec_parts = [covs[k] for k in _XSEC_TOTAL_CATEGORIES if k in covs]
    if xsec_parts:
        xsec_total = sum_cov_frac_matrices(xsec_parts)
        assert xsec_total is not None
        pack[TOTAL_XSEC] = _total_block(
            _XSEC_TOTAL_CATEGORIES, cov_frac=xsec_total, abs_mc=nominal_mc
        )
    return pack


def export_category_syst_summary(
    out_npz: str,
    var_configs: Sequence[Any],
    *,
    flux_npz: Any,
    g4_npz: Any,
    cosmics_npz: Any,
    detector_npz: Any = None,
    wiremod_npz: Any = None,
    dent_npz: Any = None,
    sce_npz: Any = None,
    mcstat_npz: Any = None,
    genie_blob: Optional[Mapping] = None,
    syst_disk_root: Optional[str] = None,
    include_flat: bool = True,
    nominal_mc_by_var: Optional[Mapping[str, np.ndarray]] = None,
    nominal_mc_signal_by_var: Optional[Mapping[str, np.ndarray]] = None,
    nominal_mc_total_by_var: Optional[Mapping[str, np.ndarray]] = None,
) -> Dict[str, Any]:
    """Write ``category_syst_summary.npz`` and companion manifest JSON."""
    out_npz = os.path.abspath(out_npz)
    os.makedirs(os.path.dirname(out_npz), exist_ok=True)
    manifest_path = category_summary_manifest_path(out_npz)
    nominal_mc_by_var = nominal_mc_by_var or {}
    nominal_mc_signal_by_var = nominal_mc_signal_by_var or {}
    nominal_mc_total_by_var = nominal_mc_total_by_var or {}

    by_var: Dict[str, Any] = {}
    skipped: list[str] = []
    for var_config in var_configs:
        vsn = var_config.var_save_name
        try:
            by_var[vsn] = build_variable_pack(
                var_config,
                flux_npz=flux_npz,
                g4_npz=g4_npz,
                cosmics_npz=cosmics_npz,
                detector_npz=detector_npz,
                wiremod_npz=wiremod_npz,
                dent_npz=dent_npz,
                sce_npz=sce_npz,
                mcstat_npz=mcstat_npz,
                genie_blob=genie_blob,
                include_flat=include_flat,
                nominal_mc=nominal_mc_by_var.get(vsn),
                nominal_mc_signal=nominal_mc_signal_by_var.get(vsn),
                nominal_mc_total=nominal_mc_total_by_var.get(vsn),
            )
        except Exception as ex:
            skipped.append(f"{vsn}: {ex}")

    np.savez_compressed(out_npz, **by_var)
    manifest = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "syst_disk_root": normalized_root(syst_disk_root) if syst_disk_root else None,
        "npz_path": out_npz,
        "variables": sorted(by_var.keys()),
        "skipped": skipped,
        "category_keys": list(CATEGORY_KEYS),
        "total_keys": list(TOTAL_KEYS),
        "cosmics_kind": "selected_rate",
        "flat_frac_unc_pct": {"pot": POT_FRAC_UNC_PCT, "ntargets": NTARGETS_FRAC_UNC_PCT},
        "matrix_keys_per_category": ["cov_frac", "corr", "cov", "frac_unc_pct"],
        "rate_rebase_to_total_categories": list(RATE_REBASE_TO_TOTAL_CATEGORIES),
        "note_total_rate": (
            "total_rate rebases flux/g4/mcstat/genie_rate from signal CV onto "
            "total-selected CV when nominal_mc_signal and nominal_mc_total are "
            "supplied at export; categories[] keep unrebased source fracs. "
            "total_xsec is never rebased."
        ),
        "usage": {
            "load": "load_category_syst_summary(npz_path) or load_category_syst_summary(syst_disk_root=...)",
            "cov_frac": "category_cov_frac(summary, var_save_name, category_key)",
            "corr": "pack['categories'][key]['corr']",
            "cov": "pack['categories'][key]['cov'] (when nominal MC supplied at export)",
            "frac_unc_pct": "category_frac_unc_pct(summary, var_save_name, category_key)",
            "total": "total_frac_unc_pct(summary, var_save_name, kind='xsec'|'rate')",
        },
    }
    with open(manifest_path, "w", encoding="utf-8") as mf:
        json.dump(manifest, mf, indent=2)
        mf.write("\n")

    return manifest


def rebuild_total_rate_in_summary_pack(
    pack: MutableMapping[str, Any],
    var_config: Any,
    *,
    n_signal: np.ndarray,
    n_total: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Rewrite ``pack[TOTAL_RATE]`` from unrebased ``categories`` + overlay counts.

    Leaves ``categories`` and ``total_xsec`` untouched. Returns the new
    ``total_rate`` block.
    """
    cats = pack.get("categories") or {}
    cat_covs = {
        k: np.asarray(cats[k]["cov_frac"], dtype=np.float64)
        for k in _RATE_TOTAL_CATEGORIES
        if k in cats and isinstance(cats[k], Mapping) and "cov_frac" in cats[k]
    }
    rate_total = assemble_total_rate_cov_frac(
        cat_covs, n_signal=n_signal, n_total=n_total, rebase=True
    )
    if rate_total is None:
        raise ValueError(
            "No rate category matrices to assemble for %r"
            % (getattr(var_config, "var_save_name", "?"),)
        )
    blk = _category_block(rate_total, var_config, nominal_mc=n_total)
    blk["category_keys"] = np.array(list(_RATE_TOTAL_CATEGORIES), dtype=object)
    pack[TOTAL_RATE] = blk
    return blk


def load_category_syst_summary(
    npz_path: Optional[str] = None,
    *,
    syst_disk_root: Optional[str] = None,
) -> Dict[str, Any]:
    """Load summary NPZ; returns ``{'manifest', 'by_var', 'npz_path'}``."""
    if npz_path is None:
        root = syst_disk_root or os.environ.get(SYST_DISK_ENV)
        if not root:
            raise ValueError(
                "Pass npz_path or syst_disk_root, or set %s" % SYST_DISK_ENV
            )
        npz_path = category_summary_npz_path(root)
    npz_path = os.path.abspath(npz_path)
    blob = np.load(npz_path, allow_pickle=True)
    by_var = {k: blob[k].item() for k in blob.files}
    manifest_path = category_summary_manifest_path(npz_path)
    manifest: Dict[str, Any] = {}
    if os.path.isfile(manifest_path):
        with open(manifest_path, encoding="utf-8") as mf:
            manifest = json.load(mf)
    return {"npz_path": npz_path, "manifest": manifest, "by_var": by_var}


def _var_pack(summary: Mapping[str, Any], var_save_name: str) -> Mapping[str, Any]:
    by_var = summary["by_var"]
    if var_save_name not in by_var:
        raise KeyError(
            "Variable %r not in category summary (have: %s)"
            % (var_save_name, ", ".join(sorted(by_var.keys())))
        )
    return by_var[var_save_name]


def category_cov_frac(
    summary: Mapping[str, Any], var_save_name: str, category_key: str
) -> np.ndarray:
    pack = _var_pack(summary, var_save_name)
    return np.asarray(pack["categories"][category_key]["cov_frac"], dtype=np.float64)


def category_frac_unc_pct(
    summary: Mapping[str, Any], var_save_name: str, category_key: str
) -> np.ndarray:
    pack = _var_pack(summary, var_save_name)
    return np.asarray(pack["categories"][category_key]["frac_unc_pct"], dtype=np.float64)


def total_cov_frac(
    summary: Mapping[str, Any], var_save_name: str, *, kind: str = "rate"
) -> np.ndarray:
    key = TOTAL_XSEC if kind == "xsec" else TOTAL_RATE
    pack = _var_pack(summary, var_save_name)
    if key not in pack:
        raise KeyError("No %r entry for variable %r" % (key, var_save_name))
    return np.asarray(pack[key]["cov_frac"], dtype=np.float64)


def total_frac_unc_pct(
    summary: Mapping[str, Any], var_save_name: str, *, kind: str = "rate"
) -> np.ndarray:
    key = TOTAL_XSEC if kind == "xsec" else TOTAL_RATE
    pack = _var_pack(summary, var_save_name)
    if key not in pack:
        raise KeyError("No %r entry for variable %r" % (key, var_save_name))
    return np.asarray(pack[key]["frac_unc_pct"], dtype=np.float64)


def combine_categories_cov_frac(
    summary: Mapping[str, Any],
    var_save_name: str,
    category_keys: Sequence[str],
) -> np.ndarray:
    """Sum category fractional covariances (independent sources)."""
    mats = [category_cov_frac(summary, var_save_name, k) for k in category_keys]
    out = sum_cov_frac_matrices(mats)
    if out is None:
        raise ValueError("No category matrices to combine")
    return out
