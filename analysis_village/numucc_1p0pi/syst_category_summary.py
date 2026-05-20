"""Export / load per-category fractional systematics for plotting and unfolding.

Written by ``systematics-summary.ipynb`` into ``<SYST_DISK_ROOT>/CategorySummary/``.
Each variable entry stores fractional covariance matrices and per-bin uncertainties in percent,
using the same conventions as the summary breakdown plots (cosmics = contamination-scaled
``SelectedRate``; GENIE rate and xsec totals available separately).

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

from analysis_village.numucc_1p0pi.syst_disk_layout import (
    SYST_DISK_ENV,
    category_summary_dir,
    category_summary_manifest_path,
    category_summary_npz_path,
    normalized_root,
)

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

# Flat correlated terms (percent, matching ``frac_unc_pct`` scale in summary plots).
POT_FRAC_UNC_PCT = 2.0
NTARGETS_FRAC_UNC_PCT = 1.0

_RATE_TOTAL_CATEGORIES = (
    CAT_FLUX,
    CAT_G4,
    CAT_MCSTAT,
    CAT_DETECTOR,
    CAT_COSMICS,
    CAT_GENIE_RATE,
)
_XSEC_TOTAL_CATEGORIES = (
    CAT_FLUX,
    CAT_G4,
    CAT_MCSTAT,
    CAT_DETECTOR,
    CAT_COSMICS,
    CAT_GENIE_XSEC,
)


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
    return genie_blob[var_name]


def genie_knob_covs(gd: Optional[Mapping]) -> Optional[Dict[str, Any]]:
    if gd is None:
        return None
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
    u = float(frac_unc_pct_val) / 100.0
    return np.diag(np.full(nbins, u * u, dtype=np.float64))


def _category_block(cov_frac: np.ndarray, var_config: Any) -> Dict[str, np.ndarray]:
    cov_frac = np.asarray(cov_frac, dtype=np.float64)
    w = frac_weights_for_plot(cov_frac, var_config)
    return {
        "cov_frac": cov_frac,
        "frac_unc_pct": w,
        "integrated_frac_variance": np.array(integrated_rate_frac_variance(cov_frac)),
    }


def build_category_cov_frac(
    vsn: str,
    *,
    flux_npz: Any,
    g4_npz: Any,
    cosmics_npz: Any,
    detector_npz: Any = None,
    mcstat_npz: Any = None,
    genie_blob: Optional[Mapping] = None,
) -> Dict[str, np.ndarray]:
    """Return ``{category_key: cov_frac}`` for one variable slug."""
    covs: Dict[str, np.ndarray] = {
        CAT_FLUX: np.asarray(flux_npz[vsn].item()["flux"]["cov_frac"], dtype=np.float64),
        CAT_G4: np.asarray(g4_npz[vsn].item()["G4"]["cov_frac"], dtype=np.float64),
    }
    if mcstat_npz is not None and vsn in mcstat_npz:
        covs[CAT_MCSTAT] = np.asarray(
            mcstat_npz[vsn].item()["MCstat"]["cov_frac"], dtype=np.float64
        )
    if detector_npz is not None:
        det_item = dict(detector_npz)["detector"].item()
        if vsn in det_item:
            covs[CAT_DETECTOR] = np.asarray(det_item[vsn]["cov_frac"], dtype=np.float64)
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
    mcstat_npz: Any = None,
    genie_blob: Optional[Mapping] = None,
    include_flat: bool = True,
) -> Dict[str, Any]:
    """Full export payload for one ``VariableConfig``."""
    vsn = var_config.var_save_name
    nbins = len(var_config.bin_centers)
    covs = build_category_cov_frac(
        vsn,
        flux_npz=flux_npz,
        g4_npz=g4_npz,
        cosmics_npz=cosmics_npz,
        detector_npz=detector_npz,
        mcstat_npz=mcstat_npz,
        genie_blob=genie_blob,
    )
    if include_flat:
        covs[CAT_POT] = _flat_cov_frac(nbins, POT_FRAC_UNC_PCT)
        covs[CAT_NTARGETS] = _flat_cov_frac(nbins, NTARGETS_FRAC_UNC_PCT)

    categories: Dict[str, Dict[str, np.ndarray]] = {
        k: _category_block(c, var_config) for k, c in covs.items()
    }

    def _total_block(keys: Sequence[str]) -> Dict[str, np.ndarray]:
        parts = [covs[k] for k in keys if k in covs]
        if not parts:
            return {}
        total = sum_cov_frac_matrices(parts)
        assert total is not None
        blk = _category_block(total, var_config)
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
    rate_blk = _total_block(_RATE_TOTAL_CATEGORIES)
    if rate_blk:
        pack[TOTAL_RATE] = rate_blk
    xsec_blk = _total_block(_XSEC_TOTAL_CATEGORIES)
    if xsec_blk:
        pack[TOTAL_XSEC] = xsec_blk
    return pack


def export_category_syst_summary(
    out_npz: str,
    var_configs: Sequence[Any],
    *,
    flux_npz: Any,
    g4_npz: Any,
    cosmics_npz: Any,
    detector_npz: Any = None,
    mcstat_npz: Any = None,
    genie_blob: Optional[Mapping] = None,
    syst_disk_root: Optional[str] = None,
    include_flat: bool = True,
) -> Dict[str, Any]:
    """Write ``category_syst_summary.npz`` and companion manifest JSON."""
    out_npz = os.path.abspath(out_npz)
    os.makedirs(os.path.dirname(out_npz), exist_ok=True)
    manifest_path = category_summary_manifest_path(out_npz)

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
                mcstat_npz=mcstat_npz,
                genie_blob=genie_blob,
                include_flat=include_flat,
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
        "usage": {
            "load": "load_category_syst_summary(npz_path) or load_category_syst_summary(syst_disk_root=...)",
            "cov_frac": "category_cov_frac(summary, var_save_name, category_key)",
            "frac_unc_pct": "category_frac_unc_pct(summary, var_save_name, category_key)",
            "total": "total_frac_unc_pct(summary, var_save_name, kind='xsec'|'rate')",
        },
    }
    with open(manifest_path, "w", encoding="utf-8") as mf:
        json.dump(manifest, mf, indent=2)
        mf.write("\n")

    return manifest


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
    summary: Mapping[str, Any], var_save_name: str, *, kind: str = "xsec"
) -> np.ndarray:
    key = TOTAL_XSEC if kind == "xsec" else TOTAL_RATE
    pack = _var_pack(summary, var_save_name)
    if key not in pack:
        raise KeyError("No %r entry for variable %r" % (key, var_save_name))
    return np.asarray(pack[key]["cov_frac"], dtype=np.float64)


def total_frac_unc_pct(
    summary: Mapping[str, Any], var_save_name: str, *, kind: str = "xsec"
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
