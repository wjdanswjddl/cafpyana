#!/usr/bin/env python3
"""
Batch-process BNB gsimple ROOT flux files like ``gsimple.ipynb``, without loading
all rays into memory.

The notebook concatenates every file's arrays, which can exhaust RAM on multi-thousand-file
sets. This script processes ``--batch-size`` files at a time, adds weighted energy
histograms (and ray counts) into accumulators, then divides by ``sum(POT)`` once at the end.

Physics matches the notebook's intended normalization: spectra are ``sum(wgt in bin) /
(area * TOTAL_POT)`` with ``TOTAL_POT = sum_k POT_k`` over processed files.

**Outputs** (when ``--out-dir`` and/or ``--save-npz`` is set):

* ``flux_histograms.npz`` — arrays ``hist_weighted_zwf`` and ``flux_density_zwf`` with shape
  ``(n_z, n_window, n_flavor, n_energy)`` for every combination of projection plane ``z``,
  rectangular flux window, and neutrino flavor; plus ``n_rays_zwf``, bin edges, and window
  metadata. Ray-traced FV sub-volumes add ``hist_weighted_raytrace_wvf`` /
  ``flux_density_raytrace_wvf`` with shape ``(n_raytrace_volume, n_flavor, n_energy)``.
* ``flux_histograms_meta.json`` — run configuration and window list.
* ``raytrace_flux_histograms.npz`` — same ray-traced FV volume spectra as
  ``gsimple_raytrace_batch.py`` (written next to ``flux_histograms.npz`` for cross-checks).
* ``slab_zconfig_flux_histograms.npz`` — z-slab flux on several uniform $z$ grids
  (``--slab-z-ns``), with volume-averaged spectra for each ray-traced FV comparison volume.
* ``slab_zconfig_flux_histograms_meta.json`` — slab grid definitions and compare-volume mapping.
* ``processed_files.txt`` / ``skipped_files.json`` — per-file success list and skip reasons
  (unless ``--no-manifest``).
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path

import matplotlib

# Batch-friendly default; override with MPLBACKEND or GSIMPLE_MPL_BACKEND (e.g. TkAgg) for plt.show().
matplotlib.use(os.environ.get("GSIMPLE_MPL_BACKEND", os.environ.get("MPLBACKEND", "Agg")))
import matplotlib.pyplot as plt
import numpy as np
import uproot
from tqdm.auto import tqdm

# ---------------------------------------------------------------------------
# Defaults (same as gsimple.ipynb)
# ---------------------------------------------------------------------------

DEFAULT_GSIMPLE_DIR = (
    "/cvmfs/sbnd.osgstorage.org/pnfs/fnal.gov/usr/sbnd/persistent/stash/"
    "fluxFiles/bnb/BooNEtoGSimple/configK-v1/july2023/neutrinoMode/"
)

Z_POSITIONS_CM = np.linspace(0.0, 500.0, 26)

SBND_VOLUME = dict(z_range_cm=(0.0, 500.0), y_range=(-200.0, 200.0), x_range=(-200.0, 200.0))
FLUX_WINDOW_XY = dict(x_range=(-200.0, 200.0), y_range=(-200.0, 200.0))

FACE_Z_CM = 0.0
DOWNSTREAM_FACE_Z_CM = FACE_Z_CM + 500.0

PDG_BY_FLAVOR = {
    "nue": 12,
    "nuebar": -12,
    "numu": 14,
    "numubar": -14,
}
FLAVORS = list(PDG_BY_FLAVOR.keys())
FLAVOR_LATEX = {
    "nue": r"$\nu_e$",
    "nuebar": r"$\bar{\nu}_e$",
    "numu": r"$\nu_\mu$",
    "numubar": r"$\bar{\nu}_\mu$",
}

FACE_WINDOWS = {
    r"AV ($\pm$200, $\pm$200)": dict(x_range=(-200.0, 200.0), y_range=(-200.0, 200.0)),
    r"FV ($\pm$190, $\pm$190)": dict(x_range=(-190.0, 190.0), y_range=(-190.0, 190.0)),
    r"FV 2 ($\pm$190, [-190,100])": dict(x_range=(-190.0, 190.0), y_range=(-190.0, 100.0)),
}

# Rectangular (x, y) flux windows [cm], in fixed order for the z × window × flavor tensor.
# Index 0 matches SBND_AV / default dk2nu-style flux window (same geometry as ``FLUX_WINDOW_XY``).
FLUX_WINDOW_DEFS: list[tuple[str, tuple[float, float], tuple[float, float]]] = [
    ("AV_pm200_cm", (-200.0, 200.0), (-200.0, 200.0)),
    ("FV_pm190_cm", (-190.0, 190.0), (-190.0, 190.0)),
    ("FV2_pm190_y-190_100_cm", (-190.0, 190.0), (-190.0, 100.0)),
]

# Map FACE_WINDOWS plot labels to FLUX_WINDOW_DEFS index (geometries match in same order).
_FACE_LABELS_ORDERED = list(FACE_WINDOWS.keys())
assert len(_FACE_LABELS_ORDERED) == len(FLUX_WINDOW_DEFS)
FACE_LABEL_TO_WINDOW_IDX = {lbl: i for i, lbl in enumerate(_FACE_LABELS_ORDERED)}

CONFIGS = {
    "SBND volume (avg over z)": ("volume", None),
    "upstream face": ("plane", (FACE_Z_CM, -200.0, 200.0, -200.0, 200.0)),
    "downstream face": ("plane", (DOWNSTREAM_FACE_Z_CM, -200.0, 200.0, -200.0, 200.0)),
}

# Ray-traced 3D volumes (canonical definitions in raytrace_volume_defs.py / gsimple_raytrace.ipynb)
from raytrace_volume_defs import (
    BIN_CENTERS,
    E_BINS,
    FV_SPLIT_BOXES,
    FV_SPLIT_TRUNCY_BOXES,
    N_E_BINS,
    RAYTRACE_VOLUME_DEFS,
    RAYTRACE_VOLUME_LABEL,
    SLAB_COMPARE_Z_RANGE_CM,
    SLAB_Z_RANGE_FV_CM,
    mask_xy_in_volume_at_z,
    plane_area_m2_at_z,
)

# Default z-slab grids for convergence studies (uniform centers on FV z extent).
DEFAULT_SLAB_Z_NS = [5, 11, 21, 26, 51, 101, 201]


def build_slab_z_grid_defs(
    n_slabs_list: list[int],
    *,
    z_range_fv_cm: tuple[float, float] = SLAB_Z_RANGE_FV_CM,
    include_batch_grid: bool = True,
) -> list[dict]:
    """Build slab z-grid configs: uniform linspace grids plus optional batch default grid."""
    z_lo, z_hi = z_range_fv_cm
    defs: list[dict] = []
    for n in n_slabs_list:
        n = int(n)
        if n < 2:
            raise ValueError(f"slab z grid needs n >= 2, got {n}")
        z_cm = np.linspace(z_lo, z_hi, n)
        defs.append(
            {
                "key": f"z{z_lo:g}_{z_hi:g}_n{n}",
                "z_cm": z_cm,
                "z_range_cm": z_range_fv_cm,
                "dz_cm": float((z_hi - z_lo) / (n - 1)),
                "n_z": n,
                "batch_default": False,
            }
        )
    if include_batch_grid:
        z_cm = np.linspace(0.0, 500.0, 26)
        defs.append(
            {
                "key": "z0_500_n26_batch",
                "z_cm": z_cm,
                "z_range_cm": (0.0, 500.0),
                "dz_cm": 20.0,
                "n_z": 26,
                "batch_default": True,
            }
        )
    return defs


def pack_hist_compare_czvf_from_union(
    hist_union_czvf: np.ndarray,
    slab_iz_maps: list[np.ndarray],
) -> np.ndarray:
    """Slice union (n_z, n_compare_vol, f, E) into padded (n_config, n_z_max, n_cv, f, E)."""
    n_cfg = len(slab_iz_maps)
    n_z_max = max(len(iz) for iz in slab_iz_maps)
    n_cv, n_f, n_e = hist_union_czvf.shape[1:]
    hist_czvf = np.zeros((n_cfg, n_z_max, n_cv, n_f, n_e))
    for icfg, iz_map in enumerate(slab_iz_maps):
        hist_czvf[icfg, : len(iz_map)] = hist_union_czvf[iz_map]
    return hist_czvf


def finalize_slab_zconfig_arrays(
    hist_czwf: np.ndarray,
    hist_compare_czvf: np.ndarray,
    slab_defs: list[dict],
    *,
    total_pot: float,
    window_area_m2_arr: np.ndarray,
    compare_volume_defs: list[tuple[str, list[dict]]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return (flux_density_czwf, flux_density_mean_cvwf, z_positions_cm_cz).

    flux_density_czwf: (n_config, n_z_max, n_window, n_flavor, n_E) — legacy rectangular windows
    flux_density_mean_cvwf: (n_config, n_compare_vol, n_flavor, n_E) — z-mean with raytrace geometry
    z_positions_cm_cz: (n_config, n_z_max), NaN-padded
    """
    n_config, n_z_max, n_windows, n_flavors, n_e = hist_czwf.shape
    inv_apot = 1.0 / (window_area_m2_arr[np.newaxis, :, np.newaxis, np.newaxis] * total_pot)
    flux_czwf = hist_czwf * inv_apot

    n_cv = len(compare_volume_defs)
    flux_mean_cvwf = np.zeros((n_config, n_cv, n_flavors, n_e))

    z_positions_cm_cz = np.full((n_config, n_z_max), np.nan, dtype=np.float64)
    for icfg, sdef in enumerate(slab_defs):
        nz = len(sdef["z_cm"])
        z_positions_cm_cz[icfg, :nz] = sdef["z_cm"]
        z_cm = sdef["z_cm"]
        for iv, (vkey, boxes) in enumerate(compare_volume_defs):
            z_lo, z_hi = SLAB_COMPARE_Z_RANGE_CM[vkey]
            iz = np.where((z_cm >= z_lo) & (z_cm <= z_hi))[0]
            if iz.size == 0:
                raise ValueError(
                    f"No z slabs of config {sdef['key']!r} in [{z_lo}, {z_hi}] for {vkey}"
                )
            plane_flux = []
            for iz in iz:
                z = float(z_cm[iz])
                area_m2 = plane_area_m2_at_z(z, boxes)
                if area_m2 <= 0.0:
                    continue
                plane_flux.append(hist_compare_czvf[icfg, iz, iv] / (area_m2 * total_pot))
            if not plane_flux:
                raise ValueError(f"No active z planes for {vkey} in config {sdef['key']!r}")
            flux_mean_cvwf[icfg, iv] = np.mean(np.stack(plane_flux, axis=0), axis=0)

    return flux_czwf, flux_mean_cvwf, z_positions_cm_cz


def write_slab_zconfig_npz(
    path: Path,
    *,
    total_pot: float,
    slab_defs: list[dict],
    hist_czwf: np.ndarray,
    flux_density_czwf: np.ndarray,
    flux_density_mean_cvwf: np.ndarray,
    z_positions_cm_cz: np.ndarray,
    window_area_m2_arr: np.ndarray,
) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    config_keys = np.array([d["key"] for d in slab_defs], dtype=object)
    config_n_z = np.array([len(d["z_cm"]) for d in slab_defs], dtype=np.int32)
    config_dz_cm = np.array([d["dz_cm"] for d in slab_defs], dtype=np.float64)
    config_batch_default = np.array([bool(d.get("batch_default", False)) for d in slab_defs], dtype=bool)
    compare_volume_keys = np.array([v for v, _ in RAYTRACE_VOLUME_DEFS], dtype=object)
    window_keys = np.array([t[0] for t in FLUX_WINDOW_DEFS], dtype=object)

    np.savez_compressed(
        path,
        total_pot=np.array(total_pot),
        E_BINS=E_BINS,
        BIN_CENTERS=BIN_CENTERS,
        flavors=np.array(FLAVORS, dtype=object),
        window_keys=window_keys,
        window_area_m2=window_area_m2_arr,
        slab_config_keys=config_keys,
        slab_config_n_z=config_n_z,
        slab_config_dz_cm=config_dz_cm,
        slab_config_batch_default=config_batch_default,
        z_positions_cm_cz=z_positions_cm_cz,
        compare_volume_keys=compare_volume_keys,
        hist_weighted_czwf=hist_czwf,
        flux_density_czwf=flux_density_czwf,
        flux_density_mean_cvwf=flux_density_mean_cvwf,
    )


def _mplstyle() -> None:
    style = Path(__file__).resolve().parent.parent / "numucc_1p0pi" / "presentation.mplstyle"
    try:
        plt.style.use(style)
    except OSError:
        pass


def _get_pot(f: uproot.reading.ReadOnlyDirectory) -> float:
    meta_key = next((k for k in f.keys() if "meta" in k.lower()), None)
    if meta_key is None:
        raise RuntimeError("No meta tree found; cannot read POT.")
    meta = f[meta_key]
    pot_key = next((k for k in meta.keys() if "proton" in k.lower()), None)
    if pot_key is None:
        raise RuntimeError(f"No 'protons' branch found in meta tree; keys = {list(meta.keys())}")
    return float(meta[pot_key].array(library="np").sum())


def load_gsimple(filename: str) -> dict:
    """Load one file's ray arrays (cm for vertices). See gsimple.ipynb."""
    with uproot.open(filename) as f:
        pot = _get_pot(f)
        tree = f["flux"]
        keys = tree.keys()
        if "entry/pdg" in keys:
            prefix = "entry/"
        elif "pdg" in keys:
            prefix = ""
        else:
            raise RuntimeError(f"Cannot find 'pdg' branch in {filename}; keys = {list(keys)}")

        branches = [
            prefix + "pdg",
            prefix + "wgt",
            prefix + "vtxx",
            prefix + "vtxy",
            prefix + "vtxz",
            prefix + "px",
            prefix + "py",
            prefix + "pz",
            prefix + "E",
        ]
        arr = tree.arrays(branches, library="np")

    return {
        "pot": pot,
        "pdg": arr[prefix + "pdg"],
        "wgt": arr[prefix + "wgt"],
        "vtxx": arr[prefix + "vtxx"] * 100.0,
        "vtxy": arr[prefix + "vtxy"] * 100.0,
        "vtxz": arr[prefix + "vtxz"] * 100.0,
        "px": arr[prefix + "px"],
        "py": arr[prefix + "py"],
        "pz": arr[prefix + "pz"],
        "E": arr[prefix + "E"],
    }


def window_area_m2(x_range: tuple[float, float], y_range: tuple[float, float]) -> float:
    return (x_range[1] - x_range[0]) * 1.0e-2 * (y_range[1] - y_range[0]) * 1.0e-2


def project_to_z(data: dict, z_target_cm: float) -> tuple[np.ndarray, np.ndarray]:
    slope_x = data["px"] / data["pz"]
    slope_y = data["py"] / data["pz"]
    dz = z_target_cm - data["vtxz"]
    x_cm = data["vtxx"] + slope_x * dz
    y_cm = data["vtxy"] + slope_y * dz
    return x_cm, y_cm


def mask_xy(
    x_cm: np.ndarray,
    y_cm: np.ndarray,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
) -> np.ndarray:
    return (
        (x_cm >= x_range[0])
        & (x_cm < x_range[1])
        & (y_cm >= y_range[0])
        & (y_cm < y_range[1])
    )


def add_weighted_e_hist(
    data: dict,
    flavor: str,
    z_cm: float,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    out: np.ndarray,
) -> int:
    """Add np.histogram(..., weights=wgt) for rays in the window; return ray count."""
    x_cm, y_cm = project_to_z(data, z_cm)
    pdg = PDG_BY_FLAVOR[flavor]
    m = (data["pdg"] == pdg) & mask_xy(x_cm, y_cm, x_range, y_range)
    h, _ = np.histogram(data["E"][m], bins=E_BINS, weights=data["wgt"][m])
    out += h
    return int(np.sum(m))


def build_z_union_maps(
    z_slabs_default: list[float] | np.ndarray,
    slab_defs: list[dict],
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """
    Merge default and slab-sweep z grids into one sorted union.

    Returns (z_union_cm, iz_default, slab_iz_maps) where iz_default indexes z_union
  for the legacy 0–500 cm / n=26 grid, and each slab_iz_maps[icfg] indexes planes
    for that config — each physical z is projected only once per ROOT file.
    """
    pieces = [np.asarray(z_slabs_default, dtype=np.float64)]
    for sdef in slab_defs:
        pieces.append(np.asarray(sdef["z_cm"], dtype=np.float64))
    z_union = np.unique(np.concatenate(pieces))
    z_union.sort()

    def _indices(z_query: np.ndarray) -> np.ndarray:
        z_query = np.asarray(z_query, dtype=np.float64)
        iz = np.searchsorted(z_union, z_query)
        if not np.all(np.isclose(z_union[iz], z_query, rtol=0, atol=1e-9)):
            raise ValueError("z grid values not found in union (non-linspace z?)")
        return iz

    iz_default = _indices(np.asarray(z_slabs_default, dtype=np.float64))
    slab_iz_maps = [_indices(np.asarray(sdef["z_cm"], dtype=np.float64)) for sdef in slab_defs]
    return z_union, iz_default, slab_iz_maps


def pack_hist_czwf_from_union(
    hist_union_zwf: np.ndarray,
    slab_iz_maps: list[np.ndarray],
) -> np.ndarray:
    """Slice union histogram (n_z, w, f, E) into padded (n_config, n_z_max, w, f, E)."""
    n_cfg = len(slab_iz_maps)
    n_z_max = max(len(iz) for iz in slab_iz_maps)
    n_w, n_f, n_e = hist_union_zwf.shape[1:]
    hist_czwf = np.zeros((n_cfg, n_z_max, n_w, n_f, n_e))
    for icfg, iz_map in enumerate(slab_iz_maps):
        hist_czwf[icfg, : len(iz_map)] = hist_union_zwf[iz_map]
    return hist_czwf


def accumulate_z_histograms_for_file(
    data: dict,
    z_union_cm: np.ndarray,
    hist_zwf: np.ndarray,
    n_rays_zwf: np.ndarray,
    *,
    hist_compare_czvf: np.ndarray | None = None,
    n_rays_compare_czvf: np.ndarray | None = None,
    compare_volume_defs: list[tuple[str, list[dict]]] | None = None,
) -> None:
    """Project each z plane once per file; fill rectangular windows and optional FV slab volumes."""
    e = data["E"]
    wgt = data["wgt"]
    pdg = data["pdg"]
    tally_compare = hist_compare_czvf is not None
    if tally_compare:
        if compare_volume_defs is None or n_rays_compare_czvf is None:
            raise ValueError("compare_volume_defs and n_rays_compare_czvf required for FV slab tally")
    for iz, z in enumerate(z_union_cm):
        z_f = float(z)
        x_cm, y_cm = project_to_z(data, z_f)
        for iw, (_name, xr, yr) in enumerate(FLUX_WINDOW_DEFS):
            in_window = mask_xy(x_cm, y_cm, xr, yr)
            for iflav, flav in enumerate(FLAVORS):
                m = (pdg == PDG_BY_FLAVOR[flav]) & in_window
                if not np.any(m):
                    continue
                h, _ = np.histogram(e[m], bins=E_BINS, weights=wgt[m])
                hist_zwf[iz, iw, iflav] += h
                n_rays_zwf[iz, iw, iflav] += int(np.count_nonzero(m))
        if tally_compare:
            for ivol, (_vname, boxes) in enumerate(compare_volume_defs):
                in_vol = mask_xy_in_volume_at_z(x_cm, y_cm, z_f, boxes)
                for iflav, flav in enumerate(FLAVORS):
                    m = (pdg == PDG_BY_FLAVOR[flav]) & in_vol
                    if not np.any(m):
                        continue
                    h, _ = np.histogram(e[m], bins=E_BINS, weights=wgt[m])
                    hist_compare_czvf[iz, ivol, iflav] += h
                    n_rays_compare_czvf[iz, ivol, iflav] += int(np.count_nonzero(m))


def _slab_t(v: np.ndarray, d: np.ndarray, lo: float, hi: float) -> tuple[np.ndarray, np.ndarray]:
    with np.errstate(divide="ignore", invalid="ignore"):
        t1 = (lo - v) / d
        t2 = (hi - v) / d
    t_near = np.minimum(t1, t2)
    t_far = np.maximum(t1, t2)
    parallel = d == 0
    if np.any(parallel):
        inside = (v >= lo) & (v <= hi)
        t_near = np.where(parallel, np.where(inside, -np.inf, np.inf), t_near)
        t_far = np.where(parallel, np.where(inside, np.inf, -np.inf), t_far)
    return t_near, t_far


def path_length_box(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    dx: np.ndarray,
    dy: np.ndarray,
    dz: np.ndarray,
    box: dict,
) -> np.ndarray:
    """Path length [cm] through one axis-aligned box, clipped to forward ray (t >= 0)."""
    tx_n, tx_f = _slab_t(vx, dx, *box["x_range"])
    ty_n, ty_f = _slab_t(vy, dy, *box["y_range"])
    tz_n, tz_f = _slab_t(vz, dz, *box["z_range"])
    t_enter = np.maximum.reduce([tx_n, ty_n, tz_n, np.zeros_like(vx)])
    t_exit = np.minimum.reduce([tx_f, ty_f, tz_f])
    return np.maximum(0.0, t_exit - t_enter)


def path_length_volume(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    dx: np.ndarray,
    dy: np.ndarray,
    dz: np.ndarray,
    boxes: list[dict],
) -> np.ndarray:
    L = np.zeros_like(vx)
    for box in boxes:
        L += path_length_box(vx, vy, vz, dx, dy, dz, box)
    return L


def volume_cm3(boxes: list[dict]) -> float:
    v = 0.0
    for b in boxes:
        dx = b["x_range"][1] - b["x_range"][0]
        dy = b["y_range"][1] - b["y_range"][0]
        dz = b["z_range"][1] - b["z_range"][0]
        v += dx * dy * dz
    return v


def add_raytrace_weighted_hist(
    data: dict,
    flavor: str,
    path_lengths: np.ndarray,
    out: np.ndarray,
) -> int:
    """Add sum(wgt * L) per energy bin for rays with L > 0; return ray count."""
    pdg = PDG_BY_FLAVOR[flavor]
    m = (data["pdg"] == pdg) & (path_lengths > 0)
    h, _ = np.histogram(data["E"][m], bins=E_BINS, weights=data["wgt"][m] * path_lengths[m])
    out += h
    return int(np.sum(m))


def integrated_flux(spectrum: np.ndarray) -> float:
    return float(np.sum(spectrum))


def write_raytrace_npz(
    path: Path,
    *,
    total_pot: float,
    hist_wvf: np.ndarray,
    n_rays_wvf: np.ndarray,
    flux_density_wvf: np.ndarray,
    volume_cm3_arr: np.ndarray,
) -> None:
    """Write ray-traced volume histograms in the same schema as ``gsimple_raytrace_batch.py``."""
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    volume_keys = np.array([t[0] for t in RAYTRACE_VOLUME_DEFS], dtype=object)
    np.savez_compressed(
        path,
        total_pot=np.array(total_pot),
        E_BINS=E_BINS,
        BIN_CENTERS=BIN_CENTERS,
        volume_keys=volume_keys,
        volume_cm3=volume_cm3_arr,
        flavors=np.array(FLAVORS, dtype=object),
        hist_weighted_wvf=hist_wvf,
        n_rays_wvf=n_rays_wvf,
        flux_density_wvf=flux_density_wvf,
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--gsimple-dir",
        type=str,
        default=DEFAULT_GSIMPLE_DIR,
        help="Directory containing gsimple *.root files",
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=500,
        help="Number of ROOT files to load before freeing their arrays (default: 500)",
    )
    p.add_argument(
        "--n-files",
        type=int,
        default=None,
        help="Process at most this many files (default: all)",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="If set, save figures under this directory instead of plt.show()",
    )
    p.add_argument(
        "--save-npz",
        type=Path,
        default=None,
        help="Write z×window×flavor histograms (+metadata) to this .npz path. "
        "If omitted but --out-dir is set, defaults to OUT_DIR/flux_histograms.npz",
    )
    p.add_argument(
        "--save-raytrace-npz",
        type=Path,
        default=None,
        help="Write ray-traced FV volume histograms (gsimple_raytrace_batch schema) to this path. "
        "If omitted but flux_histograms.npz is written, defaults to "
        "OUT_DIR/raytrace_flux_histograms.npz",
    )
    p.add_argument(
        "--no-manifest",
        action="store_true",
        help="Do not write processed_files.txt / skipped_files.json next to the histogram .npz",
    )
    p.add_argument(
        "--slab-z-ns",
        type=int,
        nargs="+",
        default=None,
        help=(
            "Uniform z-slab grids on [10,450] cm with this many centers each "
            f"(default: {DEFAULT_SLAB_Z_NS}). Written to slab_zconfig_flux_histograms.npz."
        ),
    )
    p.add_argument(
        "--no-batch-slab-grid",
        action="store_true",
        help="Do not add the legacy linspace(0,500,26) grid to the slab z-config sweep.",
    )
    p.add_argument(
        "--no-slab-zconfig-npz",
        action="store_true",
        help="Skip writing slab_zconfig_flux_histograms.npz (saves memory/time).",
    )
    p.add_argument(
        "--save-slab-zconfig-npz",
        type=Path,
        default=None,
        help="Path for slab z-config histograms. Default: OUT_DIR/slab_zconfig_flux_histograms.npz",
    )
    p.add_argument("--explore-only", action="store_true", help="Print one file's structure and exit")
    p.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip diagnostic PNG figures at end of run (recommended for batch jobs)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    _mplstyle()

    all_files = sorted(glob.glob(os.path.join(args.gsimple_dir, "*.root")))
    if not all_files:
        raise SystemExit(f"No *.root files under {args.gsimple_dir!r}")

    if args.explore_only:
        sample = all_files[0]
        print(f"Found {len(all_files)} files")
        print("Sample:", sample)
        with uproot.open(sample) as f:
            print("Keys:", f.keys())
            print("flux branches:", list(f["flux"].keys()))
        return

    if args.n_files is not None:
        all_files = all_files[: args.n_files]

    z_lo, z_hi = SBND_VOLUME["z_range_cm"]
    z_slabs = list(Z_POSITIONS_CM)
    n_z_slabs = len(z_slabs)
    z_volume_indices = [iz for iz, z in enumerate(z_slabs) if z_lo <= z <= z_hi]
    IW_AV = 0  # AV ±200 cm; same geometry as ``FLUX_WINDOW_XY`` / SBND AV plots

    total_pot = 0.0
    n_files_ok = 0

    n_windows = len(FLUX_WINDOW_DEFS)
    n_flavors = len(FLAVORS)
    n_raytrace_volumes = len(RAYTRACE_VOLUME_DEFS)

    write_slab_zconfig = not args.no_slab_zconfig_npz
    slab_defs: list[dict] = []
    slab_iz_maps: list[np.ndarray] = []
    if write_slab_zconfig:
        z_ns = args.slab_z_ns if args.slab_z_ns is not None else DEFAULT_SLAB_Z_NS
        slab_defs = build_slab_z_grid_defs(
            z_ns,
            include_batch_grid=not args.no_batch_slab_grid,
        )

    z_union, iz_default, slab_iz_maps_sweep = build_z_union_maps(z_slabs, slab_defs)
    if write_slab_zconfig:
        slab_iz_maps = slab_iz_maps_sweep
        n_redundant = n_z_slabs + sum(len(d["z_cm"]) for d in slab_defs)
        print(
            f"Slab z-config sweep: {len(slab_defs)} grids | "
            f"unique z planes per file: {len(z_union)} (deduped from {n_redundant})"
        )

    n_z_union = len(z_union)
    hist_union_zwf = np.zeros((n_z_union, n_windows, n_flavors, N_E_BINS))
    n_rays_union_zwf = np.zeros((n_z_union, n_windows, n_flavors), dtype=np.int64)

    n_compare_volumes = len(RAYTRACE_VOLUME_DEFS)
    hist_compare_union_czvf = (
        np.zeros((n_z_union, n_compare_volumes, n_flavors, N_E_BINS))
        if write_slab_zconfig
        else None
    )
    n_rays_compare_union_czvf = (
        np.zeros((n_z_union, n_compare_volumes, n_flavors), dtype=np.int64)
        if write_slab_zconfig
        else None
    )

    # volume × flavor × energy: raw sum of wgt * path_length [cm] per E bin
    hist_raytrace_wvf = np.zeros((n_raytrace_volumes, n_flavors, N_E_BINS))
    n_rays_raytrace_wvf = np.zeros((n_raytrace_volumes, n_flavors), dtype=np.int64)

    ok_files: list[str] = []
    skipped: list[dict[str, str]] = []

    def process_one_file(fpath: str) -> None:
        nonlocal total_pot, n_files_ok
        d = load_gsimple(fpath)
        pot = float(d["pot"])

        accumulate_z_histograms_for_file(
            d,
            z_union,
            hist_union_zwf,
            n_rays_union_zwf,
            hist_compare_czvf=hist_compare_union_czvf,
            n_rays_compare_czvf=n_rays_compare_union_czvf,
            compare_volume_defs=RAYTRACE_VOLUME_DEFS if write_slab_zconfig else None,
        )

        p_mag = np.sqrt(d["px"] ** 2 + d["py"] ** 2 + d["pz"] ** 2)
        dx = d["px"] / p_mag
        dy = d["py"] / p_mag
        dz = d["pz"] / p_mag
        for ivol, (_vname, boxes) in enumerate(RAYTRACE_VOLUME_DEFS):
            path_L = path_length_volume(d["vtxx"], d["vtxy"], d["vtxz"], dx, dy, dz, boxes)
            for iflav, flav in enumerate(FLAVORS):
                n_rays_raytrace_wvf[ivol, iflav] += add_raytrace_weighted_hist(
                    d, flav, path_L, hist_raytrace_wvf[ivol, iflav]
                )

        total_pot += pot
        n_files_ok += 1
        ok_files.append(fpath)

    # Batch loop: load at most batch_size files, accumulate, then drop references
    batch: list[str] = []
    for fpath in tqdm(all_files, desc="gsimple files"):
        batch.append(fpath)
        if len(batch) < args.batch_size:
            continue
        for fp in batch:
            try:
                process_one_file(fp)
            except Exception as e:
                tqdm.write(f"  [skip] {os.path.basename(fp)}: {e}")
                skipped.append({"path": fp, "error": str(e)})
        batch.clear()

    for fp in batch:
        try:
            process_one_file(fp)
        except Exception as e:
            tqdm.write(f"  [skip] {os.path.basename(fp)}: {e}")
            skipped.append({"path": fp, "error": str(e)})

    if n_files_ok == 0:
        raise SystemExit("No files were loaded successfully.")

    print(f"Files read OK: {n_files_ok} / {len(all_files)}")
    print(f"Total simulated POT (sum over files): {total_pot:.6g}")

    hist_zwf = hist_union_zwf[iz_default]
    n_rays_zwf = n_rays_union_zwf[iz_default]

    # --- finalize spectra (/ m^2 / POT / bin) ---
    def to_flux(hist: np.ndarray, area_m2: float) -> np.ndarray:
        return hist / (area_m2 * total_pot)

    area_sbnd = window_area_m2(SBND_VOLUME["x_range"], SBND_VOLUME["y_range"])
    area_flux_win = window_area_m2(FLUX_WINDOW_XY["x_range"], FLUX_WINDOW_XY["y_range"])
    area_face_cfg = window_area_m2((-200.0, 200.0), (-200.0, 200.0))

    sbnd_volume_spectra: dict[str, np.ndarray] = {}
    for iflav, flav in enumerate(FLAVORS):
        plane_specs = np.stack(
            [to_flux(hist_zwf[iz, IW_AV, iflav], area_sbnd) for iz in z_volume_indices],
            axis=0,
        )
        sbnd_volume_spectra[flav] = np.mean(plane_specs, axis=0)

    raytrace_volume_cm3 = np.array([volume_cm3(boxes) for _, boxes in RAYTRACE_VOLUME_DEFS])
    inv_vpot = 1.0e4 / (
        raytrace_volume_cm3[:, np.newaxis, np.newaxis] * total_pot
    )
    flux_density_raytrace_wvf = hist_raytrace_wvf * inv_vpot
    raytrace_volume_spectra: dict[str, dict[str, np.ndarray]] = {
        vname: {} for vname, _ in RAYTRACE_VOLUME_DEFS
    }
    for ivol, (vname, _boxes) in enumerate(RAYTRACE_VOLUME_DEFS):
        for iflav, flav in enumerate(FLAVORS):
            raytrace_volume_spectra[vname][flav] = flux_density_raytrace_wvf[ivol, iflav]

    z_slab_spectra = np.zeros((n_flavors, n_z_slabs, N_E_BINS))
    z_slab_integ = np.zeros((n_flavors, n_z_slabs))
    for iflav, flav in enumerate(FLAVORS):
        for iz in range(n_z_slabs):
            spec = to_flux(hist_zwf[iz, IW_AV, iflav], area_flux_win)
            z_slab_spectra[iflav, iz] = spec
            z_slab_integ[iflav, iz] = integrated_flux(spec)

    iz_face = z_slabs.index(FACE_Z_CM) if FACE_Z_CM in z_slabs else 0
    face_spectra_final: dict[str, dict[str, tuple[np.ndarray, int]]] = {}
    for iflav, flav in enumerate(FLAVORS):
        face_spectra_final[flav] = {}
        for lbl, win in FACE_WINDOWS.items():
            iw = FACE_LABEL_TO_WINDOW_IDX[lbl]
            a = window_area_m2(win["x_range"], win["y_range"])
            spec = to_flux(hist_zwf[iz_face, iw, iflav], a)
            face_spectra_final[flav][lbl] = (spec, int(n_rays_zwf[iz_face, iw, iflav]))

    histogram_npz = args.save_npz
    if histogram_npz is None and args.out_dir is not None:
        histogram_npz = args.out_dir / "flux_histograms.npz"

    window_area_m2_arr = np.array(
        [window_area_m2(t[1], t[2]) for t in FLUX_WINDOW_DEFS],
        dtype=np.float64,
    )

    if histogram_npz is not None:
        histogram_npz = histogram_npz.resolve()
        histogram_npz.parent.mkdir(parents=True, exist_ok=True)

        window_keys = np.array([t[0] for t in FLUX_WINDOW_DEFS], dtype=object)
        wx = np.array([t[1] for t in FLUX_WINDOW_DEFS], dtype=np.float64)
        wy = np.array([t[2] for t in FLUX_WINDOW_DEFS], dtype=np.float64)
        # flux density φ: (z, window, flavor, Ebin) in /m²/POT/bin
        inv_apot = 1.0 / (window_area_m2_arr[np.newaxis, :, np.newaxis, np.newaxis] * total_pot)
        flux_density_zwf = hist_zwf * inv_apot

        meta = {
            "gsimple_dir": args.gsimple_dir,
            "n_files_attempted": len(all_files),
            "n_files_ok": n_files_ok,
            "n_skipped": len(skipped),
            "batch_size": args.batch_size,
            "window_definitions": [
                {"key": k, "x_range_cm": list(xr), "y_range_cm": list(yr)}
                for k, xr, yr in FLUX_WINDOW_DEFS
            ],
            "flavors": FLAVORS,
            "pdg_by_flavor": PDG_BY_FLAVOR,
            "raytrace_volume_definitions": [
                {
                    "key": vname,
                    "label": RAYTRACE_VOLUME_LABEL[vname],
                    "volume_cm3": volume_cm3(boxes),
                    "boxes_cm": [
                        {
                            "x_range": list(b["x_range"]),
                            "y_range": list(b["y_range"]),
                            "z_range": list(b["z_range"]),
                        }
                        for b in boxes
                    ],
                }
                for vname, boxes in RAYTRACE_VOLUME_DEFS
            ],
        }
        (histogram_npz.parent / "flux_histograms_meta.json").write_text(
            json.dumps(meta, indent=2), encoding="utf-8"
        )
        print(f"Wrote {histogram_npz.parent / 'flux_histograms_meta.json'}")

        raytrace_volume_keys = np.array([t[0] for t in RAYTRACE_VOLUME_DEFS], dtype=object)

        np.savez_compressed(
            histogram_npz,
            total_pot=np.array(total_pot),
            E_BINS=E_BINS,
            BIN_CENTERS=BIN_CENTERS,
            z_positions_cm=np.array(z_slabs, dtype=np.float64),
            window_keys=window_keys,
            window_x_range_cm=wx,
            window_y_range_cm=wy,
            window_area_m2=window_area_m2_arr,
            flavors=np.array(FLAVORS, dtype=object),
            hist_weighted_zwf=hist_zwf,
            n_rays_zwf=n_rays_zwf,
            flux_density_zwf=flux_density_zwf,
            raytrace_volume_keys=raytrace_volume_keys,
            raytrace_volume_cm3=raytrace_volume_cm3,
            hist_weighted_raytrace_wvf=hist_raytrace_wvf,
            n_rays_raytrace_wvf=n_rays_raytrace_wvf,
            flux_density_raytrace_wvf=flux_density_raytrace_wvf,
        )
        print(f"Wrote {histogram_npz}")

        raytrace_npz = args.save_raytrace_npz
        if raytrace_npz is None:
            raytrace_npz = histogram_npz.parent / "raytrace_flux_histograms.npz"
        write_raytrace_npz(
            raytrace_npz,
            total_pot=total_pot,
            hist_wvf=hist_raytrace_wvf,
            n_rays_wvf=n_rays_raytrace_wvf,
            flux_density_wvf=flux_density_raytrace_wvf,
            volume_cm3_arr=raytrace_volume_cm3,
        )
        print(f"Wrote {raytrace_npz}")

        raytrace_meta = {
            "source": "gsimple_batch.py",
            "paired_flux_histograms_npz": str(histogram_npz),
            "n_files_ok": n_files_ok,
            "batch_size": args.batch_size,
            "volume_definitions": meta["raytrace_volume_definitions"],
        }
        raytrace_meta_path = raytrace_npz.parent / "raytrace_flux_histograms_meta.json"
        raytrace_meta_path.write_text(json.dumps(raytrace_meta, indent=2), encoding="utf-8")
        print(f"Wrote {raytrace_meta_path}")

        if not args.no_manifest:
            manifest_dir = histogram_npz.parent
            proc_txt = manifest_dir / "processed_files.txt"
            proc_txt.write_text("\n".join(ok_files) + ("\n" if ok_files else ""), encoding="utf-8")
            print(f"Wrote {proc_txt} ({len(ok_files)} paths)")
            skip_path = manifest_dir / "skipped_files.json"
            skip_path.write_text(json.dumps(skipped, indent=2), encoding="utf-8")
            print(f"Wrote {skip_path} ({len(skipped)} entries)")

    if write_slab_zconfig and slab_defs:
        if histogram_npz is not None:
            out_parent = histogram_npz.parent
        elif args.out_dir is not None:
            out_parent = args.out_dir.resolve()
            out_parent.mkdir(parents=True, exist_ok=True)
        else:
            raise SystemExit(
                "slab z-config output requires --out-dir or --save-npz / --save-slab-zconfig-npz"
            )
        hist_czwf = pack_hist_czwf_from_union(hist_union_zwf, slab_iz_maps)
        hist_compare_czvf = pack_hist_compare_czvf_from_union(
            hist_compare_union_czvf, slab_iz_maps
        )
        flux_czwf, flux_mean_cvwf, z_pos_cz = finalize_slab_zconfig_arrays(
            hist_czwf,
            hist_compare_czvf,
            slab_defs,
            total_pot=total_pot,
            window_area_m2_arr=window_area_m2_arr,
            compare_volume_defs=RAYTRACE_VOLUME_DEFS,
        )
        slab_zconfig_npz = args.save_slab_zconfig_npz
        if slab_zconfig_npz is None:
            slab_zconfig_npz = out_parent / "slab_zconfig_flux_histograms.npz"
        write_slab_zconfig_npz(
            slab_zconfig_npz,
            total_pot=total_pot,
            slab_defs=slab_defs,
            hist_czwf=hist_czwf,
            flux_density_czwf=flux_czwf,
            flux_density_mean_cvwf=flux_mean_cvwf,
            z_positions_cm_cz=z_pos_cz,
            window_area_m2_arr=window_area_m2_arr,
        )
        print(f"Wrote {slab_zconfig_npz}")

        slab_zconfig_meta = {
            "source": "gsimple_batch.py",
            "paired_flux_histograms_npz": str(histogram_npz) if histogram_npz else None,
            "n_files_ok": n_files_ok,
            "batch_size": args.batch_size,
            "slab_z_ns": args.slab_z_ns if args.slab_z_ns is not None else DEFAULT_SLAB_Z_NS,
            "slab_z_range_fv_cm": list(SLAB_Z_RANGE_FV_CM),
            "slab_config_definitions": [
                {
                    "key": d["key"],
                    "n_z": d["n_z"],
                    "dz_cm": d["dz_cm"],
                    "z_range_cm": list(d["z_range_cm"]),
                    "z_positions_cm": np.asarray(d["z_cm"]).tolist(),
                    "batch_default": bool(d.get("batch_default", False)),
                }
                for d in slab_defs
            ],
            "compare_volume_definitions": [
                {
                    "raytrace_volume_key": vkey,
                    "z_range_cm": list(SLAB_COMPARE_Z_RANGE_CM[vkey]),
                    "boxes_cm": [
                        {
                            "x_range": list(b["x_range"]),
                            "y_range": list(b["y_range"]),
                            "z_range": list(b["z_range"]),
                        }
                        for b in boxes
                    ],
                }
                for vkey, boxes in RAYTRACE_VOLUME_DEFS
            ],
        }
        slab_meta_path = slab_zconfig_npz.parent / "slab_zconfig_flux_histograms_meta.json"
        slab_meta_path.write_text(json.dumps(slab_zconfig_meta, indent=2), encoding="utf-8")
        print(f"Wrote {slab_meta_path}")

        if not args.no_manifest and histogram_npz is None:
            proc_txt = out_parent / "processed_files.txt"
            proc_txt.write_text("\n".join(ok_files) + ("\n" if ok_files else ""), encoding="utf-8")
            print(f"Wrote {proc_txt} ({len(ok_files)} paths)")

    if args.no_plots:
        return

    # --- plots (mirror notebook) ---
    def _save_or_show(fig: plt.Figure, name: str) -> None:
        if args.out_dir:
            args.out_dir.mkdir(parents=True, exist_ok=True)
            out = args.out_dir / f"{name}.png"
            fig.savefig(out, dpi=150)
            plt.close(fig)
            print(f"Saved {out}")
        else:
            plt.show()

    # SBND volume-averaged spectrum
    fig, ax = plt.subplots()
    for flav in FLAVORS:
        spectrum = sbnd_volume_spectra[flav]
        integ = integrated_flux(spectrum)
        ax.step(
            BIN_CENTERS,
            spectrum,
            where="mid",
            label=rf"{FLAVOR_LATEX[flav]} ($\int\phi\,dE$ = {integ:.3e})",
        )
    ax.set_xlim(0, 4)
    ax.set_xlabel("Neutrino energy [GeV]")
    ax.set_ylabel(r"$\phi$ [/m$^2$/POT/50MeV]")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.4)
    ax.set_yscale("log")
    fig.tight_layout()
    _save_or_show(fig, "sbnd_volume_mean_spectrum")

    for vname, _boxes in RAYTRACE_VOLUME_DEFS:
        fig, ax = plt.subplots()
        for flav in FLAVORS:
            spectrum = raytrace_volume_spectra[vname][flav]
            integ = integrated_flux(spectrum)
            ax.step(
                BIN_CENTERS,
                spectrum,
                where="mid",
                label=rf"{FLAVOR_LATEX[flav]} ($\int\phi\,dE$ = {integ:.3e})",
            )
        ax.set_xlim(0, 4)
        ax.set_xlabel("Neutrino energy [GeV]")
        ax.set_ylabel(r"$\phi$ [/m$^2$/POT/50MeV]")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=10)
        ax.set_title(f"Ray-traced volume mean — {vname}\n{RAYTRACE_VOLUME_LABEL[vname]}")
        fig.tight_layout()
        _save_or_show(fig, f"raytrace_volume_{vname}")

    for flav in FLAVORS:
        fig, ax = plt.subplots()
        for vname, _boxes in RAYTRACE_VOLUME_DEFS:
            spec = raytrace_volume_spectra[vname][flav]
            ivol = next(i for i, (vn, _) in enumerate(RAYTRACE_VOLUME_DEFS) if vn == vname)
            n_rays = int(n_rays_raytrace_wvf[ivol, FLAVORS.index(flav)])
            integ = integrated_flux(spec)
            ax.step(
                BIN_CENTERS,
                spec,
                where="mid",
                label=rf"{vname}, {n_rays:,} rays ($\int$={integ:.2e})",
            )
        ax.set_xlim(0, 4)
        ax.set_xlabel("Neutrino energy [GeV]")
        ax.set_ylabel(r"$\phi$ [/m$^2$/POT/50MeV]")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=9)
        ax.set_title(f"{FLAVOR_LATEX[flav]} — ray-traced FV volumes")
        fig.tight_layout()
        _save_or_show(fig, f"raytrace_volumes_compare_{flav}")

    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=z_slabs[0], vmax=z_slabs[-1])
    for iflav, flav in enumerate(FLAVORS):
        fig, ax = plt.subplots()
        for iz, z in enumerate(z_slabs):
            ax.step(
                BIN_CENTERS,
                z_slab_spectra[iflav, iz],
                where="mid",
                color=cmap(norm(z)),
                alpha=0.85,
                linewidth=1.2,
                label=f"z={z:.0f} cm",
            )
        ax.set_xlim(0, 4)
        ax.set_xlabel("Neutrino Energy [GeV]")
        ax.set_ylabel(r"$\phi$ [/m$^2$/POT]")
        ax.set_title(f"{FLAVOR_LATEX[flav]}")
        ax.grid(True, alpha=0.4)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label="z [cm]")
        fig.tight_layout()
        _save_or_show(fig, f"z_slabs_{flav}")

    fig, ax = plt.subplots()
    for iflav, flav in enumerate(FLAVORS):
        ax.plot(z_slabs, z_slab_integ[iflav], marker="o", markersize=4, label=FLAVOR_LATEX[flav])
    ax.set_xlabel("z slab center [cm]")
    ax.set_ylabel(r"$\int \phi\,dE$ [/m$^2$/POT]")
    ax.set_title("Integrated flux through the flux window vs. z")
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    _save_or_show(fig, "integrated_flux_vs_z")

    for flav in FLAVORS:
        fig, ax = plt.subplots()
        for label, (spec, n_rays) in face_spectra_final[flav].items():
            integ = integrated_flux(spec)
            ax.step(
                BIN_CENTERS,
                spec,
                where="mid",
                label=rf"{label}, {n_rays} rays ($\int$={integ:.2e})",
            )
        ax.set_xlim(0, 4)
        ax.set_xlabel("Neutrino energy [GeV]")
        ax.set_ylabel(r"$\phi$ [/m$^2$/POT]")
        ax.set_title(f"{FLAVOR_LATEX[flav]}")
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=9, loc="upper right")
        fig.tight_layout()
        _save_or_show(fig, f"face_windows_{flav}")

    cfg_plane_z_idx: dict[str, int] = {}
    for lbl, (kind, payload) in CONFIGS.items():
        if kind == "plane":
            zm = float(payload[0])
            cfg_plane_z_idx[lbl] = min(range(n_z_slabs), key=lambda i: abs(z_slabs[i] - zm))

    for iflav, flav in enumerate(FLAVORS):
        fig, ax = plt.subplots()
        for label, (kind, payload) in CONFIGS.items():
            if kind == "volume":
                plane_specs = np.stack(
                    [to_flux(hist_zwf[iz, IW_AV, iflav], area_sbnd) for iz in z_volume_indices],
                    axis=0,
                )
                spec = np.mean(plane_specs, axis=0)
                n_meta = n_vol_planes
            else:
                izp = cfg_plane_z_idx[label]
                spec = to_flux(hist_zwf[izp, IW_AV, iflav], area_face_cfg)
                n_meta = int(n_rays_zwf[izp, IW_AV, iflav])
            integ = integrated_flux(spec)
            ax.step(
                BIN_CENTERS,
                spec,
                where="mid",
                label=rf"{label}, {n_meta} planes/rays ($\int$={integ:.2e})",
            )
        ax.set_xlim(0, 4)
        ax.set_xlabel("Neutrino Energy [GeV]")
        ax.set_ylabel(r"$\phi$ [/m$^2$/POT]")
        ax.set_title(f"{FLAVOR_LATEX[flav]}")
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=9)
        fig.tight_layout()
        _save_or_show(fig, f"configs_compare_{flav}")


if __name__ == "__main__":
    main()
