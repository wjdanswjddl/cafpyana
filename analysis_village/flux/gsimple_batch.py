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
  metadata.
* ``flux_histograms_meta.json`` — run configuration and window list.
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

E_BINS = np.linspace(0.0, 4.0, 81)
BIN_CENTERS = 0.5 * (E_BINS[:-1] + E_BINS[1:])
N_E_BINS = len(BIN_CENTERS)

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


def integrated_flux(spectrum: np.ndarray) -> float:
    return float(np.sum(spectrum))


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
        "--no-manifest",
        action="store_true",
        help="Do not write processed_files.txt / skipped_files.json next to the histogram .npz",
    )
    p.add_argument("--explore-only", action="store_true", help="Print one file's structure and exit")
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
    n_vol_planes = len(z_volume_indices)
    IW_AV = 0  # AV ±200 cm; same geometry as ``FLUX_WINDOW_XY`` / SBND AV plots

    total_pot = 0.0
    n_files_ok = 0

    n_windows = len(FLUX_WINDOW_DEFS)
    n_flavors = len(FLAVORS)

    # z × window × flavor × energy: raw sum of weights in each E bin (same as np.histogram weights).
    hist_zwf = np.zeros((n_z_slabs, n_windows, n_flavors, N_E_BINS))
    n_rays_zwf = np.zeros((n_z_slabs, n_windows, n_flavors), dtype=np.int64)

    ok_files: list[str] = []
    skipped: list[dict[str, str]] = []

    def process_one_file(fpath: str) -> None:
        nonlocal total_pot, n_files_ok
        d = load_gsimple(fpath)
        pot = float(d["pot"])

        for iz, z in enumerate(z_slabs):
            for iw, (_name, xr, yr) in enumerate(FLUX_WINDOW_DEFS):
                for iflav, flav in enumerate(FLAVORS):
                    n_rays_zwf[iz, iw, iflav] += add_weighted_e_hist(
                        d, flav, z, xr, yr, hist_zwf[iz, iw, iflav]
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

    if histogram_npz is not None:
        histogram_npz = histogram_npz.resolve()
        histogram_npz.parent.mkdir(parents=True, exist_ok=True)

        window_keys = np.array([t[0] for t in FLUX_WINDOW_DEFS], dtype=object)
        wx = np.array([t[1] for t in FLUX_WINDOW_DEFS], dtype=np.float64)
        wy = np.array([t[2] for t in FLUX_WINDOW_DEFS], dtype=np.float64)
        window_area_m2_arr = np.array(
            [window_area_m2(t[1], t[2]) for t in FLUX_WINDOW_DEFS],
            dtype=np.float64,
        )
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
        }
        (histogram_npz.parent / "flux_histograms_meta.json").write_text(
            json.dumps(meta, indent=2), encoding="utf-8"
        )
        print(f"Wrote {histogram_npz.parent / 'flux_histograms_meta.json'}")

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
        )
        print(f"Wrote {histogram_npz}")

        if not args.no_manifest:
            manifest_dir = histogram_npz.parent
            proc_txt = manifest_dir / "processed_files.txt"
            proc_txt.write_text("\n".join(ok_files) + ("\n" if ok_files else ""), encoding="utf-8")
            print(f"Wrote {proc_txt} ({len(ok_files)} paths)")
            skip_path = manifest_dir / "skipped_files.json"
            skip_path.write_text(json.dumps(skipped, indent=2), encoding="utf-8")
            print(f"Wrote {skip_path} ({len(skipped)} entries)")

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
