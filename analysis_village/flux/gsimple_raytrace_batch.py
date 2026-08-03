#!/usr/bin/env python3
"""
Batch-process BNB gsimple ROOT flux files for ray-traced volume averages, like
``gsimple_raytrace.ipynb``, without loading all rays into memory.

The notebook concatenates every file's arrays, which can exhaust RAM on multi-thousand-file
sets. This script processes ``--batch-size`` files at a time, accumulates
``sum(wgt * path_length)`` per energy bin for each FV sub-volume, then divides by
``V * sum(POT)`` once at the end.

Physics matches ``gsimple_raytrace.ipynb``:

    phi(E) = 1e4 * sum(wgt * L_cm) / (V_cm3 * TOTAL_POT)

with ``TOTAL_POT = sum_k POT_k`` over processed files.

**Outputs** (when ``--out-dir`` and/or ``--save-npz`` is set):

* ``raytrace_flux_histograms.npz`` — ``hist_weighted_wvf``, ``flux_density_wvf`` with shape
  ``(n_volume, n_flavor, n_energy)``; ``n_rays_wvf``, bin edges, volume metadata.
* ``raytrace_flux_histograms_meta.json`` — run configuration and volume box definitions.
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

matplotlib.use(os.environ.get("GSIMPLE_MPL_BACKEND", os.environ.get("MPLBACKEND", "Agg")))
import matplotlib.pyplot as plt
import numpy as np
import uproot
from tqdm.auto import tqdm

# ---------------------------------------------------------------------------
# Defaults (same as gsimple_raytrace.ipynb)
# ---------------------------------------------------------------------------

DEFAULT_GSIMPLE_DIR = (
    "/cvmfs/sbnd.osgstorage.org/pnfs/fnal.gov/usr/sbnd/persistent/stash/"
    "fluxFiles/bnb/BooNEtoGSimple/configK-v1/july2023/neutrinoMode/"
)

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

from raytrace_volume_defs import (
    BIN_CENTERS,
    E_BINS,
    E_MAX_GEV,
    FV_SPLIT_BOXES,
    FV_SPLIT_TRUNCY_BOXES,
    N_E_BINS,
    RAYTRACE_VOLUME_DEFS,
    RAYTRACE_VOLUME_LABEL,
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
    """Load one file's ray arrays (cm for vertices). See gsimple_raytrace.ipynb."""
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
        v += (
            (b["x_range"][1] - b["x_range"][0])
            * (b["y_range"][1] - b["y_range"][0])
            * (b["z_range"][1] - b["z_range"][0])
        )
    return v


def add_raytrace_weighted_hist(
    data: dict,
    flavor: str,
    path_lengths: np.ndarray,
    out: np.ndarray,
) -> int:
    pdg = PDG_BY_FLAVOR[flavor]
    m = (data["pdg"] == pdg) & (path_lengths > 0)
    h, _ = np.histogram(data["E"][m], bins=E_BINS, weights=data["wgt"][m] * path_lengths[m])
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
        help="Write volume×flavor histograms to this .npz path. "
        "If omitted but --out-dir is set, defaults to OUT_DIR/raytrace_flux_histograms.npz",
    )
    p.add_argument(
        "--no-manifest",
        action="store_true",
        help="Do not write processed_files.txt / skipped_files.json next to the .npz",
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

    n_volumes = len(RAYTRACE_VOLUME_DEFS)
    n_flavors = len(FLAVORS)

    total_pot = 0.0
    n_files_ok = 0
    hist_wvf = np.zeros((n_volumes, n_flavors, N_E_BINS))
    n_rays_wvf = np.zeros((n_volumes, n_flavors), dtype=np.int64)

    ok_files: list[str] = []
    skipped: list[dict[str, str]] = []

    def process_one_file(fpath: str) -> None:
        nonlocal total_pot, n_files_ok
        d = load_gsimple(fpath)
        pot = float(d["pot"])

        p_mag = np.sqrt(d["px"] ** 2 + d["py"] ** 2 + d["pz"] ** 2)
        dx = d["px"] / p_mag
        dy = d["py"] / p_mag
        dz = d["pz"] / p_mag

        for ivol, (_vname, boxes) in enumerate(RAYTRACE_VOLUME_DEFS):
            path_L = path_length_volume(d["vtxx"], d["vtxy"], d["vtxz"], dx, dy, dz, boxes)
            for iflav, flav in enumerate(FLAVORS):
                n_rays_wvf[ivol, iflav] += add_raytrace_weighted_hist(
                    d, flav, path_L, hist_wvf[ivol, iflav]
                )

        total_pot += pot
        n_files_ok += 1
        ok_files.append(fpath)

    batch: list[str] = []
    for fpath in tqdm(all_files, desc="gsimple raytrace files"):
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

    volume_cm3_arr = np.array([volume_cm3(boxes) for _, boxes in RAYTRACE_VOLUME_DEFS])
    inv_vpot = 1.0e4 / (volume_cm3_arr[:, np.newaxis, np.newaxis] * total_pot)
    flux_density_wvf = hist_wvf * inv_vpot

    volume_spectra: dict[str, dict[str, np.ndarray]] = {
        vname: {} for vname, _ in RAYTRACE_VOLUME_DEFS
    }
    for ivol, (vname, _boxes) in enumerate(RAYTRACE_VOLUME_DEFS):
        for iflav, flav in enumerate(FLAVORS):
            volume_spectra[vname][flav] = flux_density_wvf[ivol, iflav]

    histogram_npz = args.save_npz
    if histogram_npz is None and args.out_dir is not None:
        histogram_npz = args.out_dir / "raytrace_flux_histograms.npz"

    if histogram_npz is not None:
        histogram_npz = histogram_npz.resolve()
        histogram_npz.parent.mkdir(parents=True, exist_ok=True)

        volume_keys = np.array([t[0] for t in RAYTRACE_VOLUME_DEFS], dtype=object)
        meta = {
            "gsimple_dir": args.gsimple_dir,
            "n_files_attempted": len(all_files),
            "n_files_ok": n_files_ok,
            "n_skipped": len(skipped),
            "batch_size": args.batch_size,
            "flavors": FLAVORS,
            "pdg_by_flavor": PDG_BY_FLAVOR,
            "volume_definitions": [
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
        meta_path = histogram_npz.parent / "raytrace_flux_histograms_meta.json"
        meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        print(f"Wrote {meta_path}")

        np.savez_compressed(
            histogram_npz,
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
        print(f"Wrote {histogram_npz}")

        if not args.no_manifest:
            manifest_dir = histogram_npz.parent
            proc_txt = manifest_dir / "processed_files.txt"
            proc_txt.write_text("\n".join(ok_files) + ("\n" if ok_files else ""), encoding="utf-8")
            print(f"Wrote {proc_txt} ({len(ok_files)} paths)")
            skip_path = manifest_dir / "skipped_files.json"
            skip_path.write_text(json.dumps(skipped, indent=2), encoding="utf-8")
            print(f"Wrote {skip_path} ({len(skipped)} entries)")

    def _save_or_show(fig: plt.Figure, name: str) -> None:
        if args.out_dir:
            args.out_dir.mkdir(parents=True, exist_ok=True)
            out = args.out_dir / f"{name}.png"
            fig.savefig(out, dpi=150)
            plt.close(fig)
            print(f"Saved {out}")
        else:
            plt.show()

    for vname, _boxes in RAYTRACE_VOLUME_DEFS:
        fig, ax = plt.subplots()
        for flav in FLAVORS:
            spectrum = volume_spectra[vname][flav]
            integ = integrated_flux(spectrum)
            ax.step(
                BIN_CENTERS,
                spectrum,
                where="mid",
                label=rf"{FLAVOR_LATEX[flav]} ($\int\phi\,dE$ = {integ:.3e})",
            )
        ax.set_xlim(0, E_MAX_GEV)
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
        for ivol, (vname, _boxes) in enumerate(RAYTRACE_VOLUME_DEFS):
            spec = volume_spectra[vname][flav]
            n_rays = int(n_rays_wvf[ivol, FLAVORS.index(flav)])
            integ = integrated_flux(spec)
            ax.step(
                BIN_CENTERS,
                spec,
                where="mid",
                label=rf"{vname}, {n_rays:,} rays ($\int$={integ:.2e})",
            )
        ax.set_xlim(0, E_MAX_GEV)
        ax.set_xlabel("Neutrino energy [GeV]")
        ax.set_ylabel(r"$\phi$ [/m$^2$/POT/50MeV]")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=9)
        ax.set_title(f"{FLAVOR_LATEX[flav]} — ray-traced FV volumes")
        fig.tight_layout()
        _save_or_show(fig, f"raytrace_volumes_compare_{flav}")


if __name__ == "__main__":
    main()
