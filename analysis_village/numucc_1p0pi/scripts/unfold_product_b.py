#!/usr/bin/env python3
"""PRL Product B Wiener-SVD unfold (script twin of ``notebooks/unfolding.ipynb``).

Loads:

* ``PRL/response_matrices/response_matrices.npz`` (eff, response, MC rates, data−bkg)
* ``PRL/systematics/productB_sel_mup`` CategorySummary ``total_xsec``
  (nominal since 2026-09-29: GENIE = ``GENIE_slim_v3`` = base × FSI_v3 with the MEC
  knobs swapped for their May versions; DENT = rolling 80% w=3 + Gaussian σ=1 on σ/N;
  Detector = WireMod nested + that DENT + dE/dx smear26). Raw (unsmoothed) DENT:
  ``productB_sel_mup__dent_raw`` → ``productB_sel_mup__genie_FSIv3_MEC_May`` /
  ``PRL/unfolded_dent_raw``. Former GENIE (``GENIE_slim_both``, FSI v1×v3):
  ``productB_sel_mup__FSI_v1v3`` / ``PRL/unfolded_FSI_v1v3``.
* Gen1 ray-trace flux → ``xsec_unit``
* MC rates scaled by ``MC_POT_FIX`` by default (recorded MC POT is high)

Runs MC Asimov closure + data unfold (``C_type=2``, ``Norm_type=0``), writes
results/plots under ``PRL/unfolded/``. Per variable: syst-source breakdown plus
response / fractional-covariance / ``A_c`` heatmaps.

Legacy (GENIE_slim_v3 only): ``productB_sel_mup_legacy_slim_v3`` /
``unfolded_legacy_slim_v3``.
"""

from __future__ import annotations

import argparse
import json
import pickle
import shutil
import sys
import warnings
from datetime import datetime, timezone
from os import makedirs
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np

_REPO = Path(__file__).resolve().parents[3]
_SCRIPTS = Path(__file__).resolve().parent
for p in (_REPO, _SCRIPTS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

warnings.filterwarnings("ignore", category=FutureWarning)

from analysis_village.numucc_1p0pi.constants import M_AR, MC_POT_FIX, N_A, RHO  # noqa: E402
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
)
from analysis_village.numucc_1p0pi.syst_category_summary import (  # noqa: E402
    load_category_syst_summary,
    total_cov_frac,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import (  # noqa: E402
    category_summary_npz_path,
)
from analysis_village.numucc_1p0pi.utils import (  # noqa: E402
    cov_from_fraccov,
    get_chi2,
    get_integrated_flux,
    plot_unfolded_result,
    save_unfold_diagnostics,
)
from analysis_village.unfolding.wienersvd import WienerSVD  # noqa: E402
from analysis_village.flux.raytrace_volume_defs import (  # noqa: E402
    FV_SPLIT_TRUNCY_BOXES,
    RAYTRACE_VOLUME_LABEL,
)
PRL_ROOT = Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL")
DEFAULT_RESP = PRL_ROOT / "response_matrices" / "response_matrices.npz"
DEFAULT_SYST = PRL_ROOT / "systematics" / "productB_sel_mup"
DEFAULT_OUT = PRL_ROOT / "unfolded"
FLUX_FILE = Path("/exp/sbnd/data/users/munjung/flux/SBND_gsimple_raytrace/Gen1.root")

C_TYPE = 2
NORM_TYPE = 0.0


def pack_unfold_results(unfold: dict, var_config) -> dict:
    """Bin-wise unfold, covariance diagonals, and per-bin-width errors."""
    bins = var_config.bins
    bin_widths = np.diff(bins)
    if len(bins) == 2:
        bin_widths = np.array([1.0])
    unfolded = np.asarray(unfold["unfold"], dtype=float)
    return {
        "bins": np.asarray(bins, dtype=float),
        "bin_centers": np.asarray(var_config.bin_centers, dtype=float),
        "bin_widths": bin_widths,
        "unfold": unfolded,
        "unfold_per_bin_width": unfolded / bin_widths,
        "stat_err": np.sqrt(np.maximum(np.diag(unfold["StatUnfoldCov"]), 0.0)),
        "syst_err": np.sqrt(np.maximum(np.diag(unfold["SystUnfoldCov"]), 0.0)),
        "total_err": np.sqrt(np.maximum(np.diag(unfold["UnfoldCov"]), 0.0)),
        "stat_err_per_bin_width": np.sqrt(np.maximum(np.diag(unfold["StatUnfoldCov"]), 0.0)) / bin_widths,
        "syst_err_per_bin_width": np.sqrt(np.maximum(np.diag(unfold["SystUnfoldCov"]), 0.0)) / bin_widths,
        "total_err_per_bin_width": np.sqrt(np.maximum(np.diag(unfold["UnfoldCov"]), 0.0)) / bin_widths,
        "AddSmear": np.asarray(unfold["AddSmear"], dtype=float),
        "UnfoldCov": np.asarray(unfold["UnfoldCov"], dtype=float),
        "StatUnfoldCov": np.asarray(unfold["StatUnfoldCov"], dtype=float),
        "SystUnfoldCov": np.asarray(unfold["SystUnfoldCov"], dtype=float),
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--response-npz", default=str(DEFAULT_RESP))
    p.add_argument("--syst-root", default=str(DEFAULT_SYST))
    p.add_argument("--out-dir", default=str(DEFAULT_OUT))
    p.add_argument("--flux-file", default=str(FLUX_FILE))
    p.add_argument("--skip-closure", action="store_true")
    p.add_argument(
        "--mc-pot-fix",
        dest="mc_pot_fix",
        action="store_true",
        default=True,
        help="Recorded neutrino-MC POT is high: multiply nevts_allmc and "
        "nevts_sel_reco by MC_POT_FIX (same as overlay --mc-pot-fix). Data and "
        "response matrix unchanged. Default on.",
    )
    p.add_argument(
        "--no-mc-pot-fix",
        dest="mc_pot_fix",
        action="store_false",
        help="Disable the MC_POT_FIX scale on neutrino-MC rates.",
    )
    return p.parse_args(argv)


def asimov_covariance(nevts_sel_reco: np.ndarray, xsec_unit: float) -> np.ndarray:
    n = np.asarray(nevts_sel_reco, dtype=float)
    frac = np.where(n > 0, 1.0 / n, 0.0)
    return cov_from_fraccov(np.diag(frac), n) * (xsec_unit ** 2)


def _fix_integrated_n_sel_data(
    variables: Dict[str, Dict[str, Any]],
) -> None:
    """``integrated`` may lack overlay counts; use total from a CORE differential var."""
    integ = variables.get("integrated")
    if integ is None:
        return
    n_sel = np.asarray(integ.get("n_sel_data", []), dtype=float)
    if n_sel.size and float(np.sum(n_sel)) > 0:
        return
    for ref in ("muon-p", "tki-del_Tp", "proton-p"):
        if ref not in variables:
            continue
        tot = float(np.sum(np.asarray(variables[ref]["n_sel_data"], dtype=float)))
        if tot <= 0:
            continue
        integ["n_sel_data"] = np.array([tot], dtype=float)
        if "n_data" in variables[ref]:
            integ["n_data"] = np.array(
                [float(np.sum(variables[ref]["n_data"]))], dtype=float
            )
        if "n_mc_bkg" in variables[ref]:
            integ["n_mc_bkg"] = np.array(
                [float(np.sum(variables[ref]["n_mc_bkg"]))], dtype=float
            )
        print(f"  filled integrated n_sel_data={tot:.1f} from {ref}", flush=True)
        return
    print("  WARN: could not fill integrated n_sel_data", flush=True)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    resp_npz = Path(args.response_npz)
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "plots"
    makedirs(fig_dir, exist_ok=True)

    cat_npz = Path(category_summary_npz_path(str(args.syst_root)))
    flux_file = Path(args.flux_file)

    if not resp_npz.is_file():
        raise FileNotFoundError(resp_npz)
    if not cat_npz.is_file():
        raise FileNotFoundError(cat_npz)
    if not flux_file.is_file():
        raise FileNotFoundError(flux_file)

    print("response :", resp_npz)
    print("category :", cat_npz)
    print("flux     :", flux_file)
    print("out      :", out_dir)

    blob = np.load(resp_npz, allow_pickle=True)
    meta = json.loads(str(blob["meta_json"][0]))
    data_tot_pot = float(meta["data_tot_pot"])

    variables: Dict[str, Dict[str, Any]] = {}
    vsns = sorted({k.split("::", 1)[0] for k in blob.files if "::" in k})
    for vsn in vsns:
        pack: Dict[str, Any] = {}
        for key in blob.files:
            if key.startswith(vsn + "::"):
                pack[key.split("::", 1)[1]] = blob[key]
        variables[vsn] = pack
    _fix_integrated_n_sel_data(variables)

    mc_pot_fix_factor = MC_POT_FIX if args.mc_pot_fix else None
    if mc_pot_fix_factor is not None:
        # Overlay convention: recorded MC POT 3.2% high → MC rates at data POT
        # were 3.2% low. Scale truth and selected-signal MC; R and data stay put.
        for pack in variables.values():
            for key in ("nevts_allmc", "nevts_sel_reco"):
                if key in pack:
                    pack[key] = np.asarray(pack[key], dtype=float) * mc_pot_fix_factor
        print(
            f"MC POT fix ×{mc_pot_fix_factor}: scaled nevts_allmc and "
            "nevts_sel_reco (data / response unchanged)",
            flush=True,
        )

    summary = load_category_syst_summary(str(cat_npz))
    vc_by = {vc.var_save_name: vc for vc in CORE_SELECTED_EVT_VARIABLE_CONFIGS}

    # --- xsec unit ---
    print("Fiducial volume (FV_split_truncY):", RAYTRACE_VOLUME_LABEL["FV_split_truncY"])
    v_sbnd = 0.0
    for i, box in enumerate(FV_SPLIT_TRUNCY_BOXES):
        dx = box["x_range"][1] - box["x_range"][0]
        dy = box["y_range"][1] - box["y_range"][0]
        dz = box["z_range"][1] - box["z_range"][0]
        v_box = dx * dy * dz
        v_sbnd += v_box
        print(f"  slab {i+1}: -> {v_box:.4e} cm3")
    print(f"V_SBND = {v_sbnd:.6e} cm3")

    integrated_flux_per_pot = get_integrated_flux(str(flux_file), plot=False)
    integrated_flux = integrated_flux_per_pot * data_tot_pot
    n_targets = (RHO * v_sbnd / M_AR) * N_A
    xsec_unit = 1.0 / (integrated_flux * n_targets)
    flux_info = {
        "flux_file": str(flux_file),
        "flux_fv": "FV_split_truncY",
        "flux_fv_label": RAYTRACE_VOLUME_LABEL["FV_split_truncY"],
        "volume_cm3": float(v_sbnd),
        "integrated_flux_per_pot": float(integrated_flux_per_pot),
        "data_tot_pot": float(data_tot_pot),
        "integrated_flux": float(integrated_flux),
        "n_targets": float(n_targets),
        "xsec_unit": float(xsec_unit),
    }
    print(f"xsec_unit = {xsec_unit:.6e} cm2/nucleon")

    # --- MC Asimov closure ---
    closure_results: Dict[str, Any] = {}
    if not args.skip_closure:
        for vsn in vsns:
            vc = vc_by.get(vsn)
            if vc is None:
                continue
            pack = variables[vsn]
            response = np.asarray(pack["response"], dtype=float)
            nevts_allmc = np.asarray(pack["nevts_allmc"], dtype=float)
            nevts_sel_reco = np.asarray(pack["nevts_sel_reco"], dtype=float)
            model = np.atleast_1d(np.asarray(nevts_allmc, dtype=float) * xsec_unit)
            measured = np.atleast_1d(np.asarray(nevts_sel_reco, dtype=float) * xsec_unit)
            cov = asimov_covariance(nevts_sel_reco, xsec_unit)
            unfold = WienerSVD(
                response, model, measured, cov, C_TYPE, NORM_TYPE, stat_scaling=xsec_unit
            )
            ac_model = np.atleast_1d(
                np.asarray(unfold["AddSmear"], dtype=float) @ model
            )
            u = np.atleast_1d(np.asarray(unfold["unfold"], dtype=float))
            ucov = np.asarray(unfold["UnfoldCov"], dtype=float)
            mask = (u > 0) & (ac_model > 0)
            ndof = int(np.count_nonzero(mask))
            if ndof >= 1:
                chi2, p_val = get_chi2(u[mask], ac_model[mask], ucov[np.ix_(mask, mask)])
            else:
                chi2, p_val = float("nan"), float("nan")
            print(f"[closure] {vsn}: chi2/ndof = {chi2:.3f}/{ndof}  p={p_val:.3g}")
            # plot_unfolded_result expects models[name] = [array, color]
            plot_unfolded_result(
                unfold,
                measured,
                {
                    "GENIE (smeared truth)": [ac_model, "C0"],
                    "GENIE (truth)": [model, "C1"],
                },
                vc,
                plot=False,
                save_fig=True,
                save_name=str(fig_dir / f"{vsn}-closure"),
                closure_test=True,
            )
            closure_results[vsn] = {
                "chi2": float(chi2),
                "ndof": int(ndof),
                "p_val": float(p_val),
            }
        print("closure done", flush=True)

    # --- data unfold ---
    data_results: Dict[str, Any] = {}
    ingredients: Dict[str, Any] = {}
    for vsn in vsns:
        vc = vc_by.get(vsn)
        if vc is None:
            continue
        pack = variables[vsn]
        response = np.asarray(pack["response"], dtype=float)
        nevts_allmc = np.asarray(pack["nevts_allmc"], dtype=float)
        nevts_sel_reco = np.asarray(pack["nevts_sel_reco"], dtype=float)
        n_sel_data = np.asarray(pack["n_sel_data"], dtype=float)
        if float(np.sum(n_sel_data)) <= 0:
            print(f"[data] skip {vsn}: n_sel_data empty")
            continue

        model = np.atleast_1d(np.asarray(nevts_allmc, dtype=float) * xsec_unit)
        measured = np.atleast_1d(np.asarray(n_sel_data, dtype=float) * xsec_unit)
        try:
            frac_cov = total_cov_frac(summary, vsn, kind="xsec")
        except KeyError as ex:
            print(f"[data] skip {vsn}: no total_xsec ({ex})")
            continue

        Covariance = cov_from_fraccov(frac_cov, nevts_sel_reco) * (xsec_unit ** 2)
        unfold = WienerSVD(
            response,
            model,
            measured,
            Covariance,
            C_TYPE,
            NORM_TYPE,
            stat_scaling=xsec_unit,
        )
        ac_model = np.atleast_1d(np.asarray(unfold["AddSmear"], dtype=float) @ model)
        u = np.atleast_1d(np.asarray(unfold["unfold"], dtype=float))
        ucov = np.asarray(unfold["UnfoldCov"], dtype=float)
        mask = (u > 0) & (ac_model > 0)
        ndof = int(np.count_nonzero(mask))
        if ndof >= 1:
            chi2, p_val = get_chi2(u[mask], ac_model[mask], ucov[np.ix_(mask, mask)])
        else:
            chi2, p_val = float("nan"), float("nan")

        print(f"[data] {vsn}: chi2/ndof = {chi2:.3f}/{ndof}  p={p_val:.3g}")
        # plot_unfolded_result expects models[name] = [array, color]
        # Pass unsmeared truth; plot_unfolded_result applies A_c internally.
        plot_unfolded_result(
            unfold,
            measured,
            {"GENIE AR23": [model, "C0"]},
            vc,
            plot=False,
            save_fig=True,
            save_name=str(fig_dir / f"{vsn}-data_unfold"),
            data=True,
        )
        save_unfold_diagnostics(
            vc,
            summary,
            response,
            frac_cov,
            np.asarray(unfold["AddSmear"], dtype=float),
            fig_dir,
            plot=False,
        )

        packed = pack_unfold_results(unfold, vc)
        packed["chi2_vs_genie"] = float(chi2)
        packed["ndof_vs_genie"] = int(ndof)
        packed["measured"] = measured
        packed["model"] = model
        packed["response"] = response
        packed["Covariance_input"] = Covariance
        packed["frac_cov_total_xsec"] = np.asarray(frac_cov, dtype=float)
        packed["n_sel_data"] = n_sel_data
        packed["nevts_sel_reco"] = nevts_sel_reco
        packed["nevts_allmc"] = nevts_allmc
        packed["xsec_unit"] = float(xsec_unit)
        data_results[vsn] = packed
        ingredients[vsn] = {
            "bins": np.asarray(pack["bins"], dtype=float),
            "response": response,
            "eff": np.asarray(pack["eff"], dtype=float),
            "reco_vs_true": np.asarray(pack["reco_vs_true"], dtype=float),
            "n_data": np.asarray(pack.get("n_data", np.zeros_like(n_sel_data)), dtype=float),
            "n_mc_bkg": np.asarray(pack.get("n_mc_bkg", np.zeros_like(n_sel_data)), dtype=float),
            "n_sel_data": n_sel_data,
            "nevts_allmc": nevts_allmc,
            "nevts_sel_reco": nevts_sel_reco,
            "measured": measured,
            "model": model,
            "Covariance_input": Covariance,
            "frac_cov_total_xsec": np.asarray(frac_cov, dtype=float),
        }

    print("data unfold done for", list(data_results), flush=True)

    # --- save ---
    flux_dest = out_dir / "Gen1_flux.root"
    if not flux_dest.exists() or flux_dest.stat().st_size != flux_file.stat().st_size:
        shutil.copy2(flux_file, flux_dest)
        print("copied flux ->", flux_dest)
    else:
        print("flux already present:", flux_dest)

    flux_json = out_dir / "flux_info.json"
    with open(flux_json, "w") as f:
        json.dump(flux_info, f, indent=2)

    release = {
        "meta": {
            "schema": "prl_productB_unfolded_v1",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "response_npz": str(resp_npz),
            "category_summary_npz": str(cat_npz),
            "c_type": C_TYPE,
            "norm_type": NORM_TYPE,
            "xsec_unit": float(xsec_unit),
            "data_tot_pot": float(data_tot_pot),
            "mc_pot_fix_factor": mc_pot_fix_factor,
            "flux_info": flux_info,
            "response_meta": meta,
            "closure": {
                k: {"chi2": v["chi2"], "ndof": v["ndof"]}
                for k, v in closure_results.items()
            },
            "TODO_genie_flat_version_correction_1p03": (
                "Keep investigating ~3% GENIE-flat vs production residual "
                "(notebook VERSION_CORRECTION=1.03)."
            ),
        },
        "ingredients": ingredients,
        "results": data_results,
    }

    pkl_path = out_dir / "unfolding_ingredients_and_results.pkl"
    with open(pkl_path, "wb") as f:
        pickle.dump(release, f, protocol=pickle.HIGHEST_PROTOCOL)
    print("wrote", pkl_path)

    savez: Dict[str, Any] = {
        "meta_json": np.array([json.dumps(release["meta"], default=str)], dtype=object),
        "xsec_unit": np.float64(xsec_unit),
        "data_tot_pot": np.float64(data_tot_pot),
    }
    for vsn, res in data_results.items():
        for key in (
            "bins",
            "bin_centers",
            "bin_widths",
            "unfold",
            "unfold_per_bin_width",
            "stat_err",
            "syst_err",
            "total_err",
            "stat_err_per_bin_width",
            "syst_err_per_bin_width",
            "total_err_per_bin_width",
            "AddSmear",
            "UnfoldCov",
            "StatUnfoldCov",
            "SystUnfoldCov",
            "measured",
            "model",
            "response",
            "Covariance_input",
            "frac_cov_total_xsec",
            "n_sel_data",
            "nevts_sel_reco",
            "nevts_allmc",
        ):
            if key in res:
                savez[f"{vsn}::{key}"] = np.asarray(res[key])
        savez[f"{vsn}::chi2_vs_genie"] = np.float64(res.get("chi2_vs_genie", np.nan))
        savez[f"{vsn}::ndof_vs_genie"] = np.int32(res.get("ndof_vs_genie", -1))

    npz_path = out_dir / "unfolding_ingredients_and_results.npz"
    np.savez_compressed(npz_path, **savez)
    print("wrote", npz_path)

    per_var_dir = out_dir / "by_variable"
    makedirs(per_var_dir, exist_ok=True)
    for vsn, res in data_results.items():
        p = per_var_dir / f"{vsn}.npz"
        np.savez_compressed(
            p,
            bins=res["bins"],
            bin_centers=res["bin_centers"],
            unfold=res["unfold"],
            UnfoldCov=res["UnfoldCov"],
            StatUnfoldCov=res["StatUnfoldCov"],
            SystUnfoldCov=res["SystUnfoldCov"],
            AddSmear=res["AddSmear"],
            measured=res["measured"],
            model=res["model"],
            response=res["response"],
            xsec_unit=np.float64(xsec_unit),
            chi2_vs_genie=np.float64(res.get("chi2_vs_genie", np.nan)),
            ndof_vs_genie=np.int32(res.get("ndof_vs_genie", -1)),
        )
        print("wrote", p)

    manifest = {
        "schema": "prl_productB_unfolded_v1",
        "created_utc": release["meta"]["created_utc"],
        "pkl": str(pkl_path),
        "npz": str(npz_path),
        "flux_info": str(flux_json),
        "flux_file_copy": str(flux_dest),
        "plots_dir": str(fig_dir),
        "variables": sorted(data_results),
        "closure_chi2": {
            k: {"chi2": v["chi2"], "ndof": v["ndof"]} for k, v in closure_results.items()
        },
        "data_chi2_vs_genie": {
            k: {"chi2": v.get("chi2_vs_genie"), "ndof": v.get("ndof_vs_genie")}
            for k, v in data_results.items()
        },
    }
    man_path = out_dir / "unfolded_manifest.json"
    with open(man_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print("wrote", man_path)
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
