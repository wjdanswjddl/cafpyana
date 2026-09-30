#!/usr/bin/env python3
"""Plot unfolded data vs generator predictions (eye-check).

Matches ``generator_comparison.ipynb`` / ``generator_predictions.ipynb``:
topology signal mask + ``40 × fScaleFactor × Weight``, GiBUU ``/ n_runs``.

That phase space is the **extracted** cross section
(``signal_truth_fv='none'`` on ``nevts_allmc``): selection / per_tpc end FV
are unfolded out via the response, not part of the generator prediction.

Requires a response+unfold rebuilt with ``signal_truth_fv='none'`` for
``nevts_allmc`` (see ``response_matrices_product_b.py``).

  unset PYTHONPATH
  envs/venv_py310_cafpyana/bin/python \\
    analysis_village/numucc_1p0pi/scripts/plot_generator_predictions_eye.py
"""

from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
_SCRIPTS = Path(__file__).resolve().parent
for p in (_REPO, _SCRIPTS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
)
from analysis_village.numucc_1p0pi.utils import dpi, fig_ext, plot_unfolded_result  # noqa: E402
from run_detector_diag_test_unfold_compare import (  # noqa: E402
    GENIE_AR23_LIANG_BUGFIX_FLAT,
    GENIE_AR23_LIANG_BUGFIX_LABEL,
    GIBUU_FLAT,
    NEUT_FLAT,
    NOMINAL_UNFOLD,
    TEST_UNFOLD,
    _gibuu_n_runs,
    _load_flat_spectra,
    _resolve_genie_ar23_nominal_flat,
)

OUT_DIR = TEST_UNFOLD / "plots_generator_predictions"
VARS = [
    "muon-p",
    "muon-dir_z",
    "proton-p",
    "proton-dir_z",
    "tki-del_Tp",
    "tki-del_alpha",
    "tki-del_phi",
]


def _unfold_dict(res: dict) -> dict:
    return {
        "unfold": np.asarray(res["unfold"], float),
        "AddSmear": np.asarray(res["AddSmear"], float),
        "UnfoldCov": np.asarray(res["UnfoldCov"], float),
        "StatUnfoldCov": np.asarray(res["StatUnfoldCov"], float),
        "SystUnfoldCov": np.asarray(res["SystUnfoldCov"], float),
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    vc_by = {vc.var_save_name: vc for vc in CORE_SELECTED_EVT_VARIABLE_CONFIGS}

    with open(NOMINAL_UNFOLD / "unfolding_ingredients_and_results.pkl", "rb") as f:
        nom = pickle.load(f)
    results = nom["results"]

    n_runs = _gibuu_n_runs(GIBUU_FLAT)
    ar23_flat_path = _resolve_genie_ar23_nominal_flat()
    print(
        f"Loading flats (40×fSF×W topology; GiBUU /{n_runs})...\n"
        f"  nominal AR23: {ar23_flat_path}\n"
        f"  {GENIE_AR23_LIANG_BUGFIX_LABEL}: {GENIE_AR23_LIANG_BUGFIX_FLAT}",
        flush=True,
    )
    genie_flat = _load_flat_spectra(ar23_flat_path, 1.0)
    liang_flat = _load_flat_spectra(GENIE_AR23_LIANG_BUGFIX_FLAT, 1.0)
    neut_flat = _load_flat_spectra(NEUT_FLAT, 1.0)
    gibuu_flat = _load_flat_spectra(GIBUU_FLAT, 1.0 / float(n_runs))

    pb_model = float(np.sum(results["muon-p"]["model"]))
    print(
        "Integrated σ muon-p: "
        f"AR23_flat={np.sum(genie_flat['muon-p']):.4e}  "
        f"Liang={np.sum(liang_flat['muon-p']):.4e}  "
        f"NEUT={np.sum(neut_flat['muon-p']):.4e}  "
        f"GiBUU={np.sum(gibuu_flat['muon-p']):.4e}  "
        f"unfold_model={pb_model:.4e}",
        flush=True,
    )

    manifest = {
        "out_dir": str(OUT_DIR),
        "nominal_unfold": str(NOMINAL_UNFOLD),
        "paths": {
            "genie_ar23_nominal": ar23_flat_path,
            "genie_liang_bugfix": GENIE_AR23_LIANG_BUGFIX_FLAT,
            "neut": NEUT_FLAT,
            "gibuu": GIBUU_FLAT,
        },
        "gibuu_n_runs": n_runs,
        "note": (
            "Weights = 40×fSF×W, topology mask (generator_comparison.ipynb). "
            "GiBUU / n_runs. Phase space = unfold signal_truth_fv='none'. "
            "Nominal AR23 flat = pre–Liang-fix Nieves QE splines; "
            f"{GENIE_AR23_LIANG_BUGFIX_LABEL!r} = munjung SBND_gen1. "
            "plot_unfolded_result applies A_C."
        ),
        "plots": [],
    }

    for vsn in VARS:
        if vsn not in results or vsn not in vc_by:
            print(f"skip {vsn}", flush=True)
            continue
        res = results[vsn]
        vc = vc_by[vsn]
        models = {
            "GENIE AR23 (prod model)": [np.asarray(res["model"], float), "C0"],
            "GENIE AR23 flat": [np.asarray(genie_flat[vsn], float), "C3"],
            GENIE_AR23_LIANG_BUGFIX_LABEL: [np.asarray(liang_flat[vsn], float), "C4"],
            "NEUT 6.1.4": [np.asarray(neut_flat[vsn], float), "C2"],
            f"GiBUU 2025 /{n_runs}": [np.asarray(gibuu_flat[vsn], float), "C1"],
        }
        save_name = str(OUT_DIR / f"{vsn}-xsec_generators")
        plot_unfolded_result(
            _unfold_dict(res),
            np.asarray(res["measured"], float),
            models,
            vc,
            xsec_unit=float(res.get("xsec_unit", nom["meta"]["xsec_unit"])),
            save_fig=True,
            save_name=save_name,
            textloc=[0.55, 0.92],
            approval="",
            data=True,
            closure_test=False,
        )
        print("wrote", save_name + fig_ext, flush=True)
        manifest["plots"].append(save_name + fig_ext)

    man_path = OUT_DIR / "manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2))
    print("wrote", man_path)
    print("OUT_DIR =", OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
