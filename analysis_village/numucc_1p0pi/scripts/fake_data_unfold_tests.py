#!/usr/bin/env python3
"""Stream Product B MC + GiBUU dfs for unfolding fake-data tests.

Reweight tests match the archived notebook (``FakeDataWeights`` + the same
``FAKE_DATA_TEST_SPECS``). Product B updates:

* Response / syst / ``xsec_unit`` come from the live unfold
  (``unfolding.ipynb``): Gen1 flux × ``FV_split_truncY`` × data POT,
  ``MC_POT_FIX`` from ``constants``, ``signal_truth_fv='none'`` for the truth model.
* Injected-truth weights are applied on the unfold signal
  (``IsNuInFV_NumuCC_1p0pi(..., 'none')``), not ``topo_categ==1`` / ``per_tpc``.
* GiBUU reco is ``sel_mup`` selected events; GiBUU truth is ``sel_mup`` ``mcnu``
  (this campaign's ``sel_all`` dfs have no ``mcnu``). Weights are
  ``genweight``; no flat ×1000 and no GENIE ``MC_POT_FIX``.
* ``mec_test`` / ``qe_test`` are GENIE-only: CAF ``mc.genie_mode`` is 0/10 for
  GENIE QE/MEC and 1 / 35|36 for GiBUU. Do not apply those reweights to GiBUU.
"""
from __future__ import annotations

import argparse
import gc
import glob
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

_REPO = Path(__file__).resolve().parents[3]
_SCRIPTS = Path(__file__).resolve().parent
for p in (_REPO, _SCRIPTS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

from pyanalib.split_df_helpers_new import get_n_split, load_dfs  # noqa: E402

from analysis_village.numucc_1p0pi.categories import IsNuInFV_NumuCC_1p0pi  # noqa: E402
from analysis_village.numucc_1p0pi.fake_data_test_configs import (  # noqa: E402
    FakeDataWeights,
    format_bump_test_label,
)
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
)
from analysis_village.numucc_1p0pi.utils import get_clipped_evts  # noqa: E402
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402
import data_mc_overlay_products as dmo  # noqa: E402
from response_matrices_product_b import _prepare_evt_mcnu  # noqa: E402

SIGNAL_TRUTH_FV_UNFOLD = "none"
TABLE_VAR_NAMES = (
    "muon-p",
    "muon-dir_z",
    "proton-p",
    "proton-dir_z",
    "tki-del_Tp",
    "tki-del_alpha",
    "tki-del_phi",
)
FAKE_DATA_TEST_SPECS: List[Tuple[str, Optional[str], Dict[str, Any]]] = [
    ("mec_test", "MEC Scale", {"scale_factor": 0.5}),
    ("qe_test", "QE Scale", {"scale_factor": 1.2}),
    ("q2_test_alpha_0.3", r"$Q^2$ tilt", {"scale_factor": 0.3}),
    ("costh_weight_scale_0.7", "Forward Muon Scale", {"scale_factor": 0.7}),
    ("proton_P_tilt_alpha_0.3", "Proton Momentum Tilt", {"scale_factor": 0.3}),
    ("bump_area0p5_center_bin", None, {"bump_area_bin_fraction": 0.5}),
    ("bump_area1_center_bin", None, {"bump_area_bin_fraction": 1.0}),
    ("bump_area1p5_center_bin", None, {"bump_area_bin_fraction": 1.5}),
]
GIBUU_SEL_MUP_DIR = (
    "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"
    "2026_09_09_142843__sel_mup-mc-GiBUU"
)
GIBUU_SEL_ALL_DIR = (
    "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"
    "2026_09_09_133316__sel_all-mc-CV-GiBUU"
)
DEFAULT_MC_LIST = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/unfolded/"
    "fake_data_tests/mc_files_used.txt"
)
DEFAULT_OUT = Path(
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/unfolded/fake_data_tests"
)


def table_var_configs() -> List[VariableConfig]:
    want = set(TABLE_VAR_NAMES)
    return [vc for vc in CORE_SELECTED_EVT_VARIABLE_CONFIGS if vc.var_save_name in want]


def test_label(var_config: VariableConfig, test_name: str, test_lbl, kwargs: Dict[str, Any]) -> str:
    if test_name.startswith("bump_"):
        return format_bump_test_label(
            var_config,
            kwargs["bump_area_bin_fraction"],
            bump_pos=kwargs.get("bump_pos"),
        )
    return str(test_lbl)


def _list_df_files(sample_dir: str, filename_str: str) -> List[str]:
    return sorted(glob.glob(str(Path(sample_dir) / f"*{filename_str}*.df")))


def _read_file_list(path: str) -> List[str]:
    files = [ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip()]
    return [p for p in files if Path(p).is_file()]


def _store_keys(fp: str) -> set:
    with pd.HDFStore(fp, "r") as st:
        return {k.lstrip("/") for k in st.keys()}


def _load_available(fp: str, keys: Sequence[str]) -> Dict[str, pd.DataFrame]:
    n_split = int(get_n_split(fp))
    have = _store_keys(fp)
    out: Dict[str, pd.DataFrame] = {}
    for key in keys:
        frames = []
        for i in range(n_split):
            name = f"{key}_{i}"
            if name in have:
                frames.append(pd.read_hdf(fp, key=name))
        if frames:
            out[key] = frames[0] if len(frames) == 1 else pd.concat(frames, axis=0, sort=False)
    return out


def _genweight(df: pd.DataFrame, n: int) -> np.ndarray:
    """GiBUU uses ``genweight``; Product B GENIE stores 0 → fall back to ones."""
    if n == 0:
        return np.zeros(0, dtype=np.float64)
    w = None
    try:
        if "mc" in df.columns.get_level_values(0):
            cand = np.asarray(df.mc.genweight, dtype=np.float64)
            if cand.shape[0] == n:
                w = cand
    except Exception:
        w = None
    if w is None:
        for c in df.columns:
            if "genweight" in str(c).lower():
                cand = np.asarray(df[c], dtype=np.float64)
                if cand.shape[0] == n:
                    w = cand
                    break
    if w is None:
        return np.ones(n, dtype=np.float64)
    w = np.where(np.isfinite(w), w, 0.0)
    if not np.any(w > 0):
        return np.ones(n, dtype=np.float64)
    return w


def _empty_var(nbins: int, test_names: Sequence[str]) -> Dict[str, Any]:
    z = np.zeros(nbins, dtype=np.float64)
    return {
        "n_bkg_reco": z.copy(),
        "tests": {tn: {"n_allsel_reco": z.copy(), "n_allmc": z.copy()} for tn in test_names},
    }


def _hist(df: pd.DataFrame, col, bins: np.ndarray, weights: np.ndarray, vsn: str) -> np.ndarray:
    if len(df) == 0:
        return np.zeros(len(bins) - 1, dtype=np.float64)
    tmp = df.copy()
    tmp["pot_weight"] = np.asarray(weights, dtype=np.float64)
    var, w = get_clipped_evts(tmp, col, bins, var_save_name=vsn)
    h, _ = np.histogram(var, bins=bins, weights=w)
    return np.asarray(h, dtype=np.float64)


def _accumulate_genie_one(fp: str, var_specs: Sequence[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    try:
        dfs = _load_available(fp, ["hdr", "evt", "mcnu"])
        if "hdr" not in dfs or "evt" not in dfs or "mcnu" not in dfs:
            return None
        pot = float(dfs["hdr"]["pot"].sum())
        evt, mcnu = _prepare_evt_mcnu(dfs["evt"], dfs["mcnu"])
        del dfs
        mcnu_sig = mcnu[IsNuInFV_NumuCC_1p0pi(mcnu, signal_truth_fv=SIGNAL_TRUTH_FV_UNFOLD)].copy()
        mcnu_sig.loc[:, "topo_categ"] = 1
        w_evt0 = _genweight(evt, len(evt))
        w_nu0 = _genweight(mcnu_sig, len(mcnu_sig))
        evt = evt.copy()
        evt["pot_weight"] = w_evt0
        mcnu_sig["pot_weight"] = w_nu0
        test_names = [t[0] for t in FAKE_DATA_TEST_SPECS]
        vc_map = _vc_by_name()
        out: Dict[str, Any] = {"pot": pot, "by_var": {}}
        for spec in var_specs:
            vsn = spec["var_save_name"]
            bins = np.asarray(spec["bins"], dtype=float)
            vc = vc_map[vsn]
            bkg = evt[evt.topo_categ != 1]
            pack = _empty_var(len(bins) - 1, test_names)
            pack["n_bkg_reco"] = _hist(bkg, spec["reco_col"], bins, _genweight(bkg, len(bkg)), vsn)
            fake_w = FakeDataWeights(evt, mcnu_sig, vc)
            for test_name, _, kwargs in FAKE_DATA_TEST_SPECS:
                # GENIE CAF codes only (QE=0, MEC=10). GiBUU is not reweighted here.
                w_evt, w_nu = fake_w.get_weights(test_name, generator="genie", **kwargs)
                pack["tests"][test_name]["n_allsel_reco"] = _hist(
                    evt, spec["reco_col"], bins, w_evt0 * np.asarray(w_evt, dtype=float), vsn
                )
                pack["tests"][test_name]["n_allmc"] = _hist(
                    mcnu_sig, spec["nu_col"], bins, w_nu0 * np.asarray(w_nu, dtype=float), vsn
                )
            out["by_var"][vsn] = pack
        del evt, mcnu, mcnu_sig
        gc.collect()
        return out
    except Exception as exc:
        print(f"  WARN GENIE skip {fp}: {exc}", flush=True)
        return None


def _accumulate_gibuu_one(fp: str, var_specs: Sequence[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    try:
        dfs = _load_available(fp, ["hdr", "evt", "mcnu"])
        if "hdr" not in dfs:
            return None
        pot = float(dfs["hdr"]["pot"].sum())
        out: Dict[str, Any] = {"pot": pot, "by_var": {}, "has_mcnu": "mcnu" in dfs, "has_evt": "evt" in dfs}
        if "evt" not in dfs or "mcnu" not in dfs:
            return out
        evt, mcnu = _prepare_evt_mcnu(dfs["evt"], dfs["mcnu"])
        del dfs
        mcnu_sig = mcnu[IsNuInFV_NumuCC_1p0pi(mcnu, signal_truth_fv=SIGNAL_TRUTH_FV_UNFOLD)]
        w_evt = _genweight(evt, len(evt))
        w_nu = _genweight(mcnu_sig, len(mcnu_sig))
        bkg = evt[evt.topo_categ != 1]
        w_bkg = _genweight(bkg, len(bkg))
        for spec in var_specs:
            vsn = spec["var_save_name"]
            bins = np.asarray(spec["bins"], dtype=float)
            out["by_var"][vsn] = {
                "n_allsel_reco": _hist(evt, spec["reco_col"], bins, w_evt, vsn),
                "n_bkg_reco": _hist(bkg, spec["reco_col"], bins, w_bkg, vsn),
                "n_allmc": _hist(mcnu_sig, spec["nu_col"], bins, w_nu, vsn),
            }
        del evt, mcnu, mcnu_sig
        gc.collect()
        return out
    except Exception as exc:
        print(f"  WARN GiBUU skip {fp}: {exc}", flush=True)
        return None


def _merge_genie(acc: Dict[str, Any], one: Dict[str, Any], var_specs, test_names) -> None:
    acc["pot"] += float(one["pot"])
    acc["n_files"] += 1
    for spec in var_specs:
        vsn = spec["var_save_name"]
        src = one["by_var"][vsn]
        dst = acc["by_var"][vsn]
        dst["n_bkg_reco"] += src["n_bkg_reco"]
        for tn in test_names:
            dst["tests"][tn]["n_allsel_reco"] += src["tests"][tn]["n_allsel_reco"]
            dst["tests"][tn]["n_allmc"] += src["tests"][tn]["n_allmc"]


def _merge_gibuu(acc: Dict[str, Any], one: Dict[str, Any], var_specs) -> None:
    acc["pot"] += float(one["pot"])
    acc["n_files"] += 1
    if not one.get("by_var"):
        acc["n_no_mcnu"] += 1
        return
    acc["n_with_mcnu"] += 1
    for spec in var_specs:
        vsn = spec["var_save_name"]
        src = one["by_var"][vsn]
        dst = acc["by_var"][vsn]
        for k in ("n_allsel_reco", "n_bkg_reco", "n_allmc"):
            dst[k] += src[k]


def _vc_by_name() -> Dict[str, VariableConfig]:
    return {vc.var_save_name: vc for vc in CORE_SELECTED_EVT_VARIABLE_CONFIGS}


def _var_specs_payload(var_configs: Sequence[VariableConfig]) -> List[Dict[str, Any]]:
    # Picklable (no VariableConfig objects) so ProcessPool workers can rebuild.
    return [
        {
            "var_save_name": vc.var_save_name,
            "bins": np.asarray(vc.bins, dtype=float),
            "reco_col": vc.var_evt_reco_col,
            "truth_col": vc.var_evt_truth_col,
            "nu_col": vc.var_nu_col,
        }
        for vc in var_configs
    ]


def _init_genie_acc(var_specs, test_names) -> Dict[str, Any]:
    return {
        "pot": 0.0,
        "n_files": 0,
        "n_skip": 0,
        "by_var": {spec["var_save_name"]: _empty_var(len(spec["bins"]) - 1, test_names) for spec in var_specs},
    }


def _init_gibuu_acc(var_specs) -> Dict[str, Any]:
    zpack = {
        spec["var_save_name"]: {
            "n_allsel_reco": np.zeros(len(spec["bins"]) - 1, dtype=np.float64),
            "n_bkg_reco": np.zeros(len(spec["bins"]) - 1, dtype=np.float64),
            "n_allmc": np.zeros(len(spec["bins"]) - 1, dtype=np.float64),
        }
        for spec in var_specs
    }
    return {"pot": 0.0, "n_files": 0, "n_with_mcnu": 0, "n_no_mcnu": 0, "by_var": zpack}


def _run_pool(fn, files: Sequence[str], var_specs, workers: int, desc: str):
    if workers <= 1:
        for fp in tqdm(files, desc=desc):
            yield fn(fp, var_specs)
        return
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(fn, fp, var_specs) for fp in files]
        for fut in tqdm(as_completed(futs), total=len(futs), desc=desc):
            yield fut.result()


def accumulate(
    *,
    genie_files: Sequence[str],
    gibuu_mup_files: Sequence[str],
    gibuu_all_files: Sequence[str],
    var_configs: Optional[Sequence[VariableConfig]] = None,
    workers: int = 16,
) -> Dict[str, Any]:
    var_configs = list(var_configs or table_var_configs())
    var_specs = _var_specs_payload(var_configs)
    test_names = [t[0] for t in FAKE_DATA_TEST_SPECS]

    genie = _init_genie_acc(var_specs, test_names)
    for one in _run_pool(_accumulate_genie_one, genie_files, var_specs, workers, "GENIE reweight"):
        if one is None:
            genie["n_skip"] += 1
            continue
        _merge_genie(genie, one, var_specs, test_names)

    gibuu = _init_gibuu_acc(var_specs)
    for one in _run_pool(_accumulate_gibuu_one, gibuu_mup_files, var_specs, workers, "GiBUU sel_mup"):
        if one is None:
            continue
        _merge_gibuu(gibuu, one, var_specs)

    all_pot = 0.0
    n_all = 0
    for fp in tqdm(gibuu_all_files, desc="GiBUU sel_all POT"):
        try:
            dfs = _load_available(fp, ["hdr"])
            if "hdr" in dfs:
                all_pot += float(dfs["hdr"]["pot"].sum())
                n_all += 1
        except Exception:
            continue

    return {
        "genie": genie,
        "gibuu": gibuu,
        "gibuu_sel_all_pot": all_pot,
        "gibuu_sel_all_n_files": n_all,
        "test_names": test_names,
        "variables": [vc.var_save_name for vc in var_configs],
    }


def save_acc(acc: Dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: Dict[str, Any] = {
        "gibuu_sel_all_pot": np.float64(acc["gibuu_sel_all_pot"]),
        "gibuu_sel_all_n_files": np.int64(acc["gibuu_sel_all_n_files"]),
        "test_names": np.array(acc["test_names"], dtype=object),
        "variables": np.array(acc["variables"], dtype=object),
        "genie_meta": np.array(
            [{"pot": acc["genie"]["pot"], "n_files": acc["genie"]["n_files"], "n_skip": acc["genie"]["n_skip"]}],
            dtype=object,
        ),
        "gibuu_meta": np.array(
            [
                {
                    "pot": acc["gibuu"]["pot"],
                    "n_files": acc["gibuu"]["n_files"],
                    "n_with_mcnu": acc["gibuu"]["n_with_mcnu"],
                    "n_no_mcnu": acc["gibuu"]["n_no_mcnu"],
                }
            ],
            dtype=object,
        ),
    }
    for vsn, pack in acc["genie"]["by_var"].items():
        payload[f"genie::{vsn}::n_bkg_reco"] = pack["n_bkg_reco"]
        for tn, tpack in pack["tests"].items():
            payload[f"genie::{vsn}::{tn}::n_allsel_reco"] = tpack["n_allsel_reco"]
            payload[f"genie::{vsn}::{tn}::n_allmc"] = tpack["n_allmc"]
    for vsn, pack in acc["gibuu"]["by_var"].items():
        for k in ("n_allsel_reco", "n_bkg_reco", "n_allmc"):
            payload[f"gibuu::{vsn}::{k}"] = pack[k]
    np.savez_compressed(path, **payload)
    print("wrote", path)


def load_acc(path: Path) -> Dict[str, Any]:
    blob = np.load(path, allow_pickle=True)
    test_names = [str(x) for x in blob["test_names"]]
    variables = [str(x) for x in blob["variables"]]
    gmeta = blob["genie_meta"][0]
    bmeta = blob["gibuu_meta"][0]
    genie = {
        "pot": float(gmeta["pot"]),
        "n_files": int(gmeta["n_files"]),
        "n_skip": int(gmeta["n_skip"]),
        "by_var": {},
    }
    gibuu = {
        "pot": float(bmeta["pot"]),
        "n_files": int(bmeta["n_files"]),
        "n_with_mcnu": int(bmeta["n_with_mcnu"]),
        "n_no_mcnu": int(bmeta["n_no_mcnu"]),
        "by_var": {},
    }
    for vsn in variables:
        genie["by_var"][vsn] = {
            "n_bkg_reco": np.asarray(blob[f"genie::{vsn}::n_bkg_reco"], dtype=float),
            "tests": {},
        }
        for tn in test_names:
            genie["by_var"][vsn]["tests"][tn] = {
                "n_allsel_reco": np.asarray(blob[f"genie::{vsn}::{tn}::n_allsel_reco"], dtype=float),
                "n_allmc": np.asarray(blob[f"genie::{vsn}::{tn}::n_allmc"], dtype=float),
            }
        gibuu["by_var"][vsn] = {
            k: np.asarray(blob[f"gibuu::{vsn}::{k}"], dtype=float)
            for k in ("n_allsel_reco", "n_bkg_reco", "n_allmc")
        }
    return {
        "genie": genie,
        "gibuu": gibuu,
        "gibuu_sel_all_pot": float(blob["gibuu_sel_all_pot"]),
        "gibuu_sel_all_n_files": int(blob["gibuu_sel_all_n_files"]),
        "test_names": test_names,
        "variables": variables,
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mc-list", default=DEFAULT_MC_LIST)
    p.add_argument("--mc-dir", default=dmo.PRODUCT_B_MC_DIR)
    p.add_argument("--mc-fn", default=dmo.PRODUCT_B_MC_FN)
    p.add_argument("--gibuu-mup-dir", default=GIBUU_SEL_MUP_DIR)
    p.add_argument("--gibuu-all-dir", default=GIBUU_SEL_ALL_DIR)
    p.add_argument("--out-dir", default=str(DEFAULT_OUT))
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--max-files", type=int, default=0)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.mc_list and Path(args.mc_list).is_file():
        genie_files = _read_file_list(args.mc_list)
        print(f"GENIE from list: {len(genie_files)}  {args.mc_list}")
    else:
        genie_files = _list_df_files(args.mc_dir, args.mc_fn)
        print(f"GENIE from glob: {len(genie_files)}  {args.mc_dir}")
    gibuu_mup = _list_df_files(args.gibuu_mup_dir, "sel_mup-mc-GiBUU")
    gibuu_all = _list_df_files(args.gibuu_all_dir, "sel_all-mc-CV-GiBUU")
    if args.max_files > 0:
        genie_files = genie_files[: args.max_files]
        gibuu_mup = gibuu_mup[: args.max_files]
        gibuu_all = gibuu_all[: args.max_files]
    print(f"GiBUU sel_mup={len(gibuu_mup)}  sel_all={len(gibuu_all)}")
    acc = accumulate(
        genie_files=genie_files,
        gibuu_mup_files=gibuu_mup,
        gibuu_all_files=gibuu_all,
        workers=args.workers,
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    save_acc(acc, out_dir / "fake_data_histograms.npz")
    meta = {
        "n_genie": acc["genie"]["n_files"],
        "n_genie_skip": acc["genie"]["n_skip"],
        "genie_pot": acc["genie"]["pot"],
        "n_gibuu_mup": acc["gibuu"]["n_files"],
        "n_gibuu_mup_mcnu": acc["gibuu"]["n_with_mcnu"],
        "gibuu_mup_pot": acc["gibuu"]["pot"],
        "gibuu_sel_all_pot": acc["gibuu_sel_all_pot"],
        "gibuu_sel_all_n_files": acc["gibuu_sel_all_n_files"],
    }
    (out_dir / "accumulate_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
