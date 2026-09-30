#!/usr/bin/env python3
"""Build stacked constraint joints, check diagonal blocks vs PRL, run two tests.

May 14 GENIE chunks and May 15 pair joints are **read-only**. All outputs go under
``PRL/constraint_validation_mu_p_dirz/`` (never into May campaign directories).

Recipe
------
* Reconstruct May MEC interpolator + SBN_v1 MEC stacked universes by concatenating
  the four constraint variables from the May 14 marginal chunks (same universe index).
* Build a Product B–diagonal hybrid joint: block-diagonal PRL ``flux`` + ``g4`` +
  ``genie_rate`` (overlay μ_nu) plus May MEC **cross-variable** blocks only.
  Detector/cosmics/POT/ntargets stay block-diagonal from CategorySummary at load time.
* **Gate:** diagonal blocks of that constructed joint vs PRL CategorySummary (same
  categories the constraint loader uses). May-archive per-knob rate is diagnostic only.
* If the PRL check passes, write ``PRL/constraint_validation_mu_p_dirz/`` and run
  (muon-p, muon-dir_z) → proton-p and → proton-dir_z.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
from datetime import datetime, timezone
from glob import glob
from pathlib import Path
from typing import Any, Dict, Mapping, Sequence, Tuple

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from pyanalib.covariance import corr_from_fraccov, cov_from_fraccov, fraccov_from_cov  # noqa: E402

from analysis_village.numucc_1p0pi.dataset_locations import (  # noqa: E402
    prl_overlay_counts_npz,
    prl_syst_disk_root,
)
from analysis_village.numucc_1p0pi.scripts.conditional_constraint_validation import (  # noqa: E402
    chi2_and_pull,
    data_mc_chi2_ndof,
    load_overlay_constraint_hists,
    sigma_yy_with_data_poisson_on_y,
)
from analysis_village.numucc_1p0pi.scripts.conditional_constraint_validation import (  # noqa: E402
    conditional_gaussian_update,
)
from analysis_village.numucc_1p0pi.syst_category_summary import (  # noqa: E402
    CAT_COSMICS,
    CAT_DETECTOR,
    CAT_FLUX,
    CAT_G4,
    CAT_GENIE_RATE,
    CAT_NTARGETS,
    CAT_POT,
    category_cov_frac,
    load_category_syst_summary,
    total_cov_frac,
)
from analysis_village.numucc_1p0pi.scripts._genie_pkl_compat import (  # noqa: E402
    load_cov_mat_dict_compat,
)
from analysis_village.numucc_1p0pi.syst_cc_joint_multisim_common import (  # noqa: E402
    JOINT_CC_EXTRAS_CATEGORIES,
    STACK_SLUG,
    joint_stack_meta,
)
from analysis_village.numucc_1p0pi.syst_disk_cc_layout import (  # noqa: E402
    FILE_JOINT_GENIE_COMBINED,
    joint_genie_out_dir,
    joint_multisim_category_inner_key,
    joint_multisim_category_npz_basename,
    joint_multisim_category_out_dir,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import (  # noqa: E402
    FILE_GENIE,
    SUB_GENIE,
    category_summary_npz_path,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402

# ---------------------------------------------------------------------------
# Protected source trees (never write / delete / truncate)
# ---------------------------------------------------------------------------
PROTECTED_MAY_ROOTS = (
    Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/genie_syst-chunked-20260514"),
    Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/joint_genie_cc-chunked-20260513"),
    Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/joint_genie_cc-chunked-20260514"),
    Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/joint_genie_cc-chunked-20260515"),
    Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final-archive"),
    Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-genie-final"),
)

MAY_AR23P_CHUNKS = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "genie_syst-chunked-20260514/chunks/genie__Ar23p__*.pkl"
)
MAY_MEC_CHUNKS = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "genie_syst-chunked-20260514/chunks/genie__MEC__*.pkl"
)
MAY_ARCHIVE_PKL = Path(
    "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/"
    "systematics-final-archive/GENIE/cov_mat_dict.pkl"
)

STACK_VARS = (
    VariableConfig.muon_momentum(),
    VariableConfig.muon_direction(),
    VariableConfig.proton_momentum(),
    VariableConfig.proton_direction(),
)
STACK_NAMES = tuple(v.var_save_name for v in STACK_VARS)
STACK_SIZES = tuple(int(len(v.bin_centers)) for v in STACK_VARS)
STACK_OFFSETS = tuple(int(s) for s in np.cumsum((0,) + STACK_SIZES[:-1]))

MAY_INTERP_KNOBS = tuple(
    f"MECq0q3InterpWeighting_SuSAv2To{_to}_q0binned_MECResponse_q0bin{_b}"
    for _to in ("Valenica", "Martini")
    for _b in range(4)
)
MAY_SBNV1_MEC_KNOBS = (
    "GENIEReWeight_SBN_v1_multisim_NormCCMEC",
    "GENIEReWeight_SBN_v1_multisim_NormNCMEC",
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",
)
MAY_MEC_KNOBS = MAY_INTERP_KNOBS + MAY_SBNV1_MEC_KNOBS

# PRL neutrino + extras used by the stacked joint loader (no MCstat).
PRL_JOINT_CATS = (CAT_FLUX, CAT_G4, CAT_GENIE_RATE, CAT_DETECTOR, CAT_COSMICS, CAT_POT, CAT_NTARGETS)
PRL_NEUTRINO_CATS = (CAT_FLUX, CAT_G4, CAT_GENIE_RATE)

# Relative RMS of (diag_joint - diag_ref) / max(diag_ref, floor) for a pass.
DIAG_MATCH_TOL = 0.12


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _assert_not_protected(path: Path) -> None:
    resolved = path.expanduser().resolve()
    for root in PROTECTED_MAY_ROOTS:
        try:
            resolved.relative_to(root.resolve())
        except ValueError:
            continue
        raise RuntimeError("refusing to write under protected May tree: %s" % resolved)


def _slices() -> Dict[str, slice]:
    return {
        name: slice(off, off + sz)
        for name, off, sz in zip(STACK_NAMES, STACK_OFFSETS, STACK_SIZES)
    }


def _cov_from_univ(univ: np.ndarray, cv: np.ndarray) -> Dict[str, np.ndarray]:
    u = np.asarray(univ, dtype=np.float64)
    c = np.asarray(cv, dtype=np.float64).reshape(-1)
    if u.ndim != 2 or u.shape[1] != c.size:
        raise ValueError("univ shape %s vs cv %s" % (u.shape, c.shape))
    n_u = int(u.shape[0])
    d = u - c
    cov = (d.T @ d) / max(n_u, 1)
    frac = fraccov_from_cov(cov, c)
    corr = corr_from_fraccov(frac)
    return {"cov": cov, "cov_frac": frac, "corr": corr, "n_univ": n_u}


def _sum_rate_acc(paths: Sequence[str], knobs: Sequence[str], var_names: Sequence[str]) -> Dict[str, Dict[str, Any]]:
    """Sum May 14 ``rate_univ_cv[knob][var]`` across chunk pickles (read-only)."""
    acc: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {
        k: {v: None for v in var_names} for k in knobs  # type: ignore[misc]
    }
    n_files = 0
    for p in paths:
        with open(p, "rb") as f:
            blob = pickle.load(f)
        rate = blob.get("rate_univ_cv") or {}
        n_files += 1
        for kn in knobs:
            if kn not in rate:
                continue
            for vn in var_names:
                slot = rate[kn].get(vn)
                if not isinstance(slot, dict) or "univ" not in slot:
                    continue
                u = np.asarray(slot["univ"], dtype=np.float64)
                c = np.asarray(slot["cv"], dtype=np.float64).reshape(-1)
                cur = acc[kn][vn]
                if cur is None:
                    acc[kn][vn] = {"univ": u.copy(), "cv": c.copy()}
                else:
                    if cur["univ"].shape != u.shape:
                        raise ValueError(
                            "shape mismatch %s %s: %s vs %s" % (kn, vn, cur["univ"].shape, u.shape)
                        )
                    cur["univ"] += u
                    cur["cv"] += c
    out: Dict[str, Dict[str, Any]] = {}
    for kn in knobs:
        missing = [v for v in var_names if acc[kn][v] is None]
        if missing:
            continue
        us = [acc[kn][v]["univ"] for v in var_names]
        cs = [acc[kn][v]["cv"] for v in var_names]
        u_j = np.hstack(us)
        c_j = np.concatenate(cs)
        pack = _cov_from_univ(u_j, c_j)
        pack["univ"] = u_j
        pack["cv"] = c_j
        pack["n_files"] = n_files
        out[kn] = pack
    return out


def _zero_diag_blocks(mat: np.ndarray) -> np.ndarray:
    """Keep only cross-variable blocks of a stacked matrix."""
    out = np.array(mat, dtype=np.float64, copy=True)
    for sl in _slices().values():
        out[sl, sl] = 0.0
    return out


def _extract_blocks(mat: np.ndarray) -> Dict[str, np.ndarray]:
    sl = _slices()
    return {name: np.array(mat[sl[name], sl[name]], copy=True) for name in STACK_NAMES}


def _rel_rms_diag(a: np.ndarray, b: np.ndarray) -> float:
    da = np.diag(np.asarray(a, dtype=float))
    db = np.diag(np.asarray(b, dtype=float))
    denom = np.maximum(np.abs(db), 1e-30)
    return float(np.sqrt(np.mean(((da - db) / denom) ** 2)))


def _max_abs_rel_diag(a: np.ndarray, b: np.ndarray) -> float:
    da = np.diag(np.asarray(a, dtype=float))
    db = np.diag(np.asarray(b, dtype=float))
    denom = np.maximum(np.abs(db), 1e-30)
    return float(np.max(np.abs((da - db) / denom)))


def _diag_cmp(a: np.ndarray, b: np.ndarray, tol: float) -> Dict[str, Any]:
    rel = _rel_rms_diag(a, b)
    mx = _max_abs_rel_diag(a, b)
    return {"rel_rms_diag": rel, "max_abs_rel_diag": mx, "pass": rel <= tol}


def _blkdiag_frac(fracs: Sequence[np.ndarray]) -> np.ndarray:
    n = int(sum(m.shape[0] for m in fracs))
    out = np.zeros((n, n), dtype=np.float64)
    o = 0
    for m in fracs:
        k = int(m.shape[0])
        out[o : o + k, o : o + k] = np.asarray(m, dtype=np.float64)
        o += k
    return out


def _write_stacked_npz(
    path: Path,
    inner_key: str,
    frac: np.ndarray,
    cov: np.ndarray,
    *,
    extra_meta: Mapping[str, Any] | None = None,
    by_knob: Mapping[str, Mapping[str, np.ndarray]] | None = None,
) -> None:
    _assert_not_protected(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = joint_stack_meta(STACK_VARS, bkgd_subtract=False)
    if extra_meta:
        meta = {**meta, **dict(extra_meta)}
    cell: Dict[str, Any] = {
        inner_key: {
            "cov_frac": np.asarray(frac, dtype=float),
            "cov": np.asarray(cov, dtype=float),
            "corr": corr_from_fraccov(frac),
        },
        "meta": meta,
    }
    if by_knob:
        cell["%s_by_knob" % inner_key] = {
            kn: {
                "cov_frac": np.asarray(p["cov_frac"], dtype=float),
                "cov": np.asarray(p["cov"], dtype=float),
                "corr": np.asarray(p.get("corr", corr_from_fraccov(p["cov_frac"])), dtype=float),
            }
            for kn, p in by_knob.items()
        }
    np.savez_compressed(path, **{STACK_SLUG: np.array(cell, dtype=object)})


def _overlay_vectors(counts_npz: Path) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    mu_nu = {}
    mu_tot = {}
    n_data = {}
    for v in STACK_VARS:
        h = load_overlay_constraint_hists(counts_npz, v)
        mu_nu[v.var_save_name] = np.asarray(h["mu_nu"], dtype=float)
        mu_tot[v.var_save_name] = np.asarray(h["mu_total"], dtype=float)
        n_data[v.var_save_name] = np.asarray(h["n_data"], dtype=float)
        if mu_nu[v.var_save_name].size != len(v.bin_centers):
            raise ValueError(
                "overlay %s has %d bins, VariableConfig has %d"
                % (v.var_save_name, mu_nu[v.var_save_name].size, len(v.bin_centers))
            )
    return mu_nu, mu_tot, n_data


def default_out_root() -> Path:
    return Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/constraint_validation_mu_p_dirz")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-root", type=Path, default=default_out_root())
    p.add_argument("--diag-tol", type=float, default=DIAG_MATCH_TOL)
    p.add_argument("--skip-tests", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out_root = args.out_root.expanduser().resolve()
    _assert_not_protected(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "README.txt").write_text(
        "Conditional-constraint test: muon-p + muon-dir_z constrain proton-p / proton-dir_z.\n"
        "May 14/15 GENIE sources are read-only; this directory is the only write target.\n"
    )

    prl_root = prl_syst_disk_root("B")
    counts_npz = prl_overlay_counts_npz("B")
    mu_nu, mu_tot, n_data = _overlay_vectors(counts_npz)
    mu_nu_stack = np.concatenate([mu_nu[n] for n in STACK_NAMES])
    mu_tot_stack = np.concatenate([mu_tot[n] for n in STACK_NAMES])

    summary_path = category_summary_npz_path(str(prl_root))
    cat_sum = load_category_syst_summary(summary_path)

    prl_abs: Dict[str, Dict[str, np.ndarray]] = {n: {} for n in STACK_NAMES}
    prl_frac: Dict[str, Dict[str, np.ndarray]] = {n: {} for n in STACK_NAMES}
    for name, vc in zip(STACK_NAMES, STACK_VARS):
        for cat in PRL_JOINT_CATS:
            frac = category_cov_frac(cat_sum, name, cat)
            prl_frac[name][cat] = frac
            prl_abs[name][cat] = cov_from_fraccov(frac, mu_tot[name] if cat in (CAT_DETECTOR, CAT_COSMICS, CAT_POT, CAT_NTARGETS) else mu_nu[name])
        prl_frac[name]["total_rate"] = total_cov_frac(cat_sum, name, kind="rate")
        prl_abs[name]["total_rate"] = cov_from_fraccov(prl_frac[name]["total_rate"], mu_tot[name])

    # ---- May 14 stacked joints (read-only globs) ----
    ar23p_paths = sorted(glob(MAY_AR23P_CHUNKS))
    mec_paths = sorted(glob(MAY_MEC_CHUNKS))
    if len(ar23p_paths) < 1 or len(mec_paths) < 1:
        raise SystemExit("May 14 GENIE chunks missing; refusing to continue (sources must stay intact).")
    print("[may] Ar23p chunks %d  MEC chunks %d (read-only)" % (len(ar23p_paths), len(mec_paths)), flush=True)
    may_by_knob = {}
    may_by_knob.update(_sum_rate_acc(ar23p_paths, MAY_INTERP_KNOBS, STACK_NAMES))
    may_by_knob.update(_sum_rate_acc(mec_paths, MAY_SBNV1_MEC_KNOBS, STACK_NAMES))
    missing_knobs = [k for k in MAY_MEC_KNOBS if k not in may_by_knob]
    if missing_knobs:
        raise SystemExit("May chunks missing knobs: %s" % missing_knobs)

    may_frac = np.sum([may_by_knob[k]["cov_frac"] for k in MAY_MEC_KNOBS], axis=0)
    may_cov_univ = np.sum([may_by_knob[k]["cov"] for k in MAY_MEC_KNOBS], axis=0)
    # Scale May *fractional* cov to overlay neutrino-MC counts (Product B recipe).
    may_cov_overlay = cov_from_fraccov(may_frac, mu_nu_stack)
    may_blocks = _extract_blocks(may_cov_overlay)

    # ---- May archive per-knob rate (diagnostic; not a save gate) ----
    print("[may] loading archive pkl (read-only)", flush=True)
    archive = load_cov_mat_dict_compat(MAY_ARCHIVE_PKL)
    archive_cmp: Dict[str, Any] = {}
    for name in STACK_NAMES:
        row = archive[name]
        acc = None
        for kn in MAY_MEC_KNOBS:
            if kn + "_rate" in row:
                mat = np.asarray(row[kn + "_rate"], dtype=float)
            elif kn in row:
                mat = np.asarray(row[kn], dtype=float)
            else:
                raise KeyError("archive %s missing %s or %s_rate" % (name, kn, kn))
            abs_m = cov_from_fraccov(mat, mu_nu[name]) if mat.shape == (mu_nu[name].size, mu_nu[name].size) else mat
            acc = abs_m if acc is None else acc + abs_m
        archive_cmp[name] = _diag_cmp(may_blocks[name], acc, args.diag_tol)
        print("[diag] May stacked vs archive %s  rel_rms=%.4f  max=%.4f  (informational)" % (
            name, archive_cmp[name]["rel_rms_diag"], archive_cmp[name]["max_abs_rel_diag"]
        ), flush=True)

    hybrid_genie_abs = _blkdiag_from_prl(prl_abs, CAT_GENIE_RATE) + _zero_diag_blocks(may_cov_overlay)
    hybrid_flux_abs = _blkdiag_from_prl(prl_abs, CAT_FLUX)
    hybrid_g4_abs = _blkdiag_from_prl(prl_abs, CAT_G4)
    extras_abs = None
    for cat in JOINT_CC_EXTRAS_CATEGORIES:
        extras_abs = _blkdiag_from_prl(prl_abs, cat) if extras_abs is None else extras_abs + _blkdiag_from_prl(prl_abs, cat)

    prl_pkl = Path(prl_root) / SUB_GENIE / FILE_GENIE
    prl_genie_row = load_cov_mat_dict_compat(prl_pkl) if prl_pkl.is_file() else None

    prl_cmp: Dict[str, Any] = {}
    prl_full_cmp: Dict[str, Any] = {}
    prl_pkl_cmp: Dict[str, Any] = {}
    prl_ok = True
    sl = _slices()
    for name in STACK_NAMES:
        prl_nu = prl_abs[name][CAT_FLUX] + prl_abs[name][CAT_G4] + prl_abs[name][CAT_GENIE_RATE]
        hyb_nu = (
            hybrid_flux_abs[sl[name], sl[name]]
            + hybrid_g4_abs[sl[name], sl[name]]
            + hybrid_genie_abs[sl[name], sl[name]]
        )
        prl_cmp[name] = _diag_cmp(hyb_nu, prl_nu, args.diag_tol)
        prl_full = prl_nu
        for cat in JOINT_CC_EXTRAS_CATEGORIES:
            prl_full = prl_full + prl_abs[name][cat]
        hyb_full = hyb_nu + extras_abs[sl[name], sl[name]]
        prl_full_cmp[name] = _diag_cmp(hyb_full, prl_full, args.diag_tol)
        prl_ok = prl_ok and prl_cmp[name]["pass"] and prl_full_cmp[name]["pass"]
        print("[check] hybrid vs PRL CategorySummary %s  neutrino rel_rms=%.4e  full rel_rms=%.4e  %s" % (
            name,
            prl_cmp[name]["rel_rms_diag"],
            prl_full_cmp[name]["rel_rms_diag"],
            "PASS" if prl_cmp[name]["pass"] and prl_full_cmp[name]["pass"] else "FAIL",
        ), flush=True)
        if prl_genie_row is not None and name in prl_genie_row and "genie_rate" in prl_genie_row[name]:
            pkl_frac = np.asarray(prl_genie_row[name]["genie_rate"], dtype=float)
            pkl_abs = cov_from_fraccov(pkl_frac, mu_nu[name])
            prl_pkl_cmp[name] = _diag_cmp(hybrid_genie_abs[sl[name], sl[name]], pkl_abs, args.diag_tol)
            prl_ok = prl_ok and prl_pkl_cmp[name]["pass"]
            print("[check] hybrid GENIE vs PRL GENIE/cov_mat_dict.pkl %s  rel_rms=%.4e  %s" % (
                name, prl_pkl_cmp[name]["rel_rms_diag"], "PASS" if prl_pkl_cmp[name]["pass"] else "FAIL"
            ), flush=True)

    report = {
        "created_utc": _utc(),
        "protected_may_roots": [str(p) for p in PROTECTED_MAY_ROOTS],
        "may_ar23p_chunks": len(ar23p_paths),
        "may_mec_chunks": len(mec_paths),
        "stack_names": list(STACK_NAMES),
        "stack_sizes": list(STACK_SIZES),
        "diag_tol": args.diag_tol,
        "may_vs_archive_informational": archive_cmp,
        "hybrid_neutrino_vs_prl": prl_cmp,
        "hybrid_full_vs_prl": prl_full_cmp,
        "hybrid_genie_vs_prl_pkl": prl_pkl_cmp,
        "prl_diag_ok": prl_ok,
        "overlay_counts_npz": str(counts_npz),
        "prl_syst_root": str(prl_root),
        "may_mec_knobs": list(MAY_MEC_KNOBS),
        "joint_extras_categories": list(JOINT_CC_EXTRAS_CATEGORIES),
    }
    (out_root / "diagonal_check.json").write_text(json.dumps(report, indent=2) + "\n")

    if not prl_ok:
        print("[stop] PRL diagonal check failed; not writing JointCC / not running tests", flush=True)
        print(json.dumps(report, indent=2))
        return 2

    # ---- write hybrid JointCC (new PRL subdirectory only) ----
    cc_root = out_root / "JointCC"
    _assert_not_protected(cc_root)
    flux_frac = fraccov_from_cov(hybrid_flux_abs, mu_nu_stack)
    g4_frac = fraccov_from_cov(hybrid_g4_abs, mu_nu_stack)
    genie_frac = fraccov_from_cov(hybrid_genie_abs, mu_nu_stack)
    _write_stacked_npz(
        Path(joint_multisim_category_out_dir(str(cc_root), "Flux")) / joint_multisim_category_npz_basename("Flux"),
        joint_multisim_category_inner_key("Flux"),
        flux_frac,
        hybrid_flux_abs,
        extra_meta={"recipe": "PRL flux block-diagonal; no joint off-diag (campaign still mapping)"},
    )
    _write_stacked_npz(
        Path(joint_multisim_category_out_dir(str(cc_root), "G4")) / joint_multisim_category_npz_basename("G4"),
        joint_multisim_category_inner_key("G4"),
        g4_frac,
        hybrid_g4_abs,
        extra_meta={"recipe": "PRL G4 block-diagonal; no joint off-diag (campaign still mapping)"},
    )
    _write_stacked_npz(
        Path(joint_genie_out_dir(str(cc_root))) / FILE_JOINT_GENIE_COMBINED,
        "JointGenie",
        genie_frac,
        hybrid_genie_abs,
        extra_meta={
            "recipe": "PRL genie_rate block-diagonal + May MEC interpolator/SBN_v1 cross-variable joints",
            "may_chunk_globs": [MAY_AR23P_CHUNKS, MAY_MEC_CHUNKS],
        },
        by_knob={k: {"cov_frac": may_by_knob[k]["cov_frac"], "cov": cov_from_fraccov(may_by_knob[k]["cov_frac"], mu_nu_stack), "corr": may_by_knob[k]["corr"]} for k in MAY_MEC_KNOBS},
    )
    np.savez_compressed(
        out_root / "may_mec_stacked.npz",
        cov_frac=may_frac,
        cov_overlay=may_cov_overlay,
        cov_univ_native=may_cov_univ,
        mu_nu=mu_nu_stack,
        var_save_names=np.array(STACK_NAMES, dtype=object),
        n_bins_per_var=np.array(STACK_SIZES),
    )

    from analysis_village.numucc_1p0pi.cc_joint_cov import (  # noqa: PLC0415
        build_joint_multi_covariance_abs,
        split_joint_multi_sigma,
    )

    os.environ["NUMUCC_SYST_DISK_ROOT"] = str(prl_root)
    os.environ["NUMUCC_SYST_DISK_CC_ROOT"] = str(cc_root)

    # Loader-level diagonal check: permute stacked joint to each test layout and
    # compare XX / YY variable blocks to PRL CategorySummary.
    tests = (
        ("proton_p", VariableConfig.proton_momentum()),
        ("proton_dir_z", VariableConfig.proton_direction()),
    )
    ys = (VariableConfig.muon_momentum(), VariableConfig.muon_direction())
    loader_cmp: Dict[str, Any] = {}
    loader_ok = True
    for tag, var_X in tests:
        hx = load_overlay_constraint_hists(counts_npz, var_X)
        hy = [load_overlay_constraint_hists(counts_npz, vy) for vy in ys]
        mu_X = hx["mu_total"]
        mu_nu_X = hx["mu_nu"]
        mu_Ys = tuple(h["mu_total"] for h in hy)
        mu_nu_Ys = tuple(h["mu_nu"] for h in hy)
        sigma_joint = build_joint_multi_covariance_abs(
            var_X,
            ys,
            mu_X,
            mu_Ys,
            syst_cc_root=str(cc_root),
            syst_marginal_root=str(prl_root),
            mu_nu_X=mu_nu_X,
            mu_nu_Ys=mu_nu_Ys,
        )
        sigma_XX, _sigma_XY, sigma_YY, _ = split_joint_multi_sigma(
            sigma_joint, len(mu_X), [len(m) for m in mu_Ys]
        )
        prl_X = prl_abs[var_X.var_save_name][CAT_FLUX] + prl_abs[var_X.var_save_name][CAT_G4] + prl_abs[var_X.var_save_name][CAT_GENIE_RATE]
        for cat in JOINT_CC_EXTRAS_CATEGORIES:
            prl_X = prl_X + prl_abs[var_X.var_save_name][cat]
        rec = {"XX": _diag_cmp(sigma_XX, prl_X, args.diag_tol), "YY": {}}
        y_off = 0
        for vy, mu_y in zip(ys, mu_Ys):
            n_y = int(mu_y.size)
            blk = sigma_YY[y_off : y_off + n_y, y_off : y_off + n_y]
            prl_y = prl_abs[vy.var_save_name][CAT_FLUX] + prl_abs[vy.var_save_name][CAT_G4] + prl_abs[vy.var_save_name][CAT_GENIE_RATE]
            for cat in JOINT_CC_EXTRAS_CATEGORIES:
                prl_y = prl_y + prl_abs[vy.var_save_name][cat]
            rec["YY"][vy.var_save_name] = _diag_cmp(blk, prl_y, args.diag_tol)
            y_off += n_y
        loader_ok = loader_ok and rec["XX"]["pass"] and all(v["pass"] for v in rec["YY"].values())
        loader_cmp[tag] = rec
        print(
            "[check] loaded joint %s  XX(%s) rel_rms=%.4e  YY %s  %s"
            % (
                tag,
                var_X.var_save_name,
                rec["XX"]["rel_rms_diag"],
                {k: "%.4e" % v["rel_rms_diag"] for k, v in rec["YY"].items()},
                "PASS" if rec["XX"]["pass"] and all(v["pass"] for v in rec["YY"].values()) else "FAIL",
            ),
            flush=True,
        )
    report["loader_diag_vs_prl"] = loader_cmp
    report["loader_diag_ok"] = loader_ok
    (out_root / "diagonal_check.json").write_text(json.dumps(report, indent=2) + "\n")
    if not loader_ok:
        print("[stop] loaded-joint diagonal blocks do not match PRL; not running tests", flush=True)
        print(json.dumps(report, indent=2))
        return 2

    if args.skip_tests:
        print("[ok] matrices saved under %s (tests skipped)" % out_root, flush=True)
        return 0

    # ---- run the two constraint tests ----
    run_summaries = []
    for tag, var_X in tests:
        run_dir = out_root / "runs" / tag
        run_dir.mkdir(parents=True, exist_ok=True)
        hx = load_overlay_constraint_hists(counts_npz, var_X)
        hy = [load_overlay_constraint_hists(counts_npz, vy) for vy in ys]
        n_X = hx["n_data"]
        mu_X = hx["mu_total"]
        mu_nu_X = hx["mu_nu"]
        n_Ys = tuple(h["n_data"] for h in hy)
        mu_Ys = tuple(h["mu_total"] for h in hy)
        mu_nu_Ys = tuple(h["mu_nu"] for h in hy)
        sigma_joint = build_joint_multi_covariance_abs(
            var_X,
            ys,
            mu_X,
            mu_Ys,
            syst_cc_root=str(cc_root),
            syst_marginal_root=str(prl_root),
            mu_nu_X=mu_nu_X,
            mu_nu_Ys=mu_nu_Ys,
        )
        sigma_XX, sigma_XY, sigma_YY, _ = split_joint_multi_sigma(
            sigma_joint, len(mu_X), [len(m) for m in mu_Ys]
        )
        n_Y = np.concatenate(n_Ys)
        mu_Y = np.concatenate(mu_Ys)
        sigma_YY_eff = sigma_yy_with_data_poisson_on_y(sigma_YY, n_Y)
        mu_X_c, sigma_XX_c, _gain = conditional_gaussian_update(
            mu_X, mu_Y, n_Y, sigma_XX, sigma_XY, sigma_YY_eff
        )
        chi2_pre, ndof_pre, chi2_pre_per = data_mc_chi2_ndof(n_X, mu_X, sigma_XX)
        chi2_post, pval, ndof, pull = chi2_and_pull(n_X, mu_X_c, sigma_XX_c)
        rec = {
            "target": var_X.var_save_name,
            "constrain_with": [v.var_save_name for v in ys],
            "chi2_prefit": float(chi2_pre),
            "chi2_prefit_per_ndof": float(chi2_pre_per),
            "ndof_prefit": int(ndof_pre),
            "chi2_postfit": float(chi2_post),
            "chi2_postfit_per_ndof": float(chi2_post / max(ndof, 1)),
            "p_value": float(pval),
            "ndof": int(ndof),
            "n_bins_X": int(n_X.size),
            "n_bins_Y": [int(m.size) for m in mu_Ys],
            "frobenius_joint": float(np.linalg.norm(sigma_joint, ord="fro")),
            "pull": np.asarray(pull, dtype=float).tolist(),
        }
        (run_dir / "summary.json").write_text(json.dumps(rec, indent=2) + "\n")
        np.savez_compressed(
            run_dir / "joint_blocks.npz",
            sigma_joint=sigma_joint,
            sigma_XX=sigma_XX,
            sigma_XY=sigma_XY,
            sigma_YY=sigma_YY,
            mu_X=mu_X,
            mu_X_conditional=mu_X_c,
            n_X=n_X,
        )
        run_summaries.append(rec)
        print(
            "[test] %s  chi2/ndof pre=%.3f post=%.3f  p=%.4g"
            % (tag, rec["chi2_prefit_per_ndof"], rec["chi2_postfit_per_ndof"], rec["p_value"]),
            flush=True,
        )

    (out_root / "test_summaries.json").write_text(json.dumps(run_summaries, indent=2) + "\n")
    print("[ok] wrote %s" % out_root, flush=True)
    return 0


def _blkdiag_from_prl(
    prl_abs: Mapping[str, Mapping[str, np.ndarray]],
    cat: str,
) -> np.ndarray:
    mats = [prl_abs[n][cat] for n in STACK_NAMES]
    n = int(sum(m.shape[0] for m in mats))
    out = np.zeros((n, n), dtype=np.float64)
    o = 0
    for m in mats:
        k = int(m.shape[0])
        out[o : o + k, o : o + k] = m
        o += k
    return out


if __name__ == "__main__":
    raise SystemExit(main())
