"""Stack per-variable ``rate_univ_cv`` histograms into a 49-bin joint.

May 14 / Sep 12 GENIE chunk trees are **read-only**. Universe index *u* is the
same event-weight universe for every variable in a chunk, so concatenating
``[muon-p, muon-dir_z, proton-p, proton-dir_z]`` at fixed *u* is a valid
inclusive joint (same recipe as ``constraint_joint_prl_check._sum_rate_acc``).

Covariance matches :func:`pyanalib.covariance.get_covariance_matrix` (mean of
outer products, including *n_univ* = 1 unisims).
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from glob import glob
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Sequence, Tuple

import numpy as np

from analysis_village.numucc_1p0pi.scripts._genie_pkl_compat import _NumpyCompatUnpickler
from analysis_village.numucc_1p0pi.scripts.constraint_joint_prl_check import (
    MAY_AR23P_CHUNKS,
    MAY_INTERP_KNOBS,
    MAY_MEC_CHUNKS,
    MAY_SBNV1_MEC_KNOBS,
    STACK_NAMES,
)
from pyanalib.covariance import corr_from_fraccov, cov_from_fraccov

SEP_AR23P_CHUNKS = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "genie_syst-chunked-sel_mup_20260912_Ar23p/chunks/genie__Ar23p__*.pkl"
)
SEP_INTERP_KNOBS = tuple(
    f"MECq0q3InterpWeighting_SBN_v3_SuSATo{_to}_MECResponse_q0bin{_b}"
    for _to in ("Val", "Mar")
    for _b in range(4)
)
# Sep interpolator name → May interpolator name (Product B MEC map).
INTERP_SEP_TO_MAY: Dict[str, str] = dict(zip(SEP_INTERP_KNOBS, MAY_INTERP_KNOBS))

_FRAC_CV_EPS = 1e-12


def _load_pickle(path: str) -> dict:
    with open(path, "rb") as f:
        return _NumpyCompatUnpickler(f).load()


def cov_from_univ(univ: np.ndarray, cv: np.ndarray) -> Dict[str, np.ndarray]:
    """Vectorized :func:`pyanalib.covariance.get_covariance_matrix`."""
    u = np.asarray(univ, dtype=np.float64)
    c = np.asarray(cv, dtype=np.float64).reshape(-1)
    if u.ndim != 2 or u.shape[1] != c.size:
        raise ValueError("univ shape %s vs cv %s" % (u.shape, c.shape))
    n_u = max(int(u.shape[0]), 1)
    d = u - c
    cov = (d.T @ d) / n_u
    scale = np.maximum(c, _FRAC_CV_EPS)
    df = d / scale
    frac = (df.T @ df) / n_u
    return {
        "cov": cov,
        "cov_frac": frac,
        "corr": corr_from_fraccov(frac),
        "n_univ": int(u.shape[0]),
        "cv": c,
        "univ": u,
    }


def _extract_knob_vars(
    path: str, knobs: Sequence[str], var_names: Sequence[str]
) -> Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]]:
    blob = _load_pickle(path)
    rate = blob.get("rate_univ_cv") or {}
    out: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]] = {}
    for kn in knobs:
        slot = rate.get(kn)
        if not isinstance(slot, dict):
            continue
        for vn in var_names:
            cell = slot.get(vn)
            if not isinstance(cell, dict) or "univ" not in cell:
                continue
            u = np.asarray(cell["univ"], dtype=np.float64)
            c = np.asarray(cell["cv"], dtype=np.float64).reshape(-1)
            out[(kn, vn)] = (u, c)
    return out


def _extract_knob_vars_star(
    args: Tuple[str, Sequence[str], Sequence[str]],
) -> Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]]:
    path, knobs, var_names = args
    return _extract_knob_vars(path, knobs, var_names)


def _add_into(
    acc: Dict[str, Dict[str, Dict[str, np.ndarray]]],
    extracted: Mapping[Tuple[str, str], Tuple[np.ndarray, np.ndarray]],
) -> None:
    for (kn, vn), (u, c) in extracted.items():
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


def sum_rate_univ_cv(
    paths: Sequence[str],
    knobs: Sequence[str],
    var_names: Sequence[str] = STACK_NAMES,
    *,
    workers: int = 1,
) -> Dict[str, Dict[str, Any]]:
    """Sum chunk ``rate_univ_cv`` and stack ``var_names`` at the same universe index."""
    acc: Dict[str, Dict[str, Any]] = {k: {v: None for v in var_names} for k in knobs}
    n_files = 0
    path_list = [str(p) for p in paths]
    knobs_t = tuple(knobs)
    vars_t = tuple(var_names)
    jobs = [(p, knobs_t, vars_t) for p in path_list]
    n_jobs = len(jobs)
    pool: ProcessPoolExecutor | None = None
    if workers > 1:
        try:
            pool = ProcessPoolExecutor(max_workers=int(workers))
            extracted_iter = pool.map(_extract_knob_vars_star, jobs, chunksize=8)
        except Exception:
            pool = None
            extracted_iter = (_extract_knob_vars_star(j) for j in jobs)
    else:
        extracted_iter = (_extract_knob_vars_star(j) for j in jobs)
    try:
        for extracted in extracted_iter:
            _add_into(acc, extracted)
            n_files += 1
            if n_jobs >= 200 and n_files % 200 == 0:
                print("[univ-stack] %d / %d chunks" % (n_files, n_jobs), flush=True)
    finally:
        if pool is not None:
            pool.shutdown(wait=True)
    out: Dict[str, Dict[str, Any]] = {}
    for kn in knobs:
        missing = [v for v in var_names if acc[kn][v] is None]
        if missing:
            continue
        u_j = np.hstack([acc[kn][v]["univ"] for v in var_names])
        c_j = np.concatenate([acc[kn][v]["cv"] for v in var_names])
        pack = cov_from_univ(u_j, c_j)
        pack["n_files"] = n_files
        out[kn] = pack
    return out


def sum_frac(by_knob: Mapping[str, Mapping[str, np.ndarray]], knobs: Iterable[str]) -> np.ndarray:
    mats = [np.asarray(by_knob[k]["cov_frac"], dtype=np.float64) for k in knobs if k in by_knob]
    if not mats:
        raise ValueError("no knobs to sum: %s" % list(knobs))
    return np.sum(mats, axis=0)


def overlay_pack(
    frac: np.ndarray, mu_nu: np.ndarray
) -> Dict[str, np.ndarray]:
    cov = cov_from_fraccov(frac, mu_nu)
    return {"cov_frac": np.asarray(frac, dtype=float), "cov": cov, "corr": corr_from_fraccov(frac)}


def psd_clip(mat: np.ndarray) -> Tuple[np.ndarray, float]:
    s = 0.5 * (np.asarray(mat, dtype=np.float64) + np.asarray(mat, dtype=np.float64).T)
    w, v = np.linalg.eigh(s)
    min_eig = float(w.min()) if w.size else 0.0
    if min_eig < 0:
        w = np.clip(w, 0.0, None)
        s = (v * w) @ v.T
        s = 0.5 * (s + s.T)
    return s, min_eig


def load_may_mec_joints(*, workers: int = 1) -> Dict[str, Dict[str, Any]]:
    ar23p = sorted(glob(MAY_AR23P_CHUNKS))
    mec = sorted(glob(MAY_MEC_CHUNKS))
    if not ar23p or not mec:
        raise RuntimeError("May 14 GENIE chunks missing (Ar23p=%d MEC=%d)" % (len(ar23p), len(mec)))
    out = {}
    out.update(sum_rate_univ_cv(ar23p, MAY_INTERP_KNOBS, STACK_NAMES, workers=workers))
    out.update(sum_rate_univ_cv(mec, MAY_SBNV1_MEC_KNOBS, STACK_NAMES, workers=workers))
    return out


def load_sep_interp_joints(*, workers: int = 8) -> Dict[str, Dict[str, Any]]:
    paths = sorted(glob(SEP_AR23P_CHUNKS))
    if not paths:
        raise RuntimeError("Sep 12 Ar23p GENIE chunks missing: %s" % SEP_AR23P_CHUNKS)
    return sum_rate_univ_cv(paths, SEP_INTERP_KNOBS, STACK_NAMES, workers=workers)


def product_b_mec_frac_delta(
    may_by_knob: Mapping[str, Mapping[str, Any]],
    sep_by_knob: Mapping[str, Mapping[str, Any]],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(delta_frac, may_interp_frac, sep_interp_frac)`` for the 8 interpolators.

    Constraint variables only swap interpolators (SBN_v1 MEC has no Sep per-knob
    matrix on muon/proton kinematics).
    """
    missing_may = [k for k in MAY_INTERP_KNOBS if k not in may_by_knob]
    missing_sep = [k for k in SEP_INTERP_KNOBS if k not in sep_by_knob]
    if missing_may or missing_sep:
        raise RuntimeError("missing interpolator knobs may=%s sep=%s" % (missing_may, missing_sep))
    may_frac = sum_frac(may_by_knob, MAY_INTERP_KNOBS)
    sep_frac = sum_frac(sep_by_knob, SEP_INTERP_KNOBS)
    return may_frac - sep_frac, may_frac, sep_frac
