#!/usr/bin/env python3
"""Reduce phase: merge joint GENIE chunk pickles → ``syst_disk_CC/JointGenie/joint_genie_combined.npz``.

For each GENIE knob group with ``nu__joint_cc_genie__<GROUP>__*.pkl`` files:

1. :func:`get_systematics_genie.merge_genie_chunk_pickles` sums rate accumulators across ``.df`` stems.
2. Within the group, per preset *pair_slug*, knob-level covariances are built with
   :func:`pyanalib.covariance.get_covariance_matrix` and combined as independent knobs via
   :func:`syst_multisim_common.combine_indep_knob_cov_packs` (**sum ``cov_frac``**, rebuild
   ``cov``; same recipe as marginal GENIE ``chunk-merge`` for the rate path).
3. Per group, ``cov_frac`` / ``corr`` are built from that group's merged absolute covariance
   and nominal CV (each GENIE bundle may use a different ``.df`` tree in
   ``dataset_locations.GENIE_GROUP_GLOBS``). **Fractional** matrices are **summed** across
   groups (same recipe as :mod:`syst_genie_aggregate` for marginal GENIE). Stored ``cov`` /
   ``corr`` are reconstructed from the summed ``cov_frac`` and the **first** group's CV so
   the NPZ cells stay self-consistent; :func:`cc_joint_cov.build_joint_covariance_abs` still
   rescales ``cov_frac`` with the analysis ``mu`` stack. Per-reweight-knob matrices are stored
   under ``JointGenie_by_knob`` in the output NPZ.

**Note:** This pipeline uses the GENIE **rate** universe histograms (aligned with
``genie_rate`` in the marginal pickle), not the response-matrix **xsec** tensors in the
``genie`` key of ``cov_mat_dict.pkl``.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import defaultdict
from os import path
from typing import Sequence

import numpy as np
from tqdm import tqdm

_REPO_ROOT = path.abspath(path.join(path.dirname(__file__), "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from pyanalib.covariance import cov_from_fraccov, get_covariance_matrix  # noqa: E402

from analysis_village.numucc_1p0pi.dataset_locations import GENIE_GROUP_ORDER, _genie_glob_map  # noqa: E402
from analysis_village.numucc_1p0pi.scripts.get_systematics_genie import RATE_ACC_KEY, merge_genie_chunk_pickles  # noqa: E402
from analysis_village.numucc_1p0pi.syst_cc_joint_multisim_common import (  # noqa: E402
    JOINT_CC_GENIE_CHUNK_GLOB,
    default_kinematic_joint_pairs,
    joint_meta,
    save_joint_genie_combined_npz,
)
from analysis_village.numucc_1p0pi.syst_disk_cc_layout import normalized_root  # noqa: E402
from analysis_village.numucc_1p0pi.syst_multisim_common import combine_indep_knob_cov_packs  # noqa: E402


def collect_joint_genie_chunks(chunks_dir: str) -> list[str]:
    root = path.abspath(path.expanduser(chunks_dir.rstrip(os.sep)))
    return sorted(glob.glob(path.join(root, JOINT_CC_GENIE_CHUNK_GLOB)))


def collect_joint_genie_chunks_many(chunks_dirs: Sequence[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for d in chunks_dirs:
        d = str(d).strip()
        if not d:
            continue
        for p in collect_joint_genie_chunks(d):
            if p not in seen:
                seen.add(p)
                out.append(p)
    return sorted(out)


def _parse_group_from_filename(fp: str) -> str | None:
    base = path.basename(fp)
    if base.endswith(".pkl"):
        if base.startswith("nu__joint_cc_genie__"):
            rest = base[len("nu__joint_cc_genie__") : -len(".pkl")]
        elif base.startswith("nu__joint_genie__"):
            rest = base[len("nu__joint_genie__") : -len(".pkl")]
        else:
            return None
    else:
        return None
    for g in sorted(GENIE_GROUP_ORDER, key=lambda s: -len(s)):
        prefix = g + "__"
        if rest.startswith(prefix):
            return g
    return None


def _ordered_groups_present(paths: list[str], stage: str) -> list[str]:
    gmap = _genie_glob_map(stage)
    have = {_parse_group_from_filename(p) for p in paths}
    have.discard(None)
    out: list[str] = []
    seen: set[str] = set()
    for g in GENIE_GROUP_ORDER:
        if g in gmap and g in have:
            out.append(g)
            seen.add(g)
    for g in sorted(have - seen):
        out.append(g)
    return out


def _triplet_from_gcm(mat: dict) -> dict[str, np.ndarray]:
    return {
        "cov": np.asarray(mat["cov"], dtype=np.float64),
        "cov_frac": np.asarray(mat["cov_frac"], dtype=np.float64),
        "corr": np.asarray(mat["corr"], dtype=np.float64),
    }


def _per_pair_by_knob_one_group(merged: dict) -> dict[str, dict[str, dict[str, np.ndarray]]]:
    """Per pair_slug, per GENIE knob: covariance from merged rate universes (one knob group)."""
    rate = merged.get(RATE_ACC_KEY) or {}
    out: dict[str, dict[str, dict[str, np.ndarray]]] = {}
    for knob, blk in sorted(rate.items()):
        if not isinstance(blk, dict):
            continue
        for pair_slug, pack in blk.items():
            univ = np.asarray(pack["univ"], dtype=float)
            cv = np.asarray(pack["cv"], dtype=float).reshape(-1)
            if univ.shape[0] < 2:
                continue
            out.setdefault(pair_slug, {})[knob] = _triplet_from_gcm(get_covariance_matrix(univ, cv))
    return out


def _per_pair_cov_one_group(merged: dict) -> dict[str, dict[str, np.ndarray]]:
    """Per pair_slug: combine all knobs in one merged GENIE group into absolute cov + CV."""
    rate = merged.get(RATE_ACC_KEY) or {}
    pair_slugs: set[str] = set()
    for _knob, d in rate.items():
        if isinstance(d, dict):
            pair_slugs.update(d.keys())
    out: dict[str, dict[str, np.ndarray]] = {}
    for pair_slug in sorted(pair_slugs):
        packs = []
        cv_ref: np.ndarray | None = None
        for knob in sorted(rate.keys()):
            blk = rate[knob]
            if not isinstance(blk, dict):
                continue
            pack = blk.get(pair_slug)
            if pack is None:
                continue
            univ = np.asarray(pack["univ"], dtype=float)
            cv = np.asarray(pack["cv"], dtype=float).reshape(-1)
            if univ.shape[0] < 2:
                continue
            if cv_ref is None:
                cv_ref = cv.copy()
            elif cv_ref.shape != cv.shape:
                raise ValueError(
                    "[cc-joint-genie-agg] CV shape mismatch for pair %s knob %s: %s vs %s"
                    % (pair_slug, knob, cv_ref.shape, cv.shape)
                )
            packs.append(get_covariance_matrix(univ, cv))
        if not packs or cv_ref is None:
            continue
        comb = combine_indep_knob_cov_packs(packs, cv_ref)
        out[pair_slug] = {"cov": np.asarray(comb["cov"], dtype=float), "cv": cv_ref.copy()}
    return out


def _finalize_pair(pack: dict[str, np.ndarray]) -> dict:
    cov = np.asarray(pack["cov"], dtype=float)
    cv = np.asarray(pack["cv"], dtype=float).reshape(-1)
    safe = np.outer(np.maximum(cv, 1e-18), np.maximum(cv, 1e-18))
    cov_frac = np.divide(cov, safe, out=np.zeros_like(cov), where=safe > 0)
    d = np.sqrt(np.maximum(np.diag(cov), 0.0))
    outer = np.outer(np.maximum(d, 1e-18), np.maximum(d, 1e-18))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = np.where(outer > 0, cov / outer, 0.0)
    np.fill_diagonal(corr, 1.0)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    return {"cov": cov, "cov_frac": cov_frac, "corr": corr}


def run_joint_genie_aggregate(
    chunks_dirs: Sequence[str],
    syst_disk_cc_root: str,
    *,
    mc_df_stage: str = "final",
) -> str:
    roots = [path.abspath(path.expanduser(str(d).strip())) for d in chunks_dirs if str(d).strip()]
    if not roots:
        raise SystemExit("[cc-joint-genie-agg] no chunk directories")
    paths = collect_joint_genie_chunks_many(roots)
    if not paths:
        raise RuntimeError("[cc-joint-genie-agg] no nu__joint_cc_genie__*.pkl under %s" % roots)

    groups = _ordered_groups_present(paths, mc_df_stage)
    if not groups:
        raise RuntimeError("[cc-joint-genie-agg] could not parse GENIE groups from filenames")

    by_group: dict[str, list[str]] = defaultdict(list)
    for p in paths:
        g = _parse_group_from_filename(p)
        if g:
            by_group[g].append(p)

    # Sum fractional covariances across groups (each group may use a different GENIE .df
    # bundle with its own nominal CV — same pattern as syst_genie_aggregate).
    cov_frac_sum_by_pair: dict[str, np.ndarray | None] = {}
    cv_ref_by_pair: dict[str, np.ndarray | None] = {}
    by_knob_by_pair: dict[str, dict[str, dict[str, np.ndarray]]] = defaultdict(dict)

    for grp in tqdm(groups, desc="GENIE groups (joint)"):
        gpaths = sorted(by_group.get(grp, ()))
        if not gpaths:
            continue
        merged = merge_genie_chunk_pickles(gpaths)
        per_knob = _per_pair_by_knob_one_group(merged)
        for pair_slug, kmap in per_knob.items():
            for knob, trip in kmap.items():
                if knob in by_knob_by_pair[pair_slug]:
                    raise ValueError(
                        "[cc-joint-genie-agg] duplicate knob %r for pair %r across GENIE groups "
                        "(unexpected — each knob should belong to a single group)"
                        % (knob, pair_slug)
                    )
                by_knob_by_pair[pair_slug][knob] = trip
        per = _per_pair_cov_one_group(merged)
        for pair_slug, ppack in per.items():
            fin_g = _finalize_pair({"cov": ppack["cov"], "cv": ppack["cv"]})
            cf = fin_g["cov_frac"]
            if cov_frac_sum_by_pair.get(pair_slug) is None:
                cov_frac_sum_by_pair[pair_slug] = cf.copy()
                cv_ref_by_pair[pair_slug] = np.asarray(ppack["cv"], dtype=np.float64).reshape(-1).copy()
            else:
                cov_frac_sum_by_pair[pair_slug] = cov_frac_sum_by_pair[pair_slug] + cf

    per_out: dict[str, dict] = {}
    for pair_slug, cfsum in cov_frac_sum_by_pair.items():
        if cfsum is None or cv_ref_by_pair.get(pair_slug) is None:
            continue
        cv_ref = cv_ref_by_pair[pair_slug]
        cov_stored = cov_from_fraccov(cfsum, cv_ref)
        fin = _finalize_pair({"cov": cov_stored, "cv": cv_ref})
        nx = ny = None
        ntot = int(cv_ref.size)
        var_x_name = var_y_name = ""
        for slug, vx, vy in default_kinematic_joint_pairs():
            if slug == pair_slug:
                nx = len(vx.bin_centers)
                ny = len(vy.bin_centers)
                var_x_name = vx.var_save_name
                var_y_name = vy.var_save_name
                if nx + ny != ntot:
                    raise ValueError(
                        "pair %s: n_X+n_Y=%d+%d != cv len %d" % (pair_slug, nx, ny, ntot)
                    )
                break
        if nx is None:
            raise ValueError("[cc-joint-genie-agg] unknown pair_slug %r" % pair_slug)
        per_out[pair_slug] = {
            **fin,
            "meta": {**joint_meta(nx, ny, var_x_name, var_y_name), "source": "GENIE_rate_joint", "n_genie_groups": len(groups)},
            "by_knob": dict(by_knob_by_pair.get(pair_slug, {})),
        }

    if not per_out:
        raise RuntimeError("[cc-joint-genie-agg] built zero pair covariances")
    out_path = save_joint_genie_combined_npz(per_out, syst_disk_cc_root)
    print("[cc-joint-genie-agg] wrote", out_path)
    summ = path.join(normalized_root(syst_disk_cc_root), "JointGenie", "joint_genie_aggregate_summary.json")
    os.makedirs(path.dirname(summ), exist_ok=True)
    with open(summ, "w") as f:
        json.dump(
            {
                "chunks": paths,
                "groups": groups,
                "pair_slugs": sorted(per_out.keys()),
                "genie_knobs_by_pair": {slug: sorted((per_out[slug].get("by_knob") or {}).keys()) for slug in per_out},
            },
            f,
            indent=2,
        )
    print("[cc-joint-genie-agg] wrote", summ)
    return out_path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--chunks_dir",
        dest="chunks_dirs",
        nargs="+",
        required=True,
        metavar="DIR",
        help="Directories containing nu__joint_cc_genie__*.pkl",
    )
    p.add_argument("--syst-disk-cc-root", dest="syst_disk_cc_root", required=True)
    p.add_argument("--mc-df-stage", choices=("final", "sel_all"), default="final")
    return p.parse_args()


def main():
    args = parse_args()
    run_joint_genie_aggregate(args.chunks_dirs, args.syst_disk_cc_root, mc_df_stage=args.mc_df_stage)


if __name__ == "__main__":
    main()
