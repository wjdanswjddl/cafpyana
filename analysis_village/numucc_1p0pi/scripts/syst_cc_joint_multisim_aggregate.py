#!/usr/bin/env python3
"""Reduce phase: merge joint multisim chunk pickles → per-category ``syst_disk_CC`` trees.

Merges ``nu__joint_cc__*.pkl`` from :mod:`syst_cc_joint_multisim_chunk`, builds per-systematic
covariance matrices with :func:`pyanalib.covariance.get_covariance_matrix`, combines Flux/G4
knobs independently via :func:`syst_multisim_common.combine_indep_knob_cov_packs`, then writes
**one NPZ per category** under ``JointMCstat/``, ``JointFlux/``, and ``JointG4/`` (parallel
directories). :mod:`cc_joint_cov` sums fractional covariances across those files when
building the total joint multisim block (legacy ``JointMultisim/joint_multisim_combined.npz``
is still supported if present).

Legacy map shards named ``nu__joint__*.pkl`` are **not** merged (pair coverage may be stale);
re-run the map phase to produce ``nu__joint_cc__*`` shards.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import sys
from collections.abc import Sequence
from os import path

import numpy as np
from tqdm import tqdm

sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))

from pyanalib.covariance import cov_from_fraccov, get_covariance_matrix  # noqa: E402

from analysis_village.numucc_1p0pi.syst_cc_joint_multisim_common import (  # noqa: E402
    JOINT_CC_MULTISIM_CHUNK_GLOB,
    default_kinematic_joint_pairs,
    joint_meta,
    save_joint_multisim_category_npz,
)
from analysis_village.numucc_1p0pi.syst_disk_cc_layout import normalized_root  # noqa: E402
from analysis_village.numucc_1p0pi.syst_multisim_common import (  # noqa: E402
    combine_indep_knob_cov_packs,
    knob_nested_syst_block,
    parse_neutrino_syst_type_csv,
    syst_acc_bucket_nonempty,
)


def collect_joint_chunks(chunks_dir: str) -> list[str]:
    root = path.abspath(path.expanduser(chunks_dir.rstrip(os.sep)))
    paths = sorted(glob.glob(path.join(root, JOINT_CC_MULTISIM_CHUNK_GLOB)))
    for sub in ("Combined", "MCstat", "Flux", "G4"):
        d = path.join(root, sub)
        if path.isdir(d):
            for p in sorted(glob.glob(path.join(d, JOINT_CC_MULTISIM_CHUNK_GLOB))):
                if p not in paths:
                    paths.append(p)
    return sorted(set(paths))


def collect_joint_chunks_many(chunks_dirs: Sequence[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for d in chunks_dirs:
        d = str(d).strip()
        if not d:
            continue
        for p in collect_joint_chunks(d):
            if p not in seen:
                seen.add(p)
                out.append(p)
    return sorted(out)


def _merge_joint_pack(dst: dict, pair_slug: str, univ: np.ndarray, cv: np.ndarray, err: str) -> None:
    u = np.asarray(univ, dtype=float)
    c = np.asarray(cv, dtype=float)
    if pair_slug not in dst:
        dst[pair_slug] = {"univ_events": u.copy(), "cv_events": c.copy()}
    else:
        cur = dst[pair_slug]
        if cur["univ_events"].shape != u.shape:
            raise ValueError("joint merge shape mismatch %s: %s vs %s" % (err, cur["univ_events"].shape, u.shape))
        cur["univ_events"] += u
        cur["cv_events"] += c


def _merge_flux_or_g4_joint_block(merged_block: dict, raw_block: dict, label: str) -> None:
    if not raw_block:
        return
    if knob_nested_syst_block(raw_block):
        if merged_block and not knob_nested_syst_block(merged_block):
            raise ValueError("[cc-joint-agg] cannot mix flat %s with knob-nested chunks" % label)
        for knob, kb in raw_block.items():
            tgt = merged_block.setdefault(knob, {})
            for pair_slug, pack in kb.items():
                _merge_joint_pack(tgt, pair_slug, pack["univ_events"], pack["cv_events"], "%s/%s/%s" % (label, knob, pair_slug))
        return
    if merged_block and knob_nested_syst_block(merged_block):
        raise ValueError("[cc-joint-agg] cannot mix knob-nested %s with flat chunks" % label)
    for pair_slug, pack in raw_block.items():
        _merge_joint_pack(merged_block, pair_slug, pack["univ_events"], pack["cv_events"], "%s/%s" % (label, pair_slug))


def merge_joint_chunks(paths: list[str], syst_types: tuple[str, ...]) -> dict:
    merged = None
    for fp in tqdm(paths, desc="merge joint nu chunks"):
        with open(fp, "rb") as f:
            d = pickle.load(f)
        if merged is None:
            merged = {"syst": {sn: {} for sn in syst_types}, "meta": []}
        merged["meta"].append(
            {
                "df_file": d.get("df_file"),
                "pairs": d.get("pairs"),
                "syst_names_computed": d.get("syst_names_computed"),
            }
        )
        raw_syst = d.get("syst") or {}
        for sn in syst_types:
            block = raw_syst.get(sn, {})
            if sn in ("G4", "Flux") and knob_nested_syst_block(block):
                _merge_flux_or_g4_joint_block(merged["syst"][sn], block, sn)
                continue
            for pair_slug, pack in block.items():
                _merge_joint_pack(
                    merged["syst"][sn],
                    pair_slug,
                    pack["univ_events"],
                    pack["cv_events"],
                    "%s/%s" % (sn, pair_slug),
                )
    if merged is None:
        raise RuntimeError("no joint chunks merged")
    return merged


def _triplet(ret: dict) -> dict[str, np.ndarray]:
    return {
        "cov": np.asarray(ret["cov"], dtype=np.float64),
        "cov_frac": np.asarray(ret["cov_frac"], dtype=np.float64),
        "corr": np.asarray(ret["corr"], dtype=np.float64),
    }


def _cov_for_flat_joint_block(cat_block: dict, pair_slug: str) -> dict | None:
    pack = cat_block.get(pair_slug)
    if pack is None:
        return None
    univ = np.asarray(pack["univ_events"], dtype=float)
    cv = np.asarray(pack["cv_events"], dtype=float)
    if univ.size == 0 or cv.size == 0:
        return None
    return get_covariance_matrix(univ, cv)


def _pair_layout_meta(pair_slug: str) -> tuple[int, int, str, str]:
    for slug, vx, vy in default_kinematic_joint_pairs():
        if slug == pair_slug:
            return len(vx.bin_centers), len(vy.bin_centers), vx.var_save_name, vy.var_save_name
    raise ValueError("[cc-joint-agg] unknown pair_slug %r" % pair_slug)


def joint_covariance_by_category_from_merged(
    merged: dict, syst_types: tuple[str, ...]
) -> dict[str, dict[str, dict]]:
    """``category -> pair_slug -> pack`` with single-category ``cov`` / ``cov_frac`` / ``corr`` / ``meta``."""
    pair_slugs: set[str] = set()
    for sn in syst_types:
        blk = merged["syst"].get(sn) or {}
        if knob_nested_syst_block(blk):
            for kb in blk.values():
                pair_slugs.update(kb.keys())
        else:
            pair_slugs.update(blk.keys())

    out: dict[str, dict[str, dict]] = {sn: {} for sn in syst_types}
    for pair_slug in sorted(pair_slugs):
        nx, ny, var_x_name, var_y_name = _pair_layout_meta(pair_slug)
        ntot = nx + ny
        for sn in syst_types:
            blk = merged["syst"].get(sn) or {}
            if not syst_acc_bucket_nonempty(sn, blk):
                continue
            ret = None
            by_knob_payload: dict[str, dict[str, np.ndarray]] | None = None
            cv_flat: np.ndarray | None = None
            if sn in ("Flux", "G4") and knob_nested_syst_block(blk):
                knob_rets = []
                knob_triplets: dict[str, dict[str, np.ndarray]] = {}
                cv_ref_sn = None
                for knob in sorted(blk.keys()):
                    pack = blk[knob].get(pair_slug)
                    if pack is None:
                        continue
                    univ = np.asarray(pack["univ_events"], dtype=float)
                    cv = np.asarray(pack["cv_events"], dtype=float).reshape(-1)
                    if univ.shape[0] < 2:
                        continue
                    if cv_ref_sn is None:
                        cv_ref_sn = cv.copy()
                    elif cv_ref_sn.shape != cv.shape:
                        raise ValueError(
                            "[cc-joint-agg] CV shape mismatch %s pair %s knob %s: %s vs %s"
                            % (sn, pair_slug, knob, cv_ref_sn.shape, cv.shape)
                        )
                    one = get_covariance_matrix(univ, cv)
                    knob_rets.append(one)
                    knob_triplets[knob] = _triplet(one)
                if not knob_rets or cv_ref_sn is None:
                    continue
                ret = combine_indep_knob_cov_packs(knob_rets, cv_ref_sn)
                cv_flat = np.asarray(cv_ref_sn, dtype=np.float64).reshape(-1)
                by_knob_payload = knob_triplets
            else:
                ret = _cov_for_flat_joint_block(blk, pair_slug)
                if ret is None:
                    continue
                pack = blk.get(pair_slug)
                if pack is None:
                    continue
                cv_flat = np.asarray(pack["cv_events"], dtype=np.float64).reshape(-1)
            if cv_flat.size != ntot:
                raise ValueError(
                    "[cc-joint-agg] %s pair %s: cv len %d != n_X+n_Y=%d+%d"
                    % (sn, pair_slug, int(cv_flat.size), nx, ny)
                )
            cov_sum = cov_from_fraccov(np.asarray(ret["cov_frac"], dtype=np.float64), cv_flat)
            d = np.sqrt(np.maximum(np.diag(cov_sum), 0.0))
            outer = np.outer(np.maximum(d, 1e-18), np.maximum(d, 1e-18))
            with np.errstate(divide="ignore", invalid="ignore"):
                corr = np.where(outer > 0, cov_sum / outer, 0.0)
            np.fill_diagonal(corr, 1.0)
            corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
            pack_out: dict = {
                "cov": cov_sum,
                "cov_frac": np.asarray(ret["cov_frac"], dtype=np.float64),
                "corr": corr,
                "meta": {
                    **joint_meta(nx, ny, var_x_name, var_y_name),
                    "joint_multisim_category": sn,
                },
            }
            if by_knob_payload:
                pack_out["by_knob"] = by_knob_payload
            out[sn][pair_slug] = pack_out
    return out


def run_joint_aggregate(
    chunks_dirs: Sequence[str],
    syst_disk_cc_root: str,
    syst_types: tuple[str, ...] | None = None,
) -> None:
    roots = [path.abspath(path.expanduser(str(d).strip())) for d in chunks_dirs if str(d).strip()]
    if not roots:
        raise SystemExit("[cc-joint-agg] no chunk directories")
    st = syst_types if syst_types is not None else parse_neutrino_syst_type_csv(None)
    ck = collect_joint_chunks_many(roots)
    if not ck:
        raise RuntimeError("[cc-joint-agg] no nu__joint_cc__*.pkl under %s" % roots)
    print("[cc-joint-agg] merging %d chunk(s) syst_types=%s" % (len(ck), ",".join(st)))
    merged = merge_joint_chunks(ck, st)
    empty = [sn for sn in st if not syst_acc_bucket_nonempty(sn, merged["syst"].get(sn))]
    if empty:
        raise RuntimeError("[cc-joint-agg] empty merged data for: %s" % (empty,))
    by_cat = joint_covariance_by_category_from_merged(merged, st)
    written: list[str] = []
    all_pairs: set[str] = set()
    for sn in st:
        per_pair = by_cat.get(sn) or {}
        if not per_pair:
            continue
        all_pairs.update(per_pair.keys())
        outp = save_joint_multisim_category_npz(per_pair, syst_disk_cc_root, sn)
        written.append(outp)
        print("[cc-joint-agg] wrote", outp)
    if not written:
        raise RuntimeError("[cc-joint-agg] built zero category covariances")
    summ = path.join(normalized_root(syst_disk_cc_root), "joint_multisim_aggregate_summary.json")
    os.makedirs(path.dirname(summ), exist_ok=True)
    with open(summ, "w") as f:
        json.dump(
            {
                "chunks": merged["meta"],
                "pair_slugs": sorted(all_pairs),
                "syst_types": list(st),
                "npz_outputs": written,
            },
            f,
            indent=2,
        )
    print("[cc-joint-agg] wrote", summ)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--chunks_dir",
        dest="chunks_dirs",
        nargs="+",
        required=True,
        metavar="DIR",
        help="One or more directories containing nu__joint_cc__*.pkl (flat or Combined/…).",
    )
    p.add_argument(
        "--syst-disk-cc-root",
        dest="syst_disk_cc_root",
        required=True,
        help="Output root for syst_disk_CC (JointFlux/, JointG4/, JointMCstat/).",
    )
    p.add_argument(
        "--syst-types",
        default=None,
        metavar="CSV",
        help="Subset of MCstat,Flux,G4 (default: all three).",
    )
    return p.parse_args()


def main():
    args = parse_args()
    st = parse_neutrino_syst_type_csv(args.syst_types)
    run_joint_aggregate(args.chunks_dirs, args.syst_disk_cc_root, syst_types=st)


if __name__ == "__main__":
    main()
