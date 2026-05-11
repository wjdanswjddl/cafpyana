#!/usr/bin/env python
"""Reduce phase for chunked Flux / G4 / MCstat covariances (neutrino multisim only).

* Input: ``nu__*.pkl`` from ``syst_multisim_chunk.py``. Pass one or more ``--chunks_dir`` paths
  (each is a **chunks** directory): under every root, pickles are collected from ``MCstat/``,
  ``Flux/``, ``G4/``, ``Combined/`` if present, and from the root itself (legacy flat layout).
  ``run_syst_multisim_chunked.sh`` passes ``multisim_syst-chunked-*/chunks``,
  ``g4_syst-chunked-*/chunks``, and ``flux_syst-chunked-*/chunks`` together.
* Sum ``univ_events`` / ``cv_events`` across shards, build per-systematic covariances,
  write ``MCstat/``, ``Flux/``, ``G4/`` under ``--syst-disk-root``.

Cosmic background uncertainties are **not** produced here; use
``run_syst_cosmics_chunked.sh`` → ``syst_cosmics_aggregate.py`` (writes ``Cosmics/``).
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

import matplotlib

matplotlib.use("Agg")
import numpy as np
from tqdm import tqdm

sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))

from pyanalib.covariance import get_covariance_matrix
from analysis_village.numucc_1p0pi.syst_disk_layout import category_out_dir, normalized_root
from analysis_village.numucc_1p0pi.syst_multisim_common import (
    NEUTRINO_SYST_ORDER,
    build_var_configs,
    combine_indep_knob_cov_packs,
    count_merged_knob_nested_var_slots,
    knob_nested_syst_block,
    save_neutrino_multisim_npzs,
)
from analysis_village.numucc_1p0pi.utils import (
    plot_heatmap,
    plot_univ_hists,
)

# turn off performance warnings
import warnings
import pandas as pd
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


def collect_nu_chunks(chunks_dir: str) -> list[str]:
    """Paths to ``nu__*.pkl`` under typed subdirs and/or legacy flat ``chunks_dir``."""
    root = path.abspath(path.expanduser(chunks_dir.rstrip(os.sep)))
    paths: list[str] = []
    seen: set[str] = set()
    for sub in ("MCstat", "Flux", "G4", "Combined"):
        d = path.join(root, sub)
        if path.isdir(d):
            for p in glob.glob(path.join(d, "nu__*.pkl")):
                if p not in seen:
                    seen.add(p)
                    paths.append(p)
    for p in glob.glob(path.join(root, "nu__*.pkl")):
        if p not in seen:
            seen.add(p)
            paths.append(p)
    return sorted(paths)


def collect_nu_chunks_many(chunks_dirs: Sequence[str]) -> list[str]:
    """Union of :func:`collect_nu_chunks` over several roots (deduped, sorted)."""
    seen: set[str] = set()
    out: list[str] = []
    for d in chunks_dirs:
        d = d.strip()
        if not d:
            continue
        for p in collect_nu_chunks(d):
            if p not in seen:
                seen.add(p)
                out.append(p)
    return sorted(out)


def _detect_input_stage(paths: list[str]) -> str:
    """Peek at the first pickle to decide whether chunks are final-stage or sel_all."""
    if not paths:
        raise RuntimeError("[multisim-agg] no chunks to inspect")
    with open(paths[0], "rb") as f:
        d = pickle.load(f)
    stage = d.get("input_stage")
    if stage in ("final", "sel_all"):
        return stage
    # Legacy chunks: no input_stage field → treat as final.
    return "final"


def _write_covariance_manifest(
    out_dir: str,
    syst_dict: dict,
    merged_meta: list,
    mc_df_stage: str,
    var_set: str,
) -> None:
    """Small JSON sidecar listing variables and NPZ outputs for downstream tools."""
    neutrino_npz = [
        "mcstat_syst_dict.npz",
        "flux_syst_dict.npz",
        "g4_syst_dict.npz",
    ]
    var_names: set[str] = set()
    for _cat, block in (syst_dict or {}).items():
        if isinstance(block, dict):
            var_names.update(block.keys())
    manifest = {
        "schema": "numucc_multisim_covariance_v1",
        "description": "Fractional and absolute covariance packs per category; NPZs use numpy.savez_compressed",
        "mc_df_stage": mc_df_stage,
        "var_set": var_set,
        "categories_present": sorted([k for k, v in (syst_dict or {}).items() if v]),
        "variables": sorted(var_names),
        "neutrino_multisim_npz": [
            {"file": n, "role": "per-variable dict → inner syst matrices"} for n in neutrino_npz
        ],
        "map_shard_metadata": merged_meta,
    }
    outp = path.join(out_dir, "covariance_manifest.json")
    with open(outp, "w") as f:
        json.dump(manifest, f, indent=2)
    print("[multisim-agg] wrote", outp)


def _merge_pack_into(dst: dict, vsn: str, pack: dict, sn: str, err_label: str) -> None:
    u = np.asarray(pack["univ_events"], dtype=float)
    c = np.asarray(pack["cv_events"], dtype=float)
    stage_key = pack.get("stage_key", "")
    if vsn not in dst:
        dst[vsn] = {
            "univ_events": u.copy(),
            "cv_events": c.copy(),
            "stage_key": stage_key,
        }
    else:
        mu = dst[vsn]
        if mu["univ_events"].shape != u.shape:
            raise ValueError(
                "shape mismatch %s %s: %s vs %s" % (sn, err_label, mu["univ_events"].shape, u.shape)
            )
        mu["univ_events"] += u
        mu["cv_events"] += c


def _merge_flux_or_g4_block(merged_block: dict, raw_block: dict, label: str) -> None:
    """Merge one chunk's Flux or G4 block (flat ``{var: pack}`` or nested ``{knob: {var: pack}}``)."""
    if not raw_block:
        return
    if knob_nested_syst_block(raw_block):
        if merged_block and not knob_nested_syst_block(merged_block):
            raise ValueError(
                "[multisim-agg] cannot mix flat %s chunks with knob-nested chunks in one aggregate"
                % label
            )
        for knob, kb in raw_block.items():
            tgt = merged_block.setdefault(knob, {})
            for vsn, pack in kb.items():
                _merge_pack_into(tgt, vsn, pack, label, "%s/%s" % (knob, vsn))
        return
    if merged_block and knob_nested_syst_block(merged_block):
        raise ValueError(
            "[multisim-agg] cannot mix knob-nested %s chunks with flat bundled chunks in one aggregate"
            % label
        )
    for vsn, pack in raw_block.items():
        _merge_pack_into(merged_block, vsn, pack, label, vsn)


def merge_nu_chunks(paths: list[str]) -> dict:
    merged = None
    for fp in tqdm(paths, desc="merge nu chunks"):
        with open(fp, "rb") as f:
            d = pickle.load(f)
        if merged is None:
            merged = {"syst": {sn: {} for sn in NEUTRINO_SYST_ORDER}, "meta": []}
        merged["meta"].append(
            {
                "df_file": d.get("df_file"),
                "splits": d.get("splits_processed"),
                "syst_names_computed": d.get("syst_names_computed"),
                "input_stage": d.get("input_stage", "final"),
                "g4_mode": d.get("g4_mode"),
                "flux_mode": d.get("flux_mode"),
                "flux_knob_groups": d.get("flux_knob_groups"),
            }
        )
        raw_syst = d.get("syst") or {}
        for sn in NEUTRINO_SYST_ORDER:
            block = raw_syst.get(sn, {})
            if sn in ("G4", "Flux"):
                _merge_flux_or_g4_block(merged["syst"][sn], block, sn)
                continue
            for vsn, pack in block.items():
                _merge_pack_into(merged["syst"][sn], vsn, pack, sn, vsn)
    if merged is None:
        raise RuntimeError("no chunks merged")
    n_tot = 0
    for sn in NEUTRINO_SYST_ORDER:
        b = merged["syst"][sn]
        if sn in ("G4", "Flux") and knob_nested_syst_block(b):
            n_tot += count_merged_knob_nested_var_slots(b)
        else:
            n_tot += len(b)
    if n_tot == 0:
        hint = (
            "Every chunk's ``syst`` block was empty (no ``univ_events`` keys). "
            "Typical causes: (1) stale ``nu__*.pkl`` from an older failed run — remove them under "
            "``chunks`` and re-run the map phase; (2) ``--n-universe`` larger than stored "
            "``univ_*`` columns (fixed in syst_multisim_chunk for new chunks); (3) all variables "
            "skipped in the chunk script (see map logs). First chunk file: %s"
        ) % (paths[0],)
        raise RuntimeError("[multisim-agg] merged zero systematic variables. " + hint)
    return merged


def _syst_plot_key(sn: str):
    return ("mc", sn) if sn in ("Flux", "G4") else sn


def _covariance_nested_knob_block(
    merged_block: dict,
    syst_dict_sn: dict,
    sn: str,
    vc_by: dict,
    cat_dir: str,
    save_fig: bool,
    sk,
    tag: str,
) -> None:
    """Fill ``syst_dict_sn`` for Flux or G4 when ``merged_block`` is ``{knob: {var: pack}}``."""
    var_names = set()
    for kb in merged_block.values():
        var_names.update(kb.keys())
    for vsn in tqdm(sorted(var_names), desc="cov %s" % sn):
        vc = vc_by.get(vsn)
        if vc is None:
            continue
        knob_rets = []
        cv_ref = None
        univ_first = None
        cv_first = None
        for knob in sorted(merged_block.keys()):
            pack = merged_block[knob].get(vsn)
            if pack is None:
                continue
            univ = np.asarray(pack["univ_events"], dtype=float)
            cv = np.asarray(pack["cv_events"], dtype=float)
            if cv_ref is None:
                cv_ref = cv.copy()
            if univ_first is None:
                univ_first, cv_first = univ, cv
            knob_rets.append(get_covariance_matrix(univ, cv))
        if not knob_rets or cv_ref is None:
            continue
        ret = combine_indep_knob_cov_packs(knob_rets, cv_ref)
        syst_dict_sn[vsn] = ret
        if save_fig and univ_first is not None and cv_first is not None:
            plot_univ_hists(
                univ_first,
                cv_first,
                sk,
                vc,
                plot=False,
                save_fig=True,
                save_name=path.join(cat_dir, "{}-{}-universes".format(vsn, tag)),
            )
            for matrix_type in ["cov", "cov_frac", "corr"]:
                plot_heatmap(
                    ret[matrix_type],
                    vc.bins,
                    plot_labels=[
                        vc.var_labels[1],
                        vc.var_labels[1],
                        matrix_type.capitalize(),
                    ],
                    plot=False,
                    save_fig=True,
                    save_name=path.join(
                        cat_dir,
                        "{}-{}-{}".format(vsn, tag, matrix_type),
                    ),
                )


def covariance_dict_from_merged(
    merged: dict,
    var_configs: list,
    syst_disk_root: str,
    save_fig: bool,
) -> dict:
    """Build per-category covariance dict; plots go under ``<root>/MCstat|Flux|G4/``.

    ``syst_disk_root`` is normalized to an absolute path so figures and NPZs are not
    written relative to the process working directory.
    """
    root = normalized_root(syst_disk_root)
    syst_dict = {sn: {} for sn in NEUTRINO_SYST_ORDER}
    vc_by = {v.var_save_name: v for v in var_configs}
    for sn in NEUTRINO_SYST_ORDER:
        cat_dir = category_out_dir(root, sn)
        os.makedirs(cat_dir, exist_ok=True)
        sk = _syst_plot_key(sn)
        tag = sn
        cat_block = merged["syst"][sn]
        if sn in ("G4", "Flux") and knob_nested_syst_block(cat_block):
            _covariance_nested_knob_block(
                cat_block, syst_dict[sn], sn, vc_by, cat_dir, save_fig, sk, tag
            )
            continue

        for vsn, pack in tqdm(list(merged["syst"][sn].items()), desc="cov %s" % sn):
            vc = vc_by.get(vsn)
            if vc is None:
                continue
            univ = np.asarray(pack["univ_events"], dtype=float)
            cv = np.asarray(pack["cv_events"], dtype=float)
            ret = get_covariance_matrix(univ, cv)
            syst_dict[sn][vsn] = ret
            if save_fig:
                # ``plot=False``: batch/Agg backend; ``save_name`` has no extension here
                # because ``plot_univ_hists`` / ``plot_heatmap`` append ``fig_ext`` (.png).
                plot_univ_hists(
                    univ,
                    cv,
                    sk,
                    vc,
                    plot=False,
                    save_fig=True,
                    save_name=path.join(cat_dir, "{}-{}-universes".format(vsn, tag)),
                )
                for matrix_type in ["cov", "cov_frac", "corr"]:
                    plot_heatmap(
                        ret[matrix_type],
                        vc.bins,
                        plot_labels=[
                            vc.var_labels[1],
                            vc.var_labels[1],
                            matrix_type.capitalize(),
                        ],
                        plot=False,
                        save_fig=True,
                        save_name=path.join(
                            cat_dir,
                            "{}-{}-{}".format(vsn, tag, matrix_type),
                        ),
                    )
    return syst_dict


def run_syst_multisim_aggregate(
    chunks_dirs: Sequence[str],
    syst_disk_root: str,
    mc_df_stage: str = "final",
    var_set: str = "final",
    no_plots: bool = False,
    no_legacy_npz: bool = False,
) -> None:
    roots = [
        os.path.abspath(os.path.expanduser(str(d).rstrip(os.sep)))
        for d in chunks_dirs
        if str(d).strip()
    ]
    if not roots:
        raise SystemExit("[multisim-agg] no non-empty --chunks_dir paths provided")
    root = normalized_root(syst_disk_root)
    os.makedirs(root, exist_ok=True)
    save_fig = not no_plots

    print("[multisim-agg] chunks_dir(s)=%s  syst-disk-root=%s" % ("; ".join(roots), root))

    ck = collect_nu_chunks_many(roots)
    if not ck:
        raise RuntimeError(
            "[multisim-agg] no nu__*.pkl under any of: %s (expected per-root MCstat/, Flux/, "
            "G4/, Combined/ and/or legacy nu__*.pkl in each root)" % "; ".join(roots)
        )
    detected_stage = _detect_input_stage(ck)
    # Match the chunk var registry to the input stage so the aggregator knows the
    # bins for every var_save_name that the chunks produced.
    effective_var_set = "sel_all" if detected_stage == "sel_all" else var_set
    if effective_var_set != var_set:
        print(
            "[multisim-agg] input_stage=%s detected; overriding --var-set %r → %r"
            % (detected_stage, var_set, effective_var_set)
        )
    var_configs = build_var_configs(effective_var_set)
    print(
        "[multisim-agg] merging %d nu chunk(s) (input_stage=%s var_set=%s)"
        % (len(ck), detected_stage, effective_var_set)
    )
    merged = merge_nu_chunks(ck)
    syst_dict = covariance_dict_from_merged(merged, var_configs, root, save_fig)

    if not no_legacy_npz:
        save_neutrino_multisim_npzs(syst_dict, root)
        print("[multisim-agg] neutrino multisim NPZs ->", root)

    _write_covariance_manifest(
        root,
        syst_dict,
        merged["meta"],
        mc_df_stage,
        effective_var_set,
    )

    summ = path.join(root, "syst_multisim_aggregate_summary.txt")
    with open(summ, "w") as f:
        f.write("# merged nu chunk metadata\n")
        for m in merged["meta"]:
            f.write("%s\n" % m)
    print("[multisim-agg] wrote", summ)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--chunks_dir",
        dest="chunks_dirs",
        nargs="+",
        required=True,
        metavar="DIR",
        help="One or more chunk directories (each may contain MCstat/, Flux/, G4/, Combined/ or "
        "flat nu__*.pkl). The multisim driver passes multisim, g4, and flux chunk roots.",
    )
    p.add_argument(
        "--syst-disk-root",
        "--out_dir",
        "--out-dir",
        dest="syst_disk_root",
        required=True,
        help="Syst disk root (see analysis_village.numucc_1p0pi.syst_disk_layout): writes "
        "MCstat/, Flux/, and G4/ subfolders for neutrino multisim.",
    )
    p.add_argument("--mc-df-stage", choices=("final", "sel_all"), default="final")
    p.add_argument(
        "--var-set",
        choices=("final", "intermediate", "both", "sel_all"),
        default="final",
        help="Variable catalogue. ``sel_all`` selects pipeline-walker cut variables + "
        "final variables; auto-applied when sel_all chunk pickles are detected.",
    )
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--no-legacy-npz", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    run_syst_multisim_aggregate(
        chunks_dirs=args.chunks_dirs,
        syst_disk_root=args.syst_disk_root,
        mc_df_stage=args.mc_df_stage,
        var_set=args.var_set,
        no_plots=args.no_plots,
        no_legacy_npz=args.no_legacy_npz,
    )


if __name__ == "__main__":
    main()
