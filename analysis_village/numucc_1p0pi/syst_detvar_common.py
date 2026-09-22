"""Shared helpers for detector unisim matching + WireMod / DENT Product-B processing.

Matching is always at **sel_all** (before selection), delegated to
``dent_match_common_events``. Product B notebooks then walk the selection
pipeline on those matched files to fill final-stage histograms.

Also covers:

* asserting matched-file key equality
* batched hist accumulation (no full DF concat)
* WireMod total envelope (calo ± + efield) → unisim packs
* DENT CV-vs-var unisim packs
* combined Detector NPZ (WireMod YZ / XTXW + DENT + total)
* summary frac-unc comparison plots
"""
from __future__ import annotations

import gc
import json
import pickle
import re
from datetime import datetime
from os import makedirs, path
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

from pyanalib.covariance import get_covariance_matrix
from pyanalib.split_df_helpers_new import get_n_split

from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import FILE_DETECTOR, SUB_DETECTOR, normalized_root
from analysis_village.numucc_1p0pi.syst_histcounts import unisim_cov_from_cv_and_var
from analysis_village.numucc_1p0pi.syst_multisim_common import combine_indep_knob_cov_packs
from analysis_village.numucc_1p0pi.selection_framework import multicol_get_series
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig

EventKey = Tuple[float, int, int, int]  # (E / nuE, run, subrun, evt)

CALO_PARAMS = ("ccal", "alpha", "beta", "R")
WIREMOD_EFIELD_UNIV = "efield"
# Calo-only tables (legacy name kept for callers that only need ± calo stems).
WIREMOD_CALO_UNIVERSES = ("cv",) + tuple(f"{c}_{s}" for c in CALO_PARAMS for s in ("p", "m"))
# Full walk set: CV + eight calo ± + efield redo (``evt_efield`` / ``trk_efield``).
WIREMOD_UNIVERSES = WIREMOD_CALO_UNIVERSES + (WIREMOD_EFIELD_UNIV,)
# All WireMod universes (in-file cv + calo ± + efield) vs the external matched CV.
WIREMOD_ENVELOPE_SHIFTED = tuple(WIREMOD_UNIVERSES)
WIREMOD_KNOB_TAGS = {"YZ": "wiremod_yz", "XTXW": "wiremod_xtxw"}

# WireMod updatecalo products live in ``chi2_*_new`` (calo) / ``chi2_*_new_efield``.
_WIREMOD_CHI2_PLANES = ("I0", "I1", "I2")
_WIREMOD_CHI2_QUANTS = ("chi2_muon", "chi2_proton")


def wiremod_component_shifted_univs() -> Dict[str, Tuple[str, ...]]:
    """Per-component envelopes for inspection (not the final total).

    Keys: each calo parameter (± pair) and ``efield`` alone.
    """
    out: Dict[str, Tuple[str, ...]] = {
        p: (f"{p}_p", f"{p}_m") for p in CALO_PARAMS
    }
    out[WIREMOD_EFIELD_UNIV] = (WIREMOD_EFIELD_UNIV,)
    return out

SOURCE_COLORS = {
    "wiremod_yz": "#1f77b4",
    "wiremod_xtxw": "#ff7f0e",
    "DENT": "#2ca02c",
    "Detector": "k",
}
SOURCE_DISPLAY = {
    "wiremod_yz": "WireMod YZ",
    "wiremod_xtxw": "WireMod XTXW",
    "DENT": "DENT",
    "Detector": "Detector (total)",
}


def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def glob_matched_dfs(
    search_dir: Path | str,
    *,
    filename_str: str = "sel_mup",
    matched_suffix: str = "_matched",
) -> List[str]:
    """List matched DFs under *search_dir*.

    ``matched_suffix`` is the tag before ``.df`` (``_matched`` or ``_matched_hs``
    for the high-stats DENT campaign). Prefer names containing *filename_str*.
    """
    root = Path(search_dir)
    suf = matched_suffix if matched_suffix.startswith("_") else f"_{matched_suffix}"
    pattern = f"*{suf}.df"
    all_matched = sorted({str(p) for p in root.rglob(pattern) if p.is_file()})
    # Avoid picking ``*_matched_hs.df`` when asking for plain ``_matched``.
    if suf == "_matched":
        all_matched = [f for f in all_matched if "_matched_hs.df" not in Path(f).name]
    pref = [f for f in all_matched if filename_str in Path(f).name]
    return pref if pref else all_matched


def _matched_file_index(fpath: str | Path) -> Optional[int]:
    """Parse ``..._123_matched.df`` / ``..._123_matched_hs.df`` → 123."""
    m = re.search(r"_(\d+)_matched(?:_hs)?\.df$", Path(fpath).name)
    return int(m.group(1)) if m else None


def _per_evt_col(evt_df: pd.DataFrame, col_tuple) -> np.ndarray:
    try:
        return multicol_get_series(evt_df, col_tuple).to_numpy(dtype=float)
    except Exception:
        return np.array([], dtype=float)


def default_final_var_defs() -> Dict[str, dict]:
    """``{var_save_name: {label, bins, extract}}`` for final-selected kinematics."""
    configs = with_final_selected_evt_variables(list(CORE_SELECTED_EVT_VARIABLE_CONFIGS))
    names = {c.var_save_name for c in configs}
    if "integrated" not in names:
        configs = list(configs) + [VariableConfig.all_events()]

    out: Dict[str, dict] = {}
    for vc in configs:
        xlab = vc.var_labels[0] if getattr(vc, "var_labels", None) else vc.var_save_name
        if vc.var_save_name == "integrated":
            out["integrated"] = {
                "label": xlab,
                "bins": np.asarray(vc.bins, dtype=float),
                "extract": lambda df: np.full(len(df), 500.0, dtype=float),
                "var_config": vc,
            }
            continue
        col = vc.var_evt_reco_col
        out[vc.var_save_name] = {
            "label": xlab,
            "bins": np.asarray(vc.bins, dtype=float),
            "extract": lambda df, col=col: _per_evt_col(df, col),
            "var_config": vc,
        }
    return out


def collect_meta_event_keys(files: Sequence[str], *, max_files: Optional[int] = None) -> Set[EventKey]:
    """Union of ``(E, run, subrun, evt)`` from ``meta`` tables in matched/unmatched files."""
    from analysis_village.numucc_1p0pi.scripts.wiremod_match_common_events import (
        collect_event_keys,
    )

    use = list(files) if max_files is None else list(files)[: int(max_files)]
    return collect_event_keys(use)


def collect_sel_all_event_keys(files: Sequence[str], *, max_files: Optional[int] = None) -> Set[EventKey]:
    """Union of ``(E, run, subrun, evt)`` from sel_all ``hdr``+``evt`` tables."""
    from analysis_village.numucc_1p0pi.scripts.dent_match_common_events import (
        collect_sel_all_event_keys as _collect,
    )

    use = list(files) if max_files is None else list(files)[: int(max_files)]
    return _collect(use)


def assert_variations_matched(
    files_by_variation: Mapping[str, Sequence[str]],
    *,
    max_files_per_var: int = 50,
    format: str = "sel_mup",
    common_keys: Optional[Set[EventKey]] = None,
    common_keys_path: Optional[Path | str] = None,
) -> Set[EventKey]:
    """Assert matched products are consistent across variations.

    Detector matching builds a **global** common-key set, then filters each input
    file independently. Same file index on CV vs DENT need **not** share keys
    (CAF list orderings differ). Therefore:

    * If ``common_keys`` / ``common_keys_path`` is given — assert that keys from
      a sample of each variation's files are a **subset** of that set.
    * Else — fall back to equality of key unions (only reliable with a full scan
      or when file lists are truly pairwise-matched).
    """
    if len(files_by_variation) < 2:
        raise ValueError("need at least two variations to assert matching")
    key_fn = collect_sel_all_event_keys if format == "sel_all" else collect_meta_event_keys

    ref_keys: Optional[Set[EventKey]] = common_keys
    if ref_keys is None and common_keys_path is not None:
        path = Path(common_keys_path)
        if path.is_file():
            with open(path, "rb") as fh:
                loaded = pickle.load(fh)
            ref_keys = set(loaded)
            log(f"loaded common keys ({len(ref_keys):,}) from {path}")

    keysets: Dict[str, Set[EventKey]] = {}
    for lab, files in files_by_variation.items():
        if not files:
            raise AssertionError(f"{lab}: no matched files")
        use = list(files) if max_files_per_var is None else list(files)[: int(max_files_per_var)]
        keysets[lab] = key_fn(use, max_files=None)
        log(
            f"  {lab}: {len(files)} files, scan {len(use)} → "
            f"{len(keysets[lab]):,} event keys"
        )

    if ref_keys is not None:
        for lab, keys in keysets.items():
            extra = keys - ref_keys
            if extra:
                raise AssertionError(
                    f"{lab}: {len(extra):,} sampled keys not in common-key set "
                    f"(sample={len(keys):,}, common={len(ref_keys):,})"
                )
            if not keys:
                raise AssertionError(f"{lab}: empty key set in sampled matched files")
        log(
            f"assert OK: {len(files_by_variation)} variations ⊆ common keys "
            f"({len(ref_keys):,})"
        )
        return ref_keys

    # Equality fallback (pairwise / full-scan campaigns).
    ref_lab = next(iter(keysets))
    ref = keysets[ref_lab]
    for lab, keys in keysets.items():
        if keys != ref:
            only_ref = len(ref - keys)
            only_lab = len(keys - ref)
            raise AssertionError(
                f"matched-key mismatch: {ref_lab} vs {lab} "
                f"(only in {ref_lab}={only_ref}, only in {lab}={only_lab}). "
                f"For global-match campaigns pass common_keys_path="
                f"dent_common_keys_sel_all.pkl"
            )
    log(f"assert OK: {len(files_by_variation)} variations share {len(ref):,} event keys")
    return ref


def run_wiremod_match(
    variations: Mapping[str, str],
    *,
    filename_str: str = "sel_mup",
    phase: str = "all",
    common_keys_pkl: Optional[str] = None,
    summary_csv: Optional[str] = None,
    max_files: Optional[int] = None,
) -> int:
    """Run ``wiremod_match_common_events.main`` for sel_mup-style meta matching."""
    from analysis_village.numucc_1p0pi.scripts import wiremod_match_common_events as wm

    argv: List[str] = ["--phase", phase, "--filename-str", filename_str]
    for name, directory in variations.items():
        argv += ["--variation", name, str(directory)]
    if common_keys_pkl:
        argv += ["--common-keys-pkl", str(common_keys_pkl)]
    if summary_csv:
        argv += ["--summary-csv", str(summary_csv)]
    if max_files is not None:
        argv += ["--max-files", str(int(max_files))]
    return int(wm.main(argv))


def run_dent_match(
    variations: Mapping[str, str],
    *,
    fmt: str = "sel_all",
    filename_str: Optional[str] = None,
    summary_csv: Optional[str] = None,
    matched_out_dir: Optional[str] = None,
    matched_suffix: Optional[str] = None,
    max_files: Optional[int] = None,
    common_keys_pkl: Optional[str] = None,
    phase: Optional[str] = None,
    n_workers: Optional[int] = None,
    file_timeout: Optional[float] = None,
    meta_retries: Optional[int] = None,
    skip_existing_matched: bool = False,
) -> int:
    """Run ``dent_match_common_events.main`` (sel_all or sel_mup)."""
    from analysis_village.numucc_1p0pi.scripts import dent_match_common_events as dm

    filename_str = filename_str or fmt
    argv: List[str] = ["--format", fmt, "--filename-str", filename_str]
    for name, directory in variations.items():
        argv += ["--variation", name, str(directory)]
    if summary_csv:
        argv += ["--summary-csv", str(summary_csv)]
    if matched_out_dir:
        argv += ["--matched-out-dir", str(matched_out_dir)]
    if matched_suffix:
        argv += ["--matched-suffix", matched_suffix]
    if max_files is not None:
        argv += ["--max-files", str(int(max_files))]
    if common_keys_pkl:
        argv += ["--common-keys-pkl", str(common_keys_pkl)]
    if phase:
        argv += ["--phase", str(phase)]
    if n_workers is not None:
        argv += ["--n-workers", str(int(n_workers))]
    if file_timeout is not None:
        argv += ["--file-timeout", str(float(file_timeout))]
    if meta_retries is not None:
        argv += ["--meta-retries", str(int(meta_retries))]
    if skip_existing_matched:
        argv += ["--skip-existing-matched"]
    return int(dm.main(argv))


# ---------------------------------------------------------------------------
# POT + hist accumulation
# ---------------------------------------------------------------------------

def _series_col(df: pd.DataFrame, name: str):
    if name in df.columns:
        return df[name]
    for c in df.columns:
        if isinstance(c, tuple) and len(c) and c[0] == name:
            return df[c]
    return None


def sum_pot_from_matched_file(fpath: str) -> float:
    tot = 0.0
    try:
        with pd.HDFStore(fpath, "r") as store:
            for key in store.keys():
                name = key.lstrip("/")
                if not (name.startswith("hdr") or name.startswith("meta")):
                    continue
                df = store[name]
                pot = _series_col(df, "pot")
                if pot is None:
                    continue
                fis = _series_col(df, "first_in_subrun")
                if fis is not None:
                    tot += float(np.asarray(pot[fis == 1], dtype=float).sum())
                else:
                    tot += float(np.asarray(pot, dtype=float).sum())
    except Exception as ex:
        log(f"  pot read failed for {path.basename(fpath)}: {ex}")
    return tot


def sum_pot_from_matched_files(matched_files: Sequence[str]) -> float:
    return float(sum(sum_pot_from_matched_file(f) for f in matched_files))


def pot_scales_to_cv(pot_by_variation: Mapping[str, float], reference: str = "CV") -> Dict[str, float]:
    pot_ref = float(pot_by_variation.get(reference, 0.0))
    scales: Dict[str, float] = {}
    for lab, pot in pot_by_variation.items():
        pot = float(pot)
        if lab == reference:
            scales[lab] = 1.0
        elif pot > 0 and pot_ref > 0:
            scales[lab] = pot_ref / pot
        else:
            scales[lab] = 0.0
    return scales


def apply_pot_scales_to_hists(all_hists: dict, pot_scales: Mapping[str, float]) -> dict:
    for lab, scale in pot_scales.items():
        if lab not in all_hists or abs(float(scale) - 1.0) < 1e-15:
            continue
        for univ, var_hists in all_hists[lab].items():
            for var_name in list(var_hists.keys()):
                var_hists[var_name] = np.asarray(var_hists[var_name], dtype=float) * float(scale)
    return all_hists


def _fill_hists_from_split(df, var_defs, hists, weight: float = 1.0) -> int:
    if df is None or len(df) == 0:
        return 0
    w = float(weight)
    for var_name, cfg in var_defs.items():
        try:
            vals = cfg["extract"](df)
        except Exception:
            continue
        vals = np.asarray(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        bins = np.asarray(cfg["bins"], dtype=float)
        eps = (bins[-1] - bins[0]) * 1e-9
        vals = np.clip(vals, bins[0], bins[-1] - eps)
        counts, _ = np.histogram(vals, bins=bins)
        hists[var_name] += counts * w
    return len(df)


def accumulate_matched_sel_all_final(
    matched_files: Sequence[str],
    *,
    final_var_defs: Optional[Mapping[str, dict]] = None,
    include_cut_stage: bool = False,
) -> Tuple[dict, float]:
    """Batched hist fill from matched **sel_all** files via the selection pipeline.

    Loads ``evt_i`` / ``trk_i`` / ``hdr_i`` only (never concatenates across files),
    walks cuts with :func:`syst_pipeline_walker.walk_pipeline`, and fills Product **B**
    final kinematics at ``FINAL_STAGE_KEY``.  Matching must already have been done at
    sel_all so efficiency differences across variations are retained.
    """
    payload = accumulate_matched_sel_all_products(
        matched_files,
        final_var_defs=final_var_defs,
        include_cut_stage=include_cut_stage,
    )
    hists = {**payload["hists_cut"], **payload["hists_final"]} if include_cut_stage else payload["hists_final"]
    return hists, float(payload["pot"])


def accumulate_matched_sel_all_products(
    matched_files: Sequence[str],
    *,
    final_var_defs: Optional[Mapping[str, dict]] = None,
    cut_var_defs: Optional[Mapping[str, dict]] = None,
    include_cut_stage: bool = True,
    mu_p_candidate_kwargs: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Walk matched sel_all files; return cut-stage + final hist products + POT.

    Returns
    -------
    dict with keys:
      ``hists_cut``, ``hists_final``, ``pot``, ``cut_var_names``, ``final_var_names``
    """
    from analysis_village.numucc_1p0pi.scripts import dent_compare as dc

    final_defs = (
        dict(final_var_defs)
        if final_var_defs is not None
        else dict(dc.build_final_var_defs())
    )
    cut_defs = (
        dict(cut_var_defs)
        if cut_var_defs is not None
        else (dc.build_sel_all_var_defs() if include_cut_stage else {})
    )
    var_defs = {**cut_defs, **final_defs}
    hists = {v: np.zeros(len(cfg["bins"]) - 1, dtype=float) for v, cfg in var_defs.items()}
    stage_specs = dc._stage_specs_by_key() if cut_defs else {}
    summary = dc.SampleSummary(variation="", n_matched_files=len(matched_files))
    tot_pot = 0.0
    pid_kw = dict(mu_p_candidate_kwargs) if mu_p_candidate_kwargs else None

    for fpath in tqdm(matched_files, desc="matched sel_all walk"):
        pot = dc.process_sel_all_file(
            fpath,
            hists=hists,
            var_defs=var_defs,
            summary=summary,
            stage_specs=stage_specs,
            keyed_maps=None,
            final_var_defs=final_defs,
            mu_p_candidate_kwargs=pid_kw,
        )
        tot_pot += float(pot)
        gc.collect()

    hists_cut = {k: hists[k] for k in cut_defs if k in hists}
    hists_final = {k: hists[k] for k in final_defs if k in hists}
    return {
        "hists_cut": hists_cut,
        "hists_final": hists_final,
        "pot": tot_pot,
        "cut_var_names": list(cut_defs.keys()),
        "final_var_names": list(final_defs.keys()),
        "stages": {k: sm.__dict__ for k, sm in summary.stages.items()},
        "mu_p_candidate_kwargs": dict(pid_kw) if pid_kw else {},
    }


def accumulate_matched_sel_all_cv_vs_var(
    matched_files_cv: Sequence[str],
    matched_files_var: Sequence[str],
    *,
    final_var_defs: Optional[Mapping[str, dict]] = None,
    include_cut_stage: bool = True,
) -> Tuple[dict, dict, float, float]:
    """DENT: walk matched sel_all CV and variation → combined cut+final hist dicts."""
    p_cv = accumulate_matched_sel_all_products(
        matched_files_cv,
        final_var_defs=final_var_defs,
        include_cut_stage=include_cut_stage,
    )
    p_var = accumulate_matched_sel_all_products(
        matched_files_var,
        final_var_defs=final_var_defs,
        include_cut_stage=include_cut_stage,
    )
    h_cv = {**p_cv["hists_cut"], **p_cv["hists_final"]}
    h_var = {**p_var["hists_cut"], **p_var["hists_final"]}
    return h_cv, h_var, float(p_cv["pot"]), float(p_var["pot"])


def accumulate_matched_sel_all_cv_vs_var_products(
    matched_files_cv: Sequence[str],
    matched_files_var: Sequence[str],
    *,
    final_var_defs: Optional[Mapping[str, dict]] = None,
    include_cut_stage: bool = True,
) -> dict:
    """DENT dual-product walk: returns structured CV/VAR cut+final payloads."""
    p_cv = accumulate_matched_sel_all_products(
        matched_files_cv,
        final_var_defs=final_var_defs,
        include_cut_stage=include_cut_stage,
    )
    p_var = accumulate_matched_sel_all_products(
        matched_files_var,
        final_var_defs=final_var_defs,
        include_cut_stage=include_cut_stage,
    )
    return {"cv": p_cv, "var": p_var}


# ---------------------------------------------------------------------------
# WireMod: require calo universes + walk each universe
# ---------------------------------------------------------------------------

def list_hdf_stems(fpath: str) -> Set[str]:
    with pd.HDFStore(fpath, "r") as store:
        names = {k.lstrip("/") for k in store.keys()}
    stems: Set[str] = set()
    for n in names:
        if "_" not in n:
            stems.add(n)
            continue
        # strip trailing _<int> split index
        head, _, tail = n.rpartition("_")
        if tail.isdigit() and head:
            stems.add(head)
        else:
            stems.add(n)
    return stems


def assert_wiremod_calo_universes(
    matched_files: Sequence[str],
    *,
    sample_files: int = 3,
    require_efield: bool = True,
) -> None:
    """Abort unless matched files contain ``evt_cv``, all calo ±, and (by default) efield."""
    if not matched_files:
        raise RuntimeError("WireMod: no matched files")
    required = set(WIREMOD_CALO_UNIVERSES)
    if require_efield:
        required.add(WIREMOD_EFIELD_UNIV)
    missing_any: Dict[str, Set[str]] = {}
    for fpath in list(matched_files)[: max(int(sample_files), 1)]:
        stems = list_hdf_stems(fpath)
        evt_univs = {s[len("evt_"):] for s in stems if s.startswith("evt_")}
        # also accept bare evt as cv only for detection — still require calo
        have = set(evt_univs)
        if "evt" in stems and "cv" not in have:
            have.add("cv")
        missing = required - have
        if missing:
            missing_any[path.basename(fpath)] = missing
    if missing_any:
        detail = "; ".join(f"{f}: missing {sorted(m)}" for f, m in missing_any.items())
        need = "evt_cv + evt_<ccal|alpha|beta|R>_{p,m}"
        if require_efield:
            need += " + evt_efield"
        raise RuntimeError(
            f"WireMod requires variation tables ({need}) in matched sel_all files. {detail}"
        )
    log(
        f"WireMod universes OK ({len(required)} univs"
        f"{', +efield' if require_efield else ''}) on "
        f"{min(len(matched_files), sample_files)} file(s)"
    )


assert_wiremod_universes = assert_wiremod_calo_universes


def apply_wiremod_chi2_variation(trk: Optional[pd.DataFrame], univ: str) -> Optional[pd.DataFrame]:
    """Overwrite standard ``chi2_*`` columns with WireMod variation products.

    Selection / PID use ``chi2_muon`` / ``chi2_proton``. WireMod updatecalo stores
    the varied values in ``chi2_*_new`` (calo univ tables, including in-file cv)
    and ``chi2_*_new_efield`` (efield table). Without this remap every universe
    walks identical standard χ² and the envelope collapses.
    """
    if trk is None or len(trk) == 0:
        return trk
    suf = "_new_efield" if univ == WIREMOD_EFIELD_UNIV else "_new"
    trk = trk.copy()
    n_applied = 0
    for plane in _WIREMOD_CHI2_PLANES:
        for quant in _WIREMOD_CHI2_QUANTS:
            src = ("pfp", "trk", "chi2pid", plane, f"{quant}{suf}", "")
            dst = ("pfp", "trk", "chi2pid", plane, quant, "")
            if src in trk.columns and dst in trk.columns:
                trk[dst] = trk[src]
                n_applied += 1
    if n_applied == 0:
        raise RuntimeError(
            f"WireMod chi2 remap failed for univ={univ!r}: expected columns "
            f"*{{chi2_muon,chi2_proton}}{suf} on planes {_WIREMOD_CHI2_PLANES}"
        )
    return trk


def _load_univ_sel_all_tables(fpath: str, split_i: int, univ: str):
    """Load (evt, trk, hdr) for one universe from a matched sel_all(+calo) file."""
    hdr = None
    for hk in (f"hdr_{split_i}", f"hdr_{univ}_{split_i}"):
        try:
            hdr = pd.read_hdf(fpath, key=hk)
            break
        except Exception:
            continue
    evt = None
    for ek in (f"evt_{univ}_{split_i}", f"evt_{split_i}" if univ == "cv" else None):
        if ek is None:
            continue
        try:
            evt = pd.read_hdf(fpath, key=ek)
            break
        except Exception:
            continue
    trk = None
    for tk in (f"trk_{univ}_{split_i}", f"trk_{split_i}", f"trk_cv_{split_i}"):
        try:
            trk = pd.read_hdf(fpath, key=tk)
            break
        except Exception:
            continue
    if trk is not None:
        trk = apply_wiremod_chi2_variation(trk, univ)
    return evt, trk, hdr


def _walk_fill_state(
    evt,
    trk,
    hdr,
    *,
    hists: dict,
    var_defs: Mapping[str, dict],
    stage_specs: dict,
    final_defs: Mapping[str, dict],
    mu_p_candidate_kwargs: Optional[Mapping[str, Any]] = None,
) -> float:
    """Fill cut+final hists by walking one in-memory sel_all split (DENT-compatible)."""
    from pyanalib.variable_calculator import add_reco_cc1p0pi_tki_evtdf
    from analysis_village.numucc_1p0pi.event_selection_batch_core import (
        attach_intrinsic_weights,
        ensure_phi_and_kinematics_cols,
        hdr_chunk_pot,
    )
    from analysis_village.numucc_1p0pi.syst_pipeline_walker import (
        FINAL_STAGE_KEY,
        get_var_series,
        histogram_var,
        walk_pipeline,
    )

    if evt is None or len(evt) == 0:
        return float(hdr_chunk_pot(hdr))
    pot = float(hdr_chunk_pot(hdr))
    attach_intrinsic_weights(evt, trk, "mc", use_mc_genweight=False)
    evt, _ = ensure_phi_and_kinematics_cols(evt, trk, None)
    state = {"evt": evt, "trk": trk, "hdr": hdr, "mcnu": None}
    if mu_p_candidate_kwargs:
        state["_mu_p_candidate_kwargs"] = dict(mu_p_candidate_kwargs)
    for stage_key, cur in walk_pipeline(state, sample="mc"):
        cur_evt = cur.get("evt")
        for var_name, vc, target in stage_specs.get(stage_key, []):
            if var_name not in hists:
                continue
            got = get_var_series(cur, vc, target)
            if got is None:
                continue
            vals, _ = got
            hists[var_name] += histogram_var(vals, var_defs[var_name]["bins"])
        if stage_key == FINAL_STAGE_KEY and final_defs and cur_evt is not None and len(cur_evt) > 0:
            pe = add_reco_cc1p0pi_tki_evtdf(cur_evt)
            for var_name, cfg in final_defs.items():
                if var_name not in hists:
                    continue
                try:
                    vals = cfg["extract"](pe)
                except Exception:
                    continue
                if len(vals) == 0:
                    continue
                hists[var_name] += histogram_var(vals, cfg["bins"])
    return pot


def _filter_tables_drop_entries(evt, trk, hdr, drop_entries: set):
    """Drop rows whose (__ntuple, entry) is in *drop_entries* (cross-file artkey dups)."""
    if not drop_entries:
        return evt, trk, hdr

    def _mask(df):
        if df is None or len(df) == 0:
            return None if df is None else df
        df_r = df.reset_index()
        if "__ntuple" not in df_r.columns or "entry" not in df_r.columns:
            return df
        keep = np.fromiter(
            (
                (int(nt), int(en)) not in drop_entries
                for nt, en in zip(df_r["__ntuple"], df_r["entry"])
            ),
            dtype=bool,
            count=len(df_r),
        )
        if keep.all():
            return df
        if not keep.any():
            return df.iloc[0:0]
        return df.iloc[np.flatnonzero(keep)]

    return _mask(evt), _mask(trk), _mask(hdr)


def accumulate_wiremod_matched_products(
    matched_files: Sequence[str],
    *,
    universes: Sequence[str] = WIREMOD_UNIVERSES,
    final_var_defs: Optional[Mapping[str, dict]] = None,
    include_cut_stage: bool = True,
    drop_map: Optional[Mapping[str, Mapping[int, set]]] = None,
    mu_p_candidate_kwargs: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Walk matched WireMod sel_all(+calo/efield) files for every universe.

    Default *universes* are CV + eight calo ± + ``efield``. Requires those
    tables (aborts via :func:`assert_wiremod_calo_universes`).  Returns
    ``{univ: {hists_cut, hists_final}, pot, ...}``.

    *drop_map* (optional): ``{abspath: {split_i: set[(__ntuple, entry)]}}`` from
    :mod:`dedupe_matched_artkeys` — skips cross-file duplicate artkeys.

    *mu_p_candidate_kwargs* (optional): forwarded to ``get_mu_p_candidate``
    (e.g. ``{"mu_chi2mu_th": 25}``). Default selection uses ``MU_CHI2MU_TH=30``.
    """
    from analysis_village.numucc_1p0pi.scripts import dent_compare as dc

    assert_wiremod_calo_universes(
        matched_files, require_efield=(WIREMOD_EFIELD_UNIV in set(universes))
    )
    final_defs = (
        dict(final_var_defs)
        if final_var_defs is not None
        else dict(dc.build_final_var_defs())
    )
    cut_defs = dc.build_sel_all_var_defs() if include_cut_stage else {}
    var_defs = {**cut_defs, **final_defs}
    stage_specs = dc._stage_specs_by_key() if cut_defs else {}
    pid_kw = dict(mu_p_candidate_kwargs) if mu_p_candidate_kwargs else None

    per_univ = {
        u: {v: np.zeros(len(cfg["bins"]) - 1, dtype=float) for v, cfg in var_defs.items()}
        for u in universes
    }
    tot_pot = 0.0
    n_filled = 0

    for fpath in tqdm(matched_files, desc="WireMod sel_all+calo walk"):
        try:
            n_split = get_n_split(fpath)
        except Exception as ex:
            log(f"  skip {path.basename(fpath)}: {ex}")
            continue
        abs_f = path.abspath(fpath)
        file_drops = (drop_map or {}).get(abs_f) or {}
        # POT once per file (from hdr); informational when envelopes are unscaled
        tot_pot += sum_pot_from_matched_file(fpath)
        for i in range(n_split):
            drop_entries = file_drops.get(i) or file_drops.get(str(i)) or set()
            for univ in universes:
                evt, trk, hdr = _load_univ_sel_all_tables(fpath, i, univ)
                if evt is None:
                    continue
                if drop_entries:
                    evt, trk, hdr = _filter_tables_drop_entries(evt, trk, hdr, drop_entries)
                    if evt is None or len(evt) == 0:
                        continue
                # Embedded trk1/trk2 without a trk table: cannot do true sel_all walk.
                if trk is None:
                    top = set(evt.columns.get_level_values(0).unique()) if hasattr(evt.columns, "get_level_values") else set()
                    if "trk1" in top or "trk2" in top:
                        raise RuntimeError(
                            f"WireMod file {path.basename(fpath)} has evt_{univ} with "
                            "embedded trk1/trk2 but no sel_all trk table. Remake as "
                            "sel_all+updatecalo (evt_* + trk_*/trk) so the full "
                            "selection walk can run."
                        )
                    raise RuntimeError(
                        f"WireMod: missing trk table for univ={univ} in {path.basename(fpath)}"
                    )
                _walk_fill_state(
                    evt,
                    trk,
                    hdr,
                    hists=per_univ[univ],
                    var_defs=var_defs,
                    stage_specs=stage_specs,
                    final_defs=final_defs,
                    mu_p_candidate_kwargs=pid_kw,
                )
                n_filled += 1
                del evt, trk, hdr
            gc.collect()

    if n_filled == 0:
        raise RuntimeError("WireMod walk filled zero universe splits — check matched file keys")

    out_univ = {}
    for u, h in per_univ.items():
        out_univ[u] = {
            "hists_cut": {k: h[k] for k in cut_defs if k in h},
            "hists_final": {k: h[k] for k in final_defs if k in h},
        }
    return {
        "by_universe": out_univ,
        "pot": tot_pot,
        "cut_var_names": list(cut_defs.keys()),
        "final_var_names": list(final_defs.keys()),
        "universes": list(universes),
        "mu_p_candidate_kwargs": dict(pid_kw) if pid_kw else {},
    }


def wiremod_geometry_hists_for_envelope(by_universe: Mapping[str, dict], *, product: str) -> dict:
    """``{univ: {var: hist}}`` for :func:`build_wiremod_detector_dict` (*product* = cut|final)."""
    key = "hists_cut" if product == "cut" else "hists_final"
    return {u: dict(payload[key]) for u, payload in by_universe.items()}


def accumulate_matched_universes(
    matched_files: Sequence[str],
    universes: Sequence[str],
    var_defs: Mapping[str, dict],
    *,
    pot_scale: float = 1.0,
    key_fmt: str = "evt_{univ}_{i}",
    df_prepare=None,
) -> Tuple[dict, dict, float]:
    """Legacy hist fill from ``evt_<univ>_<split>`` without a selection walk."""
    hists = {
        u: {v: np.zeros(len(cfg["bins"]) - 1, dtype=float) for v, cfg in var_defs.items()}
        for u in universes
    }
    n_evts = {u: 0 for u in universes}
    tot_pot = 0.0

    for fpath in tqdm(matched_files, desc="matched files"):
        tot_pot += sum_pot_from_matched_file(fpath)
        try:
            n_split = get_n_split(fpath)
        except Exception as ex:
            log(f"  skip {path.basename(fpath)}: {ex}")
            continue
        for i in range(n_split):
            for univ in universes:
                key = key_fmt.format(univ=univ, i=i)
                try:
                    df = pd.read_hdf(fpath, key=key)
                except Exception:
                    continue
                if df_prepare is not None:
                    try:
                        df = df_prepare(df)
                    except Exception:
                        del df
                        continue
                n = _fill_hists_from_split(df, var_defs, hists[univ], weight=pot_scale)
                n_evts[univ] += n
                del df
            gc.collect()
    return hists, n_evts, tot_pot


def accumulate_matched_cv_vs_var(
    matched_files_cv: Sequence[str],
    matched_files_var: Sequence[str],
    var_defs: Mapping[str, dict],
    *,
    cv_key_fmt: str = "evt_cv_{i}",
    var_key_fmt: str = "evt_{i}",
    df_prepare=None,
) -> Tuple[dict, dict, float, float]:
    """DENT-style: CV and variation as separate file lists (often one univ table each).

    Tries ``evt_cv_i`` / ``evt_i`` then falls back to ``evt_0``-style keys present in the file.
    """
    def _accum(files, key_fmt, label):
        hists = {v: np.zeros(len(cfg["bins"]) - 1, dtype=float) for v, cfg in var_defs.items()}
        n_evts = 0
        tot_pot = 0.0
        for fpath in tqdm(files, desc=f"{label} files"):
            tot_pot += sum_pot_from_matched_file(fpath)
            try:
                n_split = get_n_split(fpath)
            except Exception as ex:
                log(f"  skip {path.basename(fpath)}: {ex}")
                continue
            with pd.HDFStore(fpath, "r") as store:
                keys = set(k.lstrip("/") for k in store.keys())
            for i in range(n_split):
                candidates = [
                    key_fmt.format(i=i),
                    f"evt_{i}",
                    f"evt_cv_{i}",
                    "evt",
                ]
                df = None
                for key in candidates:
                    if key not in keys and f"/{key}" not in {f"/{k}" for k in keys}:
                        # allow missing from keys set check
                        pass
                    try:
                        df = pd.read_hdf(fpath, key=key)
                        break
                    except Exception:
                        continue
                if df is None:
                    continue
                if df_prepare is not None:
                    try:
                        df = df_prepare(df)
                    except Exception:
                        del df
                        continue
                n_evts += _fill_hists_from_split(df, var_defs, hists, weight=1.0)
                del df
            gc.collect()
        return hists, n_evts, tot_pot

    h_cv, n_cv, pot_cv = _accum(matched_files_cv, cv_key_fmt, "CV")
    h_var, n_var, pot_var = _accum(matched_files_var, var_key_fmt, "VAR")
    return h_cv, h_var, pot_cv, pot_var


# ---------------------------------------------------------------------------
# Cache I/O
# ---------------------------------------------------------------------------

def save_hist_cache(cache_path: Path | str, payload: dict) -> Path:
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    # strip non-pickleable extract callables from var_defs
    safe = dict(payload)
    if "var_defs" in safe:
        safe["var_defs"] = {
            k: {"label": v.get("label"), "bins": np.asarray(v["bins"])}
            for k, v in safe["var_defs"].items()
        }
    with open(cache_path, "wb") as fh:
        pickle.dump(safe, fh, protocol=pickle.HIGHEST_PROTOCOL)
    log(f"saved hist cache → {cache_path}")
    return cache_path


def load_hist_cache(cache_path: Path | str) -> dict:
    with open(cache_path, "rb") as fh:
        payload = pickle.load(fh)
    log(f"loaded hist cache ← {cache_path}")
    return payload


def load_or_build_cache(cache_path: Path | str, builder, *, force: bool = False) -> dict:
    cache_path = Path(cache_path)
    if cache_path.is_file() and not force:
        return load_hist_cache(cache_path)
    payload = builder()
    save_hist_cache(cache_path, payload)
    return payload


# ---------------------------------------------------------------------------
# WireMod envelope unisim
# ---------------------------------------------------------------------------

def sanitize_matrix_pack(ret: Mapping[str, Any]) -> dict:
    return {
        k: np.nan_to_num(np.asarray(ret[k], dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
        for k in ("cov", "cov_frac", "corr")
    }


def max_envelope_univ_counts(n_cv, hists, var_name, shifted_univs) -> np.ndarray:
    """Per-bin counts used to encode ``max_u |n_u - n_cv|`` in a two-univ pack.

    Plots show the **actual** min/max among universes. The uncertainty on each
    bin is the larger absolute deviation from CV (either direction). Returns
    ``n_cv + delta`` so ``|n_var - n_cv| = delta`` in the unisim cov.
    """
    n_cv = np.asarray(n_cv, dtype=float)
    keys = [u for u in shifted_univs if u in hists and var_name in hists[u]]
    if not keys:
        return n_cv.copy()
    stacked = np.stack([np.asarray(hists[u][var_name], dtype=float) for u in keys], axis=0)
    delta = np.maximum(np.abs(stacked.max(axis=0) - n_cv), np.abs(stacked.min(axis=0) - n_cv))
    return n_cv + delta


def envelope_delta_counts(n_cv, hists, var_name, shifted_univs) -> np.ndarray:
    """Per-bin ``max_u |n_u - n_cv|`` (uncertainty magnitude, not a plot band)."""
    n_cv = np.asarray(n_cv, dtype=float)
    keys = [u for u in shifted_univs if u in hists and var_name in hists[u]]
    if not keys:
        return np.zeros_like(n_cv)
    stacked = np.stack([np.asarray(hists[u][var_name], dtype=float) for u in keys], axis=0)
    return np.maximum(np.abs(stacked.max(axis=0) - n_cv), np.abs(stacked.min(axis=0) - n_cv))


def envelope_univ_lo_hi(hists, var_name, shifted_univs) -> tuple[np.ndarray, np.ndarray]:
    """Per-bin actual min/max counts among *shifted_univs* (no symmetrization)."""
    keys = [u for u in shifted_univs if u in hists and var_name in hists[u]]
    if not keys:
        raise KeyError(f"no universes in {list(shifted_univs)} for {var_name}")
    stacked = np.stack([np.asarray(hists[u][var_name], dtype=float) for u in keys], axis=0)
    return stacked.min(axis=0), stacked.max(axis=0)


def cov_pack_two_universe(n_cv, n_var) -> dict:
    ret = get_covariance_matrix(
        np.asarray([n_var], dtype=float),
        np.asarray(n_cv, dtype=float),
    )
    return sanitize_matrix_pack(ret)


def build_wiremod_detector_dict(
    all_hists: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    var_names: Sequence[str],
    *,
    wiremod_labels: Sequence[str] = ("YZ", "XTXW"),
    knob_tags: Optional[Mapping[str, str]] = None,
    shifted_univs: Optional[Sequence[str]] = None,
    cv_hists: Optional[Mapping[str, np.ndarray]] = None,
) -> dict:
    """WireMod-only detector dict (per-geometry + combined WireMod total).

    Default *shifted_univs* is the **total** envelope: in-file cv + calo ± + efield.
    *cv_hists* is the **external matched CV sample** ``{var: counts}`` (Sep-4).
    Per-bin unc = ``max_u |n_u - n_cv| / n_cv``. If *cv_hists* is omitted, falls
    back to each geometry's in-file ``cv``.
    """
    knob_tags = dict(knob_tags or WIREMOD_KNOB_TAGS)
    shifted = list(shifted_univs) if shifted_univs is not None else list(WIREMOD_ENVELOPE_SHIFTED)

    detector_dict: Dict[str, Any] = {"detector": {}}
    for lab in wiremod_labels:
        detector_dict[f"detector-{knob_tags[lab]}"] = {}
    detector_by_wiremod: Dict[str, dict] = {}

    for var_name in var_names:
        per_tag = {}
        packs = []
        n_cv_ext = None
        if cv_hists is not None and var_name in cv_hists:
            n_cv_ext = np.asarray(cv_hists[var_name], dtype=float)
            if float(n_cv_ext.sum()) <= 0:
                n_cv_ext = None
        for lab in wiremod_labels:
            if lab not in all_hists:
                continue
            hists = all_hists[lab]
            if n_cv_ext is not None:
                n_cv = n_cv_ext
            else:
                if "cv" not in hists or var_name not in hists["cv"]:
                    continue
                n_cv = np.asarray(hists["cv"][var_name], dtype=float)
                if float(n_cv.sum()) <= 0:
                    continue
            if not any(u in hists and var_name in hists[u] for u in shifted):
                continue
            n_var = max_envelope_univ_counts(n_cv, hists, var_name, shifted)
            pack = cov_pack_two_universe(n_cv, n_var)
            tag = knob_tags[lab]
            per_tag[tag] = pack
            detector_dict[f"detector-{tag}"][var_name] = pack
            packs.append(pack)
        if not packs:
            continue
        if n_cv_ext is not None:
            n_cv_ref = n_cv_ext
        else:
            lab0 = next(lab for lab in wiremod_labels if lab in all_hists and "cv" in all_hists[lab])
            n_cv_ref = np.asarray(all_hists[lab0]["cv"][var_name], dtype=float)
        combined = combine_indep_knob_cov_packs(packs, n_cv_ref)
        detector_dict["detector"][var_name] = sanitize_matrix_pack(combined)
        detector_by_wiremod[var_name] = per_tag

    if detector_by_wiremod:
        detector_dict["detector_by_wiremod"] = detector_by_wiremod
    return detector_dict


def build_dent_detector_dict(
    h_cv: Mapping[str, np.ndarray],
    h_dent: Mapping[str, np.ndarray],
    var_names: Sequence[str],
) -> dict:
    """Pure unisim DENT packs under ``detector-DENT`` (+ single-knob ``detector``)."""
    out: Dict[str, Any] = {"detector": {}, "detector-DENT": {}, "detector_by_wiremod": {}}
    for var_name in var_names:
        if var_name not in h_cv or var_name not in h_dent:
            continue
        n_cv = np.asarray(h_cv[var_name], dtype=float)
        n_var = np.asarray(h_dent[var_name], dtype=float)
        if n_cv.shape != n_var.shape or float(n_cv.sum()) <= 0:
            continue
        pack = sanitize_matrix_pack(unisim_cov_from_cv_and_var(n_cv, n_var))
        out["detector-DENT"][var_name] = pack
        out["detector"][var_name] = pack
        out["detector_by_wiremod"][var_name] = {"DENT": pack}
    return out


def combine_wiremod_dent_detector_dict(
    wiremod_dict: Mapping[str, Any],
    dent_dict: Mapping[str, Any],
) -> dict:
    """Merge WireMod YZ/XTXW + DENT into one NPZ payload with combined Detector total."""
    out: Dict[str, Any] = {
        "detector": {},
        "detector_by_wiremod": {},
    }
    # copy wiremod knobs
    for key, val in wiremod_dict.items():
        if key in ("detector", "detector_by_wiremod"):
            continue
        if isinstance(val, dict):
            out[key] = dict(val)
    # copy DENT knob
    if "detector-DENT" in dent_dict:
        out["detector-DENT"] = dict(dent_dict["detector-DENT"])

    vars_wm = set(wiremod_dict.get("detector", {}).keys())
    vars_dent = set(dent_dict.get("detector-DENT", {}).keys())
    for vsn in sorted(vars_wm | vars_dent):
        packs = []
        cv_ref = None
        by_knob = {}
        wm_by = wiremod_dict.get("detector_by_wiremod", {}).get(vsn, {})
        for tag, pack in wm_by.items():
            by_knob[tag] = pack
            packs.append(pack)
            if cv_ref is None and "cov" in pack:
                # recover a CV scale from diag if needed — use ones; combine uses cov_frac
                nb = np.asarray(pack["cov_frac"]).shape[0]
                cv_ref = np.ones(nb, dtype=float)
        dent_pack = dent_dict.get("detector-DENT", {}).get(vsn)
        if dent_pack is not None:
            by_knob["DENT"] = dent_pack
            packs.append(dent_pack)
            if cv_ref is None:
                nb = np.asarray(dent_pack["cov_frac"]).shape[0]
                cv_ref = np.ones(nb, dtype=float)
        if not packs or cv_ref is None:
            continue
        out["detector"][vsn] = sanitize_matrix_pack(
            combine_indep_knob_cov_packs(packs, cv_ref)
        )
        out["detector_by_wiremod"][vsn] = by_knob
    return out


def load_detector_dict_npz(path: Path | str) -> dict:
    """Load ``detector_syst_dict.npz`` into plain nested dicts."""
    z = np.load(path, allow_pickle=True)
    out = {}
    for k in z.files:
        val = z[k]
        out[k] = val.item() if hasattr(val, "item") else val
    log(f"loaded detector npz ← {path}")
    return out


def save_detector_npz(detector_dict: Mapping[str, Any], out_path: Path | str, *, manifest: Optional[dict] = None) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_path, **detector_dict)
    log(f"wrote {out_path}")
    if manifest is not None:
        man_path = out_path.parent.parent / "detector_covariance_manifest.json"
        if out_path.parent.name == SUB_DETECTOR:
            man_path = out_path.parent.parent / "detector_covariance_manifest.json"
        else:
            man_path = out_path.with_name(out_path.stem + "_manifest.json")
        manifest = dict(manifest)
        manifest["npz_path"] = str(out_path)
        manifest["produced"] = datetime.now().isoformat(timespec="seconds")
        with open(man_path, "w") as fh:
            json.dump(manifest, fh, indent=2)
            fh.write("\n")
        log(f"wrote {man_path}")
    return out_path


# ---------------------------------------------------------------------------
# Summary plots (DENT vs WireMod)
# ---------------------------------------------------------------------------

def frac_unc_pct_from_pack(pack: Mapping[str, Any]) -> np.ndarray:
    cf = np.asarray(pack["cov_frac"], dtype=float)
    return 100.0 * np.sqrt(np.maximum(np.diag(cf), 0.0))


def plot_detector_source_frac_unc(
    detector_dict: Mapping[str, Any],
    var_config,
    *,
    ax=None,
    source_order: Sequence[str] = ("wiremod_yz", "wiremod_xtxw", "DENT", "Detector"),
    save_path: Optional[Path | str] = None,
    dpi: int = 140,
):
    """Overlay WireMod YZ / XTXW / DENT / total fractional unc for one variable."""
    vsn = var_config.var_save_name
    show = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
    else:
        fig = ax.figure

    centers = np.asarray(var_config.bin_centers, dtype=float)
    bins = np.asarray(var_config.bins, dtype=float)
    is_int = vsn == "integrated" or len(centers) == 1

    for src in source_order:
        if src == "Detector":
            pack = detector_dict.get("detector", {}).get(vsn)
        else:
            pack = detector_dict.get(f"detector-{src}", {}).get(vsn)
        if pack is None:
            continue
        w = frac_unc_pct_from_pack(pack)
        if is_int:
            w = np.full_like(w, float(w[0]))
        ax.hist(
            centers, bins=bins, weights=w, histtype="step", linewidth=1.8,
            color=SOURCE_COLORS.get(src, "gray"),
            label=SOURCE_DISPLAY.get(src, src),
        )

    if is_int:
        ax.set_xlabel("All Events")
        ax.set_xticks([centers[0]])
        ax.set_xticklabels(["All Events"])
    else:
        xlab = var_config.var_labels[1] if var_config.var_labels else vsn
        ax.set_xlabel(xlab)
    ax.set_ylabel("Uncertainty [%]")
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        log(f"wrote {save_path}")
    if show:
        plt.show()
    return fig, ax
