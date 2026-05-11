"""Shared helpers for Flux / G4 / MCstat covariance workflows (monolithic or chunked)."""

from __future__ import annotations

import os
import numpy as np
import pandas as pd

from analysis_village.numucc_1p0pi.syst_disk_layout import (
    FILE_COSMICS,
    FILE_FLUX,
    FILE_G4,
    FILE_MCSTAT,
    SUB_COSMICS,
    SUB_FLUX,
    SUB_G4,
    SUB_MCSTAT,
    normalized_root,
)
from analysis_village.numucc_1p0pi.selection_framework import multicol_resolve_column_key
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    INTERMEDIATE_CUT_SYST_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)


def g4_mc_knob_names():
    """Geant4 reinteraction knob names (same list as ``makedf.g4syst.g4_systematics``)."""
    from makedf.g4syst import g4_systematics

    return tuple(str(k) for k in g4_systematics)


def flux_mc_knob_names(group_spec: str = "all") -> tuple[str, ...]:
    """BNB flux regen knob names under ``(mc, <knob>, univ_i)`` (see ``makedf.bnbsyst``).

    * ``all`` / empty / ``*`` → ``makedf.bnbsyst.regen_systematics`` (full list for
      all-knobs-in-one-table builds).
    * ``beam``, ``hadron``, ``xsec``, or comma-separated combinations → union of the
      corresponding ``BNB_FLUX_GROUPS`` knob lists (deduped, order preserved).
    """
    from makedf import bnbsyst

    s = (group_spec or "all").strip().lower()
    if s in ("", "all", "*"):
        return tuple(str(k) for k in bnbsyst.regen_systematics)
    out: list[str] = []
    for g in [x.strip().lower() for x in s.split(",") if x.strip()]:
        if g in bnbsyst.BNB_FLUX_GROUPS:
            out.extend(str(k) for k in bnbsyst.BNB_FLUX_GROUPS[g][0])
        else:
            raise ValueError(
                "unknown flux knob group %r (use all, beam, hadron, xsec, or comma-separated)"
                % (g,)
            )
    seen: set[str] = set()
    uniq: list[str] = []
    for k in out:
        if k not in seen:
            seen.add(k)
            uniq.append(k)
    return tuple(uniq)


def drop_bad_g4_weights(mc_evt_df, max_wgt=1e3, n_univ=100):
    """Drop evt rows with any bundled ``mc.G4`` universe weight above ``max_wgt`` (no-op if missing)."""
    df = mc_evt_df
    try:
        col = df.mc.G4
    except (AttributeError, KeyError):
        return df
    bad_mask = pd.Series(False, index=df.index)
    for i in range(n_univ):
        try:
            var = col["univ_{}".format(i)]
        except (KeyError, TypeError):
            return df
        bad_mask = bad_mask | (var > max_wgt)
    if not bad_mask.any():
        return df
    return df.loc[~bad_mask]


def drop_bad_flux_knob_weights(mc_evt_df, knobs=None, max_wgt=1e3, n_univ=100):
    """Drop rows where any listed ``(mc, flux_knob, univ_i)`` weight exceeds ``max_wgt``."""
    df = mc_evt_df
    if knobs is None:
        knobs = flux_mc_knob_names()
    if not knobs:
        return df
    bad_mask = pd.Series(False, index=df.index)
    for knob in knobs:
        sk = ("mc", knob)
        for i in range(int(n_univ)):
            probe = sk + ("univ_{}".format(i),)
            key = multicol_resolve_column_key(df, probe)
            if key is None:
                break
            try:
                w = df.loc[:, key]
                bad_mask = bad_mask | (w > max_wgt)
            except Exception:
                break
    if not bad_mask.any():
        return df
    return df.loc[~bad_mask]


def drop_bad_g4_knob_weights(mc_evt_df, knobs=None, max_wgt=1e3, n_univ=100):
    """Drop rows where any listed ``(mc, knob, univ_i)`` weight exceeds ``max_wgt``."""
    df = mc_evt_df
    if knobs is None:
        knobs = g4_mc_knob_names()
    if not knobs:
        return df
    bad_mask = pd.Series(False, index=df.index)
    for knob in knobs:
        sk = ("mc", knob)
        for i in range(int(n_univ)):
            probe = sk + ("univ_{}".format(i),)
            key = multicol_resolve_column_key(df, probe)
            if key is None:
                break
            try:
                w = df.loc[:, key]
                bad_mask = bad_mask | (w > max_wgt)
            except Exception:
                break
    if not bad_mask.any():
        return df
    return df.loc[~bad_mask]


def knob_nested_syst_block(block) -> bool:
    """True when Flux/G4 chunk payload is ``{knob_name: {var_save_name: pack}}``."""
    if not isinstance(block, dict) or not block:
        return False
    v0 = next(iter(block.values()))
    return isinstance(v0, dict) and "univ_events" not in v0


# Backward-compatible aliases (same shape for G4 and Flux knob modes)
g4_block_is_nested_knobs = knob_nested_syst_block
flux_block_is_nested_knobs = knob_nested_syst_block


def count_merged_knob_nested_var_slots(block) -> int:
    """Number of (knob × variable) histogram slots under a nested Flux/G4 merge block."""
    if not block:
        return 0
    if knob_nested_syst_block(block):
        return sum(len(kb) for kb in block.values())
    return len(block)


count_merged_g4_var_slots = count_merged_knob_nested_var_slots


def combine_indep_knob_cov_packs(packs: list, cv_events: np.ndarray) -> dict:
    """Treat per-knob covariances as independent: sum absolute ``cov``, rebuild ``cov_frac`` / ``corr``."""
    if not packs:
        raise ValueError("combine_indep_knob_cov_packs: empty packs")
    cov = np.sum([np.asarray(p["cov"], dtype=float) for p in packs], axis=0)
    mu = np.asarray(cv_events, dtype=float).reshape(-1)
    safe = np.outer(np.maximum(mu, 1e-18), np.maximum(mu, 1e-18))
    cov_frac = np.divide(cov, safe, out=np.zeros_like(cov), where=safe > 0)
    d = np.sqrt(np.maximum(np.diag(cov), 0.0))
    outer = np.outer(np.maximum(d, 1e-18), np.maximum(d, 1e-18))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = np.where(outer > 0, cov / outer, 0.0)
    np.fill_diagonal(corr, 1.0)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    return {"cov": cov, "cov_frac": cov_frac, "corr": corr}


combine_indep_g4_knob_cov_packs = combine_indep_knob_cov_packs


def syst_acc_bucket_nonempty(category: str, block) -> bool:
    """Whether a systematic accumulator has any data (flat ``{var: pack}`` or knob-nested)."""
    if not block:
        return False
    if category in ("Flux", "G4") and knob_nested_syst_block(block):
        return any(bool(inner) for inner in block.values())
    return True


def build_var_configs(var_set: str):
    final_extra = [
        VariableConfig.vertex_x(),
        VariableConfig.vertex_y(),
        VariableConfig.vertex_z(),
        VariableConfig.muon_direction_x(),
        VariableConfig.muon_direction_y(),
        VariableConfig.proton_direction_x(),
        VariableConfig.proton_direction_y(),
        VariableConfig.opening_angle(),
    ]
    final_list = with_final_selected_evt_variables(list(CORE_SELECTED_EVT_VARIABLE_CONFIGS) + final_extra)
    if var_set == "final":
        return final_list
    if var_set == "intermediate":
        return list(INTERMEDIATE_CUT_SYST_VARIABLE_CONFIGS)
    if var_set == "sel_all":
        # Cut-stage variables produced by the pipeline-walker chunk path PLUS the
        # final-selected variables (deduped). Avoids circular import by lazy-loading
        # ``syst_pipeline_walker``.
        from analysis_village.numucc_1p0pi.syst_pipeline_walker import (
            CUT_STAGE_VAR_SPECS,
        )
        seen: set[str] = set()
        merged = []
        for spec in CUT_STAGE_VAR_SPECS:
            vc = spec.var_config
            if vc.var_save_name in seen:
                continue
            seen.add(vc.var_save_name)
            merged.append(vc)
        for vc in final_list:
            if vc.var_save_name in seen:
                continue
            seen.add(vc.var_save_name)
            merged.append(vc)
        return merged
    seen = set()
    merged = []
    for vc in list(INTERMEDIATE_CUT_SYST_VARIABLE_CONFIGS) + final_list:
        if vc.var_save_name in seen:
            continue
        seen.add(vc.var_save_name)
        merged.append(vc)
    return merged


NEUTRINO_SYST_ORDER = ("MCstat", "Flux", "G4")


def syst_key_for_name(sname: str):
    return ("mc", sname) if sname in ("Flux", "G4") else sname


def legacy_npz_wrap(inner_key: str, ret_by_var: dict):
    """Layout expected by ``selected_events.get_syst_unc`` / ``utils.get_syst_unc``."""
    return {vn: np.array({inner_key: ret}, dtype=object) for vn, ret in ret_by_var.items()}


def save_neutrino_multisim_npzs(syst_dict: dict, syst_disk_root: str) -> None:
    """Write MCstat / Flux / G4 NPZs under ``<syst_disk_root>/<MCstat|Flux|G4>/``."""
    root = normalized_root(syst_disk_root)
    spec = (
        ("MCstat", SUB_MCSTAT, FILE_MCSTAT, "MCstat"),
        ("Flux", SUB_FLUX, FILE_FLUX, "flux"),
        ("G4", SUB_G4, FILE_G4, "G4"),
    )
    for dict_key, subdir, fname, inner_key in spec:
        block = syst_dict.get(dict_key)
        if not block:
            continue
        d = os.path.join(root, subdir)
        os.makedirs(d, exist_ok=True)
        payload = legacy_npz_wrap(inner_key, block)
        np.savez_compressed(os.path.join(d, fname), **payload)


def save_cosmics_legacy_npz(syst_dict: dict, syst_disk_root: str) -> None:
    """Write ``Cosmics/cosmics_syst_dict.npz`` under ``syst_disk_root`` if the block is present."""
    block = syst_dict.get("cosmics")
    if not block:
        return
    root = normalized_root(syst_disk_root)
    d = os.path.join(root, SUB_COSMICS)
    os.makedirs(d, exist_ok=True)
    payload = legacy_npz_wrap("Cosmics", block)
    np.savez_compressed(os.path.join(d, FILE_COSMICS), **payload)


def save_legacy_category_npzs(syst_dict: dict, syst_disk_root: str) -> None:
    """Write multisim + cosmics NPZs into a ``syst_disk_layout`` tree."""
    save_neutrino_multisim_npzs(syst_dict, syst_disk_root)
    save_cosmics_legacy_npz(syst_dict, syst_disk_root)
