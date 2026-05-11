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
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    INTERMEDIATE_CUT_SYST_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)


def drop_bad_g4_weights(mc_evt_df, max_wgt=1e3, n_univ=100):
    """Drop evt rows with any G4 universe weight above ``max_wgt`` (no-op if columns missing)."""
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
