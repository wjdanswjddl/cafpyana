"""Shared variable registry for cosmics chunk + aggregate + ``get_systematics_cosmics``."""
from __future__ import annotations

from typing import Any, List, Sequence

from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig


def build_variable_configs(arg_vars: Sequence[str] | None) -> List[Any]:
    registry = {
        "integrated": VariableConfig.all_events,
        "vertex_x": VariableConfig.vertex_x,
        "vertex_y": VariableConfig.vertex_y,
        "vertex_z": VariableConfig.vertex_z,
        "muon-p": VariableConfig.muon_momentum,
        "muon-dir_z": VariableConfig.muon_direction,
        "muon-dir_x": VariableConfig.muon_direction_x,
        "muon-dir_y": VariableConfig.muon_direction_y,
        "muon-dir_phi": VariableConfig.muon_direction_phi,
        "proton-p": VariableConfig.proton_momentum,
        "proton-dir_z": VariableConfig.proton_direction,
        "proton-dir_x": VariableConfig.proton_direction_x,
        "proton-dir_y": VariableConfig.proton_direction_y,
        "proton-dir_phi": VariableConfig.proton_direction_phi,
        "muon-end_x": VariableConfig.muon_end_x,
        "muon-end_y": VariableConfig.muon_end_y,
        "muon-end_z": VariableConfig.muon_end_z,
        "opening_angle": VariableConfig.opening_angle,
        "tki-del_alpha": VariableConfig.tki_del_alpha,
        "tki-del_phi": VariableConfig.tki_del_phi,
        "tki-del_Tp": VariableConfig.tki_del_Tp,
        "tki-del_p": VariableConfig.tki_del_p,
        "tki-del_Tp_x": VariableConfig.tki_del_Tp_x,
        "tki-del_Tp_y": VariableConfig.tki_del_Tp_y,
    }
    if arg_vars:
        out = []
        for name in arg_vars:
            key = name.strip()
            if key not in registry:
                raise ValueError("Unknown variable key '%s'. Choices: %s" % (key, sorted(registry)))
            out.append(registry[key]())
        return out
    return with_final_selected_evt_variables(
        list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)
        + [
            VariableConfig.vertex_x(),
            VariableConfig.vertex_y(),
            VariableConfig.vertex_z(),
            VariableConfig.muon_direction_x(),
            VariableConfig.muon_direction_y(),
            VariableConfig.proton_direction_x(),
            VariableConfig.proton_direction_y(),
            VariableConfig.opening_angle(),
        ]
    )
