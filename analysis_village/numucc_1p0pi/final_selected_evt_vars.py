"""Canonical kinematic VariableConfig sets for final-sample plots/systematics.

* :data:`CORE_SELECTED_EVT_VARIABLE_CONFIGS` — baseline distributions (integrated + μ/p
  kinematics + TKI) used across notebooks and syst scripts.
* :data:`FINAL_SELECTED_EVT_VARIABLE_CONFIGS` — extra φ / vertex / endpoint components;
  merged via :func:`with_final_selected_evt_variables` without duplicating ``var_save_name``.
"""

from __future__ import annotations

from typing import List, Sequence

from analysis_village.numucc_1p0pi.variable_configs import VariableConfig

__all__ = (
    "CORE_SELECTED_EVT_VARIABLE_CONFIGS",
    "FINAL_SELECTED_EVT_VARIABLE_CONFIGS",
    "INTERMEDIATE_CUT_SYST_VARIABLE_CONFIGS",
    "with_final_selected_evt_variables",
)

# Variables defined on loose / pre-final evt dfs (nu score, multiplicity, vertex).
# Use with MC dfs from ``get_ana_dfs("systs", systs_mc_df_tag="-sel_all-wgts", ...)``.
INTERMEDIATE_CUT_SYST_VARIABLE_CONFIGS: tuple[VariableConfig, ...] = (
    VariableConfig.all_events(),
    VariableConfig.nu_score(),
    VariableConfig.n_trks(),
    VariableConfig.vertex_x(),
    VariableConfig.vertex_y(),
    VariableConfig.vertex_z(),
)

CORE_SELECTED_EVT_VARIABLE_CONFIGS: tuple[VariableConfig, ...] = (
    VariableConfig.all_events(),
    VariableConfig.muon_momentum(),
    VariableConfig.muon_direction(),
    VariableConfig.proton_momentum(),
    VariableConfig.proton_direction(),
    VariableConfig.tki_del_Tp(),
    VariableConfig.tki_del_Tp_x(),
    VariableConfig.tki_del_Tp_y(),
    VariableConfig.tki_del_p(),
    VariableConfig.tki_del_alpha(),
    VariableConfig.tki_del_phi(),
)

FINAL_SELECTED_EVT_VARIABLE_CONFIGS: tuple[VariableConfig, ...] = (
    VariableConfig.muon_direction_phi(),
    VariableConfig.proton_direction_phi(),
    VariableConfig.muon_direction_x(),
    VariableConfig.muon_direction_y(),
    VariableConfig.proton_direction_x(),
    VariableConfig.proton_direction_y(),
    VariableConfig.vertex_x(),
    VariableConfig.vertex_y(),
    VariableConfig.vertex_z(),
    VariableConfig.muon_end_x(),
    VariableConfig.muon_end_y(),
    VariableConfig.muon_end_z(),
)


def with_final_selected_evt_variables(configs: Sequence[VariableConfig]) -> List[VariableConfig]:
    """``configs`` plus any entry from :data:`FINAL_SELECTED_EVT_VARIABLE_CONFIGS` not already present."""
    seen = {c.var_save_name for c in configs}
    out: List[VariableConfig] = list(configs)
    for vc in FINAL_SELECTED_EVT_VARIABLE_CONFIGS:
        if vc.var_save_name not in seen:
            out.append(vc)
            seen.add(vc.var_save_name)
    return out
