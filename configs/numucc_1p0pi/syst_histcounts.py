# Systematic histogram counts on the grid (BNB / dirt MC with weights, or unisim samples).
#
# Walks the numuCC 1p0pi event selection on each CAF and writes long-format bin
# counts under HDF key ``syst_hists`` (+ ``hdr`` for POT/gates). No heavy weight
# event tables — see ``analysis_village.numucc_1p0pi.syst_histcounts``.
#
# Modes via env ``SYST_HIST_MODE`` (default: genie):
#   all         — all GENIE knobs (incl. Ar23p) + Flux + G4 in one CAF pass
#                 (stored per knob, e.g. CoulombCCQE / NormCCMEC — not per mode)
#                 + slim products: slim_multisim (true multisim only) and slim
#                   (× Gaussian throws of multisigma/morph; weights ≥ 0)
#   genie       — GENIE_KNOB_GROUP=<CCQE|MEC|…> selects which knobs to load;
#                 unset = all groups. Output still one hist stack per knob
#                 (+ slim_multisim + slim by default)
#   genie_slim  — getsyst slim load; histcounts still build slim_multisim + slim
#   flux        — BNB flux multisim knobs (rate only) + Flux_slim_* by default
#   g4          — Geant4 reinteraction knobs (rate only) + G4_slim_* by default
#   nowgt       — no weights; CV counts only (WireMod, DENT, intime, offbeam)
#
# Optional: ``SYST_HIST_SAMPLE=mc|dirt|intime|offbeam|data``
# Optional: ``SYST_HIST_EXCLUDE_AR23P=1`` to drop Ar23p knobs from all/genie.
# Optional: ``SYST_HIST_EXCLUDE_SLIM=1`` to skip slim product histcounts.
#
# Notebook: do not sum slim / slim_multisim frac-cov with the sum of per-knob
# multisim frac-covs (double-counting). Prefer ``slim`` as the combined envelope
# (already includes slim_multisim + thrown ±σ/morph); or use per-knob sum.
#
# VariableConfig freeze: each job HDF includes key ``var_configs``; grid submit
# also writes ``variable_configs.json`` next to outputs. Notebooks must load that
# snapshot for binning/labels (do not import live ``variable_configs.py``).
## Examples:
#   SYST_HIST_MODE=all python run_df_maker.py \
#     -c configs/numucc_1p0pi/syst_histcounts.py -l mc.list -o hist_all -ngrid 500
#
#   SYST_HIST_MODE=genie GENIE_KNOB_GROUP=CCQE python run_df_maker.py \
#     -c configs/numucc_1p0pi/syst_histcounts.py -l mc.list -o hist_genie_CCQE -ngrid 500
#
#   SYST_HIST_MODE=flux python run_df_maker.py \
#     -c configs/numucc_1p0pi/syst_histcounts.py -l mc.list -o hist_flux -ngrid 500
#
#   SYST_HIST_MODE=nowgt SYST_HIST_SAMPLE=offbeam python run_df_maker.py \
#     -c configs/numucc_1p0pi/syst_histcounts.py -l offbeam.list -o hist_offbeam -ngrid 200
import os

from analysis_village.numucc_1p0pi.makedf.makedf import build_syst_histcounts_config

_mode = os.environ.get("SYST_HIST_MODE", "genie").strip().lower()
_sample = os.environ.get("SYST_HIST_SAMPLE", "mc").strip().lower() or "mc"
_group = os.environ.get("GENIE_KNOB_GROUP", "").strip() or None

DFS, ARGS, NAMES = build_syst_histcounts_config(
    mode=_mode,
    group_filter=_group,
    sample=_sample,
)
