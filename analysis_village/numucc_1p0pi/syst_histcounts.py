"""Grid-side systematic histogram counts (selection + xsec variables).

Production path
---------------
CAF → :func:`fill_syst_histcounts` (walk event selection, fill per-knob /
per-universe bins) → long-format ``DataFrame`` stored as HDF key ``syst_hists``
by ``run_df_maker`` → notebook loads + sums across files → covariance.

Each job also stores a frozen VariableConfig table under HDF key ``var_configs``
(and grid submit writes ``variable_configs.json`` next to outputs). Plotting must
read that snapshot for bin edges/labels — not the live ``variable_configs.py``.

Sparse packing always emits a trailing-bin shape sentinel so unpack recovers the
true ``nbins`` even when the last bins are empty.

GENIE **rate** vs **xsec** (critical)
-------------------------------------
* **rate** — reweight the selected reco spectrum (signal + Δbackground). Stored as
  ``rate_cv`` / ``rate_univ`` rows.
* **xsec** — does **not** reweight selected yield the same way. We store additive
  response-matrix tensors (same recipe as ``scripts/get_systematics_genie.py``):

  - ``xsec_nevts_allmc`` — CV generated signal (truth bins)
  - ``xsec_cv_sel_reco`` / ``xsec_cv_allsel_reco`` — CV selected reco (bkg-sub / all)
  - ``xsec_signal_sel_truth[u]``, ``xsec_signal_allmc[u]`` — for
    ``ε_u = sel_truth(u) / allmc(u)``
  - ``xsec_reco_vs_true[u]`` — migration ``histogram2d(truth, reco)`` with univ weights
  - ``xsec_bg_cv``, ``xsec_bg_univ[u]`` — background reco yields

  Universe xsec spectrum (notebook / :func:`finalize_genie_xsec_univ`)::

      R_u = response(reco_vs_true_u, ε_u)
      N_u = R_u @ N_gen^CV  +  (bg_u - bg_CV)

  so pure rate/normalization knobs that cancel between selected and generated
  are suppressed on the xsec path.

Flux / G4 are **rate-only** multisims (same fractional cov enters rate and xsec
totals downstream). Unisim samples (WireMod, DENT, intime, offbeam) store only
``rate_cv`` (no universe weights on the CAF).

"""
from __future__ import annotations

import gc
import json
import os
from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from pyanalib.covariance import get_covariance_matrix
from pyanalib.variable_calculator import (
    add_mc_cc1p0pi_tki_mcnu,
    add_reco_cc1p0pi_tki_evtdf,
    add_truth_cc1p0pi_tki_evtdf,
)

from analysis_village.numucc_1p0pi.categories import get_genie_category, get_topo_category
from analysis_village.numucc_1p0pi.evt_derived_kinematics import (
    ensure_derived_trk_kinematics_cols,
    ensure_mc_level_phi_mcnu,
)
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.syst_pipeline_walker import (
    CUT_STAGE_RATE_ONLY_SLUGS,
    CUT_STAGE_VAR_SPECS,
    FINAL_STAGE_KEY,
    get_var_series,
    histogram_var,
    walk_pipeline,
)
from analysis_village.numucc_1p0pi.utils import (
    genie_univ_weight_series,
    get_clipped_evts,
    get_response_matrix,
    get_univ_rates,
    signal_hists,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig

# Long-format column schema for HDF ``syst_hists`` tables.
HIST_DF_COLUMNS: Tuple[str, ...] = (
    "kind",
    "family",
    "knob",
    "var",
    "univ",
    "bin0",
    "bin1",
    "value",
)

SystName = Tuple[str, str]

# kind tags
KIND_RATE_CV = "rate_cv"
KIND_RATE_UNIV = "rate_univ"
KIND_XSEC_NEVTS_ALLMC = "xsec_nevts_allmc"
KIND_XSEC_CV_SEL_RECO = "xsec_cv_sel_reco"
KIND_XSEC_CV_ALLSEL_RECO = "xsec_cv_allsel_reco"
KIND_XSEC_BG_CV = "xsec_bg_cv"
KIND_XSEC_SIGNAL_ALLMC = "xsec_signal_allmc"
KIND_XSEC_SIGNAL_SEL_TRUTH = "xsec_signal_sel_truth"
KIND_XSEC_BG_UNIV = "xsec_bg_univ"
KIND_XSEC_RECO_VS_TRUE = "xsec_reco_vs_true"


def empty_histcounts_df() -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series(dtype=object if c in ("kind", "family", "knob", "var") else float)
                         for c in HIST_DF_COLUMNS})


def final_var_configs() -> List[VariableConfig]:
    return with_final_selected_evt_variables(list(CORE_SELECTED_EVT_VARIABLE_CONFIGS))


# ---------------------------------------------------------------------------
# Frozen VariableConfig snapshot (written with job outputs; used for plotting)
# ---------------------------------------------------------------------------
VAR_CONFIG_SNAPSHOT_SCHEMA = 1
VAR_CONFIG_SNAPSHOT_JSON_NAME = "variable_configs.json"


@dataclass
class FrozenVarConfig:
    """Plotting-only VariableConfig freeze (bins + labels); no live module import."""

    var_save_name: str
    var_plot_name: str
    var_labels: List[str]
    bins: np.ndarray
    xsec_label: str = ""
    category_syst_var_save_name: Optional[str] = None
    var_evt_reco_col: Optional[Tuple[Any, ...]] = None
    var_evt_truth_col: Optional[Tuple[Any, ...]] = None
    var_nu_col: Optional[Tuple[Any, ...]] = None

    def __post_init__(self) -> None:
        self.bins = np.asarray(self.bins, dtype=float)
        self.var_labels = [str(x) for x in self.var_labels]
        self.bin_centers = (self.bins[:-1] + self.bins[1:]) / 2.0

    @property
    def n_bins(self) -> int:
        return max(int(len(self.bins) - 1), 0)


def histcounts_var_configs() -> List[VariableConfig]:
    """All VariableConfigs used when filling histcounts (cut-stage + final)."""
    by_name: Dict[str, VariableConfig] = {}
    for spec in CUT_STAGE_VAR_SPECS:
        by_name[spec.var_config.var_save_name] = spec.var_config
    for vc in final_var_configs():
        by_name[vc.var_save_name] = vc
    return list(by_name.values())


def _col_to_jsonable(col: Any) -> Any:
    if col is None:
        return None
    if isinstance(col, tuple):
        return list(col)
    if isinstance(col, list):
        return list(col)
    return col


def var_config_to_record(vc: Any) -> Dict[str, Any]:
    return {
        "var_save_name": str(vc.var_save_name),
        "var_plot_name": str(getattr(vc, "var_plot_name", vc.var_save_name)),
        "var_labels": [str(x) for x in getattr(vc, "var_labels", [])],
        "bins": np.asarray(vc.bins, dtype=float).tolist(),
        "xsec_label": str(getattr(vc, "xsec_label", "") or ""),
        "category_syst_var_save_name": getattr(vc, "category_syst_var_save_name", None),
        "var_evt_reco_col": _col_to_jsonable(getattr(vc, "var_evt_reco_col", None)),
        "var_evt_truth_col": _col_to_jsonable(getattr(vc, "var_evt_truth_col", None)),
        "var_nu_col": _col_to_jsonable(getattr(vc, "var_nu_col", None)),
    }


def record_to_frozen_var_config(rec: Mapping[str, Any]) -> FrozenVarConfig:
    def _tup(x: Any) -> Optional[Tuple[Any, ...]]:
        if x is None:
            return None
        if isinstance(x, (list, tuple)):
            return tuple(x)
        return (x,)

    return FrozenVarConfig(
        var_save_name=str(rec["var_save_name"]),
        var_plot_name=str(rec.get("var_plot_name", rec["var_save_name"])),
        var_labels=list(rec.get("var_labels") or []),
        bins=np.asarray(rec["bins"], dtype=float),
        xsec_label=str(rec.get("xsec_label") or ""),
        category_syst_var_save_name=rec.get("category_syst_var_save_name"),
        var_evt_reco_col=_tup(rec.get("var_evt_reco_col")),
        var_evt_truth_col=_tup(rec.get("var_evt_truth_col")),
        var_nu_col=_tup(rec.get("var_nu_col")),
    )


def var_configs_to_snapshot_dict(configs: Sequence[Any]) -> Dict[str, Any]:
    return {
        "schema_version": VAR_CONFIG_SNAPSHOT_SCHEMA,
        "variables": [var_config_to_record(vc) for vc in configs],
    }


def frozen_var_configs_from_snapshot(
    snapshot: Mapping[str, Any],
) -> Dict[str, FrozenVarConfig]:
    out: Dict[str, FrozenVarConfig] = {}
    for rec in snapshot.get("variables") or []:
        fc = record_to_frozen_var_config(rec)
        out[fc.var_save_name] = fc
    return out


def nbins_by_var_from_configs(configs: Mapping[str, Any]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for name, vc in configs.items():
        bins = getattr(vc, "bins", None)
        if bins is None:
            continue
        out[str(name)] = max(int(len(np.asarray(bins)) - 1), 0)
    return out


def histcounts_var_configs_df(configs: Optional[Sequence[Any]] = None) -> pd.DataFrame:
    """Serialize VariableConfigs to a one-row-per-var DataFrame (HDF ``var_configs``)."""
    if configs is None:
        configs = histcounts_var_configs()
    rows = []
    for vc in configs:
        rec = var_config_to_record(vc)
        rows.append(
            {
                "var_save_name": rec["var_save_name"],
                "var_plot_name": rec["var_plot_name"],
                "xsec_label": rec["xsec_label"],
                "category_syst_var_save_name": rec["category_syst_var_save_name"],
                "bins_json": json.dumps(rec["bins"]),
                "var_labels_json": json.dumps(rec["var_labels"]),
                "var_evt_reco_col_json": json.dumps(rec["var_evt_reco_col"]),
                "var_evt_truth_col_json": json.dumps(rec["var_evt_truth_col"]),
                "var_nu_col_json": json.dumps(rec["var_nu_col"]),
                "n_bins": max(len(rec["bins"]) - 1, 0),
            }
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "var_save_name",
                "var_plot_name",
                "xsec_label",
                "category_syst_var_save_name",
                "bins_json",
                "var_labels_json",
                "var_evt_reco_col_json",
                "var_evt_truth_col_json",
                "var_nu_col_json",
                "n_bins",
            ]
        )
    return pd.DataFrame(rows)


def frozen_var_configs_from_df(df: pd.DataFrame) -> Dict[str, FrozenVarConfig]:
    if df is None or len(df) == 0:
        return {}
    out: Dict[str, FrozenVarConfig] = {}
    for _, row in df.iterrows():
        rec = {
            "var_save_name": row["var_save_name"],
            "var_plot_name": row.get("var_plot_name", row["var_save_name"]),
            "xsec_label": row.get("xsec_label", ""),
            "category_syst_var_save_name": row.get("category_syst_var_save_name"),
            "bins": json.loads(row["bins_json"]),
            "var_labels": json.loads(row.get("var_labels_json") or "[]"),
            "var_evt_reco_col": json.loads(row.get("var_evt_reco_col_json") or "null"),
            "var_evt_truth_col": json.loads(row.get("var_evt_truth_col_json") or "null"),
            "var_nu_col": json.loads(row.get("var_nu_col_json") or "null"),
        }
        fc = record_to_frozen_var_config(rec)
        out[fc.var_save_name] = fc
    return out


def write_var_config_snapshot_json(
    path: str,
    configs: Optional[Sequence[Any]] = None,
) -> str:
    """Write frozen VariableConfig JSON (call at job submission / next to outputs)."""
    if configs is None:
        configs = histcounts_var_configs()
    path = os.path.abspath(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(var_configs_to_snapshot_dict(configs), fh, indent=2, sort_keys=True)
        fh.write("\n")
    return path


def load_var_config_snapshot_json(path: str) -> Dict[str, FrozenVarConfig]:
    with open(path, "r", encoding="utf-8") as fh:
        return frozen_var_configs_from_snapshot(json.load(fh))


def load_var_configs_from_df_file(path: str) -> Dict[str, FrozenVarConfig]:
    """Read ``var_configs_*`` from one ``run_df_maker`` HDF, else sidecar JSON."""
    with pd.HDFStore(path, mode="r") as store:
        keys = [k.lstrip("/") for k in store.keys()]
        vc_keys = sorted(k for k in keys if k.startswith("var_configs"))
        for k in vc_keys:
            cfg = frozen_var_configs_from_df(store[k])
            if cfg:
                return cfg
    sidecar = os.path.join(os.path.dirname(path), VAR_CONFIG_SNAPSHOT_JSON_NAME)
    if os.path.isfile(sidecar):
        return load_var_config_snapshot_json(sidecar)
    return {}


def load_var_configs_from_glob(paths: Sequence[str]) -> Dict[str, FrozenVarConfig]:
    """First non-empty snapshot among histcount files (or their output-dir JSON)."""
    for p in paths:
        try:
            cfg = load_var_configs_from_df_file(p)
        except Exception:
            cfg = {}
        if cfg:
            return cfg
    seen_dirs = set()
    for p in paths:
        d = os.path.dirname(os.path.abspath(p))
        if d in seen_dirs:
            continue
        seen_dirs.add(d)
        sidecar = os.path.join(d, VAR_CONFIG_SNAPSHOT_JSON_NAME)
        if os.path.isfile(sidecar):
            try:
                return load_var_config_snapshot_json(sidecar)
            except Exception:
                pass
    return {}


def _prefix_mcnu_columns(mc_nu_df: pd.DataFrame) -> pd.DataFrame:
    if mc_nu_df is None or len(mc_nu_df) == 0:
        return mc_nu_df
    if not isinstance(mc_nu_df.columns, pd.MultiIndex):
        return mc_nu_df
    try:
        first_level = mc_nu_df.columns.get_level_values(0)
        need_prefix = not np.all(first_level == "mc")
    except Exception:
        need_prefix = True
    if need_prefix:
        out = mc_nu_df.copy()
        out.columns = pd.MultiIndex.from_tuples(
            [tuple(["mc"] + list(c)) for c in mc_nu_df.columns]
        )
        return out
    return mc_nu_df


def _annotate_topo(evt_df: pd.DataFrame, mc_nu_df: Optional[pd.DataFrame]) -> None:
    if evt_df is not None and len(evt_df) > 0 and "topo_categ" not in evt_df.columns:
        evt_df.loc[:, "topo_categ"] = get_topo_category(evt_df)
        try:
            evt_df.loc[:, "genie_categ"] = get_genie_category(evt_df)
        except Exception:
            pass
    if mc_nu_df is not None and len(mc_nu_df) > 0 and "topo_categ" not in mc_nu_df.columns:
        mc_nu_df.loc[:, "topo_categ"] = get_topo_category(mc_nu_df)
        try:
            mc_nu_df.loc[:, "genie_categ"] = get_genie_category(mc_nu_df)
        except Exception:
            pass


def _align_evt_mcnu(evt_df: pd.DataFrame, mc_nu_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if evt_df is None or len(evt_df) == 0:
        empty = mc_nu_df.iloc[0:0] if mc_nu_df is not None else evt_df
        return evt_df, empty
    if mc_nu_df is None or len(mc_nu_df) == 0:
        return evt_df, mc_nu_df
    ix = evt_df.index.intersection(mc_nu_df.index)
    return evt_df.loc[ix], mc_nu_df.loc[ix]


def discover_syst_names_on_df(
    evt_df: pd.DataFrame,
    *,
    family: str,
    knob_names: Optional[Sequence[str]] = None,
) -> List[SystName]:
    """Return ``(mc, <knob>)`` blocks present on ``evt_df``.

    Keys are **individual CAF weight knobs** (e.g. ``..._CoulombCCQE``,
    ``..._NormCCMEC``), never mode/group labels (``CCQE``, ``MEC``, …).

    * If ``knob_names`` is given, keep those present under ``mc`` (may include
      ``slim`` / ``slim_multisim`` / ``Flux_slim`` / …).
    * Else auto-detect knobs under ``mc`` that have ``univ_*`` / ``ps1`` / ``morph``.
      Auto-detect returns per-knob blocks; slim product names are included when present.
    """
    if evt_df is None or len(evt_df) == 0 or not isinstance(evt_df.columns, pd.MultiIndex):
        return []
    try:
        mc = evt_df.mc
    except (AttributeError, KeyError):
        return []
    available = set(mc.columns.get_level_values(0).unique())
    if knob_names is not None:
        return [("mc", str(k)) for k in knob_names if str(k) in available]

    out: List[SystName] = []
    for knob in sorted(available):
        if knob in (None, ""):
            continue
        try:
            leaves = {
                str(c[0]) if isinstance(c, tuple) else str(c) for c in mc[knob].columns
            }
        except Exception:
            continue
        if leaves & {"ps1", "ms1", "morph"} or any(x.startswith("univ_") for x in leaves):
            out.append(("mc", str(knob)))
    return out


def _mc_knob_leaf_set(df: pd.DataFrame, knob: str) -> set:
    try:
        block = df.mc[knob]
    except Exception:
        return set()
    leaves = set()
    cols = block.columns
    if isinstance(cols, pd.MultiIndex):
        for c in cols:
            if isinstance(c, tuple):
                for part in c:
                    if part:
                        leaves.add(str(part))
                        break
            else:
                leaves.add(str(c))
    else:
        leaves = {str(c) for c in cols}
    return leaves


def _is_true_multisim_knob(df: pd.DataFrame, knob: str) -> bool:
    """True if knob has ``univ_*`` leaves (CAF type-0 multisim), not only ±σ/morph."""
    leaves = _mc_knob_leaf_set(df, knob)
    return any(x.startswith("univ_") for x in leaves)


def _is_multisigma_or_morph_knob(df: pd.DataFrame, knob: str) -> bool:
    leaves = _mc_knob_leaf_set(df, knob)
    return bool(leaves & {"ps1", "ms1", "morph"}) and not any(
        x.startswith("univ_") for x in leaves
    )


def _resolve_mc_weight_col(df: pd.DataFrame, knob: str, leaf: str):
    """Return full MultiIndex column key for ``mc.<knob>.<leaf>`` if present."""
    from analysis_village.numucc_1p0pi.selection_framework import multicol_resolve_column_key

    probe = ("mc", knob, leaf)
    return multicol_resolve_column_key(df, probe)


def slim_product_names(family: str) -> Tuple[str, str]:
    """Return ``(slim_multisim_name, slim_name)`` for a weight family.

    GENIE uses the short names ``slim_multisim`` / ``slim`` (user-facing).
    Flux/G4 are prefixed so they do not collide on a shared ``mode=all`` frame.
    """
    fam = str(family).upper()
    if fam == "GENIE":
        return "slim_multisim", "slim"
    if fam == "FLUX":
        return "Flux_slim_multisim", "Flux_slim"
    if fam == "G4":
        return "G4_slim_multisim", "G4_slim"
    return f"{family}_slim_multisim", f"{family}_slim"


def multisim_throw_seed(tag: str, knob: str, univ_i: int) -> int:
    """Stable 32-bit seed for a (product, knob, universe) Gaussian throw."""
    import zlib

    payload = f"{tag}|{knob}|{int(univ_i)}".encode("utf-8")
    return int(zlib.crc32(payload) & 0xFFFFFFFF)


def _dst_weight_col(sample_key: tuple, product_name: str, leaf: str) -> tuple:
    nlevels = len(sample_key) if isinstance(sample_key, tuple) else 2
    dst = tuple(["mc", product_name, leaf] + [""] * max(0, nlevels - 3))
    if len(dst) != nlevels:
        dst = tuple(list(sample_key[:1]) + [product_name, leaf] + [""] * max(0, nlevels - 3))
        if len(dst) < nlevels:
            dst = dst + tuple([""] * (nlevels - len(dst)))
        elif len(dst) > nlevels:
            dst = dst[:nlevels]
    return dst


def _clip_physical_wgt(w: np.ndarray, *, hi: Optional[float] = 10.0) -> np.ndarray:
    w = np.nan_to_num(np.asarray(w, dtype=np.float64), nan=1.0, posinf=1.0, neginf=0.0)
    if hi is None:
        return np.clip(w, 0.0, None)
    return np.clip(w, 0.0, float(hi))


def _throw_multisigma_factor(ps1: np.ndarray, z: float) -> np.ndarray:
    """Scalar-``z`` multisigma throw: ``1 + (ps1 - 1) * z``, clipped ``≥ 0``.

    Matches ``multisigma_to_multisim.ipynb`` / historical ``getsyst`` slim
    (``isigma == 0`` uses ``ps1`` only; one ``N(0,1)`` draw per universe, shared
    by all events — not per-event).
    """
    w = 1.0 + (np.asarray(ps1, dtype=np.float64) - 1.0) * float(z)
    return _clip_physical_wgt(w, hi=None)


def _throw_morph_factor(morph: np.ndarray, z: float) -> np.ndarray:
    """Morph throw: ``1 + (morph - 1) * 2 * |z|``, clipped ``≥ 0``.

    The factor ``2`` matches the notebook / historical ``getsyst`` recipe
    (morph CAF weight is treated as a half-σ template; ``2*|z|`` restores a
    unit-Gaussian envelope). Weights stay non-negative.
    """
    w = 1.0 + (np.asarray(morph, dtype=np.float64) - 1.0) * 2.0 * abs(float(z))
    return _clip_physical_wgt(w, hi=None)


def attach_slim_multisim_product(
    df: pd.DataFrame,
    *,
    product_name: str,
    n_univ: int,
    knob_names: Optional[Sequence[str]] = None,
    wgt_clip_hi: float = 10.0,
    exclude_names: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Attach ``mc.<product_name>.univ_i`` = ∏ true-multisim knobs (CAF type 0).

    Physical rules:
    * Only knobs with ``univ_*`` leaves enter the product.
    * Multisigma (``ps*``/``ms*``) and morph (``morph``) are **never** multiplied in
      here — they have no shared CAF universe index.
    * Every factor is clipped to ``[0, wgt_clip_hi]`` before multiplying.

    No-op if ``product_name`` already has ``univ_0`` on ``df``.
    """
    if df is None or len(df) == 0 or n_univ <= 0:
        return df
    if _resolve_mc_weight_col(df, product_name, "univ_0") is not None:
        return df

    try:
        available = list(df.mc.columns.get_level_values(0).unique())
    except Exception:
        return df

    skip = set(exclude_names or ())
    skip.update({product_name, "slim", "slim_multisim", "GENIE", "Flux", "G4"})
    if knob_names is None:
        candidates = [str(k) for k in available if k not in skip and k not in (None, "")]
    else:
        candidates = [str(k) for k in knob_names if str(k) in available and str(k) not in skip]

    multisim_knobs = [k for k in candidates if _is_true_multisim_knob(df, k)]
    if not multisim_knobs:
        return df

    sample_key = _resolve_mc_weight_col(df, multisim_knobs[0], "univ_0")
    if sample_key is None:
        return df

    out = df
    for u in range(int(n_univ)):
        leaf = f"univ_{u}"
        prod = np.ones(len(out), dtype=np.float64)
        any_factor = False
        for knob in multisim_knobs:
            key = _resolve_mc_weight_col(out, knob, leaf)
            if key is None:
                continue
            w = _clip_physical_wgt(out.loc[:, key], hi=wgt_clip_hi)
            prod *= w
            any_factor = True
        if not any_factor:
            continue
        prod = _clip_physical_wgt(prod, hi=None)
        out.loc[:, _dst_weight_col(sample_key, product_name, leaf)] = prod
    return out


def attach_slim_full_product(
    df: pd.DataFrame,
    *,
    product_name: str,
    n_univ: int,
    knob_names: Optional[Sequence[str]] = None,
    slim_multisim_name: Optional[str] = None,
    seed_tag: str = "slim",
    wgt_clip_hi: float = 10.0,
    exclude_names: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Attach ``mc.<product_name>.univ_i`` = slim_multisim × multisigma/morph throws.

    Recipe (``multisigma_to_multisim.ipynb`` / historical ``getsyst`` slim):
    * Start from the true-multisim product (or 1 if absent).
    * Multisigma knobs (``ps1``): one scalar ``z ~ N(0,1)`` per universe;
      ``wgt = max(0, 1 + (ps1 - 1) * z)``.
    * Morph knobs: one scalar ``z ~ N(0,1)`` per universe;
      ``wgt = max(0, 1 + (morph - 1) * 2 * |z|)``.
    * Only ``ps1`` (1σ) is used for multisigma — not ``ps2``/``ps3``.

    No-op if ``product_name`` already has ``univ_0`` on ``df``.
    """
    if df is None or len(df) == 0 or n_univ <= 0:
        return df
    if _resolve_mc_weight_col(df, product_name, "univ_0") is not None:
        return df

    try:
        available = list(df.mc.columns.get_level_values(0).unique())
    except Exception:
        return df

    skip = set(exclude_names or ())
    skip.update(
        {
            product_name,
            slim_multisim_name or "",
            "slim",
            "slim_multisim",
            "GENIE",
            "Flux",
            "G4",
        }
    )
    if knob_names is None:
        candidates = [str(k) for k in available if k not in skip and k not in (None, "")]
    else:
        candidates = [str(k) for k in knob_names if str(k) in available and str(k) not in skip]

    ms_knobs = [k for k in candidates if _is_multisigma_or_morph_knob(df, k)]
    # Prefer existing slim_multisim product as the multisim base.
    base_name = slim_multisim_name
    if base_name and _resolve_mc_weight_col(df, base_name, "univ_0") is None:
        base_name = None
    if base_name is None:
        # Fall back: product of true multisim knobs inline (same as attach_slim_multisim).
        ms_true = [k for k in candidates if _is_true_multisim_knob(df, k)]
    else:
        ms_true = []

    # Need a sample key for column padding.
    sample_key = None
    if base_name:
        sample_key = _resolve_mc_weight_col(df, base_name, "univ_0")
    if sample_key is None:
        for k in candidates:
            sample_key = _resolve_mc_weight_col(df, k, "univ_0") or _resolve_mc_weight_col(
                df, k, "ps1"
            ) or _resolve_mc_weight_col(df, k, "morph")
            if sample_key is not None:
                break
    if sample_key is None:
        return df

    # Pre-resolve ±σ / morph source columns.
    ms_sources: List[Tuple[str, str, object]] = []  # (knob, kind, colkey)
    for k in ms_knobs:
        leaves = _mc_knob_leaf_set(df, k)
        if "ps1" in leaves:
            key = _resolve_mc_weight_col(df, k, "ps1")
            if key is not None:
                ms_sources.append((k, "multisigma", key))
        elif "morph" in leaves:
            key = _resolve_mc_weight_col(df, k, "morph")
            if key is not None:
                ms_sources.append((k, "morph", key))

    if base_name is None and not ms_true and not ms_sources:
        return df

    out = df
    rng = np.random.RandomState(0)  # re-seeded per throw
    for u in range(int(n_univ)):
        leaf = f"univ_{u}"
        if base_name is not None:
            bkey = _resolve_mc_weight_col(out, base_name, leaf)
            if bkey is None:
                prod = np.ones(len(out), dtype=np.float64)
            else:
                prod = _clip_physical_wgt(out.loc[:, bkey], hi=wgt_clip_hi)
        else:
            prod = np.ones(len(out), dtype=np.float64)
            for knob in ms_true:
                key = _resolve_mc_weight_col(out, knob, leaf)
                if key is None:
                    continue
                prod *= _clip_physical_wgt(out.loc[:, key], hi=wgt_clip_hi)

        for knob, kind, src_key in ms_sources:
            rng.seed(multisim_throw_seed(seed_tag, knob, u))
            z = float(rng.normal(0.0, 1.0))
            src = _clip_physical_wgt(out.loc[:, src_key], hi=None)
            if kind == "multisigma":
                prod *= _throw_multisigma_factor(src, z)
            else:
                prod *= _throw_morph_factor(src, z)

        prod = _clip_physical_wgt(prod, hi=None)
        out.loc[:, _dst_weight_col(sample_key, product_name, leaf)] = prod
    return out


def attach_family_slim_products(
    df: pd.DataFrame,
    *,
    family: str,
    n_univ: int,
    knob_names: Optional[Sequence[str]] = None,
    wgt_clip_hi: float = 10.0,
) -> Tuple[pd.DataFrame, List[str]]:
    """Attach ``slim_multisim`` + ``slim`` products for ``family``. Returns ``(df, names)``."""
    name_ms, name_full = slim_product_names(family)
    exclude = [name_ms, name_full]
    # Legacy bundled name from getsyst(slim=True, slimname=family)
    legacy = {"GENIE": "GENIE", "FLUX": "Flux", "G4": "G4"}.get(str(family).upper())
    if legacy:
        exclude.append(legacy)

    out = attach_slim_multisim_product(
        df,
        product_name=name_ms,
        n_univ=n_univ,
        knob_names=knob_names,
        wgt_clip_hi=wgt_clip_hi,
        exclude_names=exclude,
    )
    # If per-knob multisim was slimmed away already, copy legacy GENIE/Flux/G4 → slim_multisim.
    if (
        _resolve_mc_weight_col(out, name_ms, "univ_0") is None
        and legacy
        and _resolve_mc_weight_col(out, legacy, "univ_0") is not None
    ):
        sample = _resolve_mc_weight_col(out, legacy, "univ_0")
        for u in range(int(n_univ)):
            leaf = f"univ_{u}"
            src = _resolve_mc_weight_col(out, legacy, leaf)
            if src is None or sample is None:
                continue
            out.loc[:, _dst_weight_col(sample, name_ms, leaf)] = _clip_physical_wgt(
                out.loc[:, src], hi=wgt_clip_hi
            )

    out = attach_slim_full_product(
        out,
        product_name=name_full,
        n_univ=n_univ,
        knob_names=knob_names,
        slim_multisim_name=name_ms,
        seed_tag=name_full,
        wgt_clip_hi=wgt_clip_hi,
        exclude_names=exclude,
    )
    attached = []
    if _resolve_mc_weight_col(out, name_ms, "univ_0") is not None:
        attached.append(name_ms)
    if _resolve_mc_weight_col(out, name_full, "univ_0") is not None:
        attached.append(name_full)
    return out, attached


def clip_nonnegative_weight_leaves(
    df: pd.DataFrame,
    syst_name: SystName,
    *,
    wgt_clip_hi: float = 20.0,
) -> None:
    """In-place: force ``ps*`` / ``ms*`` / ``morph`` / ``univ_*`` under a knob ≥ 0."""
    if df is None or len(df) == 0:
        return
    key = tuple(syst_name)
    try:
        block_cols = [c for c in df.columns if isinstance(c, tuple) and c[: len(key)] == key]
    except Exception:
        return
    for c in block_cols:
        leaf = next((str(p) for p in c[len(key) :] if p), "")
        if not (
            leaf.startswith("univ_")
            or leaf.startswith("ps")
            or leaf.startswith("ms")
            or leaf == "morph"
            or leaf == "cv"
        ):
            continue
        vals = np.asarray(df.loc[:, c], dtype=np.float64)
        df.loc[:, c] = np.clip(np.nan_to_num(vals, nan=1.0, posinf=1.0, neginf=0.0), 0.0, wgt_clip_hi)


def _infer_n_univ_block(block: pd.DataFrame) -> int:
    cols = block.columns
    leaves = []
    if isinstance(cols, pd.MultiIndex):
        for c in cols:
            if isinstance(c, tuple):
                for part in c:
                    if part:
                        leaves.append(str(part))
                        break
            else:
                leaves.append(str(c))
    else:
        leaves = [str(c) for c in cols]
    n = 0
    while f"univ_{n}" in leaves:
        n += 1
    if n > 0:
        return n
    leaf_set = set(leaves)
    if "ps1" in leaf_set and "ms1" in leaf_set:
        return 2
    if "ps1" in leaf_set or "morph" in leaf_set:
        return 1
    return 0


def _ensure_univ_aliases(evt_df: pd.DataFrame, mc_nu_df: Optional[pd.DataFrame], syst_name: SystName) -> int:
    """Map multisigma/morph leaves → ``univ_*`` (mutates frames). Returns n_univ.

    Non-multisim knobs (``ps1``/``ms1``/``morph``) are treated as discrete universes for
    covariance; their weights are clipped to ``≥ 0`` (never negative reweights).
    If a ``cv`` leaf exists, multisigma universes are ``ps|ms / cv``.
    """
    key = tuple(syst_name)
    try:
        block = evt_df.loc[:, key]
    except Exception:
        return 0
    n = _infer_n_univ_block(block)
    if n > 0 and any(
        (isinstance(c, tuple) and str(c[0]).startswith("univ_"))
        or (not isinstance(c, tuple) and str(c).startswith("univ_"))
        for c in block.columns
    ):
        clip_nonnegative_weight_leaves(evt_df, syst_name)
        if mc_nu_df is not None:
            clip_nonnegative_weight_leaves(mc_nu_df, syst_name)
        return n

    leaves = set()
    for c in block.columns:
        if isinstance(c, tuple):
            for part in c:
                if part:
                    leaves.add(str(part))
                    break
        else:
            leaves.add(str(c))

    def _copy_leaf(df: pd.DataFrame, src: str, dst: str, *, divide_by_cv: bool = False) -> None:
        if df is None or len(df) == 0:
            return
        # Resolve to a *full* MultiIndex column key (partial tuples can look "in"
        # columns but cannot be used to create new leaves).
        src_key = None
        cv_key = None
        for c in df.columns:
            if not isinstance(c, tuple) or c[: len(key)] != key:
                continue
            leaf = next((str(p) for p in c[len(key) :] if p), "")
            if leaf == src and src_key is None:
                src_key = c
            if leaf == "cv" and cv_key is None:
                cv_key = c
        if src_key is None:
            return
        dst_key = tuple(list(key) + [dst] + [""] * (len(src_key) - len(key) - 1))
        vals = np.asarray(df.loc[:, src_key], dtype=np.float64)
        # Physical: event weights cannot be negative (esp. ±σ / morph unisim).
        vals = np.nan_to_num(vals, nan=1.0, posinf=1.0, neginf=0.0)
        if divide_by_cv and cv_key is not None:
            cv = np.asarray(df.loc[:, cv_key], dtype=np.float64)
            cv = np.nan_to_num(cv, nan=1.0, posinf=1.0, neginf=1.0)
            cv = np.where(cv == 0.0, 1.0, cv)
            vals = vals / cv
        vals = np.clip(vals, 0.0, None)
        df.loc[:, dst_key] = vals

    # Match get_systematics_genie: GENIE_MULTISIGMA_DIVIDE_BY_CV (default on).
    _div_env = os.environ.get("GENIE_MULTISIGMA_DIVIDE_BY_CV", "1").strip().lower()
    _div_on = _div_env not in ("0", "false", "no", "off")
    div_cv = "cv" in leaves and _div_on
    if "ps1" in leaves and "ms1" in leaves:
        for src, dst in (("ps1", "univ_0"), ("ms1", "univ_1")):
            _copy_leaf(evt_df, src, dst, divide_by_cv=div_cv)
            if mc_nu_df is not None:
                _copy_leaf(mc_nu_df, src, dst, divide_by_cv=div_cv)
        clip_nonnegative_weight_leaves(evt_df, syst_name)
        if mc_nu_df is not None:
            clip_nonnegative_weight_leaves(mc_nu_df, syst_name)
        return 2
    if "ps1" in leaves:
        _copy_leaf(evt_df, "ps1", "univ_0", divide_by_cv=div_cv)
        if mc_nu_df is not None:
            _copy_leaf(mc_nu_df, "ps1", "univ_0", divide_by_cv=div_cv)
        clip_nonnegative_weight_leaves(evt_df, syst_name)
        if mc_nu_df is not None:
            clip_nonnegative_weight_leaves(mc_nu_df, syst_name)
        return 1
    if "morph" in leaves:
        _copy_leaf(evt_df, "morph", "univ_0")
        if mc_nu_df is not None:
            _copy_leaf(mc_nu_df, "morph", "univ_0")
        clip_nonnegative_weight_leaves(evt_df, syst_name)
        if mc_nu_df is not None:
            clip_nonnegative_weight_leaves(mc_nu_df, syst_name)
        return 1
    return n


def _empty_xsec_acc(n_univ: int, nb: int) -> Dict[str, np.ndarray]:
    return {
        "nevts_allmc": np.zeros(nb, dtype=np.float64),
        "cv_sel_reco": np.zeros(nb, dtype=np.float64),
        "cv_allsel_reco": np.zeros(nb, dtype=np.float64),
        "bg_cv": np.zeros(nb, dtype=np.float64),
        "reco_vs_true": np.zeros((n_univ, nb, nb), dtype=np.float64),
        "signal_allmc": np.zeros((n_univ, nb), dtype=np.float64),
        "signal_sel_truth": np.zeros((n_univ, nb), dtype=np.float64),
        "bg_univ": np.zeros((n_univ, nb), dtype=np.float64),
    }


def _accumulate_xsec_tensors(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    var_config: VariableConfig,
    syst_name: SystName,
    n_univ: int,
    acc: MutableMapping[str, np.ndarray],
) -> None:
    """Additive xsec tensors — identical math to ``get_systematics_genie.accumulate_xsec_path_chunk``."""
    from analysis_village.numucc_1p0pi.categories import topology_list

    bins = var_config.bins
    nb = len(bins) - 1
    evtdf_signal = mc_evt_df[mc_evt_df.topo_categ == 1]
    nudf_signal = mc_nu_df[mc_nu_df.topo_categ == 1]
    evtdf_div_topo = [mc_evt_df[mc_evt_df.topo_categ == mode] for mode in topology_list]

    ret = signal_hists(mc_evt_df, mc_nu_df, var_config, return_data=True, plot=False)
    nevts_allmc = ret["nevts_allmc"]
    if nevts_allmc is None:
        return
    acc["nevts_allmc"] += np.asarray(nevts_allmc, dtype=np.float64)
    acc["cv_sel_reco"] += np.asarray(ret["nevts_sel_reco"], dtype=np.float64)
    acc["cv_allsel_reco"] += np.asarray(ret["nevts_allsel_reco"], dtype=np.float64)

    wblock_evt = evtdf_signal[syst_name]
    wblock_nu = nudf_signal[syst_name]
    for uidx in range(n_univ):
        w_evt_univ = np.clip(
            np.asarray(genie_univ_weight_series(wblock_evt, uidx), dtype=np.float64), 0.0, 20.0
        )
        w_nu_univ = np.clip(
            np.asarray(genie_univ_weight_series(wblock_nu, uidx), dtype=np.float64), 0.0, 20.0
        )
        if nb == 1:
            reco_vs_true = np.array([[1.0]], dtype=np.float64)
        else:
            reco_vs_true, _, _ = np.histogram2d(
                ret["var_sel_truth"],
                ret["var_sel_reco"],
                weights=ret["wgt_sel_truth"] * w_evt_univ,
                bins=bins,
            )
        acc["reco_vs_true"][uidx] += reco_vs_true
        sam, _ = np.histogram(ret["var_allmc"], weights=ret["wgt_allmc"] * w_nu_univ, bins=bins)
        sst, _ = np.histogram(
            ret["var_sel_truth"], weights=ret["wgt_sel_truth"] * w_evt_univ, bins=bins
        )
        acc["signal_allmc"][uidx] += sam
        acc["signal_sel_truth"][uidx] += sst

    for this_evtdf in evtdf_div_topo[1:]:
        if this_evtdf is None or len(this_evtdf) == 0:
            continue
        var, wgt = get_clipped_evts(
            this_evtdf,
            var_config.var_evt_reco_col,
            bins,
            var_save_name=var_config.var_save_name,
        )
        acc["bg_cv"] += np.histogram(var, bins=bins, weights=wgt)[0].astype(np.float64)
        wblock_bg = this_evtdf[syst_name]
        for uidx in range(n_univ):
            uw = np.asarray(genie_univ_weight_series(wblock_bg, uidx), dtype=np.float64).copy()
            uw[np.isnan(uw)] = 1.0
            uw = np.clip(uw, 0.0, 20.0)
            acc["bg_univ"][uidx] += np.histogram(var, bins=bins, weights=wgt * uw)[0].astype(
                np.float64
            )


def finalize_genie_xsec_univ(
    acc: Mapping[str, np.ndarray],
    *,
    xsec_unit: float = 1.0,
) -> np.ndarray:
    """Universe xsec spectra from stored tensors (response × CV gen + Δbg)."""
    nevts_allmc = np.asarray(acc["nevts_allmc"], dtype=np.float64)
    nb = int(nevts_allmc.shape[0])
    n_univ = int(acc["reco_vs_true"].shape[0])
    scale = float(xsec_unit)
    rows: List[np.ndarray] = []
    for uidx in range(n_univ):
        if nb == 1:
            reco_vs_true = np.array([[1.0]], dtype=np.float64)
        else:
            reco_vs_true = acc["reco_vs_true"][uidx]
        denom = acc["signal_allmc"][uidx]
        eff = np.divide(
            acc["signal_sel_truth"][uidx],
            denom,
            out=np.zeros(nb, dtype=np.float64),
            where=denom != 0,
        )
        response = get_response_matrix(reco_vs_true, eff)
        signal_univ = response @ nevts_allmc
        signal_univ = signal_univ + (acc["bg_univ"][uidx] - acc["bg_cv"])
        signal_univ *= scale
        rows.append(signal_univ)
    return np.asarray(rows, dtype=np.float64)


def finalize_genie_xsec_cv(
    acc: Mapping[str, np.ndarray],
    *,
    xsec_unit: float = 1.0,
    bkgd_subtract: bool = True,
) -> np.ndarray:
    scale = float(xsec_unit)
    base = acc["cv_sel_reco"] if bkgd_subtract else acc["cv_allsel_reco"]
    return np.asarray(base, dtype=np.float64) * scale


# ---------------------------------------------------------------------------
# Pack / unpack long DataFrame
# ---------------------------------------------------------------------------
def _rows_1d(
    kind: str,
    family: str,
    knob: str,
    var: str,
    univ: int,
    arr: np.ndarray,
    *,
    sparse: bool = True,
) -> List[dict]:
    """Emit long-format rows. With ``sparse=True`` (default), skip exact zeros.

    Unpack reconstructs dense arrays from ``bin0`` / ``univ`` indices, so zeros are
    optional. Sparse packing is essential for GENIE xsec migration matrices.

    Always emit a shape sentinel at ``bin0 = nbins - 1`` (value 0 if that bin is
    empty) so trailing empty bins survive ``max(bin0)+1`` unpack.
    """
    a = np.asarray(arr, dtype=np.float64).reshape(-1)
    rows = []
    for i in range(len(a)):
        v = float(a[i])
        if sparse and v == 0.0:
            continue
        rows.append(
            {
                "kind": kind,
                "family": family,
                "knob": knob,
                "var": var,
                "univ": int(univ),
                "bin0": int(i),
                "bin1": -1,
                "value": v,
            }
        )
    if len(a) > 0:
        last = int(len(a) - 1)
        if not rows or max(r["bin0"] for r in rows) < last:
            rows.append(
                {
                    "kind": kind,
                    "family": family,
                    "knob": knob,
                    "var": var,
                    "univ": int(univ),
                    "bin0": last,
                    "bin1": -1,
                    "value": 0.0,
                }
            )
    return rows


def _rows_2d(
    kind: str,
    family: str,
    knob: str,
    var: str,
    univ: int,
    mat: np.ndarray,
    *,
    sparse: bool = True,
) -> List[dict]:
    m = np.asarray(mat, dtype=np.float64)
    rows = []
    nz = np.argwhere(m != 0.0) if sparse else np.indices(m.shape).reshape(2, -1).T
    for i, j in nz:
        rows.append(
            {
                "kind": kind,
                "family": family,
                "knob": knob,
                "var": var,
                "univ": int(univ),
                "bin0": int(i),
                "bin1": int(j),
                "value": float(m[i, j]),
            }
        )
    # Always keep true (n_reco, n_true) via a corner sentinel when needed.
    if m.size > 0:
        last0 = int(m.shape[0] - 1)
        last1 = int(m.shape[1] - 1)
        max0 = max((r["bin0"] for r in rows), default=-1)
        max1 = max((r["bin1"] for r in rows), default=-1)
        if max0 < last0 or max1 < last1:
            rows.append(
                {
                    "kind": kind,
                    "family": family,
                    "knob": knob,
                    "var": var,
                    "univ": int(univ),
                    "bin0": last0,
                    "bin1": last1,
                    "value": 0.0,
                }
            )
    return rows


def _resolve_nbins(
    slug: str,
    observed: int,
    nbins_by_var: Optional[Mapping[str, int]],
) -> int:
    nb = max(int(observed), 0)
    if nbins_by_var is not None and slug in nbins_by_var:
        nb = max(nb, int(nbins_by_var[slug]))
    return nb


def _pad_1d(arr: np.ndarray, nb: int) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64).reshape(-1)
    if len(a) == nb:
        return a
    if len(a) > nb:
        return a[:nb].copy()
    out = np.zeros(nb, dtype=np.float64)
    out[: len(a)] = a
    return out


def _pad_2d_univ(arr: np.ndarray, nb: int) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    if a.ndim != 2:
        a = a.reshape(1, -1)
    n_univ, n_old = a.shape
    if n_old == nb:
        return a
    out = np.zeros((n_univ, nb), dtype=np.float64)
    n_copy = min(n_old, nb)
    out[:, :n_copy] = a[:, :n_copy]
    return out


def pack_blob_to_df(blob: Mapping[str, Any]) -> pd.DataFrame:
    """Flatten rate + optional xsec accumulators into the long ``syst_hists`` schema."""
    rows: List[dict] = []
    family = str(blob.get("family", "unknown"))
    rate = blob.get("rate", {})
    for knob, vars_d in rate.items():
        for slug, pack in vars_d.items():
            rows.extend(_rows_1d(KIND_RATE_CV, family, knob, slug, -1, pack["cv"]))
            univ = np.asarray(pack["univ"], dtype=np.float64)
            for u in range(univ.shape[0]):
                rows.extend(_rows_1d(KIND_RATE_UNIV, family, knob, slug, u, univ[u]))

    xsec = blob.get("xsec", {})
    for knob, vars_d in xsec.items():
        for slug, acc in vars_d.items():
            rows.extend(
                _rows_1d(KIND_XSEC_NEVTS_ALLMC, family, knob, slug, -1, acc["nevts_allmc"])
            )
            rows.extend(
                _rows_1d(KIND_XSEC_CV_SEL_RECO, family, knob, slug, -1, acc["cv_sel_reco"])
            )
            rows.extend(
                _rows_1d(KIND_XSEC_CV_ALLSEL_RECO, family, knob, slug, -1, acc["cv_allsel_reco"])
            )
            rows.extend(_rows_1d(KIND_XSEC_BG_CV, family, knob, slug, -1, acc["bg_cv"]))
            n_univ = int(acc["reco_vs_true"].shape[0])
            for u in range(n_univ):
                rows.extend(
                    _rows_1d(KIND_XSEC_SIGNAL_ALLMC, family, knob, slug, u, acc["signal_allmc"][u])
                )
                rows.extend(
                    _rows_1d(
                        KIND_XSEC_SIGNAL_SEL_TRUTH, family, knob, slug, u, acc["signal_sel_truth"][u]
                    )
                )
                rows.extend(_rows_1d(KIND_XSEC_BG_UNIV, family, knob, slug, u, acc["bg_univ"][u]))
                rows.extend(
                    _rows_2d(
                        KIND_XSEC_RECO_VS_TRUE, family, knob, slug, u, acc["reco_vs_true"][u]
                    )
                )

    if not rows:
        return empty_histcounts_df()
    return pd.DataFrame(rows, columns=list(HIST_DF_COLUMNS))


def sum_histcounts_dfs(dfs: Iterable[pd.DataFrame]) -> pd.DataFrame:
    """Sum ``value`` over identical keys (disjoint CAF / grid chunks)."""
    parts = [d for d in dfs if d is not None and len(d) > 0]
    if not parts:
        return empty_histcounts_df()
    cat = pd.concat(parts, ignore_index=True)
    keys = ["kind", "family", "knob", "var", "univ", "bin0", "bin1"]
    return cat.groupby(keys, as_index=False, sort=False)["value"].sum()


def _apply_hist_filters(
    df: pd.DataFrame,
    *,
    family: Optional[str] = None,
    knob: Optional[str] = None,
    var: Optional[str] = None,
    vars: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Subset long histcounts by family / knob / variable slug(s)."""
    sub = df
    if family is not None:
        sub = sub[sub["family"] == family]
    if knob is not None:
        sub = sub[sub["knob"] == knob]
    if vars is not None:
        keep = {str(v) for v in vars}
        sub = sub[sub["var"].astype(str).isin(keep)]
    elif var is not None:
        sub = sub[sub["var"] == var]
    return sub


def unpack_rate_from_df(
    df: pd.DataFrame,
    *,
    family: Optional[str] = None,
    knob: Optional[str] = None,
    var: Optional[str] = None,
    vars: Optional[Sequence[str]] = None,
    nbins_by_var: Optional[Mapping[str, int]] = None,
) -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
    """``out[knob][var] = {cv, univ}`` from summed histcounts.

    ``nbins_by_var`` pads to the frozen VariableConfig length (safety net for
    older files that lacked a trailing-bin shape sentinel).
    Pass ``vars`` to keep only final observables (skips cut-stage histograms).
    """
    if df is None or len(df) == 0:
        return {}
    sub = _apply_hist_filters(df, family=family, knob=knob, var=var, vars=vars)
    out: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    cv_rows = sub[sub["kind"] == KIND_RATE_CV]
    univ_rows = sub[sub["kind"] == KIND_RATE_UNIV]
    for (kn, slug), g in cv_rows.groupby(["knob", "var"]):
        g2 = g.sort_values("bin0")
        slug_s = str(slug)
        nb = _resolve_nbins(slug_s, int(g2["bin0"].max()) + 1, nbins_by_var)
        cv = np.zeros(nb, dtype=np.float64)
        idx = g2["bin0"].astype(int).values
        ok = idx < nb
        cv[idx[ok]] = g2["value"].astype(float).values[ok]
        out.setdefault(str(kn), {})[slug_s] = {
            "cv": cv,
            "univ": np.zeros((0, nb), dtype=np.float64),
        }
    for (kn, slug), g in univ_rows.groupby(["knob", "var"]):
        g2 = g.sort_values(["univ", "bin0"])
        slug_s = str(slug)
        n_univ = int(g2["univ"].max()) + 1
        nb = _resolve_nbins(slug_s, int(g2["bin0"].max()) + 1, nbins_by_var)
        univ = np.zeros((n_univ, nb), dtype=np.float64)
        u_idx = g2["univ"].astype(int).values
        b_idx = g2["bin0"].astype(int).values
        ok = b_idx < nb
        univ[u_idx[ok], b_idx[ok]] = g2["value"].astype(float).values[ok]
        slot = out.setdefault(str(kn), {}).setdefault(
            slug_s, {"cv": np.zeros(nb, dtype=np.float64), "univ": univ}
        )
        slot["univ"] = univ
        if slot["cv"].shape[0] != nb:
            slot["cv"] = _pad_1d(slot["cv"], nb)
    return out


def unpack_xsec_from_df(
    df: pd.DataFrame,
    *,
    family: Optional[str] = None,
    knob: Optional[str] = None,
    var: Optional[str] = None,
    vars: Optional[Sequence[str]] = None,
    nbins_by_var: Optional[Mapping[str, int]] = None,
) -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
    """``out[knob][var] = xsec accumulator dict`` ready for :func:`finalize_genie_xsec_univ`.

    Pass ``vars`` to keep only selected observables (skips cut-stage histograms).
    """
    if df is None or len(df) == 0:
        return {}
    sub = _apply_hist_filters(df, family=family, knob=knob, var=var, vars=vars)
    xkinds = {
        KIND_XSEC_NEVTS_ALLMC,
        KIND_XSEC_CV_SEL_RECO,
        KIND_XSEC_CV_ALLSEL_RECO,
        KIND_XSEC_BG_CV,
        KIND_XSEC_SIGNAL_ALLMC,
        KIND_XSEC_SIGNAL_SEL_TRUTH,
        KIND_XSEC_BG_UNIV,
        KIND_XSEC_RECO_VS_TRUE,
    }
    sub = sub[sub["kind"].isin(xkinds)]
    if len(sub) == 0:
        return {}

    out: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    for (kn, slug), g in sub.groupby(["knob", "var"]):
        slug_s = str(slug)
        obs = int(
            max(
                g["bin0"].max(),
                g.loc[g["bin1"] >= 0, "bin1"].max()
                if (g["bin1"] >= 0).any()
                else g["bin0"].max(),
            )
        ) + 1
        nb = _resolve_nbins(slug_s, obs, nbins_by_var)
        univ_max = g.loc[g["univ"] >= 0, "univ"].max() if (g["univ"] >= 0).any() else -1
        n_univ = int(univ_max) + 1 if univ_max >= 0 else 0
        acc = _empty_xsec_acc(max(n_univ, 1), nb)
        if n_univ == 0:
            acc["reco_vs_true"] = np.zeros((0, nb, nb), dtype=np.float64)
            acc["signal_allmc"] = np.zeros((0, nb), dtype=np.float64)
            acc["signal_sel_truth"] = np.zeros((0, nb), dtype=np.float64)
            acc["bg_univ"] = np.zeros((0, nb), dtype=np.float64)

        def _fill_1d(kind: str, target: np.ndarray, univ_axis: bool) -> None:
            gg = g[g["kind"] == kind]
            if len(gg) == 0:
                return
            if univ_axis:
                target[gg["univ"].astype(int).values, gg["bin0"].astype(int).values] = (
                    gg["value"].astype(float).values
                )
            else:
                target[gg["bin0"].astype(int).values] = gg["value"].astype(float).values

        _fill_1d(KIND_XSEC_NEVTS_ALLMC, acc["nevts_allmc"], False)
        _fill_1d(KIND_XSEC_CV_SEL_RECO, acc["cv_sel_reco"], False)
        _fill_1d(KIND_XSEC_CV_ALLSEL_RECO, acc["cv_allsel_reco"], False)
        _fill_1d(KIND_XSEC_BG_CV, acc["bg_cv"], False)
        if n_univ > 0:
            _fill_1d(KIND_XSEC_SIGNAL_ALLMC, acc["signal_allmc"], True)
            _fill_1d(KIND_XSEC_SIGNAL_SEL_TRUTH, acc["signal_sel_truth"], True)
            _fill_1d(KIND_XSEC_BG_UNIV, acc["bg_univ"], True)
            gg = g[g["kind"] == KIND_XSEC_RECO_VS_TRUE]
            if len(gg) > 0:
                acc["reco_vs_true"][
                    gg["univ"].astype(int).values,
                    gg["bin0"].astype(int).values,
                    gg["bin1"].astype(int).values,
                ] = gg["value"].astype(float).values
        out.setdefault(str(kn), {})[str(slug)] = acc
    return out


def load_syst_hists_from_df_file(
    path: str,
    *,
    vars: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Read and sum all ``syst_hists_*`` splits from one ``run_df_maker`` output.

    Optional ``vars`` drops other slugs immediately after read (HDF Fixed format
    still loads the full table — this only shrinks the in-memory working set).
    """
    parts: List[pd.DataFrame] = []
    with pd.HDFStore(path, mode="r") as store:
        keys = [k.lstrip("/") for k in store.keys()]
        hist_keys = sorted(k for k in keys if k.startswith("syst_hists"))
        for k in hist_keys:
            parts.append(store[k])
    out = sum_histcounts_dfs(parts)
    if vars is not None and len(out) > 0:
        keep = {str(v) for v in vars}
        out = out[out["var"].astype(str).isin(keep)].reset_index(drop=True)
    return out


def sum_histcounts_paths_streaming(
    paths: Sequence[str],
    *,
    progress_every: int = 10,
    progress_cb: Optional[Callable[[int, int, str], None]] = None,
) -> pd.DataFrame:
    """Sum histcounts across files **one file at a time** (never materialize all).

    Peak RAM is ~2× one file (current accumulator + next file during concat),
    not ``N_files ×`` one file. Prefer :func:`stream_sum_rate_dense` for Flux/G4
    rate-only campaigns (smaller resident set after unpack).
    """
    acc: Optional[pd.DataFrame] = None
    n_paths = len(paths)
    for i, p in enumerate(paths, start=1):
        df = load_syst_hists_from_df_file(p)
        if acc is None:
            acc = df
        else:
            acc = sum_histcounts_dfs([acc, df])
            del df
            gc.collect()
        if progress_cb is not None:
            progress_cb(i, n_paths, p)
        elif progress_every > 0 and (i % progress_every == 0 or i == n_paths):
            n_rows = 0 if acc is None else len(acc)
            print(
                "[histcounts-stream] %d/%d  rows=%d  %s"
                % (i, n_paths, n_rows, os.path.basename(p)),
                flush=True,
            )
    return empty_histcounts_df() if acc is None else acc


def add_rate_dense_inplace(
    acc: MutableMapping[str, Dict[str, Dict[str, np.ndarray]]],
    rate: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
) -> None:
    """``acc[knob][var][{cv,univ}] += rate[...]`` with shape padding."""

    def _fit_univ(a: np.ndarray, n_univ: int, nb: int) -> np.ndarray:
        a = np.asarray(a, dtype=np.float64)
        if a.size == 0:
            return np.zeros((n_univ, nb), dtype=np.float64)
        if a.ndim != 2:
            a = a.reshape(1, -1)
        a = _pad_2d_univ(a, nb)
        if a.shape[0] == n_univ:
            return a
        out = np.zeros((n_univ, nb), dtype=np.float64)
        n_copy = min(int(a.shape[0]), n_univ)
        out[:n_copy, :] = a[:n_copy, :]
        return out

    for knob, vars_d in rate.items():
        slot_k = acc.setdefault(str(knob), {})
        for slug, pack in vars_d.items():
            cv = np.asarray(pack["cv"], dtype=np.float64).reshape(-1)
            univ = np.asarray(pack["univ"], dtype=np.float64)
            if slug not in slot_k:
                if univ.size:
                    univ = univ.reshape(-1, univ.shape[-1]) if univ.ndim == 1 else univ
                    slot_k[str(slug)] = {"cv": cv.copy(), "univ": univ.copy()}
                else:
                    slot_k[str(slug)] = {
                        "cv": cv.copy(),
                        "univ": np.zeros((0, cv.shape[0]), dtype=np.float64),
                    }
                continue
            cur = slot_k[str(slug)]
            nb = max(int(cur["cv"].shape[0]), int(cv.shape[0]))
            cur["cv"] = _pad_1d(cur["cv"], nb) + _pad_1d(cv, nb)
            if univ.size == 0:
                if cur["univ"].size:
                    cur["univ"] = _pad_2d_univ(cur["univ"], nb)
                continue
            if univ.ndim != 2:
                univ = univ.reshape(1, -1)
            n_univ = max(
                int(cur["univ"].shape[0]) if cur["univ"].size else 0,
                int(univ.shape[0]),
            )
            cur["univ"] = _fit_univ(cur["univ"], n_univ, nb) + _fit_univ(univ, n_univ, nb)


def stream_sum_rate_dense(
    paths: Sequence[str],
    *,
    family: Optional[str] = None,
    nbins_by_var: Optional[Mapping[str, int]] = None,
    vars: Optional[Sequence[str]] = None,
    progress_every: int = 10,
) -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
    """Stream files → dense rate accumulators ``out[knob][var] = {cv, univ}``.

    Intended for Flux / G4 (rate-only). Loads one long hist DF at a time, unpacks,
    adds into numpy arrays, then drops the DF. Pass ``vars`` to skip cut-stage
    histograms (much faster cov build).
    """
    acc: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    n_paths = len(paths)
    for i, p in enumerate(paths, start=1):
        df = load_syst_hists_from_df_file(p, vars=vars)
        rate = unpack_rate_from_df(
            df, family=family, vars=vars, nbins_by_var=nbins_by_var
        )
        del df
        add_rate_dense_inplace(acc, rate)
        del rate
        gc.collect()
        if progress_every > 0 and (i % progress_every == 0 or i == n_paths):
            n_knob = len(acc)
            n_var = len(next(iter(acc.values()), {}))
            print(
                "[histcounts-stream-dense] %d/%d  knobs=%d  vars~=%d  %s"
                % (i, n_paths, n_knob, n_var, os.path.basename(p)),
                flush=True,
            )
    return acc


def rate_dense_to_hist_df(
    rate: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    *,
    family: str,
) -> pd.DataFrame:
    """Pack dense rate accumulators back into the long ``syst_hists`` schema."""
    return pack_blob_to_df({"family": family, "rate": rate, "xsec": {}})


def add_xsec_dense_inplace(
    acc: MutableMapping[str, Dict[str, Dict[str, np.ndarray]]],
    xsec: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
) -> None:
    """Add unpacked xsec accumulators into ``acc[knob][var]`` (GENIE)."""

    def _add_arr(dst: np.ndarray, src: np.ndarray) -> np.ndarray:
        a = np.asarray(dst, dtype=np.float64)
        b = np.asarray(src, dtype=np.float64)
        if a.shape == b.shape:
            return a + b
        # Pad trailing dims to the max shape (bin / univ growth across files).
        out_shape = tuple(max(x, y) for x, y in zip(a.shape, b.shape))
        if len(a.shape) != len(b.shape):
            raise ValueError("xsec array rank mismatch: %s vs %s" % (a.shape, b.shape))
        aa = np.zeros(out_shape, dtype=np.float64)
        bb = np.zeros(out_shape, dtype=np.float64)
        aa[tuple(slice(0, s) for s in a.shape)] = a
        bb[tuple(slice(0, s) for s in b.shape)] = b
        return aa + bb

    for knob, vars_d in xsec.items():
        slot_k = acc.setdefault(str(knob), {})
        for slug, pack in vars_d.items():
            if slug not in slot_k:
                slot_k[str(slug)] = {
                    k: np.array(v, dtype=np.float64, copy=True) for k, v in pack.items()
                }
                continue
            cur = slot_k[str(slug)]
            for key, arr in pack.items():
                if key not in cur:
                    cur[key] = np.array(arr, dtype=np.float64, copy=True)
                else:
                    cur[key] = _add_arr(cur[key], arr)


def stream_sum_genie_dense(
    paths: Sequence[str],
    *,
    nbins_by_var: Optional[Mapping[str, int]] = None,
    vars: Optional[Sequence[str]] = None,
    progress_every: int = 10,
    max_files: int = 0,
) -> Tuple[Dict[str, Dict[str, Dict[str, np.ndarray]]], Dict[str, Dict[str, Dict[str, np.ndarray]]]]:
    """Stream GENIE histcount files → dense ``(rate, xsec)`` (one file at a time).

    Peak RAM ≈ one file (~10 GiB for Spring MC GENIE) + dense accumulators.
    Never materializes all long DataFrames together.

    Pass ``vars`` (e.g. final selected observables only) to skip cut-stage
    histograms — unpack/cov become much cheaper; HDF Fixed still reads full file.
    """
    files = list(paths)
    if max_files > 0:
        files = files[: int(max_files)]
    rate_acc: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    xsec_acc: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    n_paths = len(files)
    for i, p in enumerate(files, start=1):
        df = load_syst_hists_from_df_file(p, vars=vars)
        rate = unpack_rate_from_df(
            df, family="GENIE", vars=vars, nbins_by_var=nbins_by_var
        )
        xsec = unpack_xsec_from_df(
            df, family="GENIE", vars=vars, nbins_by_var=nbins_by_var
        )
        del df
        add_rate_dense_inplace(rate_acc, rate)
        add_xsec_dense_inplace(xsec_acc, xsec)
        del rate, xsec
        gc.collect()
        if progress_every > 0 and (i % progress_every == 0 or i == n_paths):
            print(
                "[genie-stream-dense] %d/%d  rate_knobs=%d  xsec_knobs=%d  %s"
                % (
                    i,
                    n_paths,
                    len(rate_acc),
                    len(xsec_acc),
                    os.path.basename(p),
                ),
                flush=True,
            )
    return rate_acc, xsec_acc


def build_genie_covs_from_dense(
    rate: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    xsec: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    *,
    xsec_unit: float = 1.0,
    slim_skip: Optional[Sequence[str]] = None,
) -> Dict[str, Dict[str, Any]]:
    """Build per-var GENIE cov packs from dense rate/xsec (same keys as notebook helper)."""
    skip = set(slim_skip or ())
    by_var: Dict[str, Dict[str, Any]] = {}
    rate_packs_by_var: Dict[str, List[Any]] = {}
    xsec_packs_by_var: Dict[str, List[Any]] = {}
    cv_rate_by_var: Dict[str, np.ndarray] = {}
    cv_xsec_by_var: Dict[str, np.ndarray] = {}

    for knob, vars_d in rate.items():
        for slug, pack in vars_d.items():
            univ = np.asarray(pack["univ"], dtype=float)
            if univ.size == 0:
                continue
            rp = rate_cov_from_univ_cv(univ, pack["cv"])
            by_var.setdefault(str(slug), {})[f"{knob}_rate"] = rp
            if knob not in skip:
                rate_packs_by_var.setdefault(str(slug), []).append(rp)
            cv_rate_by_var[str(slug)] = np.asarray(pack["cv"], dtype=float)

    for knob, vars_d in xsec.items():
        for slug, acc in vars_d.items():
            univ = finalize_genie_xsec_univ(acc, xsec_unit=xsec_unit)
            cv = finalize_genie_xsec_cv(acc, xsec_unit=xsec_unit, bkgd_subtract=True)
            if univ.size == 0:
                continue
            xp = rate_cov_from_univ_cv(univ, cv)
            by_var.setdefault(str(slug), {})[str(knob)] = xp
            if knob not in skip:
                xsec_packs_by_var.setdefault(str(slug), []).append(xp)
            cv_xsec_by_var[str(slug)] = np.asarray(cv, dtype=float)

    for slug in set(rate_packs_by_var) | set(xsec_packs_by_var):
        if slug in rate_packs_by_var:
            by_var.setdefault(slug, {})["genie_rate"] = combine_indep_knob_frac_covs(
                rate_packs_by_var[slug], cv_rate_by_var[slug]
            )
        if slug in xsec_packs_by_var:
            by_var.setdefault(slug, {})["genie"] = combine_indep_knob_frac_covs(
                xsec_packs_by_var[slug], cv_xsec_by_var[slug]
            )
    return by_var


def load_syst_hists_from_glob(paths: Sequence[str]) -> pd.DataFrame:
    """Sum histcounts across files via streaming (never loads all files at once)."""
    return sum_histcounts_paths_streaming(paths)


def rate_cov_from_univ_cv(
    univ: np.ndarray,
    cv: np.ndarray,
) -> Dict[str, np.ndarray]:
    return get_covariance_matrix(np.asarray(univ, dtype=float), np.asarray(cv, dtype=float))


def unisim_cov_from_cv_and_var(
    n_cv: np.ndarray,
    n_var: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Single-universe unisim (WireMod / DENT / cosmics template)."""
    univ = np.asarray(n_var, dtype=float).reshape(1, -1)
    cv = np.asarray(n_cv, dtype=float).reshape(-1)
    if univ.shape[1] != cv.shape[0]:
        raise ValueError(
            "unisim CV/var bin mismatch: cv has %d bins, var has %d"
            % (cv.shape[0], univ.shape[1])
        )
    return get_covariance_matrix(univ, cv)


def scale_cov_by_contamination(
    cov_pack: Mapping[str, np.ndarray],
    frac: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Cosmic template → selected-rate: ``cov_sel = cov_tmpl * outer(f, f)``."""
    f = np.asarray(frac, dtype=float).reshape(-1)
    outer = np.outer(f, f)
    cov = np.asarray(cov_pack["cov"], dtype=float) * outer
    cov_frac = np.asarray(cov_pack["cov_frac"], dtype=float) * outer
    corr = np.zeros_like(cov)
    eps = 1e-12
    for i in range(cov.shape[0]):
        for j in range(cov.shape[1]):
            di = max(float(cov[i, i]), 0.0)
            dj = max(float(cov[j, j]), 0.0)
            denom = np.sqrt(di * dj)
            corr[i, j] = (cov[i, j] / denom) if denom > eps else 0.0
    return {"cov": cov, "cov_frac": cov_frac, "corr": corr}


def combine_indep_knob_frac_covs(
    packs: Sequence[Mapping[str, np.ndarray]],
    cv: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Sum fractional covariances (independent knobs), rebuild absolute cov from ``cv``."""
    if not packs:
        nb = len(cv)
        z = np.zeros((nb, nb), dtype=float)
        return {"cov": z, "cov_frac": z, "corr": z}
    cov_frac = np.sum([np.asarray(p["cov_frac"], dtype=float) for p in packs], axis=0)
    v = np.asarray(cv, dtype=float).reshape(-1)
    cov = cov_frac * np.outer(v, v)
    corr = np.zeros_like(cov_frac)
    eps = 1e-12
    for i in range(cov_frac.shape[0]):
        for j in range(cov_frac.shape[1]):
            di = max(float(cov_frac[i, i]), 0.0)
            dj = max(float(cov_frac[j, j]), 0.0)
            denom = np.sqrt(di * dj)
            corr[i, j] = (cov_frac[i, j] / denom) if denom > eps else 0.0
    return {"cov": cov, "cov_frac": cov_frac, "corr": corr}


# ---------------------------------------------------------------------------
# Fill from in-memory evt/trk/mcnu (called by makedf wrappers)
# ---------------------------------------------------------------------------
def fill_syst_histcounts(
    evt: pd.DataFrame,
    trk: pd.DataFrame,
    mcnu: Optional[pd.DataFrame],
    *,
    family: str,
    syst_names: Sequence[SystName],
    sample: str = "mc",
    do_xsec: bool = False,
    bkgd_subtract: bool = True,
) -> pd.DataFrame:
    """Walk selection; fill rate (± xsec tensors for GENIE); return long DataFrame."""
    if evt is None or len(evt) == 0:
        return empty_histcounts_df()

    evt = evt.copy()
    trk = trk.copy() if trk is not None else None
    mcnu = mcnu.copy() if mcnu is not None else None
    if mcnu is not None:
        mcnu = _prefix_mcnu_columns(mcnu)
        mcnu = ensure_mc_level_phi_mcnu(mcnu)
    _annotate_topo(evt, mcnu)

    blob: Dict[str, Any] = {"family": family, "rate": {}, "xsec": {}}
    cut_by_stage: Dict[str, List[Any]] = {}
    for spec in CUT_STAGE_VAR_SPECS:
        cut_by_stage.setdefault(spec.stage_key, []).append(spec)
    final_vcs = final_var_configs()

    # Resolve n_univ per knob (with aliases)
    knob_n_univ: Dict[str, int] = {}
    live_systs: List[SystName] = []
    for sn in syst_names:
        n_u = _ensure_univ_aliases(evt, mcnu, sn)
        if n_u <= 0 and family.upper() != "UNISIM":
            continue
        knob_n_univ[sn[1]] = max(n_u, 0)
        live_systs.append(sn)

    if family.upper() == "UNISIM":
        # No weights: still walk and fill CV-only rate histograms under knob="cv".
        live_systs = [("mc", "cv")]
        knob_n_univ["cv"] = 0

    state0: Dict[str, Any] = {"evt": evt, "trk": trk, "hdr": None, "mcnu": None}
    for stage_key, post_state in walk_pipeline(state0, sample=sample):
        post_evt = post_state.get("evt")
        if post_evt is None or len(post_evt) == 0:
            continue
        pe = post_evt
        pn = None
        if mcnu is not None:
            pe, pn = _align_evt_mcnu(pe, mcnu)
            if len(pe) == 0:
                continue
            post_state = dict(post_state)
            post_state["evt"] = pe

        # Cut-stage rate (simple reweight; no xsec)
        for spec in cut_by_stage.get(stage_key, ()):
            extracted = get_var_series(post_state, spec.var_config, spec.target)
            if extracted is None:
                continue
            values, evt_idx = extracted
            bins = np.asarray(spec.var_config.bins, dtype=float)
            nb = len(bins) - 1
            slug = spec.var_config.var_save_name
            cv_h = histogram_var(values, bins)
            for sn in live_systs:
                knob = sn[1]
                n_univ = knob_n_univ.get(knob, 0)
                rate_slot = blob["rate"].setdefault(knob, {}).setdefault(
                    slug,
                    {
                        "univ": np.zeros((max(n_univ, 0), nb), dtype=np.float64),
                        "cv": np.zeros(nb, dtype=np.float64),
                    },
                )
                rate_slot["cv"] += cv_h
                if n_univ <= 0 or family.upper() == "UNISIM":
                    continue
                try:
                    wblock = pe[sn]
                except Exception:
                    continue
                for u in range(n_univ):
                    w = np.asarray(genie_univ_weight_series(wblock, u), dtype=np.float64)
                    w = np.clip(np.nan_to_num(w, nan=1.0), 0.0, 20.0)
                    rate_slot["univ"][u] += histogram_var(values, bins, weights=w[evt_idx])

        if stage_key != FINAL_STAGE_KEY:
            continue

        # Final stage: derived kinematics + full rate (bkg-sub) + optional xsec tensors
        pe = ensure_derived_trk_kinematics_cols(pe)
        pe = add_reco_cc1p0pi_tki_evtdf(pe)
        pe = add_truth_cc1p0pi_tki_evtdf(pe)
        if pn is not None and len(pn) > 0:
            pn = add_mc_cc1p0pi_tki_mcnu(pn)
        _annotate_topo(pe, pn)

        for vc in final_vcs:
            slug = vc.var_save_name
            bins = np.asarray(vc.bins, dtype=float)
            nb = len(bins) - 1

            if family.upper() == "UNISIM" or not live_systs or knob_n_univ.get(live_systs[0][1], 0) == 0:
                # CV-only selected reco counts (no bkg subtraction weights)
                try:
                    var, wgt = get_clipped_evts(pe, vc.var_evt_reco_col, bins, var_save_name=slug)
                    cv_h, _ = np.histogram(var, bins=bins, weights=wgt)
                except Exception:
                    cv_h = np.zeros(nb, dtype=np.float64)
                blob["rate"].setdefault("cv", {}).setdefault(
                    slug, {"univ": np.zeros((0, nb), dtype=np.float64), "cv": np.zeros(nb)}
                )
                blob["rate"]["cv"][slug]["cv"] += np.asarray(cv_h, dtype=np.float64)
                continue

            for sn in live_systs:
                knob = sn[1]
                n_univ = knob_n_univ[knob]
                if n_univ <= 0:
                    continue
                try:
                    univ_r, cv_r = get_univ_rates(
                        cov_type="rate",
                        syst_type="GENIE" if family.upper() == "GENIE" else family.upper(),
                        evtdf=pe,
                        nudf=pn,
                        var_config=vc,
                        syst_name=sn,
                        n_univ=n_univ,
                        bkgd_subtract=bkgd_subtract,
                        plot=False,
                    )
                except Exception:
                    continue
                slot = blob["rate"].setdefault(knob, {}).setdefault(
                    slug,
                    {
                        "univ": np.zeros_like(univ_r, dtype=np.float64),
                        "cv": np.zeros_like(cv_r, dtype=np.float64),
                    },
                )
                slot["univ"] += np.asarray(univ_r, dtype=np.float64)
                slot["cv"] += np.asarray(cv_r, dtype=np.float64)

                if do_xsec and pn is not None and slug not in CUT_STAGE_RATE_ONLY_SLUGS:
                    xslot = blob["xsec"].setdefault(knob, {})
                    if slug not in xslot:
                        xslot[slug] = _empty_xsec_acc(n_univ, nb)
                    try:
                        _accumulate_xsec_tensors(pe, pn, vc, sn, n_univ, xslot[slug])
                    except Exception:
                        pass

    return pack_blob_to_df(blob)
