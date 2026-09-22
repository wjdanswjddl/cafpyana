"""Helpers for GENIE systematics inspection / plotting notebooks.

Expects ``cov_mat_dict.pkl`` layout written by ``syst_genie_aggregate``::

    {var: {"genie": cov_frac_xsec, "genie_rate": cov_frac_rate,
           "<knob>": cov_frac_xsec, "<knob>_rate": cov_frac_rate, ...}}
"""
from __future__ import annotations

import csv
import pickle
import re
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def _text_color_for(value: float, cmap_name: str, vmin: float, vmax: float) -> str:
    """Black/white label color from colormap luminance at ``value``."""
    if vmax <= vmin:
        return "black"
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    rgba = plt.get_cmap(cmap_name)(norm(value))
    lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return "black" if lum > 0.5 else "white"

# ---------------------------------------------------------------------------
# Paths / constants
# ---------------------------------------------------------------------------

MODE_ORDER: Tuple[str, ...] = ("CCQE", "MEC", "RES", "nonRES", "DIS", "Other", "Ar23p")
PHYSICS_MODES: Tuple[str, ...] = ("CCQE", "MEC", "RES", "nonRES", "DIS", "Other")  # no Ar23p
# Distributed Ar23p: Other split into Other (COH/NCEL) + FSI (rest of Other + leftover Ar23p)
DISTRIBUTED_MODE_ORDER: Tuple[str, ...] = (
    "CCQE", "MEC", "RES", "nonRES", "DIS", "Other", "FSI",
)

MODE_COLORS: Dict[str, str] = {
    "CCQE": "#1f77b4",
    "MEC": "#ff7f0e",
    "RES": "#2ca02c",
    "nonRES": "#17becf",
    "DIS": "#d62728",
    "Other": "#9467bd",
    "FSI": "#e377c2",
    "Ar23p": "#8c564b",
    "Total": "k",
}

_STRIP_PREFIXES: Tuple[str, ...] = (
    "GENIEReWeight_SBN_v1_multisim_",
    "GENIEReWeight_SBN_v1_multisigma_",
    "GENIEReWeight_SBN_v3_",
    "GENIEReWeight_",
    "CCQETemplateReweight_SBN_v3_",
    "CCQETemplateReweight_SBNNuSyst_",
    "CCQETemplateReweight_",
    "QEInterference_SBN_v3_",
    "QEInterference_",
    "ZExpPCAWeighter_SBN_v3_",
    "ZExpPCAWeighter_SBNNuSyst_",
    "ZExpPCAWeighter_",
    "MECq0q3InterpWeighting_SBN_v3_",
    "MECq0q3InterpWeighting_",
    "CCQEXSecCorr_SBN_v3_",
    "CCQEXSecCorr_",
)

TOTAL_KEYS = frozenset({"genie", "genie_rate", "genie_ar23", "genie_ar23_rate"})


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

def load_cov_mat_dict(path: Path | str) -> dict:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix == ".pkl":
        with open(path, "rb") as f:
            return pickle.load(f)
    if path.suffix == ".npz":
        return _npz_to_cov_mat_dict(path)
    raise ValueError(f"unsupported cov file: {path}")


def _npz_to_cov_mat_dict(path: Path) -> dict:
    """Legacy ``genie-*_syst_dict.npz`` → same flat cov_frac layout as aggregate pkl."""
    z = np.load(path, allow_pickle=True)
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for slug in z.files:
        cell = z[slug].item()
        if not isinstance(cell, dict):
            continue
        row: Dict[str, np.ndarray] = {}
        tot_x = tot_r = None
        for kn, pack in cell.items():
            if not isinstance(pack, dict):
                continue
            for kind, suffix in (("xsec", ""), ("rate", "_rate")):
                if kind not in pack:
                    continue
                sub = pack[kind]
                cf = np.asarray(sub["cov_frac"] if isinstance(sub, dict) else sub, dtype=np.float64)
                row[kn + suffix] = cf
                if kind == "xsec":
                    tot_x = cf.copy() if tot_x is None else tot_x + cf
                else:
                    tot_r = cf.copy() if tot_r is None else tot_r + cf
        if tot_x is not None:
            row["genie"] = tot_x
        if tot_r is not None:
            row["genie_rate"] = tot_r
        out[slug] = row
    return out


def load_groups(sources: Mapping[str, Path | str]) -> Dict[str, dict]:
    """``{mode: cov_mat_dict}`` from ``{mode: path}``."""
    out = {}
    for mode, path in sources.items():
        out[mode] = load_cov_mat_dict(path)
        print(f"loaded {mode}: {path}")
    return out


def _knob_to_mode_lookup() -> Dict[str, str]:
    """Exact knob name → physics mode from ``GENIE_KNOB_GROUPS`` (+ Ar23p / VecFF).

    Skips composite packs (``FSI_compare``, ``slim``) that re-list knobs already
    owned by a physics mode — those would otherwise overwrite CCQE/MEC/… with
    ``FSI_compare``. ``FSI_v1_N`` / ``FSI_v3_N`` map to ``FSI`` only via
    ``setdefault`` so they do not steal knobs already claimed by ``Other``/Ar23p.
    """
    from makedf.geniesyst import GENIE_KNOB_GROUPS, ar23p_genie_systematics

    skip = {"slim", "FSI_compare"}
    pack_dest = {
        "VecFF": "CCQE",
        "ZExp": "CCQE",
        "FSI_v1_N": "FSI",
        "FSI_v3_N": "FSI",
    }
    rev: Dict[str, str] = {}
    for mode, knobs in GENIE_KNOB_GROUPS.items():
        if mode in skip:
            continue
        dest = pack_dest.get(mode, mode)
        soft = mode in ("FSI_v1_N", "FSI_v3_N", "Ar23p")
        for kn in knobs:
            name = str(kn)
            if soft:
                rev.setdefault(name, dest)
            else:
                rev[name] = dest
    for kn in ar23p_genie_systematics:
        rev.setdefault(str(kn), "Ar23p")
    return rev


def assign_knob_to_mode(knob: str, lookup: Optional[Mapping[str, str]] = None) -> str:
    """Assign a knob from a combined ``cov_mat_dict`` to a GENIE mode bucket.

    Order: exact ``GENIE_KNOB_GROUPS`` / Ar23p match → EDepFSI heuristics →
    name-based fallback (ZExp/QE → CCQE, MEC → MEC, …) → ``Other``.
    """
    kn = str(knob)
    table = lookup if lookup is not None else _knob_to_mode_lookup()
    if kn in table:
        return table[kn]
    # rate twin of an exact match
    if kn.endswith("_rate"):
        base = kn[: -len("_rate")]
        if base in table:
            return table[base]

    kn_u = kn.upper()
    if "EDEPFSI" in kn_u:
        if "MEC" in kn_u:
            return "MEC"
        if "QE" in kn_u or "VECFF" in kn_u or "COULOMB" in kn_u:
            return "CCQE"
        return "Other"  # EDepFSI FSI dials (MFP / Fr*)
    if "ZEXP" in kn_u:
        return "CCQE"
    if "MEC" in kn_u:
        return "MEC"
    if "NONRES" in kn_u:
        return "nonRES"
    if "RES" in kn_u:
        return "RES"
    if "DIS" in kn_u:
        return "DIS"
    if "QE" in kn_u or "CRPA" in kn_u or "SF_Q0" in kn_u or "_SF_" in kn_u:
        return "CCQE"
    return "Other"


def split_cov_mat_dict_by_mode(cov: Mapping[str, dict]) -> Dict[str, dict]:
    """Split a combined (all-modes) ``cov_mat_dict`` into ``{mode: cov_mat_dict}``.

    Each mode's per-variable row keeps only that mode's knobs and rebuilds
    ``genie`` / ``genie_rate`` as the sum of those knobs (so mode totals match
    the knobs shown in section 1).
    """
    lookup = _knob_to_mode_lookup()
    modes: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}

    for slug, row in cov.items():
        # Map each non-total key to a mode; totals are rebuilt per mode.
        by_mode_keys: Dict[str, List[str]] = {}
        for key in row:
            if key in TOTAL_KEYS:
                continue
            base = key[: -len("_rate")] if key.endswith("_rate") else key
            mode = assign_knob_to_mode(base, lookup)
            by_mode_keys.setdefault(mode, []).append(key)

        for mode, keys in by_mode_keys.items():
            new_row: Dict[str, np.ndarray] = {
                k: np.asarray(row[k], dtype=np.float64) for k in keys
            }
            tot_x = tot_r = None
            for k, mat in new_row.items():
                if k.endswith("_rate"):
                    tot_r = mat.copy() if tot_r is None else tot_r + mat
                else:
                    tot_x = mat.copy() if tot_x is None else tot_x + mat
            if tot_x is not None:
                new_row["genie"] = tot_x
            if tot_r is not None:
                new_row["genie_rate"] = tot_r
            modes.setdefault(mode, {})[slug] = new_row

    return modes


def load_combined_as_mode_groups(path: Path | str) -> Dict[str, dict]:
    """Load a combined Product-B ``cov_mat_dict.pkl`` and split into mode groups.

    This is the preferred Product-B inspect source: one file already merges
    CCQE/MEC/RES/… plus both nominal and ``EDepFSI_*`` dials.
    """
    path = Path(path)
    cov = load_cov_mat_dict(path)
    groups = split_cov_mat_dict_by_mode(cov)
    print(f"loaded combined→modes from {path}")
    for mode in MODE_ORDER:
        if mode in groups:
            n_kn = sum(
                1
                for k in groups[mode].get("integrated", {})
                if k not in TOTAL_KEYS and not str(k).endswith("_rate")
            )
            print(f"  {mode:7s} vars={len(groups[mode])}  xsec_knobs@integrated={n_kn}")
    extra = sorted(set(groups) - set(MODE_ORDER))
    for mode in extra:
        print(f"  {mode:7s} vars={len(groups[mode])}")
    return groups


def merge_cov_mat_dict_into_groups(
    groups: MutableMapping[str, dict],
    extra: Mapping[str, dict],
    *,
    mode_for_knob: Optional[Callable[[str], str]] = None,
) -> int:
    """Add knobs from ``extra`` into ``groups`` (in place). Returns #keys added.

    Used to fold a parallel single-knob product (e.g. VecFF) or EDepFSI-only
    knobs into New-era mode groups without rebuilding the combined pickle.

    Shape gate: a key is accepted only if its matrix matches an existing row for
    the same ``slug`` in *any* loaded mode (or the destination row). This stops
    May-era EDepFSI matrices (e.g. 50-bin ``vertex_y``) from creating a new
    incompatible slug in a mode that did not already have that variable.
    """
    lookup = _knob_to_mode_lookup()
    assign: Callable[[str], str] = mode_for_knob or (lambda kn: assign_knob_to_mode(kn, lookup))

    def _ref_shape(slug: str, dest_row: Mapping[str, np.ndarray]) -> Optional[Tuple[int, ...]]:
        sample = next(
            (np.asarray(v).shape for k, v in dest_row.items() if k not in TOTAL_KEYS),
            None,
        )
        if sample is not None:
            return sample
        for cov in groups.values():
            row = cov.get(slug)
            if not row:
                continue
            for k, v in row.items():
                if k not in TOTAL_KEYS:
                    return tuple(np.asarray(v).shape)
        return None

    n_added = 0
    n_skipped = 0
    for slug, row in extra.items():
        for key, mat in row.items():
            if key in TOTAL_KEYS:
                continue
            base = key[: -len("_rate")] if key.endswith("_rate") else key
            mode = assign(base)
            dest_row = groups.setdefault(mode, {}).setdefault(slug, {})
            if key in dest_row:
                continue
            arr = np.asarray(mat, dtype=np.float64)
            ref = _ref_shape(slug, dest_row)
            if ref is not None and tuple(arr.shape) != ref:
                n_skipped += 1
                continue
            if ref is None:
                # No New-era coverage for this slug yet — do not introduce
                # Old-era-only binning into the inspect set.
                n_skipped += 1
                continue
            dest_row[key] = arr
            n_added += 1
            tot_key = "genie_rate" if key.endswith("_rate") else "genie"
            if tot_key in dest_row:
                dest_row[tot_key] = dest_row[tot_key] + arr
            else:
                dest_row[tot_key] = arr.copy()
    if n_skipped:
        print(f"  merge skip: {n_skipped} keys (binning mismatch vs existing mode row)")
    return n_added


def sum_cov_fracs(mats: Iterable[np.ndarray]) -> Optional[np.ndarray]:
    """Sum fractional cov matrices; skip shape mismatches with a warning."""
    tot = None
    skipped = 0
    for m in mats:
        arr = np.asarray(m, dtype=np.float64)
        if tot is None:
            tot = arr.copy()
            continue
        if arr.shape != tot.shape:
            skipped += 1
            continue
        tot = tot + arr
    if skipped:
        print(f"  sum_cov_fracs: skipped {skipped} matrices with mismatched shape")
    return tot


# ---------------------------------------------------------------------------
# Knob / matrix accessors
# ---------------------------------------------------------------------------

def short_knob_label(knob: str) -> str:
    kn = str(knob)
    for pfx in _STRIP_PREFIXES:
        if kn.startswith(pfx):
            kn = kn[len(pfx) :]
            break
    for tok in ("_multisim", "_multisigma", "_SBN_v1", "_SBN_v3", "_SBNNuSyst"):
        kn = kn.replace(tok, "")
    while "__" in kn:
        kn = kn.replace("__", "_")
    return kn.strip("_").replace("_", " ")


def display_mode_name(mode: str) -> str:
    """Plot label for a physics mode (CCQE → QE)."""
    return "QE" if mode == "CCQE" else str(mode)


def family_knob_label(knob: str) -> str:
    """Collapse binned / dial / q0bin knobs into one family label.

    Examples: ``ZExp_*_b0..bN`` → ``ZExp``; ``*_dial_N`` → base dial name;
    Martini/Valencia ``*_q0binN`` → ``MEC Martini`` / ``MEC Valencia``.
    """
    kn = str(knob)
    if "MvA" in kn and re.search(r"_b\d+$", kn):
        return "MvA"
    if "ZExp" in kn and re.search(r"_b\d+$", kn):
        return "ZExp"
    m = re.search(r"_dial_\d+$", kn)
    if m:
        return short_knob_label(kn[: m.start()])
    m = re.search(r"_q0bin\d+$", kn)
    if m:
        prefix = kn[: m.start()]
        if "Martini" in prefix:
            return "MEC Martini"
        if "Valencia" in prefix or "Valenica" in prefix:
            return "MEC Valencia"
        return short_knob_label(prefix)
    return short_knob_label(kn)


def _display_knob_path(lab: str) -> str:
    """Rewrite ``CCQE`` → ``QE`` in compound knob / mode labels."""
    s = str(lab)
    if s.startswith("CCQE/") or s.startswith("CCQE "):
        s = "QE" + s[4:]
    s = s.replace("/CCQE/", "/QE/")
    s = s.replace(" CCQE ", " QE ")
    s = s.replace("_CCQE", "_QE")
    if s.endswith(" CCQE"):
        s = s[:-5] + " QE"
    if s.endswith("/CCQE"):
        s = s[:-5] + "/QE"
    # leftover token as whole word
    s = re.sub(r"\bCCQE\b", "QE", s)
    return s


def iter_knob_cov_fracs(row: Mapping[str, np.ndarray], kind: str) -> Iterable[Tuple[str, np.ndarray]]:
    """Yield ``(knob_name, cov_frac)`` for rate or xsec from one variable row."""
    for key, mat in row.items():
        if key in TOTAL_KEYS:
            continue
        if kind == "rate":
            if not key.endswith("_rate"):
                continue
            yield key[: -len("_rate")], np.asarray(mat, dtype=np.float64)
        else:
            if key.endswith("_rate"):
                continue
            yield key, np.asarray(mat, dtype=np.float64)


def mode_total_cov_frac(row: Mapping[str, np.ndarray], kind: str) -> Optional[np.ndarray]:
    key = "genie_rate" if kind == "rate" else "genie"
    if key in row:
        return np.asarray(row[key], dtype=np.float64)
    # fallback: sum knobs
    mats = [m for _k, m in iter_knob_cov_fracs(row, kind)]
    if not mats:
        return None
    tot = mats[0].copy()
    for m in mats[1:]:
        tot = tot + m
    return tot


def frac_unc_pct(cov_frac: np.ndarray) -> np.ndarray:
    c = np.asarray(cov_frac, dtype=np.float64)
    return 100.0 * np.sqrt(np.maximum(np.diag(c), 0.0))


def integrated_frac_unc_pct(cov_frac: np.ndarray) -> float:
    """Scalar ranking metric: frac. unc. on the integrated (first) bin."""
    u = frac_unc_pct(cov_frac)
    return float(u[0]) if u.size else 0.0


def corr_from_cov_frac(cov_frac: np.ndarray) -> np.ndarray:
    c = np.asarray(cov_frac, dtype=np.float64)
    d = np.sqrt(np.maximum(np.diag(c), 0.0))
    denom = np.outer(np.maximum(d, 1e-18), np.maximum(d, 1e-18))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = np.where(denom > 0, c / denom, 0.0)
    np.fill_diagonal(corr, 1.0)
    return np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)


# ---------------------------------------------------------------------------
# Ar23p redistribution
# ---------------------------------------------------------------------------

def is_coh_or_ncel_knob(knob: str) -> bool:
    kn_u = str(knob).upper()
    return "COH" in kn_u or "NCEL" in kn_u


def assign_other_mode_knob_to_bucket(knob: str) -> str:
    """Split the GENIE ``Other`` mode product: COH/NCEL → Other, else → FSI."""
    return "Other" if is_coh_or_ncel_knob(knob) else "FSI"


def assign_ar23p_knob_to_mode(knob: str) -> str:
    """Map an Ar23p knob onto CCQE / MEC / FSI.

    Rules (case-insensitive):
      * name contains ``ZExp`` / ``Zexp`` → CCQE
      * name contains ``MEC`` → MEC
      * name contains ``QE`` → CCQE
      * else → FSI  (non-CCQE/MEC Ar23p knobs)
    """
    kn = str(knob)
    kn_u = kn.upper()
    if "ZEXP" in kn_u:
        return "CCQE"
    if "MEC" in kn_u:
        return "MEC"
    if "QE" in kn_u:
        return "CCQE"
    return "FSI"


def _add_cov(out: MutableMapping[str, np.ndarray], dest: str, mat: np.ndarray) -> None:
    if dest not in out:
        out[dest] = np.asarray(mat, dtype=np.float64).copy()
    else:
        out[dest] = out[dest] + np.asarray(mat, dtype=np.float64)


def mode_totals_ar23p_standalone(
    groups: Mapping[str, dict],
    slug: str,
    kind: str,
) -> Dict[str, np.ndarray]:
    """Per-mode totals with Ar23p kept as its own mode."""
    out: Dict[str, np.ndarray] = {}
    for mode, cov in groups.items():
        if slug not in cov:
            continue
        tot = mode_total_cov_frac(cov[slug], kind)
        if tot is not None:
            out[mode] = tot
    return out


def mode_totals_ar23p_distributed(
    groups: Mapping[str, dict],
    slug: str,
    kind: str,
) -> Dict[str, np.ndarray]:
    """Per-mode totals with Ar23p + Other-mode knobs redistributed.

    * CCQE / MEC / RES / nonRES / DIS: keep each mode's total as-is.
    * Other mode knobs: ``*COH`` / ``*NCEL`` → ``Other``; remaining → ``FSI``.
    * Ar23p knobs: QE/ZExp → CCQE; MEC → MEC; else → ``FSI``.
    """
    out: Dict[str, np.ndarray] = {}
    for mode in ("CCQE", "MEC", "RES", "nonRES", "DIS"):
        if mode not in groups or slug not in groups[mode]:
            continue
        tot = mode_total_cov_frac(groups[mode][slug], kind)
        if tot is not None:
            out[mode] = tot.copy()

    if "Other" in groups and slug in groups["Other"]:
        for knob, mat in iter_knob_cov_fracs(groups["Other"][slug], kind):
            _add_cov(out, assign_other_mode_knob_to_bucket(knob), mat)

    if "Ar23p" in groups and slug in groups["Ar23p"]:
        for knob, mat in iter_knob_cov_fracs(groups["Ar23p"][slug], kind):
            _add_cov(out, assign_ar23p_knob_to_mode(knob), mat)

    return out


# ---------------------------------------------------------------------------
# Collectors for plots
# ---------------------------------------------------------------------------

def knob_parts_for_mode(
    groups: Mapping[str, dict],
    mode: str,
    slug: str,
    kind: str,
    *,
    collapse_binned: bool = False,
) -> Tuple[Dict[str, np.ndarray], Optional[np.ndarray]]:
    """``{short_label: cov_frac}`` + mode total for one mode."""
    if mode not in groups or slug not in groups[mode]:
        return {}, None
    row = groups[mode][slug]
    parts: Dict[str, np.ndarray] = {}
    for knob, mat in iter_knob_cov_fracs(row, kind):
        lab = family_knob_label(knob) if collapse_binned else short_knob_label(knob)
        lab = _display_knob_path(lab)
        parts[lab] = mat.copy() if lab not in parts else parts[lab] + mat
    return parts, mode_total_cov_frac(row, kind)


def all_knob_parts(
    groups: Mapping[str, dict],
    slug: str,
    kind: str,
    *,
    distribute_ar23p: bool = False,
    collapse_binned: bool = True,
) -> Dict[str, np.ndarray]:
    """``{mode/label: cov_frac}`` across all modes (for top-N ranking).

    With ``collapse_binned=True`` (default), binned dial families are summed
    into one entry (e.g. all ``ZExp_b*`` → ``…/ZExp``).
    """
    parts: Dict[str, np.ndarray] = {}
    label_fn = family_knob_label if collapse_binned else short_knob_label
    for mode, cov in groups.items():
        if slug not in cov:
            continue
        if distribute_ar23p and mode == "Other":
            for knob, mat in iter_knob_cov_fracs(cov[slug], kind):
                dest = assign_other_mode_knob_to_bucket(knob)
                lab = _display_knob_path(f"{dest}/{label_fn(knob)}")
                parts[lab] = mat.copy() if lab not in parts else parts[lab] + mat
            continue
        if mode == "Ar23p" and distribute_ar23p:
            for knob, mat in iter_knob_cov_fracs(cov[slug], kind):
                dest = assign_ar23p_knob_to_mode(knob)
                lab = _display_knob_path(f"{dest}/{label_fn(knob)}")
                parts[lab] = mat.copy() if lab not in parts else parts[lab] + mat
            continue
        for knob, mat in iter_knob_cov_fracs(cov[slug], kind):
            lab = _display_knob_path(f"{mode}/{label_fn(knob)}")
            parts[lab] = mat.copy() if lab not in parts else parts[lab] + mat
    return parts


def strip_mode_prefix(lab: str) -> str:
    """``QE/ZExp`` → ``ZExp`` (legend without mode)."""
    s = _display_knob_path(lab)
    if "/" in s:
        return s.split("/", 1)[1]
    return s


# Nominal ↔ EDepFSI twin dials (same GENIE parameter, retired SBNNuSyst naming).
# Keep SBN_v1 / group-product counterparts; do not stack these EDepFSI twins.
EDEPFSI_TWIN_BASES: Tuple[str, ...] = ("NormCCMEC", "DecayAngMEC")
_RETIRED_EDEPFSI_TWIN_RE = re.compile(
    r"EDepFSI_(?:Norm\w*MEC|DecayAngMEC|VecFFCCQEshape|CoulombCCQE)"
)


def is_edepfsi_knob(knob: str) -> bool:
    return "EDepFSI" in str(knob)


def is_retired_edepfsi_twin(knob: str) -> bool:
    """True for EDepFSI copies of dials already covered by SBN_v1 / mode groups.

    Retire: ``EDepFSI_*MEC*`` (NormCCMEC, NormNCMEC, DecayAngMEC), plus the
    same-class QE twins ``EDepFSI_VecFFCCQEshape`` and ``EDepFSI_CoulombCCQE``.
    Keep EDepFSI FSI π/N dials.
    """
    kn = str(knob)
    if kn.endswith("_rate"):
        kn = kn[: -len("_rate")]
    return bool(_RETIRED_EDEPFSI_TWIN_RE.search(kn))


def knob_matches_twin_base(knob: str, base: str, *, edepfsi: bool) -> bool:
    """Match ``…_NormCCMEC`` or ``…_EDepFSI_NormCCMEC`` (and ``*_rate``)."""
    kn = str(knob)
    if kn.endswith("_rate"):
        kn = kn[: -len("_rate")]
    if edepfsi:
        return kn.endswith(f"EDepFSI_{base}") or f"EDepFSI_{base}" in kn
    if "EDepFSI" in kn:
        return False
    return kn.endswith(f"_{base}") or kn.endswith(base)


def find_twin_cov_frac(
    row: Mapping[str, np.ndarray],
    kind: str,
    base: str,
    *,
    edepfsi: bool,
) -> Optional[np.ndarray]:
    """Return the first matching knob cov_frac for ``base`` in ``row``."""
    for kn, mat in iter_knob_cov_fracs(row, kind):
        if knob_matches_twin_base(kn, base, edepfsi=edepfsi):
            return np.asarray(mat, dtype=np.float64)
    return None


def collect_edepfsi_twin_parts(
    groups: Mapping[str, dict],
    slug: str,
    kind: str,
    base: str,
) -> Dict[str, np.ndarray]:
    """``{NormCCMEC: mat, EDepFSI NormCCMEC: mat}`` summed over modes if needed.

    Keys are plot labels (no mode prefix). Missing side is omitted.
    """
    out: Dict[str, np.ndarray] = {}
    label_nom = str(base)
    label_edep = f"EDepFSI {base}"
    for _mode, cov in groups.items():
        row = cov.get(slug)
        if not row:
            continue
        for edep, lab in ((False, label_nom), (True, label_edep)):
            mat = find_twin_cov_frac(row, kind, base, edepfsi=edep)
            if mat is None:
                continue
            if lab in out:
                if out[lab].shape != mat.shape:
                    continue
                out[lab] = out[lab] + mat
            else:
                out[lab] = mat.copy()
    return out


def twin_parts_from_cov_mat_dict(
    cov: Mapping[str, dict],
    slug: str,
    kind: str,
    base: str,
    *,
    mode_filter: Optional[Callable[[str], bool]] = None,
) -> Dict[str, np.ndarray]:
    """Pull nominal + EDepFSI twin mats for ``base`` from a combined cov dict.

    If the dict is already mode-split (``{mode: {slug: row}}``), pass that via
    :func:`collect_edepfsi_twin_parts` instead. This helper expects the flat
    aggregate layout ``{slug: {knob: mat}}`` **or** mode-split groups.
    """
    # Mode-split?
    sample = next(iter(cov.values()), None)
    if isinstance(sample, dict) and sample and isinstance(next(iter(sample.values()), None), dict):
        groups = {
            m: g
            for m, g in cov.items()
            if mode_filter is None or mode_filter(str(m))
        }
        return collect_edepfsi_twin_parts(groups, slug, kind, base)

    row = cov.get(slug)
    if not isinstance(row, dict):
        return {}
    out: Dict[str, np.ndarray] = {}
    for edep, lab in ((False, str(base)), (True, f"EDepFSI {base}")):
        mat = find_twin_cov_frac(row, kind, base, edepfsi=edep)
        if mat is not None:
            out[lab] = mat.copy()
    return out


def top_n_knobs_by_integrated(
    parts: Mapping[str, np.ndarray],
    integrated_parts: Mapping[str, np.ndarray],
    n: int = 10,
) -> List[str]:
    """Order keys by integrated frac. unc. of ``integrated_parts``; keep top ``n``."""
    scores = {k: integrated_frac_unc_pct(m) for k, m in integrated_parts.items()}
    return sorted(scores, key=scores.get, reverse=True)[: max(int(n), 0)]


def contribution_rows(
    parts: Mapping[str, np.ndarray],
    total: Optional[np.ndarray],
    *,
    kind: str,
    slug: str,
    ordered_keys: Optional[Sequence[str]] = None,
    use_integrated_score: bool = True,
) -> List[Dict[str, object]]:
    """Rows for CSV: per-knob (and Total) fractional-unc % used in legends."""
    keys = list(ordered_keys) if ordered_keys is not None else list(parts.keys())
    rows: List[Dict[str, object]] = []
    for lab in keys:
        if lab not in parts:
            continue
        mat = parts[lab]
        score = (
            integrated_frac_unc_pct(mat)
            if use_integrated_score
            else float(np.mean(frac_unc_pct(mat)))
        )
        rows.append(
            {
                "kind": kind,
                "variable": slug,
                "knob": lab,
                "frac_unc_pct": score,
            }
        )
    if total is not None:
        wtot = frac_unc_pct(total)
        score_tot = float(wtot[0]) if use_integrated_score and wtot.size else float(np.mean(wtot))
        rows.append(
            {
                "kind": kind,
                "variable": slug,
                "knob": "Total",
                "frac_unc_pct": score_tot,
            }
        )
    return rows


def write_contribution_csv(path: Path | str, rows: Sequence[Mapping[str, object]]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["kind", "variable", "knob", "frac_unc_pct"]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})
    return path


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _centers_bins(vc, nbin: int):
    if vc is not None and len(getattr(vc, "bin_centers", [])) == nbin:
        centers = np.asarray(vc.bin_centers, float)
        bins = np.asarray(vc.bins, float) if getattr(vc, "bins", None) is not None else None
        if bins is not None and len(bins) == nbin + 1:
            return centers, bins
    return np.arange(nbin, dtype=float), np.arange(nbin + 1, dtype=float) - 0.5


def _bin_range_labels(edges: np.ndarray) -> List[str]:
    return [f"[{edges[i]:.3g}, {edges[i + 1]:.3g})" for i in range(len(edges) - 1)]


def _draw_matrix_on_ax(
    ax,
    matrix: np.ndarray,
    bins: np.ndarray,
    *,
    cmap: str,
    title: str,
    xlab: str,
    ylab: str,
    vmin=None,
    vmax=None,
    annotate: bool = True,
    cbar_scale_only: bool = False,
    show_cbar_label: bool = False,
    cbar_label: str = "",
    tick_labels: Optional[Sequence[str]] = None,
):
    nbins = len(bins)
    if not (nbins - 1 == matrix.shape[0] == matrix.shape[1]):
        raise ValueError(
            f"binning mismatch: len(bins)-1={nbins - 1} vs matrix {matrix.shape} "
            f"(title={title!r})"
        )
    unif = np.linspace(0.0, float(nbins - 1), nbins)
    extent = [unif[0], unif[-1], unif[0], unif[-1]]
    tick_pos = (unif[:-1] + unif[1:]) / 2
    if tick_labels is not None:
        tick_lab = list(tick_labels)
    else:
        tick_lab = _bin_range_labels(np.asarray(bins, float))

    flat = matrix[~np.isnan(matrix)]
    if vmin is None and vmax is None and flat.size:
        vmin = float(np.min(flat))
        vmax = float(np.max(flat))
        if vmin == vmax:
            vmax = vmin + 1e-12
    kw = dict(extent=extent, origin="lower", cmap=cmap, aspect="equal")
    if vmin is not None:
        kw["vmin"] = vmin
    if vmax is not None:
        kw["vmax"] = vmax
    im = ax.imshow(matrix, **kw)

    exponent = 0
    if flat.size > 0:
        nonzero = flat[np.isfinite(flat) & (flat != 0)]
        if nonzero.size:
            example = float(np.max(np.abs(nonzero)))
            logv = np.log10(example)
            if np.isfinite(logv):
                exponent = int(np.floor(logv))
    cbar = plt.colorbar(im, ax=ax, shrink=0.75, fraction=0.046, pad=0.04)
    if cbar_scale_only and exponent not in (0, -1):
        cbar.set_label(f"[10$^{{{exponent}}}$]", fontsize=12)
        cbar.ax.yaxis.set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, _, e=exponent: f"{x / 10**e:.2f}")
        )
    elif show_cbar_label:
        if exponent not in (0, -1):
            cbar.set_label(cbar_label + f" [10$^{{{exponent}}}$]", fontsize=12)
            cbar.ax.yaxis.set_major_formatter(
                mpl.ticker.FuncFormatter(lambda x, _, e=exponent: f"{x / 10**e:.2f}")
            )
        else:
            cbar.set_label(cbar_label, fontsize=12)
    elif cbar_scale_only and exponent in (0, -1):
        # no text label; still format ticks without scale factor
        pass
    else:
        # unlabeled colorbar (e.g. correlation)
        pass

    v0 = float(vmin) if vmin is not None else 0.0
    v1 = float(vmax) if vmax is not None else 1.0
    if annotate:
        for i in range(nbins - 1):
            for j in range(nbins - 1):
                value = matrix[i, j]
                if np.isnan(value):
                    continue
                sig = value / 10**exponent if exponent != -1 else value
                ax.text(
                    j + 0.5, i + 0.5, f"{sig:.2f}",
                    ha="center", va="center",
                    color=_text_color_for(float(value), cmap, v0, v1), fontsize=8,
                )

    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lab, rotation=45, ha="right")
    ax.set_yticks(tick_pos)
    ax.set_yticklabels(tick_lab)
    ax.set_xlabel(xlab, fontsize=14)
    ax.set_ylabel(ylab, fontsize=14)
    ax.set_title(title, fontsize=14)


def plot_frac_unc_breakdown(
    ax,
    parts: Mapping[str, np.ndarray],
    total: Optional[np.ndarray],
    *,
    vc=None,
    title: str = "",
    kind: str = "xsec",
    top_n: Optional[int] = None,
    ordered_keys: Optional[Sequence[str]] = None,
    ylabel: Optional[str] = None,
    show_pct_in_legend: bool = True,
    show_title: bool = False,
    legend_knob_only: bool = False,
    colors: Optional[Mapping[str, object]] = None,
):
    """Step histograms of fractional unc [%] for selected parts + Total.

    ``colors``: optional ``{part_key: matplotlib color}``; falls back to tab20
    by draw order when a key is missing.
    """
    if not parts and total is None:
        if show_title:
            ax.set_title(title + " (empty)")
        return

    if ordered_keys is not None:
        show = [k for k in ordered_keys if k in parts]
    else:
        scores = {k: float(np.mean(frac_unc_pct(m))) for k, m in parts.items()}
        ordered = sorted(scores, key=scores.get, reverse=True)
        show = ordered if top_n is None else ordered[: max(int(top_n), 0)]

    ref = total if total is not None else next(iter(parts.values()))
    nbin = int(np.asarray(ref).shape[0])
    centers, bins = _centers_bins(vc, nbin)
    is_int = getattr(vc, "var_save_name", "") == "integrated" or nbin == 1
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(show), 1)))

    for i, lab in enumerate(show):
        w = frac_unc_pct(parts[lab])
        if is_int:
            w = np.full_like(w, float(w[0]))
        score = integrated_frac_unc_pct(parts[lab]) if is_int else float(np.mean(w))
        disp = strip_mode_prefix(lab) if legend_knob_only else _display_knob_path(lab)
        label = f"{disp} ({score:.2f}%)" if show_pct_in_legend else disp
        color = (
            colors[lab]
            if colors is not None and lab in colors
            else cmap[i % len(cmap)]
        )
        ax.hist(
            centers, bins=bins, weights=w, histtype="step", linewidth=1.3,
            color=color,
            label=label,
        )

    if total is not None:
        wtot = frac_unc_pct(total)
        if is_int:
            wtot = np.full_like(wtot, float(wtot[0]))
        tot_lab = (
            f"Total ({float(np.mean(wtot)):.2f}%)"
            if show_pct_in_legend
            else "Total"
        )
        ax.hist(
            centers, bins=bins, weights=wtot, histtype="step",
            linewidth=2.6, color="k",
            label=tot_lab,
        )

    xlab = vc.var_labels[1] if vc is not None and getattr(vc, "var_labels", None) else ""
    if is_int:
        ax.set_xlabel("All Events")
        ax.set_xticks([centers[0]])
        ax.set_xticklabels(["All Events"])
    else:
        ax.set_xlabel(xlab or "bin")
    ax.set_ylabel(ylabel if ylabel is not None else "Uncertainty [%]")
    if show_title and title:
        ax.set_title(title)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="best", framealpha=0.9)


def plot_mode_frac_unc(
    ax,
    mode_totals: Mapping[str, np.ndarray],
    *,
    vc=None,
    title: str = "",
    kind: str = "xsec",
    mode_order: Sequence[str] = MODE_ORDER,
    show_title: bool = False,
    ylabel: Optional[str] = None,
):
    """Per-mode + Total fractional uncertainty curves."""
    if not mode_totals:
        if show_title:
            ax.set_title(title + " (empty)")
        return
    ref = next(iter(mode_totals.values()))
    nbin = int(np.asarray(ref).shape[0])
    centers, bins = _centers_bins(vc, nbin)
    is_int = getattr(vc, "var_save_name", "") == "integrated" or nbin == 1
    grand = sum_cov_fracs(list(mode_totals.values()))

    for mode in mode_order:
        if mode not in mode_totals:
            continue
        w = frac_unc_pct(mode_totals[mode])
        if is_int:
            w = np.full_like(w, float(w[0]))
        ax.hist(
            centers, bins=bins, weights=w, histtype="step", linewidth=1.8,
            color=MODE_COLORS.get(mode, "gray"),
            label=display_mode_name(mode),
        )
    if grand is not None:
        wtot = frac_unc_pct(grand)
        if is_int:
            wtot = np.full_like(wtot, float(wtot[0]))
        ax.hist(
            centers, bins=bins, weights=wtot, histtype="step",
            linewidth=2.6, color="k", label="Total",
        )

    xlab = vc.var_labels[1] if vc is not None and getattr(vc, "var_labels", None) else ""
    if is_int:
        ax.set_xlabel("All Events")
        ax.set_xticks([centers[0]])
        ax.set_xticklabels(["All Events"])
    else:
        ax.set_xlabel(xlab or "bin")
    ax.set_ylabel(ylabel if ylabel is not None else "Uncertainty [%]")
    if show_title and title:
        ax.set_title(title)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="best")


def show_cov_corr_heatmaps(
    cov_frac: np.ndarray,
    vc,
    *,
    kind: str = "xsec",
    suptitle: Optional[str] = None,
    save_path: Optional[Path | str] = None,
    dpi: int = 140,
    show: bool = True,
    title_prefix: str = "",  # kept for call-site compatibility; unused
):
    """Frac. cov (viridis) + correlation as one wide figure (1×2 subplots)."""
    cf = np.asarray(cov_frac, dtype=np.float64)
    bins = np.asarray(vc.bins, float)
    if cf.shape[0] != len(bins) - 1 or cf.shape[0] != cf.shape[1]:
        print(
            f"skip cov/corr heatmap for {getattr(vc, 'var_save_name', '?')}: "
            f"matrix {cf.shape} vs bins→{len(bins) - 1} "
            f"(suptitle={suptitle!r})"
        )
        return None
    corr = corr_from_cov_frac(cf)
    is_int = getattr(vc, "var_save_name", "") == "integrated" or len(bins) == 2
    if is_int:
        xlab = "All Events"
        tick_labels = ["All Events"]
    else:
        xlab = vc.var_labels[1] if vc.var_labels else vc.var_save_name
        tick_labels = None

    fig, axes = plt.subplots(1, 2, figsize=(18.0, 7.5))
    _draw_matrix_on_ax(
        axes[0], cf, bins,
        cmap="viridis",
        title="Fractional Covariance",
        xlab=xlab, ylab=xlab,
        cbar_scale_only=True,
        tick_labels=tick_labels,
    )
    _draw_matrix_on_ax(
        axes[1], corr, bins,
        cmap="bwr",
        title="Correlation",
        xlab=xlab, ylab=xlab,
        vmin=-1, vmax=1,
        cbar_scale_only=False,
        show_cbar_label=False,
        tick_labels=tick_labels,
    )
    if suptitle:
        fig.suptitle(suptitle, fontsize=14, y=1.02)
    fig.tight_layout()
    if save_path is not None:
        out = Path(save_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print("wrote", out)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig
