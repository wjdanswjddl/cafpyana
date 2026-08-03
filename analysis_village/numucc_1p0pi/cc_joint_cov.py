"""Load **joint** (cross-variable) systematic covariances from ``syst_disk_CC``.

This module exposes two complementary builders:

* :func:`build_joint_covariance_abs` (legacy) — single ``(var_X, var_Y)`` pair, stacked bin
  layout ``[X; Y]`` of size ``n_X + n_Y``.

* :func:`build_joint_multi_covariance_abs` (new) — arbitrary ``var_X`` constrained by an
  ordered list ``var_Ys = [Y_1, Y_2, …]``. Stacked bins ``[X; Y_1; Y_2; …]`` of size
  ``n_X + Σ_i n_{Y_i}``. Required NPZ pairs are auto-discovered: ``(X, Y_i)`` for the
  diagonal/X-Y cross blocks, and ``(Y_i, Y_j)`` for the Y-Y cross blocks. Pairs are read
  from ``JointMCstat/`` ``JointFlux/`` ``JointG4/`` (legacy ``JointMultisim/`` accepted) plus
  optional ``JointGenie/`` and oriented with :func:`syst_cc_joint_multisim_common.joint_pair_layout`
  so the same NPZ works regardless of how the ``(A, B)`` order was stored on disk.

Per-category ``joint_*_combined.npz`` files may include optional ``Joint{Flux,G4,MCstat}_by_knob``
breakdowns. :func:`load_joint_multisim_frac_cov` sums ``cov_frac`` across selected categories
(by default all files on disk), using each category's **combined** block unless a per-category
knob list is supplied. Legacy ``JointMultisim/joint_multisim_combined.npz`` is still read when
present (indivisible single block).

Bin layout matches :mod:`analysis_village.numucc_1p0pi.syst_disk_cc_layout`: indices
``0 .. n_X-1`` are variable **X** bins, ``n_X .. n_X+n_Y-1`` are variable **Y** bins.
"""

from __future__ import annotations

import os
from collections.abc import Mapping, Sequence

import numpy as np

from analysis_village.numucc_1p0pi.dataset_locations import default_syst_disk_cc_root
from analysis_village.numucc_1p0pi.syst_cc_joint_multisim_common import (
    joint_multisim_npz_pair_key,
    joint_pair_layout,
    reorient_joint_block,
)
from analysis_village.numucc_1p0pi.syst_disk_cc_layout import (
    joint_multisim_category_inner_key,
    syst_disk_cc_paths,
)
from analysis_village.numucc_1p0pi.utils import get_syst_unc as get_syst_unc_disk
from pyanalib.covariance import cov_from_fraccov

JOINT_MULTISIM_CATEGORY_ORDER: tuple[str, ...] = ("MCstat", "Flux", "G4")

_JOINT_MULTISIM_PATH_KEYS: dict[str, str] = {
    "MCstat": "joint_multisim_mcstat",
    "Flux": "joint_multisim_flux",
    "G4": "joint_multisim_g4",
}


def resolve_syst_disk_cc_root(explicit: str | os.PathLike[str] | None) -> str:
    root = explicit or os.environ.get("NUMUCC_SYST_DISK_CC_ROOT")
    if not root:
        root = str(default_syst_disk_cc_root())
    return syst_disk_cc_paths(str(root))["root"]


_CROSS_DISK_TO_MULTISIM_NAME: dict[str, str] = {
    "flux": "Flux",
    "g4": "G4",
    "genie": "GENIE",
}
_MS_DISK_TO_JOINT_CATEGORY: dict[str, str] = {
    "mcstat": "MCstat",
    "flux": "Flux",
    "g4": "G4",
}
_JOINT_CATEGORY_TO_MARGINAL_SKIP: dict[str, str] = {
    "MCstat": "mcstat",
    "Flux": "flux",
    "G4": "g4",
}


def resolve_primary_syst_disk_joint_options(
    primary_syst_disk_keys: Sequence[str] | None,
    syst_cc_root: str | os.PathLike[str],
    *,
    joint_multisim_categories: Sequence[str] | None = None,
    joint_multisim_exclude_categories: Sequence[str] = ("MCstat",),
) -> dict[str, object]:
    """Map ``PRIMARY_SYST_DISK_KEYS`` to joint CC loader flags and legacy XY names.

    Auto-detects whether the joint NPZ path is usable for the requested keys (i.e.
    whether any joint GENIE or joint multisim NPZ exists on disk for those keys).
    The caller should use ``use_joint_multisim_cc`` from the returned dict as the
    master flag controlling the joint path — no manual ``USE_JOINT_MULTISIM_CC``
    variable is needed.

    When ``primary_syst_disk_keys`` is a subset, joint multisim loading is limited to
    the requested ``mcstat`` / ``flux`` / ``g4`` disk keys. The indivisible legacy
    ``JointMultisim/joint_multisim_combined.npz`` is only used when
    ``primary_syst_disk_keys is None`` (load everything on disk).
    """
    paths = syst_disk_cc_paths(resolve_syst_disk_cc_root(syst_cc_root))
    exc = frozenset(joint_multisim_exclude_categories)

    if primary_syst_disk_keys is None:
        pk: tuple[str, ...] | None = None
        pk_set: frozenset[str] | None = None
    else:
        pk = tuple(str(x).lower() for x in primary_syst_disk_keys)
        pk_set = frozenset(pk)

    if pk is None:
        syst_for_cross_cov: tuple[str, ...] = ("Flux", "G4", "GENIE")
    else:
        x = tuple(_CROSS_DISK_TO_MULTISIM_NAME[k] for k in pk if k in _CROSS_DISK_TO_MULTISIM_NAME)
        syst_for_cross_cov = x if x else ("Flux", "G4", "GENIE")

    if pk is None:
        include_genie = True
        include_multisim = True
        ms_order = (
            tuple(joint_multisim_categories)
            if joint_multisim_categories is not None
            else JOINT_MULTISIM_CATEGORY_ORDER
        )
    else:
        include_genie = "genie" in pk_set
        include_multisim = bool(pk_set & frozenset(_MS_DISK_TO_JOINT_CATEGORY))
        if joint_multisim_categories is not None:
            ms_order = tuple(
                c
                for c in joint_multisim_categories
                if c in JOINT_MULTISIM_CATEGORY_ORDER
                and _JOINT_CATEGORY_TO_MARGINAL_SKIP[c] in pk_set
                and c not in exc
            )
        else:
            ms_order = tuple(
                cat
                for cat in JOINT_MULTISIM_CATEGORY_ORDER
                if _JOINT_CATEGORY_TO_MARGINAL_SKIP[cat] in pk_set and cat not in exc
            )

    legacy_exists = os.path.isfile(paths["joint_multisim_legacy"])
    use_legacy = (
        legacy_exists
        and pk is None
        and joint_multisim_categories is None
        and include_multisim
    )

    if use_legacy:
        cats_effective: tuple[str, ...] = ()
    elif include_multisim:
        cats_effective = tuple(
            c
            for c in ms_order
            if c not in exc and os.path.isfile(paths[_JOINT_MULTISIM_PATH_KEYS[c]])
        )
    else:
        cats_effective = ()

    include_multisim_eff = include_multisim and (use_legacy or len(cats_effective) > 0)

    # Auto-detect: the joint path is usable if at least one joint NPZ exists for the
    # requested keys. Check GENIE NPZ existence when genie is requested.
    genie_npz_exists = os.path.isfile(paths["joint_genie_combined"])
    include_genie_eff = include_genie and genie_npz_exists
    use_joint_multisim_cc = include_multisim_eff or include_genie_eff

    skip: set[str] = set()
    if include_multisim_eff:
        if use_legacy:
            skip.update(_JOINT_CATEGORY_TO_MARGINAL_SKIP.values())
        else:
            skip.update(_JOINT_CATEGORY_TO_MARGINAL_SKIP[c] for c in cats_effective)
    if include_genie_eff:
        skip.add("genie")

    if pk is None:
        marginal: tuple[str, ...] | None = None
    else:
        marginal = tuple(k for k in pk if k not in skip)

    return {
        "use_joint_multisim_cc": use_joint_multisim_cc,
        "syst_for_cross_cov": syst_for_cross_cov,
        "syst_unc_components": primary_syst_disk_keys,
        "include_multisim": include_multisim_eff,
        "include_genie": include_genie_eff,
        "multisim_categories_effective": cats_effective,
        "use_legacy_multisim_combined": use_legacy,
        "marginal_syst_diag_components": marginal,
    }


def _multisim_frac_from_category_cell(
    cell: dict,
    category: str,
    knob_names: Sequence[str] | None,
) -> np.ndarray:
    """Fractional joint covariance for one multisim category from one NPZ *pair_slug* cell."""
    inner = joint_multisim_category_inner_key(category)
    if knob_names is None:
        return np.asarray(cell[inner]["cov_frac"], dtype=float)
    by_key = "%s_by_knob" % inner
    if by_key not in cell:
        raise KeyError(
            "category %r requested knob subset %s but %r is missing from NPZ cell (only combined %r is available)"
            % (category, list(knob_names), by_key, inner)
        )
    bk = cell[by_key]
    if not isinstance(bk, dict):
        bk = dict(bk)
    out = None
    for k in knob_names:
        if k not in bk:
            raise KeyError("knob %r not in %s (have: %s)" % (k, by_key, ", ".join(sorted(bk.keys()))))
        tri = bk[k]
        frac = np.asarray(tri["cov_frac"], dtype=float)
        out = frac if out is None else out + frac
    if out is None:
        raise ValueError("empty knob_names for category %r" % category)
    return out


def load_joint_multisim_frac_cov(
    var_X,
    var_Y,
    syst_cc_root: str | os.PathLike[str] | None = None,
    *,
    categories: Sequence[str] | None = None,
    knobs_by_category: Mapping[str, Sequence[str] | None] | None = None,
) -> tuple[np.ndarray, dict]:
    """Return fractional covariance ``(n_X+n_Y, n_X+n_Y)`` and metadata for one kinematic pair.

    * **Legacy** ``JointMultisim/joint_multisim_combined.npz``: loaded when the file exists.
      Per-category or per-knob selection is **not** supported; pass ``categories=None`` and
      ``knobs_by_category=None`` only.

    * **Per-category** ``JointMCstat/``, ``JointFlux/``, ``JointG4/``: sums ``cov_frac`` across
      categories. ``categories=None`` includes every category file that exists on disk (order:
      MCstat, Flux, G4). ``knobs_by_category`` maps a category name to ``None`` (use the NPZ
      **combined** block for that category) or to a sequence of knob names (Flux/G4 only;
      sums those knobs' ``cov_frac`` from ``Joint*_by_knob``).

    The returned matrix is **oriented to match the caller**: even if the NPZ stores the joint
    block as ``[var_Y bins; var_X bins]``, this function transparently re-permutes to
    ``[var_X bins; var_Y bins]`` (see :func:`joint_pair_layout`).
    """
    root = resolve_syst_disk_cc_root(syst_cc_root)
    paths = syst_disk_cc_paths(root)
    pair_slug, swap, n_first, n_second = joint_pair_layout(var_X, var_Y)
    legacy = paths["joint_multisim_legacy"]

    want_filter = categories is not None or (
        knobs_by_category is not None and len(knobs_by_category) > 0
    )
    if os.path.isfile(legacy):
        if want_filter:
            raise ValueError(
                "legacy %s is indivisible: use categories=None and knobs_by_category=None, "
                "or remove the legacy file and use per-category JointMCstat/JointFlux/JointG4 NPZs"
                % legacy
            )
        blob = np.load(legacy, allow_pickle=True)
        if pair_slug not in blob.files:
            raise KeyError(
                "pair_slug %r not in %s (available: %s)"
                % (pair_slug, legacy, ", ".join(sorted(blob.files)))
            )
        cell = dict(blob)[pair_slug].item()
        pack = cell["JointMultisim"]
        frac = np.asarray(pack["cov_frac"], dtype=float)
        frac = reorient_joint_block(frac, swap=swap, n_first=n_first, n_second=n_second)
        meta = dict(cell.get("meta") or {})
        meta["joint_multisim_layout"] = "legacy_JointMultisim"
        return frac, meta

    cat_order = tuple(categories) if categories is not None else JOINT_MULTISIM_CATEGORY_ORDER
    fracs: list[np.ndarray] = []
    meta: dict = {
        "joint_multisim_layout": "per_category",
        "joint_multisim_categories": [],
        "joint_multisim_knob_mode": [],
    }
    for cat in cat_order:
        if cat not in _JOINT_MULTISIM_PATH_KEYS:
            raise ValueError("unknown joint multisim category %r (allowed: %s)" % (cat, JOINT_MULTISIM_CATEGORY_ORDER))
        pth = paths[_JOINT_MULTISIM_PATH_KEYS[cat]]
        if not os.path.isfile(pth):
            if categories is not None:
                raise FileNotFoundError("joint multisim category %r required but missing: %s" % (cat, pth))
            continue
        blob = np.load(pth, allow_pickle=True)
        if pair_slug not in blob.files:
            raise KeyError(
                "pair_slug %r not in %s (available: %s)"
                % (pair_slug, pth, ", ".join(sorted(blob.files)))
            )
        cell = dict(blob)[pair_slug].item()
        knob_sel: Sequence[str] | None = None
        if knobs_by_category is not None and cat in knobs_by_category:
            knob_sel = knobs_by_category[cat]
            if knob_sel is not None and len(knob_sel) == 0:
                raise ValueError("knobs_by_category[%r] is an empty sequence; use None for combined" % cat)
        frac = _multisim_frac_from_category_cell(cell, cat, knob_sel)
        frac = reorient_joint_block(frac, swap=swap, n_first=n_first, n_second=n_second)
        fracs.append(frac)
        meta["joint_multisim_categories"].append(cat)
        meta["joint_multisim_knob_mode"].append("combined" if knob_sel is None else list(knob_sel))
    if not fracs:
        raise FileNotFoundError(
            "Joint multisim covariance not found under %s (expected legacy %s or one of "
            "JointMCstat/joint_mcstat_combined.npz, JointFlux/joint_flux_combined.npz, "
            "JointG4/joint_g4_combined.npz)"
            % (root, legacy)
        )
    total = np.sum(fracs, axis=0)
    return total, meta


def load_joint_genie_frac_cov(
    var_X,
    var_Y,
    syst_cc_root: str | os.PathLike[str] | None = None,
    *,
    knob_names: Sequence[str] | None = None,
) -> tuple[np.ndarray, dict] | None:
    """Return ``(cov_frac, meta)`` for joint GENIE (**rate** universes only), or ``None`` if absent.

    ``knob_names=None`` uses the combined ``JointGenie`` block. Otherwise sums ``cov_frac`` from
    ``JointGenie_by_knob`` for each listed reweight knob (independent-knob recipe used at write time).

    Re-orients the returned matrix to ``[var_X; var_Y]`` regardless of the on-disk layout, like
    :func:`load_joint_multisim_frac_cov`.
    """
    root = resolve_syst_disk_cc_root(syst_cc_root)
    pth = syst_disk_cc_paths(root)["joint_genie_combined"]
    if not os.path.isfile(pth):
        return None
    pair_slug, swap, n_first, n_second = joint_pair_layout(var_X, var_Y)
    blob = np.load(pth, allow_pickle=True)
    if pair_slug not in blob.files:
        return None
    cell = dict(blob)[pair_slug].item()
    meta = dict(cell.get("meta") or {})
    meta["joint_genie_cov_source"] = "rate_universes"
    if knob_names is None:
        pack = cell["JointGenie"]
        frac = np.asarray(pack["cov_frac"], dtype=float)
        frac = reorient_joint_block(frac, swap=swap, n_first=n_first, n_second=n_second)
        meta["joint_genie_knob_mode"] = "combined"
        return frac, meta
    bk_key = "JointGenie_by_knob"
    if bk_key not in cell:
        raise KeyError(
            "joint GENIE knob subset requested but %r missing from NPZ cell for pair %r" % (bk_key, pair_slug)
        )
    bk = cell[bk_key]
    if not isinstance(bk, dict):
        bk = dict(bk)
    out = None
    for k in knob_names:
        if k not in bk:
            raise KeyError(
                "GENIE knob %r not in %s (have: %s)" % (k, bk_key, ", ".join(sorted(bk.keys())))
            )
        tri = bk[k]
        frac = np.asarray(tri["cov_frac"], dtype=float)
        out = frac if out is None else out + frac
    if out is None:
        raise ValueError("empty knob_names for joint GENIE")
    out = reorient_joint_block(out, swap=swap, n_first=n_first, n_second=n_second)
    meta["joint_genie_knob_mode"] = list(knob_names)
    return out, meta


def inspect_joint_cc_disk(
    var_X,
    var_Y,
    syst_cc_root: str | os.PathLike[str] | None = None,
) -> dict:
    """Summarize which joint multisim / GENIE blocks exist on disk and list optional *by_knob* keys."""
    root = resolve_syst_disk_cc_root(syst_cc_root)
    paths = syst_disk_cc_paths(root)
    pair_slug = joint_multisim_npz_pair_key(var_X, var_Y)
    out: dict = {"root": root, "pair_slug": pair_slug, "multisim": {}, "genie": {}}
    for cat in JOINT_MULTISIM_CATEGORY_ORDER:
        pth = paths[_JOINT_MULTISIM_PATH_KEYS[cat]]
        if not os.path.isfile(pth):
            out["multisim"][cat] = {"path": pth, "present": False}
            continue
        blob = np.load(pth, allow_pickle=True)
        if pair_slug not in blob.files:
            out["multisim"][cat] = {"path": pth, "present": True, "pair_in_file": False}
            continue
        cell = dict(blob)[pair_slug].item()
        inner = joint_multisim_category_inner_key(cat)
        by_key = "%s_by_knob" % inner
        knobs = sorted(dict(cell[by_key]).keys()) if by_key in cell else []
        out["multisim"][cat] = {
            "path": pth,
            "present": True,
            "pair_in_file": True,
            "knobs": knobs,
        }
    gpath = paths["joint_genie_combined"]
    if not os.path.isfile(gpath):
        out["genie"] = {"path": gpath, "present": False}
    elif pair_slug not in np.load(gpath, allow_pickle=True).files:
        out["genie"] = {"path": gpath, "present": True, "pair_in_file": False}
    else:
        cell = dict(np.load(gpath, allow_pickle=True))[pair_slug].item()
        bk = cell.get("JointGenie_by_knob") or {}
        if not isinstance(bk, dict):
            bk = dict(bk)
        knobs = sorted(bk.keys()) if bk else []
        out["genie"] = {"path": gpath, "present": True, "pair_in_file": True, "knobs": knobs}
    return out


def build_joint_covariance_abs(
    var_X,
    var_Y,
    mu_X: np.ndarray,
    mu_Y: np.ndarray,
    syst_cc_root: str | os.PathLike[str] | None = None,
    syst_marginal_root: str | os.PathLike[str] | None = None,
    marginal_syst_components: tuple[str, ...] | None = None,
    *,
    joint_multisim_categories: Sequence[str] | None = None,
    joint_multisim_knobs_by_category: Mapping[str, Sequence[str] | None] | None = None,
    use_joint_multisim: bool = True,
    use_joint_genie: bool = True,
    joint_genie_knobs: Sequence[str] | None = None,
    marginal_genie_cov_frac_key: str = "genie_rate",
) -> np.ndarray:
    """Absolute-count joint covariance ``Σ`` of stacked bins ``[X; Y]``.

    Thin wrapper over :func:`build_joint_multi_covariance_abs` with ``var_Ys=(var_Y,)``. See
    that function for the full multi-Y semantics. Kept for backward compatibility with the
    original single-Y workflow.
    """
    sigma_multi = build_joint_multi_covariance_abs(
        var_X,
        (var_Y,),
        mu_X,
        (mu_Y,),
        syst_cc_root=syst_cc_root,
        syst_marginal_root=syst_marginal_root,
        marginal_syst_components=marginal_syst_components,
        joint_multisim_categories=joint_multisim_categories,
        joint_multisim_knobs_by_category=joint_multisim_knobs_by_category,
        use_joint_multisim=use_joint_multisim,
        use_joint_genie=use_joint_genie,
        joint_genie_knobs=joint_genie_knobs,
        marginal_genie_cov_frac_key=marginal_genie_cov_frac_key,
    )
    return sigma_multi


def _xx_yy_xy_from_pair_frac(
    frac_oriented: np.ndarray,
    n_X: int,
    n_Y: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split an oriented ``(n_X+n_Y)``-square fractional joint cov into ``(frac_XX, frac_XY, frac_YY)``."""
    frac_XX = frac_oriented[0:n_X, 0:n_X]
    frac_XY = frac_oriented[0:n_X, n_X : n_X + n_Y]
    frac_YY = frac_oriented[n_X : n_X + n_Y, n_X : n_X + n_Y]
    return frac_XX, frac_XY, frac_YY


def _frac_pair_oriented(
    var_A,
    var_B,
    *,
    syst_cc_root: str | os.PathLike[str] | None,
    use_joint_multisim: bool,
    joint_multisim_categories: Sequence[str] | None,
    joint_multisim_knobs_by_category: Mapping[str, Sequence[str] | None] | None,
    use_joint_genie: bool,
    joint_genie_knobs: Sequence[str] | None,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Per-pair fractional joint cov, oriented as ``[A; B]``. Returns ``(frac_multisim, frac_genie)``.

    Either entry may be ``None`` when the corresponding source is disabled or not present on
    disk. Both sources are *separately* fractional with respect to the same CV used at write
    time; consumers must rescale them with the analysis ``μ`` stack before combining (this is
    what :func:`build_joint_multi_covariance_abs` does).
    """
    paths = syst_disk_cc_paths(resolve_syst_disk_cc_root(syst_cc_root))
    p_ms_legacy = paths["joint_multisim_legacy"]
    have_ms_on_disk = os.path.isfile(p_ms_legacy) or any(
        os.path.isfile(paths[k]) for k in _JOINT_MULTISIM_PATH_KEYS.values()
    )

    frac_ms: np.ndarray | None = None
    if use_joint_multisim and have_ms_on_disk:
        frac_ms, _meta = load_joint_multisim_frac_cov(
            var_A,
            var_B,
            syst_cc_root=syst_cc_root,
            categories=joint_multisim_categories,
            knobs_by_category=joint_multisim_knobs_by_category,
        )

    frac_g: np.ndarray | None = None
    if use_joint_genie:
        pack = load_joint_genie_frac_cov(
            var_A,
            var_B,
            syst_cc_root=syst_cc_root,
            knob_names=joint_genie_knobs,
        )
        if pack is not None:
            frac_g, _gmeta = pack

    return frac_ms, frac_g


def build_joint_multi_covariance_abs(
    var_X,
    var_Ys: Sequence,
    mu_X: np.ndarray,
    mu_Ys: Sequence[np.ndarray],
    syst_cc_root: str | os.PathLike[str] | None = None,
    syst_marginal_root: str | os.PathLike[str] | None = None,
    marginal_syst_components: tuple[str, ...] | None = None,
    *,
    joint_multisim_categories: Sequence[str] | None = None,
    joint_multisim_knobs_by_category: Mapping[str, Sequence[str] | None] | None = None,
    use_joint_multisim: bool = True,
    use_joint_genie: bool = True,
    joint_genie_knobs: Sequence[str] | None = None,
    marginal_genie_cov_frac_key: str = "genie_rate",
) -> np.ndarray:
    """Absolute joint covariance for stacked vector ``[X; Y_1; Y_2; …]``.

    **This is one joint model, not a sequence of separate single-Y constraints.** The map
    pipeline writes **pairwise** NPZs (each is the covariance of stacked histograms
    ``[A; B]`` in one GENIE/multisim universe family). Here those blocks are **recombined**
    into the unique big matrix ``Σ`` for     ``[X; Y_1; Y_2; …]``: cross blocks ``Σ_{X Y_i}``
    and ``Σ_{Y_i Y_j}`` come from the corresponding pair NPZs, and ``Σ_{Y_i Y_j}`` for
    ``i≠j`` is essential so that conditioning uses **one** inverted ``Σ_{YY}`` on the
    **concatenated** data vector ``n_Y = [n_{Y_1}; n_{Y_2}; …]`` (same formula as single-Y,
    larger blocks). Pairwise shards exist because the chunk step only histograms two
    variables at a time; they are not themselves separate published constraints.

    Per the MicroBooNE data-driven model-validation prescription (Abratenko *et al.*,
    "Data-driven model validation for neutrino-nucleus cross section measurements"), the full
    joint covariance ``Σ`` is assembled from independent contributions:

    * **Joint multisim** (Flux / G4 / MCstat) — per-category per-pair fractional covariances
      from :func:`load_joint_multisim_frac_cov`. Disabled with ``use_joint_multisim=False``.

    * **Joint GENIE** (rate universes) — per-pair fractional covariance from
      :func:`load_joint_genie_frac_cov`. Disabled with ``use_joint_genie=False``.
      ``joint_genie_knobs`` selects a subset of reweight knobs from ``JointGenie_by_knob``.

    * **Marginal disk** — block-diagonal augmentation from :func:`utils.get_syst_unc` applied
      independently to each variable's diagonal block (the joint NPZ tree does **not** encode
      flux/G4/GENIE detector/cosmics/etc. cross blocks beyond what is in the joint multisim /
      GENIE files). GENIE is **omitted** from the marginal pass when joint GENIE is loaded
      so the diagonal GENIE term is not double-counted (matches the legacy single-Y
      :func:`build_joint_covariance_abs` semantics).

    Each pair ``(X, Y_i)`` and ``(Y_i, Y_j)`` is loaded from its preset *pair_slug* NPZ and
    re-oriented as ``[A; B]`` so the cross block can be placed into the correct off-diagonal
    location of the stacked ``Σ``. Per-variable diagonal blocks are averaged across all NPZ
    pairs that touch that variable (each variable typically appears in multiple pair NPZs;
    fractional covariances for the *same* variable should be identical up to MC-stat noise,
    so averaging is a deliberate stabilization).

    Parameters
    ----------
    var_X : VariableConfig
        Target / constrained variable.
    var_Ys : Sequence[VariableConfig]
        Ordered list of constraining variables. May be length 1 (single-Y, equivalent to
        :func:`build_joint_covariance_abs`).
    mu_X : np.ndarray
        Central-value MC prediction for ``X``, length ``n_X``.
    mu_Ys : Sequence[np.ndarray]
        Central-value MC predictions for each ``Y_i``, same order as ``var_Ys``.

    Returns
    -------
    sigma : np.ndarray
        ``(n_tot, n_tot)`` symmetric covariance, ``n_tot = n_X + Σ_i n_{Y_i}``. Bin layout
        ``[X bins; Y_1 bins; Y_2 bins; …]``.
    """
    from analysis_village.numucc_1p0pi.utils import SYST_UNC_DISK_KEYS

    if marginal_genie_cov_frac_key not in ("genie", "genie_rate"):
        raise ValueError("marginal_genie_cov_frac_key must be 'genie' or 'genie_rate'")
    if not use_joint_multisim and not use_joint_genie:
        raise ValueError("build_joint_multi_covariance_abs: enable at least one of use_joint_multisim, use_joint_genie")

    var_Ys = tuple(var_Ys)
    mu_Ys = tuple(np.asarray(m, dtype=float) for m in mu_Ys)
    if len(var_Ys) != len(mu_Ys):
        raise ValueError(
            "len(var_Ys)=%d does not match len(mu_Ys)=%d" % (len(var_Ys), len(mu_Ys))
        )
    if len(var_Ys) == 0:
        raise ValueError("build_joint_multi_covariance_abs: var_Ys must be non-empty")

    mu_X = np.asarray(mu_X, dtype=float)
    n_X = mu_X.size
    n_Ys = [int(m.size) for m in mu_Ys]
    n_tot = n_X + int(sum(n_Ys))

    # Block-offset table: 0 = X, i = Y_i (1-indexed into var_Ys)
    offsets = [0]
    sizes = [n_X] + list(n_Ys)
    for s in sizes[:-1]:
        offsets.append(offsets[-1] + s)
    # offsets[k] = starting row/col of block k; sizes[k] = block size

    mu_blocks = [mu_X] + list(mu_Ys)
    var_blocks = [var_X] + list(var_Ys)
    mu_stacked_all = np.concatenate(mu_blocks)
    if mu_stacked_all.size != n_tot:
        raise AssertionError("internal stacked length mismatch")

    sigma = np.zeros((n_tot, n_tot), dtype=float)

    # Track marginal-diag augmentation eligibility per pair (GENIE skipped when joint GENIE found)
    has_joint_genie_any = False
    applied_any_multisim = False

    # Build absolute cov by iterating over UNORDERED variable pairs (A,B) with A_idx <= B_idx in our stacking.
    # Each pair NPZ gives us, after orientation as [A; B]:
    #   frac_AA, frac_AB, frac_BB  (each may come from multisim, GENIE, or both, all wrt CV at write time)
    # We rescale frac_AA -> cov_AA with μ_A from analysis, similarly for AB and BB. This matches what
    # build_joint_covariance_abs has always done in the single-Y case.
    # Diagonal blocks may be visited multiple times if a variable appears in more than one pair NPZ;
    # we accumulate the *sum across pairs* per source and divide by the visit count at the end so
    # the per-variable diagonal contribution from joint NPZs stays self-consistent.
    diag_sums_ms = [np.zeros((sz, sz), dtype=float) for sz in sizes]
    diag_sums_g = [np.zeros((sz, sz), dtype=float) for sz in sizes]
    diag_counts_ms = [0 for _ in sizes]
    diag_counts_g = [0 for _ in sizes]

    n_blocks = len(sizes)
    for a in range(n_blocks):
        for b in range(a + 1, n_blocks):
            var_a = var_blocks[a]
            var_b = var_blocks[b]
            mu_a = mu_blocks[a]
            mu_b = mu_blocks[b]
            na = sizes[a]
            nb = sizes[b]
            mu_pair = np.concatenate([mu_a, mu_b])

            frac_ms, frac_g = _frac_pair_oriented(
                var_a,
                var_b,
                syst_cc_root=syst_cc_root,
                use_joint_multisim=use_joint_multisim,
                joint_multisim_categories=joint_multisim_categories,
                joint_multisim_knobs_by_category=joint_multisim_knobs_by_category,
                use_joint_genie=use_joint_genie,
                joint_genie_knobs=joint_genie_knobs,
            )

            if frac_ms is not None:
                applied_any_multisim = True
                cov_pair = cov_from_fraccov(frac_ms, mu_pair)
                cov_AA, cov_AB, cov_BB = _xx_yy_xy_from_pair_frac(cov_pair, na, nb)
                # Off-diagonal cross block placed at (a, b) and symmetrically (b, a)
                _place_block(sigma, cov_AB, offsets[a], offsets[b])
                _place_block(sigma, cov_AB.T, offsets[b], offsets[a])
                diag_sums_ms[a] += cov_AA
                diag_sums_ms[b] += cov_BB
                diag_counts_ms[a] += 1
                diag_counts_ms[b] += 1

            if frac_g is not None:
                has_joint_genie_any = True
                cov_pair_g = cov_from_fraccov(frac_g, mu_pair)
                cov_AA_g, cov_AB_g, cov_BB_g = _xx_yy_xy_from_pair_frac(cov_pair_g, na, nb)
                sigma[offsets[a] : offsets[a] + na, offsets[b] : offsets[b] + nb] += cov_AB_g
                sigma[offsets[b] : offsets[b] + nb, offsets[a] : offsets[a] + na] += cov_AB_g.T
                diag_sums_g[a] += cov_AA_g
                diag_sums_g[b] += cov_BB_g
                diag_counts_g[a] += 1
                diag_counts_g[b] += 1

    # Diagonal accumulation: average the per-variable contribution across the pair NPZs that
    # touched it. (Each variable appears in (n_blocks - 1) pairs at most.)
    for k in range(n_blocks):
        if diag_counts_ms[k] > 0:
            sigma_block = diag_sums_ms[k] / float(diag_counts_ms[k])
            sigma[offsets[k] : offsets[k] + sizes[k], offsets[k] : offsets[k] + sizes[k]] += sigma_block
        if diag_counts_g[k] > 0:
            sigma_block = diag_sums_g[k] / float(diag_counts_g[k])
            sigma[offsets[k] : offsets[k] + sizes[k], offsets[k] : offsets[k] + sizes[k]] += sigma_block

    if not applied_any_multisim and not has_joint_genie_any:
        paths = syst_disk_cc_paths(resolve_syst_disk_cc_root(syst_cc_root))
        raise FileNotFoundError(
            "No joint syst_disk_CC block found for any (X, Y_i) / (Y_i, Y_j) pair: joint multisim "
            "disabled or missing and joint GENIE absent (%s)" % paths["joint_genie_combined"]
        )

    if syst_marginal_root is None:
        return sigma

    # Marginal block-diagonal augmentation (no off-diagonal contribution by construction).
    if marginal_syst_components is None:
        skip = {"mcstat", "flux", "g4"}
        if has_joint_genie_any:
            skip = skip | {"genie"}
        marginal_syst_components = tuple(k for k in SYST_UNC_DISK_KEYS if k not in skip) + (
            "pot",
            "ntargets",
        )

    for k, vk in enumerate(var_blocks):
        _, frac_kk = get_syst_unc_disk(
            vk,
            syst_disk_root=str(syst_marginal_root),
            syst_components=marginal_syst_components,
            genie_cov_frac_key=marginal_genie_cov_frac_key,
        )
        sigma[
            offsets[k] : offsets[k] + sizes[k], offsets[k] : offsets[k] + sizes[k]
        ] += cov_from_fraccov(frac_kk, mu_blocks[k])
    return sigma


def _place_block(target: np.ndarray, block: np.ndarray, row0: int, col0: int) -> None:
    nr, nc = block.shape
    target[row0 : row0 + nr, col0 : col0 + nc] += block


def split_joint_sigma(
    sigma_joint: np.ndarray, n_x: int, n_y: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split single-Y ``Σ`` into ``(Σ_XX, Σ_XY, Σ_YY)`` given ``n_x``, ``n_y``."""
    sigma_joint = np.asarray(sigma_joint, dtype=float)
    s_xx = sigma_joint[0:n_x, 0:n_x]
    s_xy = sigma_joint[0:n_x, n_x : n_x + n_y]
    s_yy = sigma_joint[n_x : n_x + n_y, n_x : n_x + n_y]
    return s_xx, s_xy, s_yy


def split_joint_multi_sigma(
    sigma_joint: np.ndarray, n_X: int, n_Ys: Sequence[int]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[tuple[int, int]]]:
    """Split multi-Y ``Σ`` (layout ``[X; Y_1; Y_2; …]``) into ``(Σ_XX, Σ_XY, Σ_YY, slices_Y)``.

    Returns
    -------
    sigma_XX : np.ndarray
        ``(n_X, n_X)`` block.
    sigma_XY : np.ndarray
        ``(n_X, n_Y_total)`` cross block where ``n_Y_total = Σ_i n_{Y_i}``. Columns are
        stacked in the input order of ``n_Ys``.
    sigma_YY : np.ndarray
        ``(n_Y_total, n_Y_total)`` symmetric block containing all ``Y_i × Y_j`` correlations.
    slices_Y : list of (start, stop)
        Per-``Y_i`` slice within the stacked Y vector (``stop`` is exclusive). Useful for
        slicing ``n_Y`` / ``μ_Y`` per-variable plots out of the stacked conditioning vector.
    """
    sigma_joint = np.asarray(sigma_joint, dtype=float)
    n_Ys = [int(n) for n in n_Ys]
    n_Y_total = int(sum(n_Ys))
    n_tot = n_X + n_Y_total
    if sigma_joint.shape != (n_tot, n_tot):
        raise ValueError(
            "split_joint_multi_sigma: shape %s != (%d, %d) (n_X=%d, n_Y_total=%d)"
            % (sigma_joint.shape, n_tot, n_tot, n_X, n_Y_total)
        )
    s_xx = sigma_joint[0:n_X, 0:n_X]
    s_xy = sigma_joint[0:n_X, n_X : n_X + n_Y_total]
    s_yy = sigma_joint[n_X : n_X + n_Y_total, n_X : n_X + n_Y_total]
    slices_Y: list[tuple[int, int]] = []
    cur = 0
    for nyi in n_Ys:
        slices_Y.append((cur, cur + nyi))
        cur += nyi
    return s_xx, s_xy, s_yy, slices_Y


def inspect_joint_multi_cc_disk(
    var_X,
    var_Ys: Sequence,
    syst_cc_root: str | os.PathLike[str] | None = None,
) -> dict:
    """Summarize joint NPZ availability for every pair required by a multi-Y constraint.

    Returns a dict mapping ``"(varA, varB)"`` pair labels to the same per-pair structure as
    :func:`inspect_joint_cc_disk`. Handy for verifying that all needed NPZs are on disk before
    calling :func:`build_joint_multi_covariance_abs`.
    """
    var_blocks = [var_X] + list(var_Ys)
    out: dict = {}
    for a in range(len(var_blocks)):
        for b in range(a + 1, len(var_blocks)):
            label = "%s__%s" % (var_blocks[a].var_save_name, var_blocks[b].var_save_name)
            try:
                out[label] = inspect_joint_cc_disk(var_blocks[a], var_blocks[b], syst_cc_root=syst_cc_root)
            except KeyError as ex:
                out[label] = {"error": str(ex)}
    return out
