"""Shared helpers for **joint-bin** (cross-variable) multisim covariance production.

See :mod:`analysis_village.numucc_1p0pi.syst_disk_cc_layout` for directory names
(``JointMCstat/``, ``JointFlux/``, ``JointG4/``, legacy ``JointMultisim/``, ``JointGenie/``).
"""

from __future__ import annotations

import os
from typing import Sequence

import numpy as np

# -----------------------------------------------------------------------------
# Chunk pickle names (map phase under Flux/G4/MCstat work trees)
# -----------------------------------------------------------------------------
# Use ``nu__joint_cc__`` / ``nu__joint_cc_genie__`` so ``skip_existing`` in the parallel
# drivers does not reuse older ``nu__joint__*`` shards built before the kinematic-pair list
# grew (e.g. same-side pairs for multi-Y). Aggregate globs match these prefixes only.
JOINT_CC_MULTISIM_CHUNK_PREFIX = "nu__joint_cc__"
JOINT_CC_MULTISIM_CHUNK_GLOB = JOINT_CC_MULTISIM_CHUNK_PREFIX + "*.pkl"

JOINT_CC_GENIE_CHUNK_PREFIX = "nu__joint_cc_genie__"
JOINT_CC_GENIE_CHUNK_GLOB = JOINT_CC_GENIE_CHUNK_PREFIX + "*.pkl"


def joint_cc_multisim_chunk_basename(syst_tag: str, df_stem: str) -> str:
    """Basename for one joint multisim map shard, e.g. ``nu__joint_cc__Flux__merged_0001.pkl``."""
    return "%s%s__%s.pkl" % (JOINT_CC_MULTISIM_CHUNK_PREFIX, syst_tag, df_stem)


def joint_cc_genie_chunk_basename(genie_group: str, df_stem: str) -> str:
    """Basename for one joint GENIE map shard, e.g. ``nu__joint_cc_genie__MaCCQE__merged_0001.pkl``."""
    return "%s%s__%s.pkl" % (JOINT_CC_GENIE_CHUNK_PREFIX, genie_group, df_stem)

from analysis_village.numucc_1p0pi.syst_disk_cc_layout import (
    FILE_JOINT_GENIE_COMBINED,
    joint_genie_out_dir,
    joint_multisim_category_inner_key,
    joint_multisim_category_npz_basename,
    joint_multisim_category_out_dir,
    normalized_root,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig


def kinematic_pair_slug(var_X, var_Y) -> str:
    """Stable key from ``var_save_name`` values: ``var_X`` then ``var_Y`` (not the preset CSV slug)."""
    return "%s__%s" % (var_X.var_save_name, var_Y.var_save_name)


def joint_multisim_npz_pair_key(var_X, var_Y) -> str:
    """Preset *pair_slug* key shared by per-category ``joint_*_combined.npz`` files under ``syst_disk_CC``.

    Lookup is **unordered** on ``var_save_name`` — the same pair_slug is returned whether the
    caller passes ``(A, B)`` or ``(B, A)``. The internal NPZ bin layout follows the tuple
    order in :func:`default_kinematic_joint_pairs`; if you need to re-orient the joint block
    so that ``A`` comes first, use :func:`joint_pair_layout` and transpose when ``swap_needed``.
    """
    names = {var_X.var_save_name, var_Y.var_save_name}
    for slug, _vx, _vy in default_kinematic_joint_pairs():
        if {_vx.var_save_name, _vy.var_save_name} == names:
            return slug
    raise KeyError(
        "no preset joint multisim pair for var_save_names %r; extend default_kinematic_joint_pairs()"
        % (sorted(names),)
    )


def joint_pair_layout(var_A, var_B) -> tuple[str, bool, int, int]:
    """Return ``(pair_slug, swap_needed, n_first, n_second)`` for variables ``(A, B)``.

    The NPZ for ``pair_slug`` stores a ``(n_first + n_second, n_first + n_second)`` joint
    covariance with bins ordered as ``[var_X_npz bins; var_Y_npz bins]`` (see the tuple in
    :func:`default_kinematic_joint_pairs`). ``swap_needed=True`` means the on-disk layout is
    ``(B, A)``; consumers wishing to view the cov in ``[A; B]`` order must apply the
    permutation ``[range(n_B, n_B + n_A), range(0, n_B)]`` (e.g. via
    :func:`reorient_joint_block`).
    """
    name_a = var_A.var_save_name
    name_b = var_B.var_save_name
    for slug, vx, vy in default_kinematic_joint_pairs():
        if vx.var_save_name == name_a and vy.var_save_name == name_b:
            return slug, False, len(vx.bin_centers), len(vy.bin_centers)
        if vx.var_save_name == name_b and vy.var_save_name == name_a:
            return slug, True, len(vx.bin_centers), len(vy.bin_centers)
    raise KeyError(
        "no preset joint multisim pair for var_save_names %r; extend default_kinematic_joint_pairs()"
        % (sorted({name_a, name_b}),)
    )


def reorient_joint_block(block: "np.ndarray", *, swap: bool, n_first: int, n_second: int) -> "np.ndarray":
    """Re-orient a ``(n_first+n_second)``-square joint block ``[var_X_npz; var_Y_npz]`` to ``[A; B]``.

    Used together with :func:`joint_pair_layout`. When ``swap=False`` the matrix is returned
    unchanged; when ``swap=True`` the rows / cols are permuted so the **A** block (of size
    ``n_second``) sits first and the **B** block (``n_first``) sits second.
    """
    import numpy as _np

    arr = _np.asarray(block, dtype=float)
    n_tot = n_first + n_second
    if arr.shape != (n_tot, n_tot):
        raise ValueError(
            "reorient_joint_block: expected shape (%d,%d), got %s" % (n_tot, n_tot, arr.shape)
        )
    if not swap:
        return arr
    perm = _np.concatenate([_np.arange(n_first, n_first + n_second), _np.arange(0, n_first)])
    return arr[_np.ix_(perm, perm)]


def default_kinematic_joint_pairs() -> tuple[tuple[str, object, object], ...]:
    """Kinematic pairs for joint multisim covariance production.

    Each row is ``(preset_slug, var_X, var_Y)`` where ``(var_X, var_Y)`` defines the **on-disk
    bin ordering** for the pair (NPZ stacks ``[var_X bins; var_Y bins]``). The first four pairs
    correspond to the original conditional-constraint pipeline (X = proton / constrained,
    Y = muon / constraining). The last two are **same-side** pairs needed so that a multi-variable
    Y constraint can assemble the full ``Y_i × Y_j`` (and ``X_i × X_j``) cross blocks of the joint
    covariance ``Σ``. With those present, :func:`cc_joint_cov.build_joint_multi_covariance_abs`
    can build a stacked ``[X, Y_1, Y_2, …]`` covariance for arbitrary subsets of muon and proton
    kinematic observables.
    """
    return (
        ("muon_p__proton_p", VariableConfig.proton_momentum(), VariableConfig.muon_momentum()),
        (
            "muon_costheta__proton_costheta",
            VariableConfig.proton_direction(),
            VariableConfig.muon_direction(),
        ),
        ("muon_p__proton_costheta", VariableConfig.proton_direction(), VariableConfig.muon_momentum()),
        ("muon_costheta__proton_p", VariableConfig.proton_momentum(), VariableConfig.muon_direction()),
        # Same-side pairs added for multi-variable Y conditional constraints:
        # both muons → builds Σ(Y_i, Y_j) cross block; both protons → enables X-multi if desired.
        (
            "muon_p__muon_costheta",
            VariableConfig.muon_momentum(),
            VariableConfig.muon_direction(),
        ),
        (
            "proton_p__proton_costheta",
            VariableConfig.proton_momentum(),
            VariableConfig.proton_direction(),
        ),
    )


def parse_pair_slugs_csv(spec: str | None) -> tuple[str, ...] | None:
    """Parse comma-separated pair slugs, or ``None`` when ``spec`` is blank → all defaults."""
    if not spec or not str(spec).strip():
        return None
    raw = [x.strip() for x in str(spec).split(",") if x.strip()]
    allowed = {t[0] for t in default_kinematic_joint_pairs()}
    unk = [x for x in raw if x not in allowed]
    if unk:
        raise ValueError("unknown pair slug(s): %s; allowed: %s" % (unk, sorted(allowed)))
    return tuple(raw)


def select_joint_pairs(pair_slugs: Sequence[str] | None) -> tuple[tuple[str, object, object], ...]:
    allp = default_kinematic_joint_pairs()
    if pair_slugs is None:
        return allp
    want = frozenset(pair_slugs)
    return tuple(p for p in allp if p[0] in want)


def save_joint_multisim_category_npz(
    per_pair: dict[str, dict],
    syst_disk_cc_root: str,
    category: str,
) -> str:
    """Write ``Joint<Category>/<joint_<category>_combined>.npz`` with keys = *pair_slug*.

    Each cell holds the **category-only** combined covariance under ``JointFlux`` / ``JointG4`` /
    ``JointMCstat``, plus optional ``Joint*_by_knob`` for Flux/G4 knob-nested breakdowns.
    """
    root = normalized_root(syst_disk_cc_root)
    inner = joint_multisim_category_inner_key(category)
    by_knob_inner = "%s_by_knob" % inner
    outd = joint_multisim_category_out_dir(root, category)
    os.makedirs(outd, exist_ok=True)
    out_path = os.path.join(outd, joint_multisim_category_npz_basename(category))
    payload: dict[str, np.ndarray] = {}
    for pair_slug, pack in per_pair.items():
        cell = {
            inner: {
                "cov_frac": np.asarray(pack["cov_frac"], dtype=float),
                "cov": np.asarray(pack["cov"], dtype=float),
                "corr": np.asarray(pack["corr"], dtype=float),
            },
            "meta": pack.get("meta") or {},
        }
        bk = pack.get("by_knob")
        if bk:
            cell[by_knob_inner] = bk
        payload[pair_slug] = np.array(cell, dtype=object)
    np.savez_compressed(out_path, **payload)
    return out_path


def save_joint_genie_combined_npz(per_pair: dict[str, dict], syst_disk_cc_root: str) -> str:
    """Write ``JointGenie/joint_genie_combined.npz`` (same cell layout as multisim, inner key ``JointGenie``).

    When ``by_knob`` is present on a pair pack, also writes ``JointGenie_by_knob``:
    ``{knob_name: {cov, cov_frac, corr}}`` for each GENIE reweight knob (one group per knob).
    """
    root = normalized_root(syst_disk_cc_root)
    outd = joint_genie_out_dir(root)
    os.makedirs(outd, exist_ok=True)
    out_path = os.path.join(outd, FILE_JOINT_GENIE_COMBINED)
    payload: dict[str, np.ndarray] = {}
    for pair_slug, pack in per_pair.items():
        cell = {
            "JointGenie": {
                "cov_frac": np.asarray(pack["cov_frac"], dtype=float),
                "cov": np.asarray(pack["cov"], dtype=float),
                "corr": np.asarray(pack["corr"], dtype=float),
            },
            "meta": pack.get("meta") or {},
        }
        bk = pack.get("by_knob")
        if bk:
            cell["JointGenie_by_knob"] = bk
        payload[pair_slug] = np.array(cell, dtype=object)
    np.savez_compressed(out_path, **payload)
    return out_path


def joint_meta(nx: int, ny: int, var_x: str, var_y: str) -> dict:
    return {
        "n_bins_X": int(nx),
        "n_bins_Y": int(ny),
        "n_bins_total": int(nx + ny),
        "index_order": "X_bins_indices_0_to_nx_minus_1_then_Y_bins",
        "var_X": var_x,
        "var_Y": var_y,
    }
