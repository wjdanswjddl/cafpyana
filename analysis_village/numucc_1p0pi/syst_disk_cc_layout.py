"""On-disk layout for **cross-variable (joint-bin)** systematic covariances.

These live in a tree **parallel** to the marginal ``syst_disk`` tree (see
``analysis_village.numucc_1p0pi.syst_disk_layout``), rooted at
``NUMUCC_SYST_DISK_CC_ROOT`` or :func:`dataset_locations.default_syst_disk_cc_root`.

Joint multisim outputs are written **per category** under ``JointMCstat/``, ``JointFlux/``,
and ``JointG4/`` (each with its own ``joint_*_combined.npz``). A legacy single bundle
``JointMultisim/joint_multisim_combined.npz`` may still exist on older disks; :mod:`cc_joint_cov`
prefers that file when present, otherwise sums the per-category files.

Bin indexing convention for a kinematic pair ``(var_X, var_Y)`` (e.g. proton vs muon
kinematics in the conditional-constraint study):

* Global joint index ``k = 0 .. n_X + n_Y - 1`` maps **X bins first**, then **Y bins**::

      k in [0, n_X)           → variable X, bin k
      k in [n_X, n_X + n_Y) → variable Y, bin (k - n_X)

The symmetric covariance ``C`` has shape ``(n_X + n_Y, n_X + n_Y)`` with::

      C[0:n_X, 0:n_X]           = Σ_XX   (marginal over X bins)
      C[n_X:, n_X:]             = Σ_YY
      C[0:n_X, n_X:] = C[n_X:, 0:n_X].T = Σ_XY / Σ_YX

Preset **pair_slug** strings (e.g. ``muon_p__proton_costheta``) match the conditional
notebook; stacked bins are **X (proton) then Y (muon)** — see
:mod:`analysis_village.numucc_1p0pi.syst_cc_joint_multisim_common`.
"""

from __future__ import annotations

import os

SYST_DISK_CC_ENV = "NUMUCC_SYST_DISK_CC_ROOT"

# Legacy single-tree layout (still read by :mod:`cc_joint_cov` if present).
SUB_JOINT_MULTISIM = "JointMultisim"
FILE_JOINT_MULTISIM_COMBINED = "joint_multisim_combined.npz"

# Per-neutrino-multisim category (parallel to each other under ``syst_disk_CC``).
SUB_JOINT_MCSTAT = "JointMCstat"
FILE_JOINT_MCSTAT_COMBINED = "joint_mcstat_combined.npz"
SUB_JOINT_FLUX = "JointFlux"
FILE_JOINT_FLUX_COMBINED = "joint_flux_combined.npz"
SUB_JOINT_G4 = "JointG4"
FILE_JOINT_G4_COMBINED = "joint_g4_combined.npz"

SUB_JOINT_GENIE = "JointGenie"
FILE_JOINT_GENIE_COMBINED = "joint_genie_combined.npz"

_JOINT_MULTISIM_CAT: dict[str, tuple[str, str]] = {
    "MCstat": (SUB_JOINT_MCSTAT, FILE_JOINT_MCSTAT_COMBINED),
    "Flux": (SUB_JOINT_FLUX, FILE_JOINT_FLUX_COMBINED),
    "G4": (SUB_JOINT_G4, FILE_JOINT_G4_COMBINED),
}


def normalized_root(root: str) -> str:
    return os.path.abspath(os.path.expanduser(root.rstrip(os.sep)))


def syst_disk_cc_paths(root: str) -> dict[str, str]:
    """Absolute paths under a conditional-constraint syst disk root."""
    r = normalized_root(root)
    legacy_ms = os.path.join(r, SUB_JOINT_MULTISIM, FILE_JOINT_MULTISIM_COMBINED)
    return {
        "root": r,
        "joint_multisim_combined": legacy_ms,
        "joint_multisim_legacy": legacy_ms,
        "joint_multisim_mcstat": os.path.join(r, SUB_JOINT_MCSTAT, FILE_JOINT_MCSTAT_COMBINED),
        "joint_multisim_flux": os.path.join(r, SUB_JOINT_FLUX, FILE_JOINT_FLUX_COMBINED),
        "joint_multisim_g4": os.path.join(r, SUB_JOINT_G4, FILE_JOINT_G4_COMBINED),
        "joint_genie_combined": os.path.join(r, SUB_JOINT_GENIE, FILE_JOINT_GENIE_COMBINED),
    }


def joint_multisim_category_inner_key(category: str) -> str:
    """NPZ cell inner key for one neutrino multisim category (``JointFlux``, ``JointG4``, ``JointMCstat``)."""
    if category not in _JOINT_MULTISIM_CAT:
        raise ValueError("unknown joint multisim category %r" % category)
    return "Joint%s" % category


def joint_multisim_category_out_dir(root: str, category: str) -> str:
    if category not in _JOINT_MULTISIM_CAT:
        raise ValueError("unknown joint multisim category %r" % category)
    sub, _fname = _JOINT_MULTISIM_CAT[category]
    return os.path.join(normalized_root(root), sub)


def joint_multisim_category_npz_basename(category: str) -> str:
    """Basename only (e.g. ``joint_flux_combined.npz``)."""
    if category not in _JOINT_MULTISIM_CAT:
        raise ValueError("unknown joint multisim category %r" % category)
    return _JOINT_MULTISIM_CAT[category][1]


def joint_multisim_out_dir(root: str) -> str:
    """Deprecated: legacy ``JointMultisim/`` directory. Prefer :func:`joint_multisim_category_out_dir`."""
    return os.path.join(normalized_root(root), SUB_JOINT_MULTISIM)


def joint_genie_out_dir(root: str) -> str:
    return os.path.join(normalized_root(root), SUB_JOINT_GENIE)
