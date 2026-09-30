"""On-disk layout for **cross-variable (joint-bin)** systematic covariances.

These live under ``JointCC/`` on the PRL Product **B** consumer tree (see
:func:`dataset_locations.default_syst_disk_cc_root`), or ``NUMUCC_SYST_DISK_CC_ROOT``.

Default production writes one inclusive stacked NPZ cell ``stacked_mu_p``
(muon *p*, muon cosθ, proton *p*, proton cosθ) under ``JointMCstat/``,
``JointFlux/``, ``JointG4/``, and ``JointGenie/``. :mod:`cc_joint_cov` permutes
that stack to the analysis order ``[X; Y_1; Y_2; …]``. Pairwise cells from older
campaigns may still exist; the stacked loader does not stitch them.

Bin indexing for the stacked cell follows ``meta["var_save_names"]`` /
``meta["n_bins_per_var"]`` (concatenated inclusive selected-rate histograms,
``bkgd_subtract=False``).
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

# Unisim extras (detector / cosmics / POT / ntargets) — stacked rank-1 joints.
SUB_JOINT_DETECTOR = "JointDetector"
FILE_JOINT_DETECTOR_COMBINED = "joint_detector_combined.npz"
SUB_JOINT_COSMICS = "JointCosmics"
FILE_JOINT_COSMICS_COMBINED = "joint_cosmics_combined.npz"
SUB_JOINT_POT = "JointPOT"
FILE_JOINT_POT_COMBINED = "joint_pot_combined.npz"
SUB_JOINT_NTARGETS = "JointNtargets"
FILE_JOINT_NTARGETS_COMBINED = "joint_ntargets_combined.npz"

_JOINT_MULTISIM_CAT: dict[str, tuple[str, str]] = {
    "MCstat": (SUB_JOINT_MCSTAT, FILE_JOINT_MCSTAT_COMBINED),
    "Flux": (SUB_JOINT_FLUX, FILE_JOINT_FLUX_COMBINED),
    "G4": (SUB_JOINT_G4, FILE_JOINT_G4_COMBINED),
}

_JOINT_EXTRAS_CAT: dict[str, tuple[str, str]] = {
    "detector": (SUB_JOINT_DETECTOR, FILE_JOINT_DETECTOR_COMBINED),
    "cosmics": (SUB_JOINT_COSMICS, FILE_JOINT_COSMICS_COMBINED),
    "pot": (SUB_JOINT_POT, FILE_JOINT_POT_COMBINED),
    "ntargets": (SUB_JOINT_NTARGETS, FILE_JOINT_NTARGETS_COMBINED),
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
        "joint_detector": os.path.join(r, SUB_JOINT_DETECTOR, FILE_JOINT_DETECTOR_COMBINED),
        "joint_cosmics": os.path.join(r, SUB_JOINT_COSMICS, FILE_JOINT_COSMICS_COMBINED),
        "joint_pot": os.path.join(r, SUB_JOINT_POT, FILE_JOINT_POT_COMBINED),
        "joint_ntargets": os.path.join(r, SUB_JOINT_NTARGETS, FILE_JOINT_NTARGETS_COMBINED),
    }


def joint_extras_out_dir(root: str, category: str) -> str:
    if category not in _JOINT_EXTRAS_CAT:
        raise ValueError("unknown joint extras category %r" % category)
    return os.path.join(normalized_root(root), _JOINT_EXTRAS_CAT[category][0])


def joint_extras_npz_basename(category: str) -> str:
    if category not in _JOINT_EXTRAS_CAT:
        raise ValueError("unknown joint extras category %r" % category)
    return _JOINT_EXTRAS_CAT[category][1]


def joint_extras_inner_key(category: str) -> str:
    keys = {
        "detector": "JointDetector",
        "cosmics": "JointCosmics",
        "pot": "JointPOT",
        "ntargets": "JointNtargets",
    }
    if category not in keys:
        raise ValueError("unknown joint extras category %r" % category)
    return keys[category]


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
