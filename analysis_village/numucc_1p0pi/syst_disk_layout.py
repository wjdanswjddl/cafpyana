"""On-disk layout for precomputed systematic covariances consumed by ``utils.get_syst_unc``.

Set ``NUMUCC_SYST_DISK_ROOT`` to a directory with **exactly** this structure (one subdirectory per
systematic **source**):

.. code-block:: text

    <SYST_DISK_ROOT>/
      MCstat/mcstat_syst_dict.npz
      Flux/flux_syst_dict.npz
      G4/g4_syst_dict.npz
      GENIE/cov_mat_dict.pkl
      Cosmics/cosmics_syst_dict.npz
      Detector/detector_syst_dict.npz

Producer scripts write into these paths; loaders **fail** if any expected file is missing.
"""

from __future__ import annotations

import os

SYST_DISK_ENV = "NUMUCC_SYST_DISK_ROOT"

SUB_MCSTAT = "MCstat"
SUB_FLUX = "Flux"
SUB_G4 = "G4"
SUB_GENIE = "GENIE"
SUB_COSMICS = "Cosmics"
SUB_DETECTOR = "Detector"

FILE_MCSTAT = "mcstat_syst_dict.npz"
FILE_FLUX = "flux_syst_dict.npz"
FILE_G4 = "g4_syst_dict.npz"
FILE_GENIE = "cov_mat_dict.pkl"
FILE_COSMICS = "cosmics_syst_dict.npz"
FILE_DETECTOR = "detector_syst_dict.npz"


def normalized_root(root: str) -> str:
    return os.path.abspath(os.path.expanduser(root.rstrip(os.sep)))


def syst_disk_paths(root: str) -> dict[str, str]:
    """Map logical keys to absolute paths (includes ``\"root\"`` for the resolved tree root)."""
    r = normalized_root(root)
    return {
        "root": r,
        "mcstat": os.path.join(r, SUB_MCSTAT, FILE_MCSTAT),
        "flux": os.path.join(r, SUB_FLUX, FILE_FLUX),
        "g4": os.path.join(r, SUB_G4, FILE_G4),
        "genie": os.path.join(r, SUB_GENIE, FILE_GENIE),
        "cosmics": os.path.join(r, SUB_COSMICS, FILE_COSMICS),
        "detector": os.path.join(r, SUB_DETECTOR, FILE_DETECTOR),
    }


def category_out_dir(root: str, category: str) -> str:
    """Directory for a producer category (``MCstat``, ``Flux``, …)."""
    return os.path.join(normalized_root(root), category)
