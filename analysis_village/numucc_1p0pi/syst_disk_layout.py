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
      CategorySummary/category_syst_summary.npz   (from ``systematics-summary.ipynb``)

Producer scripts write into these paths; loaders **fail** if any expected file is missing.
"""

from __future__ import annotations

import os
import pickle
from typing import Any

SYST_DISK_ENV = "NUMUCC_SYST_DISK_ROOT"

SUB_MCSTAT = "MCstat"
SUB_FLUX = "Flux"
SUB_G4 = "G4"
SUB_GENIE = "GENIE"
SUB_COSMICS = "Cosmics"
SUB_DETECTOR = "Detector"
SUB_CATEGORY_SUMMARY = "CategorySummary"

FILE_MCSTAT = "mcstat_syst_dict.npz"
FILE_FLUX = "flux_syst_dict.npz"
FILE_G4 = "g4_syst_dict.npz"
FILE_GENIE = "cov_mat_dict.pkl"
FILE_COSMICS = "cosmics_syst_dict.npz"
FILE_DETECTOR = "detector_syst_dict.npz"
FILE_CATEGORY_SUMMARY = "category_syst_summary.npz"
FILE_CATEGORY_SUMMARY_MANIFEST = "category_syst_summary_manifest.json"


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
        "category_summary": category_summary_npz_path(r),
    }


def category_summary_dir(root: str) -> str:
    return os.path.join(normalized_root(root), SUB_CATEGORY_SUMMARY)


def category_summary_npz_path(root: str) -> str:
    return os.path.join(category_summary_dir(root), FILE_CATEGORY_SUMMARY)


def category_summary_manifest_path(npz_path: str) -> str:
    d = os.path.dirname(os.path.abspath(npz_path))
    return os.path.join(d, FILE_CATEGORY_SUMMARY_MANIFEST)


def category_out_dir(root: str, category: str) -> str:
    """Directory for a producer category (``MCstat``, ``Flux``, …)."""
    return os.path.join(normalized_root(root), category)


def resolve_genie_disk_path(root: str) -> str | None:
    """Return the on-disk GENIE payload path (pickle or ``savez_compressed`` sidecar)."""
    r = normalized_root(root)
    base = os.path.join(r, SUB_GENIE, FILE_GENIE)
    for candidate in (base, f"{base}.npz"):
        if os.path.isfile(candidate):
            return candidate
    return None


def load_genie_disk_payload(root: str) -> Any | None:
    """Load ``GENIE/cov_mat_dict.pkl`` (pickle) or the ``.npz`` written by multisim notebooks."""
    genie_path = resolve_genie_disk_path(root)
    if genie_path is None:
        return None
    if genie_path.endswith(".npz"):
        import numpy as np

        return np.load(genie_path, allow_pickle=True)
    try:
        with open(genie_path, "rb") as gf:
            return pickle.load(gf)
    except Exception:
        import numpy as np

        return np.load(genie_path, allow_pickle=True)
