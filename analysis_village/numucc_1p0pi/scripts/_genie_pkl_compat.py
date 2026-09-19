"""Load cov_mat_dict pickles written under either numpy 1.x or 2.x.

numpy 2 renamed ``numpy.core`` to ``numpy._core``, so pickles made there name
modules this 1.26 venv does not have. Remapping in ``find_class`` keeps the
rename local to unpickling; aliasing ``sys.modules`` instead breaks the pandas
C-extension import check.
"""
from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np


class _NumpyCompatUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str):
        try:
            return super().find_class(module, name)
        except ModuleNotFoundError:
            if not module.startswith("numpy._core"):
                raise
            return super().find_class("numpy.core" + module[len("numpy._core"):], name)


def load_cov_mat_dict_compat(path: Path | str) -> dict:
    path = Path(path)
    with open(path, "rb") as f:
        return _NumpyCompatUnpickler(f).load()
