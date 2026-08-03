#!/usr/bin/env python
"""Replace ``integrated`` in ``systematics-final/GENIE/cov_mat_dict.pkl`` from a notebook run.

Use after re-running ``systematics-genie.ipynb`` with ``today_str = "integrated"`` and the
fixed ``get_systematics_genie`` xsec path (response-matrix accumulators, not rate copy).

Example::

    python merge_integrated_genie_into_final.py \\
        --src /exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-notebook-genie-integrated/GENIE/cov_mat_dict.pkl
"""
from __future__ import annotations

import argparse
import pickle
import shutil
import sys
from datetime import datetime
from os import path
from pathlib import Path

import numpy as np

sys.path.append(path.join(path.dirname(__file__), "..", "..", ".."))

from analysis_village.numucc_1p0pi.files_config import save_fig_base_dir  # noqa: E402
from analysis_village.numucc_1p0pi.syst_disk_layout import FILE_GENIE, SUB_GENIE, category_out_dir  # noqa: E402


def merge_integrated_row(dst: dict, src: dict, *, slug: str = "integrated") -> None:
    if slug not in src:
        raise KeyError(f"source missing {slug!r}; keys={list(src)}")
    if slug not in dst:
        raise KeyError(f"destination missing {slug!r}; keys={list(dst)[:5]}...")
    new_row = {}
    for key, val in src[slug].items():
        arr = np.asarray(val, dtype=np.float64)
        if arr.ndim != 2:
            raise ValueError(f"src[{slug!r}][{key!r}] not 2D: {arr.shape}")
        new_row[key] = arr.copy()
    dst[slug] = new_row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--src",
        type=Path,
        default=Path(save_fig_base_dir)
        / "systematics-notebook-genie-integrated"
        / "GENIE"
        / FILE_GENIE,
        help="Integrated GENIE cov_mat_dict.pkl from the notebook run",
    )
    parser.add_argument(
        "--dst",
        type=Path,
        default=Path(save_fig_base_dir) / "systematics-final" / "GENIE" / FILE_GENIE,
        help="Final syst-disk GENIE pickle to update",
    )
    parser.add_argument(
        "--backup-suffix",
        default=".pre_integrated_merge",
        help="Backup dst to dst + suffix if that backup does not exist yet",
    )
    args = parser.parse_args()

    src_path = args.src.expanduser().resolve()
    dst_path = args.dst.expanduser().resolve()
    if not src_path.is_file():
        raise FileNotFoundError(src_path)
    if not dst_path.is_file():
        raise FileNotFoundError(dst_path)

    backup_path = Path(str(dst_path) + args.backup_suffix)
    if not backup_path.is_file():
        shutil.copy2(dst_path, backup_path)
        print(f"backup -> {backup_path}")
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        extra = dst_path.with_name(dst_path.name + f".backup_{ts}")
        shutil.copy2(dst_path, extra)
        print(f"backup already at {backup_path}; extra copy -> {extra}")

    with open(src_path, "rb") as f:
        src = pickle.load(f)
    with open(dst_path, "rb") as f:
        dst = pickle.load(f)

    old = dst["integrated"]
    merge_integrated_row(dst, src)
    new = dst["integrated"]

    with open(dst_path, "wb") as f:
        pickle.dump(dst, f, protocol=pickle.HIGHEST_PROTOCOL)

    def unc(row, k):
        return 100.0 * np.sqrt(max(float(np.asarray(row[k]).flat[0]), 0.0))

    print(f"merged {src_path.name} -> {dst_path}")
    print(f"  keys: {len(old)} -> {len(new)}")
    for k in ("genie", "genie_rate", "genie_ar23", "genie_ar23_rate"):
        if k in old and k in new:
            print(f"  {k}: {unc(old,k):.4f}% -> {unc(new,k):.4f}%")


if __name__ == "__main__":
    main()
