#!/usr/bin/env python3
"""Drop the non-canonical xsec block from a ``sel_all`` GENIE cov_mat_dict.

A ``sel_all`` GENIE run fills the final-stage (Product B) variables as well as
the cut-stage (Product A) ones, so it emits xsec covariances for Product B. Those
disagree with the canonical ``sel_mup`` / ``MC_DF_STAGE=final`` run (up to ~15
percentage points on ``vertex_z``), and only the ``sel_mup`` ones are used.
Removing them here keeps the file from being picked up by mistake.

Only the Product B rows are touched: their xsec entries (every key without the
``_rate`` suffix) go, their rate entries stay. Cut-stage rows are left byte-for-byte
as they were, since this file is the authoritative Product A source.

Writes a new file; the input is never modified.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

sys.path.insert(0, "/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")

from analysis_village.numucc_1p0pi.scripts._genie_pkl_compat import (  # noqa: E402
    load_cov_mat_dict_compat,
)
from analysis_village.numucc_1p0pi.syst_pipeline_walker import (  # noqa: E402
    final_stage_var_configs,
)


def prune(src: Path, dst: Path) -> None:
    cov = load_cov_mat_dict_compat(src)
    prod_b = {vc.var_save_name for vc in final_stage_var_configs() if vc.var_save_name}

    n_dropped = 0
    touched = []
    for var, row in cov.items():
        if var not in prod_b:
            continue
        drop = [k for k in row if not k.endswith("_rate")]
        if not drop:
            continue
        for k in drop:
            del row[k]
        n_dropped += len(drop)
        touched.append(var)

    with open(dst, "wb") as f:
        pickle.dump(cov, f, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"pruned {n_dropped} xsec keys across {len(touched)} Product B variables")
    print(f"  {src}  ({src.stat().st_size / 1e6:.1f} MB)")
    print(f"  -> {dst}  ({dst.stat().st_size / 1e6:.1f} MB)")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("src", type=Path)
    p.add_argument("dst", type=Path)
    p.add_argument(
        "--manifest",
        type=Path,
        help="manifest to copy alongside dst, tagged rate_only",
    )
    a = p.parse_args()

    prune(a.src, a.dst)

    if a.manifest and a.manifest.is_file():
        man = json.loads(a.manifest.read_text())
        man["output_pkl"] = str(a.dst)
        man["productB_xsec_pruned"] = True
        man["productB_xsec_note"] = (
            "sel_all-derived Product B xsec removed; canonical source is the "
            "MC_DF_STAGE=final sel_mup run (syst_disk_sel_mup_20260912_Ar23p)"
        )
        out = a.dst.with_name(a.dst.stem + "_manifest.json")
        out.write_text(json.dumps(man, indent=2) + "\n")
        print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
