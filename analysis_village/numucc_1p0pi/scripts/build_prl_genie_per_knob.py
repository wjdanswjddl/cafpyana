#!/usr/bin/env python3
"""Assemble Product-B per-knob GENIE covariances into one PRL file + mode map.

Sources (waves), in merge priority (earlier wins on key collisions)::

  1. Aug24 mode disks: CCQE, MEC, RES, nonRES, DIS, Other
  2. Sep12 Ar23p disk (full Product-B vars)
  3. VecFF product (forgotten SBN_v1 VecFFCCQEshape)
  4. May combined archive — **EDepFSI FSI π/N only** (retired QE/MEC twins skipped)

Writes under ``PRL/systematics/productB_sel_mup/GENIE/``::

  cov_mat_dict_per_knob.pkl   — flat cov_mat_dict (all knobs; genie*/totals rebuilt)
  knob_mode_map.json          — per-knob wave + production/distributed mode
  per_knob_manifest.json      — provenance / skips / counts

Distributed mode (Ar23p *distributed*): same rules as ``syst_genie_inspect`` —
Ar23p → CCQE/MEC/FSI; Other-mode → Other (COH/NCEL) or FSI; VecFF/ZExp → CCQE.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from analysis_village.numucc_1p0pi.syst_genie_inspect import (  # noqa: E402
    TOTAL_KEYS,
    assign_ar23p_knob_to_mode,
    assign_knob_to_mode,
    assign_other_mode_knob_to_bucket,
    is_retired_edepfsi_twin,
    load_cov_mat_dict,
)
from makedf.geniesyst import GENIE_KNOB_GROUPS  # noqa: E402

SYST_BASE = Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi")
ARCHIVE_GENIE = Path(
    "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/"
    "systematics-final-archive/GENIE/cov_mat_dict.pkl"
)
DEFAULT_OUT = Path(
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/systematics/"
    "productB_sel_mup/GENIE"
)

# Merge order = priority (first source keeps the key on collision).
WAVE_SOURCES: List[Tuple[str, Path]] = [
    ("Aug24_CCQE", SYST_BASE / "syst_disk_sel_mup_Aug24_CCQE/GENIE/cov_mat_dict.pkl"),
    ("Aug24_MEC", SYST_BASE / "syst_disk_sel_mup_Aug24_MEC/GENIE/cov_mat_dict.pkl"),
    ("Aug24_RES", SYST_BASE / "syst_disk_sel_mup_Aug24_RES/GENIE/cov_mat_dict.pkl"),
    ("Aug24_nonRES", SYST_BASE / "syst_disk_sel_mup_Aug24_nonRES/GENIE/cov_mat_dict.pkl"),
    ("Aug24_DIS", SYST_BASE / "syst_disk_sel_mup_Aug24_DIS/GENIE/cov_mat_dict.pkl"),
    ("Aug24_Other", SYST_BASE / "syst_disk_sel_mup_Aug24_Other/GENIE/cov_mat_dict.pkl"),
    ("Sep12_Ar23p", SYST_BASE / "syst_disk_sel_mup_20260912_Ar23p/GENIE/cov_mat_dict.pkl"),
    ("VecFF", SYST_BASE / "syst_disk_sel_mup_VecFF/GENIE/cov_mat_dict.pkl"),
    ("May_EDepFSI_FSI", ARCHIVE_GENIE),
]

# Wave tag → production-mode label for knobs that lack an exact GENIE_KNOB_GROUPS hit.
WAVE_DEFAULT_MODE = {
    "Aug24_CCQE": "CCQE",
    "Aug24_MEC": "MEC",
    "Aug24_RES": "RES",
    "Aug24_nonRES": "nonRES",
    "Aug24_DIS": "DIS",
    "Aug24_Other": "Other",
    "Sep12_Ar23p": "Ar23p",
    "VecFF": "VecFF",
    "May_EDepFSI_FSI": "EDepFSI_FSI",
}


def _base_knob(key: str) -> str:
    return key[: -len("_rate")] if key.endswith("_rate") else key


def _production_mode(knob: str, wave: str) -> str:
    """Physics bucket from GENIE_KNOB_GROUPS / inspect heuristics / wave default."""
    mode = assign_knob_to_mode(knob)
    # Wave disk is authoritative when inspect heuristics land on a soft bucket.
    if mode in ("Other", "Ar23p", "FSI", "EDepFSI_FSI") and wave in WAVE_DEFAULT_MODE:
        wmode = WAVE_DEFAULT_MODE[wave]
        if wmode == "EDepFSI_FSI":
            return "EDepFSI_FSI"
        if wmode != "May_EDepFSI_FSI":
            # Prefer exact wave label for Aug24_* / Sep12_Ar23p / VecFF.
            if wave.startswith("Aug24_") or wave in ("Sep12_Ar23p", "VecFF"):
                return wmode
    if mode == "ZExp":
        return "CCQE"
    if mode == "VecFF":
        return "CCQE"
    if "EDepFSI" in knob and mode in ("Other", "FSI"):
        return "EDepFSI_FSI"
    return mode


def _distributed_mode(knob: str, production_mode: str, wave: str) -> str:
    """Ar23p-distributed bucket used by inspect plots."""
    if wave == "Sep12_Ar23p" or production_mode == "Ar23p":
        return assign_ar23p_knob_to_mode(knob)
    if wave == "Aug24_Other" or production_mode == "Other":
        return assign_other_mode_knob_to_bucket(knob)
    if production_mode in ("VecFF", "ZExp", "CCQE"):
        return "CCQE"
    if production_mode == "EDepFSI_FSI" or ("EDepFSI" in knob and not is_retired_edepfsi_twin(knob)):
        return "FSI"
    if production_mode == "FSI":
        return "FSI"
    if production_mode in ("MEC", "RES", "nonRES", "DIS"):
        return production_mode
    return production_mode


def _accept_edepfsi_key(key: str) -> bool:
    """Keep only non-retired EDepFSI FSI dials from the May combined file."""
    base = _base_knob(key)
    if "EDepFSI" not in base:
        return False
    if is_retired_edepfsi_twin(base):
        return False
    return True


def _ref_shape_for_slug(
    combined: Mapping[str, Mapping[str, np.ndarray]], slug: str
) -> Optional[Tuple[int, ...]]:
    row = combined.get(slug)
    if not row:
        return None
    for k, v in row.items():
        if k in TOTAL_KEYS:
            continue
        return tuple(np.asarray(v).shape)
    return None


def merge_waves(
    sources: Sequence[Tuple[str, Path]],
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, dict], dict]:
    """Return (cov_mat_dict, knob_map, stats)."""
    combined: Dict[str, Dict[str, np.ndarray]] = {}
    # knob base name → metadata (shared by rate/xsec twin keys)
    knob_meta: Dict[str, dict] = {}
    stats: Dict[str, Any] = {
        "sources": [],
        "n_keys_added": 0,
        "n_keys_collision_kept_earlier": 0,
        "n_keys_shape_skip": 0,
        "n_keys_edep_filter_skip": 0,
        "collisions": [],
        "shape_skips": [],
    }

    for wave, path in sources:
        src_stat: Dict[str, Any] = {
            "wave": wave,
            "path": str(path),
            "exists": path.is_file(),
            "n_vars": 0,
            "n_keys_added": 0,
            "n_keys_collision": 0,
            "n_keys_shape_skip": 0,
            "n_keys_edep_filter_skip": 0,
        }
        if not path.is_file():
            stats["sources"].append(src_stat)
            print(f"[skip missing] {wave}: {path}")
            continue

        cov = load_cov_mat_dict(path)
        src_stat["n_vars"] = len(cov)
        print(f"[load] {wave}: {path}  vars={len(cov)}")

        for slug, row in cov.items():
            if not isinstance(row, dict):
                continue
            dest_row = combined.setdefault(slug, {})
            ref = _ref_shape_for_slug(combined, slug)

            for key, mat in row.items():
                if key in TOTAL_KEYS:
                    continue
                base = _base_knob(key)

                if wave == "May_EDepFSI_FSI" and not _accept_edepfsi_key(key):
                    src_stat["n_keys_edep_filter_skip"] += 1
                    stats["n_keys_edep_filter_skip"] += 1
                    continue

                arr = np.asarray(mat, dtype=np.float64)
                if key in dest_row:
                    src_stat["n_keys_collision"] += 1
                    stats["n_keys_collision_kept_earlier"] += 1
                    if len(stats["collisions"]) < 200:
                        stats["collisions"].append(
                            {
                                "slug": slug,
                                "key": key,
                                "kept_wave": knob_meta.get(base, {}).get("wave"),
                                "skipped_wave": wave,
                            }
                        )
                    continue

                if ref is not None and tuple(arr.shape) != ref:
                    src_stat["n_keys_shape_skip"] += 1
                    stats["n_keys_shape_skip"] += 1
                    if len(stats["shape_skips"]) < 200:
                        stats["shape_skips"].append(
                            {
                                "slug": slug,
                                "key": key,
                                "wave": wave,
                                "shape": list(arr.shape),
                                "ref_shape": list(ref),
                            }
                        )
                    continue

                # First non-total matrix for this slug sets the reference shape.
                if ref is None:
                    ref = tuple(arr.shape)

                dest_row[key] = arr
                src_stat["n_keys_added"] += 1
                stats["n_keys_added"] += 1

                if base not in knob_meta:
                    prod = _production_mode(base, wave)
                    knob_meta[base] = {
                        "knob": base,
                        "wave": wave,
                        "wave_path": str(path),
                        "production_mode": prod,
                        "distributed_mode": _distributed_mode(base, prod, wave),
                        "has_rate": key.endswith("_rate"),
                        "has_xsec": not key.endswith("_rate"),
                    }
                else:
                    if key.endswith("_rate"):
                        knob_meta[base]["has_rate"] = True
                    else:
                        knob_meta[base]["has_xsec"] = True

        stats["sources"].append(src_stat)

    # Rebuild independent-sum totals from kept knobs (do not trust per-wave genie*).
    for slug, row in combined.items():
        tot_x = tot_r = None
        for key, mat in row.items():
            if key in TOTAL_KEYS:
                continue
            arr = np.asarray(mat, dtype=np.float64)
            if key.endswith("_rate"):
                tot_r = arr.copy() if tot_r is None else tot_r + arr
            else:
                tot_x = arr.copy() if tot_x is None else tot_x + arr
        # drop any stale totals then rewrite
        for tk in list(TOTAL_KEYS):
            row.pop(tk, None)
        if tot_x is not None:
            row["genie"] = tot_x
        if tot_r is not None:
            row["genie_rate"] = tot_r

    return combined, knob_meta, stats


def write_outputs(
    out_dir: Path,
    combined: Mapping[str, dict],
    knob_meta: Mapping[str, dict],
    stats: dict,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    pkl_path = out_dir / "cov_mat_dict_per_knob.pkl"
    map_path = out_dir / "knob_mode_map.json"
    man_path = out_dir / "per_knob_manifest.json"

    with open(pkl_path, "wb") as f:
        pickle.dump(dict(combined), f, protocol=pickle.HIGHEST_PROTOCOL)

    # Stable sorted map
    knobs_sorted = [knob_meta[k] for k in sorted(knob_meta)]
    by_dist: Dict[str, List[str]] = {}
    by_prod: Dict[str, List[str]] = {}
    by_wave: Dict[str, List[str]] = {}
    for m in knobs_sorted:
        by_dist.setdefault(m["distributed_mode"], []).append(m["knob"])
        by_prod.setdefault(m["production_mode"], []).append(m["knob"])
        by_wave.setdefault(m["wave"], []).append(m["knob"])

    map_doc = {
        "schema": "numucc_genie_knob_mode_map_v1",
        "created": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "notes": {
            "distributed_mode": (
                "Ar23p knobs → CCQE (QE/ZExp) / MEC / FSI; "
                "Other-mode knobs → Other (COH/NCEL) or FSI; "
                "VecFF/ZExp → CCQE; EDepFSI FSI π/N → FSI. "
                "Matches syst_genie_inspect.mode_totals_ar23p_distributed."
            ),
            "production_mode": (
                "Wave / GENIE_KNOB_GROUPS bucket before Ar23p redistribution "
                "(VecFF/ZExp folded into CCQE)."
            ),
            "totals_in_pkl": (
                "genie / genie_rate are independent sums of kept per-knob cov_frac "
                "(not the slim_v3 product total)."
            ),
        },
        "n_knobs": len(knobs_sorted),
        "knobs": knobs_sorted,
        "by_distributed_mode": {k: sorted(v) for k, v in sorted(by_dist.items())},
        "by_production_mode": {k: sorted(v) for k, v in sorted(by_prod.items())},
        "by_wave": {k: sorted(v) for k, v in sorted(by_wave.items())},
    }
    with open(map_path, "w") as f:
        json.dump(map_doc, f, indent=2)
        f.write("\n")

    # Integrated summary for quick QA
    integ = combined.get("integrated", {})
    integ_summary = []
    for m in knobs_sorted:
        kn = m["knob"]
        row = {"knob": kn, "distributed_mode": m["distributed_mode"], "wave": m["wave"]}
        for kind, key in (("xsec", kn), ("rate", kn + "_rate")):
            if key in integ:
                cf = np.asarray(integ[key], dtype=np.float64)
                row[f"{kind}_frac_unc_pct"] = float(100.0 * np.sqrt(max(cf[0, 0], 0.0)))
        integ_summary.append(row)

    man = {
        "schema": "numucc_genie_per_knob_manifest_v1",
        "created": map_doc["created"],
        "outputs": {
            "cov_mat_dict_per_knob_pkl": str(pkl_path),
            "knob_mode_map_json": str(map_path),
        },
        "n_vars": len(combined),
        "variables": sorted(combined.keys()),
        "n_knobs": len(knobs_sorted),
        "merge_stats": stats,
        "genie_knob_groups_keys": sorted(GENIE_KNOB_GROUPS.keys()),
        "integrated_frac_unc_pct": integ_summary,
    }
    with open(man_path, "w") as f:
        json.dump(man, f, indent=2)
        f.write("\n")

    print(f"wrote {pkl_path}  ({pkl_path.stat().st_size / 1e6:.1f} MB)")
    print(f"wrote {map_path}  n_knobs={len(knobs_sorted)}")
    print(f"wrote {man_path}")
    print("by_distributed_mode counts:")
    for k, v in sorted(by_dist.items()):
        print(f"  {k:8s} {len(v)}")
    print("by_wave counts:")
    for k, v in sorted(by_wave.items()):
        print(f"  {k:16s} {len(v)}")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    combined, knob_meta, stats = merge_waves(WAVE_SOURCES)
    write_outputs(args.out_dir, combined, knob_meta, stats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
