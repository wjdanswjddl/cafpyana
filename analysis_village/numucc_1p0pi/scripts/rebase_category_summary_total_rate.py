#!/usr/bin/env python3
"""Rewrite CategorySummary ``total_rate`` with signal→total CV rebase.

Flux / G4 / MCstat / GENIE-rate category fracs are left **unrebased** (as stored).
Only ``total_rate`` is rebuilt::

    abs = cov_from_fraccov(frac_sig, n_signal)
    frac_tot = fraccov_from_cov(abs, n_total)

using overlay ``counts_report.npz`` ``mc_signal`` / ``mc_total``. Detector,
cosmics SelectedRate, pot, and ntargets are summed without rebase.
``total_xsec`` and ``categories`` are left byte-identical (verified).

Default target: PRL Product B consumer
``productB_sel_mup``.

Usage::

  /exp/sbnd/app/users/munjung/xsec/freeze/cafpyana/envs/venv_py310_cafpyana/bin/python \\
    analysis_village/numucc_1p0pi/scripts/rebase_category_summary_total_rate.py

  # dry-run (no write):
  ... rebase_category_summary_total_rate.py --dry-run
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
_SCRIPTS = Path(__file__).resolve().parent
for p in (_REPO, _SCRIPTS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from analysis_village.numucc_1p0pi.dataset_locations import (  # noqa: E402
    PRL_OVERLAYS_ROOT,
    PRL_PRODUCT_B_DIR,
    prl_syst_disk_root,
)
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
)
from analysis_village.numucc_1p0pi.syst_category_summary import (  # noqa: E402
    TOTAL_RATE,
    TOTAL_XSEC,
    category_summary_manifest_path,
    category_summary_npz_path,
    load_category_syst_summary,
    rebuild_total_rate_in_summary_pack,
)


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _load_overlay_counts(counts_npz: Path) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Return ``{slug: (mc_signal, mc_total)}`` from a counts_report.npz."""
    z = np.load(counts_npz, allow_pickle=True)
    out: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    # Keys look like ``muon-p__mc_signal`` / ``muon-p__mc_total``.
    bases = set()
    for k in z.files:
        if k.endswith("__mc_signal"):
            bases.add(k[: -len("__mc_signal")])
    for base in bases:
        sig_k = f"{base}__mc_signal"
        tot_k = f"{base}__mc_total"
        if tot_k not in z.files:
            continue
        out[base] = (
            np.asarray(z[sig_k], dtype=np.float64),
            np.asarray(z[tot_k], dtype=np.float64),
        )
    return out


class _VC:
    """Minimal VariableConfig stand-in for pack rebuild."""

    def __init__(self, var_save_name: str, bins: np.ndarray, bin_centers: np.ndarray):
        self.var_save_name = var_save_name
        self.bins = bins
        self.bin_centers = bin_centers


def _vc_for_pack(vsn: str, pack: Mapping[str, Any], vc_by: Mapping[str, Any]):
    if vsn in vc_by:
        return vc_by[vsn]
    bins_arr = np.asarray(pack.get("bins", [0.0, 1.0]), dtype=float)
    centers = (
        0.5 * (bins_arr[:-1] + bins_arr[1:])
        if len(bins_arr) > 1
        else np.array([0.5])
    )
    return _VC(vsn, bins_arr, centers)


def main(argv: Optional[list] = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--syst-root",
        default=str(prl_syst_disk_root("B")),
        help=f"Product B syst disk (default: …/{PRL_PRODUCT_B_DIR})",
    )
    p.add_argument(
        "--counts-npz",
        default=str(
            Path(PRL_OVERLAYS_ROOT)
            / PRL_PRODUCT_B_DIR
            / "counts_report.npz"
        ),
        help="Overlay counts_report.npz with mc_signal / mc_total",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Compute + verify only; do not rewrite the NPZ",
    )
    p.add_argument(
        "--backup",
        action="store_true",
        default=True,
        help="Copy NPZ to *.pre_rate_rebase.npz before overwrite (default)",
    )
    p.add_argument(
        "--no-backup",
        action="store_false",
        dest="backup",
        help="Skip backup copy",
    )
    args = p.parse_args(argv)

    syst_root = Path(args.syst_root).resolve()
    counts_npz = Path(args.counts_npz).resolve()
    cat_npz = Path(category_summary_npz_path(str(syst_root)))
    if not cat_npz.is_file():
        print(f"missing CategorySummary: {cat_npz}", file=sys.stderr)
        return 2
    if not counts_npz.is_file():
        print(f"missing counts report: {counts_npz}", file=sys.stderr)
        return 2

    print(f"syst_root   : {syst_root}")
    print(f"category.npz: {cat_npz}")
    print(f"counts.npz  : {counts_npz}")

    counts = _load_overlay_counts(counts_npz)
    print(f"  overlay vars with signal/total: {len(counts)}")

    summary = load_category_syst_summary(str(cat_npz))
    by_var = summary["by_var"]
    vc_by = {vc.var_save_name: vc for vc in CORE_SELECTED_EVT_VARIABLE_CONFIGS}

    # Snapshot total_xsec for byte-identity check after rewrite.
    xsec_before: Dict[str, Dict[str, np.ndarray]] = {}
    for vsn, pack in by_var.items():
        if TOTAL_XSEC in pack:
            xsec_before[vsn] = {
                k: np.asarray(pack[TOTAL_XSEC][k]).copy()
                for k in pack[TOTAL_XSEC]
                if hasattr(pack[TOTAL_XSEC][k], "__array__")
                or isinstance(pack[TOTAL_XSEC][k], np.ndarray)
            }

    updated = 0
    skipped = []
    rate_delta = []
    out_packs: Dict[str, Any] = {}
    for vsn, pack in by_var.items():
        # pack may be shared; deep-ish copy of mutable tree we will edit
        pack = dict(pack)
        pack["categories"] = dict(pack.get("categories") or {})
        if TOTAL_XSEC in pack:
            pack[TOTAL_XSEC] = dict(pack[TOTAL_XSEC])
        if TOTAL_RATE in pack:
            pack[TOTAL_RATE] = dict(pack[TOTAL_RATE])

        if vsn not in counts:
            skipped.append(f"{vsn}: no overlay counts")
            out_packs[vsn] = pack
            continue
        n_sig, n_tot = counts[vsn]
        vc = _vc_for_pack(vsn, pack, vc_by)
        old_rate = None
        if TOTAL_RATE in pack and "cov_frac" in pack[TOTAL_RATE]:
            old_rate = np.asarray(pack[TOTAL_RATE]["cov_frac"], dtype=np.float64).copy()
        try:
            new_blk = rebuild_total_rate_in_summary_pack(
                pack, vc, n_signal=n_sig, n_total=n_tot
            )
        except Exception as ex:
            skipped.append(f"{vsn}: {ex}")
            out_packs[vsn] = pack
            continue
        updated += 1
        new_rate = np.asarray(new_blk["cov_frac"], dtype=np.float64)
        if old_rate is not None and old_rate.shape == new_rate.shape:
            d = float(np.nanmax(np.abs(new_rate - old_rate)))
            rate_delta.append((vsn, d))
        out_packs[vsn] = pack

    # Verify total_xsec untouched in the in-memory packs.
    xsec_ok = True
    for vsn, before in xsec_before.items():
        after = out_packs[vsn].get(TOTAL_XSEC) or {}
        for k, arr_b in before.items():
            arr_a = np.asarray(after[k])
            if arr_a.shape != arr_b.shape or not np.array_equal(arr_a, arr_b):
                print(f"ERROR: total_xsec[{vsn!r}].{k} changed", file=sys.stderr)
                xsec_ok = False
    if xsec_ok:
        print(f"  total_xsec identity: OK ({len(xsec_before)} vars)")
    else:
        print("  total_xsec identity: FAILED", file=sys.stderr)
        return 3

    rate_delta.sort(key=lambda t: -t[1])
    print(f"  total_rate rebuilt: {updated}  skipped: {len(skipped)}")
    if rate_delta:
        print("  max |Δcov_frac| (top 8):")
        for vsn, d in rate_delta[:8]:
            print(f"    {vsn:24s}  {d:.6g}")
    if skipped:
        print("  skipped:")
        for s in skipped[:20]:
            print(f"    {s}")

    if args.dry_run:
        print("dry-run: not writing")
        return 0 if xsec_ok else 3

    if args.backup:
        bak = cat_npz.with_suffix(".pre_rate_rebase.npz")
        if not bak.is_file():
            shutil.copy2(cat_npz, bak)
            print(f"  backup -> {bak}")
        else:
            print(f"  backup exists (kept): {bak}")

    # Write: object-array packing matches export / detector-swap rebuilders
    np.savez_compressed(
        cat_npz, **{k: np.array(v, dtype=object) for k, v in out_packs.items()}
    )

    # Re-load and re-check total_xsec vs pre-write snapshot
    z2 = np.load(cat_npz, allow_pickle=True)
    for vsn, before in xsec_before.items():
        after_pack = z2[vsn].item()
        after = after_pack.get(TOTAL_XSEC) or {}
        for k, arr_b in before.items():
            arr_a = np.asarray(after[k])
            if arr_a.shape != arr_b.shape or not np.array_equal(arr_a, arr_b):
                print(
                    f"ERROR post-write: total_xsec[{vsn!r}].{k} changed",
                    file=sys.stderr,
                )
                return 3
    print("  post-write total_xsec identity: OK")

    man_path = Path(category_summary_manifest_path(str(cat_npz)))
    manifest: Dict[str, Any] = {}
    if man_path.is_file():
        manifest = json.loads(man_path.read_text())
    manifest.update(
        {
            "total_rate_rebased_utc": _utc(),
            "total_rate_rebase": {
                "method": "rebase_fraccov_signal_to_total on flux/g4/mcstat/genie_rate",
                "counts_npz": str(counts_npz),
                "n_vars_updated": updated,
                "skipped": skipped,
                "total_xsec_untouched": True,
                "categories_unrebased": True,
            },
        }
    )
    man_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"  updated manifest {man_path}")
    print(f"wrote {cat_npz}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
