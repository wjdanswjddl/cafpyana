#!/usr/bin/env python3
"""Set JointGenie combined = ``GENIE_slim_v3`` (full joint, including off-diagonals).

Default (this test): combined is the FSI_compare product-universe covariance
only. The May−Sep SuSA interpolator splice is off; pass ``--mec-splice`` to
put it back::

    C' = C(slim_v3) - C_Sep(interp) + C_May(interp)

Those interpolator joints are stacked from per-variable ``rate_univ_cv``
histograms at the same universe index. SBN_v1 MEC is not swapped on these
variables (no Sep per-knob matrix).

May 14 trees are read-only. Flux / G4 / extras are not touched.
A spliced combined, if present, is copied to ``*.bak_with_may_mec_joint``.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from pyanalib.covariance import corr_from_fraccov, cov_from_fraccov  # noqa: E402

from analysis_village.numucc_1p0pi.dataset_locations import (  # noqa: E402
    default_syst_disk_cc_root,
    prl_constraint_validation_root,
    prl_overlay_counts_npz,
    prl_syst_disk_root,
)
from analysis_village.numucc_1p0pi.scripts.constraint_joint_prl_check import (  # noqa: E402
    CAT_GENIE_RATE,
    FILE_JOINT_GENIE_COMBINED,
    MAY_AR23P_CHUNKS,
    MAY_INTERP_KNOBS,
    MAY_MEC_CHUNKS,
    MAY_SBNV1_MEC_KNOBS,
    STACK_NAMES,
    _assert_not_protected,
    _rel_rms_diag,
    _slices,
    _write_stacked_npz,
)
from analysis_village.numucc_1p0pi.scripts.joint_genie_univ_stack import (  # noqa: E402
    SEP_AR23P_CHUNKS,
    SEP_INTERP_KNOBS,
    load_may_mec_joints,
    load_sep_interp_joints,
    overlay_pack,
    product_b_mec_frac_delta,
    psd_clip,
)
from analysis_village.numucc_1p0pi.syst_category_summary import (  # noqa: E402
    category_cov_frac,
    load_category_syst_summary,
    rebase_fraccov_signal_to_total,
)
from analysis_village.numucc_1p0pi.scripts.conditional_constraint_validation import (  # noqa: E402
    load_overlay_constraint_hists,
)
from analysis_village.numucc_1p0pi.syst_disk_cc_layout import joint_genie_out_dir  # noqa: E402
from analysis_village.numucc_1p0pi.syst_disk_layout import category_summary_npz_path  # noqa: E402
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402

BACKUP_SUFFIX = ".bak_before_may_mec_joint"
SPLICED_BACKUP = ".bak_with_may_mec_joint"
SLIM_KEY = "GENIE_slim_v3"
STACK_VARS = (
    VariableConfig.muon_momentum(),
    VariableConfig.muon_direction(),
    VariableConfig.proton_momentum(),
    VariableConfig.proton_direction(),
)
# muon-p (12) × proton-p (13) off-block; slim joint is ~0.56 mean |ρ|.
MIN_CROSS_RHO = 0.25


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def default_cc_roots() -> List[Path]:
    return [
        Path(default_syst_disk_cc_root()),
        prl_constraint_validation_root("mu_p_dirz") / "JointCC",
    ]


def _npz_candidates(dest: Path) -> List[Path]:
    extra = [
        dest.with_name(dest.name + BACKUP_SUFFIX),
        dest.with_name(dest.name + SPLICED_BACKUP),
        dest.with_name(dest.name + ".bak_before_slim_v3_joint"),
        dest.with_name(dest.name + ".bak_before_productB_align"),
    ]
    out = []
    for p in [dest] + extra:
        if p.is_file() and p not in out:
            out.append(p)
    return out


def _cell(path: Path) -> Dict[str, Any]:
    return np.load(path, allow_pickle=True)["stacked_mu_p"].item()


def _pack_from_cell(cell: Mapping, *, prefer_slim: bool = True) -> Optional[Dict[str, np.ndarray]]:
    bk = cell.get("JointGenie_by_knob") or {}
    if not isinstance(bk, dict):
        bk = dict(bk)
    if prefer_slim and SLIM_KEY in bk:
        slim = bk[SLIM_KEY]
        return {
            "cov_frac": np.asarray(slim["cov_frac"], dtype=float),
            "cov": np.asarray(slim["cov"], dtype=float),
            "corr": np.asarray(slim["corr"], dtype=float)
            if slim.get("corr") is not None
            else corr_from_fraccov(slim["cov_frac"]),
        }
    return None


def _load_slim(cc_roots: Sequence[Path]) -> Dict[str, np.ndarray]:
    seen: List[Path] = []
    for root in cc_roots:
        dest = Path(joint_genie_out_dir(str(root))) / FILE_JOINT_GENIE_COMBINED
        for path in _npz_candidates(dest):
            if path in seen:
                continue
            seen.append(path)
            pack = _pack_from_cell(_cell(path), prefer_slim=True)
            if pack is not None:
                print("[slim] GENIE_slim_v3 from %s" % path, flush=True)
                return pack
    raise SystemExit(
        "GENIE_slim_v3 missing from JointGenie_by_knob in %s (and backups). "
        "Re-run the FSI_compare joint GENIE map."
        % [str(p) for p in seen]
    )


def _mean_abs_rho_block(frac: np.ndarray, sl_a: slice, sl_b: slice) -> float:
    corr = corr_from_fraccov(frac)
    return float(np.mean(np.abs(corr[sl_a, sl_b])))


def _knob_overlay_map(
    by_knob: Mapping[str, Mapping[str, Any]], knobs: Sequence[str], mu_nu: np.ndarray
) -> Dict[str, Dict[str, np.ndarray]]:
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for kn in knobs:
        if kn not in by_knob:
            continue
        out[kn] = overlay_pack(np.asarray(by_knob[kn]["cov_frac"], dtype=float), mu_nu)
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cc-root",
        action="append",
        default=None,
        type=Path,
        help="JointCC root to rewrite (repeatable). Default: Product B JointCC + validation JointCC.",
    )
    p.add_argument(
        "--sep-workers",
        type=int,
        default=8,
        help="Parallel pickle readers for Sep 12 Ar23p chunks (default 8).",
    )
    p.add_argument(
        "--mec-splice",
        action="store_true",
        help="Add May−Sep SuSA interpolator joints onto slim_v3 (off by default).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cc_roots = tuple(Path(p).expanduser().resolve() for p in (args.cc_root or default_cc_roots()))
    for root in cc_roots:
        _assert_not_protected(root)

    slim = _load_slim(cc_roots)
    slim_frac = np.asarray(slim["cov_frac"], dtype=float)
    sl = _slices()
    rho_slim = _mean_abs_rho_block(slim_frac, sl["muon-p"], sl["proton-p"])
    print("[check] slim_v3 muon-p × proton-p mean |ρ| = %.3f" % rho_slim, flush=True)
    if rho_slim < MIN_CROSS_RHO:
        print("[stop] slim_v3 cross-variable |ρ| is too small; refusing to write", flush=True)
        return 2

    counts_npz = prl_overlay_counts_npz("B")
    mu_nu = {}
    mu_sig = {}
    for v in STACK_VARS:
        h = load_overlay_constraint_hists(counts_npz, v)
        mu_nu[v.var_save_name] = np.asarray(h["mu_nu"], dtype=float)
        mu_sig[v.var_save_name] = np.asarray(h["mc_signal"], dtype=float)
    mu_nu_stack = np.concatenate([mu_nu[n] for n in STACK_NAMES])

    by_knob: Dict[str, Dict[str, np.ndarray]] = {
        SLIM_KEY: overlay_pack(slim_frac, mu_nu_stack),
    }
    artifact_meta: Dict[str, Any] = {
        "may_ar23p_chunks": MAY_AR23P_CHUNKS,
        "may_mec_chunks": MAY_MEC_CHUNKS,
        "sep_ar23p_chunks": SEP_AR23P_CHUNKS,
        "may_interp_knobs": list(MAY_INTERP_KNOBS),
        "sep_interp_knobs": list(SEP_INTERP_KNOBS),
        "may_sbnv1_knobs": list(MAY_SBNV1_MEC_KNOBS),
        "sbnv1_swapped": False,
        "mec_splice": bool(args.mec_splice),
    }

    if not args.mec_splice:
        frac = slim_frac
        min_eig = 0.0
        clipped = False
        rho_mp = rho_slim
        recipe = (
            "JointGenie combined = GENIE_slim_v3 stacked joint "
            "(FSI_compare product universes, full off-diagonals). "
            "May MEC interpolators not spliced; re-run with --mec-splice to add them."
        )
        print("[info] combined = slim_v3 (MEC splice off; --mec-splice to add back)", flush=True)
    else:
        print("[may] stacking May 14 interpolator + SBN_v1 MEC universes", flush=True)
        may_by_knob = load_may_mec_joints(workers=1)
        print(
            "[may] knobs %d  n_univ interp=%s  sbnv1=%s"
            % (
                len(may_by_knob),
                [int(may_by_knob[k]["n_univ"]) for k in MAY_INTERP_KNOBS if k in may_by_knob],
                [int(may_by_knob[k]["n_univ"]) for k in MAY_SBNV1_MEC_KNOBS if k in may_by_knob],
            ),
            flush=True,
        )
        print("[sep] stacking Sep 12 Ar23p interpolator universes (workers=%d)" % args.sep_workers, flush=True)
        sep_by_knob = load_sep_interp_joints(workers=max(int(args.sep_workers), 1))
        print(
            "[sep] knobs %d  n_univ=%s  n_files=%s"
            % (
                len(sep_by_knob),
                [int(sep_by_knob[k]["n_univ"]) for k in SEP_INTERP_KNOBS if k in sep_by_knob],
                sep_by_knob[SEP_INTERP_KNOBS[0]]["n_files"] if SEP_INTERP_KNOBS[0] in sep_by_knob else None,
            ),
            flush=True,
        )
        delta, may_frac, sep_frac = product_b_mec_frac_delta(may_by_knob, sep_by_knob)
        raw = slim_frac - sep_frac + may_frac
        frac, min_eig = psd_clip(raw)
        clipped = bool(min_eig < 0)
        rho_mp = _mean_abs_rho_block(frac, sl["muon-p"], sl["proton-p"])
        rho_may = _mean_abs_rho_block(may_frac, sl["muon-p"], sl["proton-p"])
        rho_sep = _mean_abs_rho_block(sep_frac, sl["muon-p"], sl["proton-p"])
        print(
            "[check] muon-p × proton-p mean |ρ|  slim=%.3f  May_interp=%.3f  Sep_interp=%.3f  combined=%.3f"
            % (rho_slim, rho_may, rho_sep, rho_mp),
            flush=True,
        )
        print("[psd] min_eig_before_clip=%.3e  clipped=%s" % (min_eig, clipped), flush=True)
        by_knob.update(_knob_overlay_map(may_by_knob, MAY_INTERP_KNOBS + MAY_SBNV1_MEC_KNOBS, mu_nu_stack))
        by_knob.update(_knob_overlay_map(sep_by_knob, SEP_INTERP_KNOBS, mu_nu_stack))
        by_knob["MAY_INTERP_SUM"] = overlay_pack(may_frac, mu_nu_stack)
        by_knob["SEP_INTERP_SUM"] = overlay_pack(sep_frac, mu_nu_stack)
        artifact_meta.update(
            {
                "may_n_files_ar23p": int(may_by_knob[MAY_INTERP_KNOBS[0]]["n_files"]),
                "may_n_files_mec": int(may_by_knob[MAY_SBNV1_MEC_KNOBS[0]]["n_files"]),
                "sep_n_files": int(sep_by_knob[SEP_INTERP_KNOBS[0]]["n_files"]),
                "may_interp_muon_p_proton_p_mean_abs_rho": rho_may,
                "sep_interp_muon_p_proton_p_mean_abs_rho": rho_sep,
                "min_eig_before_clip": min_eig,
                "clipped": clipped,
            }
        )
        recipe = (
            "JointGenie combined = GENIE_slim_v3 + (May − Sep) SuSA interpolator joints "
            "(full 49-bin matrices from same-index universe hists). "
            "SBN_v1 MEC not swapped on constraint vars. Ar23p/VecFF extras not summed."
        )
        val_root = prl_constraint_validation_root("mu_p_dirz")
        _assert_not_protected(val_root)
        val_root.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            val_root / "may_mec_stacked_joint.npz",
            **{("%s__cov_frac" % k): np.asarray(p["cov_frac"]) for k, p in may_by_knob.items()},
            **{("%s__cv" % k): np.asarray(p["cv"]) for k, p in may_by_knob.items()},
            mu_nu=mu_nu_stack,
            var_save_names=np.array(STACK_NAMES, dtype=object),
        )
        np.savez_compressed(
            val_root / "sep_interp_stacked_joint.npz",
            **{("%s__cov_frac" % k): np.asarray(p["cov_frac"]) for k, p in sep_by_knob.items()},
            **{("%s__cv" % k): np.asarray(p["cv"]) for k, p in sep_by_knob.items()},
            mu_nu=mu_nu_stack,
            var_save_names=np.array(STACK_NAMES, dtype=object),
        )
        print("[ok] wrote stacked joints under %s" % val_root, flush=True)

    cov = cov_from_fraccov(frac, mu_nu_stack)

    prl_root = prl_syst_disk_root("B")
    cat_sum = load_category_syst_summary(category_summary_npz_path(str(prl_root)))
    print("[info] combined diag vs PRL genie_rate (rebased signal→μ_nu; informational)", flush=True)
    for name in STACK_NAMES:
        prl = rebase_fraccov_signal_to_total(
            category_cov_frac(cat_sum, name, CAT_GENIE_RATE),
            mu_sig[name],
            mu_nu[name],
        )
        rel = _rel_rms_diag(frac[sl[name], sl[name]], prl)
        print("[diag] %s  rel_rms=%.4e" % (name, rel), flush=True)

    extra_meta = {
        "recipe": recipe,
        "aligned_utc": _utc(),
        "slim_muon_p_proton_p_mean_abs_rho": rho_slim,
        "combined_muon_p_proton_p_mean_abs_rho": rho_mp,
        **artifact_meta,
    }

    written = []
    for root in cc_roots:
        out_dir = Path(joint_genie_out_dir(str(root)))
        _assert_not_protected(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        dest = out_dir / FILE_JOINT_GENIE_COMBINED
        if dest.is_file():
            dest_cell = _cell(dest)
            dest_bk = dest_cell.get("JointGenie_by_knob") or {}
            if "MAY_INTERP_SUM" in dest_bk:
                spliced_bak = dest.with_name(dest.name + SPLICED_BACKUP)
                if not spliced_bak.is_file():
                    shutil.copy2(dest, spliced_bak)
                    print("[bak] kept spliced combined %s" % spliced_bak, flush=True)
            bak = dest.with_name(dest.name + BACKUP_SUFFIX)
            if not bak.is_file():
                shutil.copy2(dest, bak)
                print("[bak] %s" % bak, flush=True)
        _write_stacked_npz(
            dest,
            "JointGenie",
            frac,
            cov,
            extra_meta=extra_meta,
            by_knob=by_knob,
        )
        (out_dir / "joint_genie_align_manifest.json").write_text(
            json.dumps(
                {
                    "created_utc": _utc(),
                    "cc_root": str(root),
                    "recipe": recipe,
                    "slim_muon_p_proton_p_mean_abs_rho": rho_slim,
                    "combined_muon_p_proton_p_mean_abs_rho": rho_mp,
                    "overlay_counts_npz": str(counts_npz),
                    "prl_syst_root": str(prl_root),
                    **artifact_meta,
                },
                indent=2,
            )
            + "\n"
        )
        written.append(str(dest))
        print("[ok] wrote %s" % dest, flush=True)

    print("[ok] JointGenie combined → %s" % written, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
