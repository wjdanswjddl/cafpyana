#!/usr/bin/env python3
"""Build stacked JointDetector / JointCosmics / JointPOT / JointNtargets.

Cosmics uses the published offbeam / intime unisim histograms (n_univ=1).
Detector consumer NPZs do not store n_cv; each independent source is a rank-1
unisim, so the +1σ vector is recovered from that source's per-variable
``cov_frac`` (intra-variable signs from the rank-1 structure; overall polarity
per variable chosen so the vector sum is non-negative). Product B DENT is the
smoothed pack (already ``outer(u_s, u_s)`` with ``u_s ≥ 0``).

Writes under each ``--out-root`` (default: Product B ``JointCC`` and the
constraint-validation JointCC). Does **not** touch May GENIE trees or Flux/G4/GENIE
chunk campaigns.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from pyanalib.covariance import corr_from_fraccov, cov_from_fraccov  # noqa: E402

from analysis_village.numucc_1p0pi.dataset_locations import (  # noqa: E402
    default_syst_disk_cc_root,
    prl_overlay_counts_npz,
    prl_syst_disk_root,
)
from analysis_village.numucc_1p0pi.syst_category_summary import (  # noqa: E402
    CAT_COSMICS,
    CAT_DETECTOR,
    CAT_NTARGETS,
    CAT_POT,
    NTARGETS_FRAC_UNC_PCT,
    POT_FRAC_UNC_PCT,
    category_cov_frac,
    load_category_syst_summary,
)
from analysis_village.numucc_1p0pi.syst_cc_joint_multisim_common import (  # noqa: E402
    STACK_SLUG,
    default_constraint_stack_variables,
    joint_stack_meta,
)
from analysis_village.numucc_1p0pi.syst_disk_cc_layout import (  # noqa: E402
    joint_extras_inner_key,
    joint_extras_npz_basename,
    joint_extras_out_dir,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import (  # noqa: E402
    FILE_COSMICS,
    FILE_DETECTOR,
    SUB_COSMICS,
    SUB_DETECTOR,
    category_summary_npz_path,
)

STACK_VARS = default_constraint_stack_variables()
STACK_NAMES = tuple(v.var_save_name for v in STACK_VARS)
# All ``detector_by_wiremod`` keys are independent (YZ/XTXW geometry + calo
# knobs + DENT + smear26). Their cov_frac sum matches the combined ``detector`` pack.


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _unisim_delta_from_rank1(cov_frac: np.ndarray) -> np.ndarray:
    c = np.asarray(cov_frac, dtype=np.float64)
    mag = np.sqrt(np.clip(np.diag(c), 0.0, None))
    if not np.any(mag > 1e-30):
        return mag
    i0 = int(np.argmax(mag))
    raw = c[:, i0]
    signs = np.ones_like(mag)
    pos = mag > 1e-30
    signs[pos] = np.sign(raw[pos])
    signs[signs == 0] = 1.0
    d = signs * mag
    if float(d.sum()) < 0.0:
        d = -d
    return d


def _write_stacked(
    out_root: Path,
    category: str,
    cov_frac: np.ndarray,
    *,
    extra_meta: Mapping[str, Any] | None = None,
    by_source: Mapping[str, np.ndarray] | None = None,
) -> Path:
    out_dir = Path(joint_extras_out_dir(str(out_root), category))
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / joint_extras_npz_basename(category)
    inner = joint_extras_inner_key(category)
    meta = joint_stack_meta(STACK_VARS, bkgd_subtract=False)
    if extra_meta:
        meta = {**meta, **dict(extra_meta)}
    frac = np.asarray(cov_frac, dtype=float)
    # Absolute cov in the NPZ is a placeholder; the loader rescales with analysis μ.
    dummy_cv = np.ones(frac.shape[0], dtype=float)
    cell: Dict[str, Any] = {
        inner: {
            "cov_frac": frac,
            "cov": cov_from_fraccov(frac, dummy_cv),
            "corr": corr_from_fraccov(frac),
        },
        "meta": meta,
    }
    if by_source:
        cell["%s_by_source" % inner] = {
            kn: {
                "cov_frac": np.asarray(cf, dtype=float),
                "cov": cov_from_fraccov(cf, dummy_cv),
                "corr": corr_from_fraccov(cf),
            }
            for kn, cf in by_source.items()
        }
    np.savez_compressed(path, **{STACK_SLUG: np.array(cell, dtype=object)})
    return path


def _detector_sources(by0: Mapping[str, Any]) -> list[str]:
    return list(by0.keys())


def _rel_rms_diag(a: np.ndarray, b: np.ndarray) -> float:
    da = np.diag(np.asarray(a, dtype=float))
    db = np.diag(np.asarray(b, dtype=float))
    denom = np.maximum(np.abs(db), 1e-30)
    return float(np.sqrt(np.mean(((da - db) / denom) ** 2)))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--out-root",
        action="append",
        default=None,
        type=Path,
        help="JointCC root to write (repeatable). Default: PRL Product B JointCC + validation JointCC.",
    )
    p.add_argument("--diag-tol", type=float, default=0.12)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    prl = prl_syst_disk_root("B")
    counts_npz = prl_overlay_counts_npz("B")
    defaults = (
        default_syst_disk_cc_root(),
        Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/constraint_validation_mu_p_dirz/JointCC"),
    )
    out_roots = tuple(args.out_root) if args.out_root else defaults

    det_z = np.load(prl / SUB_DETECTOR / FILE_DETECTOR, allow_pickle=True)
    cos_z = np.load(prl / SUB_COSMICS / FILE_COSMICS, allow_pickle=True)
    cat_sum = load_category_syst_summary(category_summary_npz_path(str(prl)))
    by = det_z["detector_by_wiremod"].item()
    leaves = _detector_sources(by[STACK_NAMES[0]])
    missing = [k for k in leaves if any(k not in by[n] for n in STACK_NAMES)]
    if missing:
        raise SystemExit("detector source %s missing on some stack vars" % missing)

    det_by: Dict[str, np.ndarray] = {}
    det_frac = None
    for src in leaves:
        deltas = []
        for name in STACK_NAMES:
            deltas.append(_unisim_delta_from_rank1(np.asarray(by[name][src]["cov_frac"], dtype=float)))
        d = np.concatenate(deltas)
        cf = np.outer(d, d)
        det_by[src] = cf
        det_frac = cf if det_frac is None else det_frac + cf

    # Cosmics SelectedRate is a rank-1 unisim; recover d from the published pack
    # (template was flattened to 100% then × contamination; some bins are sanitized).
    d_cos_parts = []
    for name in STACK_NAMES:
        cell = cos_z[name].item()
        sel_cf = np.asarray(cell["SelectedRate"]["rate"]["cov_frac"], dtype=float)
        d_cos_parts.append(_unisim_delta_from_rank1(sel_cf))
    d_cos = np.concatenate(d_cos_parts)
    cos_frac = np.outer(d_cos, d_cos)

    pot_u = float(POT_FRAC_UNC_PCT) / 100.0
    nt_u = float(NTARGETS_FRAC_UNC_PCT) / 100.0
    n_tot = int(sum(len(v.bin_centers) for v in STACK_VARS))
    pot_frac = np.outer(np.full(n_tot, pot_u), np.full(n_tot, pot_u))
    nt_frac = np.outer(np.full(n_tot, nt_u), np.full(n_tot, nt_u))

    sl = {}
    o = 0
    for v in STACK_VARS:
        n = int(len(v.bin_centers))
        sl[v.var_save_name] = slice(o, o + n)
        o += n

    report: Dict[str, Any] = {
        "created_utc": _utc(),
        "detector_sources": leaves,
        "detector_hist_source": "rank1_recovery_from_per_var_cov_frac (no n_cv on consumer NPZ)",
        "cosmics_hist_source": "rank1_recovery of SelectedRate.rate.cov_frac (offbeam/intime hists on NPZ)",
        "diag_vs_prl": {},
    }
    ok = True
    for name in STACK_NAMES:
        checks = {
            "detector": _rel_rms_diag(det_frac[sl[name], sl[name]], category_cov_frac(cat_sum, name, CAT_DETECTOR)),
            "cosmics": _rel_rms_diag(cos_frac[sl[name], sl[name]], category_cov_frac(cat_sum, name, CAT_COSMICS)),
            "pot": _rel_rms_diag(pot_frac[sl[name], sl[name]], category_cov_frac(cat_sum, name, CAT_POT)),
            "ntargets": _rel_rms_diag(nt_frac[sl[name], sl[name]], category_cov_frac(cat_sum, name, CAT_NTARGETS)),
        }
        report["diag_vs_prl"][name] = checks
        for cat, rel in checks.items():
            print("[check] %s %s rel_rms_diag=%.4e %s" % (cat, name, rel, "PASS" if rel <= args.diag_tol else "FAIL"), flush=True)
            ok = ok and rel <= args.diag_tol
    if not ok:
        print("[stop] diagonal mismatch vs PRL CategorySummary; not writing", flush=True)
        print(json.dumps(report, indent=2))
        return 2

    meta_det = {
        "recipe": "sum of independent detector unisims; +1σ vector recovered from rank-1 cov_frac",
        "sources": leaves,
        "detector_npz": str(prl / SUB_DETECTOR / FILE_DETECTOR),
    }
    meta_cos = {
        "recipe": "SelectedRate: published template is 100% unisim; joint d = stacked contamination fraction",
        "cosmics_npz": str(prl / SUB_COSMICS / FILE_COSMICS),
    }
    written = []
    for root in out_roots:
        root = Path(root).expanduser().resolve()
        root.mkdir(parents=True, exist_ok=True)
        written.append(str(_write_stacked(root, "detector", det_frac, extra_meta=meta_det, by_source=det_by)))
        written.append(str(_write_stacked(root, "cosmics", cos_frac, extra_meta=meta_cos)))
        written.append(str(_write_stacked(root, "pot", pot_frac, extra_meta={"frac_unc": pot_u})))
        written.append(str(_write_stacked(root, "ntargets", nt_frac, extra_meta={"frac_unc": nt_u})))
        (root / "unisim_extras_manifest.json").write_text(json.dumps({**report, "out_root": str(root)}, indent=2) + "\n")
        print("[ok] wrote extras under %s" % root, flush=True)
    report["written"] = written
    print(json.dumps({"written": written, "diag_vs_prl": report["diag_vs_prl"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
