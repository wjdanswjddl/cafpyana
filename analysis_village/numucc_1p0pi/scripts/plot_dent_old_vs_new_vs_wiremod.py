#!/usr/bin/env python3
"""Overlay DENT-old vs DENT-highstats vs WireMod unisim fractional uncertainties.

Keeps ``DENT-updated`` as the old baseline and adds the new highstats cache.
"""
from __future__ import annotations

import os
import pickle
import sys
from collections import defaultdict
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_REPO = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from analysis_village.numucc_1p0pi.scripts import dent_compare as dc
from analysis_village.numucc_1p0pi.syst_histcounts import (
    combine_indep_knob_frac_covs,
    unisim_cov_from_cv_and_var,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig


ROOT = "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final"
DENT_OLD_PKL = f"{ROOT}/DENT-updated/cache/dent_sel_all_hists.pkl"
DENT_NEW_PKL = f"{ROOT}/DENT-highstats/cache/dent_sel_all_hists.pkl"
WIREMOD_NPZ = f"{ROOT}/WireMod/Detector/detector_syst_dict.npz"
OUT_ROOT = f"{ROOT}/DENT-highstats-vs-old-vs-WireMod"


def _load_dent(pkl: str):
    sys.modules["__main__"].StageMetrics = dc.StageMetrics
    sys.modules["__main__"].SampleSummary = dc.SampleSummary
    with open(pkl, "rb") as fh:
        return pickle.load(fh)


def _dent_frac_by_var(payload: dict) -> dict[str, np.ndarray]:
    hists = payload["hists"]
    out = {}
    for var in hists["cv"]:
        if var not in hists["dent"]:
            continue
        pack = unisim_cov_from_cv_and_var(
            np.asarray(hists["cv"][var], dtype=float),
            np.asarray(hists["dent"][var], dtype=float),
        )
        out[var] = np.asarray(pack["frac_unc"], dtype=float) * 100.0
    return out


def _wiremod_frac_by_var(npz_path: str) -> tuple[dict, dict, dict]:
    data = np.load(npz_path, allow_pickle=True)
    yz = data["detector-wiremod_yz"].item()
    xtxw = data["detector-wiremod_xtxw"].item()
    combined = data["detector"].item()

    def _frac(entry: dict) -> dict[str, np.ndarray]:
        out = {}
        for var, pack in entry.items():
            if not isinstance(pack, dict) or "cov_frac" not in pack:
                continue
            cf = np.asarray(pack["cov_frac"], dtype=float)
            out[var] = np.sqrt(np.clip(np.diag(cf), 0, None)) * 100.0
        return out

    return _frac(yz), _frac(xtxw), _frac(combined)


def main() -> int:
    if not os.path.isfile(DENT_NEW_PKL):
        raise SystemExit(f"missing new DENT cache: {DENT_NEW_PKL}")
    if not os.path.isfile(DENT_OLD_PKL):
        raise SystemExit(f"missing old DENT cache: {DENT_OLD_PKL}")

    old_p = _load_dent(DENT_OLD_PKL)
    new_p = _load_dent(DENT_NEW_PKL)
    old_f = _dent_frac_by_var(old_p)
    new_f = _dent_frac_by_var(new_p)
    yz_f, xt_f, wm_f = _wiremod_frac_by_var(WIREMOD_NPZ)

    os.makedirs(f"{OUT_ROOT}/plots", exist_ok=True)

    common = sorted(set(old_f) & set(new_f) & set(wm_f))
    if not common:
        common = sorted(set(old_f) & set(new_f))
    print(f"n_common_vars={len(common)}", flush=True)

    integ_old = float(old_f.get("integrated", [np.nan])[0]) if "integrated" in old_f else np.nan
    integ_new = float(new_f.get("integrated", [np.nan])[0]) if "integrated" in new_f else np.nan
    print(f"integrated DENT-old={integ_old:.4f}%  DENT-highstats={integ_new:.4f}%", flush=True)

    for var in common:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        series = [
            ("WireModYZ", yz_f.get(var), "C0"),
            ("WireModXTXW", xt_f.get(var), "C1"),
            ("WireMod (quad)", wm_f.get(var), "C2"),
            ("DENT (old)", old_f.get(var), "C3"),
            ("DENT (highstats)", new_f.get(var), "C4"),
        ]
        for label, y, color in series:
            if y is None:
                continue
            x = np.arange(len(y) + 1)
            yy = np.concatenate([y, y[-1:]])
            ax.step(x, yy, where="post", label=label, color=color, lw=1.6)
        ax.set_xlabel(f"bin index ({var})")
        ax.set_ylabel("frac. unc. [%]")
        ax.set_title(f"DENT highstats vs old vs WireMod — {var}")
        ax.legend(fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(f"{OUT_ROOT}/plots/unisim_frac_unc__{var}.{ext}")
        plt.close(fig)

    with open(f"{OUT_ROOT}/PROVENANCE.txt", "w") as fh:
        fh.write(
            "\n".join(
                [
                    f"created_utc: {datetime.now(timezone.utc).isoformat()}",
                    f"DENT_old: {DENT_OLD_PKL}",
                    f"DENT_new: {DENT_NEW_PKL}",
                    f"WireMod: {WIREMOD_NPZ}",
                    f"DENT_old pots: {old_p.get('pots')}",
                    f"DENT_new pots: {new_p.get('pots')}",
                    f"integrated DENT-old: {integ_old:.6f}%",
                    f"integrated DENT-highstats: {integ_new:.6f}%",
                    f"n_common_vars: {len(common)}",
                    "",
                ]
            )
        )
    print(f"wrote plots → {OUT_ROOT}/plots", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
