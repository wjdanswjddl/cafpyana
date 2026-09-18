#!/usr/bin/env python3
"""Make WireMod dent-style envelope plots from products cache / NPZs.

Reads ``wiremod_sel_all_products.pkl`` (and optionally refreshes NPZs) then
writes Product A/B rate+ratio+unc plots under ``WireMod/plots/``.
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.dataset_locations import PLOTS_BASE
from analysis_village.numucc_1p0pi.scripts import dent_compare as dc
from analysis_village.numucc_1p0pi.syst_detvar_common import (
    WIREMOD_ENVELOPE_SHIFTED,
    build_wiremod_detector_dict,
    envelope_univ_lo_hi,
    frac_unc_pct_from_pack,
    log,
    save_detector_npz,
    wiremod_component_shifted_univs,
    wiremod_geometry_hists_for_envelope,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import (
    FILE_DETECTOR,
    FILE_DETECTOR_SEL,
    SUB_DETECTOR,
    SUB_DETECTOR_SEL,
)
from analysis_village.numucc_1p0pi.utils import dpi


def _stairs_xy(bins, y):
    """Edge-aligned step polyline for fill_between / plot (matches hist bins)."""
    bins = np.asarray(bins, dtype=float)
    y = np.asarray(y, dtype=float)
    return bins, np.append(y, y[-1])


def _envelope_lo_hi(hists, var_name, shifted_univs=WIREMOD_ENVELOPE_SHIFTED, n_cv=None):
    """Return (n_cv, lo, hi) for the **actual** min/max among WireMod universes.

    No symmetrization about CV — the band is the true envelope of the chi2_*
    variation histograms. Uncertainty (right panel / NPZs) separately uses
    ``max_u |n_u - n_cv|`` per bin.
    """
    if n_cv is None:
        n_cv = np.asarray(hists["cv"][var_name], dtype=float)
    else:
        n_cv = np.asarray(n_cv, dtype=float)
    lo, hi = envelope_univ_lo_hi(hists, var_name, shifted_univs)
    return n_cv, lo, hi


def _plot_wiremod_envelope_compare(
    vsn: str,
    all_hists: dict,
    bins,
    xlabel: str,
    dict_det: dict,
    out_dir: Path,
    *,
    tag: str,
    title: str | None = None,
    cv_hist: np.ndarray | None = None,
):
    bins = np.asarray(bins, dtype=float)
    centers = 0.5 * (bins[:-1] + bins[1:])
    fig = plt.figure(figsize=(11.0, 5.2), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 1.0], hspace=0.05, wspace=0.22)
    ax_rate = fig.add_subplot(gs[0, 0])
    ax_ratio = fig.add_subplot(gs[1, 0], sharex=ax_rate)
    ax_unc = fig.add_subplot(gs[0, 1])

    geom_styles = (
        ("YZ", "wiremod_yz", "C0", 0.28),
        ("XTXW", "wiremod_xtxw", "C1", 0.22),
    )
    n_cv_ref = None if cv_hist is None else np.asarray(cv_hist, dtype=float)
    for lab, tag_key, color, alpha in geom_styles:
        if lab not in all_hists:
            continue
        if not any(vsn in h for h in all_hists[lab].values()):
            continue
        n_cv, lo, hi = _envelope_lo_hi(all_hists[lab], vsn, n_cv=n_cv_ref)
        if n_cv_ref is None:
            n_cv_ref = n_cv
        x_lo, y_lo = _stairs_xy(bins, lo)
        _, y_hi = _stairs_xy(bins, hi)
        ax_rate.fill_between(
            x_lo,
            y_lo,
            y_hi,
            step="post",
            color=color,
            alpha=alpha,
            linewidth=0,
            label=f"{lab} envelope",
            zorder=1,
        )
        # Same envelope as rate panel, expressed as Variation / CV.
        ratio_lo = np.full_like(n_cv, np.nan, dtype=float)
        ratio_hi = np.full_like(n_cv, np.nan, dtype=float)
        ok = n_cv > 0
        ratio_lo[ok] = lo[ok] / n_cv[ok]
        ratio_hi[ok] = hi[ok] / n_cv[ok]
        if np.any(ok):
            xr, rlo = _stairs_xy(bins, np.where(ok, ratio_lo, np.nan))
            _, rhi = _stairs_xy(bins, np.where(ok, ratio_hi, np.nan))
            ax_ratio.fill_between(
                xr, rlo, rhi, step="post", color=color, alpha=alpha, linewidth=0, label=lab
            )
        pack = dict_det.get(f"detector-{tag_key}", {}).get(vsn)
        if pack is not None:
            ax_unc.hist(
                centers,
                bins=bins,
                weights=frac_unc_pct_from_pack(pack),
                histtype="step",
                lw=1.8,
                color=color,
                label=lab,
            )

    if n_cv_ref is not None:
        ax_rate.hist(
            centers,
            bins=bins,
            weights=n_cv_ref,
            histtype="step",
            lw=1.8,
            color="black",
            label="CV",
            zorder=3,
        )

    ax_rate.set_ylabel("Events")
    ax_rate.legend(fontsize=8)
    ax_rate.grid(True, alpha=0.3)
    ax_rate.tick_params(labelbottom=False)
    if len(centers) == 1:
        ax_rate.set_xlim(bins[0], bins[-1])

    ax_ratio.axhline(1.0, color="black", ls="--", lw=1.0)
    ax_ratio.set_ylabel("Variation / CV")
    ax_ratio.set_xlabel(xlabel)
    ax_ratio.legend(fontsize=8)
    ax_ratio.grid(True, alpha=0.3)

    ax_unc.set_ylabel("Uncertainty [%]")
    ax_unc.set_xlabel(xlabel)
    ax_unc.set_ylim(bottom=0)
    ax_unc.legend(fontsize=8)
    ax_unc.grid(True, alpha=0.3)
    if len(centers) == 1:
        ax_unc.set_xlim(bins[0], bins[-1])

    fig.suptitle(title if title is not None else vsn, fontsize=11)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"wiremod_{tag}__{vsn}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--out-base",
        default=os.environ.get("WIREMOD_OUT_BASE", str(Path(PLOTS_BASE) / "systematics-final" / "WireMod")),
    )
    p.add_argument("--skip-inspect", action="store_true")
    args = p.parse_args(argv)

    out_base = Path(args.out_base)
    cache_path = out_base / "cache" / "wiremod_sel_all_products.pkl"
    if not cache_path.is_file():
        raise SystemExit(f"missing products cache: {cache_path} (wait for hist walk)")

    with open(cache_path, "rb") as fh:
        payload = pickle.load(fh)
    by_geom = payload["by_geom"]
    if "cv" not in payload:
        raise SystemExit(
            f"{cache_path} lacks external CV products (cv_role envelope_baseline). Rebuild merge."
        )
    cv_prod = payload["cv"]
    cut_names = next(iter(by_geom.values()))["cut_var_names"]
    final_names = next(iter(by_geom.values()))["final_var_names"]
    log(
        f"loaded {cache_path} geoms={list(by_geom)} cv_role={payload.get('cv_role')} "
        f"pot={payload.get('pot_by_variation')}"
    )

    def _all_hists(product: str):
        return {
            lab: wiremod_geometry_hists_for_envelope(prod["by_universe"], product=product)
            for lab, prod in by_geom.items()
        }

    all_hists_a = _all_hists("cut")
    all_hists_b = _all_hists("final")
    dict_a = build_wiremod_detector_dict(
        all_hists_a,
        cut_names,
        wiremod_labels=("YZ", "XTXW"),
        shifted_univs=WIREMOD_ENVELOPE_SHIFTED,
        cv_hists=cv_prod["hists_cut"],
    )
    dict_b = build_wiremod_detector_dict(
        all_hists_b,
        final_names,
        wiremod_labels=("YZ", "XTXW"),
        shifted_univs=WIREMOD_ENVELOPE_SHIFTED,
        cv_hists=cv_prod["hists_final"],
    )

    det_a = out_base / SUB_DETECTOR_SEL
    det_b = out_base / SUB_DETECTOR
    npz_a = det_a / FILE_DETECTOR_SEL
    npz_b = det_b / FILE_DETECTOR
    save_detector_npz(
        dict_a,
        npz_a,
        manifest={
            "source": "WireMod",
            "product": "A_selection",
            "method": "maxabs_dev_unc_actual_envelope_vs_matched_cv",
            "cv_role": "envelope_baseline",
            "shifted_univs": list(WIREMOD_ENVELOPE_SHIFTED),
            "n_vars": len(dict_a.get("detector", {})),
        },
    )
    save_detector_npz(
        dict_b,
        npz_b,
        manifest={
            "source": "WireMod",
            "product": "B_measurement",
            "method": "maxabs_dev_unc_actual_envelope_vs_matched_cv",
            "cv_role": "envelope_baseline",
            "shifted_univs": list(WIREMOD_ENVELOPE_SHIFTED),
            "n_vars": len(dict_b.get("detector", {})),
        },
    )
    log(f"Product A → {npz_a}")
    log(f"Product B → {npz_b}")

    final_defs = dc.build_final_var_defs()
    cut_defs = dc.build_sel_all_var_defs()
    fig_b = out_base / "plots" / "product_B_final"
    fig_a = out_base / "plots" / "product_A_selection"

    n_b = 0
    for vsn, cfg in final_defs.items():
        _plot_wiremod_envelope_compare(
            vsn,
            all_hists_b,
            cfg["bins"],
            cfg.get("label", vsn),
            dict_b,
            fig_b,
            tag="final",
            cv_hist=cv_prod["hists_final"].get(vsn),
        )
        n_b += 1
    log(f"Product B plots: {n_b} → {fig_b}")

    n_a = 0
    for vsn, cfg in cut_defs.items():
        _plot_wiremod_envelope_compare(
            vsn,
            all_hists_a,
            cfg["bins"],
            cfg.get("label", vsn),
            dict_a,
            fig_a,
            tag="cut",
            title=cfg.get("stage_key", "") or vsn,
            cv_hist=cv_prod["hists_cut"].get(vsn),
        )
        n_a += 1
    log(f"Product A plots: {n_a} → {fig_a}")

    if not args.skip_inspect:
        inspect_dir = out_base / "inspect_envelopes"
        inspect_npz = inspect_dir / "npz"
        inspect_npz.mkdir(parents=True, exist_ok=True)
        summary = {"product": "B_measurement", "note": "inspection only", "components": {}}
        for comp, shifted in wiremod_component_shifted_univs().items():
            d = build_wiremod_detector_dict(
                all_hists_b,
                final_names,
                wiremod_labels=("YZ", "XTXW"),
                shifted_univs=shifted,
                cv_hists=cv_prod["hists_final"],
            )
            out_npz = inspect_npz / f"wiremod_envelope_{comp}_productB.npz"
            save_detector_npz(
                d,
                out_npz,
                manifest={
                    "source": "WireMod",
                    "product": "B_measurement_inspect",
                    "component": comp,
                    "method": "component_envelope",
                    "shifted_univs": list(shifted),
                    "not_for_downstream": True,
                },
            )
            row = {"shifted_univs": list(shifted)}
            for key in ("detector-wiremod_yz", "detector-wiremod_xtxw", "detector"):
                pack = d.get(key, {}).get("integrated")
                if pack is None:
                    continue
                w = frac_unc_pct_from_pack(pack)
                row[key] = float(np.asarray(w).ravel()[0]) if len(np.asarray(w).ravel()) else float("nan")
            summary["components"][comp] = row
            log(f"  inspect {comp} → {out_npz.name}")
        summary_path = inspect_dir / "component_envelope_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        log(f"wrote {summary_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
