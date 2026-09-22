"""Helpers for ``systematics-summary.ipynb``: load source disks + total breakdown plots.

Mirrors the Section-2 style of ``systematics-genie-inspect`` (per-curve uncertainty
figures + frac. cov / corr heatmaps), but categories are **systematic sources**
(Flux, G4, …) rather than GENIE interaction modes. Rate vs xsec differ only in
which GENIE pack is used (``genie_rate`` vs ``genie_xsec``).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np

from analysis_village.numucc_1p0pi.syst_category_summary import (
    NTARGETS_FRAC_UNC_PCT,
    POT_FRAC_UNC_PCT,
    cosmics_selected_rate_cov_frac,
    detector_total_cov_frac,
    genie_category_cov_frac,
    genie_knob_covs,
    genie_var_dict,
    mcstat_cov_frac,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import (
    load_genie_disk_payload,
    resolve_genie_disk_path,
    syst_disk_paths,
)
from analysis_village.numucc_1p0pi.syst_genie_inspect import (
    show_cov_corr_heatmaps,
    sum_cov_fracs,
)

# Plot order / colors for source-breakdown figures (keys match build_source_totals).
SOURCE_ORDER_RATE: Tuple[str, ...] = (
    "Flux",
    "G4",
    "MCstat",
    "Detector",
    "Cosmics",
    "GENIE",
    "Exposure",
    "Targets",
)
SOURCE_ORDER_XSEC: Tuple[str, ...] = SOURCE_ORDER_RATE

SOURCE_COLORS: Dict[str, str] = {
    "Flux": "#1f77b4",
    "G4": "#ff7f0e",
    "MCstat": "#2ca02c",
    "Detector": "#9467bd",
    "Cosmics": "#8c564b",
    "GENIE": "#d62728",
    "Exposure": "#7f7f7f",
    "Targets": "#bcbd22",
    "Total": "k",
}

SOURCE_DISPLAY: Dict[str, str] = {
    "Flux": "Flux",
    "G4": "G4",
    "MCstat": "MC stat.",
    "Detector": "Detector",
    "Cosmics": "Cosmics",
    "GENIE": "GENIE",
    "Exposure": "Exposure",
    "Targets": "Targets",
}


def _flat_cov_frac(nbins: int, frac_unc_pct_val: float) -> np.ndarray:
    u = float(frac_unc_pct_val) / 100.0
    v = u * u
    return np.full((int(nbins), int(nbins)), v, dtype=np.float64)


def _try_load_npz(path: Path | str):
    p = Path(path)
    if not p.is_file():
        return None
    return np.load(p, allow_pickle=True)


def _file_produced_str(path: Path | str | None) -> str:
    """Human-readable production time from filesystem mtime (local time)."""
    if path is None:
        return "unknown date"
    p = Path(path)
    if not p.is_file():
        return "missing"
    try:
        return datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
    except OSError:
        return "unknown date"


def _print_loaded(label: str, path: Path | str) -> None:
    print(f"Loaded {label}: {path}  (produced {_file_produced_str(path)})")


def resolve_source_root(
    source: str,
    *,
    unified_root: Path | str | None,
    source_roots: Optional[Mapping[str, Path | str]] = None,
) -> Optional[Path]:
    """Pick a disk root for *source*: explicit map entry, else unified root."""
    if source_roots:
        raw = source_roots.get(source)
        if raw:
            p = Path(raw).expanduser()
            if p.is_dir():
                return p
    if unified_root:
        p = Path(unified_root).expanduser()
        if p.is_dir():
            return p
    return None


def load_source_payloads(
    *,
    unified_root: Path | str | None = None,
    source_roots: Optional[Mapping[str, Path | str]] = None,
) -> Dict[str, Any]:
    """Load Flux/G4/MCstat/Cosmics/Detector/GENIE payloads from one or many roots.

    ``Detector`` is a **single** Product **B** source: the combined NPZ written by
    ``systematics-detector.ipynb`` (``Detector/detector_syst_dict.npz`` =
    WireMod YZ + XTXW + DENT). WireMod / DENT / SCE are not loaded separately
    here — those notebooks feed ``systematics-detector.ipynb`` first.

    ``source_roots`` keys may include ``Flux``, ``G4``, ``MCstat``, ``Cosmics``,
    ``Detector``, ``GENIE``. Missing files are skipped (value ``None``) with a
    printed note.
    """
    out: Dict[str, Any] = {
        "flux_npz": None,
        "g4_npz": None,
        "mcstat_npz": None,
        "cosmics_npz": None,
        "detector_npz": None,
        "genie_blob": None,
        "roots_used": {},
    }

    def _paths_for(label: str) -> Optional[dict]:
        root = resolve_source_root(label, unified_root=unified_root, source_roots=source_roots)
        if root is None:
            return None
        out["roots_used"][label] = str(root)
        return syst_disk_paths(str(root))

    for key, label, path_key in (
        ("flux_npz", "Flux", "flux"),
        ("g4_npz", "G4", "g4"),
        ("mcstat_npz", "MCstat", "mcstat"),
        ("cosmics_npz", "Cosmics", "cosmics"),
        ("detector_npz", "Detector", "detector"),
    ):
        paths = _paths_for(label)
        if paths is None:
            print(f"{label}: no disk root (skipped)")
            continue
        payload = _try_load_npz(paths[path_key])
        if payload is None:
            print(f"{label}: missing {paths[path_key]}")
        else:
            out[key] = payload
            _print_loaded(label, paths[path_key])

    genie_root = resolve_source_root("GENIE", unified_root=unified_root, source_roots=source_roots)
    if genie_root is not None:
        blob = load_genie_disk_payload(str(genie_root))
        gpath = resolve_genie_disk_path(str(genie_root))
        if blob is not None:
            out["genie_blob"] = blob
            out["roots_used"]["GENIE"] = str(genie_root)
            _print_loaded("GENIE", gpath)
        else:
            print(f"GENIE: missing under {genie_root}")

    return out


def _nbins_from_payloads(payloads: Mapping[str, Any], vsn: str) -> Optional[int]:
    for key, inner in (
        ("flux_npz", "flux"),
        ("g4_npz", "G4"),
        ("mcstat_npz", "MCstat"),
    ):
        z = payloads.get(key)
        if z is None or vsn not in z:
            continue
        cell = z[vsn].item()
        if isinstance(cell, dict) and inner in cell and "cov_frac" in cell[inner]:
            return int(np.asarray(cell[inner]["cov_frac"]).shape[0])
    gp = genie_knob_covs(genie_var_dict(payloads.get("genie_blob"), vsn))
    if gp:
        for tot in (gp.get("rate_total"), gp.get("xsec_total")):
            if tot is not None:
                return int(np.asarray(tot).shape[0])
    return None


def build_source_totals(
    payloads: Mapping[str, Any],
    vsn: str,
    kind: str,
    *,
    include_flat: bool = True,
) -> Dict[str, np.ndarray]:
    """``{display_source: cov_frac}`` for one variable; *kind* is ``rate`` or ``xsec``.

    Only the GENIE entry changes between rate and xsec.
    """
    if kind not in ("rate", "xsec"):
        raise ValueError(f"kind must be 'rate' or 'xsec', got {kind!r}")

    out: Dict[str, np.ndarray] = {}
    flux_npz = payloads.get("flux_npz")
    g4_npz = payloads.get("g4_npz")
    if flux_npz is not None and vsn in flux_npz:
        try:
            out["Flux"] = np.asarray(flux_npz[vsn].item()["flux"]["cov_frac"], dtype=np.float64)
        except Exception:
            pass
    if g4_npz is not None and vsn in g4_npz:
        try:
            out["G4"] = np.asarray(g4_npz[vsn].item()["G4"]["cov_frac"], dtype=np.float64)
        except Exception:
            pass

    mc = mcstat_cov_frac(payloads.get("mcstat_npz"), vsn)
    if mc is not None:
        out["MCstat"] = mc

    det = detector_total_cov_frac(
        vsn,
        detector_npz=payloads.get("detector_npz"),
    )
    if det is not None:
        out["Detector"] = det

    cc = cosmics_selected_rate_cov_frac(payloads.get("cosmics_npz"), vsn)
    if cc is not None:
        out["Cosmics"] = cc

    gp = genie_knob_covs(genie_var_dict(payloads.get("genie_blob"), vsn))
    if gp is not None:
        g = genie_category_cov_frac(gp, kind)
        if g is not None:
            out["GENIE"] = g

    if include_flat:
        nb = _nbins_from_payloads(payloads, vsn)
        if nb is None and out:
            nb = int(next(iter(out.values())).shape[0])
        if nb is not None:
            out["Exposure"] = _flat_cov_frac(nb, POT_FRAC_UNC_PCT)
            out["Targets"] = _flat_cov_frac(nb, NTARGETS_FRAC_UNC_PCT)

    return out


def plot_source_frac_unc(
    ax,
    source_totals: Mapping[str, np.ndarray],
    *,
    vc=None,
    kind: str = "rate",
    source_order: Sequence[str] = SOURCE_ORDER_RATE,
):
    """Uncertainty curves per source + Total (genie-inspect Section-2 style)."""
    from analysis_village.numucc_1p0pi.syst_genie_inspect import (
        _centers_bins,
        frac_unc_pct,
        sum_cov_fracs,
    )

    if not source_totals:
        return
    ref = next(iter(source_totals.values()))
    nbin = int(np.asarray(ref).shape[0])
    centers, bins = _centers_bins(vc, nbin)
    is_int = getattr(vc, "var_save_name", "") == "integrated" or nbin == 1
    grand = sum_cov_fracs(list(source_totals.values()))

    for src in source_order:
        if src not in source_totals:
            continue
        w = frac_unc_pct(source_totals[src])
        if is_int:
            w = np.full_like(w, float(w[0]))
        ax.hist(
            centers, bins=bins, weights=w, histtype="step", linewidth=1.8,
            color=SOURCE_COLORS.get(src, "gray"),
            label=SOURCE_DISPLAY.get(src, src),
        )
    if grand is not None:
        wtot = frac_unc_pct(grand)
        if is_int:
            wtot = np.full_like(wtot, float(wtot[0]))
        ax.hist(
            centers, bins=bins, weights=wtot, histtype="step",
            linewidth=2.6, color="k", label="Total",
        )

    xlab = vc.var_labels[1] if vc is not None and getattr(vc, "var_labels", None) else ""
    if is_int:
        ax.set_xlabel("All Events")
        ax.set_xticks([centers[0]])
        ax.set_xticklabels(["All Events"])
    else:
        ax.set_xlabel(xlab or "bin")
    ax.set_ylabel("Uncertainty [%]")
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="best")


def plot_total_source_section(
    payloads: Mapping[str, Any],
    var_configs: Sequence[Any],
    *,
    kinds: Sequence[str] = ("rate", "xsec"),
    out_dir: Optional[Path | str] = None,
    save_figs: bool = True,
    dpi: int = 140,
    include_flat: bool = True,
    plot_source_heatmaps: bool = True,
    plot_total_heatmaps: bool = True,
):
    """For each kind × variable: source unc plot; optional per-source + total heatmaps."""
    out_dir = Path(out_dir) if out_dir else None
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)

    for kind in kinds:
        order = SOURCE_ORDER_XSEC if kind == "xsec" else SOURCE_ORDER_RATE
        for vc in var_configs:
            vsn = vc.var_save_name
            totals = build_source_totals(payloads, vsn, kind, include_flat=include_flat)
            if not totals:
                print(f"skip {kind}/{vsn}: no source covs")
                continue

            fig, ax = plt.subplots(figsize=(6.4, 4.8))
            plot_source_frac_unc(ax, totals, vc=vc, kind=kind, source_order=order)
            fig.tight_layout()
            if save_figs and out_dir is not None:
                p = out_dir / f"total_sources__{kind}__{vsn}.png"
                fig.savefig(p, dpi=dpi, bbox_inches="tight")
                print("wrote", p)
            plt.show()

            grand = sum_cov_fracs(list(totals.values()))
            if plot_total_heatmaps and grand is not None:
                save = None
                if save_figs and out_dir is not None:
                    save = out_dir / f"total_covcorr__{kind}__{vsn}__Total.png"
                show_cov_corr_heatmaps(
                    grand, vc, kind=kind, suptitle="Total",
                    save_path=save, dpi=dpi,
                )

            if plot_source_heatmaps:
                for src in order:
                    if src not in totals:
                        continue
                    save = None
                    if save_figs and out_dir is not None:
                        save = out_dir / f"total_covcorr__{kind}__{vsn}__{src}.png"
                    show_cov_corr_heatmaps(
                        totals[src], vc, kind=kind,
                        suptitle=SOURCE_DISPLAY.get(src, src),
                        save_path=save, dpi=dpi,
                    )
