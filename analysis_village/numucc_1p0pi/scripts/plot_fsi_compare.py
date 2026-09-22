#!/usr/bin/env python
"""FSI compare plots: retired v1 _N vs v3 packs, and three slim totals.

Reads ``genie_syst_FSI_compare.npz`` from Product B (``final`` / sel_mup) and
Product A (``sel_all``) merge dirs produced by ``run_syst_genie_chunked.sh``
with ``GENIE_RUN_GROUPS=FSI_compare``.

Outputs under ``OUT_DIR`` (PDF + per-page PNG):
  * ``fsi_packs_<product>_<kind>.pdf`` / ``…_<var>.png`` — ``FSI_v1_N`` vs ``FSI_v3_N``
  * ``slim_totals_<product>_<kind>.pdf`` / ``…_<var>.png`` — ``GENIE_slim_{v1,v3,both}``
  * ``fsi_atomics_<product>_<kind>.pdf`` / ``…_<family>_<var>.png`` — per-knob atomics
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from makedf.geniesyst import (  # noqa: E402
    fsi_v1_n_genie_systematics,
    fsi_v3_n_genie_systematics,
)

DEFAULT_B = Path(
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "genie_syst-chunked-FSI_compare_sel_mup/merged/FSI_compare/genie_syst_FSI_compare.npz"
)
DEFAULT_A = Path(
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "genie_syst-chunked-FSI_compare_sel_all/merged/FSI_compare/genie_syst_FSI_compare.npz"
)
OUT_DIR = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/FSI_compare")

PACK_PAIR = ("FSI_v1_N", "FSI_v3_N")
SLIM_PACKS = ("GENIE_slim_v1", "GENIE_slim_v3", "GENIE_slim_both")
PACK_COLORS = {"FSI_v1_N": "C0", "FSI_v3_N": "C1"}
SLIM_COLORS = {
    "GENIE_slim_v1": "C0",
    "GENIE_slim_v3": "C1",
    "GENIE_slim_both": "C3",
}
PACK_LABELS = {
    "FSI_v1_N": "retired FSI v1 (_N pack)",
    "FSI_v3_N": "FSI v3 nucleon pack",
    "GENIE_slim_v1": "slim = base × FSI_v1",
    "GENIE_slim_v3": "slim = base × FSI_v3",
    "GENIE_slim_both": "slim = base × FSI_v1 × FSI_v3",
}


def _frac_pct(pack: dict, kind: str) -> np.ndarray | None:
    block = pack.get(kind)
    if not isinstance(block, dict) or "cov_frac" not in block:
        return None
    cov = np.asarray(block["cov_frac"], dtype=np.float64)
    return np.sqrt(np.clip(np.diag(cov), 0.0, None)) * 100.0


def _load_syst(path: Path) -> dict:
    z = np.load(path, allow_pickle=True)
    return z["syst"].item()


def _vars_for(syst: dict, knobs: list[str], kind: str) -> list[str]:
    found: set[str] = set()
    for k in knobs:
        row = syst.get(k)
        if not isinstance(row, dict):
            continue
        for var, pack in row.items():
            if _frac_pct(pack, kind) is not None:
                found.add(var)
    # Prefer common analysis vars first
    prefer = [
        "integrated",
        "tki-del_Tp",
        "tki-del_alpha",
        "tki-del_phi",
        "muon-p",
        "muon-dir_z",
        "proton-p",
        "mcs_range_diff",
        "nu_score",
        "n_trks",
        "trk_len",
        "chi2_mu",
        "chi2_p",
    ]
    ordered = [v for v in prefer if v in found]
    ordered += sorted(v for v in found if v not in ordered)
    return ordered


def _short_knob(name: str) -> str:
    for prefix in (
        "GENIEReWeight_SBN_v1_multisigma_",
        "GENIEReWeight_SBN_v3_",
        "PionAbsWeighter_SBN_v3_",
    ):
        if name.startswith(prefix):
            return name[len(prefix) :]
    return name


def _safe_slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", name).strip("_")


def _savefig(fig, pdf: PdfPages | None, png: Path | None) -> None:
    if pdf is not None:
        pdf.savefig(fig)
    if png is not None:
        png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(png, dpi=150)
    plt.close(fig)


def _plot_overlay(
    *,
    title: str,
    series: list[tuple[str, np.ndarray, str]],
    ylabel: str = r"fractional uncertainty [%]",
    pdf: PdfPages | None = None,
    png: Path | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    for label, y, color in series:
        x = np.arange(len(y))
        ax.step(x, y, where="mid", label=label, color=color, lw=1.6)
        ax.fill_between(x, 0.0, y, step="mid", color=color, alpha=0.15)
    ax.set_xlabel("bin index")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(bottom=0.0)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    _savefig(fig, pdf, png)


def write_pack_plots(syst: dict, product: str, kind: str, out_pdf: Path) -> list[Path]:
    knobs = list(PACK_PAIR)
    vars_ = _vars_for(syst, knobs, kind)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = [out_pdf]
    with PdfPages(out_pdf) as pdf:
        for var in vars_:
            series = []
            for k in knobs:
                pack = syst.get(k, {}).get(var)
                if not pack:
                    continue
                y = _frac_pct(pack, kind)
                if y is None:
                    continue
                series.append((PACK_LABELS.get(k, k), y, PACK_COLORS.get(k, "k")))
            if not series:
                continue
            png = out_pdf.with_name(f"{out_pdf.stem}_{_safe_slug(var)}.png")
            _plot_overlay(
                title=f"{product} · {kind} · {var}\nFSI_v1_N vs FSI_v3_N packs",
                series=series,
                pdf=pdf,
                png=png,
            )
            written.append(png)
    return written


def write_slim_plots(syst: dict, product: str, kind: str, out_pdf: Path) -> list[Path]:
    knobs = list(SLIM_PACKS)
    vars_ = _vars_for(syst, knobs, kind)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = [out_pdf]
    with PdfPages(out_pdf) as pdf:
        for var in vars_:
            series = []
            for k in knobs:
                pack = syst.get(k, {}).get(var)
                if not pack:
                    continue
                y = _frac_pct(pack, kind)
                if y is None:
                    continue
                series.append((PACK_LABELS.get(k, k), y, SLIM_COLORS.get(k, "k")))
            if not series:
                continue
            png = out_pdf.with_name(f"{out_pdf.stem}_{_safe_slug(var)}.png")
            _plot_overlay(
                title=f"{product} · {kind} · {var}\nslim totals (base × FSI variants)",
                series=series,
                pdf=pdf,
                png=png,
            )
            written.append(png)
    return written


def write_atomics_plots(syst: dict, product: str, kind: str, out_pdf: Path) -> list[Path]:
    """One page per atomic family×var (v1 then v3). Skip if no atomics in NPZ."""
    knobs = list(fsi_v1_n_genie_systematics) + list(fsi_v3_n_genie_systematics)
    if not any(k in syst for k in knobs):
        print(f"  [skip atomics] no atomic knobs in NPZ for {product}/{kind}")
        return []
    vars_ = _vars_for(syst, knobs, kind)
    prefer = [v for v in ("tki-del_Tp", "tki-del_alpha", "tki-del_phi", "integrated", "mcs_range_diff") if v in vars_]
    if not prefer:
        prefer = vars_[:4]
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = [out_pdf]
    with PdfPages(out_pdf) as pdf:
        for family, members, color, fam_slug in (
            ("v1 retired _N", list(fsi_v1_n_genie_systematics), "C0", "v1"),
            ("v3 nucleon FSI", list(fsi_v3_n_genie_systematics), "C1", "v3"),
        ):
            for var in prefer:
                series = []
                for k in members:
                    pack = syst.get(k, {}).get(var)
                    if not pack:
                        continue
                    y = _frac_pct(pack, kind)
                    if y is None:
                        continue
                    series.append((_short_knob(k), y, color))
                if not series:
                    continue
                n = len(series)
                ncols = 2
                nrows = int(np.ceil(n / ncols))
                fig, axes = plt.subplots(
                    nrows, ncols, figsize=(10, 2.6 * nrows), sharex=True, squeeze=False
                )
                for ax, (label, y, c) in zip(axes.ravel(), series):
                    x = np.arange(len(y))
                    ax.step(x, y, where="mid", color=c, lw=1.4)
                    ax.fill_between(x, 0.0, y, step="mid", color=c, alpha=0.2)
                    ax.set_title(label, fontsize=9)
                    ax.set_ylim(bottom=0.0)
                    ax.grid(True, alpha=0.3)
                for ax in axes.ravel()[n:]:
                    ax.axis("off")
                fig.suptitle(f"{product} · {kind} · {var} · {family} atomics", fontsize=11)
                fig.tight_layout(rect=[0, 0, 1, 0.97])
                png = out_pdf.with_name(
                    f"{out_pdf.stem}_{fam_slug}_{_safe_slug(var)}.png"
                )
                _savefig(fig, pdf, png)
                written.append(png)
    return written


def process_product(npz_path: Path, product: str, out_dir: Path) -> list[Path]:
    if not npz_path.is_file():
        print(f"[skip] missing {npz_path}")
        return []
    syst = _load_syst(npz_path)
    print(f"[load] {product}: {npz_path}  knobs={len(syst)}")
    written: list[Path] = []
    for kind in ("rate", "xsec"):
        p1 = write_pack_plots(
            syst, product, kind, out_dir / f"fsi_packs_{product}_{kind}.pdf"
        )
        p2 = write_slim_plots(
            syst, product, kind, out_dir / f"slim_totals_{product}_{kind}.pdf"
        )
        p3 = write_atomics_plots(
            syst, product, kind, out_dir / f"fsi_atomics_{product}_{kind}.pdf"
        )
        written.extend(p1 + p2 + p3)
        n_png = sum(1 for p in (p1 + p2 + p3) if p.suffix == ".png")
        msg = f"  wrote {p1[0].name}, {p2[0].name}"
        if p3:
            msg += f", {p3[0].name}"
        msg += f" (+{n_png} png)"
        print(msg)
    return written


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--product-b", type=Path, default=DEFAULT_B, help="sel_mup / final NPZ")
    p.add_argument("--product-a", type=Path, default=DEFAULT_A, help="sel_all NPZ")
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    p.add_argument(
        "--only",
        choices=("B", "A", "both"),
        default="both",
        help="Which product(s) to plot",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    if args.only in ("B", "both"):
        written += process_product(args.product_b, "ProductB_sel_mup", args.out_dir)
    if args.only in ("A", "both"):
        written += process_product(args.product_a, "ProductA_sel_all", args.out_dir)
    if not written:
        raise SystemExit("No plots written — are the merge NPZs ready?")
    n_pdf = sum(1 for p in written if p.suffix == ".pdf")
    n_png = sum(1 for p in written if p.suffix == ".png")
    print(f"Done: {n_pdf} PDFs + {n_png} PNGs under {args.out_dir}")


if __name__ == "__main__":
    main()
