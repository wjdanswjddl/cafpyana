#!/usr/bin/env python
"""May vs Sep MEC: universe overlays on CV + fractional uncertainty plots.

One pass over each chunk campaign, then write a multi-page PDF.
"""
from __future__ import annotations

import csv
import glob
import pickle
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from analysis_village.numucc_1p0pi.scripts.get_systematics_genie import (  # noqa: E402
    genie_final_var_configs,
)
from analysis_village.numucc_1p0pi.utils import get_response_matrix  # noqa: E402
from pyanalib.covariance import get_covariance_matrix  # noqa: E402

MAY_CHUNKS = sorted(
    glob.glob(
        "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
        "genie_syst-chunked-20260514/chunks/genie__Ar23p__*.pkl"
    )
)
SEP_CHUNKS = sorted(
    glob.glob(
        "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
        "genie_syst-chunked-sel_mup_20260912_Ar23p/chunks/genie__Ar23p__*.pkl"
    )
)
MAY_NPZ = Path(
    "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/"
    "systematics-genie-final/systematics-genie-Ar23p-final/genie-Ar23p_syst_dict.npz"
)
OUT_DIR = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/mec_univ_may_vs_sep")
CACHE = OUT_DIR / "agg_cache.pkl"

VARS = ["tki-del_Tp", "tki-del_alpha", "tki-del_phi"]
KNOB_PAIRS = [
    (
        "Val q0bin0",
        "MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin0",
        "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin0",
    ),
    (
        "Val q0bin1",
        "MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin1",
        "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin1",
    ),
    (
        "Val q0bin2",
        "MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin2",
        "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin2",
    ),
    (
        "Val q0bin3",
        "MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin3",
        "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin3",
    ),
    (
        "Mar q0bin0",
        "MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin0",
        "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin0",
    ),
    (
        "Mar q0bin1",
        "MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin1",
        "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin1",
    ),
    (
        "Mar q0bin2",
        "MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin2",
        "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin2",
    ),
    (
        "Mar q0bin3",
        "MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin3",
        "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin3",
    ),
]


def finalize_xsec(acc: dict) -> tuple[np.ndarray, np.ndarray]:
    N = acc["nevts_allmc"]
    n_univ = int(acc["reco_vs_true"].shape[0])
    rows = []
    for u in range(n_univ):
        rvt = acc["reco_vs_true"][u]
        allmc = acc["signal_allmc"][u]
        sel = acc["signal_sel_truth"][u]
        eff = np.divide(sel, allmc, out=np.zeros_like(sel), where=allmc != 0)
        R = get_response_matrix(rvt, eff)
        rows.append(R @ N + (acc["bg_univ"][u] - acc["bg_cv"]))
    return np.asarray(rows, dtype=np.float64), np.asarray(acc["cv_sel_reco"], dtype=np.float64)


def aggregate_campaign(chunks: list[str], knobs: list[str], vars_: list[str], tag: str) -> dict:
    """Single pass: return {(knob, var): {rate_univ, rate_cv, xsec_univ, xsec_cv, n_chunks}}."""
    rate_u: dict[tuple[str, str], np.ndarray] = {}
    rate_c: dict[tuple[str, str], np.ndarray] = {}
    accs: dict[tuple[str, str], dict] = {}
    counts: dict[tuple[str, str], int] = {}
    knob_set = set(knobs)
    var_set = set(vars_)

    for i, p in enumerate(chunks):
        if i % 100 == 0:
            print(f"  [{tag}] {i}/{len(chunks)}", flush=True)
        with open(p, "rb") as f:
            d = pickle.load(f)
        for sk in knob_set:
            if sk not in d["rate_univ_cv"]:
                continue
            for var in var_set:
                if var not in d["rate_univ_cv"][sk]:
                    continue
                key = (sk, var)
                r = d["rate_univ_cv"][sk][var]
                a = d["xsec_accumulators"][sk][var]
                if key not in rate_u:
                    rate_u[key] = np.asarray(r["univ"], dtype=np.float64).copy()
                    rate_c[key] = np.asarray(r["cv"], dtype=np.float64).copy()
                    accs[key] = {
                        k: np.zeros_like(np.asarray(v), dtype=np.float64) for k, v in a.items()
                    }
                    counts[key] = 0
                else:
                    rate_u[key] = rate_u[key] + np.asarray(r["univ"], dtype=np.float64)
                    rate_c[key] = rate_c[key] + np.asarray(r["cv"], dtype=np.float64)
                for k in accs[key]:
                    accs[key][k] += np.asarray(a[k], dtype=np.float64)
                counts[key] += 1

    out = {}
    for key in rate_u:
        xu, xc = finalize_xsec(accs[key])
        out[key] = {
            "rate_univ": rate_u[key],
            "rate_cv": rate_c[key],
            "xsec_univ": xu,
            "xsec_cv": xc,
            "n_chunks": counts[key],
        }
    print(f"  [{tag}] done, {len(out)} (knob,var) packs", flush=True)
    return out


def frac_pct(univ: np.ndarray, cv: np.ndarray) -> np.ndarray:
    pack = get_covariance_matrix(univ, cv)
    return 100.0 * np.sqrt(np.clip(np.diag(pack["cov_frac"]), 0, None))


def npz_pct(npz: dict, sk_old: str, var: str, kind: str) -> np.ndarray | None:
    if sk_old not in npz:
        return None
    o = npz[sk_old]
    if isinstance(o, np.ndarray) and o.dtype == object:
        o = o.item()
    if var not in o:
        return None
    d = o[var][kind]
    if isinstance(d, np.ndarray) and d.dtype == object:
        d = d.item()
    return 100.0 * np.sqrt(np.clip(np.diag(np.asarray(d["cov_frac"])), 0, None))


def npz_implied_univ_cv(npz: dict, sk_old: str, var: str, kind: str):
    """Recover CV and |delta| from published cov (n_univ=1 assumption for signless overlay)."""
    if sk_old not in npz:
        return None
    o = npz[sk_old]
    if isinstance(o, np.ndarray) and o.dtype == object:
        o = o.item()
    if var not in o:
        return None
    d = o[var][kind]
    if isinstance(d, np.ndarray) and d.dtype == object:
        d = d.item()
    cov = np.asarray(d["cov"], dtype=float)
    cf = np.asarray(d["cov_frac"], dtype=float)
    diag_c = np.diag(cov)
    diag_f = np.diag(cf)
    cv = np.full_like(diag_c, np.nan)
    m = diag_f > 0
    cv[m] = np.sqrt(diag_c[m] / diag_f[m])
    # for zero-frac bins, leave nan; delta = ±sqrt(cov)
    delta = np.sqrt(np.clip(diag_c, 0, None))
    return cv, delta


def centers(bins):
    b = np.asarray(bins, float)
    return 0.5 * (b[:-1] + b[1:])


def step_xy(bins, y):
    x = np.asarray(bins, float)
    yy = np.asarray(y, float)
    return x, np.r_[yy, yy[-1]]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    vcs = {vc.var_save_name: vc for vc in genie_final_var_configs()}
    may_knobs = [p[1] for p in KNOB_PAIRS]
    sep_knobs = [p[2] for p in KNOB_PAIRS]

    if CACHE.is_file():
        print(f"loading cache {CACHE}", flush=True)
        with open(CACHE, "rb") as f:
            cache = pickle.load(f)
    else:
        print(f"May chunks: {len(MAY_CHUNKS)}", flush=True)
        may = aggregate_campaign(MAY_CHUNKS, may_knobs, VARS, "may")
        print(f"Sep chunks: {len(SEP_CHUNKS)}", flush=True)
        sep = aggregate_campaign(SEP_CHUNKS, sep_knobs, VARS, "sep")
        cache = {"may": may, "sep": sep}
        with open(CACHE, "wb") as f:
            pickle.dump(cache, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"wrote {CACHE}", flush=True)

    may = cache["may"]
    sep = cache["sep"]
    npz = dict(np.load(MAY_NPZ, allow_pickle=True))

    pdf_path = OUT_DIR / "mec_univ_cv_may_vs_sep.pdf"
    summary_rows = []

    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.axis("off")
        ax.text(
            0.02,
            0.95,
            "MEC q0-interp knobs: universe vs CV and fractional uncertainty\n\n"
            "May = genie_syst-chunked-20260514 (Ar23p), n_univ=1 (ps1)\n"
            "Sep = genie_syst-chunked-sel_mup_20260912_Ar23p, n_univ=2 (ps1, ms1)\n\n"
            "Rate: selected reco yield universes\n"
            "Xsec: R_univ @ N_gen^CV (+ bg), CV = N_sel_reco\n\n"
            "x markers: published May19 NPZ fractional unc\n"
            "(NPZ sample can differ from May14 chunk aggregate).",
            transform=ax.transAxes,
            va="top",
            fontsize=12,
            family="monospace",
        )
        pdf.savefig(fig)
        plt.close(fig)

        for label, sk_m, sk_s in KNOB_PAIRS:
            for var in VARS:
                vc = vcs[var]
                bins = np.asarray(vc.bins, float)
                cx = centers(bins)
                may_d = may.get((sk_m, var))
                sep_d = sep.get((sk_s, var))
                if may_d is None or sep_d is None:
                    print("SKIP", label, var, may_d is None, sep_d is None)
                    continue

                may_r_pct = frac_pct(may_d["rate_univ"], may_d["rate_cv"])
                may_x_pct = frac_pct(may_d["xsec_univ"], may_d["xsec_cv"])
                sep_r_pct = frac_pct(sep_d["rate_univ"], sep_d["rate_cv"])
                sep_x_pct = frac_pct(sep_d["xsec_univ"], sep_d["xsec_cv"])
                npz_r = npz_pct(npz, sk_m, var, "rate")
                npz_x = npz_pct(npz, sk_m, var, "xsec")

                summary_rows.append(
                    {
                        "knob": label,
                        "var": var,
                        "may_rate_last": float(may_r_pct[-1]),
                        "may_xsec_last": float(may_x_pct[-1]),
                        "sep_rate_last": float(sep_r_pct[-1]),
                        "sep_xsec_last": float(sep_x_pct[-1]),
                        "npz_rate_last": float(npz_r[-1]) if npz_r is not None else np.nan,
                        "npz_xsec_last": float(npz_x[-1]) if npz_x is not None else np.nan,
                        "may_cv_last": float(may_d["rate_cv"][-1]),
                        "sep_cv_last": float(sep_d["rate_cv"][-1]),
                    }
                )

                fig, axes = plt.subplots(
                    3,
                    2,
                    figsize=(12, 10),
                    sharex="col",
                    gridspec_kw={"height_ratios": [1.2, 1.2, 1.0]},
                )
                fig.suptitle(f"{label}  |  {var}\n{vc.var_labels[0]}", fontsize=13)

                for col, (data, title) in enumerate(
                    [(may_d, "May rate"), (sep_d, "Sep rate")]
                ):
                    ax = axes[0, col]
                    x, ycv = step_xy(bins, data["rate_cv"])
                    ax.step(x, ycv, where="post", color="k", lw=1.8, label="CV")
                    for u in range(data["rate_univ"].shape[0]):
                        _, yu = step_xy(bins, data["rate_univ"][u])
                        ax.step(
                            x,
                            yu,
                            where="post",
                            color="C0",
                            alpha=0.55,
                            lw=1.2,
                            label=f"univ {u}",
                        )
                    ax.set_ylabel("events")
                    ax.set_title(title)
                    ax.legend(loc="best", fontsize=8)
                    ax.grid(True, alpha=0.3)

                for col, (data, title) in enumerate(
                    [(may_d, "May xsec"), (sep_d, "Sep xsec")]
                ):
                    ax = axes[1, col]
                    x, ycv = step_xy(bins, data["xsec_cv"])
                    ax.step(x, ycv, where="post", color="k", lw=1.8, label="CV")
                    for u in range(data["xsec_univ"].shape[0]):
                        _, yu = step_xy(bins, data["xsec_univ"][u])
                        ax.step(
                            x,
                            yu,
                            where="post",
                            color="C3",
                            alpha=0.55,
                            lw=1.2,
                            label=f"univ {u}",
                        )
                    ax.set_ylabel("arb. (xsec_unit=1)")
                    ax.set_title(title)
                    ax.legend(loc="best", fontsize=8)
                    ax.grid(True, alpha=0.3)

                ax = axes[2, 0]
                ax.plot(cx, may_r_pct, "o-", color="C0", label="May14 rate")
                ax.plot(cx, sep_r_pct, "s--", color="C0", alpha=0.75, label="Sep12 rate")
                if npz_r is not None and len(npz_r) == len(cx):
                    ax.plot(cx, npz_r, "x", color="navy", ms=7, label="May19 NPZ rate")
                ax.set_ylabel("frac unc [%]")
                ax.set_xlabel(vc.var_labels[1])
                ax.set_title("Rate fractional uncertainty")
                ax.legend(fontsize=7, loc="best")
                ax.grid(True, alpha=0.3)
                ax.set_ylim(bottom=0)

                ax = axes[2, 1]
                ax.plot(cx, may_x_pct, "o-", color="C3", label="May14 xsec")
                ax.plot(cx, sep_x_pct, "s--", color="C3", alpha=0.75, label="Sep12 xsec")
                if npz_x is not None and len(npz_x) == len(cx):
                    ax.plot(cx, npz_x, "x", color="darkred", ms=7, label="May19 NPZ xsec")
                ax.set_ylabel("frac unc [%]")
                ax.set_xlabel(vc.var_labels[1])
                ax.set_title("Xsec fractional uncertainty")
                ax.legend(fontsize=7, loc="best")
                ax.grid(True, alpha=0.3)
                ax.set_ylim(bottom=0)

                txt = (
                    f"last bin: May rate {may_r_pct[-1]:.2f}%  xsec {may_x_pct[-1]:.2f}%   |   "
                    f"Sep rate {sep_r_pct[-1]:.2f}%  xsec {sep_x_pct[-1]:.2f}%"
                )
                if npz_x is not None:
                    txt += f"   |   NPZ xsec {npz_x[-1]:.2f}%"
                fig.text(0.5, 0.01, txt, ha="center", fontsize=9)
                fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                pdf.savefig(fig)
                fig.savefig(OUT_DIR / f"{label.replace(' ', '_')}_{var}.png", dpi=120)
                plt.close(fig)

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis("off")
        lines = [
            "Last-bin fractional uncertainty summary\n",
            f"{'knob':12s} {'var':14s} {'MayR':7s} {'MayX':7s} {'SepR':7s} {'SepX':7s} "
            f"{'NPZr':7s} {'NPZx':7s} {'MayCV':7s} {'SepCV':7s}",
        ]
        for r in summary_rows:
            lines.append(
                f"{r['knob']:12s} {r['var']:14s} "
                f"{r['may_rate_last']:7.2f} {r['may_xsec_last']:7.2f} "
                f"{r['sep_rate_last']:7.2f} {r['sep_xsec_last']:7.2f} "
                f"{r['npz_rate_last']:7.2f} {r['npz_xsec_last']:7.2f} "
                f"{r['may_cv_last']:7.1f} {r['sep_cv_last']:7.1f}"
            )
        ax.text(
            0.02,
            0.98,
            "\n".join(lines),
            transform=ax.transAxes,
            va="top",
            fontsize=8,
            family="monospace",
        )
        pdf.savefig(fig)
        plt.close(fig)

    csv_path = OUT_DIR / "lastbin_summary.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)
    print("wrote", pdf_path)
    print("wrote", csv_path)


if __name__ == "__main__":
    main()
