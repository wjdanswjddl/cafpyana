"""Smoothers for noisy unisim fractional uncertainties (DENT Product B).

Ported verbatim from ``notebooks/dent.ipynb`` §"Exploring uncertainty smoothing"
so the *selected* shortlist there (``TH1::Smooth``, Gaussian σ=1.5, moving
average w=3, local polynomial / Savitzky–Golay w=7 p=2, rebin×2 then upsample)
can be applied reproducibly by scripts.

All smoothers act on the per-bin fractional uncertainty ``u_i = σ_i / N_i^CV``
(``|N^var − N^CV| / N^CV`` for a unisim) and return an array of the same length.
A smoothed unisim covariance is rebuilt as the fully correlated
``cov_frac = outer(u_s, u_s)`` — the same rank-1 structure as the raw unisim.

Nothing here is used by the DENT producer; the producer pack is the raw unisim.
The Product B *consumer* nominal DENT (since 2026-09-29) is ``upper80_w3``
(rolling 80% w=3 + Gaussian σ=1); the unsmoothed unisim is kept as
``productB_sel_mup__dent_raw``. See ``scripts/build_dent_smooth_test_trees.py``.
"""

from __future__ import annotations

from typing import Callable, Dict

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter

# Tunables (bin units) — the values chosen in notebooks/dent.ipynb.
GAUSS_SIGMA_BINS = 1.5
MOVING_AVG_WIDTH = 3  # odd
SAVGOL_WINDOW = 7  # odd
SAVGOL_POLY = 2
TH1_SMOOTH_NTIMES = 1  # ROOT TH1::Smooth(ntimes) default
REBIN_FACTOR = 2


def _odd(n: int) -> int:
    n = max(3, int(n))
    return n if n % 2 else n + 1


def _nan_safe_smooth(y: np.ndarray, smoother: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
    """Interpolate over non-finite bins, smooth, then mask them back to NaN.

    Fewer than 3 finite bins → all-NaN (caller decides; the tree builder keeps raw).
    """
    y = np.asarray(y, dtype=float)
    out = np.full_like(y, np.nan, dtype=float)
    m = np.isfinite(y)
    if m.sum() < 3:
        return out
    x = np.arange(len(y), dtype=float)
    y_fill = np.interp(x, x[m], y[m])
    out[:] = smoother(y_fill)
    out[~m] = np.nan
    return out


def frac_from_hists(n_cv: np.ndarray, n_var: np.ndarray, *, eps: float = 1e-12) -> np.ndarray:
    """Raw unisim fractional unc (absolute relative shift)."""
    n_cv = np.asarray(n_cv, dtype=float)
    n_var = np.asarray(n_var, dtype=float)
    return np.abs(n_var - n_cv) / np.maximum(n_cv, eps)


def smooth_moving_avg(y: np.ndarray, width: int = MOVING_AVG_WIDTH) -> np.ndarray:
    w = _odd(width)
    kernel = np.ones(w, dtype=float) / w
    return _nan_safe_smooth(y, lambda yf: np.convolve(yf, kernel, mode="same"))


def smooth_gauss(y: np.ndarray, sigma: float = GAUSS_SIGMA_BINS) -> np.ndarray:
    return _nan_safe_smooth(
        y, lambda yf: gaussian_filter1d(yf, sigma=float(sigma), mode="nearest")
    )


def smooth_savgol(y: np.ndarray, window: int = SAVGOL_WINDOW, poly: int = SAVGOL_POLY) -> np.ndarray:
    """Local polynomial (Savitzky–Golay); falls back to Gaussian for short arrays."""
    w = _odd(window)
    y = np.asarray(y, dtype=float)
    if len(y) < w:
        return smooth_gauss(y)

    def _fn(yf):
        wl = min(w, len(yf))
        if wl % 2 == 0:
            wl -= 1
        if wl < 3:
            return yf
        return savgol_filter(yf, window_length=wl, polyorder=min(poly, wl - 1))

    return _nan_safe_smooth(y, _fn)


def _median_n(vals: np.ndarray) -> float:
    return float(np.median(np.asarray(vals, dtype=float)))


def smooth_array_353qh(xx: np.ndarray, ntimes: int = 1) -> np.ndarray:
    """ROOT ``TH1::SmoothArray`` — Friedman 353QH applied twice per pass.

    Port of ``hist/hist/src/TH1.cxx`` (HBOOK ``hsmoof.F``):
      * 353 — running medians of width 3, then 5, then 3
      * Q   — quadratic interpolation on flat 3-bin plateaus
      * H   — Hanning running mean (0.25, 0.5, 0.25)
      * twice — repeat on residuals and add back (Tukey)

    If the input minimum is ≥ 0 (e.g. fractional unc), negatives are clipped to 0
    as in ROOT.
    """
    xx = np.asarray(xx, dtype=float).copy()
    nn = int(xx.size)
    if nn < 3:
        return xx
    ntimes = max(1, int(ntimes))

    yy = np.empty(nn, dtype=float)
    zz = np.empty(nn, dtype=float)
    rr = np.empty(nn, dtype=float)

    for _pass in range(ntimes):
        zz[:] = xx
        for noent in range(2):  # algorithm twice: data, then residuals
            for kk in range(3):  # medians 3, 5, 3
                yy[:] = zz
                median_type = 3 if kk != 1 else 5
                ifirst = 1 if kk != 1 else 2
                ilast = nn - 1 if kk != 1 else nn - 2
                for ii in range(ifirst, ilast):
                    zz[ii] = _median_n(yy[ii - ifirst : ii - ifirst + median_type])

                if kk == 0:  # first median 3 — special endpoints
                    hh = np.array([zz[1], zz[0], 3 * zz[1] - 2 * zz[2]], dtype=float)
                    zz[0] = _median_n(hh)
                    hh = np.array(
                        [zz[nn - 2], zz[nn - 1], 3 * zz[nn - 2] - 2 * zz[nn - 3]],
                        dtype=float,
                    )
                    zz[nn - 1] = _median_n(hh)

                if kk == 1:  # median 5 — near-edge width-3
                    zz[1] = _median_n(yy[0:3])
                    zz[nn - 2] = _median_n(yy[nn - 3 : nn])

            yy[:] = zz

            # quadratic interpolation for flat segments
            for ii in range(2, nn - 2):
                if zz[ii - 1] != zz[ii] or zz[ii] != zz[ii + 1]:
                    continue
                tmp0 = zz[ii - 2] - zz[ii]
                tmp1 = zz[ii + 2] - zz[ii]
                if tmp0 * tmp1 <= 0:
                    continue
                jk = 1 if abs(tmp1) <= abs(tmp0) else -1
                yy[ii] = -0.5 * zz[ii - 2 * jk] + zz[ii] / 0.75 + zz[ii + 2 * jk] / 6.0
                yy[ii + jk] = 0.5 * (zz[ii + 2 * jk] - zz[ii - 2 * jk]) + zz[ii]

            # Hanning running means
            for ii in range(1, nn - 1):
                zz[ii] = 0.25 * yy[ii - 1] + 0.5 * yy[ii] + 0.25 * yy[ii + 1]
            zz[0] = yy[0]
            zz[nn - 1] = yy[nn - 1]

            if noent == 0:
                rr[:] = zz
                zz[:] = xx - zz  # residuals for second pass

        xmin = float(np.min(xx))
        if xmin < 0:
            xx[:] = rr + zz
        else:
            xx[:] = np.maximum(rr + zz, 0.0)
    return xx


def smooth_th1(y: np.ndarray, ntimes: int = TH1_SMOOTH_NTIMES) -> np.ndarray:
    """Nan-safe wrapper around Friedman 353QH (ROOT ``TH1::Smooth``)."""
    return _nan_safe_smooth(y, lambda yf: smooth_array_353qh(yf, ntimes=ntimes))


def rebin_frac_unc_then_upsample(frac: np.ndarray, factor: int = REBIN_FACTOR) -> np.ndarray:
    """Average frac unc in blocks of ``factor`` bins, then expand piecewise-constant."""
    f = np.asarray(frac, dtype=float)
    n = len(f)
    fac = max(1, int(factor))
    n_pad = int(np.ceil(n / fac) * fac)
    pad = np.full(n_pad, np.nan, dtype=float)
    pad[:n] = f
    blocks = pad.reshape(-1, fac)
    with np.errstate(all="ignore"):
        means = np.nanmean(blocks, axis=1)
    return np.repeat(means, fac)[:n]


# Conservative (upward-biased) smoother: rolling high percentile, then a light Gaussian.
UPPER80_WINDOW = 5  # odd
UPPER80_WINDOW_W3 = 3  # odd; narrower envelope than UPPER80_WINDOW
UPPER80_Q = 0.8
UPPER80_GAUSS_SIGMA = 1.0


def rolling_quantile(y: np.ndarray, window: int = UPPER80_WINDOW, q: float = UPPER80_Q) -> np.ndarray:
    """Nan-aware centered rolling quantile (edges use a truncated window)."""
    y = np.asarray(y, dtype=float)
    w = _odd(window)
    half = w // 2
    out = np.empty_like(y)
    for i in range(len(y)):
        lo = max(0, i - half)
        hi = min(len(y), i + half + 1)
        out[i] = float(np.nanquantile(y[lo:hi], q))
    return out


def smooth_upper80(
    y: np.ndarray,
    window: int = UPPER80_WINDOW,
    q: float = UPPER80_Q,
    sigma: float = UPPER80_GAUSS_SIGMA,
) -> np.ndarray:
    """Rolling ``q``-quantile (default 80%) then Gaussian σ (default 1 bin).

    Biased high relative to a two-sided mean/savgol: valleys come up, spikes are
    only partly shaved. Intended as a conservative DENT diagonal.
    """

    def _fn(yf):
        qq = rolling_quantile(yf, window=window, q=q)
        if len(yf) >= 3 and sigma > 0:
            qq = gaussian_filter1d(qq, sigma=float(sigma), mode="nearest")
        return qq

    return _nan_safe_smooth(y, _fn)


# The shortlist from notebooks/dent.ipynb §"Selected smoothing comparison".
# key -> (human label, callable on u = σ/N)
SELECTED_SMOOTHERS: Dict[str, tuple] = {
    "th1": (f"TH1::Smooth (353QH×2, ntimes={TH1_SMOOTH_NTIMES})", lambda u: smooth_th1(u)),
    "gauss15": (f"Gaussian σ={GAUSS_SIGMA_BINS} bins on σ/N", lambda u: smooth_gauss(u)),
    "mavg3": (f"moving average w={_odd(MOVING_AVG_WIDTH)}", lambda u: smooth_moving_avg(u)),
    "savgol": (f"local polynomial (Savitzky–Golay w={_odd(SAVGOL_WINDOW)}, p={SAVGOL_POLY})", lambda u: smooth_savgol(u)),
    "rebin2": (f"rebin×{REBIN_FACTOR} then upsample", lambda u: rebin_frac_unc_then_upsample(u)),
}

CONSERVATIVE_SMOOTHERS: Dict[str, tuple] = {
    "upper80": (
        f"rolling {int(100 * UPPER80_Q)}% (w={_odd(UPPER80_WINDOW)}) + Gaussian σ={UPPER80_GAUSS_SIGMA:g}",
        lambda u: smooth_upper80(u),
    ),
    "upper80_w3": (
        f"rolling {int(100 * UPPER80_Q)}% (w={_odd(UPPER80_WINDOW_W3)}) + Gaussian σ={UPPER80_GAUSS_SIGMA:g}",
        lambda u: smooth_upper80(u, window=UPPER80_WINDOW_W3),
    ),
}

ALL_SMOOTHERS: Dict[str, tuple] = {**SELECTED_SMOOTHERS, **CONSERVATIVE_SMOOTHERS}


def smooth_unisim_frac(u_raw: np.ndarray, method: str) -> np.ndarray:
    """Apply a selected smoother to σ/N. Returns raw for <3 bins; clips to ≥0;
    keeps raw values where the smoother yields NaN."""
    u_raw = np.asarray(u_raw, dtype=float)
    if u_raw.size < 3:
        return u_raw.copy()
    fn = ALL_SMOOTHERS[method][1]
    u_s = np.asarray(fn(u_raw), dtype=float)
    bad = ~np.isfinite(u_s)
    u_s[bad] = u_raw[bad]
    return np.maximum(u_s, 0.0)
