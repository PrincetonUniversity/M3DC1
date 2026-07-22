from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def tanhfit(x, a0, a1, a2, a3):
    """IDL tanhfit: 0.5*c*(1+d*(1-x))*(1-tanh(b*(x-a)))."""
    xarr = np.asarray(x, dtype=float)
    t = np.tanh(float(a1) * (xarr - float(a0)))
    s = 0.5 * float(a2) * (1.0 + float(a3) * (1.0 - xarr))
    return s * (1.0 - t)


def tanhfit2(x, a0, a1, a2, a3, a4):
    """IDL tanhfit2: tanhfit plus a fitted floor value."""
    return tanhfit(x, a0, a1, a2, a3) + float(a4)


def _read_two_columns(path: Path) -> tuple[np.ndarray, np.ndarray]:
    arr = np.loadtxt(str(path), dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] < 2:
        raise ValueError(f"{path} must contain at least two numeric columns.")
    return np.asarray(arr[:, 0], dtype=float), np.asarray(arr[:, 1], dtype=float)


def _smooth_1d(values: np.ndarray, width: int) -> np.ndarray:
    n = int(width)
    if n <= 1:
        return np.asarray(values, dtype=float)
    kernel = np.ones(n, dtype=float) / float(n)
    return np.convolve(np.asarray(values, dtype=float), kernel, mode="same")


def _fit_profile(xf: np.ndarray, yf: np.ndarray, minval):
    positive = np.asarray(yf, dtype=float) > 0.0
    sigma = np.sqrt(np.maximum(np.asarray(yf, dtype=float), np.finfo(float).tiny))
    sigma = np.where(positive, sigma, 1.0)

    if minval is None:
        p0 = np.asarray([0.98, 1.0 / 0.01, float(np.max(yf)), 0.0, float(np.min(yf))], dtype=float)
        popt, _ = curve_fit(tanhfit2, xf, yf, p0=p0, sigma=sigma, absolute_sigma=True, maxfev=20000)
    else:
        p0 = np.asarray([0.98, 1.0 / 0.01, float(np.max(yf)), 0.0], dtype=float)
        popt, _ = curve_fit(
            tanhfit,
            xf,
            yf - float(minval),
            p0=p0,
            sigma=sigma,
            absolute_sigma=True,
            maxfev=20000,
        )
    return np.asarray(popt, dtype=float)


def extend_profile(
    filein: str | Path,
    *,
    psimax: float = 1.05,
    fitrange=None,
    minval=None,
    smooth=None,
    suffix: str = ".extended",
    plot: bool = True,
    ax=None,
    **plot_kwargs,
) -> np.ndarray:
    """Extend a two-column profile file using the IDL ``extend_profile.pro`` fit.

    The fit uses
    ``0.5*c*(1+d*(1-x))*(1-tanh(b*(x-a)))`` and either a fitted floor
    parameter or a supplied ``minval``. Extended data are written to
    ``filein + suffix`` using the IDL ``(2F12.6)``-style fixed-width format.
    """
    path = Path(filein)
    x, y = _read_two_columns(path)

    if fitrange is None:
        fitrange_arr = np.asarray([0.95, float(np.max(x))], dtype=float)
    else:
        vals = np.asarray(fitrange, dtype=float).reshape(-1)
        if vals.size == 1:
            fitrange_arr = np.asarray([float(vals[0]), float(np.max(x))], dtype=float)
        elif vals.size >= 2:
            fitrange_arr = np.asarray([float(vals[0]), float(vals[1])], dtype=float)
        else:
            raise ValueError("fitrange must be None, a scalar, or a two-value sequence.")

    print("Fitting points in range ", fitrange_arr)
    if float(psimax) <= float(np.max(x)):
        print("Error: profile already extends to ", psimax)
        return np.asarray([], dtype=float)

    idx = np.nonzero((x >= fitrange_arr[0]) & (x <= fitrange_arr[1]))[0]
    count = int(idx.size)
    print("Fitting to ", count, " points")
    if count == 0:
        print("Error: no data points in range")
        return np.asarray([], dtype=float)

    xf = x[idx]
    yf = y[idx]
    params = _fit_profile(xf, yf, minval)

    print("Fit center: ", params[0])
    print("Fit width: ", 1.0 / params[1])
    print("Fit height: ", params[2])
    if minval is None:
        print("Fit floor: ", params[4])
        if params[4] <= 0.0:
            print("Warning: fit asymptotes to negative values")
            print("Try setting minval = desired asymptotic value")

    print("Extending profile to ", psimax)
    n = int(x.size)
    if n < 2:
        raise ValueError("Profile must contain at least two points to determine extension spacing.")
    dx = float(x[n - 1] - x[n - 2])
    if dx == 0.0:
        raise ValueError("Last two profile x values must be distinct.")
    m = int((float(psimax) - float(x[n - 1])) / dx)

    newx = np.empty(n + m, dtype=float)
    newy = np.empty(n + m, dtype=float)
    newx[:n] = x
    newy[:n] = y
    if m > 0:
        newx[n:] = x[n - 1] + (np.arange(m, dtype=float) + 1.0) * dx
        if minval is None:
            newy[n:] = tanhfit2(newx[n:], *params)
        else:
            newy[n:] = tanhfit(newx[n:], *params) + float(minval)

    if smooth is not None:
        newy = _smooth_1d(newy, int(smooth))

    if plot:
        if ax is None:
            _, ax = plt.subplots()
        ax.plot(newx, newy, **plot_kwargs)
        ax.plot(xf, yf, linestyle="None", marker="+")
        ax.set_xlim(float(fitrange_arr[0]), float(psimax))

    fileout = Path(str(path) + str(suffix))
    with fileout.open("w", encoding="utf-8") as fh:
        for xx, yy in zip(newx, newy):
            fh.write(f"{xx:12.6f}{yy:12.6f}\n")

    return params
