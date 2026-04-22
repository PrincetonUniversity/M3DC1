from __future__ import annotations

import builtins
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .plot_legend import plot_legend
from .read_hmn import read_hmn


def _smooth_1d(x: np.ndarray, n: int) -> np.ndarray:
    if n <= 1:
        return np.asarray(x, dtype=float)
    k = np.ones(int(n), dtype=float) / float(n)
    return np.convolve(np.asarray(x, dtype=float), k, mode="same")


def plot_hmn(
    *,
    filename: str | Path = "C1.h5",
    maxn: int | None = None,
    growth: bool = False,
    outfile: str | Path | None = None,
    xrange=None,
    yrange=None,
    smooth: int | None = None,
    overplot: bool = False,
    thick: float | None = None,
    linestyle="-",
    ke: bool = False,
    me: bool = False,
    cgs: bool = False,
    mks: bool = False,
    xscale: float = 1.0,
    labelx: float = 0.9,
    nolegend: bool = False,
    print: int | None = None,
    **kwargs,
):
    """Python port of plot_hmn.pro."""
    del kwargs, ke

    meta = read_hmn(
        filename=filename,
        maxn=maxn,
        growth=growth,
        outfile=outfile,
        me=me,
        cgs=cgs,
        mks=mks,
        return_meta=True,
    )
    data = np.asarray(meta.data, dtype=float)
    time = np.asarray(meta.time, dtype=float)
    if data.size == 0 or time.size == 0:
        return None, None

    if smooth is not None:
        for n in range(data.shape[0]):
            data[n, :] = _smooth_1d(data[n, :], int(smooth))

    if print is not None:
        nprint = int(print)
        if nprint < 0 or nprint >= data.shape[0]:
            raise ValueError(f"plot_hmn(print={nprint}) is out of range for 0 <= n < {data.shape[0]}.")
        yvals = np.asarray(data[nprint, :], dtype=float)
        for i in range(0, yvals.size, 8):
            chunk = yvals[i : i + 8]
            builtins.print(" ".join(f"{value:12.6g}" for value in chunk))

    if yrange is None:
        yrange_data = data
        if meta.magnetic and data.shape[0] > 1:
            yrange_data = data[1:, :]
        ymin = float(np.nanmin(yrange_data))
        ymax = float(np.nanmax(yrange_data))
        if ymax > ymin:
            ymax = ymax + 0.1 * (ymax - ymin)
        else:
            ymax = ymax + 0.1 * max(abs(ymax), 1.0)
        yrange = [ymin, ymax]

    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
    colors = [cycle[i % len(cycle)] for i in range(data.shape[0])]

    if not overplot:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.set_xlabel(meta.xtitle)
        ax.set_ylabel(meta.ytitle)
    else:
        ax = plt.gca()
        fig = ax.figure

    lw = float(thick) if thick is not None else None
    tscaled = time * float(xscale)
    ntimes = int(data.shape[1])
    for n in range(data.shape[0]):
        c = colors[n % len(colors)] if colors else None
        ax.plot(tscaled, data[n, :], linestyle=linestyle, linewidth=lw, color=c)
        m = min(int(ntimes * float(labelx)), ntimes - 1)
        ax.text(tscaled[m], data[n, m], f"{n}", color=c if c is not None else "black")

    if xrange is not None:
        ax.set_xlim(xrange)
    else:
        xmax = ax.get_xlim()[1]
        ax.set_xlim(0.0, xmax)
    if np.isfinite(yrange[0]) and np.isfinite(yrange[1]) and yrange[0] != yrange[1]:
        ax.set_ylim(float(yrange[0]), float(yrange[1]))
    if growth:
        ax.axhline(0.0, color="black", linestyle="-", linewidth=0.8)

    if not nolegend:
        names = [f"n={i}" for i in range(data.shape[0])]
        plot_legend(names, colors=colors, linestyles=[linestyle] * data.shape[0])

    return fig, ax
