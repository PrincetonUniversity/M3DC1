from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from .get_colors import get_colors
from .lcfs import lcfs
from .read_field import read_field


def plot_shape(
    filename: str | Path | Sequence[str | Path] = "C1.h5",
    *,
    names: Sequence[str] | None = None,
    points: int = 200,
    xrange=None,
    yrange=None,
    overplot: bool = False,
    logical: bool = False,
    phi: float = 0.0,
    cgs: bool = False,
    mks: bool = True,
    iso: bool = True,
    **kwargs,
):
    """
    Python port of plot_shape.pro.
    """
    if "xrange" in kwargs or "yrange" in kwargs:
        raise TypeError("plot_shape() uses 'xrange'/'yrange' arguments, not 'xrange'/'yrange'.")

    flist = [filename] if isinstance(filename, (str, Path)) else list(filename)
    colors = get_colors(max(len(flist), 10))

    if overplot:
        ax = plt.gca()
        fig = ax.figure
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlabel("R (m)" if mks else "R")
        ax.set_ylabel("Z (m)" if mks else "Z")

    handles = []
    legend_names = None if names is None else list(names)
    for i, fn in enumerate(flist):
        meta = read_field(
            "psi",
            timeslices=-1,
            points=points,
            xrange=xrange,
            yrange=yrange,
            phi=phi,
            logical=logical,
            cgs=cgs,
            mks=mks,
            filename=fn,
            return_meta=True,
            **kwargs,
        )
        psi = np.asarray(meta.data)
        psi2d = psi[0, :, :] if psi.ndim == 3 else psi
        xv = np.asarray(meta.r, dtype=float).reshape(-1)
        zv = np.asarray(meta.z, dtype=float).reshape(-1)

        lc = lcfs(psi2d, xv, zv, filename=fn, cgs=cgs, mks=mks)
        levels = np.arange(12, dtype=float) / 10.0 * (lc.psilim - lc.flux0) + lc.flux0
        if lc.psilim < lc.flux0:
            levels = levels[::-1]

        ax.contour(xv, zv, psi2d.T, levels=levels, colors=[colors[i]], linewidths=1.0)

        if legend_names is not None and i < len(legend_names):
            handles.append(Line2D([0], [0], color=colors[i], linewidth=1.0, label=str(legend_names[i])))

    if handles:
        ax.legend(handles=handles, loc="best")
    if xrange is not None:
        ax.set_xlim(xrange)
    if yrange is not None:
        ax.set_ylim(yrange)
    if iso:
        ax.set_aspect("equal", adjustable="box")
    return fig, ax
