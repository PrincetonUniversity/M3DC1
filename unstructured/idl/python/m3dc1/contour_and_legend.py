from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


def _choose_cmap(panel: np.ndarray) -> str:
    data = np.asarray(panel, dtype=float)
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return "turbo"
    dmin = float(np.min(finite))
    dmax = float(np.max(finite))
    if dmin >= 0.0 or dmax <= 0.0:
        return "plasma"
    return "turbo"


def contour_and_legend_single(
    ax,
    panel: np.ndarray,
    xv: np.ndarray,
    yv: np.ndarray,
    *,
    contour_levels,
    vmin,
    vmax,
    fill: bool,
    colorbar: bool,
    lines: bool,
    label: str,
    title: str,
    xtitle: str | None,
    ytitle: str | None,
    cmap: str | None = None,
    colorbar_figure=None,
):
    data = np.asarray(panel, dtype=float)
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        maxval, minval, fracdiff = 0.0, 0.0, 0.0
    else:
        maxval = float(np.max(finite))
        minval = float(np.min(finite))
        scale = max(abs(maxval), abs(minval), 1e-30)
        fracdiff = abs(maxval - minval) / scale
    print("maxval, minval, fracdiff = ", maxval, minval, fracdiff)

    cmap_name = cmap if cmap else _choose_cmap(panel)
    cmap_obj = cm.get_cmap(cmap_name).copy()
    cmap_obj.set_under("black")
    cmap_obj.set_over("white")

    lev = np.asarray(contour_levels)
    if lev.ndim == 0:
        level_min = minval
        level_max = maxval
        nlevels = int(lev)
    else:
        level_min = float(np.min(lev))
        level_max = float(np.max(lev))
        nlevels = lev.size
    full_span = level_max - level_min
    if full_span > 0.0:
        level_pad = 0.01 * full_span
    else:
        level_pad = 0.01 * max(abs(level_min), abs(level_max), 1.0)
    contour_levels_use = np.linspace(level_min - level_pad, level_max + level_pad, nlevels)

    cf = None
    if fill:
        cf = ax.contourf(
            xv,
            yv,
            panel.T,
            levels=contour_levels_use,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap_obj,
            extend="both",
        )
    if lines or not fill:
        ax.contour(xv, yv, panel.T, levels=contour_levels_use, colors="k", linewidths=0.5)

    if fill and colorbar and cf is not None:
        if colorbar_figure is None:
            if label:
                plt.colorbar(cf, ax=ax, label=label)
            else:
                plt.colorbar(cf, ax=ax)
        else:
            colorbar_figure.colorbar(cf, ax=ax, label=label)

    if title:
        ax.set_title(str(title))
    if xtitle:
        ax.set_xlabel(str(xtitle))
    if ytitle:
        ax.set_ylabel(str(ytitle))


def contour_and_legend(
    z,
    x,
    y,
    *,
    label="",
    range=None,
    levels=None,
    title="",
    fill: bool = True,
    colorbar: bool = True,
    lines=False,
    overplot=False,
    cmap: str | None = None,
    xtitle: str | None = None,
    ytitle: str | None = None,
    jpeg: str | Path | None = None,
    **kwargs,
):
    """
    Lightweight Python port of contour_and_legend.pro.
    """
    zz = np.asarray(z)
    xv = np.asarray(x, dtype=float).reshape(-1)
    yv = np.asarray(y, dtype=float).reshape(-1)

    if zz.ndim == 2:
        panels = zz[None, ...]
    elif zz.ndim == 3:
        panels = zz
    else:
        raise ValueError(f"Expected 2D/3D input z, got shape {zz.shape}.")

    n = panels.shape[0]
    labels = [label] * n if isinstance(label, str) else list(label)
    titles = [title] * n if isinstance(title, str) else list(title)
    while len(labels) < n:
        labels.append("")
    while len(titles) < n:
        titles.append("")

    vmin = vmax = None
    if range is not None:
        rr = np.asarray(range, dtype=float).reshape(-1)
        if rr.size >= 2:
            vmin, vmax = float(rr[0]), float(rr[1])
    if vmin is not None and vmax is not None:
        if levels is None:
            contour_levels = np.linspace(vmin, vmax, 100)
        else:
            lev = np.asarray(levels)
            if lev.ndim == 0:
                contour_levels = np.linspace(vmin, vmax, int(lev))
            else:
                contour_levels = np.linspace(vmin, vmax, lev.size)
    else:
        contour_levels = 100 if levels is None else levels

    if n == 1:
        ax = plt.gca() if overplot else plt.figure(figsize=(7, 5)).gca()
        contour_and_legend_single(
            ax,
            panels[0],
            xv,
            yv,
            contour_levels=contour_levels,
            vmin=vmin,
            vmax=vmax,
            fill=fill,
            colorbar=colorbar,
            lines=lines,
            label=labels[0],
            title=titles[0],
            xtitle=xtitle,
            ytitle=ytitle,
            cmap=cmap,
        )
    else:
        rows = int(np.ceil(np.sqrt(n)))
        cols = int(np.ceil(n / rows))
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
        axs = np.asarray(axes).reshape(-1)
        for i in range(rows * cols):
            ax = axs[i]
            if i >= n:
                ax.axis("off")
                continue
            contour_and_legend_single(
                ax,
                panels[i],
                xv,
                yv,
                contour_levels=contour_levels,
                vmin=vmin,
                vmax=vmax,
                fill=fill,
                colorbar=colorbar,
                lines=lines,
                label=labels[i],
                title=titles[i],
                xtitle=xtitle,
                ytitle=ytitle,
                cmap=cmap,
                colorbar_figure=fig,
            )
        fig.tight_layout()

    if jpeg is not None:
        plt.savefig(str(jpeg), dpi=150)
