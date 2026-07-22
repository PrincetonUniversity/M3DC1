from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from matplotlib.figure import Figure

from .lcfs import lcfs
from .read_field import read_field


@dataclass
class ShapeContour:
    level: float
    r: np.ndarray
    z: np.ndarray


@dataclass
class ShapeResult:
    contours: list[ShapeContour]
    levels: np.ndarray
    filename: str
    symbol: str = ""
    units: str = ""


def _contours_from_field(psi: np.ndarray, r: np.ndarray, z: np.ndarray, levels: np.ndarray) -> list[ShapeContour]:
    fig = Figure()
    ax = fig.add_subplot(111)
    try:
        cs = ax.contour(r, z, psi.T, levels=levels)
        contours: list[ShapeContour] = []
        for level, segs in zip(np.asarray(cs.levels, dtype=float), cs.allsegs):
            for seg in segs:
                path = np.asarray(seg, dtype=float)
                if path.ndim != 2 or path.shape[0] == 0 or path.shape[1] < 2:
                    continue
                contours.append(
                    ShapeContour(
                        level=float(level),
                        r=np.asarray(path[:, 0], dtype=float),
                        z=np.asarray(path[:, 1], dtype=float),
                    )
                )
        return contours
    finally:
        fig.clear()


def read_shape(
    timeslice: int = 0,
    *,
    filename: str | Path | Sequence[str | Path] = "C1.h5",
    points: int = 200,
    xrange=None,
    yrange=None,
    logical: bool = False,
    phi: float = 0.0,
    cgs: bool = False,
    mks: bool = True,
    levels=None,
    **kwargs,
):
    """
    Read flux-surface shape contours using the same level logic as plot_shape().
    """
    flist = [filename] if isinstance(filename, (str, Path)) else list(filename)
    results: list[ShapeResult] = []

    for fn in flist:
        meta = read_field(
            "psi",
            timeslices=timeslice,
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
        r = np.asarray(meta.r, dtype=float).reshape(-1)
        z = np.asarray(meta.z, dtype=float).reshape(-1)

        if levels is None:
            lc = lcfs(psi2d, r, z, filename=fn, slice=timeslice, cgs=cgs, mks=mks)
            contour_levels = np.arange(22, dtype=float) / 20.0 * (lc.psilim - lc.flux0) + lc.flux0
            if lc.psilim < lc.flux0:
                contour_levels = contour_levels[::-1]
        else:
            contour_levels = np.asarray(levels, dtype=float).reshape(-1)

        results.append(
            ShapeResult(
                contours=_contours_from_field(psi2d, r, z, contour_levels),
                levels=np.asarray(contour_levels, dtype=float),
                filename=str(fn),
                symbol=str(meta.symbol),
                units=str(meta.units),
            )
        )

    if isinstance(filename, (str, Path)):
        return results[0]
    return results
