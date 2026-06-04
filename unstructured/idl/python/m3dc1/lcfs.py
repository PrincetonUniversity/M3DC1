from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .read_lcfs import read_lcfs


@dataclass
class LcfsInfo:
    psilim: float
    flux0: float
    axis: np.ndarray
    xpoint: np.ndarray


def lcfs(psi, x, z, *, axis=None, xpoint=None, flux0=None, filename="C1.h5", slice: int | None = None, **kwargs) -> LcfsInfo:
    """
    Lightweight lcfs.pro-style helper.
    """
    p = np.asarray(psi)
    if p.ndim == 3:
        p = p[0, :, :]
    xv = np.asarray(x, dtype=float).reshape(-1)
    zv = np.asarray(z, dtype=float).reshape(-1)

    try:
        meta = read_lcfs(filename=filename, slice=slice, return_meta=True, **kwargs)
        ax = np.asarray(meta.axis, dtype=float)
        xp = np.asarray(meta.xpoint, dtype=float)
        f0 = float(meta.flux0)
        ps = float(meta.psilim)
    except Exception:
        idx = np.unravel_index(int(np.argmin(p)), p.shape)
        ax = np.asarray([xv[idx[0]], zv[idx[1]]], dtype=float)
        xp = np.asarray(ax, dtype=float)
        f0 = float(np.nanmin(p))
        ps = float(np.nanmax(p))

    return LcfsInfo(psilim=ps, flux0=f0, axis=ax, xpoint=xp)
