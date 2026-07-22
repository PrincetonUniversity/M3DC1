from __future__ import annotations

from pathlib import Path

import numpy as np

from .read_parameter import read_parameter


def radius_matrix(
    x,
    z,
    t=None,
    *,
    filename: str | Path = "C1.h5",
    rzero: float | None = None,
    cgs: bool = False,
    mks: bool = False,
) -> np.ndarray:
    """Port of radius_matrix.pro.

    Returns shape [nt, nx, nz] when t is sequence, else [nx, nz].
    For slab geometry (itor=0), the effective radius is the file's rzero.
    """
    xv = np.asarray(x, dtype=float).reshape(-1)
    zv = np.asarray(z, dtype=float).reshape(-1)
    nt = 1 if t is None else int(np.atleast_1d(t).size)

    itor = int(read_parameter("itor", filename=filename))

    if int(itor) == 0:
        if rzero is None:
            rzero = float(read_parameter("rzero", filename=filename, cgs=cgs, mks=mks))
        fill_value = float(rzero)
    else:
        fill_value = np.nan

    out = np.zeros((nt, xv.size, zv.size), dtype=float)
    for k in range(nt):
        if int(itor) == 0:
            out[k, :, :] = fill_value
        else:
            out[k, :, :] = xv[:, None]
    if nt == 1:
        return out[0]
    return out
