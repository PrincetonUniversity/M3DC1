from __future__ import annotations

import numpy as np

from .dx import dx
from .dz import dz


def delstar(a: np.ndarray, x, z, toroidal: bool = False) -> np.ndarray:
    """Port of delstar.pro."""
    arr = np.asarray(a)
    xvec = np.asarray(x, dtype=float).reshape(-1)

    out = dx(dx(arr, xvec), xvec) + dz(dz(arr, z), z)
    if toroidal:
        if out.ndim == 2:
            out = out - dx(arr, xvec) / np.maximum(xvec[:, None], np.finfo(float).tiny)
        else:
            out = out - dx(arr, xvec) / np.maximum(xvec[None, :, None], np.finfo(float).tiny)
    return out
