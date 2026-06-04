from __future__ import annotations

import numpy as np


def radius_matrix(x, z, t=None) -> np.ndarray:
    """Port of radius_matrix.pro.

    Returns shape [nt, nx, nz] when t is sequence, else [nx, nz].
    """
    xv = np.asarray(x, dtype=float).reshape(-1)
    zv = np.asarray(z, dtype=float).reshape(-1)
    nt = 1 if t is None else int(np.atleast_1d(t).size)

    out = np.zeros((nt, xv.size, zv.size), dtype=float)
    for k in range(nt):
        out[k, :, :] = xv[:, None]
    if nt == 1:
        return out[0]
    return out
