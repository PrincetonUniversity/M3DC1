from __future__ import annotations

import numpy as np


def dz(a: np.ndarray, z) -> np.ndarray:
    """Port of dz.pro derivative along z-axis (axis=1 for 2D field grids)."""
    arr = np.asarray(a)
    zv = np.asarray(z, dtype=float).reshape(-1)

    if arr.ndim == 2:
        return np.gradient(arr, zv, axis=1, edge_order=2)
    if arr.ndim == 3:
        out = np.empty_like(arr)
        for k in range(arr.shape[0]):
            out[k] = np.gradient(arr[k], zv, axis=1, edge_order=2)
        return out
    raise ValueError(f"Unsupported ndim for dz: {arr.ndim}")
