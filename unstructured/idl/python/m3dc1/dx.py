from __future__ import annotations

import numpy as np


def dx(a: np.ndarray, x) -> np.ndarray:
    """Port of dx.pro derivative along x-axis (axis=0 for 2D field grids)."""
    arr = np.asarray(a)
    xv = np.asarray(x, dtype=float).reshape(-1)

    if arr.ndim == 2:
        return np.gradient(arr, xv, axis=0, edge_order=2)
    if arr.ndim == 3:
        out = np.empty_like(arr)
        for k in range(arr.shape[0]):
            out[k] = np.gradient(arr[k], xv, axis=0, edge_order=2)
        return out
    raise ValueError(f"Unsupported ndim for dx: {arr.ndim}")
