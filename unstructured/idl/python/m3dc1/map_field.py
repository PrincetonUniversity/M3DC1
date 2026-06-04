from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def _plane2(a: np.ndarray) -> np.ndarray:
    arr = np.asarray(a, dtype=float)
    if arr.ndim == 3:
        return np.asarray(arr[0], dtype=float)
    if arr.ndim == 2:
        return arr
    raise ValueError(f"Expected 2D/3D array, got shape {arr.shape}.")


def _map_field_kernel_py(
    f: np.ndarray,
    mx: np.ndarray,
    my: np.ndarray,
    val: float,
) -> tuple[np.ndarray, np.ndarray]:
    n0, n1 = f.shape
    out_mask = mx == -1.0
    interp = RegularGridInterpolator(
        (np.arange(n0, dtype=float), np.arange(n1, dtype=float)),
        f,
        bounds_error=False,
        fill_value=val,
    )
    pts = np.column_stack((mx.ravel(), my.ravel()))
    result = np.asarray(interp(pts), dtype=float).reshape(n0, n1)
    result[out_mask] = val
    return result, out_mask


def map_field(
    field: np.ndarray,
    ix: np.ndarray,
    iy: np.ndarray,
    *,
    mask=None,
    outval=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Port of map_field.pro.

    Returns (mapped_field, mask) where mask marks unmapped points.
    """
    f = _plane2(field)
    mx = _plane2(ix)
    my = _plane2(iy)
    if f.shape != mx.shape or f.shape != my.shape:
        raise ValueError(f"Shape mismatch field={f.shape}, ix={mx.shape}, iy={my.shape}")

    val = 0.0 if outval is None else float(outval)

    result, out_mask = _map_field_kernel_py(np.asarray(f, dtype=float), np.asarray(mx, dtype=float), np.asarray(my, dtype=float), val)

    if mask is not None:
        out_mask = np.logical_or(out_mask, _plane2(mask).astype(bool))

    return result, out_mask
