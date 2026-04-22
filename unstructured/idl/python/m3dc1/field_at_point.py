from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def _interp2(field2d: np.ndarray, x: np.ndarray, z: np.ndarray, x0: np.ndarray, z0: np.ndarray) -> np.ndarray:
    if x.size < 2 or z.size < 2:
        return np.full_like(x0, np.nan, dtype=float)
    interp = RegularGridInterpolator(
        (np.asarray(x, dtype=float), np.asarray(z, dtype=float)),
        np.asarray(field2d, dtype=float),
        method="nearest",
        bounds_error=False,
        fill_value=np.nan,
    )
    pts = np.column_stack([np.asarray(x0, dtype=float).reshape(-1), np.asarray(z0, dtype=float).reshape(-1)])
    return np.asarray(interp(pts), dtype=float).reshape(np.asarray(x0).shape)


def field_at_point(field, x, z, x0, z0):
    """
    IDL-like field_at_point interpolation on regular (x,z) grid.
    """
    f = np.asarray(field)
    xv = np.asarray(x, dtype=float).reshape(-1)
    zv = np.asarray(z, dtype=float).reshape(-1)
    xq, zq = np.broadcast_arrays(np.asarray(x0, dtype=float), np.asarray(z0, dtype=float))
    out_shape = xq.shape
    xqf = xq.reshape(-1)
    zqf = zq.reshape(-1)

    if f.ndim == 2:
        return _interp2(f, xv, zv, xqf, zqf).reshape(out_shape)
    if f.ndim == 3:
        out = np.zeros((f.shape[0], xqf.size), dtype=float)
        for k in range(f.shape[0]):
            out[k, :] = _interp2(np.asarray(f[k, :, :], dtype=float), xv, zv, xqf, zqf)
        return out.reshape((f.shape[0],) + out_shape)
    raise ValueError(f"field_at_point expects 2D or 3D field array, got shape {f.shape}.")
