from __future__ import annotations

import os
import numpy as np

from .dx import dx
from .dz import dz

try:
    from numba import get_num_threads, njit, prange, set_num_threads
    from numba import config as numba_config

    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False
    get_num_threads = None
    prange = range
    set_num_threads = None
    numba_config = None

    def njit(*args, **kwargs):
        if args and callable(args[0]) and len(args) == 1 and not kwargs:
            return args[0]

        def _decorator(func):
            return func

        return _decorator


def _plane2(a: np.ndarray) -> np.ndarray:
    arr = np.asarray(a, dtype=float)
    if arr.ndim == 3:
        return np.asarray(arr[0], dtype=float)
    if arr.ndim == 2:
        return arr
    raise ValueError(f"Expected 2D/3D array, got shape {arr.shape}.")


@njit(cache=True)
def _interp2_numba(a: np.ndarray, x: float, y: float) -> float:
    nx, ny = a.shape
    if x < 0.0:
        x = 0.0
    elif x > nx - 1:
        x = nx - 1.0
    if y < 0.0:
        y = 0.0
    elif y > ny - 1:
        y = ny - 1.0

    x0 = int(np.floor(x))
    y0 = int(np.floor(y))
    x1 = x0 + 1
    y1 = y0 + 1
    if x1 >= nx:
        x1 = nx - 1
    if y1 >= ny:
        y1 = ny - 1

    tx = x - x0
    ty = y - y0

    v00 = a[x0, y0]
    v10 = a[x1, y0]
    v01 = a[x0, y1]
    v11 = a[x1, y1]
    return (1.0 - tx) * (1.0 - ty) * v00 + tx * (1.0 - ty) * v10 + (1.0 - tx) * ty * v01 + tx * ty * v11


@njit(cache=True, parallel=True)
def _create_map_kernel_numba(
    rf: np.ndarray,
    zf: np.ndarray,
    rf_x: np.ndarray,
    rf_z: np.ndarray,
    zf_x: np.ndarray,
    zf_z: np.ndarray,
    mask2: np.ndarray,
    r_out: np.ndarray,
    z_out: np.ndarray,
    tol: float,
    its: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    n = rf.shape[0]
    i0 = np.full((n, n), -1.0, dtype=np.float64)
    j0 = np.full((n, n), -1.0, dtype=np.float64)
    converged = 0

    for i in prange(n):
        for j in range(n):
            best = 1.0e300
            iguess = 0.0
            jguess = 0.0
            for ii in range(n):
                for jj in range(n):
                    d = (rf[ii, jj] - r_out[i]) ** 2 + (zf[ii, jj] - z_out[j]) ** 2
                    if d < best:
                        best = d
                        iguess = float(ii)
                        jguess = float(jj)

            ok = False
            for _ in range(its):
                rval = _interp2_numba(rf, iguess, jguess)
                zval = _interp2_numba(zf, iguess, jguess)
                d2 = (rval - r_out[i]) ** 2 + (zval - z_out[j]) ** 2
                if d2 <= tol * tol:
                    i0[i, j] = iguess
                    j0[i, j] = jguess
                    converged += 1
                    ok = True
                    break

                drdi = _interp2_numba(rf_x, iguess, jguess)
                drdj = _interp2_numba(rf_z, iguess, jguess)
                dzdi = _interp2_numba(zf_x, iguess, jguess)
                dzdj = _interp2_numba(zf_z, iguess, jguess)

                dd2di = 2.0 * ((rval - r_out[i]) * drdi + (zval - z_out[j]) * dzdi)
                dd2dj = 2.0 * ((rval - r_out[i]) * drdj + (zval - z_out[j]) * dzdj)
                denom = dd2di * dd2di + dd2dj * dd2dj
                if denom <= np.finfo(np.float64).tiny:
                    break

                dl = -2.0 * d2 / denom
                di = dl * dd2di
                dj = dl * dd2dj

                if iguess + di < 0.0:
                    di = -iguess * 0.5
                elif iguess + di > n - 1:
                    di = (n - 1 - iguess) * 0.5

                if jguess + dj < 0.0:
                    dj = -jguess * 0.5
                elif jguess + dj > n - 1:
                    dj = (n - 1 - jguess) * 0.5

                for _ in range(11):
                    m = _interp2_numba(mask2, iguess + di, jguess + dj)
                    if abs(m) < 0.02:
                        break
                    di *= 0.5
                    dj *= 0.5

                iguess += di
                jguess += dj

            if not ok:
                i0[i, j] = -1.0
                j0[i, j] = -1.0

    for i in range(1, n - 1):
        for j in range(1, n - 1):
            if i0[i, j] != -1.0:
                continue
            interp = False
            if i0[i - 1, j] != -1.0 and i0[i + 1, j] != -1.0:
                i0[i, j] = 0.5 * (i0[i - 1, j] + i0[i + 1, j])
                j0[i, j] = 0.5 * (j0[i - 1, j] + j0[i + 1, j])
                interp = True
            if i0[i, j - 1] != -1.0 and i0[i, j + 1] != -1.0:
                if not interp:
                    i0[i, j] = 0.5 * (i0[i, j - 1] + i0[i, j + 1])
                    j0[i, j] = 0.5 * (j0[i, j - 1] + j0[i, j + 1])
                else:
                    i0[i, j] = 0.25 * (i0[i, j - 1] + i0[i, j + 1] + i0[i - 1, j] + i0[i + 1, j])
                    j0[i, j] = 0.25 * (j0[i, j - 1] + j0[i, j + 1] + j0[i - 1, j] + j0[i + 1, j])

    return i0, j0, converged


def create_map(
    rfield: np.ndarray,
    zfield: np.ndarray,
    *,
    r=None,
    z=None,
    ix=None,
    iy=None,
    mask=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Port of create_map.pro.

    Uses Numba JIT when available; otherwise the same decorated code runs in Python mode.
    Returns (ix, iy, r, z, mask).
    """
    rf = _plane2(rfield)
    zf = _plane2(zfield)
    if rf.shape != zf.shape:
        raise ValueError(f"rfield/zfield shape mismatch: {rf.shape} vs {zf.shape}")

    n0, n1 = rf.shape
    if n0 != n1:
        raise ValueError(f"create_map expects square grid, got {rf.shape}")
    n = n0

    if mask is None:
        mask2 = np.zeros((n, n), dtype=float)
    else:
        mask2 = _plane2(mask)
        if mask2.shape != (n, n):
            raise ValueError(f"mask shape mismatch: {mask2.shape} vs {(n, n)}")

    grid = np.arange(n, dtype=float)
    rf_x = dx(rf, grid)
    rf_z = dz(rf, grid)
    zf_x = dx(zf, grid)
    zf_z = dz(zf, grid)

    rmin = float(np.min(rf))
    rmax = float(np.max(rf))
    zmin = float(np.min(zf))
    zmax = float(np.max(zf))

    # IDL create_map always rebuilds r/z as a uniform grid from min/max.
    del r, z, ix, iy
    r_out = np.linspace(rmin, rmax, n, dtype=float)
    z_out = np.linspace(zmin, zmax, n, dtype=float)

    tol = 1.0e-4
    its = 100

    if _NUMBA_AVAILABLE and set_num_threads is not None:
        workers = 1
        w_env = os.getenv("M3DC1_NUM_WORKERS", "").strip()
        if w_env:
            try:
                workers = int(w_env)
            except ValueError:
                workers = 1
        if workers < 1:
            workers = 1
        max_threads = int(getattr(numba_config, "NUMBA_NUM_THREADS", workers))
        workers = max(1, min(int(workers), max_threads))
        set_num_threads(int(workers))

    i0, j0, converged = _create_map_kernel_numba(
        rf,
        zf,
        rf_x,
        rf_z,
        zf_x,
        zf_z,
        mask2,
        r_out,
        z_out,
        tol,
        its,
    )

    print(f"converged = {converged}")
    return i0, j0, np.asarray(r_out), np.asarray(z_out), mask2
