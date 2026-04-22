from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .read_parameter import read_parameter
from .read_mesh import MeshData

try:
    from numba import njit

    _NUMBA_AVAILABLE = True
except ImportError:
    def njit(*args, **kwargs):
        if args and callable(args[0]) and len(args) == 1 and not kwargs:
            return args[0]

        def _decorator(func):
            return func

        return _decorator

    _NUMBA_AVAILABLE = False

_MI = np.array([0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 3, 2, 1, 0], dtype=np.int64)
_NI = np.array([0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 2, 3, 4, 5], dtype=np.int64)


@dataclass
class EvalFieldResult:
    data: np.ndarray
    r: np.ndarray
    z: np.ndarray
    mask: np.ndarray
    edge_val: float


@njit(cache=True)
def _is_in_tri_numba(lp0: float, lp1: float, a: float, b: float, c: float) -> bool:
    small = (a + b + c) * 1e-4
    if lp1 < -small:
        return False
    if lp1 > c + small:
        return False
    x = 1.0 - lp1 / c
    if lp0 < -b * x - small:
        return False
    if lp0 > a * x + small:
        return False
    return True

@njit(cache=True)
def _eval_poly_numba(
    f: np.ndarray,
    lp0: float,
    lp1: float,
    lp2: float,
    theta: float,
    elm: int,
    operation: int,
) -> float:
    px = np.empty(6, dtype=np.float64)
    py = np.empty(6, dtype=np.float64)
    px[0] = 1.0
    py[0] = 1.0
    for i in range(1, 6):
        px[i] = px[i - 1] * lp0
        py[i] = py[i - 1] * lp1

    co = np.cos(theta)
    sn = np.sin(theta)
    op = operation
    op2 = (op - 1) // 10
    op1 = op - op2 * 10

    threed = f.shape[0] == 80
    out = 0.0

    for k in range(20):
        mi = _MI[k]
        ni = _NI[k]
        mi1 = mi - 1 if mi > 0 else 0
        mi2 = mi - 2 if mi > 1 else 0
        mi3 = mi - 3 if mi > 2 else 0
        ni1 = ni - 1 if ni > 0 else 0
        ni2 = ni - 2 if ni > 1 else 0
        ni3 = ni - 3 if ni > 2 else 0

        if op1 == 1:
            temp = px[mi] * py[ni]
        elif op1 == 2:
            temp = co * mi * px[mi1] * py[ni] - sn * ni * px[mi] * py[ni1]
        elif op1 == 3:
            temp = sn * mi * px[mi1] * py[ni] + co * ni * px[mi] * py[ni1]
        elif op1 == 4:
            temp = (
                co * co * mi * mi1 * px[mi2] * py[ni]
                + sn * sn * ni * ni1 * px[mi] * py[ni2]
                - 2.0 * co * sn * ni * mi * px[mi1] * py[ni1]
            )
        elif op1 == 5:
            temp = (
                co * sn * mi * mi1 * px[mi2] * py[ni]
                - co * sn * ni * ni1 * px[mi] * py[ni2]
                + (co * co - sn * sn) * ni * mi * px[mi1] * py[ni1]
            )
        elif op1 == 6:
            temp = (
                sn * sn * mi * mi1 * px[mi2] * py[ni]
                + co * co * ni * ni1 * px[mi] * py[ni2]
                + 2.0 * co * sn * ni * mi * px[mi1] * py[ni1]
            )
        elif op1 == 7:
            temp = mi * mi1 * px[mi2] * py[ni] + ni * ni1 * py[ni2] * px[mi]
        elif op1 == 8:
            temp = (
                co * mi * mi1 * mi2 * px[mi3] * py[ni]
                - sn * ni * ni1 * ni2 * py[ni3] * px[mi]
                - sn * mi * mi1 * px[mi2] * ni * py[ni1]
                + co * mi * px[mi1] * ni * ni1 * py[ni2]
            )
        elif op1 == 9:
            temp = (
                sn * mi * mi1 * mi2 * px[mi3] * py[ni]
                + co * ni * ni1 * ni2 * py[ni3] * px[mi]
                + co * mi * mi1 * px[mi2] * ni * py[ni1]
                + sn * mi * px[mi1] * ni * ni1 * py[ni2]
            )
        else:
            temp = px[mi] * py[ni]

        if op2 == 0:
            out += f[k, elm] * temp
            if threed:
                out += temp * (f[k + 20, elm] * lp2 + f[k + 40, elm] * lp2 * lp2 + f[k + 60, elm] * lp2 * lp2 * lp2)
        elif op2 == 1:
            if threed:
                out += temp * (f[k + 20, elm] + f[k + 40, elm] * lp2 * 2.0 + f[k + 60, elm] * lp2 * lp2 * 3.0)
        elif op2 == 2:
            if threed:
                out += temp * (f[k + 40, elm] * 2.0 + f[k + 60, elm] * lp2 * 6.0)

    return out

@njit(cache=True)
def _eval_field_kernel_numba(
    f: np.ndarray,
    elm_data: np.ndarray,
    xi: np.ndarray,
    yi: np.ndarray,
    result: np.ndarray,
    out_mask: np.ndarray,
    nelms: int,
    operation: int,
    p: int,
    xrange0: float,
    yrange0: float,
    xzero: float,
    zzero: float,
    dx: float,
    dy: float,
    threed: bool,
    ib: int,
    local_phi: float,
    wall_mask: bool,
    version: float,
    debug_spin: int,
) -> None:
    small = 1e-3
    nrows = elm_data.shape[0]
    for i in range(nelms):
        izone = 1.0
        localphi_i = 0.0

        if threed:
            if ib + 2 >= nrows:
                continue
            dphi = elm_data[ib + 1, i]
            phi_elm = elm_data[ib + 2, i]
            localphi_i = local_phi - phi_elm
            if localphi_i < 0.0 or localphi_i > dphi:
                continue

        a = elm_data[0, i]
        b = elm_data[1, i]
        c = elm_data[2, i]
        t = elm_data[3, i]
        x0 = elm_data[4, i]
        y0 = elm_data[5, i]
        if version >= 16.0 and ib < nrows:
            izone = elm_data[ib, i]

        co = np.cos(t)
        sn = np.sin(t)

        p1x = x0
        p1y = y0
        p2x = p1x + (b + a) * co
        p2y = p1y + (b + a) * sn
        p3x = p1x + (b * co - c * sn)
        p3y = p1y + (b * sn + c * co)

        minx = min(p1x, p2x, p3x)
        maxx = max(p1x, p2x, p3x)
        miny = min(p1y, p2y, p3y)
        maxy = max(p1y, p2y, p3y)

        if maxx < xrange0 - xzero or maxy < yrange0 - zzero:
            continue

        if p == 1:
            i0 = 0
            i1 = 0
            j0 = 0
            j1 = 0
        else:
            i0 = max(0, int(np.floor((minx - xrange0 + xzero) / dx)))
            i1 = min(p - 1, int(np.ceil((maxx - xrange0 + xzero) / dx + small)))
            j0 = max(0, int(np.floor((miny - yrange0 + zzero) / dy)))
            j1 = min(p - 1, int(np.ceil((maxy - yrange0 + zzero) / dy + small)))

        for j in range(j0, j1 + 1):
            py = yi[j] - zzero
            for ii in range(i0, i1 + 1):
                px = xi[ii] - xzero
                lp0 = (px - x0) * co + (py - y0) * sn - b
                lp1 = -(px - x0) * sn + (py - y0) * co
                lp2 = localphi_i
                if _is_in_tri_numba(lp0, lp1, a, b, c):
                    burn = 0.0
                    if debug_spin > 0:
                        for kk in range(debug_spin):
                            burn += _eval_poly_numba(f, lp0, lp1, lp2, t, i, operation)
                    if wall_mask and izone != 2.0:
                        result[ii, j] = 0.0
                    else:
                        value = _eval_poly_numba(f, lp0, lp1, lp2, t, i, operation)
                        if debug_spin > 0:
                            value = value + burn * 0.0
                        result[ii, j] = value
                        out_mask[ii, j] = 0.0

def eval_field(
    field: np.ndarray,
    mesh: MeshData,
    r=None,
    z=None,
    points: int = 100,
    operation: int = 1,
    filename: str | Path = "C1.h5",
    xrange=None,
    yrange=None,
    mask=None,
    phi: float = 0.0,
    wall_mask: bool = False,
    map_r=None,
    map_z=None,
    edge_val: float | None = None,
    debug_spin: int = 0,
) -> EvalFieldResult:
    """Evaluate per-element polynomial field on a regular grid."""
    p = int(points)
    if p <= 0:
        raise ValueError("points must be > 0")

    f = np.asarray(field, dtype=float)
    if f.ndim != 2:
        raise ValueError(f"Expected field shape [n_coeff, n_elem], got {f.shape}")
    if f.shape[0] > f.shape[1]:
        f = f.T
    f = np.ascontiguousarray(f, dtype=np.float64)

    elm_data = np.asarray(mesh.elements, dtype=float)
    if elm_data.ndim != 2:
        raise ValueError(f"Unexpected mesh elements shape: {elm_data.shape}")
    if elm_data.shape[0] > elm_data.shape[1]:
        elm_data = elm_data.T
    elm_data = np.ascontiguousarray(elm_data, dtype=np.float64)

    nelms = int(mesh.nelms)
    if nelms <= 0:
        nelms = int(elm_data.shape[1])

    version = float(read_parameter("version", filename=filename))

    if version == 0:
        xzero = float(read_parameter("xzero", filename=filename))
        zzero = float(read_parameter("zzero", filename=filename))
    else:
        xzero = 0.0
        zzero = 0.0

    period = float(mesh.period)
    local_phi = float(phi) - np.floor(float(phi) / period) * period

    threed = elm_data.shape[0] > 8

    if xrange is None or yrange is None or len(np.atleast_1d(xrange)) < 2 or len(np.atleast_1d(yrange)) < 2:
        minx = float(np.min(elm_data[4, :]))
        maxx = float(np.max(elm_data[4, :]))
        miny = float(np.min(elm_data[5, :]))
        maxy = float(np.max(elm_data[5, :]))

        ncheck = nelms if not threed else min(nelms, nelms // max(mesh.nplanes, 1))
        for i in range(ncheck):
            i_data = elm_data[:, i]
            a, b, c, t, x0, y0 = i_data[0:6]
            co = np.cos(t)
            sn = np.sin(t)
            p1 = np.array([x0, y0], dtype=float)
            p2 = p1 + np.array([(b + a) * co, (b + a) * sn], dtype=float)
            p3 = p1 + np.array([b * co - c * sn, b * sn + c * co], dtype=float)

            minx = min(minx, p2[0], p3[0])
            maxx = max(maxx, p2[0], p3[0])
            miny = min(miny, p2[1], p3[1])
            maxy = max(maxy, p2[1], p3[1])

        if xrange is None or len(np.atleast_1d(xrange)) < 2:
            xrange = np.array([minx + xzero, maxx + xzero], dtype=float)
        else:
            xrange = np.asarray(xrange, dtype=float)
        if yrange is None or len(np.atleast_1d(yrange)) < 2:
            yrange = np.array([miny + zzero, maxy + zzero], dtype=float)
        else:
            yrange = np.asarray(yrange, dtype=float)
    else:
        xrange = np.asarray(xrange, dtype=float)
        yrange = np.asarray(yrange, dtype=float)

    if p == 1:
        xi = np.array([xrange[0]], dtype=float)
        yi = np.array([yrange[0]], dtype=float)
        dx = 0.0
        dy = 0.0
    else:
        xi = np.linspace(xrange[0], xrange[1], p)
        yi = np.linspace(yrange[0], yrange[1], p)
        dx = float(xi[1] - xi[0])
        dy = float(yi[1] - yi[0])

    result = np.zeros((p, p), dtype=np.float64)
    out_mask = np.ones((p, p), dtype=np.float64) if mask is None else np.asarray(mask, dtype=np.float64)
    xi = np.ascontiguousarray(xi, dtype=np.float64)
    yi = np.ascontiguousarray(yi, dtype=np.float64)
    out_mask = np.ascontiguousarray(out_mask, dtype=np.float64)

    ib = 6 if version < 15 else 7

    _eval_field_kernel_numba(
        f,
        elm_data,
        xi,
        yi,
        result,
        out_mask,
        nelms,
        int(operation),
        int(p),
        float(xrange[0]),
        float(yrange[0]),
        float(xzero),
        float(zzero),
        float(dx),
        float(dy),
        bool(threed),
        int(ib),
        float(local_phi),
        bool(wall_mask),
        float(version),
        int(max(debug_spin, 0)),
    )

    if np.max(out_mask) == 1:
        if edge_val is None:
            edge_sum = 0.0
            n_edge = 0
            for i in range(p):
                for j in range(p):
                    if out_mask[i, j] != 0:
                        continue
                    is_edge = False
                    if i > 0 and out_mask[i, j] != out_mask[i - 1, j]:
                        is_edge = True
                    if i < p - 1 and out_mask[i, j] != out_mask[i + 1, j]:
                        is_edge = True
                    if j > 0 and out_mask[i, j] != out_mask[i, j - 1]:
                        is_edge = True
                    if j < p - 1 and out_mask[i, j] != out_mask[i, j + 1]:
                        is_edge = True
                    if is_edge:
                        edge_sum += result[i, j]
                        n_edge += 1
            edge_val = edge_sum / n_edge if n_edge > 0 else 0.0

        result = result + out_mask * float(edge_val)

    return EvalFieldResult(data=result, r=xi, z=yi, mask=out_mask, edge_val=float(edge_val or 0.0))
