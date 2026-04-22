from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .field_at_point import field_at_point
from .make_label import make_label
from .read_parameter import read_parameter
from .read_mesh import read_mesh


def _fill_masked_from_neighbors(data: np.ndarray, out_mask: np.ndarray | None) -> np.ndarray:
    arr = np.asarray(data, dtype=float)
    if out_mask is None or arr.ndim != 2:
        return arr

    mask = np.asarray(out_mask)
    if mask.shape != arr.shape:
        return arr

    src = np.asarray(arr, dtype=float)
    out = np.asarray(arr, dtype=float).copy()
    n0, n1 = src.shape
    for i in range(n0):
        for j in range(n1):
            if mask[i, j] != 1:
                continue
            vals = []
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    if di == 0 and dj == 0:
                        continue
                    ii = i + di
                    jj = j + dj
                    if ii < 0 or ii >= n0 or jj < 0 or jj >= n1:
                        continue
                    v = float(src[ii, jj])
                    if mask[ii, jj] != 1 and np.isfinite(v):
                        vals.append(v)
            if vals:
                out[i, j] = float(np.mean(np.asarray(vals, dtype=float)))
    return out


def plot_mesh(
    *,
    mesh=None,
    oplot=False,
    boundary=False,
    iso: bool = False,
    logical: bool = False,
    points: int = 200,
    phi: float = 0.0,
    filename="C1.h5",
    slice=0,
    **kwargs,
):
    """
    Plot triangular mesh edges.
    """
    m = read_mesh(filename=filename, slice=slice) if mesh is None else mesh
    el = np.asarray(m.elements, dtype=float)
    if oplot:
        ax = plt.gca()
        fig = ax.figure
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlabel(make_label("R", l0=1, **kwargs))
        ax.set_ylabel(make_label("Z", l0=1, **kwargs))

    if el.ndim != 2 or el.shape[0] < 6:
        if iso:
            ax.set_aspect("equal", adjustable="box")
        return fig, ax

    try:
        version = read_parameter("version", filename=filename)
        print(f"Output version = {version}")
    except Exception:
        pass
    try:
        igeometry = int(read_parameter("igeometry", filename=filename))
    except Exception:
        igeometry = 0

    rst_meta = None
    zst_meta = None
    if igeometry > 0 and not logical:
        # Local import avoids unnecessary dependency during default plotting path.
        from .read_field import read_field

        rst_meta = read_field(
            "rst",
            timeslices=-1,
            points=int(points),
            phi=float(phi),
            filename=filename,
            logical=True,
            return_meta=True,
        )
        zst_meta = read_field(
            "zst",
            timeslices=-1,
            points=int(points),
            phi=float(phi),
            filename=filename,
            logical=True,
            return_meta=True,
        )
        rst_meta.data = _fill_masked_from_neighbors(np.asarray(rst_meta.data), rst_meta.mask)
        zst_meta.data = _fill_masked_from_neighbors(np.asarray(zst_meta.data), zst_meta.mask)

    nelms = int(m.nelms)
    threed = bool(el.shape[0] > 8)
    nplanes = int(getattr(m, "nplanes", 0) or 0)
    nelms_plot = nelms
    if threed and nplanes > 0:
        nelms_plot = max(1, nelms // nplanes)
    minr = np.array([np.inf, np.inf], dtype=float)
    maxr = np.array([-np.inf, -np.inf], dtype=float)
    for i in range(nelms):
        if threed and i >= nelms_plot:
            break
        a = float(el[0, i])
        b = float(el[1, i])
        c = float(el[2, i])
        t = float(el[3, i])
        x0 = float(el[4, i])
        y0 = float(el[5, i])
        bound = int(el[6, i]) if el.shape[0] > 6 else 0

        p1 = np.asarray([x0, y0], dtype=float)
        p2 = p1 + np.asarray([(b + a) * np.cos(t), (b + a) * np.sin(t)], dtype=float)
        p3 = p1 + np.asarray([b * np.cos(t) - c * np.sin(t), b * np.sin(t) + c * np.cos(t)], dtype=float)

        if rst_meta is not None and zst_meta is not None:
            r_pts = np.asarray([p1[0], p2[0], p3[0]], dtype=float)
            z_pts = np.asarray([p1[1], p2[1], p3[1]], dtype=float)
            rst_pts = np.asarray(
                field_at_point(np.asarray(rst_meta.data), np.asarray(rst_meta.r), np.asarray(rst_meta.z), r_pts, z_pts),
                dtype=float,
            )
            zst_pts = np.asarray(
                field_at_point(np.asarray(zst_meta.data), np.asarray(zst_meta.r), np.asarray(zst_meta.z), r_pts, z_pts),
                dtype=float,
            )
            ok = np.isfinite(rst_pts) & np.isfinite(zst_pts)
            if np.any(ok):
                if ok[0]:
                    p1 = np.asarray([rst_pts[0], zst_pts[0]], dtype=float)
                if ok[1]:
                    p2 = np.asarray([rst_pts[1], zst_pts[1]], dtype=float)
                if ok[2]:
                    p3 = np.asarray([rst_pts[2], zst_pts[2]], dtype=float)

        minr[0] = min(minr[0], p1[0], p2[0], p3[0])
        minr[1] = min(minr[1], p1[1], p2[1], p3[1])
        maxr[0] = max(maxr[0], p1[0], p2[0], p3[0])
        maxr[1] = max(maxr[1], p1[1], p2[1], p3[1])

        if not boundary:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="0.35", linewidth=0.4)
            ax.plot([p2[0], p3[0]], [p2[1], p3[1]], color="0.35", linewidth=0.4)
            ax.plot([p3[0], p1[0]], [p3[1], p1[1]], color="0.35", linewidth=0.4)
            continue

        if bound & 1:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="k", linewidth=1.1)
        if bound & 2:
            ax.plot([p2[0], p3[0]], [p2[1], p3[1]], color="k", linewidth=1.1)
        if bound & 4:
            ax.plot([p3[0], p1[0]], [p3[1], p1[1]], color="k", linewidth=1.1)

    if iso:
        ax.set_aspect("equal", adjustable="box")

    print(f"Elements: {nelms}")
    print(f"sqrt(nodes) (estimated): {np.sqrt(nelms / 2.0)}")
    if np.isfinite(minr).all() and np.isfinite(maxr).all():
        print(f"Width: {maxr[0] - minr[0]}")
        print(f"Height: {maxr[1] - minr[1]}")

    return fig, ax
