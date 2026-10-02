from __future__ import annotations

import inspect

import numpy as np
from matplotlib.figure import Figure
from matplotlib.path import Path as MatplotlibPath

from .read_field import read_field
from .read_lcfs import read_lcfs

_READ_FIELD_KW = set(inspect.signature(read_field).parameters.keys())
_READ_LCFS_KW = set(inspect.signature(read_lcfs).parameters.keys())


def get_lcfs(psi=None, x=None, z=None, *, psival=None, axis=None, filename="C1.h5", points=200, slice=0, **kwargs):
    """
    Return an approximate LCFS contour path as shape (2, N).
    """
    field_kwargs = {k: v for k, v in kwargs.items() if k in _READ_FIELD_KW}
    lcfs_kwargs = {k: v for k, v in kwargs.items() if k in _READ_LCFS_KW}
    mesh_mask = None
    if psi is None or x is None or z is None:
        pmeta = read_field(
            "psi",
            filename=filename,
            timeslices=slice,
            points=points,
            equilibrium=True,
            return_meta=True,
            **field_kwargs,
        )
        psi2d = np.asarray(pmeta.data)[0, :, :] if np.asarray(pmeta.data).ndim == 3 else np.asarray(pmeta.data)
        xv = np.asarray(pmeta.r, dtype=float).reshape(-1)
        zv = np.asarray(pmeta.z, dtype=float).reshape(-1)
        mesh_mask = pmeta.mask
    else:
        arr = np.asarray(psi)
        psi2d = arr[0, :, :] if arr.ndim == 3 else arr
        xv = np.asarray(x, dtype=float).reshape(-1)
        zv = np.asarray(z, dtype=float).reshape(-1)

    lc = None
    if psival is None or axis is None:
        lc = read_lcfs(filename=filename, slice=slice, return_meta=True, **lcfs_kwargs)
    if psival is None:
        psival = float(lc.psilim)
    if axis is None and lc is not None:
        axis = np.asarray(lc.axis, dtype=float).reshape(-1)[:2]

    contour_data = np.asarray(psi2d)
    if mesh_mask is not None:
        mask = np.asarray(mesh_mask)
        if mask.ndim == 3:
            mask = mask[0]
        if mask.shape == contour_data.shape:
            contour_data = np.ma.masked_where(mask == 1, contour_data)

    fig = Figure()
    ax = fig.add_subplot(111)
    cs = ax.contour(
        xv,
        zv,
        contour_data.T,
        levels=[float(psival)],
        corner_mask=False,
    )
    try:
        if not cs.collections or not cs.collections[0].get_paths():
            return np.zeros((2, 0), dtype=float)
        paths = cs.collections[0].get_paths()

        def _path_area(path: MatplotlibPath) -> float:
            vertices = np.asarray(path.vertices, dtype=float)
            if vertices.shape[0] < 3:
                return 0.0
            return 0.5 * abs(
                np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1))
                - np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1))
            )

        containing = []
        if axis is not None:
            axis_point = np.asarray(axis, dtype=float).reshape(-1)
            if axis_point.size >= 2 and np.all(np.isfinite(axis_point[:2])):
                containing = [path for path in paths if path.contains_point(axis_point[:2])]
        candidates = containing if containing else paths
        selected = max(candidates, key=_path_area)
        return np.asarray(selected.vertices.T, dtype=float)
    finally:
        fig.clear()
