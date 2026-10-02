from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .plot_mesh import plot_mesh
from .read_field import read_field
from .read_lcfs import read_lcfs
from .read_parameter import read_parameter


def _flux_contour_slice(filename, slice_value: int | None) -> tuple[int, bool]:
    if slice_value is not None:
        return int(slice_value), False
    eqsubtract = int(read_parameter("eqsubtract", filename=filename))
    return (-1 if eqsubtract != 0 else 0), True


def plot_flux_contour(fval=True, *, filename="C1.h5", points=200, slice: int | None = None, overplot=False, closed=0, color="k", iso: bool = False, boundary: bool = False, xscale: float = 1.0, yscale: float = 1.0, psilim: float | None = None, **kwargs):
    """
    Plot psi contour(s) at requested flux values.

    With ``fval=True``, plot 10 levels from the magnetic-axis flux through the
    limiting flux, plus one equally spaced level beyond the limiting flux.
    ``psilim`` overrides the limiting ``psibound`` scalar for automatic levels.
    With ``slice=None``, read equilibrium psi from slice ``-1`` for
    ``eqsubtract=1`` or slice ``0`` for ``eqsubtract=0``. An explicit slice
    reads total psi at that timeslice.
    """
    if isinstance(fval, (bool, np.bool_)) and not fval:
        return None, None
    read_slice, equilibrium = _flux_contour_slice(filename, slice)
    p = read_field(
        "psi",
        filename=filename,
        timeslices=read_slice,
        points=points,
        equilibrium=equilibrium,
        return_meta=True,
        **kwargs,
    )
    psi = np.asarray(p.data)
    psi2d = psi[0, :, :] if psi.ndim == 3 else psi
    x = np.asarray(p.r, dtype=float).reshape(-1)
    z = np.asarray(p.z, dtype=float).reshape(-1)
    contour_data = np.asarray(psi2d)
    if p.mask is not None:
        mesh_mask = np.asarray(p.mask)
        if mesh_mask.ndim == 3:
            mesh_mask = mesh_mask[0]
        if mesh_mask.shape == contour_data.shape:
            contour_data = np.ma.masked_where(mesh_mask == 1, contour_data)

    if isinstance(fval, (bool, np.bool_)):
        lcfs_kwargs = {key: kwargs[key] for key in ("cgs", "mks") if key in kwargs}
        lc = read_lcfs(
            filename=filename,
            slice=read_slice,
            return_meta=True,
            **lcfs_kwargs,
        )
        limit_flux = lc.psilim if psilim is None else float(psilim)
        if not np.isfinite(limit_flux):
            raise ValueError(f"psilim must be finite, got {psilim!r}.")
        vals = np.linspace(lc.flux0, limit_flux, 10)
        vals = np.append(vals, limit_flux + (limit_flux - lc.flux0) / 9.0)
    else:
        vals = np.asarray(fval, dtype=float).reshape(-1)
    if vals.size == 0:
        return None, None

    if overplot:
        axis = plt.gca()
        figure = axis.figure
    else:
        figure, axis = plt.subplots(figsize=(6, 5))
    axis.contour(
        x * float(xscale),
        z * float(yscale),
        contour_data.T,
        levels=np.sort(vals),
        colors=color,
        linewidths=0.8,
        corner_mask=False,
    )
    if boundary:
        plot_mesh(
            mesh=kwargs.get("mesh"),
            oplot=True,
            boundary=True,
            logical=bool(kwargs.get("logical", False)),
            points=points,
            phi=float(kwargs.get("phi", 0.0)),
            filename=filename,
            slice=read_slice,
            xscale=xscale,
            yscale=yscale,
        )
    if iso:
        axis.set_aspect("equal", adjustable="box")
    return figure, axis
