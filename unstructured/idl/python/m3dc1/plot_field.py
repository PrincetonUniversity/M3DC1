from __future__ import annotations

import inspect
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .contour_and_legend import contour_and_legend
from .field_at_point import field_at_point
from .flux_at_q import flux_at_q
from .flux_coord_field import flux_coord_field
from .make_label import make_label
from .nulls import nulls
from .plot_coils import plot_coils
from .plot_flux_contour import plot_flux_contour
from .plot_lcfs import plot_lcfs
from .plot_mesh import plot_mesh
from .plot_wall_regions import plot_wall_regions
from .read_field import FieldResult, read_field

_READ_FIELD_KW = set(inspect.signature(read_field).parameters.keys())


def _as_3d(arr: np.ndarray) -> np.ndarray:
    a = np.asarray(arr)
    if a.ndim == 3:
        return a
    if a.ndim == 2:
        return a[None, :, :]
    raise ValueError(f"Expected 2D/3D field data, got shape {a.shape}.")


def _save_xy(path: str | Path, xx: np.ndarray, yy: np.ndarray) -> None:
    x = np.asarray(xx, dtype=float).reshape(-1)
    y = np.asarray(yy, dtype=float).reshape(-1)
    n = min(x.size, y.size)
    np.savetxt(str(path), np.column_stack([x[:n], y[:n]]), fmt="%16.6e")


def _primary_timeslice(timeslices) -> int:
    if isinstance(timeslices, (list, tuple, np.ndarray)):
        arr = np.asarray(timeslices).reshape(-1)
        if arr.size == 0:
            return 0
        return int(arr[0])
    return int(timeslices)


def plot_field(
    name,
    timeslices=0,
    *,
    points: int = 200,
    mesh=None,
    iso: bool = False,
    mcolor=None,
    lcfs: bool = False,
    title: str | None = None,
    units: str | None = None,
    range=None,
    xrange=None,
    yrange=None,
    linear: bool = False,
    xlim=None,
    cutx=None,
    cutz=None,
    mpeg=None,
    linfac: float = 1.0,
    mask_val=None,
    boundary: bool = False,
    q_contours=None,
    overplot: bool = False,
    phi: float = 0.0,
    logical: bool = False,
    realtime=None,
    levels=None,
    colorbar: bool = True,
    phase: bool = False,
    abs: bool = False,
    real: bool = False,
    imaginary: bool = False,
    operation: int = 1,
    magcoord: bool = False,
    outfile: str | Path | None = None,
    fac=None,
    filename: str | Path = "C1.h5",
    psin: bool = False,
    coils: bool = False,
    axis: bool = False,
    wall_regions: bool = False,
    cmap: str | None = None,
    **kwargs,
):
    """
    Python port of plot_field.pro.
    """
    if timeslices is None:
        timeslices = 0
    if "x" in kwargs or "y" in kwargs:
        raise TypeError("plot_field() no longer accepts 'x'/'y' arguments.")
    if cmap is None and "cmap" in kwargs:
        cmap = kwargs.pop("cmap")
    input_r = kwargs.pop("r", None)
    input_z = kwargs.pop("z", None)
    input_symbol = kwargs.pop("symbol", None)
    input_mask = kwargs.pop("mask", None)
    input_time = kwargs.pop("time", kwargs.pop("realtime", None))
    notitle = title is None
    read_kwargs = {k: v for k, v in kwargs.items() if k in _READ_FIELD_KW}
    plot_kwargs = dict(kwargs)
    label_kwargs = {}
    if "cgs" in read_kwargs:
        label_kwargs["cgs"] = read_kwargs["cgs"]
    if "mks" in read_kwargs:
        label_kwargs["mks"] = read_kwargs["mks"]

    main_timeslices = timeslices
    primary_slice = _primary_timeslice(timeslices)

    complex_flag = bool(phase or abs or real or imaginary)
    if isinstance(name, str):
        meta = read_field(
            name,
            timeslices=main_timeslices,
            points=points,
            xrange=xrange,
            yrange=yrange,
            linear=linear,
            phi=phi,
            logical=logical,
            operation=operation,
            complex=complex_flag,
            filename=filename,
            linfac=linfac,
            fac=fac,
            return_meta=True,
            **read_kwargs,
        )
        field = np.asarray(meta.data)
        xvec = np.asarray(meta.r, dtype=float).reshape(-1)
        yvec = np.asarray(meta.z, dtype=float).reshape(-1)
        mask = meta.mask
        fieldname = str(meta.symbol)
        if units is None:
            units = str(meta.units)
        if realtime is None:
            realtime = meta.time
    elif isinstance(name, FieldResult):
        field = np.asarray(name.data)
        xvec = np.asarray(name.r, dtype=float).reshape(-1)
        yvec = np.asarray(name.z, dtype=float).reshape(-1)
        mask = name.mask if input_mask is None else input_mask
        fieldname = str(name.symbol) if input_symbol is None else str(input_symbol)
        if units is None:
            units = str(name.units)
        if realtime is None:
            realtime = name.time
    else:
        if input_r is None or input_z is None:
            raise TypeError("plot_field requires 'r' and 'z' when plotting precomputed field data.")
        field = np.asarray(name)
        if field.ndim != 2:
            raise TypeError("plot_field array input must be a 2D array.")
        xvec = np.asarray(input_r, dtype=float).reshape(-1)
        yvec = np.asarray(input_z, dtype=float).reshape(-1)
        mask = input_mask
        fieldname = "" if input_symbol is None else str(input_symbol)
        if realtime is None:
            realtime = input_time

    f3 = _as_3d(field)

    if phase:
        f3 = np.angle(f3, deg=True)
        units = "Degrees"
        fieldname = f"Phase({fieldname})"
    elif abs:
        f3 = np.abs(f3)
    elif real:
        f3 = np.real(f3)
    elif imaginary:
        f3 = np.imag(f3)

    if realtime is not None:
        print(f"time = {realtime}")

    bad = ~np.isfinite(f3)
    if np.any(bad):
        f3 = np.array(f3, copy=True)
        f3[bad] = 0.0

    if mask_val is not None and mask is not None:
        m = np.asarray(mask, dtype=float)
        for k in range(f3.shape[0]):
            f3[k, :, :] = f3[k, :, :] - m * (f3[k, :, :] - float(mask_val))

    if not (phase or abs or real or imaginary):
        f3 = np.real(f3)

    if notitle:
        title = fieldname

    if cutx is not None:
        xv = float(cutx)
        yy = yvec
        data = field_at_point(f3[0, :, :], xvec, yvec, np.full_like(yy, xv), yy)
        xtitle = make_label("Z", l0=1, **kwargs)
        if psin:
            psi = read_field("psi_norm", timeslices=main_timeslices, points=points, equilibrium=True, logical=logical, filename=filename, **read_kwargs)
            yy = np.asarray(field_at_point(np.asarray(psi)[0, :, :], xvec, yvec, np.full_like(yvec, xv), yvec), dtype=float)
        if not overplot:
            plt.figure(figsize=(7, 4))
        plt.plot(yy, np.asarray(data).reshape(-1))
        plt.title(str(title))
        plt.xlabel(str(xtitle))
        if units:
            plt.ylabel(str(units))
        if outfile is not None:
            _save_xy(outfile, yy, np.asarray(data).reshape(-1))
        ax = plt.gca()
        if iso:
            ax.set_aspect("equal", adjustable="box")
        return ax.figure, ax

    if cutz is not None:
        zv = float(cutz)
        xx = xvec
        data = field_at_point(f3[0, :, :], xvec, yvec, xx, np.full_like(xx, zv))
        xtitle = make_label("R", l0=1, **kwargs)
        if psin:
            psi = read_field("psi_norm", timeslices=main_timeslices, points=points, equilibrium=True, logical=logical, filename=filename, **read_kwargs)
            xx = np.asarray(field_at_point(np.asarray(psi)[0, :, :], xvec, yvec, xvec, np.full_like(xvec, zv)), dtype=float)
        if not overplot:
            plt.figure(figsize=(7, 4))
        plt.plot(xx, np.asarray(data).reshape(-1))
        plt.title(str(title))
        plt.xlabel(str(xtitle))
        if units:
            plt.ylabel(str(units))
        if outfile is not None:
            _save_xy(outfile, xx, np.asarray(data).reshape(-1))
        ax = plt.gca()
        if iso:
            ax.set_aspect("equal", adjustable="box")
        return ax.figure, ax

    if magcoord:
        psi = read_field("psi", timeslices=main_timeslices, points=points, equilibrium=True, logical=logical, filename=filename, **read_kwargs)
        mapped = flux_coord_field(
            f3,
            np.asarray(psi),
            xvec,
            yvec,
            primary_slice,
            fbins=points,
            tbins=points,
            filename=filename,
            slice=primary_slice,
            **plot_kwargs,
        )
        angle = np.linspace(0.0, 360.0, int(points))
        nflux = np.linspace(0.0, 1.0, int(points))
        contour_and_legend(
            mapped[0, :, :],
            angle,
            nflux,
            title=str(title),
            label=units or "",
            levels=levels,
            xtitle="Angle (Degrees)",
            ytitle="psi",
            cmap=cmap,
            range=range,
            colorbar=colorbar,
            overplot=overplot,
            **plot_kwargs,
        )
        ax = plt.gca()
        if iso:
            ax.set_aspect("equal", adjustable="box")
        return ax.figure, ax

    contour_and_legend(
        f3[0, :, :],
        xvec,
        yvec,
        title=str(title),
        label=units or "",
        levels=levels,
        cmap=cmap,
        xtitle=make_label("R", l0=1, **label_kwargs),
        ytitle=make_label("Z", l0=1, **label_kwargs),
        range=range,
        colorbar=colorbar,
        overplot=overplot,
        **plot_kwargs,
    )

    mesh_obj = None if isinstance(mesh, bool) else mesh
    show_mesh = bool(mesh) if isinstance(mesh, bool) else (mesh is not None)

    if boundary:
        plot_mesh(
            mesh=mesh_obj,
            oplot=True,
            boundary=True,
            logical=logical,
            phi=phi,
            filename=filename,
            slice=primary_slice,
            **plot_kwargs,
        )
    elif show_mesh:
        plot_mesh(
            mesh=mesh_obj,
            oplot=True,
            boundary=False,
            logical=logical,
            phi=phi,
            filename=filename,
            slice=primary_slice,
            **plot_kwargs,
        )

    if wall_regions:
        plot_wall_regions(filename=filename, slice=primary_slice, over=True, **plot_kwargs)

    if q_contours is not None:
        fval = flux_at_q(q_contours, points=points, filename=filename, **plot_kwargs)
        if np.asarray(fval).size > 0 and float(np.asarray(fval).reshape(-1)[0]) != 0.0:
            plot_flux_contour(fval, points=points, overplot=True, filename=filename, slice=primary_slice, **plot_kwargs)

    if axis:
        ax, _ = nulls(axis=True, xpoints=True, filename=filename, slice=primary_slice, **plot_kwargs)
        dx = (xvec.max() - xvec.min()) / 50.0
        dy = (yvec.max() - yvec.min()) / 50.0
        plt.plot([ax[0] - dx, ax[0] + dx], [ax[1] - dy, ax[1] + dy], color="tab:red")
        plt.plot([ax[0] - dx, ax[0] + dx], [ax[1] + dy, ax[1] - dy], color="tab:red")

    if lcfs:
        plot_lcfs(over=True, filename=filename, slice=primary_slice, points=points, **plot_kwargs)

    if coils:
        plot_coils(filename=filename, overplot=True, **plot_kwargs)

    if xlim is not None:
        plt.xlim(xlim)
    if mpeg is not None:
        plt.savefig(str(mpeg), dpi=150)
    if outfile is not None and cutx is None and cutz is None:
        np.save(str(outfile), f3)
    ax = plt.gca()
    if iso:
        ax.set_aspect("equal", adjustable="box")
    return ax.figure, ax
