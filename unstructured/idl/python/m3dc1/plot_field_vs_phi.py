from __future__ import annotations

import inspect
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from .contour_and_legend import contour_and_legend
from .field_at_point import field_at_point
from .read_field import read_field
from .read_parameter import read_parameter
from .time_name import time_name


_READ_FIELD_KW = set(inspect.signature(read_field).parameters.keys())
_READ_FIELD_EXPLICIT_KW = {
    "timeslices",
    "filename",
    "points",
    "tpoints",
    "phi",
    "complex",
    "return_meta",
}


def _mesh_planes(filename: str | Path, timeslices: int) -> np.ndarray:
    with h5py.File(str(filename), "r") as h5:
        group_name = time_name(int(timeslices))
        if group_name not in h5:
            group_name = "equilibrium"
        mg = h5[group_name]["mesh"]
        if "phi" in mg.attrs:
            return np.asarray(mg.attrs["phi"], dtype=float).reshape(-1)
    return np.asarray([], dtype=float)


def _mesh_phi_period(filename: str | Path, timeslices: int, itor: int) -> float:
    with h5py.File(str(filename), "r") as h5:
        group_name = time_name(int(timeslices))
        if group_name not in h5:
            group_name = "equilibrium"
        attrs = h5[group_name]["mesh"].attrs
        period = float(np.asarray(attrs["period"]).reshape(-1)[0])
        nperiods = int(np.asarray(attrs["nperiods"]).reshape(-1)[0])

    if nperiods > 1:
        if itor == 1:
            sector_period = 2.0 * np.pi / float(nperiods)
        else:
            rzero = float(read_parameter("rzero", filename=filename))
            sector_period = 2.0 * np.pi * rzero / float(nperiods)
        period = min(period, sector_period)
    return period


def plot_field_vs_phi(
    field: str,
    timeslices=0,
    *,
    cutx=None,
    cutz=None,
    mesh: bool = False,
    phirange=None,
    tpoints: int | None = None,
    filename: str | Path = "C1.h5",
    points: int = 200,
    title: str | None = None,
    units: str | None = None,
    **kwargs,
):
    """
    Plot a field along an R or Z cut as a function of toroidal angle.
    """
    if (cutx is None) == (cutz is None):
        raise ValueError("plot_field_vs_phi requires exactly one of cutx or cutz.")

    read_kwargs = {k: v for k, v in kwargs.items() if k in _READ_FIELD_KW and k not in _READ_FIELD_EXPLICIT_KW}
    plot_kwargs = {k: v for k, v in kwargs.items() if k not in _READ_FIELD_KW}

    filename = Path(filename)
    slice_idx = int(timeslices)
    icomplex = int(read_parameter("icomplex", filename=filename))
    itor = int(read_parameter("itor", filename=filename))

    if phirange is None:
        phi_max = _mesh_phi_period(filename, slice_idx, itor)
        if itor == 1:
            phi_max = phi_max * 180.0 / np.pi
        phis = np.linspace(0.0, phi_max, int(tpoints or 200), endpoint=True)
    else:
        pr = np.asarray(phirange, dtype=float).reshape(-1)
        if pr.size < 2:
            raise ValueError("phirange must contain two values.")
        phis = np.linspace(float(pr[0]), float(pr[1]), int(tpoints or 200), endpoint=True)

    if slice_idx == -1:
        ntor = 0
        meta = read_field(
            field,
            timeslices=-1,
            filename=filename,
            points=points,
            return_meta=True,
            **read_kwargs,
        )
        field_stack = None
        base_field = np.asarray(meta.data)
    elif icomplex == 1:
        ntor = int(read_parameter("ntor", filename=filename))
        meta = read_field(
            field,
            timeslices=slice_idx,
            filename=filename,
            points=points,
            complex=True,
            return_meta=True,
            **read_kwargs,
        )
        field_stack = None
        base_field = np.asarray(meta.data, dtype=np.complex128)
    else:
        ntor = 0
        if tpoints is None:
            tpoints = 20
            phis = np.linspace(float(phis[0]), float(phis[-1]), int(tpoints), endpoint=True)
        meta = read_field(
            field,
            timeslices=slice_idx,
            filename=filename,
            points=points,
            tpoints=int(phis.size),
            phi=phis,
            return_meta=True,
            **read_kwargs,
        )
        field_stack = np.asarray(meta.data)
        base_field = None

    x = np.asarray(meta.r, dtype=float).reshape(-1)
    z = np.asarray(meta.z, dtype=float).reshape(-1)
    n = int(x.size)
    if n < 2:
        raise ValueError("plot_field_vs_phi requires at least two grid points.")

    rrange = np.asarray([np.nanmin(x), np.nanmax(x)], dtype=float)
    zrange = np.asarray([np.nanmin(z), np.nanmax(z)], dtype=float)
    if cutz is not None:
        rr = np.linspace(float(rrange[0]), float(rrange[1]), n)
        zz = np.full(n, float(cutz), dtype=float)
        cut_name = "Z"
        cut_value = float(cutz)
        v = rr
        xtitle = "R (m)"
    else:
        rr = np.full(n, float(cutx), dtype=float)
        zz = np.linspace(float(zrange[0]), float(zrange[1]), n)
        cut_name = "R"
        cut_value = float(cutx)
        v = zz
        xtitle = "Z (m)"

    out = np.zeros((n, phis.size), dtype=float)
    for j, phi in enumerate(phis):
        if slice_idx == -1:
            ff = np.asarray(base_field, dtype=float)
        elif icomplex == 1:
            phi_phase = float(phi) * np.pi / 180.0 if itor == 1 else float(phi)
            ff = np.real(np.asarray(base_field) * np.exp(1j * float(ntor) * phi_phase))
        else:
            ff = np.asarray(field_stack[j, :, :], dtype=float)
        out[:, j] = np.asarray(field_at_point(ff, x, z, rr, zz), dtype=float).reshape(-1)

    symbol = str(meta.symbol)
    unit_text = str(meta.units if units is None else units)
    label = symbol
    if unit_text:
        label = f"{label} ({unit_text})"
    label = f"{label} at {cut_name} = {cut_value:g}"
    ytitle = "Phi (deg)" if itor == 1 else "y (m)"

    contour_and_legend(
        out,
        v,
        phis,
        xtitle=xtitle,
        ytitle=ytitle,
        title="" if title is None else str(title),
        label=label,
        **plot_kwargs,
    )
    fig = plt.gcf()
    ax = plt.gca()

    if mesh:
        planes = _mesh_planes(filename, slice_idx)
        if planes.size > 0:
            plane_vals = planes * 180.0 / np.pi if itor == 1 else planes
            for plane in plane_vals:
                ax.axhline(float(plane), color="0.55", linewidth=0.8)

    return fig, ax
