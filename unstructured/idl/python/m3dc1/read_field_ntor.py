from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .read_field import read_field
from .read_parameter import read_parameter


@dataclass
class FieldNTorResult:
    data: np.ndarray
    n: np.ndarray
    r: np.ndarray
    z: np.ndarray
    symbol: str = ""
    units: str = ""
    time: float = 0.0
    mask: np.ndarray | None = None


def read_field_ntor(
    name: str,
    timeslices=0,
    *,
    filename: str | Path = "C1.h5",
    points: int = 200,
    tpoints: int | None = None,
    ntor=None,
    xrange=None,
    yrange=None,
    linear: bool = False,
    logical: bool = False,
    phi: float = 0.0,
    operation: int = 1,
    fac: float | None = None,
    linfac: float = 1.0,
    equilibrium: bool = False,
    cgs: bool = False,
    mks: bool = False,
    wall_mask: bool = False,
    edge_val=None,
    **kwargs,
) -> FieldNTorResult:
    """
    Read a field using read_field() and return toroidal Fourier components
    on the original 2D (r, z) grid without any poloidal Fourier transform.
    """
    threed = int(read_parameter("3d", filename=filename, cgs=cgs, mks=mks))
    if threed != 1:
        raise ValueError("read_field_ntor requires a file with 3d=1.")
    if threed == 1 and tpoints is None:
        tpoints = 32

    field_kwargs = dict(kwargs)
    if threed == 1 and tpoints is not None:
        field_kwargs["tpoints"] = int(tpoints)

    meta = read_field(
        name,
        timeslices=timeslices,
        filename=filename,
        points=points,
        xrange=xrange,
        yrange=yrange,
        equilibrium=equilibrium,
        operation=operation,
        complex=False,
        fac=fac,
        linear=linear,
        linfac=linfac,
        cgs=cgs,
        mks=mks,
        phi=phi,
        wall_mask=wall_mask,
        logical=logical,
        edge_val=edge_val,
        return_meta=True,
        **field_kwargs,
    )

    field = np.asarray(meta.data)
    if field.ndim == 2:
        field = field[None, :, :]
    if field.ndim != 3:
        raise ValueError(f"read_field_ntor expects 2D or 3D field data, got shape {field.shape}.")

    data = np.asarray(field, dtype=np.complex128)
    data = np.fft.fft(data, axis=0)
    nvals = np.arange(data.shape[0], dtype=int)

    if ntor is not None:
        req = np.asarray(ntor, dtype=int).reshape(-1)
        idx = np.asarray([int(val) % int(data.shape[0]) for val in req], dtype=int)
        data = data[idx, :, :]
        nvals = req

    if float(phi) != 0.0:
        itor = int(read_parameter("itor", filename=filename, cgs=cgs, mks=mks))
        phi_phase = float(phi)
        if itor == 1:
            phi_phase = phi_phase * np.pi / 180.0
        data = data * np.exp(1j * np.asarray(nvals, dtype=float)[:, None, None] * phi_phase)

    return FieldNTorResult(
        data=np.asarray(data, dtype=np.complex128),
        n=np.asarray(nvals, dtype=int),
        r=np.asarray(meta.r, dtype=float),
        z=np.asarray(meta.z, dtype=float),
        symbol=str(meta.symbol),
        units=str(meta.units),
        time=float(meta.time),
        mask=None if meta.mask is None else np.asarray(meta.mask),
    )
