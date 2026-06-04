from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np

from .plot_field import plot_field
from .read_field import FieldResult
from .read_field_ntor import read_field_ntor

_READ_FIELD_NTOR_KW = set(inspect.signature(read_field_ntor).parameters.keys())
_PLOT_FIELD_KW = set(inspect.signature(plot_field).parameters.keys())


def plot_field_ntor(
    name: str,
    timeslices=0,
    *,
    ntor=0,
    phi: float = 0.0,
    filename: str | Path = "C1.h5",
    points: int = 200,
    tpoints: int | None = None,
    title: str | None = None,
    units: str | None = None,
    **kwargs,
):
    """
    Read toroidal Fourier components with read_field_ntor() and plot one
    selected component with plot_field().
    """
    read_kwargs = {k: v for k, v in kwargs.items() if k in _READ_FIELD_NTOR_KW}
    plot_kwargs = {k: v for k, v in kwargs.items() if k in _PLOT_FIELD_KW}

    res = read_field_ntor(
        name,
        timeslices=timeslices,
        ntor=ntor,
        phi=phi,
        filename=filename,
        points=points,
        tpoints=tpoints,
        **read_kwargs,
    )

    data = np.asarray(res.data)
    nvals = np.asarray(res.n, dtype=int).reshape(-1)
    if data.ndim != 3 or data.shape[0] != 1:
        raise ValueError("plot_field_ntor requires a single selected ntor value.")

    nsel = int(nvals[0])
    plot_title = title if title is not None else f"{name}, n={nsel}"
    field = FieldResult(
        data=np.asarray(data[0, :, :]),
        symbol=str(res.symbol),
        units=str(res.units if units is None else units),
        dimensions=np.asarray([]),
        r=np.asarray(res.r, dtype=float),
        z=np.asarray(res.z, dtype=float),
        time=float(res.time),
        mask=None if res.mask is None else np.asarray(res.mask),
    )

    return plot_field(
        field,
        timeslices=timeslices,
        title=plot_title,
        filename=filename,
        **plot_kwargs,
    )
