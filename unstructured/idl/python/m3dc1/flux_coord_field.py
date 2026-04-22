from __future__ import annotations

import numpy as np

from .field_at_point import field_at_point
from .flux_coordinates import flux_coordinates


def flux_coord_field(
    field,
    psi,
    x,
    z,
    t=None,
    *,
    fbins: int = 200,
    tbins: int = 200,
    flux=None,
    angle=None,
    nflux=None,
    filename="C1.h5",
    slice=0,
    **kwargs,
):
    """Map field to flux coordinates using flux_coordinates()."""
    fld = np.asarray(field)
    if fld.ndim == 2:
        fld = fld[None, :, :]
    xv = np.asarray(x, dtype=float).reshape(-1)
    zv = np.asarray(z, dtype=float).reshape(-1)
    fc = flux_coordinates(
        slice=slice,
        fbins=fbins,
        tbins=tbins,
        psi0=psi,
        x=xv,
        z=zv,
        filename=filename,
        **kwargs,
    )
    out = np.zeros((fld.shape[0], fc.n, fc.m), dtype=float)
    for k in range(fld.shape[0]):
        out[k, :, :] = np.asarray(field_at_point(fld[k, :, :], xv, zv, fc.r, fc.z), dtype=float).T

    if flux is not None:
        flux[:] = np.asarray(fc.psi, dtype=float)
    if angle is not None:
        angle[:] = np.asarray(fc.theta, dtype=float)
    if nflux is not None:
        nflux[:] = np.asarray(fc.psi_norm, dtype=float)
    return out
