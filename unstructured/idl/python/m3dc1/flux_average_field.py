from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .field_at_point import field_at_point


@dataclass
class FluxAverageFieldResult:
    data: np.ndarray
    flux: np.ndarray
    nflux: np.ndarray
    area: np.ndarray
    volume: np.ndarray
    r0: float

def flux_average_field(
    field,
    psi=None,
    *,
    x=None,
    z=None,
    t=None,
    bins: int | None = None,
    flux=None,
    nflux=None,
    area=None,
    volume=None,
    r0=None,
    fc=None,
    filename: str = "C1.h5",
    surface_weight: bool = False,
    integrate: bool = False,
    return_meta: bool = False,
    **kwargs,
):
    """Flux-surface average using provided flux coordinates or generated ones."""
    f = np.asarray(field)
    if f.ndim == 2:
        f = f[None, :, :]
    if f.ndim != 3:
        raise ValueError(f"flux_average_field expects 2D/3D field, got {f.shape}")

    if fc is None:
        from .flux_coordinates import flux_coordinates
        if psi is None:
            raise ValueError("psi must be provided when fc is not provided.")
        if x is None or z is None:
            raise ValueError("x and z must be provided when fc is not provided.")
        b = bins if bins is not None else f.shape[1]
        fc = flux_coordinates(psi0=psi, x=np.asarray(x), z=np.asarray(z), fbins=int(b), tbins=int(b), filename=filename, **kwargs)

    vals = np.asarray(field_at_point(f, np.asarray(x), np.asarray(z), fc.r, fc.z))
    if vals.ndim == 2:
        vals = vals[None, :, :]

    out = np.zeros((vals.shape[0], fc.n), dtype=np.complex128 if np.iscomplexobj(vals) else float)
    for j in range(fc.n):
        w = np.asarray(fc.j[:, j], dtype=float)
        den = np.sum(w)
        if abs(den) < np.finfo(float).tiny:
            continue
        out[:, j] = np.sum(vals[:, :, j] * w[None, :], axis=1) / den

    if integrate:
        ps = np.asarray(fc.psi_norm, dtype=float)
        dV = np.asarray(fc.dV_dchi, dtype=float)
        for k in range(out.shape[0]):
            acc = np.zeros(fc.n, dtype=out.dtype)
            acc[0] = out[k, 0] * dV[0] * ps[0] / 2.0
            for j in range(1, fc.n):
                acc[j] = acc[j - 1] + (out[k, j] * dV[j] + out[k, j - 1] * dV[j - 1]) * (ps[j] - ps[j - 1]) / 2.0
            out[k, :] = acc

    if surface_weight:
        out = out * np.asarray(fc.area, dtype=float)[None, :]

    result = out[0, :] if out.shape[0] == 1 else out
    if return_meta:
        return FluxAverageFieldResult(
            data=np.asarray(result),
            flux=np.asarray(fc.psi),
            nflux=np.asarray(fc.psi_norm),
            area=np.asarray(fc.area),
            volume=np.asarray(fc.V),
            r0=float(fc.r0),
        )
    return result
