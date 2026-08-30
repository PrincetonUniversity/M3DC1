from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .field_at_point import field_at_point
from .flux_coordinates import FluxCoordinates, flux_coordinates
from .read_mesh import read_mesh
from .read_parameter import read_parameter


@dataclass
class FieldSpectrumResult:
    data: np.ndarray
    fc: FluxCoordinates
    m: np.ndarray
    n: np.ndarray
    x: np.ndarray
    xname: str
    symbol: str = ""
    units: str = ""


def _fft_normalized(values: np.ndarray, axis: int) -> np.ndarray:
    """Return Fourier coefficients normalized by the sampled axis length."""
    arr = np.asarray(values, dtype=np.complex128)
    n = int(arr.shape[axis])
    if n < 1:
        return np.asarray(arr, dtype=np.complex128)
    return np.fft.fft(arr, axis=axis) / float(n)


def _map_field_to_flux(field, x, z, fc: FluxCoordinates) -> np.ndarray:
    fld = np.asarray(field)
    if fld.ndim == 2:
        fld = fld[None, :, :]
    if fld.ndim != 3:
        raise ValueError(f"field_spectrum expects 2D or 3D field data, got shape {fld.shape}.")

    xv = np.asarray(x, dtype=float).reshape(-1)
    zv = np.asarray(z, dtype=float).reshape(-1)
    out = np.zeros((fld.shape[0], fc.n, fc.m), dtype=np.complex128)
    for k in range(fld.shape[0]):
        fk = np.asarray(fld[k, :, :])
        if np.iscomplexobj(fk):
            vals = field_at_point(np.real(fk), xv, zv, fc.r, fc.z) + 1j * field_at_point(np.imag(fk), xv, zv, fc.r, fc.z)
        else:
            vals = np.asarray(field_at_point(fk, xv, zv, fc.r, fc.z), dtype=float)
        out[k, :, :] = np.asarray(vals, dtype=np.complex128).T
    return out


def field_spectrum(
    field,
    x,
    z,
    *,
    m_val=None,
    psi_norm: bool = False,
    phi_norm: bool = False,
    rho: bool = False,
    psi0=None,
    i0=None,
    fc: FluxCoordinates | None = None,
    ignore_jacobian: bool = False,
    dpsi0_dx=None,
    dpsi0_dz=None,
    filename="C1.h5",
    nperiods: int | None = None,
    **kwargs,
) -> FieldSpectrumResult:
    """
    Python port of field_spectrum.pro.
    """
    ncoord = int(bool(psi_norm)) + int(bool(phi_norm)) + int(bool(rho))
    if ncoord > 1:
        raise TypeError("field_spectrum() accepts only one of 'psi_norm', 'phi_norm', or 'rho'.")
    param_kwargs = {k: v for k, v in kwargs.items() if k in {"cgs", "mks"}}
    if nperiods is None:
        nperiods = int(read_mesh(filename=filename, slice=int(kwargs.get("slice", 0))).nperiods)
    if int(nperiods) < 1:
        nperiods = 1
    if fc is None:
        coord_kwargs = dict(kwargs)
        if "slice" not in coord_kwargs:
            coord_linear = bool(read_parameter("linear", filename=filename, **param_kwargs))
            coord_kwargs["slice"] = -1 if coord_linear else 0
        fc = flux_coordinates(
            psi0=psi0,
            i0=i0,
            x=x,
            z=z,
            filename=filename,
            dpsi0_dx=dpsi0_dx,
            dpsi0_dz=dpsi0_dz,
            **coord_kwargs,
        )

    a = _map_field_to_flux(field, x, z, fc)
    b = np.transpose(a, (0, 2, 1)).astype(np.complex128, copy=False)
    b *= np.exp(-1j * np.asarray(fc.omega, dtype=float))[None, :, :]

    if ignore_jacobian:
        print("field_spectrum: Ignoring Jacobian")
    else:
        print("field_spectrum: Including normalized Jacobian")
        jac = np.asarray(fc.j, dtype=float)
        valid = np.isfinite(jac)
        count = np.sum(valid, axis=0)
        jac_mean = np.divide(
            np.sum(np.where(valid, jac, 0.0), axis=0),
            count,
            out=np.full(jac.shape[1], np.nan, dtype=float),
            where=count > 0,
        )
        jac_weight = np.divide(
            jac,
            jac_mean[None, :],
            out=np.full_like(jac, np.nan, dtype=float),
            where=np.isfinite(jac_mean[None, :]) & (np.abs(jac_mean[None, :]) > np.finfo(float).tiny),
        )
        b *= jac_weight[None, :, :]

    c = _fft_normalized(b, axis=1)
    m = -np.fft.fftshift(np.fft.fftfreq(int(fc.m)) * int(fc.m))
    d = np.conj(np.fft.fftshift(c, axes=1))

    nn = int(b.shape[0])
    if nn > 1:
        d = _fft_normalized(d, axis=0)
        n = np.arange(nn, dtype=int) * int(nperiods)
    else:
        d = np.asarray(d, dtype=np.complex128)
        n = np.asarray([int(read_parameter("ntor", filename=filename, **param_kwargs))], dtype=int)

    m = np.asarray(np.rint(m), dtype=int)
    poloidal_self = m == 0
    if int(fc.m) % 2 == 0:
        poloidal_self = poloidal_self | (np.abs(m) == int(fc.m) // 2)
    toroidal_self = np.zeros(nn, dtype=bool)
    toroidal_self[0] = True
    if nn % 2 == 0:
        toroidal_self[nn // 2] = True
    peak_scale = np.where(toroidal_self[:, None] & poloidal_self[None, :], 1.0, 2.0)
    d = np.asarray(d, dtype=np.complex128) * peak_scale[:, :, None]
    if m_val is not None:
        mvals = np.asarray(m_val, dtype=int).reshape(-1)
        idx = np.asarray([int(np.argmin(np.abs(m - mv))) for mv in mvals], dtype=int)
        d = d[:, idx, :]
        m = m[idx]

    if phi_norm:
        xvec = np.asarray(fc.phi_norm, dtype=float)
        xname = "phi_norm"
    elif rho:
        xvec = np.asarray(fc.rho, dtype=float)
        xname = "rho"
    else:
        xvec = np.asarray(fc.psi_norm, dtype=float)
        xname = "psi_norm"

    return FieldSpectrumResult(
        data=np.asarray(d, dtype=np.complex128),
        fc=fc,
        m=m,
        n=np.asarray(n, dtype=int),
        x=xvec,
        xname=xname,
    )
