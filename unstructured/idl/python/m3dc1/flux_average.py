from __future__ import annotations

from dataclasses import dataclass
import inspect
import numpy as np

from .dimensions import dimensions
from .flux_average_field import flux_average_field
from .flux_coordinates import flux_coordinates
from .parse_units import parse_units
from .read_field import read_field
from .read_parameter import read_parameter

_READ_FIELD_KW = set(inspect.signature(read_field).parameters.keys())


@dataclass
class FluxAverageResult:
    data: np.ndarray
    title: str
    symbol: str
    units: str
    fc: object


def flux_average(
    field,
    timeslices=0,
    *,
    psi_norm: bool = False,
    phi_norm: bool = False,
    rho: bool = False,
    psi=None,
    i0=None,
    x=None,
    z=None,
    r0=None,
    flux=None,
    nflux=None,
    area=None,
    bins=None,
    points=None,
    name=None,
    symbol=None,
    units=None,
    integrate=False,
    complex=False,
    abs=False,
    phase=False,
    stotal=False,
    fac=None,
    fc=None,
    elongation=None,
    filename="C1.h5",
    linear=False,
    linfac=1.0,
    mks=False,
    cgs=False,
    return_meta=False,
    **kwargs,
):
    """
    Python port of flux_average.pro (core behavior used by plotting routines).
    """
    if "t" in kwargs:
        if timeslices not in (0, None):
            raise TypeError("flux_average() received both 'timeslices' and legacy 't'.")
        timeslices = kwargs.pop("t")
    ncoord = int(bool(psi_norm)) + int(bool(phi_norm)) + int(bool(rho))
    if ncoord > 1:
        raise TypeError("flux_average() accepts only one of 'psi_norm', 'phi_norm', or 'rho'.")

    read_kwargs = {k: v for k, v in kwargs.items() if k in _READ_FIELD_KW}
    coord_linear = bool(read_parameter("linear", filename=filename, cgs=cgs, mks=mks))
    coord_slice = -1 if coord_linear else int(timeslices or 0)

    if fc is None:
        if points is None:
            points = 200
        nb = int(bins if bins is not None else points)
        if psi is None or x is None or z is None:
            p = read_field(
                "psi",
                filename=filename,
                timeslices=coord_slice,
                points=int(points),
                equilibrium=(coord_slice < 0),
                return_meta=True,
                cgs=cgs,
                mks=mks,
                **read_kwargs,
            )
            psi = np.asarray(p.data)[0, :, :] if np.asarray(p.data).ndim == 3 else np.asarray(p.data)
            x = np.asarray(p.r, dtype=float).reshape(-1)
            z = np.asarray(p.z, dtype=float).reshape(-1)
        fc = flux_coordinates(
            psi0=psi,
            i0=i0,
            x=x,
            z=z,
            slice=coord_slice,
            points=int(points),
            fbins=nb,
            tbins=int(points),
            filename=filename,
            **kwargs,
        )
    if phi_norm:
        xvec = np.asarray(fc.phi_norm)
    elif rho:
        xvec = np.asarray(fc.rho)
    elif psi_norm:
        xvec = np.asarray(fc.psi_norm)
    else:
        xvec = np.asarray(fc.psi)
    if flux is not None:
        flux[:] = xvec
    if nflux is not None:
        nflux[:] = np.asarray(fc.psi_norm)
    if area is not None:
        area[:] = np.asarray(fc.area)

    if isinstance(field, str):
        fkey = field.strip().lower()
        if fkey in {"safety factor", "q"}:
            vals = np.abs(np.asarray(fc.q, dtype=float))
            title = "Safety Factor"
            sym = "$q$"
            u = ""
            if return_meta:
                return FluxAverageResult(data=np.asarray(vals), title=title, symbol=sym, units=u, fc=fc)
            return vals
        if fkey == "rho":
            vals = np.asarray(fc.rho, dtype=float)
            title = "rho"
            sym = "$\\rho$"
            u = parse_units(dimensions(l0=1), cgs=cgs, mks=mks)
            if return_meta:
                return FluxAverageResult(data=np.asarray(vals), title=title, symbol=sym, units=u, fc=fc)
            return vals
        if fkey == "flux_t":
            vals = np.asarray(fc.flux_tor, dtype=float)
            if return_meta:
                return FluxAverageResult(data=np.asarray(vals), title="Toroidal Flux", symbol="$\\Phi_t$", units="", fc=fc)
            return vals
        if fkey == "flux_p":
            vals = np.asarray(fc.flux_pol, dtype=float)
            if return_meta:
                return FluxAverageResult(data=np.asarray(vals), title="Poloidal Flux", symbol="$\\Phi_p$", units="", fc=fc)
            return vals
        if fkey == "volume":
            vals = np.asarray(fc.V, dtype=float)
            u = parse_units(dimensions(l0=3), cgs=cgs, mks=mks)
            if return_meta:
                return FluxAverageResult(data=np.asarray(vals), title="Volume", symbol="$V$", units=u, fc=fc)
            return vals

        meta = read_field(
            field,
            timeslices=int(timeslices or 0),
            points=int(points or fc.m),
            linear=linear,
            complex=complex,
            abs=abs,
            phase=phase,
            fac=fac,
            linfac=linfac,
            filename=filename,
            cgs=cgs,
            mks=mks,
            return_meta=True,
            **read_kwargs,
        )
        arr = np.asarray(meta.data)
        farr = arr[0, :, :] if arr.ndim == 3 else arr
        vals = flux_average_field(
            farr,
            psi=psi,
            x=np.asarray(x),
            z=np.asarray(z),
            fc=fc,
            bins=bins,
            integrate=integrate,
            surface_weight=stotal,
            filename=filename,
            **kwargs,
        )
        title = str(meta.symbol)
        sym = str(meta.symbol) if not integrate and not stotal else str(meta.symbol)
        u = str(meta.units)
        if return_meta:
            return FluxAverageResult(data=np.asarray(vals), title=title, symbol=sym, units=u, fc=fc)
        return np.asarray(vals)

    vals = flux_average_field(
        field,
        psi=psi,
        x=np.asarray(x),
        z=np.asarray(z),
        fc=fc,
        bins=bins,
        integrate=integrate,
        surface_weight=stotal,
        filename=filename,
        **kwargs,
    )
    if return_meta:
        return FluxAverageResult(data=np.asarray(vals), title="", symbol="", units="", fc=fc)
    return np.asarray(vals)
