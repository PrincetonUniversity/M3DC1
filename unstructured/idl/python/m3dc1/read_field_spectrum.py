from __future__ import annotations

from pathlib import Path

import numpy as np

from .field_spectrum import FieldSpectrumResult, field_spectrum
from .flux_coordinates import flux_coordinates
from .read_field import read_field
from .read_parameter import read_parameter


def read_field_spectrum(
    name: str,
    timeslices=0,
    *,
    filename: str | Path = "C1.h5",
    points: int = 200,
    tpoints: int | None = None,
    ntor: int | None = None,
    xrange=None,
    yrange=None,
    linear: bool = False,
    complex: bool = True,
    logical: bool = False,
    phi: float = 0.0,
    operation: int = 1,
    fac: float | None = None,
    linfac: float = 1.0,
    equilibrium: bool = False,
    cgs: bool = False,
    mks: bool = False,
    m_val=None,
    psi_norm: bool = False,
    phi_norm: bool = False,
    rho: bool = False,
    psi0=None,
    i0=None,
    fc=None,
    ignore_jacobian: bool = False,
    dpsi0_dx=None,
    dpsi0_dz=None,
    pest: bool = False,
    boozer: bool = False,
    hamada: bool = False,
    fast: bool = False,
    fbins: int | None = None,
    tbins: int | None = None,
    psin_range=None,
    wall_mask: bool = False,
    edge_val=None,
    **kwargs,
) -> FieldSpectrumResult:
    """
    Read a field using read_field() and return its field_spectrum().
    """
    if "m_vals" in kwargs:
        if m_val is not None:
            raise TypeError("read_field_spectrum() got both 'm_val' and 'm_vals'. Use only one.")
        m_val = kwargs.pop("m_vals")

    threed = int(read_parameter("3d", filename=filename, cgs=cgs, mks=mks))
    itor = int(read_parameter("itor", filename=filename, cgs=cgs, mks=mks))
    if not boozer and not hamada and not pest and not fast:
        if itor == 1:
            pest = True
    if ntor is None:
        ntor = int(read_parameter("ntor", filename=filename, cgs=cgs, mks=mks))
    if threed == 1 and tpoints is None:
        tpoints = 32
    field_kwargs = dict(kwargs)
    if threed == 1 and tpoints is not None:
        field_kwargs["tpoints"] = int(tpoints)

    try:
        meta = read_field(
            name,
            timeslices=timeslices,
            filename=filename,
            points=points,
            xrange=xrange,
            yrange=yrange,
            equilibrium=equilibrium,
            operation=operation,
            complex=complex,
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
    except KeyError:
        if not complex:
            raise
        print(f"read_field_spectrum: complex read for '{name}' is unavailable; falling back to real field.")
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
    x = np.asarray(meta.r, dtype=float).reshape(-1)
    z = np.asarray(meta.z, dtype=float).reshape(-1)
    timeslice_idx = -1 if timeslices is None else int(np.asarray(timeslices).reshape(-1)[0])

    select_m = None if threed == 1 else m_val

    coord_linear = bool(read_parameter("linear", filename=filename, cgs=cgs, mks=mks))
    coord_slice = -1 if coord_linear else timeslice_idx
    if fc is None:
        fc = flux_coordinates(
            psi0=psi0,
            i0=i0,
            x=x,
            z=z,
            filename=filename,
            dpsi0_dx=dpsi0_dx,
            dpsi0_dz=dpsi0_dz,
            slice=coord_slice,
            pest=pest,
            boozer=boozer,
            hamada=hamada,
            fast=fast,
            fbins=fbins,
            tbins=tbins,
            psin_range=psin_range,
            cgs=cgs,
            mks=mks,
        )

    spec = field_spectrum(
        field,
        x,
        z,
        m_val=select_m,
        psi_norm=psi_norm,
        phi_norm=phi_norm,
        rho=rho,
        psi0=psi0,
        i0=i0,
        fc=fc,
        ignore_jacobian=ignore_jacobian,
        dpsi0_dx=dpsi0_dx,
        dpsi0_dz=dpsi0_dz,
        filename=filename,
        pest=pest,
        boozer=boozer,
        hamada=hamada,
        fast=fast,
        fbins=fbins,
        tbins=tbins,
        psin_range=psin_range,
        cgs=cgs,
        mks=mks,
    )
    if threed == 1:
        nvals = np.asarray(spec.n, dtype=int).reshape(-1)
        if nvals.size > 0:
            target_ntor = int(ntor)
            idx = target_ntor % int(nvals.size)
            data = np.asarray(spec.data, dtype=np.complex128)[idx : idx + 1, ...]
            mvals = np.asarray(spec.m, dtype=int).reshape(-1)

            spec.data = data
            spec.n = np.asarray([target_ntor], dtype=int)

            if m_val is not None:
                req_m = np.asarray(m_val, dtype=int).reshape(-1)
                midx = np.asarray([int(np.argmin(np.abs(mvals - mv))) for mv in req_m], dtype=int)
                spec.data = np.asarray(spec.data)[:, midx, :]
                spec.m = mvals[midx]
    spec.symbol = str(meta.symbol)
    spec.units = str(meta.units)
    return spec
