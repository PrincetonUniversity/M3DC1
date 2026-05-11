from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import h5py
import numpy as np

from .convert_units import convert_units
from .create_map import create_map
from .dimensions import dimensions
from .a_bracket import a_bracket
from .delstar import delstar
from .eval_field import eval_field
from .field_data import field_data
from .flux_average_field import flux_average_field
from .get_slice_time import get_slice_time
from .hdf5_file_test import hdf5_file_test
from .dx import dx
from .dz import dz
from .laplacian import laplacian
from .map_field import map_field
from .parse_units import parse_units
from .radius_matrix import radius_matrix
from .read_lcfs import read_lcfs
from .read_mesh import read_mesh
from .read_parameter import read_parameter
from .read_gamma import read_gamma
from .s_bracket import s_bracket
from .time_name import time_name


@dataclass
class FieldResult:
    data: np.ndarray
    symbol: str
    units: str
    dimensions: np.ndarray
    r: np.ndarray
    z: np.ndarray
    time: float
    mask: np.ndarray | None = None


def _find_field_name_case_insensitive(fields: h5py.Group, name: str) -> str | None:
    n = name.lower()
    for k in fields.keys():
        if k.lower() == n:
            return k
    return None


def _read_field_dataset(time_group: h5py.Group, name: str) -> np.ndarray:
    obj = time_group["fields"][name]
    if isinstance(obj, h5py.Dataset):
        arr = np.asarray(obj[()])
    elif isinstance(obj, h5py.Group):
        for key in ("_data", "data", "value"):
            if key in obj and isinstance(obj[key], h5py.Dataset):
                arr = np.asarray(obj[key][()])
                break
        else:
            raise KeyError(f"Unsupported field storage layout for '{name}'.")
    else:
        raise KeyError(f"Unsupported field object type for '{name}'.")

    arr = np.asarray(arr, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"Expected field coefficients to be 2D, got {arr.shape}.")
    # IDL uses [n_coeff, n_elem]
    if arr.shape[0] > arr.shape[1]:
        arr = arr.T
    return arr


def _read_primitive(
    field_name: str,
    *,
    filename: str | Path,
    slice_idx: int,
    points: int,
    xrange,
    yrange,
    operation: int,
    phi,
    wall_mask: bool,
    logical: bool,
    map_r,
    map_z,
    edge_val,
    cgs: bool,
    mks: bool,
):
    print("  reading real field")
    mesh = read_mesh(filename=filename, slice=slice_idx)
    with h5py.File(str(filename), "r") as h5:
        gname = time_name(slice_idx)
        if gname not in h5:
            raise KeyError(f"Missing group '{gname}'.")
        tg = h5[gname]
        if "fields" not in tg or not isinstance(tg["fields"], h5py.Group):
            raise KeyError(f"Missing fields group in '{gname}'.")
        match = _find_field_name_case_insensitive(tg["fields"], field_name)
        if match is None:
            raise KeyError(field_name)
        coeff = _read_field_dataset(tg, match)

    itor = int(read_parameter("itor", filename=filename))
    phi_rad = float(phi)
    if itor == 1:
        phi_rad = float(phi) * np.pi / 180.0

    ev = eval_field(
        coeff,
        mesh,
        points=points,
        operation=operation,
        filename=filename,
        xrange=xrange,
        yrange=yrange,
        phi=phi_rad,
        wall_mask=wall_mask,
    )

    symbol, d = field_data(match, itor=itor, filename=filename)

    data = convert_units(ev.data, d, filename=filename, cgs=cgs, mks=mks)
    r = convert_units(ev.r, dimensions(l0=1), filename=filename, cgs=cgs, mks=mks)
    z = convert_units(ev.z, dimensions(l0=1), filename=filename, cgs=cgs, mks=mks)
    out_mask = np.asarray(ev.mask)

    igeometry = int(read_parameter("igeometry", filename=filename))
    if igeometry > 0 and not logical:
        cache_key = (
            str(filename),
            int(points),
            tuple(np.asarray(data).shape),
            float(phi),
            bool(cgs),
            bool(mks),
        )

        mr_cache = getattr(_read_primitive, "_mr_cache", None)
        mz_cache = getattr(_read_primitive, "_mz_cache", None)
        r_cache = getattr(_read_primitive, "_r_cache", None)
        z_cache = getattr(_read_primitive, "_z_cache", None)
        map_cache_key = getattr(_read_primitive, "_map_cache_key", None)

        if map_r is None and map_z is None and map_cache_key == cache_key:
            mr = None if mr_cache is None else np.asarray(mr_cache, dtype=float)
            mz = None if mz_cache is None else np.asarray(mz_cache, dtype=float)
            if r_cache is not None:
                r = np.asarray(r_cache, dtype=float)
            if z_cache is not None:
                z = np.asarray(z_cache, dtype=float)
        else:
            mr = None if map_r is None else np.asarray(map_r, dtype=float)
            mz = None if map_z is None else np.asarray(map_z, dtype=float)
        if mr is not None and mr.ndim == 3:
            mr = mr[0]
        if mz is not None and mz.ndim == 3:
            mz = mz[0]

        if mr is None or mz is None or mr.shape != np.asarray(data).shape or mz.shape != np.asarray(data).shape:
            rst = read_field(
                "rst",
                timeslices=-1,
                filename=filename,
                points=points,
                phi=phi,
                logical=True,
                cgs=cgs,
                mks=mks,
                return_meta=True,
            )
            zst = read_field(
                "zst",
                timeslices=-1,
                filename=filename,
                points=points,
                phi=phi,
                logical=True,
                cgs=cgs,
                mks=mks,
                return_meta=True,
            )
            mr, mz, r, z, _ = create_map(
                np.asarray(rst.data),
                np.asarray(zst.data),
                r=np.asarray(r),
                z=np.asarray(z),
                mask=np.asarray(rst.mask) if rst.mask is not None else np.asarray(out_mask),
            )
            _read_primitive._mr_cache = np.asarray(mr, dtype=float)
            _read_primitive._mz_cache = np.asarray(mz, dtype=float)
            _read_primitive._r_cache = np.asarray(r, dtype=float)
            _read_primitive._z_cache = np.asarray(z, dtype=float)
            _read_primitive._map_cache_key = cache_key
        elif map_r is None and map_z is None:
            _read_primitive._mr_cache = np.asarray(mr, dtype=float)
            _read_primitive._mz_cache = np.asarray(mz, dtype=float)
            _read_primitive._r_cache = np.asarray(r, dtype=float)
            _read_primitive._z_cache = np.asarray(z, dtype=float)
            _read_primitive._map_cache_key = cache_key

        data, out_mask = map_field(np.asarray(data), mr, mz, mask=out_mask, outval=edge_val)

    units = parse_units(d, cgs=cgs, mks=mks)

    return FieldResult(
        data=np.asarray(data),
        symbol=symbol,
        units=units,
        dimensions=d,
        r=np.asarray(r),
        z=np.asarray(z),
        time=0.0,
        mask=np.asarray(out_mask),
    )


def read_field(
    name: str,
    timeslices=0,
    *,
    mesh=None,
    filename: str | Path | Sequence[str | Path] = "C1.h5",
    points: int = 200,
    mask=None,
    xrange=None,
    yrange=None,
    equilibrium: bool = False,
    h_symmetry: bool = False,
    v_symmetry: bool = False,
    diff: bool = False,
    operation: int = 1,
    complex: bool = False,
    fac: float | None = None,
    linear: bool = False,
    last: bool = False,
    average: bool = False,
    linfac: float = 1.0,
    dpsi=None,
    symbol=None,
    units=None,
    cgs: bool = False,
    mks: bool = False,
    real: bool = False,
    imaginary: bool = False,
    edge_val=None,
    phi: float = 0.0,
    time=None,
    abs: bool = False,
    phase: bool = False,
    dimensions_out=None,
    flux_average: bool = False,
    rvector: bool = False,
    zvector: bool = False,
    yvector: bool = False,
    taverage: bool = False,
    sum: bool = False,
    tpoints=None,
    logical: bool = False,
    map_r=None,
    map_z=None,
    is_nonlinear=None,
    outval=None,
    wall_mask: bool = False,
    return_meta: bool = False,
):
    """Python port of read_field.pro (core recursive behavior)."""
    del mesh, h_symmetry, v_symmetry, dpsi, dimensions_out, flux_average, is_nonlinear, outval

    if isinstance(filename, (list, tuple)):
        files = list(filename)
    else:
        files = [filename]
    if timeslices is None:
        timeslices = 0
    if isinstance(timeslices, (list, tuple, np.ndarray)):
        times = [int(t0) for t0 in np.asarray(timeslices).reshape(-1)]
    else:
        times = [int(timeslices)]

    # average/sum across files or timeslices
    if average or sum or (len(files) > 1 and not diff):
        nfiles = len(files)
        ntimes = len(times)
        if nfiles > 1:
            n = nfiles
            if ntimes == 1:
                times = times * n
        elif ntimes > 1:
            n = ntimes
            files = files * n
        else:
            n = 1

        fac_vals = np.asarray([1.0] if fac is None else np.asarray(fac, dtype=float).reshape(-1), dtype=float)
        if fac_vals.size < n:
            fac_vals = np.repeat(fac_vals[:1], n)
        linfac_vals = np.asarray(np.asarray(linfac, dtype=float).reshape(-1), dtype=float)
        if linfac_vals.size < n:
            linfac_vals = np.repeat(linfac_vals[:1], n)

        vals = []
        for i in range(n):
            r = read_field(
                name,
                timeslices=times[i],
                filename=files[i],
                points=points,
                xrange=xrange,
                yrange=yrange,
                equilibrium=equilibrium,
                operation=operation,
                complex=complex,
                fac=float(fac_vals[i]),
                linear=linear,
                last=last,
                linfac=float(linfac_vals[i]),
                cgs=cgs,
                mks=mks,
                phi=phi,
                wall_mask=wall_mask,
                logical=logical,
                map_r=map_r,
                map_z=map_z,
                edge_val=edge_val,
                return_meta=True,
            )
            vals.append(r)
        base = vals[0]
        stack = np.stack([np.asarray(v.data) for v in vals], axis=0)
        data = np.sum(stack, axis=0)
        if average:
            data = data / max(len(vals), 1)
        out = FieldResult(
            data=data,
            symbol=base.symbol,
            units=base.units,
            dimensions=base.dimensions,
            r=base.r,
            z=base.z,
            time=base.time,
            mask=base.mask,
        )
        return out if return_meta else out.data

    # diff across files or timeslices (alternating signs like IDL)
    if diff:
        nfiles = len(files)
        ntimes = len(times)
        if nfiles > 1:
            n = nfiles
            if ntimes == 1:
                times = times * n
        elif ntimes > 1:
            n = ntimes
            files = files * n
        else:
            n = 1

        vals = []
        for i in range(n):
            vals.append(
                read_field(
                    name,
                    timeslices=times[i],
                    filename=files[i],
                    points=points,
                    xrange=xrange,
                    yrange=yrange,
                    equilibrium=equilibrium,
                    operation=operation,
                    complex=complex,
                    fac=fac,
                    linear=linear,
                    last=last,
                    linfac=linfac,
                    cgs=cgs,
                    mks=mks,
                    phi=phi,
                    wall_mask=wall_mask,
                    logical=logical,
                    map_r=map_r,
                    map_z=map_z,
                    edge_val=edge_val,
                    return_meta=True,
                )
            )
        base = vals[0]
        data = np.zeros_like(np.asarray(base.data))
        for i, v in enumerate(vals):
            sgn = -1.0 if (i % 2 == 1) else 1.0
            data = data + sgn * np.asarray(v.data)
        out = FieldResult(
            data=data,
            symbol=base.symbol,
            units=base.units,
            dimensions=base.dimensions,
            r=base.r,
            z=base.z,
            time=base.time,
            mask=base.mask,
        )
        return out if return_meta else out.data

    filename = files[0]
    if not hdf5_file_test(filename):
        raise ValueError(f"Invalid HDF5 file: {filename}")

    nt = int(read_parameter("ntime", filename=filename))
    itor = int(read_parameter("itor", filename=filename))
    ntor = int(read_parameter("ntor", filename=filename))
    icomplex = int(read_parameter("icomplex", filename=filename))
    isubeq = int(read_parameter("eqsubtract", filename=filename))
    ilin = int(read_parameter("linear", filename=filename))
    if complex and icomplex == 0:
        complex = False

    slice_idx = int(timeslices)
    if last:
        slice_idx = nt - 1
    if equilibrium:
        if isubeq == 1:
            slice_idx = -1
        if ilin == 1:
            slice_idx = -1
        if isubeq == 1:
            linear = False

    if slice_idx >= nt:
        raise ValueError(f"There are only {nt} time slices (0..{nt-1}).")

    print("**********************************************************")
    print(f"Reading {name} at timeslice {slice_idx}")
    print(f"From file {filename}")
    print(f"Eqsubtract? {isubeq}")
    print(f"sum = {bool(sum)}")
    if fac is not None:
        print(f"fac = {fac}")
    print(
        f" linear={int(bool(linear))}; pts={int(points)};"
        f"equilibrium={int(bool(equilibrium))}; complex={int(bool(complex))}; op={int(operation)}"
    )

    phi_arr = np.asarray(phi, dtype=float).reshape(-1)
    period = 2.0 * np.pi if itor == 1 else 2.0 * np.pi * float(read_parameter("rzero", filename=filename))
    if tpoints is None:
        nphi = int(phi_arr.size)
    else:
        nphi = int(tpoints)
    if phi_arr.size == nphi:
        phis = phi_arr.astype(float, copy=False)
    else:
        phis = np.linspace(0.0, period, nphi, endpoint=False)
        if itor == 1:
            phis = phis * 180.0 / np.pi
    print("phi0 = ", phis)

    if nphi > 1:
        if complex:
            base = read_field(
                name,
                timeslices=slice_idx,
                filename=filename,
                points=points,
                xrange=xrange,
                yrange=yrange,
                equilibrium=equilibrium,
                operation=operation,
                complex=True,
                linear=True,
                cgs=cgs,
                mks=mks,
                phi=float(phis[0]),
                wall_mask=wall_mask,
                logical=logical,
                map_r=map_r,
                map_z=map_z,
                edge_val=edge_val,
                return_meta=True,
            )
            base_data = np.asarray(base.data, dtype=np.complex128)
            stacks = np.zeros((nphi,) + base_data.shape, dtype=float)
            stacks[0, ...] = np.real(base_data)
            if itor == 1:
                phase_vals = (np.asarray(phis, dtype=float) - float(phis[0])) * np.pi / 180.0
            else:
                phase_vals = np.asarray(phis, dtype=float) - float(phis[0])
            for i in range(1, nphi):
                stacks[i, ...] = np.real(base_data * np.exp(1j * ntor * phase_vals[i]))
            if (not linear) and slice_idx >= 0:
                eq0 = read_field(
                    name,
                    timeslices=-1,
                    filename=filename,
                    points=points,
                    xrange=xrange,
                    yrange=yrange,
                    operation=operation,
                    complex=False,
                    fac=fac,
                    cgs=cgs,
                    mks=mks,
                    phi=0.0,
                    wall_mask=wall_mask,
                    logical=logical,
                    map_r=map_r,
                    map_z=map_z,
                    edge_val=edge_val,
                    return_meta=True,
                )
                eq_data = np.asarray(eq0.data, dtype=float)
                stacks[:, ...] = stacks[:, ...] + eq_data[None, ...]
            meta = base
        else:
            stacks_list = []
            meta = None
            for p0 in phis:
                r0 = read_field(
                    name,
                    timeslices=slice_idx,
                    filename=filename,
                    points=points,
                    xrange=xrange,
                    yrange=yrange,
                    equilibrium=equilibrium,
                    operation=operation,
                    complex=False,
                    linear=linear,
                    cgs=cgs,
                    mks=mks,
                    phi=float(p0),
                    wall_mask=wall_mask,
                    logical=logical,
                    map_r=map_r,
                    map_z=map_z,
                    edge_val=edge_val,
                    return_meta=True,
                )
                meta = r0
                stacks_list.append(np.asarray(r0.data))
            assert meta is not None
            stacks = np.stack(stacks_list, axis=0)
        assert meta is not None
        out = FieldResult(
            data=stacks,
            symbol=meta.symbol,
            units=meta.units,
            dimensions=meta.dimensions,
            r=meta.r,
            z=meta.z,
            time=meta.time,
            mask=meta.mask,
        )
        return out if return_meta else out.data

    if phi_arr.size > 0:
        phi = float(phi_arr[0])

    if taverage:
        nphi = int(taverage if isinstance(taverage, int) else 16)
        if nphi <= 0:
            nphi = 16
        period = 2.0 * np.pi if itor == 1 else 2.0 * np.pi * float(read_parameter("rzero", filename=filename))
        phis = np.linspace(0.0, period, nphi, endpoint=False)
        if itor == 1:
            phis = phis * 180.0 / np.pi
        acc = None
        meta = None
        for p0 in phis:
            r0 = read_field(
                name,
                timeslices=slice_idx,
                filename=filename,
                points=points,
                xrange=xrange,
                yrange=yrange,
                equilibrium=equilibrium,
                operation=operation,
                complex=complex,
                linear=linear,
                cgs=cgs,
                mks=mks,
                phi=float(p0),
                wall_mask=wall_mask,
                logical=logical,
                map_r=map_r,
                map_z=map_z,
                edge_val=edge_val,
                return_meta=True,
            )
            meta = r0
            acc = np.asarray(r0.data) if acc is None else acc + np.asarray(r0.data)
        assert meta is not None and acc is not None
        out = FieldResult(
            data=acc / float(nphi),
            symbol=meta.symbol,
            units=meta.units,
            dimensions=meta.dimensions,
            r=meta.r,
            z=meta.z,
            time=meta.time,
            mask=meta.mask,
        )
        return out if return_meta else out.data

    # IDL block (around line 227): for eqsubtract=1 non-linear reads, build total
    # field from equilibrium plus perturbed part, then return immediately.
    if isubeq == 1 and (not linear) and slice_idx >= 0:
        data1 = None
        try:
            data1 = read_field(
                name,
                timeslices=slice_idx,
                filename=filename,
                points=points,
                xrange=xrange,
                yrange=yrange,
                equilibrium=equilibrium,
                operation=operation,
                complex=complex,
                fac=fac,
                linear=True,
                last=last,
                linfac=linfac,
                cgs=cgs,
                mks=mks,
                phi=phi,
                wall_mask=wall_mask,
                logical=logical,
                map_r=map_r,
                map_z=map_z,
                edge_val=edge_val,
                return_meta=True,
            )
        except KeyError:
            print("Perturbed field not found.")

        data0 = read_field(
            name,
            timeslices=-1,
            filename=filename,
            points=points,
            xrange=xrange,
            yrange=yrange,
            operation=operation,
            complex=False,
            fac=fac,
            cgs=cgs,
            mks=mks,
            phi=phi,
            wall_mask=wall_mask,
            logical=logical,
            map_r=map_r,
            map_z=map_z,
            edge_val=edge_val,
            return_meta=True,
        )

        if data1 is None:
            out = data0
        else:
            out = FieldResult(
                data=np.asarray(data0.data) + np.asarray(data1.data),
                symbol=data0.symbol,
                units=data0.units,
                dimensions=data0.dimensions,
                r=data0.r,
                z=data0.z,
                time=data0.time,
                mask=data0.mask,
            )
        print("**********************************************************")
        return out if return_meta else out.data

    if linear and isubeq == 0 and slice_idx >= 0 and not complex:
        print(f"  calculating perturbed part of eqsubtract=0 field. {ntor}")
        f1 = read_field(
            name,
            timeslices=slice_idx,
            filename=filename,
            points=points,
            xrange=xrange,
            yrange=yrange,
            operation=operation,
            cgs=cgs,
            mks=mks,
            phi=phi,
            wall_mask=wall_mask,
            logical=logical,
            map_r=map_r,
            map_z=map_z,
            edge_val=edge_val,
            return_meta=True,
        )
        f0 = read_field(
            name,
            timeslices=-1,
            filename=filename,
            points=points,
            xrange=xrange,
            yrange=yrange,
            operation=operation,
            cgs=cgs,
            mks=mks,
            phi=phi,
            wall_mask=wall_mask,
            logical=logical,
            map_r=map_r,
            map_z=map_z,
            edge_val=edge_val,
            return_meta=True,
        )
        out = FieldResult(
            data=np.asarray(f1.data) - np.asarray(f0.data),
            symbol=f1.symbol,
            units=f1.units,
            dimensions=f1.dimensions,
            r=f1.r,
            z=f1.z,
            time=f1.time,
            mask=f1.mask,
        )
        return out if return_meta else out.data

    if complex:
        print(f"  reading complex field. {ntor}")
        fr = read_field(
            name,
            timeslices=slice_idx,
            filename=filename,
            points=points,
            xrange=xrange,
            yrange=yrange,
            operation=operation,
            linear=True,
            cgs=cgs,
            mks=mks,
            phi=phi,
            wall_mask=wall_mask,
            return_meta=True,
        )
        fi = read_field(
            f"{name}_i",
            timeslices=slice_idx,
            filename=filename,
            points=points,
            xrange=xrange,
            yrange=yrange,
            operation=operation,
            linear=True,
            cgs=cgs,
            mks=mks,
            phi=phi,
            wall_mask=wall_mask,
            return_meta=True,
        )
        data = np.asarray(fr.data) + 1j * np.asarray(fi.data)
        if phi is not None:
            ph = float(phi)
            ph = ph * np.pi / 180.0 if itor == 1 else ph
            data = data * np.exp(1j * ntor * ph)
        out = FieldResult(
            data=data,
            symbol=fr.symbol,
            units=fr.units,
            dimensions=fr.dimensions,
            r=fr.r,
            z=fr.z,
            time=fr.time,
            mask=fr.mask,
        )
    else:
        # primitive path
        primitive = None
        primitive_error: KeyError | None = None
        try:
            primitive = _read_primitive(
                name,
                filename=filename,
                slice_idx=slice_idx,
                points=points,
                xrange=xrange,
                yrange=yrange,
                operation=operation,
                phi=phi,
                wall_mask=wall_mask,
                logical=logical,
                map_r=map_r,
                map_z=map_z,
                edge_val=edge_val,
                cgs=cgs,
                mks=mks,
            )
        except KeyError as e:
            primitive = None
            primitive_error = e

        if primitive is not None:
            out = primitive
        else:
            print("  reading composite field")
            n = str(name).strip().lower()
            if n == "zero":
                base = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                out = FieldResult(
                    data=np.zeros_like(base.data),
                    symbol="",
                    units=base.units,
                    dimensions=base.dimensions,
                    r=base.r,
                    z=base.z,
                    time=base.time,
                    mask=base.mask,
                )
            elif n in {"fp", "fp_plasma", "fp_ext"}:
                if n == "fp":
                    print("Calculating fp from f")
                elif n == "fp_plasma":
                    print("Calculating fp_plasma from f_plasma")
                else:
                    print("Calculating fp_ext from f_ext")
                base_name = "f" if n == "fp" else ("f_plasma" if n == "fp_plasma" else "f_ext")
                fbase = read_field(
                    base_name,
                    filename=filename,
                    timeslices=slice_idx,
                    points=points,
                    xrange=xrange,
                    yrange=yrange,
                    linear=linear,
                    complex=complex,
                    phi=phi,
                    cgs=cgs,
                    mks=mks,
                    return_meta=True,
                )
                data = np.asarray(fbase.data) * (1j * ntor if ntor != 0 else 0.0)
                d = fbase.dimensions + dimensions(l0=-1)
                out = FieldResult(
                    data=data,
                    symbol="$f_p$",
                    units=parse_units(d, cgs=cgs, mks=mks),
                    dimensions=d,
                    r=fbase.r,
                    z=fbase.z,
                    time=fbase.time,
                    mask=fbase.mask,
                )
            elif n in {"grad_psi"}:
                psi_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(b0=1, l0=1 + itor)
                out = FieldResult(
                    data=np.sqrt(np.asarray(psi_r.data) ** 2 + np.asarray(psi_z.data) ** 2),
                    symbol="$|\\nabla\\psi|$",
                    units=parse_units(d, cgs=cgs, mks=mks),
                    dimensions=d,
                    r=psi_r.r,
                    z=psi_r.z,
                    time=psi_r.time,
                    mask=psi_r.mask,
                )
            elif n in {"psi_norm"}:
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                lcfs = read_lcfs(filename=filename, slice=slice_idx, cgs=cgs, mks=mks, return_meta=True)
                denom = lcfs.psilim - lcfs.flux0
                if np.abs(denom) < np.finfo(float).tiny:
                    denom = np.finfo(float).tiny
                out = FieldResult(
                    data=(np.asarray(psi.data) - lcfs.flux0) / denom,
                    symbol="$\\psi_N$",
                    units="",
                    dimensions=dimensions(),
                    r=psi.r,
                    z=psi.z,
                    time=psi.time,
                    mask=psi.mask,
                )
            elif n in {"grad_psi_norm"}:
                psi_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                lcfs = read_lcfs(filename=filename, slice=slice_idx, cgs=cgs, mks=mks, return_meta=True)
                denom = float(np.abs(lcfs.psilim - lcfs.flux0))
                if denom < np.finfo(float).tiny:
                    denom = np.finfo(float).tiny
                d = dimensions(l0=-1)
                out = FieldResult(
                    data=np.sqrt(np.asarray(psi_r.data) ** 2 + np.asarray(psi_z.data) ** 2) / denom,
                    symbol="$|\\nabla\\psi_N|$",
                    units=parse_units(d, cgs=cgs, mks=mks),
                    dimensions=d,
                    r=psi_r.r,
                    z=psi_r.z,
                    time=psi_r.time,
                    mask=psi_r.mask,
                )
            elif n in {"ion temperature", "ti"}:
                ti = read_field("ti", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                out = ti
            elif n in {"vt_i", "vti"}:
                ti = read_field("ti", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(v0=1)
                out = FieldResult(
                    data=np.sqrt(np.maximum(2.0 * np.asarray(ti.data), 0.0)),
                    symbol="$v_{ti}$",
                    units=parse_units(d, cgs=cgs, mks=mks),
                    dimensions=d,
                    r=ti.r,
                    z=ti.z,
                    time=ti.time,
                    mask=ti.mask,
                )
            elif n in {"sound speed", "cs"}:
                p = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                gam = float(read_parameter("gam", filename=filename))
                d = dimensions(v0=1)
                out = FieldResult(
                    data=np.sqrt(np.maximum(gam * np.asarray(p.data) / np.maximum(np.asarray(den.data), np.finfo(float).tiny), 0.0)),
                    symbol="$c_s$",
                    units=parse_units(d, cgs=cgs, mks=mks),
                    dimensions=d,
                    r=p.r,
                    z=p.z,
                    time=p.time,
                    mask=p.mask,
                )
            elif n in {"mach", "m"}:
                cs = read_field("cs", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                v_f = read_field("v", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                chi_f = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(phi_f.r, phi_f.z) if itor == 1 else np.ones_like(np.asarray(phi_f.data), dtype=float)
                ivform = int(read_parameter("ivform", filename=filename))
                if ivform == 0:
                    v2 = (
                        s_bracket(phi_f.data, phi_f.data, phi_f.r, phi_f.z) / (rmat**2)
                        + np.asarray(v_f.data) ** 2 / (rmat**2)
                        + s_bracket(chi_f.data, chi_f.data, chi_f.r, chi_f.z)
                        + 2.0 * a_bracket(chi_f.data, phi_f.data, phi_f.r, phi_f.z) / rmat
                    )
                else:
                    v2 = (
                        (rmat**2) * s_bracket(phi_f.data, phi_f.data, phi_f.r, phi_f.z)
                        + (rmat**2) * np.asarray(v_f.data) ** 2
                        + s_bracket(chi_f.data, chi_f.data, chi_f.r, chi_f.z) / np.maximum(rmat**4, np.finfo(float).tiny)
                        + 2.0 * a_bracket(chi_f.data, phi_f.data, phi_f.r, phi_f.z) / rmat
                    )
                out = FieldResult(
                    data=np.sqrt(np.maximum(v2, 0.0)) / np.maximum(np.asarray(cs.data), np.finfo(float).tiny),
                    symbol="$M$",
                    units="",
                    dimensions=dimensions(),
                    r=cs.r,
                    z=cs.z,
                    time=cs.time,
                    mask=cs.mask,
                )
            elif n in {"ion pressure", "pi"}:
                p = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                pe = read_field("pe", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(p0=1)
                out = FieldResult(
                    data=np.asarray(p.data) - np.asarray(pe.data),
                    symbol="$p_i$",
                    units=parse_units(d, cgs=cgs, mks=mks),
                    dimensions=d,
                    r=p.r,
                    z=p.z,
                    time=p.time,
                    mask=p.mask,
                )
            elif n in {"jy", "jphi", "jy_ext", "jy_plasma"}:
                psi_name = "psi" if n in {"jy", "jphi"} else ("psi_ext" if n == "jy_ext" else "psi_plasma")
                lp = read_field(psi_name, filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=7, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    psir = read_field(psi_name, filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    rmat = radius_matrix(lp.r, lp.z)
                    data = -(np.asarray(lp.data) - np.asarray(psir.data) / rmat) / rmat
                else:
                    data = -np.asarray(lp.data)
                d = dimensions(j0=1)
                out = FieldResult(data=data, symbol="$J_{phi}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=lp.r, z=lp.z, time=lp.time, mask=lp.mask)
            elif n in {"jx", "jx_ext", "jx_plasma"}:
                i_name = "i" if n == "jx" else ("i_ext" if n == "jx_ext" else "i_plasma")
                psi_name = "psi" if n == "jx" else ("psi_ext" if n == "jx_ext" else "psi_plasma")
                fp_name = "fp" if n == "jx" else ("fp_ext" if n == "jx_ext" else "f_plasma")
                iz = read_field(i_name, filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(iz.r, iz.z) if itor == 1 else np.ones_like(np.asarray(iz.data), dtype=float)
                data = -np.asarray(iz.data) / rmat
                if ntor != 0:
                    psi_r = read_field(psi_name, filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    fp_z = read_field(fp_name, filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    if n == "jx_plasma":
                        data = data + (ntor**2) * np.asarray(fp_z.data) / np.maximum(rmat, np.finfo(float).tiny) + 1j * ntor * np.asarray(psi_r.data) / np.maximum(rmat**2, np.finfo(float).tiny)
                    else:
                        data = data + 1j * ntor * (-np.asarray(fp_z.data) / rmat + np.asarray(psi_r.data) / np.maximum(rmat**2, np.finfo(float).tiny))
                d = dimensions(j0=1)
                out = FieldResult(data=data, symbol="$J_R$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=iz.r, z=iz.z, time=iz.time, mask=iz.mask)
            elif n in {"jz", "jz_ext", "jz_plasma"}:
                i_name = "i" if n == "jz" else ("i_ext" if n == "jz_ext" else "i_plasma")
                psi_name = "psi" if n == "jz" else ("psi_ext" if n == "jz_ext" else "psi_plasma")
                fp_name = "fp" if n == "jz" else ("fp_ext" if n == "jz_ext" else "f_plasma")
                ir = read_field(i_name, filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(ir.r, ir.z) if itor == 1 else np.ones_like(np.asarray(ir.data), dtype=float)
                data = np.asarray(ir.data) / rmat
                if ntor != 0:
                    psi_z = read_field(psi_name, filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    fp_r = read_field(fp_name, filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    if n == "jz_plasma":
                        data = data - (ntor**2) * np.asarray(fp_r.data) / np.maximum(rmat, np.finfo(float).tiny) + 1j * ntor * np.asarray(psi_z.data) / np.maximum(rmat**2, np.finfo(float).tiny)
                    else:
                        data = data + 1j * ntor * (np.asarray(fp_r.data) / rmat + np.asarray(psi_z.data) / np.maximum(rmat**2, np.finfo(float).tiny))
                d = dimensions(j0=1)
                out = FieldResult(data=data, symbol="$J_Z$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=ir.r, z=ir.z, time=ir.time, mask=ir.mask)
            elif n == "jp":
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, linear=linear, wall_mask=wall_mask, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, wall_mask=wall_mask, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(psi0.r, psi0.z) if itor == 1 else np.ones_like(np.asarray(i_f.data), dtype=float)
                data = s_bracket(i_f.data, psi0.data, psi0.r, psi0.z) / np.maximum(
                    rmat * np.sqrt(np.maximum(s_bracket(psi0.data, psi0.data, psi0.r, psi0.z), np.finfo(float).tiny)),
                    np.finfo(float).tiny,
                )
                d = dimensions(j0=1)
                out = FieldResult(data=data, symbol="$J_p$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi0.r, z=psi0.z, time=psi0.time, mask=psi0.mask)
            elif n == "etaj2":
                jy = read_field("jy", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, wall_mask=wall_mask, cgs=cgs, mks=mks, return_meta=True)
                eta = read_field("eta", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, wall_mask=wall_mask, cgs=cgs, mks=mks, return_meta=True)
                eta_data = np.asarray(eta.data)
                d = dimensions(p0=1, t0=-1)
                out = FieldResult(data=eta_data * np.asarray(jy.data) ** 2, symbol="$\\eta J_y^2$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=jy.r, z=jy.z, time=jy.time, mask=jy.mask)
            elif n == "vor":
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                data = delstar(np.asarray(phi_f.data), phi_f.r, phi_f.z, toroidal=(itor == 1))
                _, d = field_data("vor", itor=itor, filename=filename)
                out = FieldResult(data=data, symbol="$\\omega$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=phi_f.r, z=phi_f.z, time=phi_f.time, mask=phi_f.mask)
            elif n == "com":
                chi_f = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                data = laplacian(np.asarray(chi_f.data), chi_f.r, chi_f.z, toroidal=(itor == 1))
                _, d = field_data("com", itor=itor, filename=filename)
                out = FieldResult(data=data, symbol="$\\nabla\\cdot v$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=chi_f.r, z=chi_f.z, time=chi_f.time, mask=chi_f.mask)
            elif n == "r2vor":
                phi_lp = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, operation=7, cgs=cgs, mks=mks, return_meta=True)
                phi_r = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, operation=2, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(phi_lp.r, phi_lp.z)
                data = np.maximum(rmat**2, np.finfo(float).tiny) * np.asarray(phi_lp.data) + 2.0 * rmat * np.asarray(phi_r.data)
                _, d = field_data("vor", itor=itor, filename=filename)
                out = FieldResult(data=data, symbol="$R^2\\omega$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=phi_lp.r, z=phi_lp.z, time=phi_lp.time, mask=phi_lp.mask)
            elif n == "helicity":
                psi_f = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                f_f = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(psi_f.r, psi_f.z)
                data = (
                    (np.asarray(psi_f.data) * np.conj(np.asarray(i_f.data)) + np.conj(np.asarray(psi_f.data)) * np.asarray(i_f.data)) / np.maximum(rmat**4, np.finfo(float).tiny)
                    - s_bracket(np.asarray(f_f.data), np.conj(np.asarray(psi_f.data)), psi_f.r, psi_f.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                    - s_bracket(np.conj(np.asarray(f_f.data)), np.asarray(psi_f.data), psi_f.r, psi_f.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                )
                d = dimensions(l0=4, b0=2)
                out = FieldResult(data=data, symbol="$H$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi_f.r, z=psi_f.z, time=psi_f.time, mask=psi_f.mask)
            elif n in {"electron temperature", "te", "te_i"}:
                pe0 = read_field("pe", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                n0 = read_field("ne", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                if linear and isubeq == 1 and slice_idx >= 0:
                    pe1 = read_field("pe", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    n1 = read_field("ne", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    data = np.asarray(pe1.data) / np.maximum(np.asarray(n0.data), np.finfo(float).tiny) - np.asarray(pe0.data) * np.asarray(n1.data) / np.maximum(np.asarray(n0.data) ** 2, np.finfo(float).tiny)
                else:
                    data = np.asarray(pe0.data) / np.maximum(np.asarray(n0.data), np.finfo(float).tiny)
                d = dimensions(temperature=1)
                out = FieldResult(data=data, symbol="$T_e$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=pe0.r, z=pe0.z, time=pe0.time, mask=pe0.mask)
            elif n == "beta":
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                p_f = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(psi.r, psi.z) if itor == 1 else np.ones_like(np.asarray(psi.data), dtype=float)
                b2 = (s_bracket(psi.data, psi.data, psi.r, psi.z) + np.asarray(i_f.data) ** 2) / np.maximum(rmat**2, np.finfo(float).tiny)
                out = FieldResult(data=2.0 * np.asarray(p_f.data) / np.maximum(b2, np.finfo(float).tiny), symbol="$\\beta$", units="", dimensions=dimensions(), r=psi.r, z=psi.z, time=psi.time, mask=psi.mask)
            elif n == "omega":
                if operation != 1:
                    print("Warning: using op on omega")
                v = read_field("v", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=operation, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ivform = int(read_parameter("ivform", filename=filename))
                if ivform == 0:
                    rmat = radius_matrix(v.r, v.z) if itor == 1 else np.ones_like(np.asarray(v.data), dtype=float)
                    data = np.asarray(v.data) / np.maximum(rmat**2, np.finfo(float).tiny)
                else:
                    data = np.asarray(v.data)
                if itor == 0:
                    rzero = float(read_parameter("rzero", filename=filename))
                    data = data / max(rzero, np.finfo(float).tiny)
                d = dimensions(t0=-1)
                out = FieldResult(data=data, symbol="$\\Omega$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=v.r, z=v.z, time=v.time, mask=v.mask)
            elif n == "pprime":
                p = read_field(
                    "p",
                    filename=filename,
                    timeslices=slice_idx,
                    points=points,
                    xrange=xrange,
                    yrange=yrange,
                    linear=linear,
                    complex=complex,
                    operation=operation,
                    phi=phi,
                    cgs=cgs,
                    mks=mks,
                    return_meta=True,
                )
                psi = read_field(
                    "psi",
                    filename=filename,
                    timeslices=slice_idx,
                    points=points,
                    xrange=xrange,
                    yrange=yrange,
                    equilibrium=True,
                    linear=linear,
                    operation=operation,
                    phi=phi,
                    cgs=cgs,
                    mks=mks,
                    return_meta=True,
                )
                num = s_bracket(np.asarray(p.data), np.asarray(psi.data), p.r, p.z)
                den = s_bracket(np.asarray(psi.data), np.asarray(psi.data), p.r, p.z)
                d = dimensions(p0=1, b0=1, l0=1 + itor)
                out = FieldResult(
                    data=num / np.maximum(den, np.finfo(float).tiny),
                    symbol="$p'$",
                    units=parse_units(d, cgs=cgs, mks=mks),
                    dimensions=d,
                    r=p.r,
                    z=p.z,
                    time=p.time,
                    mask=p.mask,
                )
            elif n == "ke":
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                v_f = read_field("v", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(den.r, den.z)
                else:
                    rmat = np.ones_like(np.asarray(den.data), dtype=float)
                ivform = int(read_parameter("ivform", filename=filename))
                if ivform == 0:
                    v2 = (
                        s_bracket(np.asarray(phi_f.data), np.asarray(phi_f.data), den.r, den.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                        + np.asarray(v_f.data) ** 2 / np.maximum(rmat**2, np.finfo(float).tiny)
                        + s_bracket(np.asarray(chi.data), np.asarray(chi.data), den.r, den.z)
                        + 2.0 * a_bracket(np.asarray(chi.data), np.asarray(phi_f.data), den.r, den.z) / np.maximum(rmat, np.finfo(float).tiny)
                    )
                else:
                    v2 = (
                        s_bracket(np.asarray(phi_f.data), np.asarray(phi_f.data), den.r, den.z)
                        + np.maximum(rmat**2, np.finfo(float).tiny) * np.asarray(v_f.data) ** 2
                        + s_bracket(np.asarray(chi.data), np.asarray(chi.data), den.r, den.z) / np.maximum(rmat**4, np.finfo(float).tiny)
                        + 2.0 * a_bracket(np.asarray(chi.data), np.asarray(phi_f.data), den.r, den.z) / np.maximum(rmat, np.finfo(float).tiny)
                    )
                d = dimensions(p0=1)
                out = FieldResult(data=0.5 * np.asarray(den.data) * v2, symbol="$Kinetic\\ Energy\\ Density$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=den.r, z=den.z, time=den.time, mask=den.mask)
            elif n == "va":
                b = read_field("b", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(v0=1)
                out = FieldResult(data=np.asarray(b.data) / np.sqrt(np.maximum(np.asarray(den.data), np.finfo(float).tiny)), symbol="$v_A$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=b.r, z=b.z, time=b.time, mask=b.mask)
            elif n == "dreicer":
                den_mks = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, mks=True, cgs=False, return_meta=True)
                te_mks = read_field("te", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, mks=True, cgs=False, return_meta=True)
                nsi = np.maximum(np.asarray(den_mks.data), np.finfo(float).tiny)
                te = np.maximum(np.asarray(te_mks.data), np.finfo(float).tiny)
                lnL = 14.9 - 0.5 * np.log(np.maximum(nsi / 1.0e20, 1.0e-300)) + np.log(np.maximum(te / 1.0e3, 1.0e-300))
                me_rest = 511875.0
                norm = 5.0244e-23
                ec = norm * nsi * lnL
                data_mks = ec * me_rest / te
                d = dimensions(potential=1, l0=-1)
                e_norm = convert_units(np.array(1.0), d, filename=filename, mks=True, cgs=False)
                out = FieldResult(data=data_mks / max(float(np.asarray(e_norm).reshape(-1)[0]), np.finfo(float).tiny), symbol="$E_{Dreicer}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=den_mks.r, z=den_mks.z, time=den_mks.time, mask=den_mks.mask)
            elif n == "delta_w" or n == "delta_w_i":
                pert = "_i" if n.endswith("_i") else ""
                jr0 = read_field("jx", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                jphi0 = read_field("jy", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                jz0 = read_field("jz", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                jr1 = read_field(f"jx{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                jphi1 = read_field(f"jy{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                jz1 = read_field(f"jz{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                br0 = read_field("bx", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                bphi0 = read_field("by", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                bz0 = read_field("bz", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                br1 = read_field(f"bx{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bphi1 = read_field(f"by{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bz1 = read_field(f"bz{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                p1 = read_field(f"p{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                pr1 = read_field(f"p{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, operation=2, cgs=cgs, mks=mks, return_meta=True)
                pz1 = read_field(f"p{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, operation=3, cgs=cgs, mks=mks, return_meta=True)
                xir = read_field(f"xi_x{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                xiphi = read_field(f"xi_y{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                xiz = read_field(f"xi_z{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(jr0.r, jr0.z) if itor == 1 else np.ones_like(np.asarray(jr0.data), dtype=float)
                pphi1 = 1j * ntor * np.asarray(p1.data)
                fr = np.asarray(jphi0.data) * np.asarray(bz1.data) - np.asarray(jz0.data) * np.asarray(bphi1.data) + np.asarray(jphi1.data) * np.asarray(bz0.data) - np.asarray(jz1.data) * np.asarray(bphi0.data) - np.asarray(pr1.data)
                fphi = np.asarray(jz0.data) * np.asarray(br1.data) - np.asarray(jr0.data) * np.asarray(bz1.data) + np.asarray(jz1.data) * np.asarray(br0.data) - np.asarray(jr1.data) * np.asarray(bz0.data) - pphi1 / np.maximum(rmat, np.finfo(float).tiny)
                fz = np.asarray(jr0.data) * np.asarray(bphi1.data) - np.asarray(jphi0.data) * np.asarray(br1.data) + np.asarray(jr1.data) * np.asarray(bphi0.data) - np.asarray(jphi1.data) * np.asarray(br0.data) - np.asarray(pz1.data)
                data = -0.5 * (np.conj(np.asarray(xir.data)) * fr + np.conj(np.asarray(xiphi.data)) * fphi + np.conj(np.asarray(xiz.data)) * fz)
                d = dimensions(p0=1)
                out = FieldResult(data=data, symbol="$\\delta W$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=jr0.r, z=jr0.z, time=jr0.time, mask=jr0.mask)
            elif n == "xi_n" or n == "xi_n_i":
                print("xi_n = normal displacement using xi = int(v dt) method")
                pert = "_i" if n.endswith("_i") else ""
                g = np.asarray(read_gamma(filename=filename), dtype=float).reshape(-1)
                gamma = float(g[0]) if g.size else 0.0
                vn = read_field(f"vn{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                data = np.asarray(vn.data) if np.isclose(gamma, 0.0) else np.asarray(vn.data) / gamma
                d = dimensions(l0=1)
                out = FieldResult(data=data, symbol="$\\xi_n$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=vn.r, z=vn.z, time=vn.time, mask=vn.mask)
            elif n == "xi_x" or n == "xi_x_i":
                print("xi_x = normal displacement using xi = int(v dt) method")
                pert = "_i" if n.endswith("_i") else ""
                vx = read_field(f"vx{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                g = np.asarray(read_gamma(filename=filename), dtype=float).reshape(-1)
                gamma = float(g[0]) if g.size else 0.0
                data = np.asarray(vx.data) / (gamma if not np.isclose(gamma, 0.0) else 1.0)
                d = dimensions(l0=1)
                out = FieldResult(data=data, symbol="$\\xi_R$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=vx.r, z=vx.z, time=vx.time, mask=vx.mask)
            elif n == "xi_y" or n == "xi_y_i":
                print("xi_y = toroidal displacement using xi = int(v dt) method")
                pert = "_i" if n.endswith("_i") else ""
                vy = read_field(f"vy{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                g = np.asarray(read_gamma(filename=filename), dtype=float).reshape(-1)
                gamma = float(g[0]) if g.size else 0.0
                data = np.asarray(vy.data) / (gamma if not np.isclose(gamma, 0.0) else 1.0)
                d = dimensions(l0=1)
                out = FieldResult(data=data, symbol="$\\xi_{\\phi}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=vy.r, z=vy.z, time=vy.time, mask=vy.mask)
            elif n == "xi_z" or n == "xi_z_i":
                print("xi_Z = verical displacement using xi = int(v dt) method")
                pert = "_i" if n.endswith("_i") else ""
                vz = read_field(f"vz{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                g = np.asarray(read_gamma(filename=filename), dtype=float).reshape(-1)
                gamma = float(g[0]) if g.size else 0.0
                data = np.asarray(vz.data) / (gamma if not np.isclose(gamma, 0.0) else 1.0)
                d = dimensions(l0=1)
                out = FieldResult(data=data, symbol="$\\xi_Z$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=vz.r, z=vz.z, time=vz.time, mask=vz.mask)
            elif n == "displacement" or n == "displacement_i":
                pert = "_i" if n.endswith("_i") else ""
                te1 = read_field(f"te{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                te0_r = read_field("te", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, operation=2, cgs=cgs, mks=mks, return_meta=True)
                te0_z = read_field("te", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, operation=3, cgs=cgs, mks=mks, return_meta=True)
                psi0 = read_field("psi", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                psi0_r = read_field("psi", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, operation=2, cgs=cgs, mks=mks, return_meta=True)
                psi0_z = read_field("psi", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, operation=3, cgs=cgs, mks=mks, return_meta=True)
                tprime = np.asarray(te0_r.data) * np.asarray(psi0_r.data) + np.asarray(te0_z.data) * np.asarray(psi0_z.data)
                gradpsi = np.sqrt(np.maximum(np.asarray(psi0_r.data) ** 2 + np.asarray(psi0_z.data) ** 2, np.finfo(float).tiny))
                data = -np.asarray(te1.data) / np.maximum(tprime, np.finfo(float).tiny) * gradpsi
                lc = read_lcfs(filename=filename, slice=slice_idx, cgs=cgs, mks=mks, return_meta=True)
                psis = float(lc.psilim)
                flux0 = float(lc.flux0)
                if psis < flux0:
                    data = np.where((np.abs(tprime) < 1e-6) | (np.asarray(psi0.data) < psis), 0.0, data)
                else:
                    data = np.where((np.abs(tprime) < 1e-6) | (np.asarray(psi0.data) > psis), 0.0, data)
                d = dimensions(l0=1)
                out = FieldResult(data=data, symbol="$\\xi$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=te1.r, z=te1.z, time=te1.time, mask=te1.mask)
            elif n == "overlap" or n == "overlap_i":
                pert = "_i" if n.endswith("_i") else ""
                xi = read_field(f"displacement{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi0 = read_field("psi", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                val = s_bracket(np.asarray(xi.data), np.asarray(psi0.data), xi.r, xi.z) / np.sqrt(np.maximum(s_bracket(np.asarray(psi0.data), np.asarray(psi0.data), xi.r, xi.z), np.finfo(float).tiny))
                data = np.abs(val) if n == "overlap" else np.imag(val)
                d = dimensions()
                out = FieldResult(data=data, symbol="$|d\\xi_r/dr|$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=xi.r, z=xi.z, time=xi.time, mask=xi.mask)
            elif n == "linearity" or n == "linearity_i":
                pert = "_i" if n.endswith("_i") else ""
                xi = read_field(f"displacement{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi0 = read_field("psi", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                p0 = read_field("p", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                lp = np.asarray(p0.data) / np.maximum(s_bracket(np.asarray(p0.data), np.asarray(psi0.data), xi.r, xi.z), np.finfo(float).tiny) * np.sqrt(np.maximum(s_bracket(np.asarray(psi0.data), np.asarray(psi0.data), xi.r, xi.z), np.finfo(float).tiny))
                data = np.asarray(xi.data) / np.maximum(lp, np.finfo(float).tiny)
                d = dimensions()
                out = FieldResult(data=data, symbol="$|\\xi_r|/L_p$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=xi.r, z=xi.z, time=xi.time, mask=xi.mask)
            elif n == "vp":
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(phi_f.r, phi_f.z)
                else:
                    rmat = np.ones_like(np.asarray(phi_f.data), dtype=float)
                psipsi = s_bracket(np.asarray(psi.data), np.asarray(psi.data), phi_f.r, phi_f.z)
                data = (
                    s_bracket(np.asarray(phi_f.data), np.asarray(psi.data), phi_f.r, phi_f.z) / np.maximum(rmat, np.finfo(float).tiny)
                    + a_bracket(np.asarray(chi.data), np.asarray(psi.data), phi_f.r, phi_f.z)
                ) / np.sqrt(np.maximum(psipsi, np.finfo(float).tiny))
                d = dimensions(v0=1)
                out = FieldResult(data=data, symbol="$u_p$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=phi_f.r, z=phi_f.z, time=phi_f.time, mask=phi_f.mask)
            elif n == "vs":
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                v_f = read_field("v", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(phi_f.r, phi_f.z)
                else:
                    rmat = np.ones_like(np.asarray(phi_f.data), dtype=float)
                psipsi = s_bracket(np.asarray(psi.data), np.asarray(psi.data), phi_f.r, phi_f.z)
                ivform = int(read_parameter("ivform", filename=filename))
                if ivform == 0:
                    b2 = psipsi + np.asarray(i_f.data) ** 2
                    data = -(
                        np.asarray(i_f.data) * s_bracket(np.asarray(phi_f.data), np.asarray(psi.data), phi_f.r, phi_f.z)
                        - np.asarray(v_f.data) * psipsi
                        + np.asarray(i_f.data) * a_bracket(np.asarray(chi.data), np.asarray(psi.data), phi_f.r, phi_f.z) * rmat
                    ) / (np.maximum(rmat**2, np.finfo(float).tiny) * np.sqrt(np.maximum(psipsi * b2, np.finfo(float).tiny)))
                else:
                    b2 = (psipsi + np.asarray(i_f.data) ** 2) / np.maximum(rmat**2, np.finfo(float).tiny)
                    data = -(
                        np.asarray(i_f.data) * s_bracket(np.asarray(phi_f.data), np.asarray(psi.data), phi_f.r, phi_f.z)
                        - np.asarray(v_f.data) * psipsi
                        + np.asarray(i_f.data) * a_bracket(np.asarray(chi.data), np.asarray(psi.data), phi_f.r, phi_f.z) / np.maximum(rmat**3, np.finfo(float).tiny)
                    ) / np.sqrt(np.maximum(b2 * psipsi, np.finfo(float).tiny))
                d = dimensions(v0=1)
                out = FieldResult(data=data, symbol="$u_s$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=phi_f.r, z=phi_f.z, time=phi_f.time, mask=phi_f.mask)
            elif n == "omega_ci":
                db = float(read_parameter("db", filename=filename))
                if db == 0.0:
                    db = 1.0
                b = read_field("b", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(t0=-1)
                out = FieldResult(data=np.asarray(b.data) / db, symbol="$\\omega_{ci}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=b.r, z=b.z, time=b.time, mask=b.mask)
            elif n in {"omega_*", "omega_*i", "omega_*e"}:
                db = float(read_parameter("db", filename=filename))
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                if n == "omega_*":
                    p = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    num = s_bracket(p.data, psi.data, psi.r, psi.z)
                elif n == "omega_*i":
                    p = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    pe = read_field("pe", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    num = s_bracket(np.asarray(p.data) - np.asarray(pe.data), psi.data, psi.r, psi.z)
                else:
                    pe = read_field("pe", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    num = -s_bracket(pe.data, psi.data, psi.r, psi.z)
                deno = np.maximum(s_bracket(psi.data, psi.data, psi.r, psi.z) * np.maximum(np.asarray(den.data), np.finfo(float).tiny), np.finfo(float).tiny)
                d = dimensions(t0=-1)
                sym = "$\\omega_*$" if n == "omega_*" else ("$\\omega_{*i}$" if n == "omega_*i" else "$\\omega_{*e}$")
                out = FieldResult(data=db * num / deno, symbol=sym, units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi.r, z=psi.z, time=psi.time, mask=psi.mask)
            elif n == "omega_exb":
                om = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                omi = read_field("omega_*i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(t0=-1)
                out = FieldResult(data=np.asarray(om.data) - np.asarray(omi.data), symbol="$\\omega_{ExB}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=om.r, z=om.z, time=om.time, mask=om.mask)
            elif n == "v_star":
                rho_i = read_field("rho_i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                omega_ci = read_field("omega_ci", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                p = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                pabs = np.maximum(np.abs(np.asarray(p.data)), np.finfo(float).tiny)
                kr = np.sqrt(np.maximum(s_bracket(np.asarray(p.data), np.asarray(p.data), p.r, p.z), 0.0)) / pabs
                d = dimensions(v0=1)
                out = FieldResult(data=np.asarray(rho_i.data) ** 2 * np.asarray(omega_ci.data) * kr, symbol="$u_*$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=p.r, z=p.z, time=p.time, mask=p.mask)
            elif n in {"omega_star", "omega_*"}:
                vstar = read_field("v_star", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(vstar.r, vstar.z)
                else:
                    rmat = np.ones_like(np.asarray(vstar.data), dtype=float)
                d = dimensions(t0=-1)
                out = FieldResult(
                    data=np.asarray(vstar.data) / np.maximum(rmat, np.finfo(float).tiny),
                    symbol="$\\omega_*$",
                    units=parse_units(d, cgs=cgs, mks=mks),
                    dimensions=d,
                    r=vstar.r,
                    z=vstar.z,
                    time=vstar.time,
                    mask=vstar.mask,
                )
            elif n == "bp_over_b":
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                bp2 = s_bracket(np.asarray(psi.data), np.asarray(psi.data), psi.r, psi.z)
                data = np.sqrt(np.maximum(bp2, 0.0)) / np.sqrt(np.maximum(bp2 + np.asarray(i_f.data) ** 2, np.finfo(float).tiny))
                d = dimensions()
                out = FieldResult(data=data, symbol="$|B_p|/|B|$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi.r, z=psi.z, time=psi.time, mask=psi.mask)
            elif n in {"normal curvature", "kn"}:
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                i0 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                p0 = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(psi0.r, psi0.z)
                else:
                    rmat = np.ones_like(np.asarray(psi0.data), dtype=float)
                B02 = s_bracket(np.asarray(psi0.data), np.asarray(psi0.data), psi0.r, psi0.z) / np.maximum(rmat**2, np.finfo(float).tiny) + np.asarray(i0.data) ** 2 / np.maximum(rmat**2, np.finfo(float).tiny)
                p02 = s_bracket(np.asarray(p0.data), np.asarray(p0.data), psi0.r, psi0.z)
                if not linear:
                    data = (p02 + s_bracket(np.asarray(p0.data), B02, psi0.r, psi0.z)) / (np.maximum(B02, np.finfo(float).tiny) * np.sqrt(np.maximum(p02, np.finfo(float).tiny)))
                else:
                    psi1 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    i1 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    p1 = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    if complex:
                        f1 = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=True, cgs=cgs, mks=mks, return_meta=True)
                        f1p = 1j * float(ntor) * np.asarray(f1.data)
                    else:
                        f1p = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, operation=21, cgs=cgs, mks=mks, return_meta=False)
                    bb = s_bracket(np.asarray(psi1.data), np.asarray(psi0.data), psi0.r, psi0.z) / np.maximum(rmat**2, np.finfo(float).tiny) + np.asarray(i1.data) * np.asarray(i0.data) / np.maximum(rmat**2, np.finfo(float).tiny) + a_bracket(np.asarray(psi0.data), np.asarray(f1p), psi0.r, psi0.z) / np.maximum(rmat, np.finfo(float).tiny)
                    pp = s_bracket(np.asarray(p0.data), np.asarray(p1.data), psi0.r, psi0.z)
                    data = (
                        pp
                        + s_bracket(np.asarray(p0.data), bb, psi0.r, psi0.z)
                        + 0.5 * s_bracket(np.asarray(p1.data), B02, psi0.r, psi0.z)
                        - (0.5 * pp / np.maximum(p02, np.finfo(float).tiny) + bb / np.maximum(B02, np.finfo(float).tiny)) * s_bracket(np.asarray(p0.data), B02, psi0.r, psi0.z)
                        - 2.0 * p02 * bb / np.maximum(B02, np.finfo(float).tiny)
                    ) / (np.maximum(B02, np.finfo(float).tiny) * np.sqrt(np.maximum(p02, np.finfo(float).tiny)))
                d = dimensions(l0=-1)
                out = FieldResult(data=-data, symbol="$\\kappa_n$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi0.r, z=psi0.z, time=psi0.time, mask=psi0.mask)
            elif n in {"geodesic curvature", "kg"}:
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                i0 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                p0 = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(psi0.r, psi0.z)
                else:
                    rmat = np.ones_like(np.asarray(psi0.data), dtype=float)
                B02 = s_bracket(np.asarray(psi0.data), np.asarray(psi0.data), psi0.r, psi0.z) / np.maximum(rmat**2, np.finfo(float).tiny) + np.asarray(i0.data) ** 2 / np.maximum(rmat**2, np.finfo(float).tiny)
                p02 = s_bracket(np.asarray(p0.data), np.asarray(p0.data), psi0.r, psi0.z)
                if not linear:
                    data = np.asarray(i0.data) * a_bracket(B02, np.asarray(p0.data), psi0.r, psi0.z) / np.maximum(rmat, np.finfo(float).tiny)
                else:
                    psi1 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    i1 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    p1 = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    if complex:
                        f1 = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=True, cgs=cgs, mks=mks, return_meta=True)
                        f1p = 1j * float(ntor) * np.asarray(f1.data)
                        p1p = 1j * float(ntor) * np.asarray(p1.data)
                    else:
                        f1p = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, operation=21, cgs=cgs, mks=mks, return_meta=False)
                        p1p = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, operation=21, cgs=cgs, mks=mks, return_meta=False)
                    bb = s_bracket(np.asarray(psi1.data), np.asarray(psi0.data), psi0.r, psi0.z) / np.maximum(rmat**2, np.finfo(float).tiny) + np.asarray(i1.data) * np.asarray(i0.data) / np.maximum(rmat**2, np.finfo(float).tiny) + a_bracket(np.asarray(psi0.data), np.asarray(f1p), psi0.r, psi0.z) / np.maximum(rmat, np.finfo(float).tiny)
                    pp = s_bracket(np.asarray(p0.data), np.asarray(p1.data), psi0.r, psi0.z)
                    bbp = 1j * float(ntor) * bb
                    data = (
                        np.asarray(i1.data) * a_bracket(B02, np.asarray(p0.data), psi0.r, psi0.z) / np.maximum(rmat, np.finfo(float).tiny)
                        + 2.0 * (np.asarray(i0.data) * a_bracket(bb, np.asarray(p0.data), psi0.r, psi0.z) / np.maximum(rmat, np.finfo(float).tiny) - bbp * s_bracket(np.asarray(psi0.data), np.asarray(p0.data), psi0.r, psi0.z) / np.maximum(rmat**2, np.finfo(float).tiny))
                        + (np.asarray(i0.data) * a_bracket(B02, np.asarray(p1.data), psi0.r, psi0.z) / np.maximum(rmat, np.finfo(float).tiny) + np.asarray(p1p) * s_bracket(np.asarray(psi0.data), np.asarray(p0.data), psi0.r, psi0.z) / np.maximum(rmat**2, np.finfo(float).tiny))
                        - np.asarray(i0.data) * a_bracket(B02, np.asarray(p0.data), psi0.r, psi0.z) / np.maximum(rmat, np.finfo(float).tiny) * (pp / np.maximum(p02, np.finfo(float).tiny) + 3.0 * bb / np.maximum(B02, np.finfo(float).tiny))
                    )
                d = dimensions(l0=-1)
                out = FieldResult(
                    data=data / (2.0 * np.maximum(B02, np.finfo(float).tiny) * np.sqrt(np.maximum(p02 * B02, np.finfo(float).tiny))),
                    symbol="$\\kappa_g$",
                    units=parse_units(d, cgs=cgs, mks=mks),
                    dimensions=d,
                    r=psi0.r,
                    z=psi0.z,
                    time=psi0.time,
                    mask=psi0.mask,
                )
            elif n == "k":
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(phi_f.r, phi_f.z)
                else:
                    rmat = np.ones_like(np.asarray(phi_f.data), dtype=float)
                psipsi = s_bracket(np.asarray(psi.data), np.asarray(psi.data), phi_f.r, phi_f.z)
                data = np.asarray(den.data) * (s_bracket(np.asarray(phi_f.data), np.asarray(psi.data), phi_f.r, phi_f.z) + rmat * a_bracket(np.asarray(chi.data), np.asarray(psi.data), phi_f.r, phi_f.z)) / np.maximum(psipsi, np.finfo(float).tiny)
                d = dimensions(n0=1, v0=1)
                out = FieldResult(data=data, symbol="$K$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=phi_f.r, z=phi_f.z, time=phi_f.time, mask=phi_f.mask)
            elif n in {"ideal omega", "omega_ideal"}:
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                v = read_field("v", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(phi_f.r, phi_f.z)
                else:
                    rmat = np.ones_like(np.asarray(phi_f.data), dtype=float)
                psipsi = s_bracket(np.asarray(psi.data), np.asarray(psi.data), phi_f.r, phi_f.z)
                data = (
                    np.asarray(v.data)
                    - np.asarray(i_f.data)
                    * (
                        s_bracket(np.asarray(phi_f.data), np.asarray(psi.data), phi_f.r, phi_f.z)
                        + rmat * a_bracket(np.asarray(chi.data), np.asarray(psi.data), phi_f.r, phi_f.z)
                    )
                    / np.maximum(psipsi, np.finfo(float).tiny)
                ) / np.maximum(rmat**2, np.finfo(float).tiny)
                d = dimensions(t0=-1)
                out = FieldResult(data=data, symbol="$\\omega_{ideal}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=phi_f.r, z=phi_f.z, time=phi_f.time, mask=phi_f.mask)
            elif n in {"s", "lundquist number"}:
                va = read_field("va", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                eta = read_field("eta", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                eta_data = np.asarray(eta.data)
                d = dimensions()
                out = FieldResult(
                    data=np.asarray(va.data) / np.maximum(eta_data, np.finfo(float).tiny),
                    symbol="$S$",
                    units=parse_units(d, cgs=cgs, mks=mks),
                    dimensions=d,
                    r=va.r,
                    z=va.z,
                    time=va.time,
                    mask=va.mask,
                )
            elif n == "wpar":
                u = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                v = read_field("v", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(u.r, u.z)
                else:
                    rmat = np.ones_like(np.asarray(u.data), dtype=float)
                psipsi = s_bracket(np.asarray(psi.data), np.asarray(psi.data), u.r, u.z)
                r2b2 = psipsi + np.asarray(i_f.data) ** 2
                com = laplacian(np.asarray(chi.data), u.r, u.z, toroidal=itor == 1)
                data = (
                    s_bracket(np.asarray(psi.data), a_bracket(np.asarray(u.data), np.asarray(psi.data), u.r, u.z) / np.maximum(rmat, np.finfo(float).tiny), u.r, u.z)
                    - 0.5 * rmat * a_bracket(np.asarray(u.data), psipsi / np.maximum(rmat**2, np.finfo(float).tiny), u.r, u.z)
                    - np.asarray(i_f.data) * rmat * a_bracket(np.asarray(psi.data), np.asarray(v.data) / np.maximum(rmat**2, np.finfo(float).tiny), u.r, u.z)
                    + 0.5 * np.maximum(rmat**2, np.finfo(float).tiny) * s_bracket(np.asarray(chi.data), psipsi / np.maximum(rmat**2, np.finfo(float).tiny), u.r, u.z)
                    + com * psipsi
                    - s_bracket(np.asarray(psi.data), s_bracket(np.asarray(psi.data), np.asarray(chi.data), u.r, u.z), u.r, u.z)
                ) / np.maximum(r2b2, np.finfo(float).tiny) - com / 3.0
                if itor == 1:
                    data = data + np.asarray(i_f.data) ** 2 * (
                        dx(np.asarray(chi.data), u.r) / np.maximum(rmat, np.finfo(float).tiny)
                        - dz(np.asarray(u.data), u.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                    ) / np.maximum(r2b2, np.finfo(float).tiny)
                d = dimensions(t0=-1)
                out = FieldResult(data=data, symbol="$W_{||}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=u.r, z=u.z, time=u.time, mask=u.mask)
            elif n in {"v_omega", "v_k", "v_k_n"}:
                omega_f = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                u = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(omega_f.r, omega_f.z) if itor == 1 else np.ones_like(np.asarray(omega_f.data), dtype=float)
                pol = (
                    np.maximum(rmat**2, np.finfo(float).tiny) * s_bracket(u.data, psi.data, u.r, u.z)
                    + a_bracket(chi.data, psi.data, chi.r, chi.z) / np.maximum(rmat, np.finfo(float).tiny)
                )
                spsi = np.maximum(s_bracket(psi.data, psi.data, psi.r, psi.z), np.finfo(float).tiny)
                if n == "v_omega":
                    data = np.asarray(omega_f.data) - np.asarray(i_f.data) * pol / (np.maximum(rmat**2, np.finfo(float).tiny) * spsi)
                    d = dimensions(t0=-1)
                    sym = "$\\omega_v$"
                elif n == "v_k":
                    den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    data = np.asarray(den.data) * pol / spsi
                    d = dimensions(v0=1, n0=1, b0=-1)
                    sym = "$K$"
                else:
                    data = pol / spsi
                    d = dimensions(v0=1, b0=-1)
                    sym = "$K/n$"
                out = FieldResult(data=data, symbol=sym, units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=omega_f.r, z=omega_f.z, time=omega_f.time, mask=omega_f.mask)
            elif n in {"omega_e", "ve_omega", "ve_k", "omega_perp_e", "rho_i"}:
                db = float(read_parameter("db", filename=filename))
                omega_f = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                u = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(omega_f.r, omega_f.z) if itor == 1 else np.ones_like(np.asarray(omega_f.data), dtype=float)
                spsi = np.maximum(s_bracket(psi.data, psi.data, psi.r, psi.z), np.finfo(float).tiny)
                i_over = np.asarray(i_f.data) / np.maximum(np.maximum(rmat**2, np.finfo(float).tiny) * spsi, np.finfo(float).tiny)
                pol = (
                    np.maximum(rmat**2, np.finfo(float).tiny) * s_bracket(u.data, psi.data, u.r, u.z)
                    + a_bracket(chi.data, psi.data, chi.r, chi.z) / np.maximum(rmat, np.finfo(float).tiny)
                )

                if n == "omega_e":
                    if int(read_parameter("ivform", filename=filename)) == 0:
                        data = np.asarray(omega_f.data) - (db * np.asarray(read_field("jy", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks))) / np.maximum(np.asarray(den.data) * np.maximum(rmat, np.finfo(float).tiny), np.finfo(float).tiny)
                    else:
                        data = np.asarray(omega_f.data) - (db * np.asarray(read_field("jy", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks))) / np.maximum(np.asarray(den.data), np.finfo(float).tiny)
                    d = dimensions(t0=-1)
                    sym = "$\\omega_e$"
                elif n == "ve_omega":
                    psi_lp = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=7, cgs=cgs, mks=mks, return_meta=True)
                    data = (
                        np.asarray(omega_f.data)
                        - i_over * (pol - (db / np.maximum(np.asarray(den.data), np.finfo(float).tiny)) * s_bracket(i_f.data, psi.data, psi.r, psi.z))
                        + (db / np.maximum(np.asarray(den.data), np.finfo(float).tiny))
                        * (np.asarray(psi_lp.data) - itor * dx(psi.data, psi.r) / np.maximum(rmat, np.finfo(float).tiny))
                        / np.maximum(rmat**2, np.finfo(float).tiny)
                    )
                    d = dimensions(t0=-1)
                    sym = "$\\omega_{ve}$"
                elif n == "ve_k":
                    data = (
                        1.0 / spsi
                        * (
                            np.asarray(den.data) * pol
                            - db * s_bracket(i_f.data, psi.data, psi.r, psi.z)
                        )
                    )
                    d = dimensions(v0=1, n0=1, b0=-1)
                    sym = "$K_e$"
                elif n == "omega_perp_e":
                    psi_lp = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=7, cgs=cgs, mks=mks, return_meta=True)
                    omega_e = (
                        np.asarray(omega_f.data)
                        - i_over * (pol - (db / np.maximum(np.asarray(den.data), np.finfo(float).tiny)) * s_bracket(i_f.data, psi.data, psi.r, psi.z))
                        + (db / np.maximum(np.asarray(den.data), np.finfo(float).tiny))
                        * (np.asarray(psi_lp.data) - itor * dx(psi.data, psi.r) / np.maximum(rmat, np.finfo(float).tiny))
                        / np.maximum(rmat**2, np.finfo(float).tiny)
                    )
                    bp = np.sqrt(np.maximum(spsi, np.finfo(float).tiny)) / np.maximum(rmat, np.finfo(float).tiny)
                    btot = np.sqrt(np.maximum(spsi + np.asarray(i_f.data) ** 2, np.finfo(float).tiny)) / np.maximum(rmat, np.finfo(float).tiny)
                    data = (bp / np.maximum(btot, np.finfo(float).tiny)) * omega_e
                    d = dimensions(t0=-1)
                    sym = "$\\omega_{\\perp e}$"
                else:
                    ti = read_field("ti", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    btot = read_field("b", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    data = db * np.sqrt(np.maximum(2.0 * np.asarray(ti.data), 0.0)) / np.maximum(np.asarray(btot.data), np.finfo(float).tiny)
                    d = dimensions(l0=1)
                    sym = "$\\rho_i$"
                out = FieldResult(data=data, symbol=sym, units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=omega_f.r, z=omega_f.z, time=omega_f.time, mask=omega_f.mask)
            elif n == "xbdotgradt":
                te0 = read_field("te", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(te0.r, te0.z) if itor == 1 else np.ones_like(np.asarray(te0.data), dtype=float)
                if ilin == 1:
                    te1 = read_field("te", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                    psi1 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                    data = a_bracket(np.asarray(te1.data), np.asarray(psi0.data), te0.r, te0.z) / np.maximum(rmat, np.finfo(float).tiny) + a_bracket(np.asarray(te0.data), np.asarray(psi1.data), te0.r, te0.z) / np.maximum(rmat, np.finfo(float).tiny)
                    if ntor != 0:
                        i1 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                        i0 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                        f1 = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                        data = data + 1j * ntor * (np.asarray(te0.data) * np.asarray(i1.data) + np.asarray(te1.data) * np.asarray(i0.data) - s_bracket(np.asarray(te0.data), np.asarray(f1.data), te0.r, te0.z))
                else:
                    data = 2.0 * a_bracket(np.asarray(te0.data), np.asarray(psi0.data), te0.r, te0.z) / np.maximum(rmat, np.finfo(float).tiny)
                d = dimensions(l0=1)
                out = FieldResult(data=data, symbol="$B\\cdot\\nabla T$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=te0.r, z=te0.z, time=te0.time, mask=te0.mask)
            elif n == "xbdotgradp":
                p0 = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(p0.r, p0.z) if itor == 1 else np.ones_like(np.asarray(p0.data), dtype=float)
                if ilin == 1:
                    p1 = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    psi1 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    data = a_bracket(np.asarray(p1.data), np.asarray(psi0.data), p0.r, p0.z) / np.maximum(rmat, np.finfo(float).tiny) + a_bracket(np.asarray(p0.data), np.asarray(psi1.data), p0.r, p0.z) / np.maximum(rmat, np.finfo(float).tiny)
                    if ntor != 0:
                        i0 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                        f1 = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                        data = data + 1j * ntor * (np.asarray(p1.data) * np.asarray(i0.data) / np.maximum(rmat**2, np.finfo(float).tiny) - s_bracket(np.asarray(p0.data), np.asarray(f1.data), p0.r, p0.z))
                else:
                    data = np.zeros_like(np.asarray(p0.data), dtype=np.complex128 if complex else float)
                d = dimensions(p0=1, b0=1, l0=-1)
                out = FieldResult(data=data, symbol="$B\\cdot\\nabla p$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=p0.r, z=p0.z, time=p0.time, mask=p0.mask)
            elif n == "omtest" or n == "omtest_i":
                pert = "_i" if n.endswith("_i") else ""
                v_omega = read_field(f"v_omega{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                v_k_n = read_field(f"v_k_n{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field(f"i{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(v_omega.r, v_omega.z) if itor == 1 else np.ones_like(np.asarray(v_omega.data), dtype=float)
                data = np.asarray(v_omega.data) + np.asarray(i_f.data) * np.asarray(v_k_n.data) / np.maximum(rmat**2, np.finfo(float).tiny)
                d = dimensions(t0=-1)
                out = FieldResult(data=data, symbol="$\\Omega_{test}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=v_omega.r, z=v_omega.z, time=v_omega.time, mask=v_omega.mask)
            elif n == "vpar":
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                w = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                i0 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(phi_f.r, phi_f.z) if itor == 1 else np.ones_like(np.asarray(phi_f.data), dtype=float)
                pp = s_bracket(psi0.data, psi0.data, psi0.r, psi0.z)
                b2 = (pp + np.asarray(i0.data) ** 2) / np.maximum(rmat**2, np.finfo(float).tiny)
                ivform = int(read_parameter("ivform", filename=filename))
                if ivform == 0:
                    data = (
                        (s_bracket(phi_f.data, psi0.data, phi_f.r, phi_f.z) + np.asarray(w.data) * np.asarray(i0.data)) / np.maximum(rmat**2, np.finfo(float).tiny)
                        + a_bracket(chi.data, psi0.data, chi.r, chi.z) / np.maximum(rmat, np.finfo(float).tiny)
                    ) / np.sqrt(np.maximum(b2, np.finfo(float).tiny))
                else:
                    data = (
                        s_bracket(phi_f.data, psi0.data, phi_f.r, phi_f.z)
                        + np.asarray(w.data) * np.asarray(i0.data)
                        + a_bracket(chi.data, psi0.data, chi.r, chi.z) / np.maximum(rmat**3, np.finfo(float).tiny)
                    ) / np.sqrt(np.maximum(b2, np.finfo(float).tiny))
                d = dimensions(v0=1)
                out = FieldResult(data=data, symbol="$u_{||}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=phi_f.r, z=phi_f.z, time=phi_f.time, mask=phi_f.mask)
            elif n == "vpol":
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(phi_f.r, phi_f.z) if itor == 1 else np.ones_like(np.asarray(phi_f.data), dtype=float)
                ivform = int(read_parameter("ivform", filename=filename))
                if ivform == 0:
                    data = np.zeros_like(np.asarray(phi_f.data), dtype=float)
                else:
                    data = (
                        np.maximum(rmat, np.finfo(float).tiny)
                        * (
                            s_bracket(phi_f.data, psi0.data, phi_f.r, phi_f.z)
                            + a_bracket(chi.data, psi0.data, chi.r, chi.z) / np.maximum(rmat**3, np.finfo(float).tiny)
                        )
                        / np.sqrt(np.maximum(s_bracket(psi0.data, psi0.data, psi0.r, psi0.z), np.finfo(float).tiny))
                    )
                d = dimensions(v0=1)
                out = FieldResult(data=data, symbol="$u_{pol}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=phi_f.r, z=phi_f.z, time=phi_f.time, mask=phi_f.mask)
            elif n == "vperp":
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                w = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                i0 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(phi_f.r, phi_f.z) if itor == 1 else np.ones_like(np.asarray(phi_f.data), dtype=float)
                pp = s_bracket(psi0.data, psi0.data, psi0.r, psi0.z)
                b2 = (pp + np.asarray(i0.data) ** 2) / np.maximum(rmat**2, np.finfo(float).tiny)
                ivform = int(read_parameter("ivform", filename=filename))
                if ivform == 0:
                    data = np.zeros_like(np.asarray(phi_f.data), dtype=float)
                else:
                    data = (
                        -np.asarray(i0.data) * s_bracket(phi_f.data, psi0.data, phi_f.r, phi_f.z)
                        + np.asarray(w.data) * pp
                        - np.asarray(i0.data) * a_bracket(chi.data, psi0.data, chi.r, chi.z) / np.maximum(rmat**3, np.finfo(float).tiny)
                    ) / np.sqrt(np.maximum(b2 * pp, np.finfo(float).tiny))
                d = dimensions(v0=1)
                out = FieldResult(data=data, symbol="$u_{\\perp}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=phi_f.r, z=phi_f.z, time=phi_f.time, mask=phi_f.mask)
            elif n in {"chi_perp", "chi_par", "mu_perp"}:
                den0 = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                if n == "chi_perp":
                    kappa = read_field("kappa", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    num = np.asarray(kappa.data)
                    sym = "$\\chi_\\perp$"
                elif n == "chi_par":
                    kappar = float(read_parameter("kappar", filename=filename))
                    num = np.ones_like(np.asarray(den0.data), dtype=float) * kappar
                    sym = "$\\chi_{||}$"
                else:
                    visc = read_field("visc", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    num = np.asarray(visc.data)
                    sym = "$\\mu_\\perp$"
                data = num / np.maximum(np.asarray(den0.data), np.finfo(float).tiny)
                d = dimensions(l0=2, t0=-1)
                out = FieldResult(data=data, symbol=sym, units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=den0.r, z=den0.z, time=den0.time, mask=den0.mask)
            elif n == "iota":
                minor_r = read_field("r", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(psi.r, psi.z) if itor == 1 else np.ones_like(np.asarray(psi.data), dtype=float)
                bt = np.sqrt(np.asarray(i_f.data) ** 2 / np.maximum(rmat**2, np.finfo(float).tiny))
                bp = np.sqrt(np.maximum(s_bracket(psi.data, psi.data, psi.r, psi.z), np.finfo(float).tiny) / np.maximum(rmat**2, np.finfo(float).tiny))
                data = 2.0 * np.pi * (rmat * bp) / np.maximum(np.asarray(minor_r.data), np.finfo(float).tiny) * bt
                out = FieldResult(data=data, symbol="$\\iota$", units="", dimensions=dimensions(), r=psi.r, z=psi.z, time=psi.time, mask=psi.mask)
            elif n == "vn":
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(phi_f.r, phi_f.z) if itor == 1 else np.ones_like(np.asarray(phi_f.data), dtype=float)
                if int(read_parameter("ivform", filename=filename)) == 0:
                    data = (
                        a_bracket(psi0.data, phi_f.data, phi_f.r, phi_f.z) / np.maximum(rmat, np.finfo(float).tiny)
                        + s_bracket(psi0.data, chi.data, psi0.r, psi0.z)
                    ) / np.sqrt(np.maximum(s_bracket(psi0.data, psi0.data, psi0.r, psi0.z), np.finfo(float).tiny))
                else:
                    data = (
                        a_bracket(psi0.data, phi_f.data, phi_f.r, phi_f.z) * rmat
                        + s_bracket(psi0.data, chi.data, psi0.r, psi0.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                    ) / np.sqrt(np.maximum(s_bracket(psi0.data, psi0.data, psi0.r, psi0.z), np.finfo(float).tiny))
                d = dimensions(v0=1)
                out = FieldResult(data=data, symbol="$u_n$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=phi_f.r, z=phi_f.z, time=phi_f.time, mask=phi_f.mask)
            elif n == "vx":
                chi_r = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, operation=2, cgs=cgs, mks=mks, return_meta=True)
                phi_z = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, operation=3, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(chi_r.r, chi_r.z) if itor == 1 else np.ones_like(np.asarray(chi_r.data), dtype=float)
                if int(read_parameter("ivform", filename=filename)) == 0:
                    data = -np.asarray(phi_z.data) / np.maximum(rmat, np.finfo(float).tiny) + np.asarray(chi_r.data)
                else:
                    data = -rmat * np.asarray(phi_z.data) + np.asarray(chi_r.data) / np.maximum(rmat**2, np.finfo(float).tiny)
                d = dimensions(v0=1)
                out = FieldResult(data=data, symbol="$u_R$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=chi_r.r, z=chi_r.z, time=chi_r.time, mask=chi_r.mask)
            elif n == "vz":
                chi_z = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, operation=3, cgs=cgs, mks=mks, return_meta=True)
                phi_r = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, operation=2, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(chi_z.r, chi_z.z) if itor == 1 else np.ones_like(np.asarray(chi_z.data), dtype=float)
                if int(read_parameter("ivform", filename=filename)) == 0:
                    data = np.asarray(phi_r.data) / np.maximum(rmat, np.finfo(float).tiny) + np.asarray(chi_z.data)
                else:
                    data = rmat * np.asarray(phi_r.data) + np.asarray(chi_z.data) / np.maximum(rmat**2, np.finfo(float).tiny)
                d = dimensions(v0=1)
                out = FieldResult(data=data, symbol="$u_Z$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=chi_z.r, z=chi_z.z, time=chi_z.time, mask=chi_z.mask)
            elif n == "bn":
                psi0_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=2, cgs=cgs, mks=mks, return_meta=True)
                psi0_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=3, cgs=cgs, mks=mks, return_meta=True)
                psi_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, operation=2, cgs=cgs, mks=mks, return_meta=True)
                psi_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, operation=3, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(psi_r.r, psi_r.z) if itor == 1 else np.ones_like(np.asarray(psi_r.data), dtype=float)
                data = (np.asarray(psi_z.data) * np.asarray(psi0_r.data) - np.asarray(psi_r.data) * np.asarray(psi0_z.data)) / np.maximum(rmat, np.finfo(float).tiny)
                if ntor != 0:
                    f_r = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, operation=2, cgs=cgs, mks=mks, return_meta=True)
                    f_z = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, operation=3, cgs=cgs, mks=mks, return_meta=True)
                    data = data + 1j * ntor * (np.asarray(f_r.data) * np.asarray(psi0_r.data) + np.asarray(f_z.data) * np.asarray(psi0_z.data))
                data = data / np.sqrt(np.maximum(np.asarray(psi0_r.data) ** 2 + np.asarray(psi0_z.data) ** 2, np.finfo(float).tiny))
                d = dimensions(b0=1)
                out = FieldResult(data=data, symbol="$B_n$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi_r.r, z=psi_r.z, time=psi_r.time, mask=psi_r.mask)
            elif n == "bpol":
                psi0_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=2, cgs=cgs, mks=mks, return_meta=True)
                psi0_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=3, cgs=cgs, mks=mks, return_meta=True)
                psi_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=2, cgs=cgs, mks=mks, return_meta=True)
                psi_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=3, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(psi_r.r, psi_r.z) if itor == 1 else np.ones_like(np.asarray(psi_r.data), dtype=float)
                data = (np.asarray(psi_z.data) * np.asarray(psi0_z.data) + np.asarray(psi_r.data) * np.asarray(psi0_r.data)) / np.maximum(rmat, np.finfo(float).tiny)
                if ntor != 0:
                    fp_r = read_field("fp", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=2, cgs=cgs, mks=mks, return_meta=True)
                    fp_z = read_field("fp", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=3, cgs=cgs, mks=mks, return_meta=True)
                    data = data + (np.asarray(fp_z.data) * np.asarray(psi0_r.data) - np.asarray(fp_r.data) * np.asarray(psi0_z.data))
                data = data / np.sqrt(np.maximum(np.asarray(psi0_r.data) ** 2 + np.asarray(psi0_z.data) ** 2, np.finfo(float).tiny))
                d = dimensions(b0=1)
                out = FieldResult(data=data, symbol="$B_{pol}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi_r.r, z=psi_r.z, time=psi_r.time, mask=psi_r.mask)
            elif n in {"jpar", "jpar_lin", "jb", "jpar_b"}:
                psi_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi_lp = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=7, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                i_r = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                i_z = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(psi_r.r, psi_r.z) if itor == 1 else np.ones_like(np.asarray(psi_r.data), dtype=float)
                rfac = 1j * ntor if itor == 1 else 1j * ntor / float(read_parameter("rzero", filename=filename))
                psi_rp = np.asarray(psi_r.data) * rfac
                psi_zp = np.asarray(psi_z.data) * rfac

                f_r = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                f_z = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                f_rpp = np.asarray(f_r.data) * (rfac**2)
                f_zpp = np.asarray(f_z.data) * (rfac**2)

                psi0_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=2, cgs=cgs, mks=mks, return_meta=True)
                psi0_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=3, cgs=cgs, mks=mks, return_meta=True)
                i0 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                i0_r = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=2, cgs=cgs, mks=mks, return_meta=True)
                i0_z = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=3, cgs=cgs, mks=mks, return_meta=True)

                if n == "jpar_b":
                    psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                    jphi = read_field("jphi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    i1 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    b0 = s_bracket(psi0.data, psi0.data, psi0.r, psi0.z) / np.maximum(rmat**2, np.finfo(float).tiny) + np.asarray(i0.data) ** 2 / np.maximum(rmat**2, np.finfo(float).tiny)
                    data = (
                        s_bracket(i1.data, psi0.data, psi0.r, psi0.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                        - np.asarray(jphi.data) * np.asarray(i0.data) / np.maximum(rmat**2, np.finfo(float).tiny)
                    ) / np.maximum(b0, np.finfo(float).tiny)
                    d = dimensions(j0=1, b0=-1)
                    out = FieldResult(data=data, symbol="$J_{||}/B$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi_r.r, z=psi_r.z, time=psi_r.time, mask=psi_r.mask)
                else:
                    term_j1b0 = (
                        -np.asarray(i0.data) * (np.asarray(psi_lp.data) - itor * np.asarray(psi_r.data) / np.maximum(rmat, np.finfo(float).tiny)) / np.maximum(rmat**2, np.finfo(float).tiny)
                        + (np.asarray(psi0_r.data) * (np.asarray(i_r.data) + f_rpp) + np.asarray(psi0_z.data) * (np.asarray(i_z.data) + f_zpp)) / np.maximum(rmat**2, np.finfo(float).tiny)
                        - (np.asarray(psi0_z.data) * psi_rp - np.asarray(psi0_r.data) * psi_zp) / np.maximum(rmat**3, np.finfo(float).tiny)
                    )
                    if n == "jb":
                        d = dimensions(j0=1, b0=1)
                        out = FieldResult(data=term_j1b0, symbol="$J\\cdot B$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi_r.r, z=psi_r.z, time=psi_r.time, mask=psi_r.mask)
                    else:
                        b0mag = np.sqrt(np.maximum(np.asarray(psi0_r.data) ** 2 + np.asarray(psi0_z.data) ** 2 + np.asarray(i0.data) ** 2, np.finfo(float).tiny)) / np.maximum(rmat, np.finfo(float).tiny)
                        data = term_j1b0 / np.maximum(b0mag, np.finfo(float).tiny)
                        if n == "jpar_lin":
                            i1 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                            f_rp = np.asarray(f_r.data) * rfac
                            f_zp = np.asarray(f_z.data) * rfac
                            psi0_lp = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=7, cgs=cgs, mks=mks, return_meta=True)
                            b0b1 = (
                                (np.asarray(psi0_r.data) * np.asarray(psi_r.data) + np.asarray(psi0_z.data) * np.asarray(psi_z.data) + np.asarray(i0.data) * np.asarray(i1.data)) / np.maximum(rmat**2, np.finfo(float).tiny)
                                + (np.asarray(psi0_z.data) * f_rp - np.asarray(psi0_r.data) * f_zp) / np.maximum(rmat, np.finfo(float).tiny)
                            )
                            j0b0 = (
                                -np.asarray(i0.data) * (np.asarray(psi0_lp.data) - itor * np.asarray(psi0_r.data) / np.maximum(rmat, np.finfo(float).tiny))
                                + (np.asarray(psi0_r.data) * np.asarray(i0_r.data) + np.asarray(psi0_z.data) * np.asarray(i0_z.data)) / np.maximum(rmat**2, np.finfo(float).tiny)
                            )
                            j0b1 = (
                                -np.asarray(i1.data) * (np.asarray(psi0_lp.data) - itor * np.asarray(psi0_r.data) / np.maximum(rmat, np.finfo(float).tiny)) / np.maximum(rmat**2, np.finfo(float).tiny)
                                + (np.asarray(psi_r.data) * np.asarray(i0_r.data) + np.asarray(psi_z.data) * np.asarray(i0_z.data)) / np.maximum(rmat**2, np.finfo(float).tiny)
                                + (np.asarray(i0_z.data) * f_rp - np.asarray(i0_r.data) * f_zp) / np.maximum(rmat, np.finfo(float).tiny)
                            ) / np.maximum(b0mag, np.finfo(float).tiny)
                            data = data + j0b1 - j0b0 * b0b1 / np.maximum(b0mag**3, np.finfo(float).tiny)
                        d = dimensions(j0=1)
                        out = FieldResult(data=data, symbol="$J_{||}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi_r.r, z=psi_r.z, time=psi_r.time, mask=psi_r.mask)
            elif n == "jn" or n == "jn_i":
                pert = "_i" if n.endswith("_i") else ""
                psi0_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=2, cgs=cgs, mks=mks, return_meta=True)
                psi0_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=3, cgs=cgs, mks=mks, return_meta=True)
                i_r = read_field(f"i{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=2, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                i_z = read_field(f"i{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=3, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(psi0_r.r, psi0_r.z) if itor == 1 else np.ones_like(np.asarray(psi0_r.data), dtype=float)
                psipsi = np.sqrt(np.maximum(np.asarray(psi0_r.data) ** 2 + np.asarray(psi0_z.data) ** 2, np.finfo(float).tiny))
                data = (np.asarray(i_z.data) * np.asarray(psi0_r.data) - np.asarray(i_r.data) * np.asarray(psi0_z.data)) / np.maximum(rmat * psipsi, np.finfo(float).tiny)
                if ntor != 0:
                    psi_r = read_field(f"psi{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=2, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    psi_z = read_field(f"psi{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=3, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    fp_r = read_field(f"fp{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=2, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    fp_z = read_field(f"fp{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=3, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    fp_r_data = np.asarray(fp_r.data)
                    fp_z_data = np.asarray(fp_z.data)
                    data = (
                        data
                        - 1j * ntor * (np.asarray(psi_r.data) * np.asarray(psi0_r.data) + np.asarray(psi_z.data) * np.asarray(psi0_z.data)) / np.maximum(rmat**2 * psipsi, np.finfo(float).tiny)
                        + 1j * ntor * (fp_z_data * np.asarray(psi0_r.data) - fp_r_data * np.asarray(psi0_z.data)) / np.maximum(rmat * psipsi, np.finfo(float).tiny)
                    )
                d = dimensions(j0=1)
                out = FieldResult(data=data, symbol="$J_n$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi0_r.r, z=psi0_r.z, time=psi0_r.time, mask=psi0_r.mask)
            elif n == "flux_denm":
                denm = read_field("denm", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                denm_data = np.asarray(denm.data)
                base = denm
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                num = denm_data * s_bracket(den.data, psi.data, den.r, den.z)
                deno = np.sqrt(np.maximum(s_bracket(psi.data, psi.data, psi.r, psi.z), np.finfo(float).tiny))
                d = dimensions(n0=1, l0=1, t0=-1)
                out = FieldResult(data=num / deno, symbol="$\\Gamma_{D\\nabla n}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=base.r, z=base.z, time=base.time, mask=base.mask)
            elif n == "flux_nv":
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                u = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                f_p = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, operation=11, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(den.r, den.z)
                bp2 = (
                    s_bracket(psi.data, psi.data, psi.r, psi.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                    + 2.0 * a_bracket(psi.data, f_p.data, psi.r, psi.z) / np.maximum(rmat, np.finfo(float).tiny)
                    + s_bracket(f_p.data, f_p.data, f_p.r, f_p.z)
                )
                bpdotv = (
                    np.maximum(rmat, np.finfo(float).tiny) * s_bracket(u.data, psi.data, u.r, u.z)
                    + a_bracket(chi.data, psi.data, chi.r, chi.z) / np.maximum(rmat**3, np.finfo(float).tiny)
                    + a_bracket(f_p.data, psi.data, f_p.r, f_p.z) / np.maximum(rmat, np.finfo(float).tiny)
                    - s_bracket(chi.data, f_p.data, chi.r, chi.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                ) / np.sqrt(np.maximum(bp2, np.finfo(float).tiny))
                d = dimensions(n0=1, v0=1)
                out = FieldResult(data=np.asarray(den.data) * bpdotv, symbol="$\\Gamma_{||}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=den.r, z=den.z, time=den.time, mask=den.mask)
            elif n in {"exb_flux_n", "exb_flux_pe", "exb_flux_pi"}:
                qname = "den" if n == "exb_flux_n" else ("pe" if n == "exb_flux_pe" else "pi")
                q = read_field(qname, filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ex = read_field("e_r", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ey = read_field("e_phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ez = read_field("e_z", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bx0 = read_field("bx", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                by0 = read_field("by", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bz0 = read_field("bz", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                b2 = np.asarray(bx0.data) ** 2 + np.asarray(by0.data) ** 2 + np.asarray(bz0.data) ** 2
                deno = np.maximum(b2 * np.sqrt(np.maximum(b2 - np.asarray(by0.data) ** 2, np.finfo(float).tiny)), np.finfo(float).tiny)
                data = np.conj(np.asarray(q.data)) * (
                    (np.asarray(ex.data) * np.asarray(bx0.data) + np.asarray(ey.data) * np.asarray(by0.data) + np.asarray(ez.data) * np.asarray(bz0.data)) * np.asarray(by0.data)
                    - np.asarray(ey.data) * b2
                ) / deno
                d = dimensions(n0=1, v0=1)
                out = FieldResult(data=data, symbol="$\\Gamma_{ExB}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=q.r, z=q.z, time=q.time, mask=q.mask)
            elif n == "flux_nv2":
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                u = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(den.r, den.z)
                data = -np.conj(np.asarray(den.data)) * (
                    np.maximum(rmat, np.finfo(float).tiny) * a_bracket(psi0.data, u.data, den.r, den.z)
                    + s_bracket(psi0.data, chi.data, den.r, den.z) / np.maximum(rmat, np.finfo(float).tiny)
                ) / np.sqrt(np.maximum(s_bracket(psi0.data, psi0.data, den.r, den.z), np.finfo(float).tiny))
                d = dimensions(l0=-3, t0=-1)
                out = FieldResult(data=data, symbol="$\\Gamma_n$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=den.r, z=den.z, time=den.time, mask=den.mask)
            elif n in {"parallel heat flux", "qpar"}:
                p = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                p_p = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, operation=11, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                f_p = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, operation=11, cgs=cgs, mks=mks, return_meta=True)
                kappar = float(read_parameter("kappar", filename=filename))
                rmat = radius_matrix(p.r, p.z)
                te = np.asarray(p.data) / np.maximum(np.asarray(den.data), np.finfo(float).tiny)
                te_p = np.asarray(p_p.data) / np.maximum(np.asarray(den.data), np.finfo(float).tiny)
                bp2 = (
                    s_bracket(psi.data, psi.data, p.r, p.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                    + 2.0 * a_bracket(psi.data, f_p.data, p.r, p.z) / np.maximum(rmat, np.finfo(float).tiny)
                    + s_bracket(f_p.data, f_p.data, p.r, p.z)
                )
                b2 = bp2 + np.asarray(i_f.data) ** 2 / np.maximum(rmat**2, np.finfo(float).tiny)
                bdotgradte = (
                    a_bracket(te, psi.data, p.r, p.z) / np.maximum(rmat, np.finfo(float).tiny)
                    + np.asarray(i_f.data) * te_p / np.maximum(rmat**2, np.finfo(float).tiny)
                    - s_bracket(f_p.data, te, p.r, p.z)
                )
                br = -dz(psi.data, p.z) / np.maximum(rmat, np.finfo(float).tiny) - dx(f_p.data, p.r)
                bbter = br * bdotgradte / np.maximum(b2, np.finfo(float).tiny)
                bz = dx(psi.data, p.r) / np.maximum(rmat, np.finfo(float).tiny) - dz(f_p.data, p.z)
                bbtez = bz * bdotgradte / np.maximum(b2, np.finfo(float).tiny)
                if rvector:
                    data = -kappar * bbter
                elif zvector:
                    data = -kappar * bbtez
                else:
                    data = -kappar * bdotgradte * np.sqrt(np.maximum(bp2, np.finfo(float).tiny)) / np.maximum(b2, np.finfo(float).tiny)
                d = dimensions(p0=1, v0=1)
                out = FieldResult(data=data, symbol="$q_{||}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=p.r, z=p.z, time=p.time, mask=p.mask)
            elif n in {"convective heat flux", "qcon"}:
                p = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                phi_f = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                gamma = float(read_parameter("gam", filename=filename))
                rmat = radius_matrix(p.r, p.z)
                if rvector:
                    v = -np.maximum(rmat, np.finfo(float).tiny) * dz(phi_f.data, p.z) + dx(chi.data, p.r) / np.maximum(rmat**2, np.finfo(float).tiny)
                elif zvector:
                    v = np.maximum(rmat, np.finfo(float).tiny) * dx(phi_f.data, p.r) + dz(chi.data, p.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                else:
                    omega = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    v = np.sqrt(
                        np.maximum(
                            np.maximum(rmat**2, np.finfo(float).tiny) * s_bracket(phi_f.data, phi_f.data, p.r, p.z)
                            + np.maximum(rmat**2, np.finfo(float).tiny) * np.asarray(omega.data) ** 2
                            + s_bracket(chi.data, chi.data, p.r, p.z) / np.maximum(rmat**4, np.finfo(float).tiny)
                            + 2.0 * a_bracket(chi.data, phi_f.data, p.r, p.z) / np.maximum(rmat, np.finfo(float).tiny),
                            0.0,
                        )
                    )
                data = -gamma / max(gamma - 1.0, np.finfo(float).tiny) * np.asarray(p.data) * v
                d = dimensions(p0=1, v0=1)
                out = FieldResult(data=data, symbol="$q_{con}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=p.r, z=p.z, time=p.time, mask=p.mask)
            elif n == "dbndt":
                psi0_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=2, cgs=cgs, mks=mks, return_meta=True)
                psi0_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=3, cgs=cgs, mks=mks, return_meta=True)
                i0 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                i0_r = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=2, cgs=cgs, mks=mks, return_meta=True)
                i0_z = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=3, cgs=cgs, mks=mks, return_meta=True)
                w0 = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                f_f = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                u = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                w_r = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=2, cgs=cgs, mks=mks, return_meta=True)
                w_z = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=3, cgs=cgs, mks=mks, return_meta=True)
                chi_lp = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=7, cgs=cgs, mks=mks, return_meta=True)
                chi_r = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=2, cgs=cgs, mks=mks, return_meta=True)
                chi_z = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=3, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(psi.r, psi.z) if itor == 1 else np.ones_like(np.asarray(psi.data), dtype=float)
                rfac = 1j * ntor
                data = (
                    -(
                        np.asarray(i0.data) * (np.asarray(chi_lp.data) - np.asarray(chi_r.data) / np.maximum(rmat, np.finfo(float).tiny)) / np.maximum(rmat**2, np.finfo(float).tiny)
                        + (np.asarray(i0_r.data) * np.asarray(chi_r.data) + np.asarray(i0_z.data) * np.asarray(chi_z.data)) / np.maximum(rmat**2, np.finfo(float).tiny)
                        - 2.0 * np.asarray(i0.data) * np.asarray(chi_r.data) / np.maximum(rmat**3, np.finfo(float).tiny)
                        + np.maximum(rmat**2, np.finfo(float).tiny) * np.asarray(w0.data) * delstar(rfac * np.asarray(f_f.data), psi.r, psi.z, toroidal=(itor == 1))
                        + s_bracket(np.asarray(w0.data) * np.maximum(rmat**2, np.finfo(float).tiny), rfac * np.asarray(f_f.data), psi.r, psi.z)
                    ) / np.maximum(rmat, np.finfo(float).tiny)
                    - a_bracket(i0.data, u.data, psi.r, psi.z)
                    + a_bracket(w0.data, psi.data, psi.r, psi.z)
                    + (np.asarray(w_z.data) * np.asarray(psi0_r.data) - np.asarray(w_r.data) * np.asarray(psi0_z.data))
                )
                d = dimensions(b0=1, t0=-1)
                out = FieldResult(data=data, symbol="$\\nabla\\times(V\\times B)$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi.r, z=psi.z, time=psi.time, mask=psi.mask)
            elif n == "curletaj":
                eta = read_field("eta", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                f_f = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(psi.r, psi.z) if itor == 1 else np.ones_like(np.asarray(psi.data), dtype=float)
                rfac = 1j * ntor
                data = (
                    -s_bracket(eta.data, np.asarray(i_f.data) + rfac**2 * np.asarray(f_f.data), psi.r, psi.z) / np.maximum(rmat, np.finfo(float).tiny)
                    - np.asarray(eta.data) * delstar(np.asarray(i_f.data) + rfac**2 * np.asarray(f_f.data), psi.r, psi.z, toroidal=(itor == 1)) / np.maximum(rmat, np.finfo(float).tiny)
                    + a_bracket(np.asarray(eta.data) / np.maximum(rmat**2, np.finfo(float).tiny), rfac * np.asarray(psi.data), psi.r, psi.z)
                )
                d = dimensions(b0=1, t0=-1)
                out = FieldResult(data=data, symbol="$\\nabla\\times(\\eta J)$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi.r, z=psi.z, time=psi.time, mask=psi.mask)
            elif n == "torque":
                force_phi = read_field("force_phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(force_phi.r, force_phi.z) if itor == 1 else np.ones_like(np.asarray(force_phi.data), dtype=float)
                d = dimensions(p0=1)
                out = FieldResult(data=np.asarray(force_phi.data) * rmat, symbol="$\\tau$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=force_phi.r, z=force_phi.z, time=force_phi.time, mask=force_phi.mask)
            elif n in {"jxb_x", "jybz", "jzby", "jxb_y", "jzbx", "jxbz", "jxb_z", "jxby", "jybx"}:
                bx = read_field("bx", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                by = read_field("by", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bz = read_field("bz", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                jx = read_field("jx", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                jz = read_field("jz", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                jy = read_field("jy_plasma", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                if n == "jxb_x":
                    data = np.asarray(jy.data) * np.asarray(bz.data) - np.asarray(jz.data) * np.asarray(by.data)
                elif n == "jybz":
                    data = np.asarray(jy.data) * np.asarray(bz.data)
                elif n == "jzby":
                    data = -np.asarray(jz.data) * np.asarray(by.data)
                elif n == "jxb_y":
                    data = np.asarray(jz.data) * np.conj(np.asarray(bx.data)) - np.asarray(jx.data) * np.conj(np.asarray(bz.data))
                elif n == "jzbx":
                    data = np.asarray(jz.data) * np.conj(np.asarray(bx.data))
                elif n == "jxbz":
                    data = -np.asarray(jx.data) * np.conj(np.asarray(bz.data))
                elif n == "jxb_z":
                    data = np.asarray(jx.data) * np.asarray(by.data) - np.asarray(jy.data) * np.asarray(bx.data)
                elif n == "jxby":
                    data = np.asarray(jx.data) * np.asarray(by.data)
                else:
                    data = -np.asarray(jy.data) * np.asarray(bx.data)
                d = dimensions(p0=1, l0=-1)
                out = FieldResult(data=data, symbol="$J\\times B$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=bx.r, z=bx.z, time=bx.time, mask=bx.mask)
            elif n in {"torque_b2", "torque_b1", "torque_p"}:
                rmat = None
                rfac = 1j * ntor
                if n == "torque_p":
                    p = read_field("p", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    out = FieldResult(data=-rfac * np.asarray(p.data), symbol="$\\tau_p$", units=parse_units(dimensions(p0=1), cgs=cgs, mks=mks), dimensions=dimensions(p0=1), r=p.r, z=p.z, time=p.time, mask=p.mask)
                else:
                    psi_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=2, cgs=cgs, mks=mks, return_meta=True)
                    psi_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=3, cgs=cgs, mks=mks, return_meta=True)
                    i_r = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=2, cgs=cgs, mks=mks, return_meta=True)
                    i_z = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=3, cgs=cgs, mks=mks, return_meta=True)
                    f_r = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=2, cgs=cgs, mks=mks, return_meta=True)
                    f_z = read_field("f", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, operation=3, cgs=cgs, mks=mks, return_meta=True)
                    rmat = radius_matrix(psi_r.r, psi_r.z) if itor == 1 else np.ones_like(np.asarray(psi_r.data), dtype=float)
                    if n == "torque_b2":
                        data = (
                            (np.conj(np.asarray(i_z.data)) * np.asarray(psi_r.data) - np.conj(np.asarray(i_r.data)) * np.asarray(psi_z.data)) / np.maximum(rmat, np.finfo(float).tiny)
                            - (np.conj(np.asarray(i_r.data)) * (rfac * np.asarray(f_r.data)) + np.conj(np.asarray(i_z.data)) * (rfac * np.asarray(f_z.data)))
                            + (np.asarray(i_z.data) * np.conj(np.asarray(psi_r.data)) - np.asarray(i_r.data) * np.conj(np.asarray(psi_z.data))) / np.maximum(rmat, np.finfo(float).tiny)
                            - (np.asarray(i_r.data) * np.conj(rfac * np.asarray(f_r.data)) + np.asarray(i_z.data) * np.conj(rfac * np.asarray(f_z.data)))
                        ) / 2.0
                    else:
                        psi0_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=2, cgs=cgs, mks=mks, return_meta=True)
                        psi0_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=3, cgs=cgs, mks=mks, return_meta=True)
                        i0_r = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=2, cgs=cgs, mks=mks, return_meta=True)
                        i0_z = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, operation=3, cgs=cgs, mks=mks, return_meta=True)
                        data = (
                            -(np.asarray(psi0_r.data) * (rfac * np.asarray(psi_r.data)) + np.asarray(psi0_z.data) * (rfac * np.asarray(psi_z.data))) / np.maximum(rmat**2, np.finfo(float).tiny)
                            + ((np.asarray(i_z.data) + rfac**2 * np.asarray(f_z.data)) * np.asarray(psi0_r.data) - (np.asarray(i_r.data) + rfac**2 * np.asarray(f_r.data)) * np.asarray(psi0_z.data)) / np.maximum(rmat, np.finfo(float).tiny)
                            + (np.asarray(i0_z.data) * np.asarray(psi_r.data) - np.asarray(i0_r.data) * np.asarray(psi_z.data)) / np.maximum(rmat, np.finfo(float).tiny)
                            - (np.asarray(i0_r.data) * rfac * np.asarray(f_r.data) + np.asarray(i0_z.data) * rfac * np.asarray(f_z.data))
                        )
                    d = dimensions(p0=1)
                    out = FieldResult(data=data, symbol="$\\tau_B$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi_r.r, z=psi_r.z, time=psi_r.time, mask=psi_r.mask)
            elif n in {"torque_mu", "torque_v1", "torque_vv2"}:
                rfac = 1j * ntor
                if n == "torque_mu":
                    mu = read_field("visc", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                    mu_c = read_field("visc_c", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                    mu_c_data = np.asarray(mu_c.data)
                    w = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    rmat = radius_matrix(mu.r, mu.z) if itor == 1 else np.ones_like(np.asarray(mu.data), dtype=float)
                    data = (
                        np.asarray(mu.data) * delstar(np.maximum(rmat**2, np.finfo(float).tiny) * np.asarray(w.data), mu.r, mu.z, toroidal=(itor == 1))
                        + np.maximum(rmat**2, np.finfo(float).tiny) * s_bracket(np.asarray(mu.data), np.asarray(w.data), mu.r, mu.z)
                    )
                    if ntor != 0:
                        u = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                        chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                        data = (
                            data
                            + 2.0 * mu_c_data * rfac**2 * np.asarray(w.data)
                            - 4.0 * rfac * (mu_c_data * dz(np.asarray(u.data), u.z))
                            - 2.0 * rfac * (np.asarray(mu.data) - mu_c_data) * delstar(np.asarray(chi.data), chi.r, chi.z, toroidal=(itor == 1))
                            + np.asarray(mu.data) * rfac * laplacian(np.asarray(chi.data), chi.r, chi.z, toroidal=(itor == 1)) / np.maximum(rmat**2, np.finfo(float).tiny)
                            + a_bracket(np.asarray(mu.data), rfac * np.asarray(chi.data), chi.r, chi.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                        )
                    base = mu
                elif n == "torque_v1":
                    w0 = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                    n0 = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                    u = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    w = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    rmat = radius_matrix(w0.r, w0.z) if itor == 1 else np.ones_like(np.asarray(w0.data), dtype=float)
                    data = -a_bracket(np.maximum(rmat**4, np.finfo(float).tiny) * np.asarray(n0.data) * np.asarray(w0.data), np.asarray(u.data), w0.r, w0.z) / np.maximum(rmat, np.finfo(float).tiny) - 2.0 * np.maximum(rmat**2, np.finfo(float).tiny) * np.asarray(n0.data) * np.asarray(w0.data) * rfac * np.asarray(w.data)
                    base = w0
                else:
                    psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                    n0 = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, linear=linear, cgs=cgs, mks=mks, return_meta=True)
                    u = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    w = read_field("omega", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, cgs=cgs, mks=mks, return_meta=True)
                    rmat = radius_matrix(psi0.r, psi0.z) if itor == 1 else np.ones_like(np.asarray(psi0.data), dtype=float)
                    data = (
                        -np.asarray(n0.data)
                        * np.conj(np.asarray(w.data))
                        * (
                            np.maximum(rmat**3, np.finfo(float).tiny) * a_bracket(np.asarray(psi0.data), np.asarray(u.data), psi0.r, psi0.z)
                            + s_bracket(np.asarray(psi0.data), np.asarray(chi.data), psi0.r, psi0.z)
                        )
                        / np.sqrt(np.maximum(s_bracket(np.asarray(psi0.data), np.asarray(psi0.data), psi0.r, psi0.z), np.finfo(float).tiny))
                    )
                    d = dimensions(p0=1, l0=1)
                    out = FieldResult(data=data, symbol="$\\tau_{vv2}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi0.r, z=psi0.z, time=psi0.time, mask=psi0.mask)
                    base = None
                if n in {"torque_mu", "torque_v1"}:
                    d = dimensions(p0=1)
                    out = FieldResult(data=data, symbol="$\\tau$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=base.r, z=base.z, time=base.time, mask=base.mask)
            elif n == "cole_ntv" or n == "cole_ntv_i":
                raise NotImplementedError("cole_ntv requires full flux-coordinate/flux-average machinery for strict IDL equivalence.")
            elif n in {"toroidal velocity", "vy"}:
                v = read_field("v", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(v.r, v.z) if itor == 1 else np.ones_like(np.asarray(v.data), dtype=float)
                ivform = int(read_parameter("ivform", filename=filename))
                data = np.asarray(v.data) / rmat if ivform == 0 else np.asarray(v.data) * rmat
                d = dimensions(v0=1)
                out = FieldResult(data=data, symbol="$u_{phi}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=v.r, z=v.z, time=v.time, mask=v.mask)
            elif n in {"toroidal velocity shear", "vzp"}:
                v = read_field("v", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(v.r, v.z) if itor == 1 else np.ones_like(np.asarray(v.data), dtype=float)
                ivform = int(read_parameter("ivform", filename=filename))
                vz = np.asarray(v.data) / rmat if ivform == 0 else np.asarray(v.data) * rmat
                denom = np.sqrt(np.maximum(s_bracket(psi.data, psi.data, psi.r, psi.z), np.finfo(float).tiny))
                data = s_bracket(vz, psi.data, v.r, v.z) / denom
                d = dimensions(t0=-1)
                out = FieldResult(data=data, symbol="$u_{phi}'$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=v.r, z=v.z, time=v.time, mask=v.mask)
            elif n in {"e_r", "er"}:
                u = read_field("phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                v = read_field("v", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                chi = read_field("chi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                b = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                pe = read_field("pe", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                eta = read_field("eta", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                eta_data = np.asarray(eta.data)

                if itor == 1:
                    rmat = radius_matrix(u.r, u.z)
                else:
                    rmat = np.ones_like(np.asarray(u.data), dtype=float)
                jphi = delstar(np.asarray(psi.data), u.r, u.z, toroidal=(itor == 1))
                db = float(read_parameter("db", filename=filename))
                psipsi = s_bracket(np.asarray(psi.data), np.asarray(psi.data), u.r, u.z)
                data = (
                    np.asarray(b.data) * s_bracket(np.asarray(u.data), np.asarray(psi.data), u.r, u.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                    + np.asarray(b.data) * a_bracket(np.asarray(chi.data), np.asarray(psi.data), u.r, u.z) / np.maximum(rmat, np.finfo(float).tiny)
                    - np.asarray(v.data) * psipsi / np.maximum(rmat**2, np.finfo(float).tiny)
                    - (db / np.maximum(np.asarray(den.data), np.finfo(float).tiny))
                    * (
                        s_bracket(np.asarray(psi.data), np.asarray(pe.data), u.r, u.z)
                        + (jphi * psipsi + np.asarray(b.data) * s_bracket(np.asarray(psi.data), np.asarray(b.data), u.r, u.z)) / np.maximum(rmat**2, np.finfo(float).tiny)
                    )
                    + eta_data * a_bracket(np.asarray(psi.data), np.asarray(b.data), u.r, u.z) / np.maximum(rmat, np.finfo(float).tiny)
                )
                data = -data / np.sqrt(np.maximum(psipsi, np.finfo(float).tiny))
                d = dimensions(potential=1, l0=-1)
                out = FieldResult(data=data, symbol="$E_r$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=u.r, z=u.z, time=u.time, mask=u.mask)
            elif n == "curl_e" or n == "curl_e_i":
                pert = "_i" if n.endswith("_i") else ""
                if not zvector and not yvector:
                    ephi_z = read_field(f"e_phi{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=3, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    ez_phi = read_field(f"e_r{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=11, phi=phi, cgs=cgs, mks=mks, return_meta=False)
                if not rvector and not yvector:
                    ephi_r = read_field(f"e_phi{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=2, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    er_phi = read_field(f"e_r{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=11, phi=phi, cgs=cgs, mks=mks, return_meta=False)
                if not rvector and not zvector:
                    ez_r = read_field(f"e_z{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=2, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    er_z = read_field(f"e_r{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=3, phi=phi, cgs=cgs, mks=mks, return_meta=True)

                if rvector:
                    data = np.asarray(ephi_z.data) - np.asarray(ez_phi)
                    base = ephi_z
                elif zvector:
                    data = np.asarray(er_phi) - np.asarray(ephi_r.data)
                    base = ephi_r
                elif yvector:
                    data = np.asarray(ez_r.data) - np.asarray(er_z.data)
                    base = ez_r
                else:
                    cx = np.asarray(ephi_z.data) - np.asarray(ez_phi)
                    cy = np.asarray(ez_r.data) - np.asarray(er_z.data)
                    cz = np.asarray(er_phi) - np.asarray(ephi_r.data)
                    data = np.sqrt(np.maximum(np.abs(cx) ** 2 + np.abs(cy) ** 2 + np.abs(cz) ** 2, 0.0))
                    base = ephi_z
                d = dimensions(potential=1, l0=-2)
                out = FieldResult(data=data, symbol="$\\nabla\\times E$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=base.r, z=base.z, time=base.time, mask=base.mask)
            elif n == "div_e" or n == "div_e_i":
                pert = "_i" if n.endswith("_i") else ""
                er_r = read_field(f"e_r{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=2, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ephi_phi = read_field(f"e_phi{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=11, phi=phi, cgs=cgs, mks=mks, return_meta=False)
                ez_z = read_field(f"e_z{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, operation=3, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(potential=1, l0=-2)
                data = np.asarray(er_r.data) + np.asarray(ephi_phi) + np.asarray(ez_z.data)
                out = FieldResult(data=data, symbol="$\\nabla\\cdot E$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=er_r.r, z=er_r.z, time=er_r.time, mask=er_r.mask)
            elif n in {"exb_x"}:
                ey = read_field("e_phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ez = read_field("e_z", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bx = read_field("bx", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                by = read_field("by", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bz = read_field("bz", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                b2 = np.asarray(bx.data) ** 2 + np.asarray(by.data) ** 2 + np.asarray(bz.data) ** 2
                d = dimensions(v0=1)
                out = FieldResult(data=(np.asarray(ey.data) * np.asarray(bz.data) - np.asarray(ez.data) * np.asarray(by.data)) / np.maximum(b2, np.finfo(float).tiny), symbol="$V_{ExB,R}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=bx.r, z=bx.z, time=bx.time, mask=bx.mask)
            elif n in {"exb_y"}:
                ex = read_field("e_r", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ez = read_field("e_z", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bx = read_field("bx", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                by = read_field("by", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bz = read_field("bz", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                b2 = np.asarray(bx.data) ** 2 + np.asarray(by.data) ** 2 + np.asarray(bz.data) ** 2
                d = dimensions(v0=1)
                out = FieldResult(data=(np.asarray(ez.data) * np.asarray(bx.data) - np.asarray(ex.data) * np.asarray(bz.data)) / np.maximum(b2, np.finfo(float).tiny), symbol="$V_{ExB,phi}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=bx.r, z=bx.z, time=bx.time, mask=bx.mask)
            elif n in {"exb_z"}:
                ex = read_field("e_r", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ey = read_field("e_phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bx = read_field("bx", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                by = read_field("by", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bz = read_field("bz", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                b2 = np.asarray(bx.data) ** 2 + np.asarray(by.data) ** 2 + np.asarray(bz.data) ** 2
                d = dimensions(v0=1)
                out = FieldResult(data=(np.asarray(ex.data) * np.asarray(by.data) - np.asarray(ey.data) * np.asarray(bx.data)) / np.maximum(b2, np.finfo(float).tiny), symbol="$V_{ExB,Z}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=bx.r, z=bx.z, time=bx.time, mask=bx.mask)
            elif n in {"electron density", "ne"}:
                print("WARNING: calculating ne using ne = z_ion * ni")
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                zion = float(read_parameter("z_ion", filename=filename))
                d = dimensions(n0=1)
                out = FieldResult(data=np.asarray(den.data) * zion, symbol="$n_e$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=den.r, z=den.z, time=den.time, mask=den.mask)
            elif n in {"angular momentum", "lz"}:
                v = read_field("v", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                den = read_field("den", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(n0=1, v0=1, l0=1)
                out = FieldResult(data=np.asarray(den.data) * np.asarray(v.data), symbol="$L_{phi}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=v.r, z=v.z, time=v.time, mask=v.mask)
            elif n in {"minor radius", "r"}:
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                lcfs = read_lcfs(filename=filename, slice=slice_idx, cgs=cgs, mks=mks, return_meta=True)
                rr = np.repeat(psi.r[:, None], psi.z.size, axis=1)
                zz = np.repeat(psi.z[None, :], psi.r.size, axis=0)
                d = dimensions(l0=1)
                out = FieldResult(data=np.sqrt((rr - lcfs.axis[0]) ** 2 + (zz - lcfs.axis[1]) ** 2), symbol="$r$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi.r, z=psi.z, time=psi.time, mask=psi.mask)
            elif n in {"q_cyl"}:
                rho = read_field("r", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                bp = read_field("bp", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                bt = read_field("bt", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(rho.r, rho.z)
                else:
                    rmat = np.ones_like(np.asarray(rho.data), dtype=float) * float(read_parameter("rzero", filename=filename))
                out = FieldResult(data=np.asarray(rho.data) * np.asarray(bt.data) / np.maximum(rmat * np.asarray(bp.data), np.finfo(float).tiny), symbol="$q_{cyl}$", units="", dimensions=dimensions(), r=rho.r, z=rho.z, time=rho.time, mask=rho.mask)
            elif n in {"poloidal angle", "theta"}:
                psi = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                lc = read_lcfs(filename=filename, slice=slice_idx, cgs=cgs, mks=mks, return_meta=True)
                rr = np.repeat(psi.r[:, None], psi.z.size, axis=1)
                zz = np.repeat(psi.z[None, :], psi.r.size, axis=0)
                data = np.arctan2(zz - float(lc.axis[1]), rr - float(lc.axis[0]))
                d = dimensions()
                out = FieldResult(data=data, symbol="$\\theta$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi.r, z=psi.z, time=psi.time, mask=psi.mask)
            elif n == "pest angle":
                raise NotImplementedError("pest angle requires flux-coordinate mapping (flux_coord_field) for strict IDL equivalence.")
            elif n in {"field strength", "b"}:
                b2 = read_field("b2", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(b0=1)
                out = FieldResult(data=np.sqrt(np.maximum(np.asarray(b2.data), 0.0)), symbol="$B$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=b2.r, z=b2.z, time=b2.time, mask=b2.mask)
            elif n in {"field energy"}:
                b2 = read_field("b2", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(energy=1)
                out = FieldResult(data=np.asarray(b2.data), symbol="$B^2$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=b2.r, z=b2.z, time=b2.time, mask=b2.mask)
            elif n in {"poloidal field strength"}:
                out = read_field("bp", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
            elif n in {"toroidal field strength"}:
                out = read_field("bt", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
            elif n in {"toroidal field", "by", "bt"}:
                I = read_field("I", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, complex=False, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                rmat = radius_matrix(I.r, I.z) if itor == 1 else np.ones_like(np.asarray(I.data), dtype=float)
                d = dimensions(b0=1)
                data = np.asarray(I.data) / rmat
                out = FieldResult(data=data, symbol="$B_{phi}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=I.r, z=I.z, time=I.time, mask=I.mask)
            elif n == "bx":
                psi_z = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                fp_r = read_field("fp", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=False)
                rmat = radius_matrix(psi_z.r, psi_z.z) if itor == 1 else np.ones_like(np.asarray(psi_z.data), dtype=float)
                d = dimensions(b0=1)
                data = -np.asarray(psi_z.data) / rmat - np.asarray(fp_r)
                out = FieldResult(data=data, symbol="$B_R$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi_z.r, z=psi_z.z, time=psi_z.time, mask=psi_z.mask)
            elif n == "bz":
                psi_r = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                fp_z = read_field("fp", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=False)
                rmat = radius_matrix(psi_r.r, psi_r.z) if itor == 1 else np.ones_like(np.asarray(psi_r.data), dtype=float)
                d = dimensions(b0=1)
                data = np.asarray(psi_r.data) / rmat - np.asarray(fp_z)
                out = FieldResult(data=data, symbol="$B_Z$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi_r.r, z=psi_r.z, time=psi_r.time, mask=psi_r.mask)
            elif n == "bz_plasma":
                psi_r = read_field("psi_plasma", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=2, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                fz = read_field("f_plasma", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, operation=3, linear=linear, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                fp_z = 1j * ntor * np.asarray(fz.data)
                rmat = radius_matrix(psi_r.r, psi_r.z) if itor == 1 else np.ones_like(np.asarray(psi_r.data), dtype=float)
                d = dimensions(b0=1)
                data = np.asarray(psi_r.data) / rmat - np.asarray(fp_z)
                out = FieldResult(data=data, symbol="$B_Z$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi_r.r, z=psi_r.z, time=psi_r.time, mask=psi_r.mask)
            elif n == "bp":
                bx = read_field("bx", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bz = read_field("bz", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(b0=1)
                out = FieldResult(data=np.sqrt(np.asarray(bx.data) ** 2 + np.asarray(bz.data) ** 2), symbol="$B_p$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=bx.r, z=bx.z, time=bx.time, mask=bx.mask)
            elif n == "b2":
                bx = read_field("bx", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                by = read_field("by", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bz = read_field("bz", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=linear, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                d = dimensions(b0=2)
                out = FieldResult(data=np.asarray(bx.data) ** 2 + np.asarray(by.data) ** 2 + np.asarray(bz.data) ** 2, symbol="$B^2$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=bx.r, z=bx.z, time=bx.time, mask=bx.mask)
            elif n == "db_b" or n == "db_b_i":
                psi0 = read_field("psi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                i0 = read_field("i", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                te0 = read_field("te", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, equilibrium=True, cgs=cgs, mks=mks, return_meta=True)
                pert = "_i" if n.endswith("_i") else ""
                psi = read_field(f"psi{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, cgs=cgs, mks=mks, return_meta=True)
                i_f = read_field(f"i{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, cgs=cgs, mks=mks, return_meta=True)
                f_f = read_field(f"f{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, cgs=cgs, mks=mks, return_meta=True)
                te1 = read_field(f"te{pert}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, cgs=cgs, mks=mks, return_meta=True)
                if itor == 1:
                    rmat = radius_matrix(psi0.r, psi0.z)
                else:
                    rmat = np.ones_like(np.asarray(psi0.data), dtype=float)
                gradpsi = np.sqrt(np.maximum(s_bracket(np.asarray(psi0.data), np.asarray(psi0.data), psi0.r, psi0.z), np.finfo(float).tiny))
                tprime = np.abs(s_bracket(np.asarray(te0.data), np.asarray(psi0.data), psi0.r, psi0.z))
                xi = -np.asarray(te1.data) / np.maximum(tprime, np.finfo(float).tiny) * gradpsi
                lc = read_lcfs(filename=filename, slice=slice_idx, cgs=cgs, mks=mks, return_meta=True)
                psis = float(lc.psilim)
                flux0 = float(lc.flux0)
                if psis < flux0:
                    xi = np.where((np.abs(tprime) < 1e-6) | (np.asarray(psi0.data) < psis), 0.0, xi)
                else:
                    xi = np.where((np.abs(tprime) < 1e-6) | (np.asarray(psi0.data) > psis), 0.0, xi)
                b0 = np.sqrt(np.maximum(gradpsi**2 + np.asarray(i0.data) ** 2, np.finfo(float).tiny)) / np.maximum(rmat, np.finfo(float).tiny)
                b1b0 = (
                    s_bracket(np.asarray(psi0.data), np.asarray(psi.data), psi0.r, psi0.z) / np.maximum(rmat**2, np.finfo(float).tiny)
                    + np.asarray(i0.data) * np.asarray(i_f.data) / np.maximum(rmat**2, np.finfo(float).tiny)
                    - 1j * float(ntor) * a_bracket(np.asarray(f_f.data), np.asarray(psi0.data), psi0.r, psi0.z) / np.maximum(rmat, np.finfo(float).tiny)
                ) / np.maximum(b0, np.finfo(float).tiny)
                b1 = b1b0 + xi * s_bracket(b0, np.asarray(psi0.data), psi0.r, psi0.z) / np.maximum(gradpsi, np.finfo(float).tiny)
                b2 = b1 * np.conj(b1)
                d = dimensions()
                out = FieldResult(data=np.asarray(b2) / np.maximum(b0**2, np.finfo(float).tiny), symbol="$\\delta B^2/B^2$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=psi0.r, z=psi0.z, time=psi0.time, mask=psi0.mask)
            elif n in {"pitch_angle_z", "pitch_angle_x", "dpitch_angle_z", "dpitch_angle_x"}:
                axis = "z" if n.endswith("_z") else "x"
                if linear:
                    by0 = read_field("by", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    b0 = read_field(f"b{axis}", filename=filename, timeslices=-1, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    by1 = read_field("by", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    b1 = read_field(f"b{axis}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, linear=True, complex=complex, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                    if n.startswith("dpitch"):
                        data = np.asarray(b1.data) / np.maximum(np.asarray(b0.data), np.finfo(float).tiny) - np.asarray(by1.data) / np.maximum(np.asarray(by0.data), np.finfo(float).tiny)
                    else:
                        data = np.asarray(b1.data) / np.maximum(np.asarray(by0.data), np.finfo(float).tiny) - np.asarray(by1.data) * np.asarray(b0.data) / np.maximum(np.asarray(by0.data), np.finfo(float).tiny)
                    rbase, zbase, tbase, mbase = by0.r, by0.z, by0.time, by0.mask
                else:
                    by = read_field("by", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    b = read_field(f"b{axis}", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, cgs=cgs, mks=mks, return_meta=True)
                    data = np.asarray(b.data) / np.maximum(np.asarray(by.data), np.finfo(float).tiny)
                    rbase, zbase, tbase, mbase = by.r, by.z, by.time, by.mask
                d = dimensions()
                sym = "$B_Z/B_{phi}$" if axis == "z" else "$B_R/B_{phi}$"
                out = FieldResult(data=data, symbol=sym, units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=rbase, z=zbase, time=tbase, mask=mbase)
            elif n == "epar":
                bx = read_field("bx", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                by = read_field("by", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bz = read_field("bz", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ex = read_field("e_r", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ey = read_field("e_phi", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                ez = read_field("e_z", filename=filename, timeslices=slice_idx, points=points, xrange=xrange, yrange=yrange, phi=phi, cgs=cgs, mks=mks, return_meta=True)
                bmag = np.sqrt(np.asarray(bx.data) ** 2 + np.asarray(by.data) ** 2 + np.asarray(bz.data) ** 2)
                d = dimensions(potential=1, l0=-1)
                data = (np.asarray(ex.data) * np.asarray(bx.data) + np.asarray(ey.data) * np.asarray(by.data) + np.asarray(ez.data) * np.asarray(bz.data)) / np.maximum(bmag, np.finfo(float).tiny)
                out = FieldResult(data=data, symbol="$E_{||}$", units=parse_units(d, cgs=cgs, mks=mks), dimensions=d, r=bx.r, z=bx.z, time=bx.time, mask=bx.mask)
            else:
                if primitive_error is not None:
                    raise primitive_error
                raise KeyError(f"Field '{name}' not found as primitive and composite field is not implemented.")

    # scale and post-processing
    if fac is not None:
        print(f"applying factor = {fac}")
        out.data = np.asarray(out.data) * float(fac)
    if linear and linfac is not None:
        print(f"scaling data by {linfac}")
        out.data = np.asarray(out.data) * float(linfac)
    print(f"converting units, mks, cgs= {bool(mks)} {bool(cgs)}")
    if abs:
        print("Taking absolute value of data")
        out.data = np.abs(np.asarray(out.data))
    if phase:
        print("Taking phase of data")
        out.data = np.angle(np.asarray(out.data))
    if real:
        out.data = np.real(np.asarray(out.data))
    if imaginary:
        out.data = np.imag(np.asarray(out.data))

    if slice_idx >= 0:
        out.time = float(get_slice_time(filename=filename, slice=slice_idx, cgs=cgs, mks=mks)[0])
    else:
        out.time = 0.0

    if return_meta or rvector or zvector:
        if symbol is not None:
            out.symbol = str(symbol)
        if units is not None:
            out.units = str(units)
        if rvector and zvector:
            return out.r, out.z, out.data
        if rvector:
            return out.r, out.data
        if zvector:
            return out.z, out.data
        print("Done reading field")
        print("**********************************************************")
        return out
    print("Done reading field")
    print("**********************************************************")
    return out.data
