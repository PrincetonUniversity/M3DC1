from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from .convert_units import convert_units
from .dimensions import dimensions
from .get_normalizations import get_normalizations
from .hdf5_file_test import hdf5_file_test
from .read_scalar import read_scalar


@dataclass
class HmnResult:
    data: np.ndarray
    time: np.ndarray
    title: str
    xtitle: str
    ytitle: str
    units: str
    magnetic: bool


def _find_group_case_insensitive(h5: h5py.File, name: str):
    n = name.lower()
    for k in h5.keys():
        if k.lower() == n:
            return h5[k]
    return None


def _first_dataset(obj) -> np.ndarray | None:
    if isinstance(obj, h5py.Dataset):
        return np.asarray(obj[()])
    if isinstance(obj, h5py.Group):
        for k in ("_data", "data", "value"):
            if k in obj and isinstance(obj[k], h5py.Dataset):
                return np.asarray(obj[k][()])
        for k in ("KEHARMONICS", "BHARMONICS"):
            if k in obj:
                out = _first_dataset(obj[k])
                if out is not None:
                    return out
        for k in obj.keys():
            out = _first_dataset(obj[k])
            if out is not None:
                return out
    return None


def _read_harmonics(filename: str | Path, magnetic: bool) -> np.ndarray:
    key = "bharmonics" if magnetic else "keharmonics"
    with h5py.File(str(filename), "r") as h5:
        obj = _find_group_case_insensitive(h5, key)
        if obj is None:
            raise KeyError(f"Missing '{key}' in file.")
        arr = _first_dataset(obj)
        if arr is None:
            raise KeyError(f"Could not parse dataset for '{key}'.")
    out = np.asarray(arr, dtype=float)
    if out.ndim != 2:
        raise ValueError(f"Expected 2D harmonics array for '{key}', got {out.shape}.")
    return out


def _orient_harmonics(h: np.ndarray, nt: int) -> np.ndarray:
    arr = np.asarray(h, dtype=float)
    if arr.shape[1] == nt:
        return arr
    if arr.shape[0] == nt:
        return arr.T
    d0 = abs(arr.shape[0] - nt)
    d1 = abs(arr.shape[1] - nt)
    return arr if d1 <= d0 else arr.T


def read_hmn(
    *,
    filename: str | Path = "C1.h5",
    maxn: int | None = None,
    growth: bool = False,
    outfile: str | Path | None = None,
    me: bool = False,
    cgs: bool = False,
    mks: bool = False,
    return_meta: bool = False,
):
    """Read kinetic or magnetic harmonic time traces for plot_hmn()."""
    if not hdf5_file_test(filename):
        return HmnResult(
            data=np.asarray([], dtype=float),
            time=np.asarray([], dtype=float),
            title="",
            xtitle="t",
            ytitle="",
            units="",
            magnetic=bool(me),
        ) if return_meta else np.asarray([], dtype=float)

    magnetic = bool(me)
    title = "Magnetic Energy" if magnetic else "Kinetic Energy"
    data = _read_harmonics(filename, magnetic=magnetic)

    tmeta = read_scalar("time", filename=filename, cgs=cgs, mks=mks, return_meta=True)
    time = np.asarray(tmeta.data, dtype=float).reshape(-1)
    time = time[np.isfinite(time)]
    if time.size == 0:
        return HmnResult(
            data=np.asarray([], dtype=float),
            time=np.asarray([], dtype=float),
            title=title,
            xtitle="t",
            ytitle="Growth Rate" if growth else title,
            units="",
            magnetic=magnetic,
        ) if return_meta else np.asarray([], dtype=float)

    data = _orient_harmonics(data, nt=time.size)
    ntimes = min(int(data.shape[1]), int(time.size))
    data = data[:, :ntimes]
    time = time[:ntimes]

    d = dimensions(energy=1)
    b0, n0, l0, mi = get_normalizations(filename=filename)
    data = convert_units(data, d, b0=b0, n0=n0, l0=l0, mi=mi, filename=filename)

    if maxn is None:
        maxn = int(data.shape[0])
    maxn = int(max(1, min(maxn, int(data.shape[0]))))
    data = np.asarray(data[:maxn, :], dtype=float)

    if outfile is not None:
        out = np.column_stack([time.reshape(-1), data.T])
        np.savetxt(str(outfile), out, fmt="%16.6e")

    if growth:
        tiny = np.finfo(float).tiny
        out = np.zeros_like(data, dtype=float)
        for n in range(maxn):
            out[n, :] = np.gradient(np.log(np.maximum(np.abs(data[n, :]), tiny)), time)
        data = out

    result = HmnResult(
        data=data,
        time=time,
        title=title,
        xtitle=f"t ({tmeta.units})" if tmeta.units else "t",
        ytitle="Growth Rate" if growth else title,
        units="",
        magnetic=magnetic,
    )
    if return_meta:
        return result
    return np.asarray(result.data, dtype=float)
