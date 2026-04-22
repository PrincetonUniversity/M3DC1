from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from .compensate_renorm import compensate_renorm as apply_compensate_renorm
from .convert_units import convert_units
from .dimensions import dimensions
from .hdf5_file_test import hdf5_file_test
from .make_units import make_units
from .power_spectrum import power_spectrum as calc_power_spectrum
from .read_gamma import read_gamma
from .read_parameter import read_parameter
from .read_scalar import read_scalar


@dataclass
class SignalResult:
    tdata: np.ndarray
    data: np.ndarray
    xtitle: str
    ylabel: str


def _find_dataset_recursive(group: h5py.Group, pred) -> np.ndarray | None:
    for key, item in group.items():
        if isinstance(item, h5py.Dataset) and pred(key, item):
            return np.asarray(item[()])
        if isinstance(item, h5py.Group):
            out = _find_dataset_recursive(item, pred)
            if out is not None:
                return out
    return None


def _read_signal_matrix(h5: h5py.File, signal_dir: str) -> np.ndarray:
    if signal_dir not in h5:
        raise KeyError(f"Missing group '{signal_dir}' in HDF5 file.")

    grp = h5[signal_dir]
    if not isinstance(grp, h5py.Group):
        raise TypeError(f"'{signal_dir}' exists but is not an HDF5 group.")

    path_candidates = ["value/_data", "value", "_data", "data"]
    for rel in path_candidates:
        full = f"{signal_dir}/{rel}"
        if full in h5 and isinstance(h5[full], h5py.Dataset):
            out = np.asarray(h5[full][()])
            break
    else:
        out = _find_dataset_recursive(grp, lambda _k, _ds: True)
        if out is None:
            raise KeyError(f"No dataset found under group '{signal_dir}'.")

    out = np.asarray(out, dtype=float)
    if out.ndim == 1:
        out = out[np.newaxis, :]
    if out.ndim != 2:
        raise ValueError(
            f"Expected '{signal_dir}' data to be 2D [n_signal, n_time], got shape {out.shape}."
        )
    return out


def read_signals(
    signal: str,
    filename: str | Path = "C1.h5",
    deriv: bool = False,
    power_spectrum: bool = False,
    compensate_renorm: bool = False,
    scale: bool = False,
    cgs: bool = False,
    mks: bool = False,
    return_meta: bool = False,
):
    """
    Read and process magnetic probe / flux-loop signals.

    Returns:
    - data (np.ndarray) by default
    - SignalResult when return_meta=True
    """
    if not hdf5_file_test(filename):
        raise ValueError(f"Invalid HDF5 file: {filename}")

    signal_l = signal.lower()
    if signal_l == "mag_probes":
        sigdir = "mag_probes"
        count_attr = "imag_probes"
        dim = dimensions(b0=1)
        ylabel = (
            f"dB/dt [{make_units(b0=1, t0=-1, cgs=cgs, mks=mks)}]"
            if deriv
            else f"B.n [{make_units(b0=1, cgs=cgs, mks=mks)}]"
        )
    elif signal_l == "flux_loops":
        sigdir = "flux_loops"
        count_attr = "iflux_loops"
        dim = dimensions(b0=1, l0=2)
        ylabel = (
            f"d/dt Flux [{make_units(potential=1, cgs=cgs, mks=mks)}]"
            if deriv
            else f"Flux [{make_units(b0=1, l0=2, cgs=cgs, mks=mks)}]"
        )
    else:
        raise ValueError(f"Unsupported signal '{signal}'.")

    with h5py.File(str(filename), "r") as h5:
        data = _read_signal_matrix(h5, sigdir)
        tmeta = read_scalar("time", filename=filename, cgs=cgs, mks=mks, return_meta=True)
        t = np.asarray(tmeta.data, dtype=float).reshape(-1)
        tu = tmeta.units

        if data.shape[1] != t.size:
            if data.shape[0] == t.size:
                data = data.T
            else:
                raise ValueError(
                    f"Time length mismatch: data has {data.shape[1]} samples but time has {t.size}."
                )

        n_attr = float(read_parameter(count_attr, filename=filename))
        if n_attr == 0:
            raise ValueError(f"No data for signal: {signal}")

        data = convert_units(data, dim, filename=filename, cgs=cgs, mks=mks)

        if compensate_renorm:
            for i in range(data.shape[0]):
                data[i, :] = apply_compensate_renorm(data[i, :])

        if scale:
            gamma = float(read_gamma(filename, cgs=cgs, mks=mks)[0])
            scale_vec = np.exp(t * gamma)
            for i in range(data.shape[0]):
                end = data[i, -1] if data.shape[1] else 1.0
                if end != 0:
                    data[i, :] = data[i, :] / scale_vec / end

    if deriv:
        data = np.gradient(data, t, axis=1)

    if power_spectrum:
        out = np.zeros_like(data, dtype=float)
        tdata = None
        for i in range(data.shape[0]):
            out[i, :], freq = calc_power_spectrum(data[i, :], float(np.max(t)))
            tdata = freq
        data = out
        assert tdata is not None
        xtitle = "Frequency"
    else:
        tdata = t
        xtitle = f"Time ({tu})" if tu else "Time"

    if return_meta:
        return SignalResult(
            tdata=np.asarray(tdata),
            data=np.asarray(data),
            xtitle=xtitle,
            ylabel=ylabel,
        )
    return np.asarray(data)
