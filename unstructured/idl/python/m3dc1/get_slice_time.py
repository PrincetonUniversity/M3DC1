from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from .convert_units import convert_units
from .dimensions import dimensions
from .hdf5_file_test import hdf5_file_test
from .read_parameter import read_parameter
from .time_name import time_name


def get_slice_time(
    filename: str | Path = "C1.h5",
    slice=None,
    cgs: bool = False,
    mks: bool = False,
) -> np.ndarray:
    """Return physical time(s) for one or more time slices."""
    if not hdf5_file_test(filename):
        raise ValueError(f"Invalid HDF5 file: {filename}")

    if slice is None:
        nt = int(read_parameter("ntime", filename=filename))
        timeslices = np.arange(nt, dtype=int)
    else:
        timeslices = np.atleast_1d(np.asarray(slice, dtype=int))

    out = np.zeros(timeslices.size, dtype=float)
    with h5py.File(str(filename), "r") as h5:
        for i, s in enumerate(timeslices):
            gname = time_name(int(s))
            if gname not in h5:
                raise KeyError(f"Missing group '{gname}' in file.")
            grp = h5[gname]
            if "time" in grp.attrs:
                tt = float(grp.attrs["time"])
            else:
                tt = 0.0
            out[i] = float(
                np.asarray(
                    convert_units(np.asarray([tt], dtype=float), dimensions(t0=1), filename=filename, cgs=cgs, mks=mks)
                ).reshape(-1)[0]
            )

    return out
