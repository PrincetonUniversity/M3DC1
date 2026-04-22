from __future__ import annotations

from pathlib import Path
from typing import Sequence

import h5py
import numpy as np
from .hdf5_file_test import hdf5_file_test


def _read_parameter_one(name: str, filename: str | Path) -> np.ndarray | float:
    with h5py.File(str(filename), "r") as h5:
        for key, val in h5.attrs.items():
            if key.lower() == name.lower():
                arr = np.asarray(val)
                if arr.size == 1:
                    return float(arr.reshape(-1)[0])
                return arr
    return 0.0


def read_parameter(
    name: str,
    filename: str | Path | Sequence[str | Path] = "C1.h5",
    pr: bool = False,
    mks: bool = False,
    cgs: bool = False,
):
    """
    Python port of read_parameter.pro.
    Supports:
      - scalar or list/tuple filename
      - case-insensitive root attribute lookup
      - optional print (pr)
      - optional cgs/mks conversion
    """
    if isinstance(filename, (list, tuple)):
        return np.asarray(
            [read_parameter(name, f, pr=pr, mks=mks, cgs=cgs) for f in filename], dtype=float
        )

    if not hdf5_file_test(filename):
        return 0.0

    attr = _read_parameter_one(name, filename)

    if cgs or mks:
        from .convert_units import convert_units
        from .field_data import field_data
        from .get_normalizations import get_normalizations

        itor = int(read_parameter("itor", filename=filename))
        _sym, d0 = field_data(name, itor=itor, filename=filename)
        b0, n0, l0, mi = get_normalizations(filename=filename)
        attr = convert_units(
            attr, d0, b0=b0, n0=n0, l0=l0, mi=mi, filename=filename, cgs=cgs, mks=mks
        )

    if pr:
        print(f"{name} = {attr}")

    return attr
