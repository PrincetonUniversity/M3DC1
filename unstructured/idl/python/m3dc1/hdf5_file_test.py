from __future__ import annotations

from pathlib import Path

import h5py


def hdf5_file_test(filename: str | Path) -> bool:
    """
    Python port of hdf5_file_test.pro.
    """
    path = str(filename)
    try:
        ok = h5py.is_hdf5(path)
    except OSError:
        ok = False

    if not ok:
        print(f"Error: {path} is not a valid file.")
    return ok
