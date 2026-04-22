from __future__ import annotations

from pathlib import Path

import numpy as np


def _read_ascii(path: Path) -> np.ndarray:
    arr = np.asarray(np.loadtxt(str(path), dtype=float))
    if arr.ndim == 0:
        return arr.reshape(1, 1)
    if arr.ndim == 1:
        return arr.reshape(-1, 1)
    return np.asarray(arr, dtype=float)


def read_coil_data(*, directory=".", rmp=False):
    d = Path(directory)
    coil_file = d / ("rmp_coil.dat" if rmp else "coil.dat")
    curr_file = d / ("rmp_current.dat" if rmp else "current.dat")
    if not coil_file.exists() or not curr_file.exists():
        return None
    coil = _read_ascii(coil_file)
    curr = _read_ascii(curr_file)
    n = min(coil.shape[0], curr.shape[0])
    if n == 0:
        return None
    m = int(curr.shape[1])
    print(f"n, m = {n}, {m}")

    out = np.zeros((10, n), dtype=float)
    out[0, :] = curr[:n, 0]
    out[1, :] = curr[:n, 1] if m >= 2 else curr[:n, 0]

    # Match IDL read_coil_data: the newer "field01" layout stores the
    # geometry in columns 4:9, while the older layout uses all coil columns.
    if coil.shape[1] >= 9:
        out[2:8, :] = coil[:n, 3:9].T
        out[8, :] = 1.0
        out[9, :] = 1.0
    else:
        l = min(int(coil.shape[1]), 8)
        if l > 0:
            out[2 : 2 + l, :] = coil[:n, :l].T

    print(out)
    return out
