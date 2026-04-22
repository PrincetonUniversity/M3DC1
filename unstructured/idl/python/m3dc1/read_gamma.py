from __future__ import annotations

from pathlib import Path

import numpy as np
from .read_scalar import read_scalar


def read_gamma(
    filename: str | Path = "C1.h5",
    sample_fraction: float | None = None,
    cgs: bool = False,
    mks: bool = False,
) -> np.ndarray:
    """
    Python port of IDL read_gamma.pro (single-file form).
    Returns an array for parity with IDL behavior.
    """
    t = np.asarray(read_scalar("time", filename=filename, cgs=cgs, mks=mks), dtype=float).reshape(-1)
    ke = np.asarray(read_scalar("ke", filename=filename, cgs=cgs, mks=mks), dtype=float).reshape(-1)

    if t.size < 2 or ke.size != t.size:
        return np.array([0.0], dtype=float)

    pt = 10
    if sample_fraction is not None:
        pt = int(t.size * sample_fraction)
    if t.size < pt or pt <= 0:
        return np.array([0.0], dtype=float)

    ke_safe = np.maximum(ke, np.finfo(float).tiny)
    g = np.gradient(np.log(ke_safe), t) / 2.0
    gamma = float(np.median(g[-pt:]))
    return np.array([gamma], dtype=float)
