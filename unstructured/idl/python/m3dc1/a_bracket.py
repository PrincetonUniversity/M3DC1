from __future__ import annotations

import numpy as np

from .dx import dx
from .dz import dz


def a_bracket(a: np.ndarray, b: np.ndarray, x, z) -> np.ndarray:
    """Port of a_bracket.pro: antisymmetric bracket in (x,z)."""
    aa = np.asarray(a)
    bb = np.asarray(b)

    if aa.ndim == 3 and bb.ndim == 2:
        out = np.zeros_like(aa)
        for i in range(aa.shape[0]):
            out[i] = a_bracket(aa[i], bb, x, z)
        return out
    if bb.ndim == 3 and aa.ndim == 2:
        out = np.zeros_like(bb)
        for i in range(bb.shape[0]):
            out[i] = a_bracket(aa, bb[i], x, z)
        return out
    if aa.ndim == 2 and bb.ndim == 2:
        return -dx(aa, x) * dz(bb, z) + dz(aa, z) * dx(bb, x)
    raise ValueError("Error: sizes do not conform for a_bracket")
