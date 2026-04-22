from __future__ import annotations

import numpy as np

from .dx import dx
from .dz import dz


def s_bracket(a: np.ndarray, b: np.ndarray, x, z) -> np.ndarray:
    """Port of s_bracket.pro: grad(a)·grad(b) in (x,z)."""
    aa = np.asarray(a)
    bb = np.asarray(b)

    if aa.ndim == 3 and bb.ndim == 2:
        out = np.zeros_like(aa)
        for i in range(aa.shape[0]):
            out[i] = s_bracket(aa[i], bb, x, z)
        return out
    if bb.ndim == 3 and aa.ndim == 2:
        out = np.zeros_like(bb)
        for i in range(bb.shape[0]):
            out[i] = s_bracket(aa, bb[i], x, z)
        return out
    if aa.ndim == 2 and bb.ndim == 2:
        return dx(aa, x) * dx(bb, x) + dz(aa, z) * dz(bb, z)
    raise ValueError("Error: sizes do not conform for s_bracket")
