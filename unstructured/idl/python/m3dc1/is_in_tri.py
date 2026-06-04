from __future__ import annotations

import numpy as np


def is_in_tri(localp: np.ndarray, a: float, b: float, c: float) -> bool:
    """Return True if local point is inside the element triangle."""
    lp = np.asarray(localp, dtype=float)
    small = (a + b + c) * 1e-4

    if lp[1] < -small:
        return False
    if lp[1] > c + small:
        return False

    x = 1.0 - lp[1] / c
    if lp[0] < -b * x - small:
        return False
    if lp[0] > a * x + small:
        return False

    return True
