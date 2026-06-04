from __future__ import annotations

import numpy as np

from .read_nulls import read_nulls


def find_nulls(psi=None, x=None, z=None, *, filename="C1.h5", slice=0, **kwargs):
    """
    Placeholder-compatible Python implementation. Uses scalar-tracked nulls.
    """
    axis, xpoint = read_nulls(filename=filename, slice=slice, **kwargs)
    return np.asarray(axis, dtype=float), np.asarray(xpoint, dtype=float)

