from __future__ import annotations

import numpy as np

from .find_nulls import find_nulls
from .read_nulls import read_nulls


def nulls(psi=None, x=None, z=None, *, axis=False, xpoints=False, filename="C1.h5", slice=0, **kwargs):
    try:
        ax, xp = read_nulls(filename=filename, slice=slice, **kwargs)
    except Exception:
        ax, xp = find_nulls(psi, x, z, filename=filename, slice=slice, **kwargs)
    if axis and xpoints:
        return ax, xp
    if axis:
        return ax
    if xpoints:
        return xp
    return ax, xp
