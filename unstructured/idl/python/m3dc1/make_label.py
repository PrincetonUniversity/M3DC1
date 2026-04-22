from __future__ import annotations

import numpy as np

from .dimensions import dimensions
from .parse_units import parse_units


def make_label(s: str, d: np.ndarray | None = None, cgs: bool = False, mks: bool = False, **kwargs: int) -> str:
    """
    Python port of make_label.pro using plain text labels.
    """
    dim = np.asarray(d) if d is not None else dimensions(**kwargs)
    if np.max(np.abs(dim)) == 0:
        return s
    u = parse_units(dim, cgs=cgs, mks=mks)
    return f"{s} ({u})" if u else s
