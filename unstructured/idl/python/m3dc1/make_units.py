from __future__ import annotations

from .dimensions import dimensions
from .parse_units import parse_units

def make_units(cgs: bool = False, mks: bool = False, **kwargs: int) -> str:
    """
    Python port of make_units.pro.
    Example: make_units(b0=1, t0=-1, mks=True)
    """
    return parse_units(dimensions(**kwargs), cgs=cgs, mks=mks)
