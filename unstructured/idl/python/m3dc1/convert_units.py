from __future__ import annotations

from pathlib import Path

import numpy as np

from .get_normalizations import get_normalizations


def convert_units(
    x: np.ndarray,
    d: np.ndarray,
    filename: str | Path | None = None,
    b0: float | None = None,
    n0: float | None = None,
    l0: float | None = None,
    mi: float | None = None,
    cgs: bool = False,
    mks: bool = False,
) -> np.ndarray:
    """
    Python port of convert_units.pro.
    Returns converted array (does not mutate input).
    """
    arr = np.asarray(x, dtype=float)
    dim = np.asarray(d, dtype=float)

    if arr.size == 0 or not (cgs or mks):
        return arr.copy()

    if b0 is None or n0 is None or l0 is None or mi is None:
        if filename is None:
            b0, n0, l0, mi = 1.0e4, 1.0e14, 100.0, 1.0
        else:
            b0, n0, l0, mi = get_normalizations(filename)

    if b0 == 0 or n0 == 0 or l0 == 0 or mi == 0:
        b0, n0, l0, mi = 1.0e4, 1.0e14, 100.0, 1.0

    fp = 4.0 * np.pi
    c0 = 3.0e10
    v0 = 2.18e11 * b0 / np.sqrt(mi * n0)
    t0 = l0 / v0
    temp0 = (b0**2) / (fp * n0) / (1.6022e-12)
    i0 = c0 * b0 * l0 / fp
    e0 = b0**2 * l0**3 / fp
    pot0 = l0 * v0 * b0 / c0

    val_cgs = (
        (fp ** dim[0])
        * (c0 ** dim[1])
        * (n0 ** dim[2])
        * (v0 ** dim[3])
        * (b0 ** dim[4])
        * (t0 ** dim[8])
        * (l0 ** dim[9])
        * (temp0 ** dim[5])
        * (i0 ** dim[6])
        * (e0 ** dim[7])
        * (pot0 ** dim[10])
    )

    if cgs:
        return arr * val_cgs

    # mks: cgs conversion followed by cgs->mks factors (IDL behavior)
    val_mks = (
        ((1.0e-2) ** dim[1])
        * ((1.0e6) ** dim[2])
        * ((1.0e-2) ** dim[3])
        * ((1.0e-4) ** dim[4])
        * ((1.0e-2) ** dim[9])
        * ((3.0e9) ** (-dim[6]))
        * ((1.0e-7) ** dim[7])
        * ((3.0e2) ** dim[10])
    )
    return arr * val_cgs * val_mks
