from __future__ import annotations

import numpy as np


def dimensions(
    energy: int = 0,
    eta: int = 0,
    j0: int = 0,
    p0: int = 0,
    temperature: int = 0,
    potential: int = 0,
    l0: int = 0,
    light: int = 0,
    fourpi: int = 0,
    n0: int = 0,
    v0: int = 0,
    t0: int = 0,
    b0: int = 0,
    mu0: int = 0,
) -> np.ndarray:
    """
    Python port of dimensions.pro.
    """
    d = np.zeros(11, dtype=int)

    fp = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=int)
    c = np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=int)
    den = np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype=int)
    vel = np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], dtype=int)
    mag = np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype=int)
    temp = np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype=int)
    i0 = np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype=int)
    e0 = np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype=int)
    time = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], dtype=int)
    length = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype=int)
    pot = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype=int)
    p = e0 - 3 * length

    d = d + fp * fourpi
    d = d + c * light
    d = d + den * n0
    d = d + vel * v0
    d = d + mag * b0
    d = d + time * t0
    d = d + length * l0
    d = d + temp * temperature
    d = d + i0 * j0 - 2 * length * j0
    d = d + e0 * energy
    d = d + pot * potential
    d = d + eta * (2 * length - time + fp - 2 * c)
    d = d + p0 * p
    d = d + mu0 * (p + time)

    return d
