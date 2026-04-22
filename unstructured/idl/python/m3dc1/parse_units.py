from __future__ import annotations

import numpy as np


def _fmt_unit_pow(name: str, exp: int) -> str:
    if exp == 1:
        return name
    return f"{name}^{{{exp}}}"


def parse_units(x: np.ndarray, cgs: bool = False, mks: bool = False) -> str:
    """
    Python port of parse_units.pro.
    Input x is an 11-element dimension vector.
    """
    v = np.asarray(x, dtype=int).copy()

    if cgs or mks:
        v[9] = v[9] - 3 * v[2] + v[3] + v[1]
        v[8] = v[8] - v[3] - v[1]
        v[0:4] = 0
        base = [
            "4\\pi",
            "c",
            "cm",
            "v_{A0}",
            "G" if cgs else "T",
            "eV",
            "statA" if cgs else "A",
            "erg" if cgs else "J",
            "s",
            "cm" if cgs else "m",
            "statV" if cgs else "V",
        ]
    else:
        v[0] = v[0] - v[5] - v[6] - v[7]
        v[1] = v[1] + v[6]
        v[4] = v[4] + 2 * v[5] + v[6] + 2 * v[7]
        v[9] = v[9] + v[6] + 3 * v[7]
        v[5:8] = 0
        base = [
            "4\\pi",
            "c",
            "n_{0}",
            "v_{A0}",
            "B_{0}",
            "temp",
            "curr",
            "energy",
            "\\tau_{A0}",
            "L_{0}",
            "potential",
        ]

    pos = [_fmt_unit_pow(base[i], int(v[i])) for i in range(len(v)) if v[i] > 0]
    neg = [_fmt_unit_pow(base[i], int(-v[i])) for i in range(len(v)) if v[i] < 0]

    if pos and neg:
        units = f"{' '.join(pos)} / {' '.join(neg)}"
    elif pos:
        units = " ".join(pos)
    elif neg:
        units = " ".join(_fmt_unit_pow(base[i], int(v[i])) for i in range(len(v)) if v[i] < 0)
    else:
        units = ""

    if units and any(ch in units for ch in ("_", "^", "\\")):
        return f"${units}$"
    return units
