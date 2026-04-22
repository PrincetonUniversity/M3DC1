from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .read_signals import read_signals


@dataclass
class SignalZeroPointResult:
    zero_points: list[np.ndarray]
    frequency: list[np.ndarray]

def read_signals_zeropoint(
    signal: str,
    filename: str | Path = "C1.h5",
    deriv: bool = False,
    power_spectrum: bool = False,
    compensate_renorm: bool = False,
    scale: bool = False,
    cgs: bool = False,
    mks: bool = False,
    min_sep: float = 3e-7,
    max_points: int = 100,
    return_meta: bool = False,
):
    meta = read_signals(
        signal,
        filename=filename,
        deriv=deriv,
        power_spectrum=power_spectrum,
        compensate_renorm=compensate_renorm,
        scale=scale,
        cgs=cgs,
        mks=mks,
        return_meta=True,
    )

    tdata = np.asarray(meta.tdata, dtype=float)
    data = np.asarray(meta.data, dtype=float)

    if power_spectrum:
        zero_points = [np.array([], dtype=float) for _ in range(data.shape[0])]
        frequency = [np.array([], dtype=float) for _ in range(data.shape[0])]
    else:
        zero_points = []
        frequency = []
        for i in range(data.shape[0]):
            tt = np.asarray(tdata, dtype=float).reshape(-1)
            yy = np.asarray(data[i, :], dtype=float).reshape(-1)
            n = min(tt.size, yy.size)
            if n < 2:
                zero_points.append(np.array([], dtype=float))
                frequency.append(np.array([], dtype=float))
                continue

            zeros: list[float] = []
            for j in range(1, n):
                y0 = yy[j]
                y1 = yy[j - 1]
                if y0 * y1 < 0.0:
                    zt = (tt[j] * abs(y1) + tt[j - 1] * abs(y0)) / (abs(y1) + abs(y0))
                    if len(zeros) == 0 or abs(zeros[-1] - zt) > float(min_sep):
                        zeros.append(float(zt))
                        if len(zeros) >= int(max_points):
                            break

            zp = np.asarray(zeros, dtype=float)
            if zp.size < 2:
                fr = np.array([], dtype=float)
            else:
                dt = np.abs(np.diff(zp))
                fr = np.where(dt > 0.0, 0.5 / dt, np.nan)
            zero_points.append(zp)
            frequency.append(np.asarray(fr, dtype=float))

    if return_meta:
        return SignalZeroPointResult(zero_points=zero_points, frequency=frequency)

    return zero_points, frequency
