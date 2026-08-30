from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .get_slice_time import get_slice_time
from .read_parameter import read_parameter
from .read_scalar import read_scalar


@dataclass
class LcfsResult:
    psilim: float
    flux0: float
    axis: np.ndarray
    xpoint: np.ndarray


def read_lcfs(
    filename: str | Path = "C1.h5",
    slice: int | None = None,
    last: bool = False,
    cgs: bool = False,
    mks: bool = False,
    return_meta: bool = False,
):
    """Lightweight Python port of read_lcfs.pro from scalar tracks."""
    tvec = np.asarray(read_scalar("time", filename=filename, cgs=cgs, mks=mks), dtype=float).reshape(-1)
    nscalars = tvec.size
    if nscalars == 0:
        raise ValueError("No scalar time data found.")
    ntime = int(read_parameter("ntime", filename=filename))
    if ntime <= 0:
        ntime = nscalars

    if last:
        idx = ntime - 1
    elif slice is None:
        idx = 0
    else:
        idx = int(slice)
        if idx < 0:
            idx = 0
    idx = max(0, min(idx, ntime - 1))

    t0 = np.asarray(get_slice_time(filename=filename, slice=idx, cgs=cgs, mks=mks), dtype=float).reshape(-1)
    if t0.size > 0:
        dt = np.abs(tvec - t0[0])
        dt[~np.isfinite(dt)] = np.inf
        it = int(np.argmin(dt))
        print("slice time = ", t0)
        print("time step time: ", float(tvec[it]))
        print("time slice: ", it)
        idx = it

    # Fallback-safe scalar arrays
    def _at(name: str, default: float = 0.0) -> float:
        try:
            arr = np.asarray(read_scalar(name, filename=filename, cgs=cgs, mks=mks), dtype=float).reshape(-1)
            if arr.size == 0:
                return default
            return float(arr[idx])
        except Exception:
            return default

    xmag = _at("xmag", 0.0)
    zmag = _at("zmag", 0.0)
    xnull = _at("xnull", xmag)
    znull = _at("znull", zmag)
    flux0 = _at("psimin", 0.0)
    psilim = _at("psibound", flux0)

    if return_meta:
        return LcfsResult(
            psilim=psilim,
            flux0=flux0,
            axis=np.asarray([xmag, zmag], dtype=float),
            xpoint=np.asarray([xnull, znull], dtype=float),
        )
    return psilim
