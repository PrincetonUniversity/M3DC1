from __future__ import annotations

import numpy as np

from .get_slice_time import get_slice_time
from .read_scalar import read_scalar


def read_nulls(*, filename="C1.h5", slice=0, cgs=False, mks=False, **kwargs):
    t0 = float(np.asarray(get_slice_time(filename=filename, slice=int(slice), cgs=cgs, mks=mks)).reshape(-1)[0])
    t = np.asarray(read_scalar("time", filename=filename, cgs=cgs, mks=mks), dtype=float).reshape(-1)
    idx = int(np.argmin(np.abs(t - t0))) if t.size > 0 else 0

    def _at(name: str, default: float = 0.0) -> float:
        try:
            arr = np.asarray(read_scalar(name, filename=filename, cgs=cgs, mks=mks), dtype=float).reshape(-1)
            return float(arr[idx]) if arr.size > 0 else default
        except Exception:
            return default

    axis = np.asarray([_at("xmag"), _at("zmag")], dtype=float)
    xpoint = np.asarray([_at("xnull"), _at("znull")], dtype=float)
    return axis, xpoint

