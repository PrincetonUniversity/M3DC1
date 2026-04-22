from __future__ import annotations

import numpy as np


def compensate_renorm(x: np.ndarray) -> np.ndarray:
    """
    Python port of compensate_renorm.pro.
    """
    arr = np.asarray(x, dtype=float)
    y = arr.copy()
    n = arr.size
    for k in range(n - 1):
        if arr[k] * arr[k + 1] <= 0.0:
            continue
        ratio = abs(arr[k + 1] / arr[k]) if arr[k] != 0 else np.inf
        if ratio < 1e-9:
            dx = arr[k + 2] - arr[k + 1] if k < n - 2 else 0.0
            f = arr[k + 1] * np.exp(-dx / (arr[k + 1] + dx / 2.0))
            if arr[k] != 0:
                y[: k + 1] = y[: k + 1] * f / arr[k]
    return y
