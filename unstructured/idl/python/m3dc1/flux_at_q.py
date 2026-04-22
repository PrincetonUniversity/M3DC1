from __future__ import annotations

import numpy as np

from .flux_coordinates import flux_coordinates


def _curve_from_fields(filename="C1.h5", points=200, normalized_flux=False, **kwargs):
    fc = flux_coordinates(filename=filename, points=points, **kwargs)
    q = np.abs(np.asarray(fc.q, dtype=float).reshape(-1))
    f = np.asarray(fc.psi_norm if normalized_flux else fc.psi, dtype=float).reshape(-1)
    return q, f


def flux_at_q(
    qval,
    q=None,
    *,
    normalized_flux=False,
    points=200,
    flux=None,
    unique=False,
    filename="C1.h5",
    **kwargs,
):
    """
    Return flux values where q crosses qval.
    """
    qtargets = np.asarray(qval, dtype=float).reshape(-1)
    if qtargets.size == 0:
        return np.asarray([0.0], dtype=float)

    if q is None or flux is None:
        qcurve, fcurve = _curve_from_fields(filename=filename, points=points, normalized_flux=normalized_flux, **kwargs)
    else:
        qcurve = np.asarray(q, dtype=float).reshape(-1)
        fcurve = np.asarray(flux, dtype=float).reshape(-1)

    if qcurve.size < 2 or fcurve.size < 2:
        return np.asarray([0.0], dtype=float)

    out: list[float] = []
    for qt in qtargets:
        hits: list[float] = []
        diff = qcurve - float(qt)
        for i in range(1, diff.size):
            d0 = diff[i - 1]
            d1 = diff[i]
            if d0 == 0.0:
                hits.append(float(fcurve[i - 1]))
            elif d0 * d1 <= 0.0:
                denom = (d1 - d0)
                frac = 0.0 if abs(denom) < np.finfo(float).tiny else -d0 / denom
                hits.append(float(fcurve[i - 1] + frac * (fcurve[i] - fcurve[i - 1])))
                if unique:
                    break
        out.extend(hits)

    if len(out) == 0:
        return np.asarray([0.0], dtype=float)
    return np.asarray(out, dtype=float)
