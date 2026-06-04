from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .read_field import read_field


def plot_flux_contour(fval, *, filename="C1.h5", points=200, slice=0, overplot=False, closed=0, color="k", **kwargs):
    """
    Plot psi contour(s) at requested flux values.
    """
    vals = np.asarray(fval, dtype=float).reshape(-1)
    if vals.size == 0:
        return
    p = read_field("psi", filename=filename, timeslices=slice, points=points, equilibrium=True, return_meta=True, **kwargs)
    psi = np.asarray(p.data)
    psi2d = psi[0, :, :] if psi.ndim == 3 else psi
    x = np.asarray(p.r, dtype=float).reshape(-1)
    z = np.asarray(p.z, dtype=float).reshape(-1)
    if not overplot:
        plt.figure(figsize=(6, 5))
    plt.contour(x, z, psi2d.T, levels=np.sort(vals), colors=color, linewidths=0.8)

