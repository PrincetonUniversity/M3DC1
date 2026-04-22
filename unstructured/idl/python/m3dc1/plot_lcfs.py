from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .get_lcfs import get_lcfs


def plot_lcfs(psi=None, x=None, z=None, *, psival=None, over=False, color="k", filename="C1.h5", slice=0, points=200, **kwargs):
    """
    Plot the last-closed flux surface.
    """
    xy = get_lcfs(psi, x, z, psival=psival, filename=filename, slice=slice, points=points, **kwargs)
    if xy.size == 0:
        return xy
    if over:
        plt.plot(xy[0, :], xy[1, :], color=color, linewidth=1.5)
    else:
        plt.figure(figsize=(6, 6))
        plt.plot(xy[0, :], xy[1, :], color=color, linewidth=1.5)
        plt.gca().set_aspect("equal", adjustable="box")
    return np.asarray(xy)
