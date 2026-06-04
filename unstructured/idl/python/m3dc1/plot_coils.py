from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .read_coil_data import read_coil_data


def plot_coils(*, filename=None, directory=None, overplot=False, rmp=False, **kwargs):
    if directory is None:
        directory = "."
        if filename is not None:
            directory = str(Path(str(filename)).parent)
    dat = read_coil_data(directory=directory, rmp=rmp)
    if dat is None:
        return

    n = dat.shape[1]
    if not overplot:
        plt.figure(figsize=(6, 6))

    for i in range(n):
        a2 = np.tan(np.deg2rad(dat[7, i]) + (np.pi / 2.0 if dat[7, i] != 0.0 else 0.0))
        a1 = np.tan(np.deg2rad(dat[6, i]))
        xp = np.asarray(
            [
                dat[2, i] - (dat[4, i] - dat[5, i] * a2) / 2.0,
                dat[2, i] + (dat[4, i] + dat[5, i] * a2) / 2.0,
                dat[2, i] + (dat[4, i] - dat[5, i] * a2) / 2.0,
                dat[2, i] - (dat[4, i] + dat[5, i] * a2) / 2.0,
            ],
            dtype=float,
        )
        zp = np.asarray(
            [
                dat[3, i] - (dat[5, i] + dat[4, i] * a1) / 2.0,
                dat[3, i] - (dat[5, i] - dat[4, i] * a1) / 2.0,
                dat[3, i] + (dat[5, i] + dat[4, i] * a1) / 2.0,
                dat[3, i] + (dat[5, i] - dat[4, i] * a1) / 2.0,
            ],
            dtype=float,
        )
        plt.plot(np.r_[xp, xp[0]], np.r_[zp, zp[0]], color="tab:orange", linewidth=1.0)
