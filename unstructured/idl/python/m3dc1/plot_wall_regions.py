from __future__ import annotations

import h5py
import matplotlib.pyplot as plt
import numpy as np

from .time_name import time_name


def plot_wall_regions(*, filename="C1.h5", slice=0, over=False, color="tab:green", **kwargs):
    """
    Plot wall region polylines if present in HDF5.
    """
    with h5py.File(str(filename), "r") as h5:
        gname = time_name(int(slice))
        if gname not in h5:
            return
        tg = h5[gname]
        version = float(tg.attrs.get("version", 0.0))
        if version < 40:
            return
        nregions = int(tg.attrs.get("iwall_regions", 0))
        if nregions <= 0:
            return

        if not over:
            plt.figure(figsize=(6, 6))
        for i in range(1, nregions + 1):
            rname = f"wall_region_{i:03d}"
            if rname not in tg or not isinstance(tg[rname], h5py.Group):
                continue
            rg = tg[rname]
            nplanes = int(rg.attrs.get("nplanes", 0))
            for j in range(1, nplanes + 1):
                pname = f"plane_{j:03d}"
                if pname not in rg or not isinstance(rg[pname], h5py.Group):
                    continue
                pg = rg[pname]
                if "x" not in pg or "y" not in pg:
                    continue
                x = np.asarray(pg["x"][()], dtype=float).reshape(-1)
                y = np.asarray(pg["y"][()], dtype=float).reshape(-1)
                plt.plot(x, y, color=color, linewidth=1.0)

