from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .read_signals import read_signals


def _write_ascii(outfile: str | Path, tdata: np.ndarray, data: np.ndarray) -> None:
    arr = np.column_stack([tdata, data.T])
    np.savetxt(outfile, arr, fmt="%12.4g")


def plot_signals(
    signal: str,
    filename: str | Path = "C1.h5",
    deriv: bool = False,
    power_spectrum: bool = False,
    compensate_renorm: bool = False,
    scale: bool = False,
    overplot: bool = False,
    noplot: bool = False,
    outfile: str | Path | None = None,
    cgs: bool = False,
    mks: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
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
    tdata = meta.tdata
    data = meta.data
    xtitle = meta.xtitle
    ylabel = meta.ylabel

    if outfile is not None:
        _write_ascii(outfile, tdata, data)

    if not noplot:
        if not overplot:
            plt.figure(figsize=(10, 5))
        for i in range(data.shape[0]):
            plt.plot(tdata, data[i, :], label=f"{signal} {i}")
        plt.xlabel(xtitle)
        plt.ylabel(ylabel)
        finite = np.asarray(data, dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size > 0:
            ymin = float(np.min(finite))
            ymax = float(np.max(finite))
            if ymax > ymin:
                yrange_pad = 0.1 * (ymax - ymin)
                ymin_pad = ymin - yrange_pad
                ymax_pad = ymax + yrange_pad
            else:
                yrange_pad = max(abs(ymin), 1.0) * 0.1
                ymin_pad = ymin - yrange_pad
                ymax_pad = ymax + yrange_pad
            plt.ylim(bottom=ymin_pad, top=ymax_pad)
            if (ymin_pad < 0.0) and (ymax_pad > 0.0):
                plt.axhline(0.0, color="black", linestyle="-", linewidth=0.8)
        if data.shape[0] <= 12:
            plt.legend(loc="best", fontsize=8)
        plt.tight_layout()

    return tdata, data
