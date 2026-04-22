from __future__ import annotations

from pathlib import Path

import numpy as np

from .plot_signals import plot_signals


def plot_mag_probes(
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
    # IDL structure parity: plot_mag_probes -> plot_signals('mag_probes')
    return plot_signals(
        "mag_probes",
        filename=filename,
        deriv=deriv,
        power_spectrum=power_spectrum,
        compensate_renorm=compensate_renorm,
        scale=scale,
        overplot=overplot,
        noplot=noplot,
        outfile=outfile,
        cgs=cgs,
        mks=mks,
    )
