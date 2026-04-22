from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .plot_legend import plot_legend
from .read_field_spectrum import read_field_spectrum
from .read_parameter import read_parameter


def _find_profile_crossings(values, xvals, target: float) -> np.ndarray:
    vals = np.asarray(values, dtype=float).reshape(-1)
    xx = np.asarray(xvals, dtype=float).reshape(-1)
    n = min(vals.size, xx.size)
    if n < 2:
        return np.asarray([], dtype=float)
    vals = vals[:n]
    xx = xx[:n]
    out: list[float] = []
    for i in range(n - 1):
        v0 = vals[i]
        v1 = vals[i + 1]
        x0 = xx[i]
        x1 = xx[i + 1]
        if not (np.isfinite(v0) and np.isfinite(v1) and np.isfinite(x0) and np.isfinite(x1)):
            continue
        d0 = v0 - target
        d1 = v1 - target
        if d0 == 0.0:
            out.append(float(x0))
        if d0 == 0.0 and d1 == 0.0:
            out.append(float(x1))
            continue
        if d0 * d1 < 0.0:
            frac = (target - v0) / (v1 - v0)
            out.append(float(x0 + frac * (x1 - x0)))
        elif d1 == 0.0:
            out.append(float(x1))
    if not out:
        return np.asarray([], dtype=float)
    arr = np.asarray(out, dtype=float)
    arr = np.unique(np.round(arr, decimals=12))
    return np.sort(arr)


def plot_field_spectrum(
    name: str,
    timeslices=0,
    *,
    m_val=None,
    phase: bool = False,
    overplot: bool = False,
    linestyle: str = "-",
    sqrtpsin: bool = True,
    symbol: str | None = None,
    units: str | None = None,
    filename: str | Path = "C1.h5",
    points: int = 200,
    xrange=None,
    yrange=None,
    linear: bool = False,
    complex: bool = True,
    logical: bool = False,
    phi: float = 0.0,
    operation: int = 1,
    fac: float | None = None,
    linfac: float = 1.0,
    equilibrium: bool = False,
    cgs: bool = False,
    mks: bool = False,
    psi_norm: bool = False,
    phi_norm: bool = False,
    rho: bool = False,
    psi0=None,
    i0=None,
    fc=None,
    ignore_jacobian: bool = False,
    dpsi0_dx=None,
    dpsi0_dz=None,
    pest: bool = False,
    boozer: bool = False,
    hamada: bool = False,
    fast: bool = False,
    fbins: int | None = None,
    tbins: int | None = None,
    psin_range=None,
    wall_mask: bool = False,
    edge_val=None,
    **kwargs,
):
    """
    Plot selected poloidal components from read_field_spectrum().
    """
    if "m_vals" in kwargs:
        raise TypeError("plot_field_spectrum() uses 'm_val', not 'm_vals'.")

    spec = read_field_spectrum(
        name,
        timeslices,
        filename=filename,
        points=points,
        xrange=xrange,
        yrange=yrange,
        linear=linear,
        complex=complex,
        logical=logical,
        phi=phi,
        operation=operation,
        fac=fac,
        linfac=linfac,
        equilibrium=equilibrium,
        cgs=cgs,
        mks=mks,
        m_val=m_val,
        psi_norm=psi_norm,
        phi_norm=phi_norm,
        rho=rho,
        psi0=psi0,
        i0=i0,
        fc=fc,
        ignore_jacobian=ignore_jacobian,
        dpsi0_dx=dpsi0_dx,
        dpsi0_dz=dpsi0_dz,
        pest=pest,
        boozer=boozer,
        hamada=hamada,
        fast=fast,
        fbins=fbins,
        tbins=tbins,
        psin_range=psin_range,
        wall_mask=wall_mask,
        edge_val=edge_val,
        **kwargs,
    )
    if symbol is None:
        symbol = str(spec.symbol)
    if units is None:
        units = str(spec.units)

    d = np.asarray(spec.data, dtype=np.complex128)
    m = np.asarray(spec.m, dtype=int)
    if m_val is None:
        amp = np.nanmax(np.abs(d), axis=(0, 2))
        amp = np.where(np.isfinite(amp), amp, -np.inf)
        take = min(5, int(m.size))
        idx = np.argsort(amp)[-take:][::-1]
        mvals = m[idx]
        d = d[:, idx, :]
        m = m[idx]
    else:
        mvals = np.asarray(m_val, dtype=int).reshape(-1)
    nflux = np.asarray(spec.fc.psi_norm, dtype=float)
    qprof = np.asarray(spec.fc.q, dtype=float)
    ntor = int(read_parameter("ntor", filename=filename, cgs=cgs, mks=mks))
    q_target = np.abs(mvals / float(ntor))
    ax = plt.gca() if overplot else plt.subplots(figsize=(7, 5))[1]
    fig = ax.figure
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
    colors = [cycle[i % len(cycle)] for i in range(max(len(mvals), 1))]
    label = f"{symbol} ({units})" if symbol is not None and units is not None else str(name)
    ytitle = "Phase" if phase else label
    plotted = []
    xplot = np.asarray(spec.x, dtype=float)
    if sqrtpsin and str(spec.xname) == "psi_norm":
        xplot = np.sqrt(np.maximum(xplot, 0.0))
    xlabel = "sqrt(psi_norm)" if sqrtpsin and str(spec.xname) == "psi_norm" else str(spec.xname)

    for i, mv in enumerate(mvals):
        j = int(np.argmin(np.abs(m - mv)))
        data = np.angle(d[0, j, :]) if phase else np.abs(d[0, j, :])
        plotted.append(np.asarray(data, dtype=float))
        if i == 0 and not overplot:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ytitle)
        ax.plot(xplot, data, color=colors[i], linestyle=linestyle)
        for fv in _find_profile_crossings(np.abs(qprof), xplot, float(q_target[i])):
            ax.axvline(fv, color=colors[i], linestyle="--", linewidth=0.8)

    ax.set_xlim(0.0, 1.0)
    finite = np.concatenate([p[np.isfinite(p)] for p in plotted]) if plotted else np.asarray([], dtype=float)
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
        ax.set_ylim(ymin_pad, ymax_pad)
        if (ymin_pad < 0.0) and (ymax_pad > 0.0):
            ax.axhline(0.0, color="black", linestyle="-", linewidth=0.8)

    if mvals.size > 1:
        names = [f"m = {int(v)}" for v in mvals]
        plot_legend(names, colors=colors[: mvals.size], linestyles=[linestyle] * mvals.size)

    return fig, ax
