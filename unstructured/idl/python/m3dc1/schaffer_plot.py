from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .contour_and_legend import contour_and_legend
from .field_spectrum import field_spectrum
from .flux_coordinates import flux_coordinates
from .plot_legend import plot_legend
from .read_field import read_field
from .read_field_spectrum import read_field_spectrum
from .read_parameter import read_parameter


def _as_1d(x) -> np.ndarray:
    return np.asarray([] if x is None else x, dtype=float).reshape(-1)


def _interp_indices(values, targets) -> np.ndarray:
    vals = np.asarray(values, dtype=float).reshape(-1)
    tgt = np.asarray(targets, dtype=float).reshape(-1)
    if vals.size == 0 or tgt.size == 0:
        return np.asarray([], dtype=float)
    order = np.argsort(vals)
    return np.interp(tgt, vals[order], np.arange(vals.size, dtype=float)[order])


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


def schaffer_plot(
    field,
    timeslices=-1,
    x=None,
    z=None,
    *,
    q=None,
    bins: int = 200,
    q_val=None,
    psi_val=None,
    ntor=None,
    tpoints: int | None = None,
    psi0=None,
    i0=None,
    m_val=None,
    phase: bool = False,
    overplot: bool = False,
    linestyle: str = "-",
    levels=None,
    lines: bool = False,
    outfile: str | Path | None = None,
    bmnfile: str | Path | None = None,
    bmncdf: str | Path | None = None,
    rhs: bool = False,
    reverse_q: bool = False,
    sqrtpsin: bool = True,
    profdata: bool = False,
    boozer: bool = False,
    pest: bool = False,
    hamada: bool = False,
    geo: bool = False,
    symbol: str | None = None,
    units: str | None = None,
    dpsi0_dx=None,
    dpsi0_dz=None,
    filename: str | Path = "C1.h5",
    points: int = 200,
    logical: bool = False,
    phi: float = 0.0,
    xrange=None,
    yrange=None,
    **kwargs,
):
    """
    Python port of schaffer_plot.pro.
    """
    del q, rhs, reverse_q
    print("Drawing schaffer plot")
    aux_kwargs = {k: v for k, v in kwargs.items() if k in {"cgs", "mks"}}

    if not boozer and not hamada and not geo:
        pest = True
    threed = int(read_parameter("3d", filename=filename, **aux_kwargs))
    if threed == 1 and tpoints is None:
        tpoints = 32

    spec = None
    if isinstance(field, str):
        spec = read_field_spectrum(
            field,
            timeslices,
            filename=filename,
            points=points,
            tpoints=tpoints,
            ntor=ntor,
            xrange=xrange,
            yrange=yrange,
            logical=logical,
            phi=phi,
            psi0=psi0,
            i0=i0,
            fc=None,
            dpsi0_dx=dpsi0_dx,
            dpsi0_dz=dpsi0_dz,
            pest=pest,
            boozer=boozer,
            hamada=hamada,
            fast=geo,
            fbins=bins,
            tbins=bins,
            **kwargs,
        )
        if symbol is None:
            symbol = str(spec.symbol)
        if units is None:
            units = str(spec.units)
        fc = spec.fc
        d = np.asarray(spec.data, dtype=np.complex128)
        m = np.asarray(spec.m, dtype=int)
        n = np.asarray(spec.n, dtype=int)
        nflux = np.asarray(fc.psi_norm, dtype=float)
        qprof = np.asarray(fc.q, dtype=float)
        plot_idx = 0
        ntor_plot = int(n[0]) if n.size > 0 else int(read_parameter("ntor", filename=filename, **aux_kwargs))
    else:
        field_data = np.asarray(field)
        if x is None or z is None:
            raise ValueError("schaffer_plot requires x and z when field data is passed directly.")
        xv = np.asarray(x, dtype=float).reshape(-1)
        zv = np.asarray(z, dtype=float).reshape(-1)

        if psi0 is None:
            print("READING PSI IN SCHAFFER_PLOT")
            psi_meta = read_field(
                "psi",
                timeslices=-1,
                points=points,
                xrange=xrange,
                yrange=yrange,
                logical=logical,
                phi=phi,
                filename=filename,
                return_meta=True,
                **kwargs,
            )
            psi0 = np.asarray(psi_meta.data)
            psi0 = psi0[0, :, :] if psi0.ndim == 3 else psi0
        else:
            psi0 = np.asarray(psi0)
            if psi0.ndim == 3:
                psi0 = psi0[0, :, :]

        if i0 is None:
            i_meta = read_field(
                "i",
                timeslices=-1,
                points=points,
                xrange=xrange,
                yrange=yrange,
                logical=logical,
                phi=phi,
                filename=filename,
                return_meta=True,
                **kwargs,
            )
            i0 = np.asarray(i_meta.data)
            i0 = i0[0, :, :] if i0.ndim == 3 else i0
        else:
            i0 = np.asarray(i0)
            if i0.ndim == 3:
                i0 = i0[0, :, :]

        if ntor is None:
            ntor = int(read_parameter("ntor", filename=filename, **aux_kwargs))
        coord_linear = bool(read_parameter("linear", filename=filename, **aux_kwargs))

        fc = flux_coordinates(
            psi0=psi0,
            i0=i0,
            x=xv,
            z=zv,
            tbins=bins,
            fbins=bins,
            pest=pest,
            boozer=boozer,
            hamada=hamada,
            dpsi0_dx=dpsi0_dx,
            dpsi0_dz=dpsi0_dz,
            filename=filename,
            slice=-1 if coord_linear else (-1 if timeslices == -1 else int(timeslices)),
            **aux_kwargs,
        )

        spec = field_spectrum(
            field_data,
            xv,
            zv,
            fc=fc,
            psi0=psi0,
            i0=i0,
            tbins=bins,
            fbins=bins,
            pest=pest,
            boozer=boozer,
            hamada=hamada,
            dpsi0_dx=dpsi0_dx,
            dpsi0_dz=dpsi0_dz,
            filename=filename,
            **aux_kwargs,
        )
        d = np.asarray(spec.data, dtype=np.complex128)
        fc = spec.fc
        m = np.asarray(spec.m, dtype=int)
        n = np.asarray(spec.n, dtype=int)
        nflux = np.asarray(fc.psi_norm, dtype=float)
        qprof = np.asarray(fc.q, dtype=float)
        plot_idx = 0
        ntor_plot = int(ntor)
        if n.size > 1:
            matches = np.where(n == int(ntor))[0]
            if matches.size != 1:
                print(f"schaffer_plot: error: ntor = {ntor} not found.")
                return None, None
            plot_idx = int(matches[0])
            print(f"plotting ntor = {ntor}")

    label = " "
    if symbol is not None and units is not None:
        label = f"{symbol} ({units})"

    if bmnfile is not None:
        with open(bmnfile, "w", encoding="ascii") as fh:
            fh.write(f"{int(ntor_plot):5d}{nflux.size:5d}{m.size:5d}\n")
            for val in nflux:
                fh.write(f"{float(val):13.6f}\n")
            for val in m:
                fh.write(f"{int(val):5d}\n")
            for i in range(nflux.size):
                for j in range(m.size):
                    fh.write(f"{float(np.real(d[plot_idx, j, i])):13.6f}{float(np.imag(d[plot_idx, j, i])):13.6f}\n")

    if bmncdf is not None or profdata:
        print("schaffer_plot: bmncdf/profdata output is not implemented in Python yet.")

    if m_val is not None:
        mvals = np.asarray(m_val, dtype=int).reshape(-1)
        if ntor_plot == 0:
            q_target = np.full(mvals.shape, np.nan, dtype=float)
        else:
            q_target = np.abs(mvals / float(ntor_plot))
        ax = plt.gca() if overplot else plt.subplots(figsize=(7, 5))[1]
        fig = ax.figure
        cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
        if not cycle:
            cycle = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
        colors = [cycle[i % len(cycle)] for i in range(max(len(mvals), 1))]
        ytitle = "Phase" if phase else label
        plotted = []
        xplot = np.sqrt(np.maximum(nflux, 0.0)) if sqrtpsin else nflux
        for i, mv in enumerate(mvals):
            j = int(np.argmin(np.abs(m - mv)))
            data = np.angle(d[plot_idx, j, :]) if phase else np.abs(d[plot_idx, j, :])
            plotted.append(np.asarray(data, dtype=float))
            if i == 0 and not overplot:
                ax.set_xlabel("sqrt(psi_norm)" if sqrtpsin else "psi_norm")
                ax.set_ylabel(ytitle)
            ax.plot(xplot, data, color=colors[i], linestyle=linestyle)
            if np.isfinite(q_target[i]):
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

    indices = np.asarray([], dtype=float)
    q_input = _as_1d(q_val)
    psi_input = _as_1d(psi_val)
    if psi_input.size > 0:
        indices = _interp_indices(nflux, psi_input)
    elif q_input.size > 0:
        indices = _interp_indices(np.abs(qprof), q_input)

    if indices.size > 0:
        print(np.abs(qprof[np.clip(indices.astype(int), 0, qprof.size - 1)]))
        print(np.abs(qprof[np.clip((indices + 1).astype(int), 0, qprof.size - 1)]))
        outfh = open(outfile, "w", encoding="ascii") if outfile is not None else None
        try:
            for i, idx in enumerate(indices):
                idxf = float(idx)
                j = int(np.argmin(np.abs(m - q_input[i] * ntor_plot))) if q_input.size > i else int(np.argmin(np.abs(m)))
                k = int(np.argmin(np.abs(m + q_input[i] * ntor_plot))) if q_input.size > i else j
                q_here = float(np.interp(idxf, np.arange(qprof.size, dtype=float), np.abs(qprof)))
                psi_here = float(np.interp(idxf, np.arange(nflux.size, dtype=float), nflux))
                print("q, Psi = ", q_here, psi_here)
                dj = np.interp(idxf, np.arange(nflux.size, dtype=float), d[plot_idx, j, :])
                dk = np.interp(idxf, np.arange(nflux.size, dtype=float), d[plot_idx, k, :])
                print("Resonant field: m (mag, phase) = ", int(m[j]), abs(dj), float(np.angle(dj)))
                print("Resonant field: m (mag, phase) = ", int(m[k]), abs(dk), float(np.angle(dk)))
                if outfh is not None:
                    dd = dj if qprof[0] > 0 else dk
                    mi = j if qprof[0] > 0 else k
                    dq = np.gradient(np.abs(qprof), nflux, edge_order=1)
                    dpsi = np.gradient(np.asarray(fc.psi, dtype=float), nflux, edge_order=1)
                    outfh.write(
                        f"{int(m[mi]):5d}{abs(dd):12.6f}{float(np.angle(dd)):12.6f}"
                        f"{psi_here:12.6f}{q_here:12.6f}"
                        f"{float(np.interp(idxf, np.arange(dq.size, dtype=float), dq)):12.6f}"
                        f"{float(np.interp(idxf, np.arange(fc.area.size, dtype=float), fc.area)):12.6f}"
                        f"{float(np.interp(idxf, np.arange(dpsi.size, dtype=float), dpsi)):12.6f}\n"
                    )
        finally:
            if outfh is not None:
                outfh.close()

    y = np.sqrt(np.maximum(nflux, 0.0)) if sqrtpsin else nflux
    ytitle = "sqrt(psi_norm)" if sqrtpsin else "psi_norm"
    contour_and_legend(
        np.abs(d[plot_idx, :, :]),
        m,
        y,
        levels=100 if levels is None else levels,
        lines=lines,
        label=label,
        xtitle="m",
        ytitle=ytitle,
        overplot=overplot,
    )
    ax = plt.gca()
    ax.set_xlim(-20, 20)
    ax.set_ylim(0, 1)
    ax.plot(ntor_plot * qprof, y, linestyle="--", color="black")
    return ax.figure, ax
