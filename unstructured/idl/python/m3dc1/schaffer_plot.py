from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .contour_and_legend import contour_and_legend
from .field_spectrum import field_spectrum
from .get_colors import get_colors
from .plot_legend import plot_legend
from .read_field import read_field
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


def schaffer_plot(
    field,
    x=None,
    z=None,
    timeslices=-1,
    *,
    q=None,
    bins: int = 200,
    q_val=None,
    psi_val=None,
    ntor=None,
    psi0=None,
    i0=None,
    m_val=None,
    phase: bool = False,
    overplot: bool = False,
    linestyle: str = "-",
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

    if isinstance(field, str):
        try:
            meta = read_field(
                field,
                timeslices=timeslices,
                points=points,
                xrange=xrange,
                yrange=yrange,
                complex=True,
                logical=logical,
                phi=phi,
                filename=filename,
                return_meta=True,
                **kwargs,
            )
        except KeyError:
            print(f"schaffer_plot: complex read for '{field}' is unavailable; falling back to real field.")
            meta = read_field(
                field,
                timeslices=timeslices,
                points=points,
                xrange=xrange,
                yrange=yrange,
                complex=False,
                logical=logical,
                phi=phi,
                filename=filename,
                return_meta=True,
                **kwargs,
            )
        field_data = np.asarray(meta.data)
        if x is None:
            x = np.asarray(meta.r, dtype=float).reshape(-1)
        if z is None:
            z = np.asarray(meta.z, dtype=float).reshape(-1)
        if symbol is None:
            symbol = str(meta.symbol)
        if units is None:
            units = str(meta.units)
    else:
        field_data = np.asarray(field)
        if x is None or z is None:
            raise ValueError("schaffer_plot requires x and z when field data is passed directly.")

    xv = np.asarray(x, dtype=float).reshape(-1)
    zv = np.asarray(z, dtype=float).reshape(-1)

    label = " "
    if symbol is not None and units is not None:
        label = f"{symbol} ({units})"

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
        if x is None:
            xv = np.asarray(psi_meta.r, dtype=float).reshape(-1)
        if z is None:
            zv = np.asarray(psi_meta.z, dtype=float).reshape(-1)
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

    spec = field_spectrum(
        field_data,
        xv,
        zv,
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
        slice=-1 if timeslices == -1 else int(timeslices),
        **aux_kwargs,
    )
    d = np.asarray(spec.data, dtype=np.complex128)
    fc = spec.fc
    m = np.asarray(spec.m, dtype=int)
    n = np.asarray(spec.n, dtype=int)
    nflux = np.asarray(fc.psi_norm, dtype=float)
    qprof = np.asarray(fc.q, dtype=float)
    d = d / np.asarray(fc.area, dtype=float)[None, None, :]
    d = d / float(fc.dpsi_dchi)

    if bmnfile is not None:
        with open(bmnfile, "w", encoding="ascii") as fh:
            fh.write(f"{int(ntor):5d}{nflux.size:5d}{m.size:5d}\n")
            for val in nflux:
                fh.write(f"{float(val):13.6f}\n")
            for val in m:
                fh.write(f"{int(val):5d}\n")
            for i in range(nflux.size):
                for j in range(m.size):
                    fh.write(f"{float(np.real(d[0, j, i])):13.6f}{float(np.imag(d[0, j, i])):13.6f}\n")

    if bmncdf is not None or profdata:
        print("schaffer_plot: bmncdf/profdata output is not implemented in Python yet.")

    if m_val is not None:
        mvals = np.asarray(m_val, dtype=int).reshape(-1)
        q_target = np.abs(mvals / float(ntor))
        indices = _interp_indices(np.abs(qprof), q_target)
        colors = get_colors(max(len(mvals) + 2, 10))
        ax = plt.gca() if overplot else plt.subplots(figsize=(7, 5))[1]
        fig = ax.figure
        ytitle = "Phase" if phase else label
        plotted = []
        for i, mv in enumerate(mvals):
            j = int(np.argmin(np.abs(m - mv)))
            data = np.angle(d[0, j, :]) if phase else np.abs(d[0, j, :])
            plotted.append(np.asarray(data, dtype=float))
            if i == 0 and not overplot:
                ax.set_xlabel("sqrt(psi_norm)" if sqrtpsin else "psi_norm")
                ax.set_ylabel(ytitle)
            yplot = np.sqrt(np.maximum(nflux, 0.0)) if sqrtpsin else nflux
            ax.plot(yplot, data, color=colors[i + 1], linestyle=linestyle)
            if i < indices.size:
                fv = float(np.interp(indices[i], np.arange(nflux.size, dtype=float), yplot))
                ax.axvline(fv, color=colors[i + 1], linestyle="--", linewidth=0.8)
        ax.set_xlim(0.0, 1.0)
        finite = np.concatenate([p[np.isfinite(p)] for p in plotted]) if plotted else np.asarray([], dtype=float)
        if finite.size > 0:
            eps = 1.0e-9
            ymin = float(np.min(finite))
            ymax = float(np.max(finite))
            ymin_pad = 1.1 * ymin if abs(ymin) > eps else -eps
            ymax_pad = 1.1 * ymax if abs(ymax) > eps else eps
            if np.all(finite >= -eps):
                ax.set_ylim(0.0, ymax_pad)
            elif np.all(finite <= eps):
                ax.set_ylim(ymin_pad, 0.0)
            else:
                ax.set_ylim(ymin_pad, ymax_pad)
                ax.axhline(0.0, color="black", linestyle="-", linewidth=0.8)
        if mvals.size > 1:
            names = [f"m = {int(v)}" for v in mvals]
            plot_legend(names, colors=colors[1 : 1 + mvals.size], linestyles=[linestyle] * mvals.size)
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
                j = int(np.argmin(np.abs(m - q_input[i] * ntor))) if q_input.size > i else int(np.argmin(np.abs(m)))
                k = int(np.argmin(np.abs(m + q_input[i] * ntor))) if q_input.size > i else j
                q_here = float(np.interp(idxf, np.arange(qprof.size, dtype=float), np.abs(qprof)))
                psi_here = float(np.interp(idxf, np.arange(nflux.size, dtype=float), nflux))
                print("q, Psi = ", q_here, psi_here)
                dj = np.interp(idxf, np.arange(nflux.size, dtype=float), d[0, j, :])
                dk = np.interp(idxf, np.arange(nflux.size, dtype=float), d[0, k, :])
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
    if n.size > 1:
        matches = np.where(n == int(ntor))[0]
        if matches.size != 1:
            print(f"schaffer_plot: error: ntor = {ntor} not found.")
            return None, None
        k = int(matches[0])
        print(f"plotting ntor = {ntor}")
    else:
        k = 0

    contour_and_legend(
        np.abs(d[k, :, :]),
        m,
        y,
        lines=True,
        label=label,
        xtitle="m",
        ytitle=ytitle,
        overplot=overplot,
    )
    ax = plt.gca()
    ax.set_xlim(-20, 20)
    ax.set_ylim(0, 1)
    ax.plot(ntor * qprof, y, linestyle="--", color="black")
    return ax.figure, ax
