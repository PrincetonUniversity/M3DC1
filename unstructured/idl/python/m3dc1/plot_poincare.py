from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .make_label import make_label
from .plot_coils import plot_coils
from .plot_lcfs import plot_lcfs
from .plot_mesh import plot_mesh
from .plot_wall_regions import plot_wall_regions
from .read_poincare import PoincareResult, read_poincare


def plot_poincare(
    files=None,
    *,
    directory: str | Path = ".",
    pattern: str = "out*",
    filename: str | Path = "C1.h5",
    max_files: int | None = None,
    rcol: int = 1,
    zcol: int = 2,
    phicol: int | None = 0,
    valuecol: int | None = None,
    marker: str = "o",
    markersize: float | None = None,
    linestyle: str = "None",
    color="black",
    c=None,
    cmap: str | None = None,
    colorbar: bool = False,
    label: str | None = None,
    title: str | None = None,
    xrange=None,
    yrange=None,
    xlim=None,
    ylim=None,
    iso: bool = False,
    overplot: bool = False,
    mesh=None,
    boundary: bool = False,
    lcfs: bool = False,
    coils: bool = False,
    wall_regions: bool = False,
    logical: bool = False,
    points: int = 200,
    slice: int = 0,
    phi: float = 0.0,
    xscale: float = 1.0,
    yscale: float = 1.0,
    outfile: str | Path | None = None,
    mpeg: str | Path | None = None,
    skip_empty: bool = True,
    **kwargs,
):
    """Plot Poincare points from M3D-C1 ``out*`` text files."""
    xscale_f = float(xscale)
    yscale_f = float(yscale)
    markersize_f = 0.1 if markersize is None and overplot else (0.5 if markersize is None else float(markersize))

    read_valuecol = 4 if colorbar and c is None and valuecol is None else valuecol

    if isinstance(files, PoincareResult):
        meta = files
    else:
        meta = read_poincare(
            files,
            directory=directory,
            pattern=pattern,
            rcol=rcol,
            zcol=zcol,
            phicol=phicol,
            valuecol=read_valuecol,
            max_files=max_files,
            skip_empty=skip_empty,
            return_meta=True,
        )

    per_file_color = isinstance(color, bool) and color
    point_color = None if per_file_color else color

    r_all = np.asarray(meta.r, dtype=float) * xscale_f
    z_all = np.asarray(meta.z, dtype=float) * yscale_f
    finite = np.isfinite(r_all) & np.isfinite(z_all)
    r = r_all[finite]
    z = z_all[finite]

    cdata = None
    if c is not None:
        raw_c = np.asarray(c, dtype=float).reshape(-1)
        if raw_c.size != finite.size:
            raise ValueError(f"c has {raw_c.size} values, expected {finite.size}.")
        cdata = raw_c[finite]
    elif colorbar:
        if meta.value.size != finite.size:
            raise ValueError("colorbar=True requires a value column or explicit c values.")
        cdata = np.asarray(meta.value, dtype=float)[finite]

    if outfile is not None:
        np.savetxt(str(outfile), np.column_stack([r, z]), fmt="%16.6e")

    if overplot:
        ax = plt.gca()
        fig = ax.figure
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlabel(make_label("R", l0=1, **kwargs))
        ax.set_ylabel(make_label("Z", l0=1, **kwargs))

    if title is None and not overplot:
        title = "Poincare"
    if title:
        ax.set_title(str(title))

    if r.size:
        if cdata is not None:
            artist = ax.scatter(r, z, s=float(markersize), marker=marker, c=cdata, cmap=cmap)
            if colorbar:
                fig.colorbar(artist, ax=ax, label=label or "")
        elif per_file_color:
            cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
            if not cycle:
                cycle = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
            for i, arr in enumerate(meta.data):
                a = np.asarray(arr, dtype=float)
                rr = a[:, int(rcol)] * xscale_f
                zz = a[:, int(zcol)] * yscale_f
                ok = np.isfinite(rr) & np.isfinite(zz)
                if np.any(ok):
                    ax.plot(
                        rr[ok],
                        zz[ok],
                        marker=marker,
                        markersize=markersize_f,
                        linestyle=linestyle,
                        color=cycle[i % len(cycle)],
                    )
        else:
            ax.plot(
                r,
                z,
                marker=marker,
                markersize=markersize_f,
                linestyle=linestyle,
                color=point_color,
                label=label,
            )

    mesh_obj = None if isinstance(mesh, bool) else mesh
    show_mesh = bool(mesh) if isinstance(mesh, bool) else (mesh is not None)
    if boundary:
        plot_mesh(
            mesh=mesh_obj,
            oplot=True,
            boundary=True,
            logical=logical,
            phi=phi,
            xscale=xscale_f,
            yscale=yscale_f,
            filename=filename,
            slice=slice,
            points=points,
            **kwargs,
        )
    elif show_mesh:
        plot_mesh(
            mesh=mesh_obj,
            oplot=True,
            boundary=False,
            logical=logical,
            phi=phi,
            xscale=xscale_f,
            yscale=yscale_f,
            filename=filename,
            slice=slice,
            points=points,
            **kwargs,
        )

    if wall_regions:
        plot_wall_regions(filename=filename, slice=slice, over=True, xscale=xscale_f, yscale=yscale_f, **kwargs)
    if lcfs:
        plot_lcfs(overplot=True, filename=filename, slice=slice, points=points, xscale=xscale_f, yscale=yscale_f, **kwargs)
    if coils:
        plot_coils(filename=filename, overplot=True, xscale=xscale_f, yscale=yscale_f, **kwargs)

    if label and cdata is None:
        ax.legend()
    if xrange is not None:
        ax.set_xlim(xrange)
    if yrange is not None:
        ax.set_ylim(yrange)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if iso:
        ax.set_aspect("equal", adjustable="box")
    if not overplot:
        fig.tight_layout()
    if mpeg is not None:
        fig.savefig(str(mpeg), dpi=150)
    return fig, ax
