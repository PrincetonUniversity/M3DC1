from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np

from .contour_and_legend import contour_and_legend
from .field_at_point import field_at_point
from .plot_legend import plot_legend
from .radius_matrix import radius_matrix
from .read_field import read_field
from .read_gamma import read_gamma
from .read_parameter import read_parameter


@dataclass
class EquationResult:
    terms: np.ndarray
    names: list[str]
    title: str
    r: np.ndarray
    z: np.ndarray


def _field(name: str, *, filename="C1.h5", timeslices=0, points=200, return_meta=False, **kwargs):
    return read_field(name, timeslices=timeslices, filename=filename, points=points, return_meta=return_meta, **kwargs)


def _data(name: str, *, filename="C1.h5", timeslices=0, points=200, **kwargs):
    return np.asarray(_field(name, filename=filename, timeslices=timeslices, points=points, **kwargs))


def _first_data(name: str, *, filename="C1.h5", timeslices=0, points=200, **kwargs):
    meta = _field(name, filename=filename, timeslices=timeslices, points=points, return_meta=True, **kwargs)
    data = np.asarray(meta.data)
    if data.ndim == 3:
        data = data[0]
    return data, np.asarray(meta.r, dtype=float), np.asarray(meta.z, dtype=float)


def _as2(a):
    arr = np.asarray(a)
    return arr[0] if arr.ndim == 3 else arr


def _zero_like(a):
    return np.zeros_like(np.asarray(a))


def _params(filename):
    return {
        "icomplex": int(read_parameter("icomplex", filename=filename)),
        "itor": int(read_parameter("itor", filename=filename)),
        "ntor": int(read_parameter("ntor", filename=filename)),
        "rzero": float(read_parameter("rzero", filename=filename)),
        "numvar": int(read_parameter("numvar", filename=filename)),
        "ivform": int(read_parameter("ivform", filename=filename)),
    }


def _rfac_and_radius(x, z, p, *, filename):
    if p["itor"] == 1:
        return radius_matrix(x, z, filename=filename), 1j * p["ntor"]
    return 1.0, 1j * p["ntor"] / p["rzero"]


def _require_complex(p, name):
    if p["icomplex"] != 1:
        raise NotImplementedError(f"{name} currently implements the icomplex=1 branch from IDL.")


def eqn_continuity(*, filename="C1.h5", timeslices=0, points=200, **kwargs) -> EquationResult:
    title = "Continuity Equation"
    names = ["dn/dt", "v.Grad(n)", "n Div(v)", "-D del^2(n)"]
    p = _params(filename)
    print("numvar = ", p["numvar"])
    print("itor = ", p["itor"])
    _require_complex(p, "eqn_continuity")

    gamma = read_gamma(filename=filename)
    n0, x, z = _first_data("den", filename=filename, timeslices=timeslices, points=points, equilibrium=True, **kwargs)
    n0_r = _as2(_data("den", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=2, **kwargs))
    n0_z = _as2(_data("den", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=3, **kwargs))
    w0 = _as2(_data("omega", filename=filename, timeslices=timeslices, points=points, equilibrium=True, **kwargs))

    n1 = _as2(_data("den", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, **kwargs))
    n1_lp = _as2(_data("den", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=7, **kwargs))
    n1_r = _as2(_data("den", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=2, **kwargs)) if p["itor"] == 1 else _zero_like(n1)
    u1_r = _as2(_data("phi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=2, **kwargs))
    u1_z = _as2(_data("phi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=3, **kwargs))
    w1 = _as2(_data("omega", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, **kwargs)) if p["numvar"] >= 2 else _zero_like(n1)

    if p["numvar"] == 3:
        chi1_r = _as2(_data("chi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=2, **kwargs))
        chi1_z = _as2(_data("chi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=3, **kwargs))
        chi1_lp = _as2(_data("chi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=7, **kwargs))
    else:
        chi1_r = _zero_like(n1)
        chi1_z = _zero_like(n1)
        chi1_lp = _zero_like(n1)

    r, rfac = _rfac_and_radius(x, z, p, filename=filename)
    if p["itor"] == 1:
        n1_lp = n1_lp - n1_r / r

    terms = np.zeros((4,) + n0.shape, dtype=np.complex128)
    print("defining dn/dt")
    terms[0] = n1 * gamma[0]
    print("defining v.Grad(n)")
    if p["ivform"] == 1:
        terms[1] = r * (n0_z * u1_r - n0_r * u1_z) + w0 * rfac * n1 + (n0_r * chi1_r + n0_z * chi1_z) / r**2
    else:
        terms[1] = (n0_z * u1_r - n0_r * u1_z) / r + w0 * rfac * n1 + n0_r * chi1_r + n0_z * chi1_z
    print("defining n Div(v)")
    if p["ivform"] == 1:
        terms[2] = n0 * (rfac * w1 + chi1_lp / r**2 + (p["itor"] == 1) * (-2.0 * u1_z - 2.0 * chi1_r / r**3))
    else:
        terms[2] = n0 * (rfac * w1 + chi1_lp)
    print("defining del^2(n)")
    terms[3] = -float(read_parameter("denm", filename=filename)) * n1_lp
    print("gamma = ", gamma[0])
    return EquationResult(terms=terms, names=names, title=title, r=x, z=z)


def eqn_gradshafranov(*, filename="C1.h5", timeslices=0, points=200, **kwargs) -> EquationResult:
    title = "Grad Shafranov Equation"
    names = ["del*(psi)", "R^2 p'", "FF'"]
    itor = int(read_parameter("itor", filename=filename))
    psi_r, x, z = _first_data("psi", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=2, **kwargs)
    psi_z = _as2(_data("psi", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=3, **kwargs))
    p_r = _as2(_data("p", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=2, **kwargs))
    p_z = _as2(_data("p", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=3, **kwargs))
    i_r = _as2(_data("i", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=2, **kwargs))
    i_z = _as2(_data("i", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=3, **kwargs))
    i0 = _as2(_data("i", filename=filename, timeslices=timeslices, points=points, equilibrium=True, **kwargs))
    jy = _as2(_data("jy", filename=filename, timeslices=timeslices, points=points, equilibrium=True, **kwargs))
    r = radius_matrix(x, z, filename=filename) if itor == 1 else 1.0
    den = psi_r**2 + psi_z**2

    terms = np.zeros((3,) + i0.shape, dtype=np.complex128)
    print("defining del*(psi)")
    terms[0] = -r * jy
    print("defining p'")
    terms[1] = r**2 * (p_r * psi_r + p_z * psi_z) / den
    print("defining FF'")
    terms[2] = i0 * (i_r * psi_r + i_z * psi_z) / den
    return EquationResult(terms=terms, names=names, title=title, r=x, z=z)


def eqn_momentum_x(*, filename="C1.h5", timeslices=0, points=200, **kwargs) -> EquationResult:
    return _eqn_momentum("x", filename=filename, timeslices=timeslices, points=points, **kwargs)


def eqn_momentum_y(*, filename="C1.h5", timeslices=0, points=200, **kwargs) -> EquationResult:
    return _eqn_momentum("y", filename=filename, timeslices=timeslices, points=points, **kwargs)


def eqn_momentum_z(*, filename="C1.h5", timeslices=0, points=200, **kwargs) -> EquationResult:
    return _eqn_momentum("z", filename=filename, timeslices=timeslices, points=points, **kwargs)


def _eqn_momentum(axis: str, *, filename="C1.h5", timeslices=0, points=200, **kwargs) -> EquationResult:
    titles = {"x": "Momentum Equation: R", "y": "Momentum Equation: Phi", "z": "Momentum Equation: Z"}
    names = ["dV/dt", "-JxB", "Grad(p)"]
    p = _params(filename)
    _require_complex(p, f"eqn_momentum_{axis}")
    gamma = read_gamma(filename=filename)

    den0, x, z = _first_data("den", filename=filename, timeslices=timeslices, points=points, equilibrium=True, **kwargs)
    r, rfac = _rfac_and_radius(x, z, p, filename=filename)

    psi0_r = _as2(_data("psi", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=2, **kwargs))
    psi0_z = _as2(_data("psi", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=3, **kwargs))
    psi0_gs = _as2(_data("psi", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=7, **kwargs))
    i0 = _as2(_data("i", filename=filename, timeslices=timeslices, points=points, equilibrium=True, **kwargs))
    i0_r = _as2(_data("i", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=2, **kwargs))
    i0_z = _as2(_data("i", filename=filename, timeslices=timeslices, points=points, equilibrium=True, operation=3, **kwargs))

    psi1_r = _as2(_data("psi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=2, **kwargs))
    psi1_z = _as2(_data("psi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=3, **kwargs))
    psi1_gs = _as2(_data("psi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=7, **kwargs))

    if p["numvar"] >= 2:
        i1 = _as2(_data("i", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, **kwargs))
        w1 = _as2(_data("omega", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, **kwargs))
        i1_r = _as2(_data("i", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=2, **kwargs))
        i1_z = _as2(_data("i", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=3, **kwargs))
        f1_r = _as2(_data("f", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=2, **kwargs))
        f1_z = _as2(_data("f", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=3, **kwargs))
    else:
        i1 = _zero_like(den0)
        w1 = _zero_like(den0)
        i1_r = _zero_like(den0)
        i1_z = _zero_like(den0)
        f1_r = _zero_like(den0)
        f1_z = _zero_like(den0)

    if p["numvar"] == 3:
        p1 = _as2(_data("p", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, **kwargs))
        p1_r = _as2(_data("p", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=2, **kwargs))
        p1_z = _as2(_data("p", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=3, **kwargs))
        chi1_r = _as2(_data("chi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=2, **kwargs))
        chi1_z = _as2(_data("chi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=3, **kwargs))
    else:
        p1 = _zero_like(den0)
        p1_r = _zero_like(den0)
        p1_z = _zero_like(den0)
        chi1_r = _zero_like(den0)
        chi1_z = _zero_like(den0)

    if p["itor"] == 1:
        psi0_gs = psi0_gs - psi0_r / r
        psi1_gs = psi1_gs - psi1_r / r

    bx0 = -psi0_z / r
    by0 = i0 / r
    bz0 = psi0_r / r
    jx0 = -i0_z / r
    jy0 = -psi0_gs / r
    jz0 = i0_r / r
    bx1 = -psi1_z / r - f1_r * rfac
    by1 = i1 / r
    bz1 = psi1_r / r - f1_z * rfac
    jx1 = -(i1_z + f1_z * rfac**2) / r + rfac * psi1_r / r**2
    jy1 = -psi1_gs / r
    jz1 = (i1_r + f1_r * rfac**2) / r + rfac * psi1_z / r**2

    terms = np.zeros((3,) + den0.shape, dtype=np.complex128)
    print("defining dv/dt")
    if axis == "x":
        vx1 = -r * _as2(_data("phi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=3, **kwargs)) + chi1_r / r**2 if p["ivform"] == 1 else -_as2(_data("phi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=3, **kwargs)) / r + chi1_r / r**2
        terms[0] = den0 * vx1 * gamma[0]
        print("defining JxB")
        terms[1] = -(jy1 * bz0 - jz1 * by0 + jy0 * bz1 - jz0 * by1)
        print("defining grad(p)")
        terms[2] = p1_r
    elif axis == "y":
        vy1 = w1 * r
        terms[0] = den0 * vy1 * gamma[0]
        print("defining JxB")
        terms[1] = -(jz1 * bx0 - jx1 * bz0 + jz0 * bx1 - jx0 * bz1)
        print("defining grad(p)")
        terms[2] = rfac * p1
    else:
        u1_r = _as2(_data("phi", filename=filename, timeslices=timeslices, points=points, linear=True, complex=True, operation=2, **kwargs))
        vz1 = r * u1_r + chi1_z / r**2 if p["ivform"] == 1 else u1_r / r + chi1_z / r**2
        terms[0] = den0 * vz1 * gamma[0]
        print("defining JxB")
        terms[1] = -(jx1 * by0 - jy1 * bx0 + jx0 * by1 - jy0 * bx1)
        print("defining grad(p)")
        terms[2] = p1_z

    print("gamma = ", gamma[0])
    return EquationResult(terms=terms, names=names, title=titles[axis], r=x, z=z)


_EQUATIONS: dict[str, Callable[..., EquationResult]] = {
    "continuity": eqn_continuity,
    "gradshafranov": eqn_gradshafranov,
    "momentum_x": eqn_momentum_x,
    "momentum_y": eqn_momentum_y,
    "momentum_z": eqn_momentum_z,
}


def _apply_func(values, func):
    if func is None or func == "":
        func = "real_part"
    if callable(func):
        return func(values)
    name = str(func).lower()
    if name in {"real_part", "real", "re"}:
        return np.real(values)
    if name in {"imaginary", "imag_part", "imag", "im"}:
        return np.imag(values)
    if name in {"abs", "absolute"}:
        return np.abs(values)
    raise ValueError(f"Unsupported func={func!r}. Use real_part, imaginary, abs, or a callable.")


def _finite_minmax(arr):
    vals = np.asarray(arr, dtype=float)
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        return 0.0, 1.0
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    if vmin == vmax:
        pad = max(abs(vmin), 1.0) * 0.05
        return vmin - pad, vmax + pad
    return vmin, vmax


def plot_equation(
    equation: str,
    *,
    cutz=None,
    cutr=None,
    func="real_part",
    outfile: str | Path | None = None,
    filename="C1.h5",
    timeslices=0,
    points: int = 200,
    **kwargs,
):
    """Port of plot_equation.pro.

    Returns ``(figure, axis)`` for the line plot.  A contour figure for
    ``abs(total)`` is also created, matching the IDL routine.
    """
    key = str(equation).lower()
    if key.startswith("eqn_"):
        key = key[4:]
    if key not in _EQUATIONS:
        raise ValueError(f"Unknown equation {equation!r}. Available: {', '.join(sorted(_EQUATIONS))}")

    result = _EQUATIONS[key](filename=filename, timeslices=timeslices, points=points, **kwargs)
    term = np.asarray(result.terms)
    total = np.sum(term, axis=0)
    term_plot = _apply_func(term, func)
    total_plot = _apply_func(total, func)

    if cutr is not None:
        r0 = float(cutr)
        x0 = np.full(result.z.size, r0)
        z0 = result.z
        xtitle = "Z"
        xdat = result.z
        contour_line = ("v", r0)
    else:
        zcut = 0.0 if cutz is None else float(cutz)
        x0 = result.r
        z0 = np.full(result.r.size, zcut)
        xtitle = "R"
        xdat = result.r
        contour_line = ("h", zcut)

    f = np.zeros((term.shape[0] + 1, xdat.size), dtype=float)
    for i in range(term.shape[0]):
        f[i, :] = np.asarray(field_at_point(term_plot[i], result.r, result.z, x0, z0), dtype=float)
    f[term.shape[0], :] = np.asarray(field_at_point(total_plot, result.r, result.z, x0, z0), dtype=float)

    fig, ax = plt.subplots()
    ymin, ymax = _finite_minmax(f)
    ax.plot([float(np.min(xdat)), float(np.max(xdat))], [ymin, ymax], alpha=0.0)
    ax.set_title(result.title)
    ax.set_xlabel(xtitle)

    colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    for i in range(term.shape[0]):
        ax.plot(xdat, f[i, :], color=colors[(i + 1) % len(colors)] if colors else None)
    ax.plot(xdat, f[term.shape[0], :], color=colors[0] if colors else None, linestyle="--")
    names = ["Total"] + result.names
    plot_legend(names, colors=[colors[0] if colors else None] + [colors[(i + 1) % len(colors)] if colors else None for i in range(term.shape[0])])

    contour_and_legend(np.abs(total), result.r, result.z, title="", **kwargs)
    contour_fig = plt.gcf()
    ax2 = plt.gca()
    if contour_line[0] == "h":
        ax2.plot([float(np.min(result.r)), float(np.max(result.r))], [contour_line[1], contour_line[1]], color="k")
    else:
        ax2.plot([contour_line[1], contour_line[1]], [float(np.min(result.z)), float(np.max(result.z))], color="k")

    if outfile is not None:
        out = Path(outfile)
        with out.open("w", encoding="ascii") as handle:
            header = [xtitle, *names]
            handle.write("".join(f"{h:>15s}" for h in header) + "\n")
            for k in range(xdat.size):
                dat = [xdat[k], f[term.shape[0], k], *[f[i, k] for i in range(term.shape[0])]]
                print(np.asarray(dat, dtype=float))
                handle.write("".join(f"{v:15.5G}" for v in dat) + "\n")

    return fig, ax
