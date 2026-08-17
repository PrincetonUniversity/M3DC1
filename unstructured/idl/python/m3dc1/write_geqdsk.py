from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .field_at_point import field_at_point
from .flux_average_field import flux_average_field
from .flux_coordinates import flux_coordinates
from .get_normalizations import get_normalizations
from .lcfs import lcfs
from .radius_matrix import radius_matrix
from .read_field import read_field
from .read_mesh import read_mesh
from .read_parameter import read_parameter
from .read_scalar import read_scalar
from .s_bracket import s_bracket


def _as_2d(a) -> np.ndarray:
    arr = np.asarray(a, dtype=float)
    if arr.ndim == 3:
        return arr[0, :, :]
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D field data, got shape {arr.shape}.")
    return arr


def _point_in_poly(path: np.ndarray, point: np.ndarray) -> bool:
    x = np.asarray(path[:, 0], dtype=float)
    y = np.asarray(path[:, 1], dtype=float)
    px = float(point[0])
    py = float(point[1])
    inside = False
    j = x.size - 1
    for i in range(x.size):
        dy = y[j] - y[i]
        if ((y[i] > py) != (y[j] > py)) and (
            px < (x[j] - x[i]) * (py - y[i]) / (dy if abs(dy) > np.finfo(float).tiny else np.finfo(float).tiny) + x[i]
        ):
            inside = not inside
        j = i
    return inside


def _path_at_flux(psi: np.ndarray, x: np.ndarray, z: np.ndarray, flux: float, axis: np.ndarray) -> np.ndarray:
    fig, ax = plt.subplots(figsize=(2, 2))
    try:
        cs = ax.contour(x, z, psi.T, levels=[float(flux)])
        paths: list[np.ndarray] = []
        if hasattr(cs, "collections") and cs.collections:
            for path in cs.collections[0].get_paths():
                vertices = np.asarray(path.vertices, dtype=float)
                if vertices.shape[0] > 4:
                    paths.append(vertices)
        elif hasattr(cs, "allsegs") and cs.allsegs:
            for seg in cs.allsegs[0]:
                vertices = np.asarray(seg, dtype=float)
                if vertices.shape[0] > 4:
                    paths.append(vertices)
        if not paths:
            print("Error: no points at this flux value", flux)
            return np.zeros((0, 2), dtype=float)
        for path in paths:
            if _point_in_poly(path, axis):
                return path
        return paths[int(np.argmax([p.shape[0] for p in paths]))]
    finally:
        plt.close(fig)


def _mesh_vertices(el: np.ndarray, i: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    a = float(el[0, i])
    b = float(el[1, i])
    c = float(el[2, i])
    theta = float(el[3, i])
    x0 = float(el[4, i])
    z0 = float(el[5, i])
    p1 = np.asarray([x0, z0], dtype=float)
    p2 = p1 + np.asarray([(b + a) * np.cos(theta), (b + a) * np.sin(theta)], dtype=float)
    p3 = p1 + np.asarray([b * np.cos(theta) - c * np.sin(theta), b * np.sin(theta) + c * np.cos(theta)], dtype=float)
    return p1, p2, p3


def _edge_zone(bound: int, edge: int) -> int:
    if edge == 1:
        return ((bound & 120) >> 3) + 1
    if edge == 2:
        return ((bound & 1920) >> 7) + 1
    return ((bound & 30720) >> 11) + 1


def _is_boundary_edge(bound: int, edge: int, imultiregion: int) -> bool:
    if not (bound & edge):
        return False
    zone = _edge_zone(bound, edge)
    if int(imultiregion) == 1:
        return zone in (1, 2)
    return zone != 0


def _get_boundary_path(*, filename: str | Path, timeslice: int) -> np.ndarray:
    try:
        imultiregion = int(read_parameter("imulti_region", filename=filename))
    except Exception:
        imultiregion = 0

    mesh = read_mesh(filename=filename, slice=timeslice)
    el = np.asarray(mesh.elements, dtype=float)
    if el.ndim != 2 or el.shape[0] < 7:
        return np.zeros((0, 2), dtype=float)

    edges: list[tuple[np.ndarray, np.ndarray]] = []
    for i in range(mesh.nelms):
        bound = int(el[6, i])
        p1, p2, p3 = _mesh_vertices(el, i)
        if _is_boundary_edge(bound, 1, imultiregion):
            edges.append((p1, p2))
        if _is_boundary_edge(bound, 2, imultiregion):
            edges.append((p2, p3))
        if _is_boundary_edge(bound, 4, imultiregion):
            edges.append((p3, p1))

    print("Number of boundary points: ", len(edges))
    if not edges:
        return np.zeros((0, 2), dtype=float)

    xy = _order_boundary_edges(edges)
    print("found ", xy.shape[0], " points")
    center = np.asarray([(np.nanmax(xy[:, 0]) + np.nanmin(xy[:, 0])) / 2.0, (np.nanmax(xy[:, 1]) + np.nanmin(xy[:, 1])) / 2.0])
    angle = np.arctan2(xy[:, 1] - center[1], xy[:, 0] - center[0])
    if np.nanmedian(np.gradient(angle)) < 0.0:
        print("reversing angle!")
        xy = xy[::-1, :]
        angle = angle[::-1]
    angle = np.where(angle < 0.0, angle + 2.0 * np.pi, angle)
    i0 = int(np.nanargmin(angle))
    return np.roll(xy, -i0 - 1, axis=0)


def _order_boundary_edges(edges: list[tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    tol = 1.0e-6

    def key(p: np.ndarray) -> tuple[int, int]:
        return tuple(np.rint(np.asarray(p, dtype=float) / tol).astype(np.int64))

    coords: dict[tuple[int, int], np.ndarray] = {}
    adjacency: dict[tuple[int, int], list[tuple[int, int]]] = {}
    edge_keys: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for p0, p1 in edges:
        k0 = key(p0)
        k1 = key(p1)
        coords.setdefault(k0, np.asarray(p0, dtype=float))
        coords.setdefault(k1, np.asarray(p1, dtype=float))
        adjacency.setdefault(k0, []).append(k1)
        adjacency.setdefault(k1, []).append(k0)
        edge_keys.add(tuple(sorted((k0, k1))))

    def walk(start: tuple[int, int], nxt: tuple[int, int], used: set[tuple[tuple[int, int], tuple[int, int]]]) -> list[tuple[int, int]]:
        path = [start]
        prev = start
        cur = nxt
        while True:
            ekey = tuple(sorted((prev, cur)))
            if ekey in used:
                break
            used.add(ekey)
            if cur == start:
                break
            path.append(cur)
            candidates = [k for k in adjacency.get(cur, []) if tuple(sorted((cur, k))) not in used]
            if not candidates:
                break
            non_backtracking = [k for k in candidates if k != prev]
            prev, cur = cur, (non_backtracking[0] if non_backtracking else candidates[0])
        return path

    used_global: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    paths: list[list[tuple[int, int]]] = []
    for ekey in list(edge_keys):
        if ekey in used_global:
            continue
        local_used = set(used_global)
        path = walk(ekey[0], ekey[1], local_used)
        used_global.update(local_used)
        if len(path) > 1:
            paths.append(path)

    if not paths:
        return np.asarray([edges[0][0], edges[0][1]], dtype=float)
    best = max(paths, key=len)
    xy = np.asarray([coords[k] for k in best], dtype=float)
    if xy.shape[0] < len(edges) // 2:
        print("WARNING: boundary path is not closed")
    return xy


def _reduce_points(r: np.ndarray, z: np.ndarray, limit: int = 500) -> tuple[np.ndarray, np.ndarray]:
    rr = np.asarray(r, dtype=float).reshape(-1)
    zz = np.asarray(z, dtype=float).reshape(-1)
    n = rr.size
    while n >= limit:
        print("reducing lim points...")
        n = n // 2
        rr = rr[0 : 2 * n : 2]
        zz = zz[0 : 2 * n : 2]
        print("new lim points = ", n)
    return rr, zz


def _psi_at_axis_or_extremum(psi: np.ndarray, x: np.ndarray, z: np.ndarray, axis: np.ndarray) -> tuple[float, np.ndarray]:
    val = np.asarray(field_at_point(psi, x, z, axis[0], axis[1]), dtype=float).reshape(-1)
    if val.size and np.isfinite(val[0]):
        return float(val[0]), np.asarray(axis, dtype=float)
    idx = np.unravel_index(int(np.nanargmin(psi)), psi.shape)
    return float(psi[idx]), np.asarray([x[idx[0]], z[idx[1]]], dtype=float)


def _valid_flux_value(value: float, psi: np.ndarray) -> bool:
    if not np.isfinite(value):
        return False
    pmin = float(np.nanmin(psi))
    pmax = float(np.nanmax(psi))
    tol = 1.0e-8 * max(1.0, abs(pmax - pmin))
    return pmin - tol <= float(value) <= pmax + tol


def _fallback_psilim(psi: np.ndarray, flux0: float) -> float:
    pmin = float(np.nanmin(psi))
    pmax = float(np.nanmax(psi))
    return pmax if abs(flux0 - pmin) <= abs(flux0 - pmax) else pmin


def _profile_by_psi_bins(field: np.ndarray, psi: np.ndarray, flux: np.ndarray) -> np.ndarray:
    vals = np.asarray(field, dtype=float)
    ps = np.asarray(psi, dtype=float)
    fl = np.asarray(flux, dtype=float).reshape(-1)
    out = np.full(fl.size, np.nan, dtype=float)
    if fl.size == 0:
        return out
    edges = np.empty(fl.size + 1, dtype=float)
    if fl.size == 1:
        span = max(abs(fl[0]) * 0.01, 1.0)
        edges[0] = fl[0] - span
        edges[1] = fl[0] + span
    else:
        edges[1:-1] = 0.5 * (fl[:-1] + fl[1:])
        edges[0] = fl[0] - 0.5 * (fl[1] - fl[0])
        edges[-1] = fl[-1] + 0.5 * (fl[-1] - fl[-2])
    lo = np.minimum(edges[:-1], edges[1:])
    hi = np.maximum(edges[:-1], edges[1:])
    finite = np.isfinite(vals) & np.isfinite(ps)
    for i in range(fl.size):
        mask = finite & (ps >= lo[i]) & (ps < hi[i])
        if i == fl.size - 1:
            mask = finite & (ps >= lo[i]) & (ps <= hi[i])
        if np.any(mask):
            out[i] = float(np.nanmean(vals[mask]))
    if np.any(np.isfinite(out)):
        good = np.flatnonzero(np.isfinite(out))
        bad = np.flatnonzero(~np.isfinite(out))
        if bad.size:
            out[bad] = np.interp(bad, good, out[good])
    return out


def _write_e16(handle, values) -> None:
    vals = np.asarray(values, dtype=float).reshape(-1)
    for i in range(0, vals.size, 5):
        handle.write("".join(f"{v:16.9E}" for v in vals[i : i + 5]) + "\n")


def _write_i5(handle, values) -> None:
    vals = np.asarray(values, dtype=int).reshape(-1)
    handle.write("".join(f"{v:5d}" for v in vals) + "\n")


def _write_geqdsk_header(handle, nr: int, nz: int) -> None:
    name = ["M3DC1", "03/17/", "2013    ", "#000000", "0000", ""]
    handle.write("".join(f"{s[:8]:<8}" for s in name) + f"{3:4d}{nr:4d}{nz:4d}\n")


def write_geqdsk(
    *,
    eqfile: str | Path = "geqdsk.out",
    psilim: float | None = None,
    filename: str | Path = "C1.h5",
    timeslice: int | None = None,
    points: int = 200,
    jsolver: bool = False,
    jsfile: str | Path = "jsfile",
    **kwargs,
) -> Path:
    """Write the equilibrium from a C1 HDF5 file in GEQDSK format.

    This follows the GEQDSK output path in ``write_geqdsk.pro``.  Values are
    converted to SI units before writing.
    """
    if timeslice is None:
        timeslice = int(read_parameter("ntime", filename=filename)) - 1
    eqfile = Path(eqfile)

    b0, _n0, l0, _mi = get_normalizations(filename=filename)
    bzero = float(read_parameter("bzero", filename=filename))
    rzero = float(read_parameter("rzero", filename=filename))

    meta = read_field("psi", timeslices=timeslice, equilibrium=True, points=points, filename=filename, return_meta=True, **kwargs)
    psi = _as_2d(meta.data)
    x = np.asarray(meta.r, dtype=float).reshape(-1)
    z = np.asarray(meta.z, dtype=float).reshape(-1)
    r = np.asarray(radius_matrix(x, z, filename=filename), dtype=float)

    boundary = _get_boundary_path(filename=filename, timeslice=timeslice)
    rwall = boundary[:, 0] if boundary.size else np.asarray([x[0], x[-1], x[-1], x[0]], dtype=float)
    zwall = boundary[:, 1] if boundary.size else np.asarray([z[0], z[0], z[-1], z[-1]], dtype=float)
    nwall = rwall.size

    ifixedb = int(read_parameter("ifixedb", filename=filename))
    print("ifixedb = ", ifixedb)

    lc = lcfs(psi, x, z, filename=filename, slice=timeslice, **kwargs)
    lcfs_psi = float(lc.psilim)
    flux0 = float(lc.flux0)
    axis = np.asarray(lc.axis, dtype=float).copy()
    xpoint = np.asarray(lc.xpoint, dtype=float)

    if ifixedb == 1:
        psilim_val = 0.0
        lcfs_psi = 0.0
        rlim = np.asarray(rwall, dtype=float).copy()
        zlim = np.asarray(zwall, dtype=float).copy()
    else:
        psilim_val = lcfs_psi if psilim is None else float(psilim)
        flux0_grid, axis_grid = _psi_at_axis_or_extremum(psi, x, z, axis)
        if not _valid_flux_value(flux0, psi):
            print("WARNING: LCFS flux0 is outside psi grid; using psi-grid magnetic-axis value")
            flux0 = flux0_grid
            axis = axis_grid
        if psilim is None and (not _valid_flux_value(psilim_val, psi) or abs(psilim_val - flux0) < np.finfo(float).eps):
            print("WARNING: LCFS psilim is outside psi grid; using psi-grid boundary value")
            psilim_val = _fallback_psilim(psi, flux0)
        print("lcfs_psi = ", lcfs_psi)
        print("psilim = ", psilim_val)
        lcfs_path = _path_at_flux(psi, x, z, psilim_val, axis)
        if lcfs_path.size == 0:
            print("WARNING: no contour found at psilim; using wall points as limiter boundary")
            lcfs_path = np.column_stack((rwall, zwall))
        if xpoint.size > 1 and (xpoint[0] != 0.0 or xpoint[1] != 0.0):
            if xpoint[1] < axis[1]:
                mask = lcfs_path[:, 1] > xpoint[1]
            else:
                mask = lcfs_path[:, 1] < xpoint[1]
        else:
            mask = np.ones(lcfs_path.shape[0], dtype=bool)
        rlim = lcfs_path[mask, 0]
        zlim = lcfs_path[mask, 1]
        print("lim points = ", rlim.size, int(np.count_nonzero(mask)))

    rlim, zlim = _reduce_points(rlim, zlim)
    nlim = rlim.size

    psi_r = _as_2d(read_field("psi", timeslices=timeslice, equilibrium=True, operation=2, points=points, filename=filename, **kwargs))
    psi_z = _as_2d(read_field("psi", timeslices=timeslice, equilibrium=True, operation=3, points=points, filename=filename, **kwargs))
    psi_lp = _as_2d(read_field("psi", timeslices=timeslice, equilibrium=True, operation=7, points=points, filename=filename, **kwargs))
    p0 = _as_2d(read_field("p", timeslices=timeslice, equilibrium=True, points=points, filename=filename, **kwargs))
    p0_r = _as_2d(read_field("p", timeslices=timeslice, equilibrium=True, operation=2, points=points, filename=filename, **kwargs))
    p0_z = _as_2d(read_field("p", timeslices=timeslice, equilibrium=True, operation=3, points=points, filename=filename, **kwargs))
    i0 = _as_2d(read_field("I", timeslices=timeslice, equilibrium=True, points=points, filename=filename, **kwargs))
    i0_r = _as_2d(read_field("I", timeslices=timeslice, equilibrium=True, operation=2, points=points, filename=filename, **kwargs))
    i0_z = _as_2d(read_field("I", timeslices=timeslice, equilibrium=True, operation=3, points=points, filename=filename, **kwargs))

    b2_num = s_bracket(psi, psi, x, z) + i0**2
    beta = r**2 * 2.0 * p0 / b2_num
    beta0 = float(np.nanmean(2.0 * p0 * r**2 / max((bzero * rzero) ** 2, np.finfo(float).tiny)))
    jphi = psi_lp - psi_r / np.maximum(r, np.finfo(float).tiny)
    tcur = np.asarray(read_scalar("ip", filename=filename, mks=True), dtype=float).reshape(-1)
    zip_current = float(tcur[0]) if tcur.size else 0.0
    print("current = ", zip_current)

    r2bp = psi_r**2 + psi_z**2
    r2bp = np.where(np.abs(r2bp) < np.finfo(float).tiny, np.nan, r2bp)

    flux_native = np.linspace(float(flux0), float(psilim_val), int(points), dtype=float)
    try:
        fc = flux_coordinates(
            psi0=psi,
            i0=i0,
            x=x,
            z=z,
            dpsi0_dx=psi_r,
            dpsi0_dz=psi_z,
            filename=filename,
            points=points,
            slice=timeslice,
            **kwargs,
        )
        flux_native = np.asarray(fc.psi, dtype=float)
        p = np.asarray(flux_average_field(p0, psi, x=x, z=z, fc=fc, filename=filename, **kwargs), dtype=float)
        pp = (p0_r * psi_r + p0_z * psi_z) / r2bp
        pprime = np.asarray(flux_average_field(pp, psi, x=x, z=z, fc=fc, filename=filename, **kwargs), dtype=float)
        i_prof = np.asarray(flux_average_field(i0, psi, x=x, z=z, fc=fc, filename=filename, **kwargs), dtype=float)
        ffp = i0 * (i0_r * psi_r + i0_z * psi_z) / r2bp
        ffprim = np.asarray(flux_average_field(ffp, psi, x=x, z=z, fc=fc, filename=filename, **kwargs), dtype=float)
        q = np.abs(np.asarray(fc.q, dtype=float))
    except Exception as exc:
        print("WARNING: flux-coordinate profiles failed; using psi-bin profiles")
        print(exc)
        p = np.full(flux_native.size, np.nan, dtype=float)
        pprime = np.full(flux_native.size, np.nan, dtype=float)
        i_prof = np.full(flux_native.size, np.nan, dtype=float)
        ffprim = np.full(flux_native.size, np.nan, dtype=float)
        q = np.zeros(flux_native.size, dtype=float)

    if not np.any(np.isfinite(p)) or not np.any(np.isfinite(i_prof)) or not np.any(np.isfinite(q)):
        print("WARNING: flux-coordinate profiles are non-finite; using psi-bin profiles")
        flux_native = np.linspace(float(flux0), float(psilim_val), int(points), dtype=float)
        p = _profile_by_psi_bins(p0, psi, flux_native)
        i_prof = _profile_by_psi_bins(i0, psi, flux_native)
        pprime = np.gradient(p, flux_native, edge_order=1)
        ffprim = i_prof * np.gradient(i_prof, flux_native, edge_order=1)
        q = np.zeros(flux_native.size, dtype=float)

    jb = (i0_z * psi_r - i0_r * psi_z - jphi * i0) / np.maximum(r**2, np.finfo(float).tiny)
    jdotb = _profile_by_psi_bins(jb, psi, flux_native)
    r2i = _profile_by_psi_bins(1.0 / np.maximum(r**2, np.finfo(float).tiny), psi, flux_native)
    betacent = float(np.asarray(field_at_point(beta, x, z, axis[0], axis[1])).reshape(-1)[0])

    c = 3.0e10
    p = p * b0**2 / (4.0 * np.pi) / 10.0
    p0 = p0 * b0**2 / (4.0 * np.pi) / 10.0
    psi_si = psi * b0 * l0**2 / 1.0e8
    flux = np.asarray(flux_native, dtype=float) * b0 * l0**2 / 1.0e8
    flux0_si = flux0 * b0 * l0**2 / 1.0e8
    psilim_si = psilim_val * b0 * l0**2 / 1.0e8
    pprime = pprime * b0**2 / (4.0 * np.pi) / 10.0
    pprime = pprime / (b0 * l0**2) * 1.0e8
    i_prof = i_prof * b0 * l0 / (1.0e4 * 100.0)
    ffprim = ffprim * (b0 * l0) ** 2 / (1.0e4 * 100.0) ** 2
    ffprim = ffprim / (b0 * l0**2) * 1.0e8
    bzero_si = bzero * b0 / 1.0e4
    x_si = x * l0 / 100.0
    z_si = z * l0 / 100.0
    rlim_si = rlim * l0 / 100.0
    zlim_si = zlim * l0 / 100.0
    rwall_si = rwall * l0 / 100.0
    zwall_si = zwall * l0 / 100.0
    axis_si = axis * l0 / 100.0
    rzero_si = rzero * l0 / 100.0
    jdotb = jdotb * b0**2 * c / (l0 * 4.0 * np.pi) / (1.0e4 * 3.0e5)
    r2i = r2i / l0**2 * 100.0**2

    nr = flux.size
    nz = z_si.size
    print("nr = ", nr)

    psimin = flux0_si
    if psimin > psilim_si:
        psi_si = -psi_si
        psilim_si = -psilim_si
        psimin = -psimin
        flux = -flux
        pprime = -pprime
        ffprim = -ffprim

    rdim = float(np.nanmax(x_si) - np.nanmin(x_si))
    zdim = float(np.nanmax(z_si) - np.nanmin(z_si))
    ccon = float(np.nanmin(x_si))
    zmid = float((np.nanmax(z_si) + np.nanmin(z_si)) / 2.0)
    rmag = float(axis_si[0])
    zmag = float(axis_si[1])
    rcentr = rmag
    bcentr = float(bzero_si * rzero_si / max(rcentr, np.finfo(float).tiny))
    beta_n = float(100.0 * (bzero_si * rzero_si / rmag) * beta0 / (zip_current / 1.0e6))
    xdum = 0.0

    print("rdim, zdim", rdim, zdim)
    print("rmag, zmag", rmag, zmag)
    print("bcentr = ", bcentr)
    print("beta0 = ", beta0)
    print("betacent = ", betacent)
    print("beta_n = ", beta_n)
    print("zip = ", zip_current)
    print("psimin, psilim = ", psimin, psilim_si)
    print("min, max (flux) = ", np.nanmin(flux), np.nanmax(flux))
    print("nr, nz, nlim, nwall")
    print(nr, nz, nlim, nwall)
    print("outputting to eqdsk format...")

    with eqfile.open("w", encoding="ascii") as handle:
        _write_geqdsk_header(handle, nr, nz)
        _write_e16(handle, [rdim, zdim, rcentr, ccon, zmid])
        _write_e16(handle, [rmag, zmag, psimin, psilim_si, bcentr])
        _write_e16(handle, [zip_current, psimin, beta0, rmag, betacent])
        _write_e16(handle, [zmag, beta_n, psilim_si, xdum, xdum])
        _write_e16(handle, i_prof)
        _write_e16(handle, p)
        _write_e16(handle, ffprim)
        _write_e16(handle, pprime)
        _write_e16(handle, psi_si.T.reshape(-1))
        _write_e16(handle, q)
        _write_i5(handle, [nlim, nwall])
        _write_e16(handle, np.column_stack((rlim_si, zlim_si)).reshape(-1))
        _write_e16(handle, np.column_stack((rwall_si, zwall_si)).reshape(-1))

    if jsolver:
        print("outputting to jsolver format...")
        _write_jsolver(
            Path(jsfile),
            flux=flux,
            p=p,
            pprime=pprime,
            i_prof=i_prof,
            jdotb=jdotb,
            r2i=r2i,
            rlim=rlim_si,
            zlim=zlim_si,
            axis=axis_si,
            bzero=bzero_si,
            rzero=rzero_si,
            ifixedb=ifixedb,
        )

    return eqfile


def _write_jsolver(
    path: Path,
    *,
    flux: np.ndarray,
    p: np.ndarray,
    pprime: np.ndarray,
    i_prof: np.ndarray,
    jdotb: np.ndarray,
    r2i: np.ndarray,
    rlim: np.ndarray,
    zlim: np.ndarray,
    axis: np.ndarray,
    bzero: float,
    rzero: float,
    ifixedb: int,
) -> None:
    npsit = max(0, flux.size - 1)
    kmax = max(0, rlim.size - 1)
    psimins = float(np.nanmin(flux))
    psilims = float(np.nanmax(flux))
    rlim_out = rlim[::-1] if ifixedb == 1 else rlim
    zlim_out = zlim[::-1] if ifixedb == 1 else zlim
    if ifixedb == 1:
        print("Reversing boundary points")
    ajpest2 = jdotb / np.maximum(i_prof * r2i, np.finfo(float).tiny)

    with path.open("w", encoding="ascii") as handle:
        handle.write("".join(f"{v:10d}" for v in [1, 0, 1, npsit, kmax]) + "\n")
        _write_e20(handle, [0.1140e-1, axis[0], axis[1], bzero * rzero, 0.3656e5, 0.5184e-2, 0.1472, 0.5458, 0.5, psimins, psilims])
        mu0 = 4.0 * np.pi * 1.0e-7
        _write_e20(handle, mu0 * p[:npsit])
        _write_e20(handle, mu0 * pprime[:npsit])
        _write_e20(handle, -mu0 * ajpest2[:npsit])
        _write_e20(handle, flux[:npsit] - flux[0])
        _write_e20(handle, rlim_out)
        _write_e20(handle, zlim_out)


def _write_e20(handle, values) -> None:
    vals = np.asarray(values, dtype=float).reshape(-1)
    for i in range(0, vals.size, 5):
        handle.write("".join(f"{v:20.12E}" for v in vals[i : i + 5]) + "\n")
