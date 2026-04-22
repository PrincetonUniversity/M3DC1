from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from .dx import dx
from .dz import dz
from .field_at_point import field_at_point
from .lcfs import lcfs
from .read_field import read_field
from .read_parameter import read_parameter


@dataclass
class FluxCoordinates:
    m: int
    n: int
    r: np.ndarray
    z: np.ndarray
    r0: float
    z0: float
    omega: np.ndarray
    psi: np.ndarray
    psi_norm: np.ndarray
    flux_pol: np.ndarray
    theta: np.ndarray
    j: np.ndarray
    q: np.ndarray
    area: np.ndarray
    dV_dchi: np.ndarray
    pest: bool
    boozer: bool
    hamada: bool
    V: np.ndarray
    flux_tor: np.ndarray
    phi_norm: np.ndarray
    rho: np.ndarray
    period: float
    itor: int
    current: np.ndarray
    dpsi_dchi: float


def _extract_path_at_level(psi2d: np.ndarray, x: np.ndarray, z: np.ndarray, level: float) -> np.ndarray | None:
    fig, ax = plt.subplots(figsize=(2, 2))
    cs = ax.contour(x, z, psi2d.T, levels=[float(level)])
    try:
        paths = []
        if hasattr(cs, "collections") and cs.collections:
            paths = [p.vertices for p in cs.collections[0].get_paths() if p.vertices.shape[0] > 4]
        elif hasattr(cs, "allsegs") and cs.allsegs and cs.allsegs[0]:
            paths = [np.asarray(seg, dtype=float) for seg in cs.allsegs[0] if np.asarray(seg).shape[0] > 4]
        if not paths:
            return None
        k = int(np.argmax([pv.shape[0] for pv in paths]))
        return np.asarray(paths[k], dtype=float)
    finally:
        plt.close(fig)


def _sample_ring_from_path(path: np.ndarray, axis: np.ndarray, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rr = path[:, 0] - float(axis[0])
    zz = path[:, 1] - float(axis[1])
    ang = np.mod(-np.arctan2(zz, rr), 2.0 * np.pi)  # clockwise
    rho = np.sqrt(rr**2 + zz**2)
    order = np.argsort(ang)
    ang = ang[order]
    rho = rho[order]
    if ang.size < 3:
        return None, None
    ang = np.r_[ang, ang[0] + 2.0 * np.pi]
    rho = np.r_[rho, rho[0]]
    rhot = np.interp(theta, ang, rho)
    r = float(axis[0]) + rhot * np.cos(-theta)
    z = float(axis[1]) + rhot * np.sin(-theta)
    return r, z


def _safe_nonzero_scalar(v: float) -> float:
    tiny = float(np.finfo(float).tiny)
    if not np.isfinite(v):
        return 1.0
    if abs(v) < tiny:
        return float(np.copysign(tiny, v if v != 0.0 else 1.0))
    return float(v)


def flux_coordinates(
    *,
    pest: bool = False,
    points: int | None = None,
    plot: bool = False,
    fast: bool = False,
    psi0=None,
    i0=None,
    x=None,
    z=None,
    fbins: int | None = None,
    tbins: int | None = None,
    boozer: bool = False,
    hamada: bool = False,
    njac: bool = False,
    itor: int | None = None,
    psin_range=None,
    r0: float | None = None,
    filename="C1.h5",
    dpsi0_dx=None,
    dpsi0_dz=None,
    slice: int = 0,
    **kwargs,
) -> FluxCoordinates:
    fn = filename[0] if isinstance(filename, (list, tuple)) else filename
    if itor is None:
        itor = int(read_parameter("itor", filename=fn, **kwargs))
    if r0 is None:
        r0 = float(read_parameter("rzero", filename=fn, **kwargs))
    if r0 == 0.0:
        r0 = 1.0
    period = float(2.0 * np.pi if itor == 1 else 2.0 * np.pi * r0)

    geo = not (bool(pest) or bool(boozer) or bool(hamada))
    if pest:
        print("Creating flux coordinates using PEST angle")
        fast = False
    elif boozer:
        print("Creating flux coordinates using BOOZER angle")
        fast = False
    elif hamada:
        print("Creating flux coordinates using HAMADA angle")
        fast = False
    else:
        print("Creating flux coordinates using GEOMETRIC angle")
    print(f"FAST MODE: {int(bool(fast))}")

    if points is None:
        points = 200 if x is None or z is None else int(np.sqrt(np.asarray(x).size * np.asarray(z).size))
    if fbins is None:
        fbins = int(points)
    if tbins is None:
        tbins = int(points)
    m = int(max(8, tbins))
    n = int(max(8, fbins))

    if psi0 is None or x is None or z is None:
        p = read_field("psi", filename=fn, timeslices=slice, points=points, equilibrium=True, return_meta=True, **kwargs)
        psi = np.asarray(p.data)
        psi0 = psi[0, :, :] if psi.ndim == 3 else psi
        x = np.asarray(p.r, dtype=float).reshape(-1)
        z = np.asarray(p.z, dtype=float).reshape(-1)
    else:
        psi_arr = np.asarray(psi0)
        psi0 = psi_arr[0, :, :] if psi_arr.ndim == 3 else psi_arr
        x = np.asarray(x, dtype=float).reshape(-1)
        z = np.asarray(z, dtype=float).reshape(-1)

    if dpsi0_dx is None or dpsi0_dz is None:
        if fast:
            dpsi0_dx = dx(psi0, x)
            dpsi0_dz = dz(psi0, z)
        else:
            dpx = np.asarray(
                read_field("psi", filename=fn, timeslices=slice, points=points, equilibrium=True, operation=2, return_meta=False, **kwargs)
            )
            dpz = np.asarray(
                read_field("psi", filename=fn, timeslices=slice, points=points, equilibrium=True, operation=3, return_meta=False, **kwargs)
            )
            dpsi0_dx = dpx[0, :, :] if dpx.ndim == 3 else dpx
            dpsi0_dz = dpz[0, :, :] if dpz.ndim == 3 else dpz
    else:
        dpsi0_dx = np.asarray(dpsi0_dx)
        dpsi0_dz = np.asarray(dpsi0_dz)
        if dpsi0_dx.ndim == 3:
            dpsi0_dx = dpsi0_dx[0, :, :]
        if dpsi0_dz.ndim == 3:
            dpsi0_dz = dpsi0_dz[0, :, :]

    if not fast and i0 is None:
        i_raw = np.asarray(read_field("I", filename=fn, timeslices=slice, points=points, equilibrium=True, return_meta=False, **kwargs))
        i0 = i_raw[0, :, :] if i_raw.ndim == 3 else i_raw
    elif i0 is None:
        i0 = np.zeros_like(psi0, dtype=float)
    else:
        i0 = np.asarray(i0)
        if i0.ndim == 3:
            i0 = i0[0, :, :]

    lc = lcfs(psi0, x, z, filename=fn, slice=slice, **kwargs)
    axis = np.asarray(lc.axis, dtype=float)
    flux0 = float(lc.flux0)
    psi_s = float(lc.psilim)
    dpsi = float(psi_s - flux0)
    tiny = float(np.finfo(float).tiny)
    if abs(dpsi) < tiny:
        dpsi = float(np.copysign(tiny, dpsi if dpsi != 0.0 else 1.0))

    if psin_range is None:
        psi_1d = (psi_s - flux0) * (np.arange(n, dtype=float) + 0.5) / n + flux0
        psin = (psi_1d - flux0) / dpsi
    else:
        pr = np.asarray(psin_range, dtype=float).reshape(-1)
        p0 = flux0 + (psi_s - flux0) * float(pr[0])
        p1 = flux0 + (psi_s - flux0) * float(pr[1])
        psi_1d = (p1 - p0) * (np.arange(n, dtype=float) + 0.5) / n + p0
        psin = (psi_1d - flux0) / dpsi
    dpsi_dpsin = dpsi

    theta = 2.0 * np.pi * np.arange(m, dtype=float) / m
    rpath = np.zeros((m, n), dtype=float)
    zpath = np.zeros((m, n), dtype=float)
    jac = np.zeros((m, n), dtype=float)
    omega = np.zeros((m, n), dtype=float)
    theta_sfl = np.zeros((m, n), dtype=float)
    fmetric = np.zeros(n, dtype=float)
    q = np.zeros(n, dtype=float)
    dV = np.zeros(n, dtype=float)
    V = np.zeros(n, dtype=float)
    phi = np.zeros(n, dtype=float)
    area = np.zeros(n, dtype=float)
    current = np.zeros(n, dtype=float)

    fac = -1.0 if psi_s < flux0 else 1.0
    gradpsi = np.sqrt(np.maximum(np.asarray(dpsi0_dx) ** 2 + np.asarray(dpsi0_dz) ** 2, np.finfo(float).tiny))

    prev_rho = 0.05 * max(float(np.max(x) - np.min(x)), float(np.max(z) - np.min(z)))
    for j in range(n):
        path = _extract_path_at_level(psi0, x, z, float(psi_1d[j]))
        if path is not None:
            rr, zz = _sample_ring_from_path(path, axis, theta)
            if rr is not None:
                rpath[:, j] = rr
                zpath[:, j] = zz
                prev_rho = float(np.nanmean(np.sqrt((rr - axis[0]) ** 2 + (zz - axis[1]) ** 2)))
                continue
        # fallback ring if contour extraction fails
        rho = prev_rho * (0.5 + (j + 1) / n)
        rpath[:, j] = axis[0] + rho * np.cos(-theta)
        zpath[:, j] = axis[1] + rho * np.sin(-theta)

    rp = np.asarray(rpath if itor == 1 else np.ones_like(rpath), dtype=float)
    for j in range(n):
        rx = np.r_[rpath[-1, j], rpath[:, j], rpath[0, j]]
        zx = np.r_[zpath[-1, j], zpath[:, j], zpath[0, j]]
        dr = np.gradient(rx)[1:-1]
        dzv = np.gradient(zx)[1:-1]
        dl = np.sqrt(np.maximum(dr**2 + dzv**2, np.finfo(float).tiny))

        br = -np.asarray(field_at_point(dpsi0_dz, x, z, rpath[:, j], zpath[:, j]), dtype=float) / np.maximum(rp[:, j], np.finfo(float).tiny)
        bzv = np.asarray(field_at_point(dpsi0_dx, x, z, rpath[:, j], zpath[:, j]), dtype=float) / np.maximum(rp[:, j], np.finfo(float).tiny)
        bp = np.asarray(field_at_point(gradpsi, x, z, rpath[:, j], zpath[:, j]), dtype=float) / np.maximum(rp[:, j], np.finfo(float).tiny)

        current[j] = float(np.sum(br * dr + bzv * dzv))
        area[j] = float(period * np.sum(dl * rp[:, j]))
        dV[j] = float(period * np.sum(fac * dpsi_dpsin * dl / np.maximum(bp, np.finfo(float).tiny)))
        if j == 0:
            V[j] = dV[j] * psin[j] / 2.0
        else:
            V[j] = V[j - 1] + (dV[j] + dV[j - 1]) * (psin[j] - psin[j - 1]) / 2.0

        if not fast:
            ix = np.asarray(field_at_point(i0, x, z, rpath[:, j], zpath[:, j]), dtype=float)
            if not geo:
                dthetadl = np.zeros(m, dtype=float)
                for i in range(m):
                    bpi = max(float(bp[i]), np.finfo(float).tiny)
                    rpi = max(float(rp[i, j]), np.finfo(float).tiny)
                    ixi = float(ix[i])
                    if pest:
                        dthetadl[i] = -fac * ixi / (rpi**2 * bpi)
                    elif boozer:
                        dthetadl[i] = -fac * (bpi**2 + (ixi / rpi) ** 2) / bpi
                    elif hamada:
                        dthetadl[i] = -fac / bpi
                for i in range(1, m):
                    theta_sfl[i, j] = theta_sfl[i - 1, j] + 0.5 * (dl[i] * dthetadl[i] + dl[i - 1] * dthetadl[i - 1])
                fj = float(np.sum(dthetadl * dl) / (2.0 * np.pi))
                if abs(fj) < np.finfo(float).tiny:
                    fj = 1.0
                fmetric[j] = fj
                theta_sfl[:, j] = theta_sfl[:, j] / fj

            fjr2 = -fac * ix / np.maximum(rp[:, j] ** 2 * bp, np.finfo(float).tiny)
            q[j] = float(np.sum(fjr2 * dl) / max(period, np.finfo(float).tiny))
            if j == 0:
                phi[j] = -period * q[j] * (psi_1d[j] - flux0)
            else:
                phi[j] = phi[j - 1] - period * (q[j] + q[j - 1]) * (psi_1d[j] - psi_1d[j - 1]) / 2.0

    # For non-geometric coordinates, remap rings to equal spacing in transformed angle.
    if not geo and not fast:
        for j in range(n):
            th = np.asarray(theta_sfl[:, j], dtype=float)
            rr = np.asarray(rpath[:, j], dtype=float)
            zz = np.asarray(zpath[:, j], dtype=float)
            order = np.argsort(th)
            th = th[order]
            rr = rr[order]
            zz = zz[order]
            th_ext = np.r_[th, th[0] + 2.0 * np.pi]
            rr_ext = np.r_[rr, rr[0]]
            zz_ext = np.r_[zz, zz[0]]
            rpath[:, j] = np.interp(theta, th_ext, rr_ext)
            zpath[:, j] = np.interp(theta, th_ext, zz_ext)

        rp = np.asarray(rpath if itor == 1 else np.ones_like(rpath), dtype=float)
        if pest:
            for j in range(n):
                qj = float(q[j])
                for i in range(m):
                    ival = float(field_at_point(i0, x, z, float(rpath[i, j]), float(zpath[i, j])))
                    jac[i, j] = -dpsi_dpsin * (rp[i, j] ** 2) * qj / _safe_nonzero_scalar(ival)
        elif boozer:
            b2r2 = np.asarray(i0) ** 2 + np.asarray(dpsi0_dx) ** 2 + np.asarray(dpsi0_dz) ** 2
            for j in range(n):
                fj = float(fmetric[j]) if abs(float(fmetric[j])) > np.finfo(float).tiny else 1.0
                for i in range(m):
                    bval = float(field_at_point(b2r2, x, z, float(rpath[i, j]), float(zpath[i, j])))
                    jac[i, j] = -dpsi_dpsin * fj * (rp[i, j] ** 2) / _safe_nonzero_scalar(bval)
        elif hamada:
            for j in range(n):
                fj = float(fmetric[j]) if abs(float(fmetric[j])) > np.finfo(float).tiny else 1.0
                jac[:, j] = -dpsi_dpsin * fj

    dr_dpsi = np.gradient(rpath, psin, axis=1, edge_order=1)
    dz_dpsi = np.gradient(zpath, psin, axis=1, edge_order=1)
    dr_dtheta = np.gradient(rpath, theta, axis=0, edge_order=1)
    dz_dtheta = np.gradient(zpath, theta, axis=0, edge_order=1)
    jac_test = -(dr_dpsi * dz_dtheta - dr_dtheta * dz_dpsi)
    if itor == 1:
        jac_test = jac_test * rpath
    if geo or njac:
        jac[:, :] = jac_test

    # omega = zeta - phi for non-geometric coordinates.
    if not geo and not fast:
        for j in range(n):
            qj = float(q[j])
            jac_pest_old = 1.0
            for i in range(m):
                ival = float(field_at_point(i0, x, z, float(rpath[i, j]), float(zpath[i, j])))
                jac_pest = -dpsi_dpsin * (rp[i, j] ** 2) * qj / _safe_nonzero_scalar(ival)
                if i == 0:
                    omega[i, j] = 0.0
                else:
                    dtheta = float(theta[i] - theta[i - 1])
                    t0 = 1.0 - jac[i, j] / _safe_nonzero_scalar(float(jac_pest))
                    t1 = 1.0 - jac[i - 1, j] / _safe_nonzero_scalar(float(jac_pest_old))
                    omega[i, j] = omega[i - 1, j] + dtheta * (t0 + t1) * 0.5
                jac_pest_old = jac_pest
            omega[:, j] = omega[:, j] * qj
        jac[:, :] = np.nan_to_num(jac, nan=0.0, posinf=0.0, neginf=0.0)
        omega[:, :] = np.nan_to_num(omega, nan=0.0, posinf=0.0, neginf=0.0)

    if not fast:
        den = float(phi[-1])
        if abs(den) < np.finfo(float).tiny:
            phi_norm = np.zeros_like(phi)
            rho = np.zeros_like(phi)
        else:
            phi_norm = phi / den
            rho = np.sqrt(phi_norm)
    else:
        phi_norm = np.zeros_like(phi)
        rho = np.zeros_like(phi)

    if plot:
        plt.figure(figsize=(6, 6))
        for j in range(0, n, max(1, n // 20)):
            plt.plot(np.r_[rpath[:, j], rpath[0, j]], np.r_[zpath[:, j], zpath[0, j]], color="k", linewidth=0.6)
        for i in range(0, m, max(1, m // 20)):
            plt.plot(rpath[i, :], zpath[i, :], color="tab:blue", linewidth=0.6)

    return FluxCoordinates(
        m=m,
        n=n,
        r=rpath,
        z=zpath,
        r0=float(axis[0]),
        z0=float(axis[1]),
        omega=omega,
        psi=psi_1d,
        psi_norm=psin,
        flux_pol=-period * (psi_1d - psi_1d[0]),
        theta=theta,
        j=jac,
        q=q,
        area=area,
        dV_dchi=dV,
        pest=bool(pest),
        boozer=bool(boozer),
        hamada=bool(hamada),
        V=V,
        flux_tor=phi,
        phi_norm=phi_norm,
        rho=rho,
        period=period,
        itor=int(itor),
        current=current,
        dpsi_dchi=dpsi_dpsin,
    )
