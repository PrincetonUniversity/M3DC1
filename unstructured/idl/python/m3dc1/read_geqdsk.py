from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_EXP_FIX_RE = re.compile(r"^([+-]?\d*\.\d+)([+-]\d+)$")


@dataclass
class GeqdskData:
    filename: str
    description: str
    nw: int
    nh: int
    rdim: float
    zdim: float
    rcentr: float
    rleft: float
    zmid: float
    rmaxis: float
    zmaxis: float
    simag: float
    sibry: float
    bcentr: float
    current: float
    fpol: np.ndarray
    pres: np.ndarray
    ffprim: np.ndarray
    pprime: np.ndarray
    psirz: np.ndarray
    qpsi: np.ndarray
    nbbbs: int
    limitr: int
    rbbbs: np.ndarray
    zbbbs: np.ndarray
    rlim: np.ndarray
    zlim: np.ndarray
    r: np.ndarray
    z: np.ndarray
    psin: np.ndarray


def _fortran_float(token: str) -> float:
    token = token.strip()
    try:
        return float(token)
    except ValueError:
        pass
    token = token.replace("D", "E").replace("d", "e")
    try:
        return float(token)
    except ValueError:
        pass
    # Fortran occasionally drops the "E" for 3-digit exponents, e.g. "-1.234567890-100"
    m = _EXP_FIX_RE.match(token)
    if m:
        return float(f"{m.group(1)}E{m.group(2)}")
    raise ValueError(f"Cannot parse Fortran float token: {token!r}")


def _read_block(lines: list[str], li: int, n: int) -> tuple[np.ndarray, int]:
    """Read n values written in fixed-width 16-char (5e16.9) records starting at line li."""
    nlines = -(-n // 5)  # ceil(n / 5)
    text = "".join(lines[li : li + nlines]).replace("\n", "")
    vals = [_fortran_float(text[16 * i : 16 * (i + 1)]) for i in range(n)]
    return np.asarray(vals, dtype=float), li + nlines


def read_geqdsk(filename: str | Path) -> GeqdskData:
    """Read a G-EQDSK equilibrium file (EFIT-style fixed-width format)."""
    with open(filename, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    header = lines[0].rstrip("\n")
    tokens = header.split()
    idum, nw, nh = (int(t) for t in tokens[-3:])
    description = header.strip()

    li = 1
    vals, li = _read_block(lines, li, 5)
    rdim, zdim, rcentr, rleft, zmid = vals

    vals, li = _read_block(lines, li, 5)
    rmaxis, zmaxis, simag, sibry, bcentr = vals

    vals, li = _read_block(lines, li, 5)
    current = vals[0]

    _, li = _read_block(lines, li, 5)

    fpol, li = _read_block(lines, li, nw)
    pres, li = _read_block(lines, li, nw)
    ffprim, li = _read_block(lines, li, nw)
    pprime, li = _read_block(lines, li, nw)
    psirz_flat, li = _read_block(lines, li, nw * nh)
    qpsi, li = _read_block(lines, li, nw)

    bnd_tokens = lines[li].split()
    nbbbs, limitr = int(float(bnd_tokens[0])), int(float(bnd_tokens[1]))
    li += 1

    bbbs, li = _read_block(lines, li, 2 * nbbbs)
    rbbbs, zbbbs = bbbs[0::2].copy(), bbbs[1::2].copy()

    lim, li = _read_block(lines, li, 2 * limitr)
    rlim, zlim = lim[0::2].copy(), lim[1::2].copy()

    r = rleft + rdim * np.arange(nw) / max(nw - 1, 1)
    z = zmid - zdim / 2.0 + zdim * np.arange(nh) / max(nh - 1, 1)
    psirz = psirz_flat.reshape(nh, nw)  # psirz[iz, ir]
    psin = np.linspace(0.0, 1.0, nw)

    return GeqdskData(
        filename=str(filename),
        description=description,
        nw=nw,
        nh=nh,
        rdim=rdim,
        zdim=zdim,
        rcentr=rcentr,
        rleft=rleft,
        zmid=zmid,
        rmaxis=rmaxis,
        zmaxis=zmaxis,
        simag=simag,
        sibry=sibry,
        bcentr=bcentr,
        current=current,
        fpol=fpol,
        pres=pres,
        ffprim=ffprim,
        pprime=pprime,
        psirz=psirz,
        qpsi=qpsi,
        nbbbs=nbbbs,
        limitr=limitr,
        rbbbs=rbbbs,
        zbbbs=zbbbs,
        rlim=rlim,
        zlim=zlim,
        r=r,
        z=z,
        psin=psin,
    )
