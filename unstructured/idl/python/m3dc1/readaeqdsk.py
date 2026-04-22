from __future__ import annotations

from dataclasses import dataclass
import re
from pathlib import Path
import sys

import numpy as np


_FLOAT_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?")


@dataclass
class AeqdskData:
    filename: str
    uday: str
    mfvers: tuple[str, str]
    ishot: int
    ktime1: int
    mastformat: bool
    values: dict[str, np.ndarray]
    rseps: np.ndarray
    jflag: np.ndarray
    jerror: np.ndarray
    limloc: np.ndarray
    nsilop0: int
    magpri0: int
    nfcoil0: int
    nesum0: int
    mco2v: int
    mco2r: int
    rco2r: np.ndarray
    rco2v: np.ndarray
    dco2r: np.ndarray
    dco2v: np.ndarray
    csilop: np.ndarray
    cmpr2: np.ndarray
    ccbrsp: np.ndarray
    eccurt: np.ndarray
    machine: str | None
    coil_currents_ka: np.ndarray


class _LineReader:
    def __init__(self, filename: str | Path):
        self.filename = str(filename)
        with open(self.filename, "r", encoding="utf-8", errors="replace") as fh:
            self._lines = fh.readlines()
        self._idx = 0

    def read_line(self) -> str:
        if self._idx >= len(self._lines):
            raise EOFError
        line = self._lines[self._idx].rstrip("\n")
        self._idx += 1
        return line

    def read_float_values(self, count: int) -> list[float]:
        vals: list[float] = []
        while len(vals) < count:
            line = self.read_line()
            vals.extend(_parse_float_tokens(line))
        return vals[:count]

    def read_int_values(self, count: int) -> list[int]:
        vals: list[int] = []
        while len(vals) < count:
            line = self.read_line()
            vals.extend(_parse_int_tokens(line))
        return vals[:count]


def _parse_float_tokens(line: str) -> list[float]:
    out: list[float] = []
    for tok in _FLOAT_RE.findall(line):
        out.append(float(tok.replace("D", "E").replace("d", "e")))
    return out


def _parse_int_tokens(line: str) -> list[int]:
    vals = _parse_float_tokens(line)
    return [int(round(v)) for v in vals]


def _slice_int(line: str, start: int, end: int, default: int = 0) -> int:
    txt = line[start:end].strip()
    return int(txt) if txt else default


def _slice_float(line: str, start: int, end: int, default: float = 0.0) -> float:
    txt = line[start:end].strip()
    return float(txt.replace("D", "E").replace("d", "e")) if txt else default


def _init_scalar_arrays(ktime1: int) -> dict[str, np.ndarray]:
    scalar_names = [
        "time",
        "eout",
        "rout",
        "zout",
        "doutu",
        "doutl",
        "aout",
        "vout",
        "betat",
        "otop",
        "betap",
        "ali",
        "oleft",
        "oright",
        "qsta",
        "rcurrt",
        "zcurrt",
        "qout",
        "olefs",
        "orighs",
        "otops",
        "sibdry",
        "areao",
        "wplasm",
        "elongm",
        "qqmagx",
        "terror",
        "rmagx",
        "zmagx",
        "obott",
        "obots",
        "alpha",
        "rttt",
        "dbpli",
        "delbp",
        "oring",
        "sepexp",
        "shearb",
        "xtch",
        "ytch",
        "qpsib",
        "vertn",
        "aaq1",
        "aaq2",
        "aaq3",
        "btaxp",
        "btaxv",
        "simagx",
        "seplim",
        "wbpol",
        "taumhd",
        "betapd",
        "betatd",
        "alid",
        "wplasmd",
        "taudia",
        "wbpold",
        "qmerci",
        "slantu",
        "slantl",
        "zeff",
        "zeffr",
        "tave",
        "rvsin",
        "zvsin",
        "rvsout",
        "zvsout",
        "wpdot",
        "wbdot",
        "vsurfa",
        "cjor95",
        "pp95",
        "ssep",
        "yyy2",
        "xnnc",
        "wtherm",
        "wfbeam",
        "taujd3",
        "tauthn",
        "qsiwant",
        "cjorsw",
        "cjor0",
        "ssiwant",
        "ssi95",
        "peak",
        "dminux",
        "dminlx",
        "dolubat",
        "dolubafm",
        "diludom",
        "diludomm",
        "ratsol",
        "rvsiu",
        "zvsiu",
        "rvsid",
        "zvsid",
        "rvsou",
        "zvsou",
        "rvsod",
        "zvsod",
        "condno",
        "psin32",
        "psin21",
        "rq32in",
        "rq21top",
        "chilibt",
        "tsaisq",
        "bcentr",
        "cpasma",
        "pasmat",
        "bpolav",
        "s1",
        "s2",
        "s3",
        "cdflux",
        "psiref",
        "xndnt",
        "vloopt",
        "pbinj",
        "zuperts",
        "cjor99",
        "psurfa",
        "dolubaf",
        "cj1ave",
        "rmidin",
        "rmidout",
    ]
    return {name: np.zeros(ktime1, dtype=float) for name in scalar_names}


def _append(out: list[float], value: float) -> None:
    out.append(float(value))


def _print_f1100(value: float) -> None:
    print(f"{value:12.4f}")


def _print_f0(*args) -> None:
    print(*args, file=sys.stderr)


def _select_machine_currents(ccbrsp: np.ndarray, eccurt: np.ndarray, nfcoil0: int, nesum0: int) -> tuple[str | None, np.ndarray]:
    out: list[float] = []
    machine: str | None = None

    if nfcoil0 == 18:
        machine = "DIII-D"
        _print_f0(" Assuming DIII-D")
        if np.max(np.abs(ccbrsp[:18, 0])) < 20000.0:
            _print_f0(" File data is in Amps")
            # File data is in Amps
            factors = [58.0, 58.0, 58.0, 58.0, 58.0, 58.0, 58.0, 55.0, 58.0, 55.0, 58.0, 55.0, 58.0, 55.0, 58.0, 55.0, 58.0, 55.0]
        else:
            _print_f0(" File data is in Amp-turns")
            # File data is in Amp-turns
            factors = [1.0] * 18
        order = [4, 3, 2, 1, 10, 11, 12, 13, 5, 14, 8, 17, 9, 18, 7, 16, 6, 15]
        for idx, factor in zip(order, factors):
            val = ccbrsp[idx - 1, 0] * factor / 1000.0
            _append(out, val)
            _print_f1100(val)

    elif nfcoil0 == 12:
        machine = "EAST"
        _print_f0(" Assuming EAST")
        # PF1
        for val in [ccbrsp[0, 0] / 1000.0, ccbrsp[6, 0] / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        # PF3
        for val in [ccbrsp[1, 0] / 1000.0, ccbrsp[7, 0] / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        # PF5
        for val in [ccbrsp[2, 0] / 1000.0, ccbrsp[8, 0] / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        # PF7 / PF9
        for val in [ccbrsp[3, 0] * (44.0 / 248.0) / 1000.0, ccbrsp[3, 0] * (204.0 / 248.0) / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        # PF8 / PF10
        for val in [ccbrsp[9, 0] * (44.0 / 248.0) / 1000.0, ccbrsp[9, 0] * (204.0 / 248.0) / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        # PF11 / PF12
        for val in [ccbrsp[4, 0] / 1000.0, ccbrsp[10, 0] / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        # PF13 / PF14
        for val in [ccbrsp[5, 0] / 1000.0, ccbrsp[11, 0] / 1000.0]:
            _append(out, val)
            _print_f1100(val)

    elif nfcoil0 == 52:
        machine = "NSTX"
        _print_f0(" Assuming NSTX")
        _print_f0(f" nesum0 = {nesum0}")
        # OH coils
        for factor in [122.625, 121.000, 119.875, 118.000, 118.000, 119.875, 121.000, 122.625]:
            val = eccurt[0, 0] * factor / 1000.0
            _append(out, val)
            _print_f1100(val)
        # PF and vessel coils
        entries = [
            (1, 20.0), (2, 14.0), (2, 14.0), (3, 15.0), (3, 15.0), (4, 9.0), (4, 8.0),
            (5, 12.0), (5, 12.0), (6, 12.0), (6, 12.0), (7, 9.0), (7, 8.0), (8, 15.0),
            (8, 15.0), (9, 14.0), (9, 14.0), (10, 20.0), (11, 32.0), (12, 1.0), (13, 1.0),
            (15, 1.0), (14, 1.0), (16, 48.0), (17, 48.0), (18, 1.0), (19, 1.0), (20, 1.0),
            (21, 0.5), (21, 0.5), (22, 1.0), (23, 0.5), (23, 0.5), (24, 0.3333), (24, 0.3333),
            (24, 0.3333), (25, 0.1250), (25, 0.1250), (25, 0.1250), (25, 0.1250), (25, 0.1250),
            (25, 0.1250), (25, 0.1250), (25, 0.1250), (26, 0.2), (26, 0.2), (26, 0.2),
            (26, 0.2), (26, 0.2), (27, 0.5), (27, 0.5), (28, 0.3333), (28, 0.3333), (28, 0.3333),
            (29, 1.0), (30, 1.0), (31, 1.0), (32, 1.0), (33, 0.3333), (33, 0.3333), (33, 0.3333),
            (34, 0.5), (34, 0.5), (35, 0.2), (35, 0.2), (35, 0.2), (35, 0.2), (35, 0.2),
            (36, 0.2), (36, 0.2), (36, 0.2), (36, 0.2), (36, 0.2), (37, 0.5), (37, 0.5),
            (38, 1.0), (39, 0.5), (39, 0.5), (40, 1.0), (41, 1.0), (42, 1.0), (43, 0.5),
            (43, 0.5), (44, 0.5), (44, 0.5), (45, 1.0), (46, 1.0), (47, 1.0), (48, 1.0),
            (49, 1.0), (50, 1.0), (51, 1.0), (52, 1.0),
        ]
        for idx, factor in entries:
            val = ccbrsp[idx - 1, 0] * factor / 1000.0
            _append(out, val)
            _print_f1100(val)

    elif nfcoil0 == 54:
        machine = "NSTX-U"
        _print_f0(" Assuming NSTX-U")
        # OH coils
        for _ in range(8):
            _append(out, 0.0)
            _print_f1100(0.0)
        entries = [
            (1, 64.0), (2, 32.0), (3, 20.0), (4, 14.0), (4, 14.0), (5, 15.0), (5, 15.0),
            (6, 9.0), (6, 8.0), (7, 12.0), (7, 12.0), (8, 12.0), (8, 12.0), (9, 9.0), (9, 8.0),
            (10, 15.0), (10, 15.0), (11, 14.0), (11, 14.0), (12, 20.0), (13, 32.0), (14, 64.0),
            (15, 0.5), (15, 0.5), (16, 0.5), (16, 0.5), (17, 1.0), (18, 1.0),
            (19, 0.3333), (19, 0.3333), (19, 0.3333),
        ]
        for idx, factor in entries:
            val = ccbrsp[idx - 1, 0] * factor / 1000.0
            _append(out, val)
            _print_f1100(val)
        for _ in range(16):
            val = ccbrsp[19, 0] * 0.0625 / 1000.0
            _append(out, val)
            _print_f1100(val)
        for _ in range(5):
            val = ccbrsp[20, 0] * 0.2000 / 1000.0
            _append(out, val)
            _print_f1100(val)
        for _ in range(4):
            val = ccbrsp[21, 0] * 0.2500 / 1000.0
            _append(out, val)
            _print_f1100(val)
        for _ in range(24):
            val = ccbrsp[22, 0] * 0.041667 / 1000.0
            _append(out, val)
            _print_f1100(val)
        for val in [ccbrsp[23, 0] * 0.5 / 1000.0, ccbrsp[23, 0] * 0.5 / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        for _ in range(5):
            val = ccbrsp[24, 0] * 0.2 / 1000.0
            _append(out, val)
            _print_f1100(val)
        for val in [ccbrsp[25, 0] * 0.5 / 1000.0, ccbrsp[25, 0] * 0.5 / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        for val in [ccbrsp[26, 0] * 0.3333 / 1000.0, ccbrsp[26, 0] * 0.3333 / 1000.0, ccbrsp[26, 0] * 0.3333 / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        for val in [ccbrsp[27, 0] / 1000.0, ccbrsp[28, 0] / 1000.0, ccbrsp[29, 0] / 1000.0, ccbrsp[30, 0] / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        for val in [ccbrsp[31, 0] * 0.3333 / 1000.0, ccbrsp[31, 0] * 0.3333 / 1000.0, ccbrsp[31, 0] * 0.3333 / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        for val in [ccbrsp[32, 0] * 0.5 / 1000.0, ccbrsp[32, 0] * 0.5 / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        for _ in range(5):
            val = ccbrsp[33, 0] * 0.2 / 1000.0
            _append(out, val)
            _print_f1100(val)
        for val in [ccbrsp[34, 0] * 0.5 / 1000.0, ccbrsp[34, 0] * 0.5 / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        for _ in range(24):
            val = ccbrsp[35, 0] * 0.041667 / 1000.0
            _append(out, val)
            _print_f1100(val)
        for _ in range(16):
            val = ccbrsp[36, 0] * 0.0625 / 1000.0
            _append(out, val)
            _print_f1100(val)
        for _ in range(4):
            val = ccbrsp[37, 0] * 0.25 / 1000.0
            _append(out, val)
            _print_f1100(val)
        for _ in range(5):
            val = ccbrsp[38, 0] * 0.2 / 1000.0
            _append(out, val)
            _print_f1100(val)
        for val in [ccbrsp[39, 0] * 0.3333 / 1000.0, ccbrsp[39, 0] * 0.3333 / 1000.0, ccbrsp[39, 0] * 0.3333 / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        for val in [ccbrsp[40, 0] / 1000.0, ccbrsp[41, 0] / 1000.0, ccbrsp[42, 0] * 0.5 / 1000.0, ccbrsp[42, 0] * 0.5 / 1000.0, ccbrsp[43, 0] * 0.5 / 1000.0, ccbrsp[43, 0] * 0.5 / 1000.0]:
            _append(out, val)
            _print_f1100(val)
        # Divertor / passive plates
        for idx, factor in [
            (45, 0.5), (45, 0.5), (46, 0.5), (46, 0.5), (47, 1.0), (48, 1.0), (49, 1.0), (50, 1.0),
            (51, 1.0), (52, 1.0), (53, 1.0), (54, 1.0),
        ]:
            val = ccbrsp[idx - 1, 0] * factor / 1000.0
            _append(out, val)
            _print_f1100(val)

    elif nfcoil0 == 31:
        machine = "KSTAR"
        _print_f0(" Assuming KSTAR ")
        # Retain the active branch from the Fortran source. The older commented
        # alternatives are preserved below as comments in the Python source.
        entries = [
            # PF coils
            (27, 1.0), (28, 1.0), (19, 1.0), (20, 1.0), (21, 1.0), (22, 1.0), (29, 1.0), (29, 1.0),
            (26, 1.0), (25, 1.0), (24, 1.0), (23, 1.0), (28, 1.0), (27, 1.0),
            # Outer vessel coils
            (1, 1.0), (1, 1.0), (2, 1.0), (2, 1.0), (2, 1.0), (3, 1.0), (3, 1.0), (4, 1.0), (4, 1.0),
            (5, 1.0), (5, 1.0), (5, 1.0), (5, 1.0), (6, 1.0), (6, 1.0), (7, 1.0), (7, 1.0), (8, 1.0),
            (8, 1.0), (8, 1.0), (8, 1.0), (9, 1.0), (9, 1.0), (10, 1.0), (10, 1.0), (11, 1.0), (11, 1.0),
            (11, 1.0), (12, 1.0), (12, 1.0),
            # Inner vessel coils
            (1, 1.0), (1, 1.0), (2, 1.0), (2, 1.0), (2, 1.0), (3, 1.0), (3, 1.0), (4, 1.0), (4, 1.0),
            (5, 1.0), (5, 1.0), (5, 1.0), (5, 1.0), (6, 1.0), (6, 1.0), (7, 1.0), (7, 1.0), (8, 1.0),
            (8, 1.0), (8, 1.0), (8, 1.0), (9, 1.0), (9, 1.0), (10, 1.0), (10, 1.0), (11, 1.0), (11, 1.0),
            (11, 1.0), (12, 1.0), (12, 1.0),
        ]
        for idx, factor in entries:
            val = ccbrsp[idx - 1, 0] * factor / 1000.0
            _append(out, val)
            _print_f1100(val)

    elif nfcoil0 == 101:
        machine = "MAST"
        _print_f0(" Assuming MAST")
        # Each turn in OH receives only ...
        entries = [
            (1, 1.0), (2, 1.0), (3, 1.0), (4, 1.0), (5, 1.0),
            (6, 1.0 / 3.0), (6, 1.0 / 3.0), (6, 1.0 / 3.0),
            (7, 1.0 / 3.0), (7, 1.0 / 3.0), (7, 1.0 / 3.0),
            (8, 1.0 / 3.0), (8, 1.0 / 3.0), (8, 1.0 / 3.0),
            (9, 1.0 / 3.0), (9, 1.0 / 3.0), (9, 1.0 / 3.0),
            (10, 1.0 / 3.0), (10, 1.0 / 3.0), (10, 1.0 / 3.0),
            (11, 1.0 / 3.0), (11, 1.0 / 3.0), (11, 1.0 / 3.0),
            (12, 1.0 / 2.0), (12, 1.0 / 2.0), (13, 1.0 / 2.0), (13, 1.0 / 2.0),
            (14, 1.0), (14, 1.0), (14, 1.0), (14, 1.0),
            (15, 1.0), (15, 1.0), (15, 1.0), (15, 1.0),
            (16, 1.0), (16, 1.0), (16, 1.0), (16, 1.0), (16, 1.0), (16, 1.0), (16, 1.0), (16, 1.0),
            (17, 1.0), (17, 1.0), (17, 1.0), (17, 1.0), (17, 1.0), (17, 1.0), (17, 1.0), (17, 1.0),
            (18, 1.0), (18, 1.0), (18, 1.0), (18, 1.0),
            (19, 1.0), (19, 1.0), (19, 1.0), (19, 1.0),
            (20, 1.0), (20, 1.0), (20, 1.0), (20, 1.0),
            (21, 1.0), (21, 1.0), (21, 1.0), (21, 1.0),
            (22, 1.0), (22, 1.0), (22, 1.0), (22, 1.0),
            (23, 1.0), (24, 1.0), (25, 1.0), (26, 1.0), (27, 1.0), (28, 1.0), (29, 1.0), (30, 1.0),
            (31, 1.0), (32, 1.0), (33, 1.0), (34, 1.0), (35, 1.0), (36, 1.0), (37, 1.0), (38, 1.0),
            (39, 1.0), (40, 1.0), (41, 1.0), (42, 1.0), (43, 1.0), (44, 1.0), (45, 1.0), (46, 1.0),
            (47, 1.0), (48, 1.0), (49, 1.0), (50, 1.0), (51, 1.0), (52, 1.0), (53, 1.0), (54, 1.0),
            (55, 1.0), (56, 1.0), (57, 1.0), (58, 1.0), (59, 1.0), (60, 1.0), (61, 1.0), (62, 1.0),
            (63, 1.0), (64, 1.0), (68, 1.0), (69, 1.0), (70, 1.0), (71, 1.0), (72, 1.0), (73, 1.0),
            (74, 1.0), (75, 1.0), (77, 1.0), (78, 1.0), (79, 1.0), (80, 1.0), (81, 1.0), (82, 1.0),
            (83, 1.0), (84, 1.0), (85, 1.0), (86, 1.0), (87, 1.0), (88, 1.0), (89, 1.0), (90, 1.0),
            (91, 1.0), (92, 1.0), (93, 1.0), (94, 1.0), (95, 1.0), (96, 1.0), (97, 1.0), (98, 1.0),
            (99, 1.0), (100, 1.0), (101, 1.0), (101, 1.0), (101, 1.0), (101, 1.0),
        ]
        for idx, factor in entries:
            val = ccbrsp[idx - 1, 0] * factor / 1000.0
            _append(out, val)
            _print_f1100(val)

    return machine, np.asarray(out, dtype=float)


def load_eqdsk_a(filename: str | Path) -> AeqdskData:
    # Reading EQDSK a-file
    _print_f0(f" Reading EQDSK a-file: {str(filename).strip()}")
    reader = _LineReader(filename)

    line = reader.read_line()
    uday = line[1:11]
    mfvers = (line[11:16].strip(), line[16:21].strip())

    line = reader.read_line()
    ishot = _slice_int(line, 1, 7, 0)
    ktime1 = _slice_int(line, 18, 23, 0)

    # MAST a-files do not contain the following line in which the times are listed.
    # Use content of first line to identify MAST a-files.
    mastformat = "EFIT+" in mfvers[0]
    values = _init_scalar_arrays(ktime1)
    if not mastformat:
        values["time"][0] = reader.read_float_values(1)[0]
    _print_f0(uday.strip(), mfvers[0], mfvers[1])
    _print_f0(f"  SHOT, time = {ishot} {ktime1}")

    jflag = np.zeros(ktime1, dtype=int)
    jerror = np.zeros(ktime1, dtype=int)
    limloc = np.full(ktime1, "", dtype="U4")
    rseps = np.zeros((2, ktime1), dtype=float)

    nsilop0 = 0
    magpri0 = 0
    nfcoil0 = 0
    nesum0 = 0
    mco2v = 0
    mco2r = 0

    rco2r = np.zeros((0, ktime1), dtype=float)
    rco2v = np.zeros((0, ktime1), dtype=float)
    dco2r = np.zeros((ktime1, 0), dtype=float)
    dco2v = np.zeros((ktime1, 0), dtype=float)
    csilop = np.zeros((0, ktime1), dtype=float)
    cmpr2 = np.zeros((0, ktime1), dtype=float)
    ccbrsp = np.zeros((0, ktime1), dtype=float)
    eccurt = np.zeros((0, 0), dtype=float)

    pre_mco_groups = [
        ("tsaisq", "rcencm", "bcentr", "pasmat"),
        ("cpasma", "rout", "zout", "aout"),
        ("eout", "doutu", "doutl", "vout"),
        ("rcurrt", "zcurrt", "qsta", "betat"),
        ("betap", "ali", "oleft", "oright"),
        ("otop", "obott", "qpsib", "vertn"),
    ]
    post_mco_groups = [
        ("shearb", "bpolav", "s1", "s2"),
        ("s3", "qout", "olefs", "orighs"),
        ("otops", "sibdry", "areao", "wplasm"),
        ("terror", "elongm", "qqmagx", "cdflux"),
        ("alpha", "rttt", "psiref", "xndnt"),
        ("rseps1", "zseps1", "rseps2", "zseps2"),
        ("sepexp", "obots", "btaxp", "btaxv"),
        ("aaq1", "aaq2", "aaq3", "seplim"),
        ("rmagx", "zmagx", "simagx", "taumhd"),
        ("betapd", "betatd", "wplasmd", "fluxx"),
        ("vloopt", "taudia", "qmerci", "tavem"),
    ]

    for jj in range(ktime1):
        # MAST a-files do not contain nlold,nlnew and thus format 1060 cannot be used
        lineaftertime = reader.read_line()
        lineaftertime = f"{lineaftertime.rstrip()}     0    0"

        values["time"][jj] = _slice_float(lineaftertime, 1, 8, 0.0)
        jflag[jj] = _slice_int(lineaftertime, 18, 23, 0)
        limloc[jj] = lineaftertime[40:43].strip()
        mco2v = _slice_int(lineaftertime, 44, 47, 0)
        mco2r = _slice_int(lineaftertime, 48, 51, 0)

        if mco2v > 3:
            _print_f0(f" Warning: mco2v > nco2v {mco2v}")
            mco2v = 3
        if mco2r > 2:
            _print_f0(f" Warning: mco2r > nco2r {mco2r}")
            mco2r = 2

        for names in pre_mco_groups:
            vals = reader.read_float_values(4)
            for name, value in zip(names, vals):
                if name in {"rcencm", "tavem", "fluxx"}:
                    continue
                values[name][jj] = value

        if mco2v > 0:
            if rco2v.shape[0] == 0:
                rco2v = np.zeros((mco2v, ktime1), dtype=float)
                dco2v = np.zeros((ktime1, mco2v), dtype=float)
            rco2v[:mco2v, jj] = reader.read_float_values(mco2v)
            dco2v[jj, :mco2v] = reader.read_float_values(mco2v)

        if mco2r > 0:
            if rco2r.shape[0] == 0:
                rco2r = np.zeros((mco2r, ktime1), dtype=float)
                dco2r = np.zeros((ktime1, mco2r), dtype=float)
            rco2r[:mco2r, jj] = reader.read_float_values(mco2r)
            dco2r[jj, :mco2r] = reader.read_float_values(mco2r)

        for names in post_mco_groups:
            vals = reader.read_float_values(4)
            if names == ("rseps1", "zseps1", "rseps2", "zseps2"):
                rseps[0, jj] = vals[0]
                values.setdefault("zseps", np.zeros((2, ktime1), dtype=float))
                values["zseps"][0, jj] = vals[1]
                rseps[1, jj] = vals[2]
                values["zseps"][1, jj] = vals[3]
            else:
                for name, value in zip(names, vals):
                    if name in {"rcencm", "tavem", "fluxx"}:
                        continue
                    values[name][jj] = value

        # In MAST a-files the following line uses e16.9 float formatting instead of i5
        ints = reader.read_int_values(4)
        nsilop0, magpri0, nfcoil0, nesum0 = ints[:4]
        _print_f0(f" nfcoil0 = {nfcoil0}")

        if nsilop0 > 44:
            _print_f0(f" Warning: nsilop0 > nsilop {nsilop0}")
        if magpri0 > 76:
            _print_f0(f" Warning: magpri0 > magpri {magpri0}")
        if nfcoil0 > 18:
            _print_f0(f" Warning: nfcoil0 > nfcoil {nfcoil0}")
        if nesum0 > 6:
            _print_f0(f" Warning: nesum0 > nesum {nesum0}")

        if csilop.shape[0] == 0:
            csilop = np.zeros((nsilop0, ktime1), dtype=float)
            cmpr2 = np.zeros((magpri0, ktime1), dtype=float)
            ccbrsp = np.zeros((nfcoil0, ktime1), dtype=float)
            eccurt = np.zeros((ktime1, nesum0), dtype=float)

        vals = reader.read_float_values(nsilop0 + magpri0)
        csilop[:, jj] = vals[:nsilop0]
        cmpr2[:, jj] = vals[nsilop0 : nsilop0 + magpri0]

        ccbrsp[:, jj] = reader.read_float_values(nfcoil0)

        # Do not read the following in case of MAST
        if not mastformat:
            if nesum0 > 0:
                eccurt[jj, :] = reader.read_float_values(nesum0)

            extra_groups = [
                ("pbinj", "rvsin", "zvsin", "rvsout"),
                ("zvsout", "vsurfa", "wpdot", "wbdot"),
                ("slantu", "slantl", "zuperts", "chipre"),
                ("cjor95", "pp95", "ssep", "yyy2"),
                ("xnnc", "cprof", "oring", "cjor0"),
                ("fexpan", "qqmin", "chigamt", "ssi01"),
                ("fexpvs", "sepnose", "ssi95", "rqqmin"),
                ("cjor99", "cj1ave", "rmidin", "rmidout"),
                ("psurfa", "peak", "dminux", "dminlx"),
                ("dolubaf", "dolubafm", "diludom", "diludomm"),
                ("ratsol", "rvsiu", "zvsiu", "rvsid"),
                ("zvsid", "rvsou", "zvsou", "rvsod"),
                ("zvsod", "condno", "psin32", "psin21"),
                ("rq32in", "rq21top", "chilibt", "xdum"),
            ]
            for names in extra_groups:
                vals = reader.read_float_values(4)
                for name, value in zip(names, vals):
                    if name in values:
                        values[name][jj] = value

    machine, coil_currents_ka = _select_machine_currents(ccbrsp, eccurt, nfcoil0, nesum0)

    return AeqdskData(
        filename=str(filename),
        uday=uday.strip(),
        mfvers=mfvers,
        ishot=ishot,
        ktime1=ktime1,
        mastformat=mastformat,
        values=values,
        rseps=rseps,
        jflag=jflag,
        jerror=jerror,
        limloc=limloc,
        nsilop0=nsilop0,
        magpri0=magpri0,
        nfcoil0=nfcoil0,
        nesum0=nesum0,
        mco2v=mco2v,
        mco2r=mco2r,
        rco2r=rco2r,
        rco2v=rco2v,
        dco2r=dco2r,
        dco2v=dco2v,
        csilop=csilop,
        cmpr2=cmpr2,
        ccbrsp=ccbrsp,
        eccurt=eccurt,
        machine=machine,
        coil_currents_ka=coil_currents_ka,
    )
