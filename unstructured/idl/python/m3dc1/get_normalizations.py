from __future__ import annotations

from pathlib import Path

from .hdf5_file_test import hdf5_file_test
from .read_parameter import read_parameter


def get_normalizations(filename: str | Path = "C1.h5") -> tuple[float, float, float, float]:
    """
    Python port of get_normalizations.pro.
    Returns: (b0_norm, n0_norm, l0_norm, ion_mass)
    """
    if not hdf5_file_test(filename):
        return 0.0, 0.0, 0.0, 1.0

    b0 = float(read_parameter("b0_norm", filename=filename))
    n0 = float(read_parameter("n0_norm", filename=filename))
    l0 = float(read_parameter("l0_norm", filename=filename))
    ion_mass = float(read_parameter("ion_mass", filename=filename))
    if ion_mass == 0.0:
        ion_mass = 1.0

    return b0, n0, l0, ion_mass
