from __future__ import annotations


def time_name(t: int) -> str:
    """Return HDF5 group name for a time slice."""
    if int(t) < 0:
        return "equilibrium"
    return f"time_{int(t):03d}"
