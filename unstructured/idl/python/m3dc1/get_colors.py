from __future__ import annotations

from matplotlib import pyplot as plt


def get_colors(maxcolors: int = 10) -> list:
    """
    Python port of get_colors.pro (returns a list of plotting colors).
    """
    n = max(maxcolors, 10)
    cmap = plt.get_cmap("tab20")
    return [cmap(i % cmap.N) for i in range(n)]
