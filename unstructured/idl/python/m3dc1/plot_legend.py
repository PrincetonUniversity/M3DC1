from __future__ import annotations

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D


def plot_legend(
    names: list[str],
    linestyles: list[str] | None = None,
    colors: list | None = None,
    nolegend: bool = False,
) -> None:
    """
    Python replacement for plot_legend.pro.
    """
    if nolegend or not names:
        return
    n = len(names)
    if linestyles is None:
        linestyles = ["-"] * n
    if colors is None:
        colors = [None] * n
    handles = [
        Line2D([0], [0], linestyle=linestyles[i], color=colors[i], label=names[i]) for i in range(n)
    ]
    plt.legend(handles=handles, loc="best", frameon=False)
