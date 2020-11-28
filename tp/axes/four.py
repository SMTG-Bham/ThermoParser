"""Provides pre-sized figures with four axes.

Each function returns a figure and an array of sets of axes. Those with
legend space also return a function to add a pre-positioned legend.

Functions:
    square
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tp import settings
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

style = settings.style()

def square(style=style):
    """A figure with four sets of axes in a square.

    Arguments:
        style : str, optional
            style sheet(s). Default: tp.

    Returns:
        figure
            figure.
        list
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(18.3/2.54, 16.6/2.54))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
          [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97,
                        hspace=0.15, wspace=0.3)

    return fig, ax
