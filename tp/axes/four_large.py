"""Provides pre-sized figures with four axes.

Each function returns a figure and an array of sets of axes. Those with
legend space also return a function to add a pre-positioned legend.
Designed for presentations etc.

Functions
---------

    square
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tp import settings
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

default_style = settings.large_style()

def square(style=[]):
    """A figure with four sets of axes in a square.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp-large.

    Returns
    -------

        figure
            figure.
        list
            axes.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.3, 16.56))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
          [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.13, right=0.96,
                        bottom=0.08, top=0.97,
                        hspace=0.22, wspace=0.31)

    return fig, ax
