"""Provides pre-sized figures with four axes.

Each function returns a figure and an array of sets of axes. Those with
legend space also return a function to add a pre-positioned legend.
Designed for papers.

Functions
---------

    square
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tp import settings
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

default_style = settings.style()

def square(style=[]):
    """A figure with four sets of axes in a square.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp.

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
    fig = plt.figure(figsize=(18.3/2.54, 16.6/2.54))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
          [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97,
                        hspace=0.15, wspace=0.3)

    return fig, ax

def square_legend(style=[]):
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
    fig = plt.figure(figsize=(25/2.54, 16.56/2.54))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
          [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.13, right=0.8,
                        bottom=0.08, top=0.97,
                        hspace=0.22, wspace=0.31)

    def add_legend(custom=False, *args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments
        ---------

            custom : bool, optional
                enable manual editing of handles and labels arguments.
                Default: False.
            *args, **kwargs
                passed to ax.legend.

        Returns
        -------

            legend
                legend.
        """

        if custom:
            legend = ax[0].legend(loc="center left", bbox_to_anchor=(0.97, 0),
                                  *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[0].legend(loc="center left", bbox_to_anchor=(0.97, 0),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend
