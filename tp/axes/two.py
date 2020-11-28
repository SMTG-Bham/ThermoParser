"""Provides pre-sized figures with two sets of axes.

Each function returns a figure and an array of sets of axes. Those with
legend space also return a function to add a pre-positioned legend.

Functions:
    h
    h_small_legend
    h_medium_legend
    v
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tp import settings
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

style = settings.style()

def h(style=style):
    """A figure with two sets of axes horizontally.

    Arguments:
        style : str, optional
            style sheet(s). Default: tp.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(17.5/2.54, 8.3/2.54))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.98,
                        bottom=0.12, top=0.95,
                        wspace=0.3)

    return fig, ax

def h_small_legend(style=style):
    """A figure with two sets of axes horizontally and space for a legend.

    Arguments:
        style : str, optional
            style sheet(s). Default: tp.

    Returns:
        figure
            figure.
        axes
            axes.
        function
            function to add a pre-posistioned legend.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(19/2.54, 8.3/2.54))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.91,
                        bottom=0.12, top=0.95,
                        wspace=0.3)

    def add_legend(*args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments:
            *args, **kwargs : optional
                passed to ax.legend.

        Returns:
            legend
                legend.
        """

        legend = ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5),
                              *args, **kwargs)

        return legend

    return fig, ax, add_legend

def h_medium_legend(style=style):
    """A figure with two sets of axes horizontally and space for a legend.

    Arguments:
        style : str, optional
            style sheet(s). Default: tp.

    Returns:
        figure
            figure.
        axes
            axes.
        function
            function to add a pre-posistioned legend.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(19.8/2.54, 8.3/2.54))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.88,
                        bottom=0.12, top=0.95,
                        wspace=0.3)

    def add_legend(*args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments:
            *args, **kwargs : optional
                passed to ax.legend.

        Returns:
            legend
                legend.
        """

        legend = ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5),
                              *args, **kwargs)

        return legend

    return fig, ax, add_legend

def v(style=style):
    """A figure with two sets of axes vertically.

    Arguments:
        style : str, optional
            style sheet(s). Default: tp.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(8.9/2.54, 16.6/2.54))
    grid = GridSpec(2, 1)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[1, 0])]

    plt.subplots_adjust(left=0.15, right=0.95,
                        bottom=0.06, top=0.98,
                        hspace=0.15)

    return fig, ax
