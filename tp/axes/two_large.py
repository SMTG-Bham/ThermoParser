"""Provides pre-sized figures with two sets of axes.

Each function returns a figure and an array of sets of axes. Those with
legend space also return a function to add a pre-positioned legend.
Designed for presentations etc.

Functions
---------

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

default_style = settings.large_style()

def h(style=[]):
    """A figure with two sets of axes horizontally.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp-large.

    Returns
    -------

        figure
            figure.
        axes
            axes.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(17.7, 8.28))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.96,
                        bottom=0.15, top=0.96,
                        wspace=0.3)

    return fig, ax

def h_small_legend(style=[]):
    """A figure with two sets of axes horizontally and space for a legend.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp-large.

    Returns
    -------

        figure
            figure.
        axes
            axes.
        function
            function to add a pre-posistioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(20, 8.28))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.855,
                        bottom=0.15, top=0.96,
                        wspace=0.3)

    def add_legend(*args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments
        ---------

            *args, **kwargs
                passed to ax.legend.

        Returns
        -------

            legend
                legend.
        """

        legend = ax[1].legend(loc="center left", bbox_to_anchor=(0.97, 0.5),
                              *args, **kwargs)

        return legend

    return fig, ax, add_legend

def h_medium_legend(style=[]):
    """A figure with two sets of axes horizontally and space for a legend.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp-large.

    Returns
    -------

        figure
            figure.
        axes
            axes.
        function
            function to add a pre-posistioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(22, 8.28))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.78,
                        bottom=0.15, top=0.96,
                        wspace=0.3)

    def add_legend(*args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments
        ---------

            *args, **kwargs
                passed to ax.legend.

        Returns
        -------

            legend
                legend.
        """

        legend = ax[1].legend(loc="center left", bbox_to_anchor=(0.97, 0.5),
                              *args, **kwargs)

        return legend

    return fig, ax, add_legend

def v(style=[]):
    """A figure with two sets of axes vertically.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp-large.

    Returns
    -------

        figure
            figure.
        axes
            axes.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(8.97, 16.56))
    grid = GridSpec(2, 1)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[1, 0])]

    plt.subplots_adjust(left=0.17, right=0.93,
                        bottom=0.08, top=0.98,
                        hspace=0.2)

    return fig, ax
