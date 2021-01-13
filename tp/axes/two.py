"""Provides pre-sized figures with two sets of axes.

Each function returns a figure and an array of sets of axes. Those with
legend space also return a function to add a pre-positioned legend.
Designed for papers.

Functions
---------

    h
    h_small_legend
    h_medium_legend
    v
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tp
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

default_style = tp.settings.style()

def h(style=[]):
    """A figure with two sets of axes horizontally.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp.

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
    fig = plt.figure(figsize=(17.5/2.54, 8.3/2.54))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.98,
                        bottom=0.12, top=0.95,
                        wspace=0.3)

    return fig, ax

def h_small_legend(style=[]):
    """A figure with two sets of axes horizontally and space for a legend.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp.

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
    fig = plt.figure(figsize=(19/2.54, 8.3/2.54))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.91,
                        bottom=0.12, top=0.95,
                        wspace=0.3)

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
            legend = ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5),
                                  *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend

def h_medium_legend(style=[]):
    """A figure with two sets of axes horizontally and space for a legend.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp.

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
    fig = plt.figure(figsize=(19.8/2.54, 8.3/2.54))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.88,
                        bottom=0.12, top=0.95,
                        wspace=0.3)

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
            legend = ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5),
                                  *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend

def v(style=[]):
    """A figure with two sets of axes vertically.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp.

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
    fig = plt.figure(figsize=(8.9/2.54, 16.6/2.54))
    grid = GridSpec(2, 1)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[1, 0])]

    plt.subplots_adjust(left=0.15, right=0.95,
                        bottom=0.06, top=0.98,
                        hspace=0.15)

    return fig, ax
