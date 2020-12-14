"""Provides pre-sized figures with one set of primary axes.

Each function returns a figure and a set of axes, or an array of sets of
axes for those with DoSs. Those with legend space also return a function
to add a pre-positioned legend.

Functions
---------

    plain
    colourbar
    colourbar_small_legend
    dos
    dos_colourbar
    dos_colourbar_small_legend
    dos_small_legend
    small_legend
    medium_legend
    wide
    wide_large_legend
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tp import settings
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

style = settings.style()

def plain(style=style):
    """A figure with a set of axes.

    Arguments
    ---------

        style : str, optional
            style sheet(s). Default: tp.

    Returns
    -------

        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(8.6/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.15, right=0.95,
                        bottom=0.12, top=0.95)

    return fig, ax

def colourbar(style=style):
    """A figure with a set of axes and colourbar space.

    Arguments
    ---------

        style : str, optional
            style sheet(s). Default: tp.

    Returns
    -------

        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(10.5/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.14, right=0.96,
                        bottom=0.12, top=0.95)

    return fig, ax

def colourbar_small_legend(style=style):
    """A figure with one set of axes and space for a colourbar and legend.

    Arguments
    ---------

        style : str, optional
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

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(11.8/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.12, right=0.85,
                        bottom=0.12, top=0.95)

    def add_legend(*args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments
        ---------

            *args, **kwargs : optional
                passed to ax.legend.

        Returns
        -------

            legend
                legend.
        """

        legend = ax.legend(loc="center left", bbox_to_anchor=(1.27, 0.5),
                           *args, **kwargs)

        return legend

    return fig, ax, add_legend

def dos(style=style):
    """A figure with a set of axes and a DoS-style attachment.

    Arguments
    ---------

        style : str, optional
            style sheet(s). Default: tp.

    Returns
    -------

        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(11.05/2.54, 8.3/2.54))
    grid = GridSpec(1, 4)
    ax = [fig.add_subplot(grid[0, :-1]), fig.add_subplot(grid[0, -1])]

    plt.subplots_adjust(left=0.12, right=0.95,
                        bottom=0.12, top=0.95,
                        wspace=0)

    return fig, ax

def dos_colourbar(style=style):
    """A figure with axes, DoS-style attachment and colourbar space.

    Arguments
    ---------

        style : str, optional
            style sheet(s). Default: tp.

    Returns
    -------

        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(12.2/2.54, 8.3/2.54))
    grid = GridSpec(1, 9)
    ax = [fig.add_subplot(grid[0, :6]), fig.add_subplot(grid[0, 6:])]

    plt.subplots_adjust(left=0.1, right=0.95,
                        bottom=0.12, top=0.95,
                        wspace=0)

    return fig, ax

def dos_colourbar_small_legend(style=style):
    """A figure with axes, DoS-style attachment and legend space.

    Arguments
    ---------

        style : str, optional
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

    plt.style.use(style)
    fig = plt.figure(figsize=(13.8/2.54, 8.3/2.54))
    grid = GridSpec(1, 9)
    ax = [fig.add_subplot(grid[0, :6]), fig.add_subplot(grid[0, 6:])]

    plt.subplots_adjust(left=0.1, right=0.85,
                        bottom=0.12, top=0.95,
                        wspace=0)

    def add_legend(*args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments
        ---------

            *args, **kwargs : optional
                passed to ax.legend.

        Returns
        -------

            legend
                legend.
    """

        legend = ax[1].legend(loc="center left", bbox_to_anchor=(1.75, 0.5),
                              *args, **kwargs)

        return legend

    return fig, ax, add_legend

def dos_small_legend(style=style):
    """A figure with axes, DoS-style attachment and legend space.

    Arguments
    ---------

        style : str, optional
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

    plt.style.use(style)
    fig = plt.figure(figsize=(12.5/2.54, 8.3/2.54))
    grid = GridSpec(1, 4)
    ax = [fig.add_subplot(grid[0, :-1]), fig.add_subplot(grid[0, -1])]

    plt.subplots_adjust(left=0.12, right=0.85,
                        bottom=0.12, top=0.95,
                        wspace=0)

    def add_legend(*args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments
        ---------

            *args, **kwargs : optional
                passed to ax.legend.

        Returns
        -------

            legend
                legend.
        """

        legend = ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5),
                              *args, **kwargs)

        return legend

    return fig, ax, add_legend

def small_legend(style=style):
    """A figure with one set of axes and space for a legend.

    Arguments
    ---------

        style : str, optional
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

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(10.15/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.15, right=0.83,
                        bottom=0.12, top=0.95)

    def add_legend(*args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments
        ---------

            *args, **kwargs : optional
                passed to ax.legend.

        Returns
        -------

            legend
                legend.
        """

        legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0.5),
                           *args, **kwargs)

        return legend

    return fig, ax, add_legend

def medium_legend(style=style):
    """A figure with one set of axes and space for a legend.

    Arguments
    ---------

        style : str, optional
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

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(10.3/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.12, right=0.79,
                        bottom=0.12, top=0.95)

    def add_legend(*args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments
        ---------

            *args, **kwargs : optional
                passed to ax.legend.

        Returns
        -------

            legend
                legend.
        """

        legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0.5),
                           *args, **kwargs)

        return legend

    return fig, ax, add_legend

def wide(style=style):
    """A figure with one set of axes and more space on the left.

    Arguments
    ---------

        style : str, optional
            style sheet(s). Default: tp.

    Returns
    -------

        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(9.2/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.2, right=0.95,
                        bottom=0.12, top=0.95)

    return fig, ax

def wide_large_legend(style=style):
    """A figure with one set of axes and space for a legend.

    Arguments
    ---------

        style : str, optional
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

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(11.5/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.15, right=0.75,
                        bottom=0.12, top=0.95)

    def add_legend(*args, **kwargs):
        """Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.).

        Arguments
        ---------

            *args, **kwargs : optional
                passed to ax.legend.

        Returns
        -------

            legend
                legend.
        """

        legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0.5),
                           *args, **kwargs)

        return legend

    return fig, ax, add_legend
