"""Provides pre-sized figures with one set of primary axes.

Each function returns a figure and a set of axes, or an array of sets of
axes for those with DoSs. Those with legend space also return a function
to add a pre-positioned legend. Designed for presentations etc.

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

default_style = settings.large_style()

def plain(style=[]):
    """A figure with a set of axes.

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
    fig, ax = plt.subplots(figsize=(8.6, 8.28))

    plt.subplots_adjust(left=0.17, right=0.95,
                        bottom=0.15, top=0.96)

    return fig, ax

def colourbar(style=[]):
    """A figure with a set of axes and colourbar space.

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
    fig, ax = plt.subplots(figsize=(10.9, 8.28))

    plt.subplots_adjust(left=0.18, right=0.94,
                        bottom=0.16, top=0.96)

    return fig, ax

def colourbar_small_legend(style=[]):
    """A figure with one set of axes and space for a colourbar and legend.

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
    fig, ax = plt.subplots(figsize=(14, 8.28))

    plt.subplots_adjust(left=0.14, right=0.74,
                        bottom=0.15, top=0.96)

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

        legend = ax.legend(loc="center left", bbox_to_anchor=(1.33, 0.5),
                           *args, **kwargs)

        return legend

    return fig, ax, add_legend

def dos(style=[]):
    """A figure with a set of axes and a DoS-style attachment.

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
    fig = plt.figure(figsize=(11.1, 8.28))
    grid = GridSpec(1, 4)
    ax = [fig.add_subplot(grid[0, :-1]), fig.add_subplot(grid[0, -1])]

    plt.subplots_adjust(left=0.15, right=0.96,
                        bottom=0.15, top=0.96,
                        wspace=0)

    return fig, ax

def dos_colourbar(style=[]):
    """A figure with axes, DoS-style attachment and colourbar space.

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
    fig = plt.figure(figsize=(13.5, 8.3))
    grid = GridSpec(1, 9)
    ax = [fig.add_subplot(grid[0, :6]), fig.add_subplot(grid[0, 6:])]

    plt.subplots_adjust(left=0.12, right=0.87,
                        bottom=0.15, top=0.96,
                        wspace=0)

    return fig, ax

def dos_colourbar_small_legend(style=[]):
    """A figure with axes, DoS-style attachment and legend space.

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
    fig = plt.figure(figsize=(16, 8.28))
    grid = GridSpec(1, 9)
    ax = [fig.add_subplot(grid[0, :6]), fig.add_subplot(grid[0, 6:])]

    plt.subplots_adjust(left=0.09, right=0.72,
                        bottom=0.15, top=0.96,
                        wspace=0)

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

        legend = ax[1].legend(loc="center left", bbox_to_anchor=(1.8, 0.5),
                              *args, **kwargs)

        return legend

    return fig, ax, add_legend

def dos_small_legend(style=[]):
    """A figure with axes, DoS-style attachment and legend space.

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
    fig = plt.figure(figsize=(13.45, 8.3))
    grid = GridSpec(1, 4)
    ax = [fig.add_subplot(grid[0, :-1]), fig.add_subplot(grid[0, -1])]

    plt.subplots_adjust(left=0.11, right=0.78,
                        bottom=0.15, top=0.96,
                        wspace=0)

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

        legend = ax[1].legend(loc="center left", bbox_to_anchor=(0.87, 0.5),
                              *args, **kwargs)

        return legend

    return fig, ax, add_legend

def small_legend(style=[]):
    """A figure with one set of axes and space for a legend.

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
    fig, ax = plt.subplots(figsize=(12, 8.28))

    plt.subplots_adjust(left=0.17, right=0.87,
                        bottom=0.15, top=0.96)

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

        legend = ax.legend(loc="center left", bbox_to_anchor=(0.98, 0.5),
                           *args, **kwargs)

        return legend

    return fig, ax, add_legend

def medium_legend(style=[]):
    """A figure with one set of axes and space for a legend.

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
    fig, ax = plt.subplots(figsize=(13.5, 8.28))

    plt.subplots_adjust(left=0.16, right=0.78,
                        bottom=0.15, top=0.96)

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

        legend = ax.legend(loc="center left", bbox_to_anchor=(0.98, 0.5),
                           *args, **kwargs)

        return legend

    return fig, ax, add_legend

def wide(style=[]):
    """A figure with one set of axes and more space on the left.

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
    fig, ax = plt.subplots(figsize=(9.9, 8.3))

    plt.subplots_adjust(left=0.28, right=0.96,
                        bottom=0.15, top=0.96)

    return fig, ax

def wide_large_legend(style=[]):
    """A figure with one set of axes and space for a legend.

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
    fig, ax = plt.subplots(figsize=(14.6, 8.28))

    plt.subplots_adjust(left=0.2, right=0.66,
                        bottom=0.15, top=0.96)

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

        legend = ax.legend(loc="center left", bbox_to_anchor=(0.98, 0.5),
                           *args, **kwargs)

        return legend

    return fig, ax, add_legend
