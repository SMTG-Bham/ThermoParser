"""Provides pre-sized figures with three sets of axes.

Each function returns a figure and an array of sets of axes. Those with
legend space also return a function to add a pre-positioned legend.
Designed for presentations etc.

Functions
---------

    h
    h_top_legend
    h_big_top_legend
    h_bottom_legend
    h_big_bottom_legend

    square_q1_legend
    square_q2_legend
    square_q3_legend
    square_q4_legend
"""

from math import ceil
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tp
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

default_style = tp.settings.large_style()

def h(style=[]):
    """Axes horizontally.

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
    fig = plt.figure(figsize=(27.3, 8.28))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.08, right=0.97,
                        bottom=0.15, top=0.96,
                        wspace=0.32)

    return fig, ax

def h_top_legend(style=[]):
    """Axes horizontally with a legend above.

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
    fig = plt.figure(figsize=(27.3, 10.3))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.08, right=0.97,
                        bottom=0.15, top=0.8,
                        wspace=0.32)

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

        if 'ncol' not in kwargs:
            if 'handles' in kwargs:
                kwargs['ncol'] = len(kwargs['handles'])
            elif 'labels' in kwargs:
                kwargs['ncol'] = len(kwargs['labels'])
            else:
                kwargs['ncol'] = len(ax[1].get_legend_handles_labels()[0])

        if custom:
            legend = ax[1].legend(loc="lower center", bbox_to_anchor=(0.5, 1),
                                  *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1].legend(loc="lower center", bbox_to_anchor=(0.5, 1),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend

def h_big_top_legend(style=[]):
    """Axes horizontally with a double-height legend above.

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
    fig = plt.figure(figsize=(27.3, 10.8))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.08, right=0.97,
                        bottom=0.14, top=0.76,
                        wspace=0.32)

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

        if 'ncol' not in kwargs:
            if 'handles' in kwargs:
                kwargs['ncol'] = ceil(len(kwargs['handles'])/2)
            elif 'labels' in kwargs:
                kwargs['ncol'] = ceil(len(kwargs['labels'])/2)
            else:
                kwargs['ncol'] = ceil(len(ax[1].get_legend_handles_labels()[0])
                                                                            /2)

        if custom:
            legend = ax[1].legend(loc="lower center", bbox_to_anchor=(0.5, 0.98),
                                  *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1].legend(loc="lower center", bbox_to_anchor=(0.5, 0.98),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend

def h_bottom_legend(style=[]):
    """Axes horizontally with a legend below.

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
    fig = plt.figure(figsize=(27.3, 10.3))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.08, right=0.97,
                        bottom=0.3, top=0.95,
                        wspace=0.32)

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

        if 'ncol' not in kwargs:
            if 'handles' in kwargs:
                kwargs['ncol'] = len(kwargs['handles'])
            elif 'labels' in kwargs:
                kwargs['ncol'] = len(kwargs['labels'])
            else:
                kwargs['ncol'] = len(ax[1].get_legend_handles_labels()[0])

        if custom:
            legend = ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),
                                  *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend

def h_big_bottom_legend(style=[]):
    """Axes horizontally with a double-height legend below.

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
    fig = plt.figure(figsize=(27.3, 10.8))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.08, right=0.97,
                        bottom=0.34, top=0.96,
                        wspace=0.32)

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

        if 'ncol' not in kwargs:
            if 'handles' in kwargs:
                kwargs['ncol'] = ceil(len(kwargs['handles'])/2)
            elif 'labels' in kwargs:
                kwargs['ncol'] = ceil(len(kwargs['labels'])/2)
            else:
                kwargs['ncol'] = ceil(len(ax[1].get_legend_handles_labels()[0])
                                                                            /2)

        if custom:
            legend = ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.14),
                                  *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.14),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend

def square_q1_legend(style=[]):
    """Axes in a square with a legend in the first quadrant.

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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.3, 16.56))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), None],
          [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.13, right=0.96,
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
            legend = ax[0][0].legend(loc="center", bbox_to_anchor=(1.81, 0.5),
                                     *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[0][0].legend(loc="center", bbox_to_anchor=(1.81, 0.5),
                                     handles=handles, labels=labels,
                                     *args, **kwargs)

        return legend

    return fig, ax, add_legend

def square_q2_legend(style=[]):
    """Axes in a square with a legend in the second quadrant.

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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.3, 16.56))

    grid = GridSpec(2,2)
    ax = [[None,                        fig.add_subplot(grid[0, 1])],
          [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.13, right=0.96,
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
            legend = ax[0][1].legend(loc="center", bbox_to_anchor=(-0.81, 0.5),
                                     *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[0][1].legend(loc="center", bbox_to_anchor=(-0.81, 0.5),
                                     handles=handles, labels=labels,
                                     *args, **kwargs)

        return legend

    return fig, ax, add_legend

def square_q3_legend(style=[]):
    """Axes in a square with a legend in the third quadrant.

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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.3, 16.56))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
          [None,                        fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.13, right=0.96,
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
            legend = ax[1][1].legend(loc="center", bbox_to_anchor=(-0.81, 0.5),
                                     *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1][1].legend(loc="center", bbox_to_anchor=(-0.81, 0.5),
                                     handles=handles, labels=labels,
                                     *args, **kwargs)

        return legend

    return fig, ax, add_legend

def square_q4_legend(style=[]):
    """Axes in a square with a legend in the fourth quadrant.

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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.3, 16.56))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
          [fig.add_subplot(grid[1, 0]), None]]

    plt.subplots_adjust(left=0.13, right=0.96,
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
            legend = ax[1][0].legend(loc="center", bbox_to_anchor=(1.81, 0.5),
                                     *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1][0].legend(loc="center", bbox_to_anchor=(1.81, 0.5),
                                     handles=handles, labels=labels,
                                     *args, **kwargs)

        return legend

    return fig, ax, add_legend
