"""Provides pre-sized figures with three sets of axes.

Each function returns a figure and an array of sets of axes. Those with
legend space also return a function to add a pre-positioned legend.
Designed for papers.

Functions
---------

    h
    h_top_legend
    h_bottom_legend

    square_q1_legend
    square_q2_legend
    square_q3_legend
    square_q4_legend
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tp
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

default_style = tp.settings.style()

def h(style=[]):
    """A figure with three sets of axes horizontally.

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
    fig = plt.figure(figsize=(27/2.54, 8.3/2.54))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.06, right=0.98,
                        bottom=0.12, top=0.95,
                        wspace=0.3)

    return fig, ax

def h_top_legend(style=[]):
    """A figure with three sets of axes horizontally and a legend.

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
    fig = plt.figure(figsize=(27/2.54, 9.1/2.54))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.06, right=0.98,
                        bottom=0.11, top=0.87,
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

def h_bottom_legend(style=[]):
    """A figure with three sets of axes horizontally and a legend.

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
    fig = plt.figure(figsize=(27/2.54, 9.1/2.54))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.06, right=0.98,
                        bottom=0.19, top=0.95,
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

        if 'ncol' not in kwargs:
            if 'handles' in kwargs:
                kwargs['ncol'] = len(kwargs['handles'])
            elif 'labels' in kwargs:
                kwargs['ncol'] = len(kwargs['labels'])
            else:
                kwargs['ncol'] = len(ax[1].get_legend_handles_labels()[0])

        if custom:
            legend = ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.1),
                                  *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.1),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend

def square_q1_legend(style=[]):
    """Axes in a square with a legend in the first quadrant.

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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.3/2.54, 16.6/2.54))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), None],
          [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97,
                        hspace=0.15, wspace=0.3)

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
            legend = ax[0][0].legend(loc="center", bbox_to_anchor=(1.8, 0.5),
                                     *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[0][0].legend(loc="center", bbox_to_anchor=(1.8, 0.5),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend

def square_q2_legend(style=[]):
    """Axes in a square with a legend in the second quadrant.

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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.3/2.54, 16.6/2.54))

    grid = GridSpec(2,2)
    ax = [[None,                        fig.add_subplot(grid[0, 1])],
          [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97,
                        hspace=0.15, wspace=0.3)

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
            legend = ax[0][1].legend(loc="center", bbox_to_anchor=(-0.8, 0.5),
                                     *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[0][1].legend(loc="center", bbox_to_anchor=(-0.8, 0.5),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend

def square_q3_legend(style=[]):
    """Axes in a square with a legend in the third quadrant.

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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.3/2.54, 16.6/2.54))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
          [None,                        fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97,
                        hspace=0.15, wspace=0.3)

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
            legend = ax[1][1].legend(loc="center", bbox_to_anchor=(-0.8, 0.5),
                                     *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1][1].legend(loc="center", bbox_to_anchor=(-0.8, 0.5),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend

def square_q4_legend(style=[]):
    """Axes in a square with a legend in the fourth quadrant.

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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.3/2.54, 16.6/2.54))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
          [fig.add_subplot(grid[1, 0]), None]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97,
                        hspace=0.15, wspace=0.3)

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
            legend = ax[1][0].legend(loc="center", bbox_to_anchor=(1.8, 0.5),
                                     *args, **kwargs)
        else:
            handles, labels = tp.axes.legend.consolidate(ax)
            legend = ax[1][0].legend(loc="center", bbox_to_anchor=(1.8, 0.5),
                                  handles=handles, labels=labels,
                                  *args, **kwargs)

        return legend

    return fig, ax, add_legend
