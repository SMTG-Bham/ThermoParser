"""Provides pre-sized figures with three sets of axes.

Each function returns a figure and an array of sets of axes. Those with
legend space also return a function to add a pre-positioned legend.

Functions:
    h
    h_top_legend
    h_bottom_legend
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tp import settings
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

style = settings.style()

def h(style=style):
    """A figure with three sets of axes horizontally.

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
    fig = plt.figure(figsize=(27/2.54, 8.3/2.54))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.06, right=0.98,
                        bottom=0.12, top=0.95,
                        wspace=0.3)

    return fig, ax

def h_top_legend(style=style):
    """A figure with three sets of axes horizontally and a legend.

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
    fig = plt.figure(figsize=(27/2.54, 9.1/2.54))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.06, right=0.98,
                        bottom=0.11, top=0.87,
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

        if 'ncol' not in kwargs:
            if 'handles' in kwargs:
                kwargs['ncol'] = len(kwargs['handles'])
            elif 'labels' in kwargs:
                kwargs['ncol'] = len(kwargs['labels'])
            else:
                kwargs['ncol'] = len(ax[1].get_legend_handles_labels()[0])

        legend = ax[1].legend(loc="lower center", bbox_to_anchor=(0.5, 1),
                              *args, **kwargs)

        return legend

    return fig, ax, add_legend

def h_bottom_legend(style=style):
    """A figure with three sets of axes horizontally and a legend.

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
    fig = plt.figure(figsize=(27/2.54, 9.1/2.54))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.06, right=0.98,
                        bottom=0.19, top=0.95,
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

        if 'ncol' not in kwargs:
            if 'handles' in kwargs:
                kwargs['ncol'] = len(kwargs['handles'])
            elif 'labels' in kwargs:
                kwargs['ncol'] = len(kwargs['labels'])
            else:
                kwargs['ncol'] = len(ax[1].get_legend_handles_labels()[0])

        legend = ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.1),
                              *args, **kwargs)

        return legend

    return fig, ax, add_legend
