"""Provides pre-sized axes.

The style is based somewhat on the Nature guidelines (we can all
dream!), where the axes.one figure has a width of 8.6 cm. All axes are 
square, and subsequent figures maintain the height of axes.one per plot
(single height plots have a height of 8.3 cm whilst double height plots
have a height of 16.6 cm, which more effectively maintains relative
scale, as legends etc. are normally added on the side). As per the
guidelines, font size is 8 pt.

Functions:
    one
    one_colourbar
    one_colourbar_small_legend
    one_dos
    one_dos_colourbar
    one_dos_colourbar_small_legend
    one_dos_small_legend
    one_small_legend
    one_medium_legend
    one_wide
    one_wide_large_legend

    two_h
    two_h_small_legend
    two_h_medium_legend
    two_v

    three_h
    three_h_top_legend
    three_h_bottom_legend

    four_square
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tp import settings
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

def __dir__():
    names =['one',
            'one_colourbar',
            'one_colourbar_small_legend',
            'one_dos',
            'one_dos_colourbar',
            'one_dos_colourbar_small_legend',
            'one_dos_small_legend',
            'one_small_legend',
            'one_medium_legend',
            'one_wide',
            'one_wide_large_legend',

            'two_h',
            'two_h_small_legend',
            'two_h_medium_legend',
            'two_v',

            'three_h',
            'three_h_top_legend',
            'three_h_bottom_legend',

            'four_square']

    return names

style = settings.style()

def one(style=style):
    """A figure with a set of axes.

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
    fig, ax = plt.subplots(figsize=(8.6/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.15, right=0.95,
                        bottom=0.12, top=0.95)

    return fig, ax

def one_colourbar(style=style):
    """A figure with a set of axes and colourbar space.

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
    fig, ax = plt.subplots(figsize=(10.5/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.14, right=0.96,
                        bottom=0.12, top=0.95)

    return fig, ax

def one_colourbar_small_legend(style=style):
    """A figure with one set of axes and space for a colourbar and legend.

    Suggestion:
    legend = ax.legend(loc="center left", bbox_to_anchor=(1.25, 0.5))

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
    fig, ax = plt.subplots(figsize=(11.8/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.12, right=0.85,
                        bottom=0.12, top=0.95)

    return fig, ax

def one_dos(style=style):
    """A figure with a set of axes and a DoS-style attachment.

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
    fig = plt.figure(figsize=(11.05/2.54, 8.3/2.54))
    grid = GridSpec(1, 4)
    ax = [fig.add_subplot(grid[0, :-1]), fig.add_subplot(grid[0, -1])]

    plt.subplots_adjust(left=0.12, right=0.95,
                        bottom=0.12, top=0.95,
                        wspace=0)

    return fig, ax

def one_dos_colourbar(style=style):
    """A figure with axes, DoS-style attachment and colourbar space.

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
    fig = plt.figure(figsize=(12.2/2.54, 8.3/2.54))
    grid = GridSpec(1, 9)
    ax = [fig.add_subplot(grid[0, :6]), fig.add_subplot(grid[0, 6:])]

    plt.subplots_adjust(left=0.1, right=0.95,
                        bottom=0.12, top=0.95,
                        wspace=0)

    return fig, ax

def one_dos_colourbar_small_legend(style=style):
    """A figure with axes, DoS-style attachment and legend space.

    Suggestion:
    legend = ax[1].legend(loc="center left", bbox_to_anchor=(1.75, 0.5))

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
    fig = plt.figure(figsize=(13.8/2.54, 8.3/2.54))
    grid = GridSpec(1, 9)
    ax = [fig.add_subplot(grid[0, :6]), fig.add_subplot(grid[0, 6:])]

    plt.subplots_adjust(left=0.1, right=0.85,
                        bottom=0.12, top=0.95,
                        wspace=0)

    return fig, ax

def one_dos_small_legend(style=style):
    """A figure with axes, DoS-style attachment and legend space.

    Suggestion:
    legend = ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5))

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
    fig = plt.figure(figsize=(12.5/2.54, 8.3/2.54))
    grid = GridSpec(1, 4)
    ax = [fig.add_subplot(grid[0, :-1]), fig.add_subplot(grid[0, -1])]

    plt.subplots_adjust(left=0.12, right=0.85,
                        bottom=0.12, top=0.95,
                        wspace=0)

    return fig, ax

def one_small_legend(style=style):
    """A figure with one set of axes and space for a legend.

    Suggestion:
    legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

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
    fig, ax = plt.subplots(figsize=(10.15/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.15, right=0.83,
                        bottom=0.12, top=0.95)

    return fig, ax

def one_medium_legend(style=style):
    """A figure with one set of axes and space for a legend.

    Suggestion:
    legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

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
    fig, ax = plt.subplots(figsize=(10.3/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.12, right=0.79,
                        bottom=0.12, top=0.95)

    return fig, ax

def one_wide(style=style):
    """A figure with one set of axes and more space on the left.

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
    fig, ax = plt.subplots(figsize=(9.2/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.2, right=0.95,
                        bottom=0.12, top=0.95)

    return fig, ax

def one_wide_large_legend(style=style):
    """A figure with one set of axes and space for a legend.

    Suggestion:
    legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

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
    fig, ax = plt.subplots(figsize=(11.5/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.15, right=0.75,
                        bottom=0.12, top=0.95)

    return fig, ax

def two_h(style=style):
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

def two_h_small_legend(style=style):
    """A figure with two sets of axes horizontally and space for a legend.

    Suggestion:
    legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

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
    fig = plt.figure(figsize=(19/2.54, 8.3/2.54))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.91,
                        bottom=0.12, top=0.95,
                        wspace=0.3)

    return fig, ax

def two_h_medium_legend(style=style):
    """A figure with two sets of axes horizontally and space for a legend.

    Suggestion:
    legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

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
    fig = plt.figure(figsize=(19.8/2.54, 8.3/2.54))
    grid = GridSpec(1, 2)
    ax = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]

    plt.subplots_adjust(left=0.08, right=0.88,
                        bottom=0.12, top=0.95,
                        wspace=0.3)

    return fig, ax

def two_v(style=style):
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

def three_h(style=style):
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

def three_h_top_legend(style=style):
    """A figure with three sets of axes horizontally and a legend.

    Suggestion:
    legend = ax[1].legend(loc="lower center", bbox_to_anchor=(0.5, 1),
                          ncol=number of labels, title=title)

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
    fig = plt.figure(figsize=(27/2.54, 9.1/2.54))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.06, right=0.98,
                        bottom=0.11, top=0.87,
                        wspace=0.3)

    return fig, ax

def three_h_bottom_legend(style=style):
    """A figure with three sets of axes horizontally and a legend.

    Suggestion:
    legend = ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.1),
                          ncol=number_of_labels, title=title)

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
    fig = plt.figure(figsize=(27/2.54, 9.1/2.54))
    grid = GridSpec(1, 3)
    ax = [fig.add_subplot(grid[0, 0]),
          fig.add_subplot(grid[0, 1]),
          fig.add_subplot(grid[0, 2])]

    plt.subplots_adjust(left=0.06, right=0.98,
                        bottom=0.19, top=0.95,
                        wspace=0.3)

    return fig, ax

def four_square(style=style):
    """A figure with four sets of axes in a square.

    Arguments:
        style : str, optional
            style sheet(s). Default: tp.

    Returns:
        figure
            figure.
        list
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(18.3/2.54, 16.6/2.54))

    grid = GridSpec(2,2)
    ax = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
          [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97,
                        hspace=0.15, wspace=0.3)


    return fig, ax
