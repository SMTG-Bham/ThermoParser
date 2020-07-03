"""Provides pre-sized axes.

The idea is each set of axes has a height of 12 cm, so one has a height
of 12 cm, and four_square has a height of 24 cm. Thing is, matplotlib
is in inches not cm, so four_square is like 2 ft tall. One day I'll be
bothered to resize the plots and update the style sheet to match...

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

    two_h

    three_h
    three_h_top_legend
    three_h_bottom_legend

    four_square
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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

            'two_h',

            'three_h',
            'three_h_top_legend',
            'three_h_bottom_legend',

            'four_square']

    return names

def one(style='pretty2'):
    """A figure with a set of axes.

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(13.15, 12))

    plt.subplots_adjust(left=0.22, right=0.95, bottom=0.15, top=0.95)

    return fig, ax

def one_colourbar(style='pretty2'):
    """A figure with a set of axes and colourbar space.

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(16.5, 12))

    plt.subplots_adjust(left=0.22, right=0.95, bottom=0.15, top=0.95)

    return fig, ax

def one_colourbar_small_legend(style='pretty2'):
    """A figure with one set of axes and space for a colourbar and legend.

    Suggestion:
    legend = ax.legend(loc="center left", bbox_to_anchor=(1.22, 0.5))

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(18.5, 12))

    plt.subplots_adjust(left=0.13, right=0.8, bottom=0.13, top=0.96)

    return fig, ax

def one_dos(style='pretty2'):
    """A figure with a set of axes and a DoS-style attachment.

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(15.9, 12))
    grid = GridSpec(1, 4)
    ax = ['', '']
    ax[0] = fig.add_subplot(grid[0, :-1])
    ax[1] = fig.add_subplot(grid[0, -1])

    plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.13,
                        wspace=0)

    return fig, ax

def one_dos_colourbar(style='pretty2'):
    """A figure with axes, DoS-style attachment and colourbar space.

    Suggestion:
    legend = ax[1].legend(loc="center left", bbox_to_anchor=(2, 0.5))

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(18.8, 12))
    grid = GridSpec(1, 5)
    ax = ['', '']
    ax[0] = fig.add_subplot(grid[0, :-2])
    ax[1] = fig.add_subplot(grid[0, -2])

    plt.subplots_adjust(left=0.12, right=1, top=0.95, bottom=0.13,
                        wspace=0)

    return fig, ax

def one_dos_colourbar_small_legend(style='pretty2'):
    """A figure with axes, DoS-style attachment and legend space.

    Suggestion:
    legend = ax[1].legend(loc="center left", bbox_to_anchor=(2, 0.5))

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(22, 12))
    grid = GridSpec(1, 4)
    ax = ['', '']
    ax[0] = fig.add_subplot(grid[0, :-1])
    ax[1] = fig.add_subplot(grid[0, -1])

    plt.subplots_adjust(left=0.12, right=0.72, top=0.95, bottom=0.13,
                        wspace=0)

    return fig, ax

def one_dos_small_legend(style='pretty2'):
    """A figure with axes, DoS-style attachment and legend space.

    Suggestion:
    legend = ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5))

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(20, 12))
    grid = GridSpec(1, 4)
    ax = ['', '']
    ax[0] = fig.add_subplot(grid[0, :-1])
    ax[1] = fig.add_subplot(grid[0, -1])

    plt.subplots_adjust(left=0.12, right=0.78, top=0.95, bottom=0.13,
                        wspace=0)

    return fig, ax

def one_small_legend(style='pretty2'):
    """A figure with one set of axes and space for a legend.

    Suggestion:
    legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(17, 12))

    plt.subplots_adjust(left=0.16, right=0.75, bottom=0.13, top=0.96)

    return fig, ax

def one_medium_legend(style='pretty2'):
    """A figure with one set of axes and space for a legend.

    Suggestion:
    legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig, ax = plt.subplots(figsize=(17.5, 12))

    plt.subplots_adjust(left=0.1, right=0.67, bottom=0.13, top=0.96)

    return fig, ax

def two_h(style='pretty2'):
    """A figure with two sets of axes horizontally.

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(27.2, 12))
    grid = GridSpec(1, 2)
    ax = ['', '']
    ax[0] = fig.add_subplot(grid[0, 0])
    ax[1] = fig.add_subplot(grid[0, 1])

    plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.13,
                        wspace=0.3)

    return fig, ax

def three_h(style='pretty2'):
    """A figure with three sets of axes horizontally.

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(39, 12))
    grid = GridSpec(1, 3)
    ax = ['', '', '']
    ax[0] = fig.add_subplot(grid[0, 0])
    ax[1] = fig.add_subplot(grid[0, 1])
    ax[2] = fig.add_subplot(grid[0, 2])

    plt.subplots_adjust(left=0.06, right=0.97, top=0.95, bottom=0.13,
                        wspace=0.3)

    return fig, ax

def three_h_top_legend(style='pretty2'):
    """A figure with three sets of axes horizontally and a legend.

    Suggestion:
    legend = ax[1].legend(loc="lower center", bbox_to_anchor=(0.5, 1),
                          ncol=number of labels, title=title)

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(39, 13.7))
    grid = GridSpec(1, 3)
    ax = ['', '', '']
    ax[0] = fig.add_subplot(grid[0, 0])
    ax[1] = fig.add_subplot(grid[0, 1])
    ax[2] = fig.add_subplot(grid[0, 2])

    plt.subplots_adjust(left=0.06, right=0.97, top=0.82, bottom=0.11,
                        wspace=0.3)

    return fig, ax

def three_h_bottom_legend(style='pretty2'):
    """A figure with three sets of axes horizontally and a legend.

    Suggestion:
    legend = ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.1),
                          ncol=number_of_labels, title=title)

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        axes
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(39, 13.7))
    grid = GridSpec(1, 3)
    ax = ['', '', '']
    ax[0] = fig.add_subplot(grid[0, 0])
    ax[1] = fig.add_subplot(grid[0, 1])
    ax[2] = fig.add_subplot(grid[0, 2])

    plt.subplots_adjust(left=0.06, right=0.97, top=0.95, bottom=0.24,
                        wspace=0.3)

    return fig, ax

def four_square(style='pretty2'):
    """A figure with four sets of axes in a square.

    Arguments:
        style : str, optional
            style sheet(s). Default: pretty2.

    Returns:
        figure
            figure.
        list
            axes.
    """

    plt.style.use(style)
    fig = plt.figure(figsize=(26.7, 24))

    grid = GridSpec(2,2)
    ax = [['', ''], ['', '']]
    ax[0][0] = fig.add_subplot(grid[0, 0])
    ax[0][1] = fig.add_subplot(grid[0, 1])
    ax[1][0] = fig.add_subplot(grid[1, 0])
    ax[1][1] = fig.add_subplot(grid[1, 1])

    plt.subplots_adjust(hspace=0.23, wspace=0.4,
                        left=0.13, right=0.975,
                        top=0.97, bottom=0.1)

    return fig, ax
