"""Provides pre-sized figures for presentation-sized figures.

Each function returns a figure, a set of axes or an array of sets of
axes and an add_legend function. The legend function comes with a
choice of positions and takes normal ax.legend arguments. Drawing large
figures is significantly slower than small ones.

Functions
---------

    one
    one_colourbar
    one_dos
    one_dos_colourbar
    two_h
    two_h_colourbars
    two_v
    two_v_colourbars
    three_h
    three_h_colourbars
    three_square
    three_square_colourbars
    four_square
    four_square_colourbars
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tp
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

default_style = tp.settings.large_style()

def one(style=[]):
    """A figure with a set of axes.

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
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig, ax = plt.subplots(figsize=(8.6, 8.3))

    plt.subplots_adjust(left=0.15, right=0.95,
                        bottom=0.12, top=0.95)

    names = [['in', 'inside'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'loc':            'center left',
                  'bbox_to_anchor': (1, 0.5)},
                 {'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1)},
                 {'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.12)}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'right')

    return fig, ax, add_legend

def one_colourbar(style=[]):
    """A figure with a set of axes and colourbar space.

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
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig, ax = plt.subplots(figsize=(10.5, 8.3))

    plt.subplots_adjust(left=0.14, right=0.96,
                        bottom=0.12, top=0.95)

    names = [['in', 'inside'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'loc':            'center left',
                  'bbox_to_anchor': (1.4, 0.5)},
                 {'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1)},
                 {'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.12)}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'right')

    return fig, ax, add_legend

def one_dos(style=[]):
    """A figure with a set of axes and a DoS-style attachment.

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
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(11.05, 8.3))
    grid = GridSpec(1, 4)
    ax = [fig.add_subplot(grid[0, :-1]), fig.add_subplot(grid[0, -1])]
    fig.__dict__['dos'] = True # Helps positioning colourbars

    plt.subplots_adjust(left=0.12, right=0.95,
                        bottom=0.12, top=0.95,
                        wspace=0)

    names = [['in', 'inside'],
             ['dos'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'axes':           1,
                  'loc':            'center left',
                  'bbox_to_anchor': (1, 0.5)},
                 {'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1)},
                 {'axes':           0,
                 'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.12)}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'right')

    return fig, ax, add_legend

def one_dos_colourbar(style=[]):
    """A figure with axes, DoS-style attachment and colourbar space.

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
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(12.15, 8.3))
    grid = GridSpec(1, 9)
    ax = [fig.add_subplot(grid[0, :6]), fig.add_subplot(grid[0, 6:])]
    fig.__dict__['dos'] = True # Helps positioning colourbars

    plt.subplots_adjust(left=0.1, right=0.95,
                        bottom=0.12, top=0.95,
                        wspace=0)

    names = [['in', 'inside'],
             ['dos'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'axes':           1,
                  'loc':            'center left',
                  'bbox_to_anchor': (1.75, 0.5)},
                 {'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1)},
                 {'axes':           0,
                 'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.12)}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'right')

    return fig, ax, add_legend

def two_h(style=[]):
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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.4, 8.3))
    grid = GridSpec(1, 7)
    ax = [fig.add_subplot(grid[0, :3]), fig.add_subplot(grid[0, 4:])]

    plt.subplots_adjust(left=0.08, right=0.98,
                        bottom=0.12, top=0.95)

    names = [['lefthand'],
             ['righthand', 'in', 'inside'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'axes':           1,
                  'loc':            'center left',
                  'bbox_to_anchor': (1, 0.5)},
                 {'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (1.2, 1),
                  'ncol':           8},
                 {'axes':           0,
                  'loc':            'upper center',
                  'bbox_to_anchor': (1.2, -0.12),
                  'ncol':           8}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'above')

    return fig, ax, add_legend

def two_h_colourbars(style=[]):
    """A figure with two sets of axes horizontally with colourbars.

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
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(23.4, 8.3))
    grid = GridSpec(1, 12)
    ax = [fig.add_subplot(grid[0, :5]), fig.add_subplot(grid[0, 7:])]

    plt.subplots_adjust(left=0.08, right=0.98,
                        bottom=0.12, top=0.95)

    names = [['lefthand'],
             ['righthand', 'in', 'inside'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'axes':           1,
                  'loc':            'center left',
                  'bbox_to_anchor': (1.4, 0.5)},
                 {'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (1.4, 1),
                  'ncol':           8},
                 {'axes':           0,
                  'loc':            'upper center',
                  'bbox_to_anchor': (1.4, -0.12),
                  'ncol':           8}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'above')

    return fig, ax, add_legend

def two_v(style=[]):
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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(8.6, 16.8))
    grid = GridSpec(11, 1)
    ax = [fig.add_subplot(grid[:5, 0]), fig.add_subplot(grid[6:, 0])]

    plt.subplots_adjust(left=0.15, right=0.95,
                        bottom=0.06, top=0.98)

    names = [['top', 'in', 'inside'],
             ['bottom'],
             ['out', 'outside', 'above'],
             ['below']]
    locations = [{'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1),
                  'ncol':           4},
                 {'axes':           1,
                  'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.12),
                  'ncol':           4}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'above')

    return fig, ax, add_legend

def two_v_colourbars(style=[]):
    """A figure with two sets of axes vertically with colourbars.

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
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(10.7, 16.8))
    grid = GridSpec(11, 1)
    ax = [fig.add_subplot(grid[:5, 0]), fig.add_subplot(grid[6:, 0])]

    plt.subplots_adjust(left=0.15, right=0.95,
                        bottom=0.06, top=0.98)

    names = [['top', 'in', 'inside'],
             ['bottom'],
             ['out', 'outside', 'above'],
             ['below']]
    locations = [{'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1),
                  'ncol':           4},
                 {'axes':           1,
                  'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.12),
                  'ncol':           4}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'above')

    return fig, ax, add_legend

def three_h(style=[]):
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
        function
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(28.7, 8.3))
    grid = GridSpec(1, 38)
    ax = [fig.add_subplot(grid[0, :10]),
          fig.add_subplot(grid[0, 14:24]),
          fig.add_subplot(grid[0, 28:])]

    plt.subplots_adjust(left=0.06, right=0.98,
                        bottom=0.12, top=0.95)

    names = [['lefthand'],
             ['middle', 'centre', 'center'],
             ['righthand', 'in', 'inside'],
             ['right'],
             ['above', 'outside', 'out'],
             ['below']]
    locations = [{'axes':           2,
                  'loc':            'center left',
                  'bbox_to_anchor': (1.29, 0.5)},
                 {'axes':           1,
                  'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1),
                  'ncol':           8},
                 {'axes':           1,
                  'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.12),
                  'ncol':           8}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'above')

    return fig, ax, add_legend

def three_h_colourbars(style=[]):
    """A figure with three sets of axes horizontally and colourbars.

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
            function to add a pre-positioned legend.
    """

    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(36.2, 8.3))
    grid = GridSpec(1, 19)
    ax = [fig.add_subplot(grid[0, :5]),
          fig.add_subplot(grid[0, 7:12]),
          fig.add_subplot(grid[0, 14:])]

    plt.subplots_adjust(left=0.06, right=0.98,
                        bottom=0.12, top=0.95)

    names = [['lefthand'],
             ['middle', 'centre', 'center'],
             ['righthand', 'in', 'inside'],
             ['right'],
             ['above', 'outside', 'out'],
             ['below']]
    locations = [{'axes':           2,
                  'loc':            'center left',
                  'bbox_to_anchor': (1.29, 0.5)},
                 {'axes':           1,
                  'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1),
                  'ncol':           8},
                 {'axes':           1,
                  'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.12),
                  'ncol':           8}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'above')

    return fig, ax, add_legend

def three_square(style=[], blank=2):
    """Axes in a square with one missing.

    The legend is placed in the blank space by default.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp.
        blank : int, optional
            empty quadrant (numbered left to right then top to bottom).

    Returns
    -------

        figure
            figure.
        list
            axes.
        function
            function to add a pre-positioned legend.
    """

    assert blank in [1, 2, 3, 4]
    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(18.4, 16.6))

    grid = GridSpec(11, 7)
    ax = [[None if blank == 1 else fig.add_subplot(grid[:5, :3]),
           None if blank == 2 else fig.add_subplot(grid[:5, 4:])],
          [None if blank == 3 else fig.add_subplot(grid[6:, :3]),
           None if blank == 4 else fig.add_subplot(grid[6:, 4:])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97)

    names = [['topleft'],
             ['topright'],
             ['bottomleft'],
             ['bottomright'],
             ['blank', 'empty', 'out', 'outside', str(blank)]]
    axes = [1, 0, 3, 2]
    bbox = [(-0.91,0.5), (1.91,0.5), (-0.91, 0.5), (1.91,0.5)]
    locations = [{'axes':           axes[blank-1],
                  'loc':            'center',
                  'bbox_to_anchor': bbox[blank-1]}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, str(blank))

    return fig, ax, add_legend

def three_square_colourbars(style=[], blank=2):
    """Axes with colourbars in a square with one missing.

    The legend is placed in the blank space by default.

    Arguments
    ---------

        style : str or array, optional
            style sheet(s). Default: tp.
        blank : int, optional
            empty quadrant (numbered left to right then top to bottom).

    Returns
    -------

        figure
            figure.
        list
            axes.
        function
            function to add a pre-positioned legend.
    """

    assert blank in [1, 2, 3, 4]
    if isinstance(style, str): style=[style]
    default_style.extend(style)
    plt.style.use(default_style)
    fig = plt.figure(figsize=(23.4, 16.6))

    grid = GridSpec(11, 12)
    ax = [[None if blank == 1 else fig.add_subplot(grid[:5, :5]),
           None if blank == 2 else fig.add_subplot(grid[:5, 7:])],
          [None if blank == 3 else fig.add_subplot(grid[6:, :5]),
           None if blank == 4 else fig.add_subplot(grid[6:, 7:])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97)

    names = [['topleft'],
             ['topright'],
             ['bottomleft'],
             ['bottomright'],
             ['blank', 'empty', 'out', 'outside', str(blank)]]
    axes = [1, 0, 3, 2]
    bbox = [(-1.31, 0.5), (2.31,0.5), (-1.31, 0.5), (2.31,0.5)]
    locations = [{'axes':           axes[blank-1],
                  'loc':            'center',
                  'bbox_to_anchor': bbox[blank-1]}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, str(blank))

    return fig, ax, add_legend

def four_square(style=[]):
    """A figure with four sets of axes in a square.

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
    fig = plt.figure(figsize=(18.5, 16.6))

    grid = GridSpec(11,7)
    ax = [[fig.add_subplot(grid[:5, :3]), fig.add_subplot(grid[:5, 4:])],
          [fig.add_subplot(grid[6:, :3]), fig.add_subplot(grid[6:, 4:])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97)

    names = [['topleft'],
             ['topright', 'in', 'inside'],
             ['bottomleft'],
             ['bottomright'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'axes':           1,
                  'loc':            'center left',
                  'bbox_to_anchor': (1, -0.12)},
                 {'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (1.2, 1),
                  'ncol':           8},
                 {'axes':           2,
                  'loc':            'upper center',
                  'bbox_to_anchor': (1.2, -0.12),
                  'ncol':           8}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'above')

    return fig, ax, add_legend

def four_square_colourbars(style=[]):
    """A figure with four sets of axes with colourbars in a square.

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
    fig = plt.figure(figsize=(23.4, 16.6))

    grid = GridSpec(11,12)
    ax = [[fig.add_subplot(grid[:5, :5]), fig.add_subplot(grid[:5, 7:])],
          [fig.add_subplot(grid[6:, :5]), fig.add_subplot(grid[6:, 7:])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97)

    names = [['topleft'],
             ['topright', 'in', 'inside'],
             ['bottomleft'],
             ['bottomright'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'axes':           1,
                  'loc':            'center left',
                  'bbox_to_anchor': (1, -0.12)},
                 {'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (1.4, 1),
                  'ncol':           8},
                 {'axes':           2,
                  'loc':            'upper center',
                  'bbox_to_anchor': (1.4, -0.12),
                  'ncol':           8}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'above')

    return fig, ax, add_legend
