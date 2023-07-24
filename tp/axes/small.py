"""Provides pre-sized figures for paper-sized figures.

Each function returns a figure, a set of axes or an array of sets of
axes and an add_legend function. The legend function comes with a
choice of positions and takes normal ax.legend arguments.

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

default_style = tp.settings.style()

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
    fig, ax = plt.subplots(figsize=(8.6/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.15, right=0.95,
                        bottom=0.12, top=0.95)

    names = [['in', 'inside'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'loc':            'center left',
                  'bbox_to_anchor': (1, 0.5)},
                 {'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1),
                  'ncol':           4},
                 {'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.1),
                  'ncol':           4}]
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
    fig, ax = plt.subplots(figsize=(10.5/2.54, 8.3/2.54))

    plt.subplots_adjust(left=0.14, right=0.96,
                        bottom=0.12, top=0.95)

    names = [['in', 'inside'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'loc':            'center left',
                  'bbox_to_anchor': (1.27, 0.5)},
                 {'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1),
                  'ncol':           4},
                 {'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.1),
                  'ncol':           4}]
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
    fig = plt.figure(figsize=(11.05/2.54, 8.3/2.54))
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
                  'bbox_to_anchor': (0.5, 1),
                  'ncol':           4},
                 {'axes':           0,
                  'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.1),
                  'ncol':           4}]
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
    fig = plt.figure(figsize=(12.2/2.54, 8.3/2.54))
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
                  'bbox_to_anchor': (1.6, 0.5)},
                 {'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (0.5, 1),
                  'ncol':           4},
                 {'axes':           0,
                  'loc':            'upper center',
                  'bbox_to_anchor': (0.5, -0.1),
                  'ncol':           4}]
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
    fig = plt.figure(figsize=(17.6/2.54, 8.3/2.54))
    grid = GridSpec(1, 9)
    ax = [fig.add_subplot(grid[0, :4]), fig.add_subplot(grid[0, 5:])]

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
                  'bbox_to_anchor': (1.15, 1),
                  'ncol':           8},
                 {'axes':           0,
                  'loc':            'upper center',
                  'bbox_to_anchor': (1.15, -0.1),
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
    fig = plt.figure(figsize=(22/2.54, 8.3/2.54))
    grid = GridSpec(1, 9)
    ax = [fig.add_subplot(grid[0, :4]), fig.add_subplot(grid[0, 5:])]

    plt.subplots_adjust(left=0.08, right=0.98,
                        bottom=0.12, top=0.95)

    names = [['lefthand'],
             ['righthand', 'in', 'inside'],
             ['out', 'outside', 'right'],
             ['above'],
             ['below']]
    locations = [{'axes':           1,
                  'loc':            'center left',
                  'bbox_to_anchor': (1.27, 0.5)},
                 {'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (1.4, 1),
                  'ncol':           8},
                 {'axes':           0,
                  'loc':            'upper center',
                  'bbox_to_anchor': (1.4, -0.1),
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
    fig = plt.figure(figsize=(8.9/2.54, 16.8/2.54))
    grid = GridSpec(15, 1)
    ax = [fig.add_subplot(grid[:7, 0]), fig.add_subplot(grid[8:, 0])]

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
                  'bbox_to_anchor': (0.5, -0.1),
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
    fig = plt.figure(figsize=(11.1/2.54, 16.8/2.54))
    grid = GridSpec(15, 1)
    ax = [fig.add_subplot(grid[:7, 0]), fig.add_subplot(grid[8:, 0])]

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
                  'bbox_to_anchor': (0.5, -0.1),
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
    fig = plt.figure(figsize=(27/2.54, 8.3/2.54))
    grid = GridSpec(1, 14)
    ax = [fig.add_subplot(grid[0, :4]),
          fig.add_subplot(grid[0, 5:9]),
          fig.add_subplot(grid[0, 10:])]

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
                  'bbox_to_anchor': (0.5, -0.1),
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
    fig = plt.figure(figsize=(33.8/2.54, 8.3/2.54))
    grid = GridSpec(1, 14)
    ax = [fig.add_subplot(grid[0, :4]),
          fig.add_subplot(grid[0, 5:9]),
          fig.add_subplot(grid[0, 10:])]

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
                  'bbox_to_anchor': (0.5, -0.1),
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
    fig = plt.figure(figsize=(18.3/2.54, 16.6/2.54))

    grid = GridSpec(15, 9)
    ax = [[None if blank == 1 else fig.add_subplot(grid[:7, :4]),
           None if blank == 2 else fig.add_subplot(grid[:7, 5:])],
          [None if blank == 3 else fig.add_subplot(grid[8:, :4]),
           None if blank == 4 else fig.add_subplot(grid[8:, 5:])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97)

    names = [['topleft'],
             ['topright'],
             ['bottomleft'],
             ['bottomright'],
             ['blank', 'empty', 'out', 'outside', str(blank)]]
    axes = [1, 0, 3, 2]
    bbox = [(-0.8,0.5), (1.8,0.5), (-0.8, 0.5), (1.8,0.5)]
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
    fig = plt.figure(figsize=(22.7/2.54, 16.6/2.54))

    grid = GridSpec(15, 9)
    ax = [[None if blank == 1 else fig.add_subplot(grid[:7, :4]),
           None if blank == 2 else fig.add_subplot(grid[:7, 5:])],
          [None if blank == 3 else fig.add_subplot(grid[8:, :4]),
           None if blank == 4 else fig.add_subplot(grid[8:, 5:])]]

    plt.subplots_adjust(left=0.1, right=0.97,
                        bottom=0.07, top=0.97)

    names = [['topleft'],
             ['topright'],
             ['bottomleft'],
             ['bottomright'],
             ['blank', 'empty', 'out', 'outside', str(blank)]]
    axes = [1, 0, 3, 2]
    bbox = [(-1.13,0.5), (2.13,0.5), (-1.13, 0.5), (2.13,0.5)]
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
    fig = plt.figure(figsize=(18.2/2.54, 16.6/2.54))

    grid = GridSpec(15,9)
    ax = [[fig.add_subplot(grid[:7, :4]), fig.add_subplot(grid[:7, 5:])],
          [fig.add_subplot(grid[8:, :4]), fig.add_subplot(grid[8:, 5:])]]

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
                  'bbox_to_anchor': (1, -0.11)},
                 {'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (1.15, 1),
                  'ncol':           8},
                 {'axes':           2,
                  'loc':            'upper center',
                  'bbox_to_anchor': (1.15, -0.1),
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
    fig = plt.figure(figsize=(22.7/2.54, 16.6/2.54))

    grid = GridSpec(15,9)
    ax = [[fig.add_subplot(grid[:7, :4]), fig.add_subplot(grid[:7, 5:])],
          [fig.add_subplot(grid[8:, :4]), fig.add_subplot(grid[8:, 5:])]]

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
                  'bbox_to_anchor': (1, -0.11)},
                 {'axes':           0,
                  'loc':            'lower center',
                  'bbox_to_anchor': (1.35, 1),
                  'ncol':           8},
                 {'axes':           2,
                  'loc':            'upper center',
                  'bbox_to_anchor': (1.35, -0.1),
                  'ncol':           8}]
    add_legend = tp.axes.legend.add_add_legend(ax, locations, names, 'above')

    return fig, ax, add_legend
