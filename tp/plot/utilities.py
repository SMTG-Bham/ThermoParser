"""Utilities to aid the plotting scripts."""

#Functions
#---------
#
#    colour_scale:
#        sorts colour limits and colourbar format.
#    parse_colours
#        returns a colour map from a range of input types.
#    scale_to_axis
#        scales data to axes limits
#    set_locators:
#        set locators one-liner.
#"""

import matplotlib as mpl
import numpy as np
import tp

def colour_scale(c, name, cmap, cmin=None, cmax=None, cscale=None,
                 unoccupied='grey'):
    """Formats the colour scale for phono3py quantities.

    Attempts to set the limits to maximise visibility of the most
    important data and determines if the colourbar should be extended.

    Arguments
    ---------

        c : array-like
            colour data.
        name : str
            name of colour variable.
        cmap : colormap
            colourmap.

        cmin : float, optional
            colour scale minimum. Default: display 99 % data.
        cmax : float, optional
            colour scale maximum. Default: display 99.9 % data.
        cscale : str, optional
            override colour scale (linear/ log). Default: None.
        unoccupied : str, optional
            if the colour variable is occupation, values below 1 are
            coloured in this colour. If set to None, or cmin is set,
            this feature is turned off. Default: grey.

    Returns
    -------

        cnorm
            colour normalisation object.
        extend
            which direction to extend the colourbar in.
    """

    extend = [False, False]
    if name in ['frequency', 'heat_capacity']:
        if cmin is None:
            cmin = np.nanmin(c)
        elif cmin > np.nanmin(c):
            extend[0] = True
        if cmax is None:
            cmax = np.nanmax(c)
        elif cmax < np.nanmax(c):
            extend[1] = True
        cnorm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)

    elif name == 'occupation':
        csort = np.ravel(np.ma.masked_invalid(np.ma.masked_equal(c, 0)).compressed())
        csort = csort[csort.argsort()]
        clim = csort[int(round(len(csort)*99.9/100 - 1, 0))]
        if cmin is None:
            cmin = np.nanmin(c)
            if unoccupied is not None and cmin < 1:
                cmin = 1
                extend[0] = True
                cmap.set_under(unoccupied)
        elif cmin > np.nanmin(c):
            extend[0] = True
        if cmax is None:
            cmax = clim
            extend[1] = True
        elif cmax < np.nanmax(c):
            extend[1] = True
        if cmax/cmin < 10:
            cnorm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
        else:
            cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    else:
        csort = np.ravel(np.ma.masked_invalid(np.ma.masked_equal(c, 0)).compressed())
        csort = csort[csort.argsort()]
        clim = [csort[int(round(len(csort)/100 - 1, 0))],
                csort[int(round(len(csort)*99/100 - 1, 0))]]
        if cmin is None:
            cmin = clim[0]
        elif cmin > clim[0]:
            extend[0] = True
        if cmax is None:
            cmax = clim[1]
        elif cmax < clim[1]:
            extend[1] = True
        cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    if cscale is not None:
        if cscale == 'linear':
            cnorm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
        elif cscale == 'log':
            cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    if extend[0] and extend[1]:
        extend = 'both'
    elif extend[0] and not extend[1]:
        extend = 'min'
    elif not extend[0] and extend[1]:
        extend = 'max'
    elif not extend[0] and not extend[1]:
        extend = 'neither'

    return cnorm, extend

def parse_colours(colour):
    """Parses a range of inputs into a colourmap

    Argument
    --------

        colour : colourmap or str or array-like or dict
            colourmap or colourmap name or #rrggbb highlight colour or
            highlight, min, max colours in that order, or dictionary
            with cmid and cmin and/or cmax keys.

    Returns
    -------

        colourmap
            colourmap
    """

    from copy import copy

    try:
        try:
            cmap = copy(mpl.cm.get_cmap(colour))
        except AttributeError:
            cmap = copy(mpl.colormaps[colour])
    except ValueError:
        if isinstance(colour, mpl.colors.ListedColormap):
            cmap = copy(colour)
        elif isinstance(colour, str):
            cmap = tp.plot.colour.uniform(colour)
        elif isinstance(colour, list):
            try: # rgb
                cmap = tp.plot.colour.uniform(colour)
            except ValueError: # hex or named
                cmap = tp.plot.colour.uniform(*colour)
        elif isinstance(colour, dict):
            cmap = tp.plot.colour.uniform(**colour)
        else:
            raise Exception('colour must be a colourmap, colourmap '
                            'name, single highlight colour or '
                            'highlight, min, max colours in that '
                            'order, or a dictionary with cmid and '
                            'cmin and/or cmax keys. Colour format must '
                            'be hex rgb (array) or a named colour '
                            'recognised by matplotlib. ')

    return cmap

def scale_to_axis(ax, data, exclude=[], scale=None, axis='y'):
    """Scale data to fit an axis.

    Assumes data is linear, but still works for linear and log axes.

    Arguments
    ---------

        ax : axes
            axes to fit to.
        data : array-like or dict
            data to scale.

        exclude : array-like, optional
            keys to exclude from scaling. Default: none.
        scale : array-like, optional
            force scaling limits to [min, max]. Default: axis limits.
        axis : bool, optional
            axis to scale to. Default: y.

    Returns
    -------

        data
            scaled data.
    """

    if scale is not None:
        assert isinstance(scale, list) and len(list(scale)) == 2, \
               'if specified, scale must be a two element array-like [min,max]'
    assert axis in ['x', 'y'], 'axis must be x or y'

    if not isinstance(data, dict):
        data = {'qwerfv': data}
    dmax = 0
    for key in data:
        if key not in exclude:
            ymax = float(np.amax(data[key]))
            if ymax > dmax: dmax = ymax

    axes = {'x': ax.get_xlim(), 'y': ax.get_ylim()}
    axscale = {'x': ax.get_xaxis().get_scale(), 'y': ax.get_yaxis().get_scale()}
    if scale is None:
        scale = axes[axis]
    log = axscale[axis] == 'log'
    if log:
        scale = np.log10(scale)

    for key in data:
        if key not in exclude:
            data[key] = np.multiply(data[key], (scale[1]-scale[0]) / dmax)
            data[key] = np.add(scale[0], data[key])
            if log:
                data[key] = np.power(10, data[key])

    if 'qwerfv' in data:
        data = data['qwerfv']

    return data

def set_locators(ax, x=None, y=None, dos=False):
    """Set locators quickly.

    If log, sets scale as well.

    Arguments
    ---------

        ax : axes
            axes to format
        x : str, optional
            locator on x axis. Accepts linear, log or null.
            Default: do nothing.
        y : str, optional
            locator on y axis. Accepts linear, log or null.
            Default: do nothing.
        dos : bool, optional
            removes axes ticks and ticklabels and y axis label. Runs
            first, so ticks can be reinstated on x, e.g. for waterfall
            plots using x='log'. Default: False.

    Returns
    -------

        None
            formats ax directly.
    """

    assert isinstance(dos, bool), 'dos must be True or False.'

    from tp.settings import locator
    if dos:
        ax.xaxis.set_major_locator(locator()['null'])
        ax.xaxis.set_minor_locator(locator()['null'])
        ax.yaxis.set_major_locator(locator()['null'])
        ax.yaxis.set_minor_locator(locator()['null'])
        ax.set_ylabel('')

    if x is not None and x != 'null':
        ax.set_xscale(x)
    if y is not None and y != 'null':
        ax.set_yscale(y)
    for a, scale in zip([ax.xaxis, ax.yaxis], [x, y]):
        if scale is not None:
            if scale == 'linear':
                a.set_major_locator(locator()['major'])
                a.set_minor_locator(locator()['minor'])
            elif scale == 'log':
                a.set_major_locator(locator()['log'])
            elif scale == 'null':
                a.set_major_locator(locator()['null'])
                a.set_minor_locator(locator()['null'])
            else:
                raise Exception('x and y must be "linear", "log" or '
                                '"null" if specified.')

    return
