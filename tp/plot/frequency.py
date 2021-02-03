"""Functions which plot phonon frequency on the x-axis.

Each function takes a set of axes and data dictionary as its main input,
and and an ``invert`` argument to plot sideways by a phonon dispersion,
and a main argument, which determines if axes limits are set and
sometimes how scaling is handled, which helps if plotting multiple
quantities on the same axes.

Functions
---------

    add_dos:
        phonon dos.
    add_cum_kappa:
        cumulative kappa vs frequency.
    add_waterfall:
        scatter plots of various values.
    add_projected_waterfall:
        waterfall, but with a second quantity projected.
    add_density:
        density of phonon modes for a property vs frequency.

    format_waterfall
        formatting for the waterfall and density plots
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import tp
import warnings
import yaml
from scipy.stats import gaussian_kde
from tp.data import resolve

warnings.filterwarnings('ignore', module='matplotlib')

try:
    filename = '{}/.config/tprc.yaml'.format(os.path.expanduser("~"))
    with open(filename, 'r') as f:
        conf = yaml.safe_load(f)
except Exception:
    conf = None

def add_dos(ax, data, total=False, main=True, invert=False, scale=False,
            colour='tab10', fill=True, fillalpha=0.2, line=True, **kwargs):
    """Adds a phonon density of states (DoS) to a set of axes.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            DoS data.

        total : bool, optional
            plot total DoS. Default: False

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            plot on y axis. Default: False.
        scale : bool, optional
            if main, scale to percent. If not main, scale to axis
            limits. Default: False.

        colour : dict or list or str or colourmap, optional
            RGB colours per atom as a dictionary or a list in POSCAR
            order. Can instead provide a colourmap or colourmap name.
            Default: tab10.
        fill : bool, optional
            fill below lines. Default: True.
        fillalpha : float, optional
            fill alpha scaled to 0-1. Default: 0.2.
        line : bool, optional
            plot lines. Default: True.

        kwargs
            keyword arguments passed to matplotlib.pyplot.fill_between.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                rasterized: False

    Returns
    -------

        None
            adds plot directly to ax.
    """

    # defaults

    defkwargs = {'rasterized': False}

    if conf is None or 'dos_kwargs' not in conf or conf['dos_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['dos_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['total', 'main', 'invert', 'scale', 'fill', 'line'],
                           [ total,   main,   invert,   scale,   fill,   line]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)
    assert fill or line, 'fill or line or both must be True.'
    assert isinstance(fillalpha, (float, int)) and fillalpha >= 0 \
                                               and fillalpha <= 1, \
           'fillalpha must be a float/ integer between 0 and 1.'

    # data scaling

    data = dict(data)
    f = data['frequency']
    del data['frequency']
    if 'meta' in data: del data['meta']
    if not total:
        del data['total']

    if scale:
        axscale = [0, 100] if main else None
        axis = 'x' if invert else 'y'
        data = tp.plot.utilities.scale_to_axis(ax, data, scale=axscale,
                                               axis=axis)
    if total:
        totaldata = data['total']
        del data['total']

    # colours
    # Tries to read the colour as a colourmap name, then colourmap
    # object, then list of colours, then dictionary, and converts it
    # into a dictionary if necessary.

    if not isinstance(colour, dict):
        if not isinstance(colour, list):
            try:
                cmap = plt.cm.get_cmap(colour)(np.linspace(0,1,len(data)))
            except Exception:
                cmap = colour(np.linspace(0,1,len(data)))
        else:
            cmap = colour
        colour = {}
        for i, c in enumerate(data):
            colour[c] = cmap[i]

    if total and 'total' not in colour:
        colour['total'] = '#000000'
    fillcolour = {}
    for c in colour:
        if isinstance(colour[c], str):
            colour[c] = tp.plot.colour.rgb2array(colour[c])
        if fill:
            fillcolour[c] = list(colour[c])
            fillcolour[c][3] = fillalpha

    # plotting

    if total:
        if invert:
            if fill:
                ax.fill_between(totaldata, f, facecolor=fillcolour['total'], linewidth=0)
            if line:
                ax.plot(totaldata, f, label='Total', color=colour['total'], **kwargs)
        else:
            if fill:
                ax.fill_between(f, totaldata, facecolor=fillcolour['total'], linewidth=0)
            if line:
                ax.plot(f, totaldata, label='Total', color=colour['total'], **kwargs)

    for key in data:
        if invert:
            if fill:
                ax.fill_between(data[key], f, facecolor=fillcolour[key], linewidth=0)
            if line:
                ax.plot(data[key], f, label=key, color=colour[key], **kwargs)
        else:
            if fill:
                ax.fill_between(f, data[key], facecolor=fillcolour[key], linewidth=0)
            if line:
                ax.plot(f, data[key], label=key, color=colour[key], **kwargs)

    # axes formatting

    if main:
        if invert:
            axlabels = tp.settings.inverted_labels()
            pad = mpl.rcParams['xtick.major.pad']
            ax.set_xlabel(axlabels['dos'], labelpad=pad)
            ax.tick_params(axis='y', labelleft=False)
            ax.set_xlim(left=0)
            tp.plot.utilities.set_locators(ax, x='null', y='linear')
        else:
            axlabels = tp.settings.labels()
            pad = mpl.rcParams['ytick.major.pad']
            ax.set_xlabel(axlabels['frequency'])
            ax.set_ylabel(axlabels['dos'], labelpad=pad)
            ax.set_ylim(bottom=0)
            tp.plot.utilities.set_locators(ax, y='null', x='linear')

    return

def add_cum_kappa(ax, data, temperature=300, direction='avg', main=True,
                  invert=False, scale=False, colour='#000000', fill=False,
                  fillcolour=0.2, line=True, **kwargs):
    """Cumulates and plots kappa against frequency.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            Phono3py data including mode_kappa and frequency.

        temperature : float, optional
            temperature in K (finds nearest). Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            plot on y axis. Default: False.
        scale : bool, optional
            if main, scale to percent. If not main, scale to axis
            limits. Default: True.

        colour : str, optional
            RGB line colour. Default: #000000.
        fill : bool, optional
            fill below lines. Default: False.
        fillcolour : int or str, optional
            if a float from 0-1 and colour in #RRGGBB format, sets
            fill colour opacity. Otherwise treats it as a colour.
            Default: 0.2.
        line : bool, optional
            plot lines. Default: True.

        kwargs
            keyword arguments passed to matplotlib.pyplot.fill_between
            if filled or matplotlib.pyplot.plot otherwise.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                label:      $\mathregular{\kappa_l}$
                rasterized: False

    Returns
    -------

        None
            adds plot directly to ax.
    """

    # defaults

    defkwargs = {'label':      '$\mathregular{\kappa_l}$',
                 'rasterized': False}

    if conf is None or 'frequency_cum_kappa_kwargs' not in conf or \
       conf['frequency_cum_kappa_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['frequency_cum_kappa_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['main', 'invert', 'scale', 'fill', 'line'],
                           [ main,   invert,   scale,   fill,   line]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)
    assert fill or line, 'fill or line or both must be True.'

    # data formatting and calculation

    data = tp.resolve.resolve(data, 'mode_kappa', temperature=temperature,
                                                  direction=direction)
    k = np.ravel(data['mode_kappa'])
    f = np.ravel(data['frequency'])

    f, k = tp.calculate.cumulate(f, k)
    np.savetxt('cumkappa-frequency-{:.0f}K-{}.dat'.format(
               data['meta']['temperature'], direction),
               np.transpose([f, k]), header='Frequency(THz) k_l(Wm-1K-1)')
    f = np.append(f, 100*f[-1])
    k = np.append(k, k[-1])

    if scale:
        axscale = [0, 100] if main else None
        axis = 'x' if invert else 'y'
        k = tp.plot.utilities.scale_to_axis(ax, k, scale=axscale, axis=axis)

    # colour
    # Tries to read the colour as an rgb code, then alpha value.

    if fill:
        try:
            fillcolour2 = tp.plot.colour.rgb2array(colour, fillcolour)
        except Exception:
            if isinstance(colour, list) and \
               isinstance(fillcolour, (float, int)) and fillcolour >= 0 and \
               fillcolour <= 1:
                fillcolour2 = list(colour)
                if len(colour) == 3:
                    fillcolour2.append(fillcolour)
                elif len(colour) == 4:
                    fillcolour2[3] = fillcolour
            else:
                fillcolour2 = fillcolour
        if not line: colour = fillcolour2

        # plotting

        if invert:
            ax.plot(k, f, color=colour, **kwargs)
            ax.fill_between(k, f, facecolor=fillcolour2)
        else:
            ax.plot(f, k, color=colour, **kwargs)
            ax.fill_between(f, k, facecolor=fillcolour2)
    else:
        if invert:
            ax.plot(k, f, color=colour, **kwargs)
        else:
            ax.plot(f, k, color=colour, **kwargs)

    # axes formatting

    if main:
        if invert:
            axlabels = tp.settings.inverted_labels()
            ax.set_xlabel(axlabels['cumulative_kappa'])
            ax.tick_params(axis='y', labelleft=False)
            ax.set_xlim(0, k[-2])
            ax.set_ylim(0, f[-2])

        else:
            axlabels = tp.settings.labels()
            ax.set_ylabel(axlabels['cumulative_kappa'])
            ax.set_xlabel(axlabels['frequency'])
            ax.set_ylim(0, k[-2])
            ax.set_xlim(0, f[-2])

        tp.plot.utilities.set_locators(ax, x='linear', y='linear')

    return

def add_waterfall(ax, data, quantity, xquantity='frequency', temperature=300,
                  direction='avg', main=True, invert=False, colour='viridis',
                  **kwargs):
    """Adds a waterfall plot of quantities against frequency.

    Has an option to change the x-quantity.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            data including frequency and quantity.
        quantity : str
            y-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path or
            mode_kappa.

        xquantity : str, optional
            x-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path or
            mode_kappa. Default: frequency.
        temperature : float, optional
            temperature in K. Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            invert x- and y-axes. Default: False.

        colour : colourmap or str or array-like, optional
            colourmap or colourmap name or list of colours (one for
            each band or one for each point) or a single colour.
            Default: viridis.

        kwargs
            keyword arguments passed to matplotlib.pyplot.scatter.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                alpha:      0.3
                s:          3
                edgecolors: black
                linewidth:  0
                marker:     '.'
                rasterized: True

    Returns
    -------

        None
            adds plot directly to ax.
    """

    # defaults

    defkwargs = {'alpha':      0.3,
                 's':          3,
                 'edgecolors': 'black',
                 'linewidth':  0,
                 'marker':     '.',
                 'rasterized': True}

    if conf is None or 'waterfall_kwargs' not in conf or \
       conf['waterfall_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['waterfall_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['main', 'invert'],
                           [ main,   invert]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)

    # data formatting

    if quantity == 'kappa': quantity = 'mode_kappa'
    if xquantity == 'kappa': xquantity = 'mode_kappa'
    tnames = tp.settings.to_tp()
    if invert: quantity, xquantity = xquantity, quantity
    quantity = tnames[quantity] if quantity in tnames else quantity
    xquantity = tnames[xquantity] if xquantity in tnames else xquantity

    data = tp.data.resolve.resolve(data, [quantity, xquantity],
                                   temperature, direction)
    x = np.ravel(data[xquantity])
    y = np.abs(np.ravel(data[quantity]))

    # colour
    # Tries to read as a colourmap name or colourmap object or list of
    # colours (of varying formats), one per band, and assigns
    # appropriately. Otherwise leaves as is, which is appropriate for 
    # a single colour or colour per point.

    s = np.shape(data[xquantity])
    try:
        colour = mpl.cm.get_cmap(colour)
        colours = [colour(i) for i in np.linspace(0, 1, s[1])]
        colours = np.tile(colours, (s[0], 1))
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            colours = [colour(i) for i in np.linspace(0, 1, s[1])]
            colours = np.tile(colours, (s[0], 1))
        elif isinstance(colour, str) and colour == 'skelton':
            colour = tp.plot.colour.skelton()
            colours = [colour(i) for i in np.linspace(0, 1, s[1])]
            colours = np.tile(colours, (s[0], 1))
        elif (isinstance(colour, list) or isinstance(colour, np.ndarray) and
              len(colour) == s[1]):
            if np.ndim(colour) == 1:
                colours = np.tile(colour, s[0])
            elif np.ndim(colour) == 2:
                colours = np.tile(colour, (s[0], 1))
        else:
            colours = colour

    # plotting

    ax.scatter(x, y, facecolor=colours, **kwargs)

    # axes formatting

    if main:
        format_waterfall(ax, {xquantity: x, quantity: y}, quantity, xquantity)
        if invert:
            axlabels = tp.settings.inverted_labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.tick_params(axis='y', labelleft=False)
        else:
            axlabels = tp.settings.labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.set_ylabel(axlabels[quantity])

    return

def add_projected_waterfall(ax, data, quantity, projected,
                            xquantity='frequency', temperature=300,
                            direction='avg', main=True, invert=False,
                            colour='viridis', cmin=None, cmax=None,
                            cscale=None, unoccupied='grey', **kwargs):
    """Adds a waterfall plot against frequency with a colour axis.

    Has an option to change the x-quantity.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            data including frequency and quantity.
        quantity : str
            y-axis quantity. Accepts gamma, group_velocity, gv_by_gv,
            heat_capacity, lifetime, mean_free_path or mode_kappa.
        projected : str
            colour-axis quantity. Accepts gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path,
            mode_kappa or occupation.

        xquantity : str, optional
            x-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path or
            mode_kappa. Default: frequency.
        temperature : float, optional
            temperature in K. Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            plot frequency on y-axis. Default: False.

        colour : colormap or str, optional
            colourmap or colourmap name. Default: viridis.
        cmin : float, optional
            colour scale minimum. Default: display 99 % data.
        cmax : float, optional
            colour scale maximum. Default: display 99.9 % data.
        cscale : str optional
            override colour scale (linear/ log). Default: None.
        unoccupied : str or array-like, optional
            if the colour variable is occupation, values below 1 are
            coloured in this colour. If set to None, or cmin is set,
            this feature is turned off. Default: grey.

        kwargs
            keyword arguments passed to matplotlib.pyplot.scatter.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                alpha:      0.3
                s:          3
                edgecolors: black
                linewidth:  0
                marker:     .
                rasterized: True

    Returns
    -------

        colorbar
            colourbar for projected data.
    """

    from copy import copy

    # defaults

    defkwargs = {'alpha':      0.3,
                 's':          3,
                 'edgecolors': 'black',
                 'linewidth':  0,
                 'marker':     '.',
                 'rasterized': True}

    if conf is None or 'projected_waterall_kwargs' not in conf or \
       conf['projected_waterfall_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['projected_waterfall_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['main', 'invert'],
                           [ main,   invert]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)

    # data formatting

    if quantity == 'kappa': quantity = 'mode_kappa'
    if xquantity == 'kappa': xquantity = 'mode_kappa'
    if projected == 'kappa': projected = 'mode_kappa'
    tnames = tp.settings.to_tp()
    axlabels = tp.settings.labels()
    if invert: quantity, xquantity = xquantity, quantity
    quantity = tnames[quantity] if quantity in tnames else quantity
    xquantity = tnames[xquantity] if xquantity in tnames else xquantity
    projected = tnames[projected] if projected in tnames else projected

    data = tp.data.resolve.resolve(data, [quantity, xquantity, projected],
                                   temperature, direction)
    x = np.ravel(data[xquantity])
    y = np.abs(np.ravel(data[quantity]))
    c = np.abs(np.ravel(data[projected]))

    # colour
    # Reads a colourmap or colourmap name.

    try:
        cmap = copy(mpl.cm.get_cmap(colour))
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            cmap = copy(colour)
        elif isinstance(colour, str) and colour == 'skelton':
            cmap = tp.plot.colour.skelton()
        else:
            raise Exception('Unrecognised colour argument. '
                            'Expected a colourmap or colourmap name.')

    cnorm, extend = tp.plot.utilities.colour_scale(c, projected, cmap, cmin,
                                                   cmax, cscale, unoccupied)

    # plotting

    scat = ax.scatter(x, y, c=c, norm=cnorm, cmap=cmap, **kwargs)

    # axes formatting

    cbar = plt.colorbar(scat, extend=extend)
    cbar.set_alpha(1)
    cbar.set_label(axlabels[projected])
    tp.plot.utilities.set_locators(cbar.ax, y=cbar.ax.yaxis.get_scale())
    cbar.draw_all()

    if main:
        format_waterfall(ax, {xquantity: x, quantity: y}, quantity, xquantity)
        if invert:
            axlabels = tp.settings.inverted_labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.tick_params(axis='y', labelleft=False)
        else:
            axlabels = tp.settings.labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.set_ylabel(axlabels[quantity])

    return cbar

def add_density(ax, data, quantity, xquantity='frequency', temperature=300,
                  direction='avg', main=True, invert=False, colour='Blues',
                  **kwargs):
    """Adds a density plot of quantities against frequency.

    Has an option to change the x-quantity.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            data including frequency and quantity.
        quantity : str
            y-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path or
            mode_kappa.

        xquantity : str, optional
            x-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path or
            mode_kappa. Default: frequency.
        temperature : float, optional
            temperature in K. Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            invert x- and y-axes. Default: False.

        colour : colourmap or str or array-like, optional
            colourmap or colourmap name. A single #rrggbb colour can be
            given to generate a custom uniform colourmap.
            Default: Blues.

        kwargs
            keyword arguments passed to matplotlib.pyplot.scatter.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                s:          1  
                rasterized: True

    Returns
    -------

        None
            adds plot directly to ax.
    """

    # defaults

    defkwargs = {'s':      1,
                 'rasterized': True}

    if conf is None or 'density_kwargs' not in conf or \
       conf['density_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['density_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['main', 'invert'],
                           [ main,   invert]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)

    # data formatting

    if quantity == 'kappa': quantity = 'mode_kappa'
    if xquantity == 'kappa': xquantity = 'mode_kappa'
    tnames = tp.settings.to_tp()
    if invert: quantity, xquantity = xquantity, quantity
    quantity = tnames[quantity] if quantity in tnames else quantity
    xquantity = tnames[xquantity] if xquantity in tnames else xquantity

    data = tp.data.resolve.resolve(data, [quantity, xquantity],
                                   temperature, direction)
    x = np.ravel(data[xquantity])
    y = np.abs(np.ravel(data[quantity]))
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x_dens, y_dens, z_dens = x[idx], y[idx], z[idx]

    # colour
    # Tries to read as a colourmap or colourmap name, or uses a single
    # #rrggbb colour as the highlight colour for a tp.plot.colour.uniform.

    try:
        colours = mpl.cm.get_cmap(colour)
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            colours = colour
        else:
            try:
                colours = tp.plot.colour.uniform(colour)
            except Exception:
                raise Exception('colour must be a colourmap, colourmap'
                                'name or single #rrggbb colour.')

    # plotting

    ax.scatter(x_dens, y_dens, c=z_dens, cmap=colour, **kwargs)

    # axes formatting

    if main:
        format_waterfall(ax, {xquantity: x_dens, quantity: y_dens}, quantity, xquantity)
        if invert:
            axlabels = tp.settings.inverted_labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.tick_params(axis='y', labelleft=False)
        else:
            axlabels = tp.settings.labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.set_ylabel(axlabels[quantity])

    return



def format_waterfall(ax, data, yquantity, xquantity='frequency',
                     temperature=None, direction=None):
    """Formats axes for waterfall plots.

    Arguments
    ---------

        ax : axes
            axes to format.
        data : array-like
            data.
        yquantity : str
            y quantity name.

        xquantity : str, optional
            x quantity name. Default: frequency.
        temperature : float, optional
            temperature in K. Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.

    Returns
    -------

        None
            formats ax directly.
    """

    data = tp.data.resolve.resolve(data, [xquantity, yquantity], temperature,
                                   direction)
    data2 = {'x': np.ravel(data[xquantity]), 'y': np.ravel(data[yquantity])}

    limit = {'x': ax.set_xlim,   'y': ax.set_ylim}
    scale = {'x': ax.set_xscale, 'y': ax.set_yscale}
    loc = {}
    lim = {}

    for axis, quantity in zip(['x', 'y'], [xquantity, yquantity]):
        if quantity in ['frequency', 'heat_capacity']:
            lim[axis] = [np.amin(data2[axis]), np.amax(data2[axis])]
            limit[axis](lim[axis][0], lim[axis][1])
            loc[axis] = 'linear'
        else:
            data2[axis] = np.abs(data2[axis])
            sort = np.ma.masked_invalid(np.ma.masked_equal(
                                        data2[axis], 0)).compressed()
            sort = sort[sort.argsort()]
            lim[axis] = [sort[int(round(len(sort)/100, 0))],
                         sort[int(round(len(sort)*99.9/100, 0))]]
            limit[axis](lim[axis][0], lim[axis][1])
            loc[axis] = 'log'

    # dummy data to enable axis formatting before plotting
    ax.scatter(lim['x'], lim['y'], alpha=0)
    tp.plot.utilities.set_locators(ax, x=loc['x'], y=loc['y'])

    return
