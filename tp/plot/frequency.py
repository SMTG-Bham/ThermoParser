"""Functions which plot phonon frequency on the x-axis.

Functions:
    add_dos:
        phonon dos.
    add_cum_kappa:
        cumulative kappa vs frequency.
    add_waterfall:
        scatter plots of various values.
    add_projected_waterfall:
        waterfall, but with a second quantity projected.

    format_waterfall
        formatting for the waterfall plots
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tp
import warnings
from tp.data import resolve

def add_dos(ax, data, colour, total=False, main=True, invert=False,
            scale=False, fill=True, fillalpha=20, line=False, **kwargs):
    """Adds a phonon density of states (DoS) to a set of axes.

    Arguments:
        ax : axes
            axes to plot on.
        data : dict
            DoS data.
        colour : dict
            RGB colours.

        total : bool, optional
            plot total DoS. Default: False

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            plot on y axis. Default: False.
        scale : bool, optional
            if main, scale to percent. If not main, scale to axis
            limits. Default: False.

        fill : bool, optional
            fill below lines. Default: True.
        fillalpha : int, optional
            fill alpha in %. Ignored if alpha specified in colour.
            Default: 20.
        line : bool, optional
            plot lines. Default: False.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.fill_between.
            Defaults: {'rasterized': False}
    """

    # defaults

    defkwargs = {'rasterized': False}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # data scaling

    exclude = ['frequency', 'meta']
    if not total: exclude.append('total')
    if scale:
        axscale = [0, 100] if main else None
        axis = 'x' if invert else 'y'
        data = tp.plot.utilities.scale_to_axis(ax, data, exclude, axscale,axis)

    # colours

    if total and 'total' not in colour:
        colour['total'] = '#000000'
    fillcolour = {}
    for c in colour:
        if fill:
            fillcolour[c] = colour[c][:7] + str(fillalpha)
        else:
            colour[c][:7] + '00'
        if not line: colour[c] = fillcolour[c]

    # plotting

    if total:
        exclude.append('total')
        if invert:
            ax.fill_between(data['total'], data['frequency'], label='Total',
                            facecolor=fillcolour['total'],
                            edgecolor=colour['total'], **kwargs)
        else:
            ax.fill_between(data['frequency'], data['total'], label='Total',
                            facecolor=fillcolour['total'],
                            edgecolor=colour['total'], **kwargs)

    for key in data:
        if key not in exclude:
            if invert:
                ax.fill_between(data[key], data['frequency'], label=key,
                                facecolor=fillcolour[key],
                                edgecolor=colour[key], **kwargs)
            else:
                ax.fill_between(data['frequency'], data[key], label=key,
                                facecolor=fillcolour[key],
                                edgecolor=colour[key], **kwargs)

    # axes formatting

    if main:
        if invert:
            axlabels = tp.settings.inverted_labels()
            ax.set_ylabel(axlabels['frequency'])
            ax.set_xlabel(axlabels['dos'])
            ax.set_xlim(left=0)
            tp.plot.utilities.set_locators(ax, x='null', y='linear')
        else:
            axlabels = tp.settings.labels()
            ax.set_xlabel(axlabels['frequency'])
            ax.set_ylabel(axlabels['dos'])
            ax.set_ylim(bottom=0)
            tp.plot.utilities.set_locators(ax, y='null', x='linear')

    return

def add_cum_kappa(ax, data, temperature=300, direction='avg', main=False,
                  invert=False, scale=False, colour='#000000', fill=False,
                  fillcolour=20, line=True, **kwargs):
    """Cumulates and plots kappa against frequency.

    Arguments:
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
            set ticks, labels, limits. Default: False.
        invert : bool, optional
            plot on y axis. Default: False.
        scale : bool, optional
            if main, scale to percent. If not main, scale to axis
            limits. Default: True.

        colour : str, optional
            #RRGGBB line colour. Default: #000000.
        fill : bool, optional
            fill below lines. Default: False.
        fillcolour : int or str, optional
            if an integer from 0-99, and colour in #RRGGBB format,
            applies alpha in % to colour. Otherwise treats it as a
            colour. Default: 20.
        line : bool, optional
            plot lines. Default: True.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.fill_between
            if filled or matplotlib.pyplot.plot otherwise.
            Defaults: {'label':      '$\mathregular{\kappa_l}$',
                       'rasterized': False}
    """

    # defaults

    defkwargs = {'label':      '$\mathregular{\kappa_l}$',
                 'rasterized': False}

    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # data formatting and calculation

    data = tp.resolve.resolve(data, 'mode_kappa', temperature=temperature,
                                                  direction=direction)
    k = np.ravel(data['mode_kappa'])
    f = np.ravel(data['frequency'])

    f, k = tp.calculate.cumulate(f, k)
    np.savetxt('cumkappa-frequency-{:.0f}K-{}.dat'.format(temperature,
                                                          direction),
               np.transpose([f, k]), header='Frequency(THz) k_l(Wm-1K-1)')
    f = np.append(f, 100*f[-1])
    k = np.append(k, k[-1])

    if scale:
        axscale = [0, 100] if main else None
        axis = 'x' if invert else 'y'
        k = tp.plot.utilities.scale_to_axis(ax, k, scale=axscale, axis=axis)

    # colour

    if fill:
        if colour.startswith('#') and isinstance(fillcolour, int) and \
           len(colour) == 7:
            if fillcolour >= 0 and fillcolour <= 9:
                fillcolour = colour + '0' + str(fillcolour)
            elif fillcolour >= 10 and fillcolour <= 99:
                fillcolour = colour + str(fillcolour)
            else:
                raise Exception('Expected alpha value between 0 and 99')
        elif isinstance(fillcolour, int):
            warnings.warn('integer fill colours ignored when colour format is '
                          'not #RRGGBB.')
            fillcolour = colour
        if not line: colour = fillcolour

        # plotting

        if invert:
            ax.fill_between(k, f, facecolor=fillcolour, edgecolor=colour,
                            **kwargs)
        else:
            ax.fill_between(f, k, facecolor=fillcolour, edgecolor=colour,
                            **kwargs)
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
            ax.set_ylabel(axlabels['frequency'])
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

    Arguments:
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

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.scatter.
            Defaults: {'alpha':      0.3,
                       'linewidth':  0,
                       'marker':     '.',
                       'rasterized': True,
                       's':          1}
    """

    # defaults

    defkwargs = {'alpha':      0.3,
                 'linewidth':  0,
                 'marker':     '.',
                 'rasterized': True,
                 's':          1}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # data formatting

    tnames = tp.settings.to_tp()
    if invert: quantity, xquantity = xquantity, quantity
    quantity = tnames[quantity] if quantity in tnames else quantity
    xquantity = tnames[xquantity] if xquantity in tnames else xquantity
    data = tp.data.resolve.resolve(data, [quantity, xquantity],
                                   temperature, direction)
    x = np.ravel(data[xquantity])
    y = np.abs(np.ravel(data[quantity]))

    # colour

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

    Arguments:
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
        cscale : str, optional
            override colour scale (linear/ log). Default: None.
        unoccupied : str, optional
            if the colour variable is occupation, values below 1 are
            coloured in this colour. If set to None, or cmin is set,
            this feature is turned off. Default: grey.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.scatter.
            Defaults: {'alpha':      0.3,
                       'linewidth':  0,
                       'marker':     '.',
                       'rasterized': True,
                       's':          1}

    Returns:
        colorbar
            colourbar for projected data.
    """

    from copy import copy

    # defaults

    defkwargs = {'alpha':      0.3,
                 'linewidth':  0,
                 'marker':     '.',
                 'rasterized': True,
                 's':          1}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # data formatting

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
        else:
            axlabels = tp.settings.labels()
        ax.set_xlabel(axlabels[xquantity])
        ax.set_ylabel(axlabels[quantity])

    return cbar

def format_waterfall(ax, data, yquantity, xquantity='frequency',
                     temperature=None, direction=None):
    """Formats axes for waterfall plots.

    Arguments:
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
    ax.plot(lim['x'], lim['y'], alpha=0)
    tp.plot.utilities.set_locators(ax, x=loc['x'], y=loc['y'])

    return
