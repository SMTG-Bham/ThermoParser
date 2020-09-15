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
    scale_to_axes
        scales data to axes limits
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tp
from tp.data import resolve

def add_dos(ax, data, colour, total=False, main=False, invert=False,
            scale=True, fill=True, fillalpha=20, line=False, rasterise=False,
            **kwargs):
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
            set ticks, labels, limits. Default: False.
        invert : bool, optional
            plot on y axis. Default: False.
        scale : bool, optional
            if main, scale to percent. If not main, scale to axis
            limits. Default: True.

        fill : bool, optional
            fill below lines. Default: True.
        fillalpha : int, optional
            fill alpha in %. Ignored if alpha specified in colour.
            Default: 20.
        line : bool, optional
            plot lines. Default: False.
        rasterise : bool, optional
            rasterise plot. Default: False.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.fill_between.

    Returns:
        axes
            axes with DoS.
    """

    # data scaling

    exempt = ['frequency', 'meta']
    dmax = 0
    for key in data:
        if key not in exempt:
            ymax = float(np.amax(data[key]))
            if ymax > dmax: dmax = ymax

    if scale: data = scale_to_axes(ax, data, dmax, main, invert)

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
        if invert:
            ax.fill_between(data['total'], data['frequency'], label='Total',
                            facecolor=fillcolour['total'],
                            edgecolor=colour['total'], rasterized=rasterise,
                            **kwargs)
        else:
            ax.fill_between(data['frequency'], data['total'], label='Total',
                            facecolor=fillcolour['total'],
                            edgecolor=colour['total'], rasterized=rasterise,
                            **kwargs)

    exempt = ['frequency', 'meta', 'total']
    for key in data:
        if key not in exempt:
            if invert:
                ax.fill_between(data[key], data['frequency'], label=key,
                                facecolor=fillcolour[key],
                                edgecolor=colour[key], rasterized=rasterise,
                                **kwargs)
            else:
                ax.fill_between(data['frequency'], data[key], label=key,
                                facecolor=fillcolour[key],
                                edgecolor=colour[key], rasterized=rasterise,
                                **kwargs)

    # axes formatting

    if main:
        axlabels = tp.settings.labels()
        if invert:
            ax.set_ylabel(axlabels['frequency'])
            ax.set_xlabel(axlabels['dos'])

            ax.set_xlim(left=0, right=dmax)

            ax.tick_params(axis='x', which='both', top=False, bottom=False)
            ax.yaxis.set_major_locator(tp.settings.locator()['major'])
            ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])
        else:
            ax.set_xlabel(axlabels['frequency'])
            ax.set_ylabel(axlabels['dos'])

            ax.set_ylim(bottom=0, top=dmax)

            ax.tick_params(axis='y', which='both', top=False, bottom=False)
            ax.xaxis.set_major_locator(tp.settings.locator()['major'])
            ax.xaxis.set_minor_locator(tp.settings.locator()['minor'])

    return ax

def add_cum_kappa(ax, data, temperature=300, direction='avg', legend='\kappa_l',
                  main=False, invert=False, scale=False, colour='#000000',
                  fill=False, fillalpha=20, line=True, rasterise=False,
                  **kwargs):
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
        legend : str, optional
            legend entry, accepts maths. Default: \kappa_l.

        main : bool, optional
            set ticks, labels, limits. Default: False.
        invert : bool, optional
            plot on y axis. Default: False.
        scale : bool, optional
            if main, scale to percent. If not main, scale to axis
            limits. Default: True.

        colour : str, optional
            RGB line colour. Default: #000000.
        fill : bool, optional
            fill below lines. Default: False.
        fillalpha : int, optional
            fill alpha in %. Default: 20.
        line : bool, optional
            plot lines. Default: True.
        rasterise : bool, optional
            rasterise plot. Default: False.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.fill_between.

    Returns:
        axes
            axes with cumulated kappa.
    """

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

    if scale: k = scale_to_axes(ax, k, k[-1], main, invert)

    # colour

    fillcolour = colour + str(fillalpha) if fill else colour + '00'
    if not line: colour = fillcolour

    # plotting

    if invert:
        ax.fill_between(k, f, label='$\mathregular{{{}}}$'.format(legend),
                        facecolor=fillcolour, edgecolor=colour,
                        rasterized=rasterise, **kwargs)
    else:
        ax.fill_between(f, k, label='$\mathregular{{{}}}$'.format(legend),
                        facecolor=fillcolour, edgecolor=colour,
                        rasterized=rasterise, **kwargs)

    # axes formatting

    if main:
        axlabels = tp.settings.labels()
        if invert:
            ax.set_xlabel(axlabels['cumulative_kappa'])
            ax.set_ylabel(axlabels['frequency'])

            ax.set_xlim(0, k[-2])
            ax.set_ylim(0, f[-2])

        else:
            ax.set_ylabel(axlabels['cumulative_kappa'])
            ax.set_xlabel(axlabels['frequency'])

            ax.set_ylim(0, k[-2])
            ax.set_xlim(0, f[-2])

        ax.xaxis.set_major_locator(tp.settings.locator()['major'])
        ax.xaxis.set_minor_locator(tp.settings.locator()['minor'])
        ax.yaxis.set_major_locator(tp.settings.locator()['major'])
        ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])

    return ax

def add_waterfall(ax, data, quantity, xquantity='frequency', temperature=300,
                  direction='avg', main=True, invert=False, colour='viridis',
                  alpha=0.3, markersize=1, marker='.', linewidth=0,
                  rasterise=True, **kwargs):
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
        alpha : float, optional
            colour alpha from 0-1. Default: 0.3.
        markersize : float, optional
            marker size in points. Default: 1.
        marker : str, optional
            marker type. Default: '.'.
        linewidth : float, optional
            marker edge linewidth. Default: 0.
        rasterise : bool, optional
            rasterise plot. Default: True.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.scatter.

    Returns:
        axes
            axes with waterfall.
    """

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

    ax.scatter(x, y, s=markersize, marker=marker, facecolor=colours,
               alpha=alpha, linewidth=linewidth, rasterized=rasterise,
               **kwargs)

    # axes formatting

    if main:
        format_waterfall(ax, x, xquantity, 'x')
        format_waterfall(ax, y, quantity, 'y')

    return ax

def add_projected_waterfall(ax, data, quantity, projected,
                            xquantity='frequency', temperature=300,
                            direction='avg', main=True, invert=False,
                            colour='viridis', alpha=0.3, cmin=None, cmax=None,
                            cscale=None, unoccupied='grey', markersize=1,
                            marker='.', linewidth=0, rasterise=True, **kwargs):
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
        alpha : float
            colour alpha from 0-1. Default: 0.3.
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
        markersize : float, optional
            marker size in points. Default: 1.
        marker : str, optional
            marker type. Default: '.'.
        linewidth : float, optional
            marker edge linewidth. Default: 0.
        rasterise : bool, optional
            rasterise plot. Default: True.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.scatter.

    Returns:
        axes
            axes with projected waterfall.
        colorbar
            colourbar for projected data.
    """

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
        cmap = mpl.cm.get_cmap(colour)
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            cmap = colour
        elif isinstance(colour, str) and colour == 'skelton':
            cmap = tp.plot.colour.skelton()
        else:
            raise Exception('Unrecognised colour argument. '
                            'Expected a colourmap or colourmap name.')

    cnorm, extend = tp.plot.utilities.colour_scale(c, projected, cmap, cmin,
                                                   cmax, cscale, unoccupied)

    # plotting

    scat = ax.scatter(x, y, s=markersize, marker=marker, c=c, norm=cnorm,
                      cmap=cmap, alpha=alpha, linewidth=linewidth,
                      rasterized=rasterise, **kwargs)

    # axes formatting

    cbar = plt.colorbar(scat, extend=extend)
    cbar.set_alpha(1)
    cbar.set_label(axlabels[projected])
    cbar.draw_all()

    if main:
        format_waterfall(ax, x, xquantity, 'x')
        format_waterfall(ax, y, quantity, 'y')

    return ax, cbar

def format_waterfall(ax, data, quantity, axis):
    """Formats axes for waterfall plots

    Arguments:
        ax : axes
            axes to format.
        data : array-like
            data for the axis in question.
        quantity : str
            quantity name.
        axis : str
            axis (x or y).
    """

    axes =  {'x': ax.xaxis,      'y': ax.yaxis}
    label = {'x': ax.set_xlabel, 'y': ax.set_ylabel}
    limit = {'x': ax.set_xlim,   'y': ax.set_ylim}
    scale = {'x': ax.set_xscale, 'y': ax.set_yscale}

    axlabels = tp.settings.labels()
    label[axis](axlabels[quantity])

    if quantity in ['frequency', 'heat_capacity']:
        limit[axis](np.amin(data), np.amax(data))
        axes[axis].set_major_locator(tp.settings.locator()['major'])
        axes[axis].set_minor_locator(tp.settings.locator()['minor'])
    else:
        sort = np.ma.masked_invalid(np.ma.masked_equal(data, 0)).compressed()
        sort = sort[sort.argsort()]
        limit[axis](sort[int(round(len(sort)/100, 0))],
                    sort[int(round(len(sort)*99.9/100, 0))])
        scale[axis]('log')
        axes[axis].set_major_locator(tp.settings.locator()['log'])

    return

def scale_to_axes(ax, data, dmax, main, invert):
    """Scale data to fit axes.

    Works for linear and log axes.

    Arguments:
        ax : axes
            axes to fit to.
        data : array-like or dict
            data to scale.
        dmax
            maximum data value.
        main : bool
            scale to percent.
        invert : bool
            inverted axes.

    Returns:
        data
            scaled data.
    """

    if not isinstance(data, dict):
        data = {'qwerfv': data}
    if invert:
        ymin = 0 if main else ax.get_xlim()[0]
        ymax = 100 if main else ax.get_xlim()[1]
        log = ax.get_xaxis().get_scale()
    else:
        ymin = 0 if main else ax.get_ylim()[0]
        ymax = 100 if main else ax.get_ylim()[1]
        log = ax.get_yaxis().get_scale()
    if log == 'log':
        ymin = np.log10(ymin)
        ymax = np.log10(ymax)

    exempt = ['frequency', 'meta']
    for key in data:
        if key not in exempt:
            data[key] = np.multiply(data[key], (ymax-ymin) / dmax)
            data[key] = np.add(ymin, data[key])
            if log == 'log':
                data[key] = np.power(10, data[key])

    if 'qwerfv' in data:
        data = data['qwerfv']

    return data
