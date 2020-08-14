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
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tp.settings
from tp.data import resolve

def add_dos(ax, data, colour, total=False, main=False, invert=False,
            scale=True, fill=True, fillalpha=20, line=False, rasterise=False):
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

    Returns:
        axes
            axes with DoS.
    """

    exempt = ['x', 'meta', 'total']
    if total:
        dmax = max(data['total'])
    else:
        dmax = 0
        for key in data:
            if key not in exempt:
                kmax = float(np.amax(data[key]))
                if kmax > dmax: dmax = kmax

    # Scaling

    if scale:
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

        for key in data:
            if key not in exempt:
                data[key] = np.multiply(data[key], (ymax-ymin) / dmax)
                data[key] = np.add(ymin, data[key])
                if log == 'log':
                    data[key] = np.power(10, data[key])
        dmax = ymax

    # Plot

    fillcolour = {}
    for c in colour:
        if fill:
            fillcolour[c] = colour[c][:7] + str(fillalpha)
        else:
            colour[c][:7] + '00'
        if not line: colour[c] = fillcolour[c]

    if total:
        if invert:
            ax.fill_between(data['total'], data['x'], label='Total',
                            facecolor=fillcolour['total'],
                            edgecolor=colour['total'], rasterized=rasterise)
        else:
            ax.fill_between(data['x'], data['total'], label='Total',
                            facecolor=fillcolour['total'],
                            edgecolor=colour['total'], rasterized=rasterise)

    for key in data:
        if key not in exempt:
            if invert:
                ax.fill_between(data[key], data['x'], label=key,
                                facecolor=fillcolour[key],
                                edgecolor=colour[key], rasterized=rasterise)
            else:
                ax.fill_between(data['x'], data[key], label=key,
                                facecolor=fillcolour[key],
                                edgecolor=colour[key], rasterized=rasterise)

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
                  fill=False, fillalpha=20, line=True, rasterise=False):
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

    Returns:
        axes
            axes with cumulated kappa.
    """

    from tp.calculate import cumulate

    data = resolve.resolve(data, 'mode_kappa', temperature=temperature,
                                               direction=direction)
    k = np.ravel(data['mode_kappa'])
    f = np.ravel(data['frequency'])

    f, k = cumulate(f, k)
    np.savetxt('cumkappa-frequency-{:.0f}K-{}.dat'.format(temperature, direction),
            np.transpose([f, k]), header='Frequency(THz) k_l(Wm-1K-1)')
    f = np.append(f, 100*f[-1])
    k = np.append(k, k[-1])

    if scale:
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

        k = np.multiply(k, (ymax-ymin) / k[-1])
        k = np.add(k, ymin)
        if log == 'log':
            k = np.power(10, k)

    fillcolour = colour + str(fillalpha) if fill else colour + '00'
    if not line: colour = fillcolour
    if invert:
        ax.fill_between(k, f, label='$\mathregular{{{}}}$'.format(legend),
                        facecolor=fillcolour, edgecolor=colour,
                        rasterized=rasterise)
    else:
        ax.fill_between(f, k, label='$\mathregular{{{}}}$'.format(legend),
                        facecolor=fillcolour, edgecolor=colour,
                        rasterized=rasterise)

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

def add_waterfall(ax, data, quantity, temperature=300, direction='avg',
                  main=True, invert=False, colour='viridis', alpha=0.3,
                  rasterise=True):
    """Adds a waterfall plot of quantities against frequency.

    Arguments:
        ax : axes
            axes to plot on.
        data : dict
            data including frequency and quantity.
        quantity : str
            y-axis quantity. Accepts gamma, group_velocity, gv_by_gv,
            heat_capacity, lifetime, mean_free_path or mode_kappa.

        temperature : float, optional
            temperature in K. Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            plot frequency on y-axis. Default: False.

        colour : colourmap or str or array-like, optional
            colourmap or colourmap name or list of colours (one for
            each band or one for each point) or a single colour.
            Default: viridis.
        alpha : float
            colour alpha from 0-1. Default: 0.3.
        rasterise : bool, optional
            rasterise plot. Default: True.

    Returns:
        axes
            axes with waterfall.
    """

    tnames = tp.settings.to_tp()
    quantity = tnames[quantity] if quantity in tnames else quantity
    data = tp.data.resolve.resolve(data, quantity, temperature, direction)
    f = np.ravel(data['frequency'])
    q = np.abs(np.ravel(data[quantity]))

    fs = np.shape(data['frequency'])
    try:
        colour = mpl.cm.get_cmap(colour)
        colours = [colour(i) for i in np.linspace(0, 1, fs[1])]
        colours = np.tile(colours, (fs[0], 1))
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            colours = [colour(i) for i in np.linspace(0, 1, fs[1])]
            colours = np.tile(colours, (fs[0], 1))
        elif isinstance(colour, str) and colour == 'skelton':
            from tp import plot
            colour = plot.colour.skelton()
            colours = [colour(i) for i in np.linspace(0, 1, fs[1])]
            colours = np.tile(colours, (fs[0], 1))
        elif (isinstance(colour, list) or isinstance(colour, np.ndarray) and
              len(colour) == fs[1]):
            if np.ndim(colour) == 1:
                colours = np.tile(colour, fs[0])
            elif np.ndim(colour) == 2:
                colours = np.tile(colour, (fs[0], 1))
        else:
            colours = colour

    if invert:
        ax.scatter(q, f, s=1, marker='.', facecolor=colours, alpha=alpha,
                   linewidth=0, rasterized=rasterise)
    else:
        ax.scatter(f, q, s=1, marker='.', facecolor=colours, alpha=alpha,
                   linewidth=0, rasterized=rasterise)

    if main:
        qsort = np.ma.masked_invalid(q)
        qsort = qsort[qsort.argsort()]
        axlabels = tp.settings.labels()
        if invert:
            ax.set_xlabel(axlabels[quantity])
            ax.set_ylabel(axlabels['frequency'])

            ax.set_xlim(qsort[int(round(len(qsort)/100, 0))],
                        qsort[int(round(len(qsort)*99.9/100, 0))])
            ax.set_ylim(0, np.amax(f))

            ax.set_xscale('log')
            ax.xaxis.set_major_locator(tp.settings.locator()['log'])
            ax.yaxis.set_major_locator(tp.settings.locator()['major'])
            ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])

        else:
            ax.set_ylabel(axlabels[quantity])
            ax.set_xlabel(axlabels['frequency'])

            ax.set_ylim(qsort[int(round(len(qsort)/100, 0))],
                        qsort[int(round(len(qsort)*99.9/100, 0))])
            ax.set_xlim(0, np.amax(f))

            ax.set_yscale('log')
            ax.yaxis.set_major_locator(tp.settings.locator()['log'])
            ax.xaxis.set_major_locator(tp.settings.locator()['major'])
            ax.xaxis.set_minor_locator(tp.settings.locator()['minor'])

    return ax

def add_projected_waterfall(ax, data, quantity, projected, temperature=300,
                            direction='avg', main=True, invert=False,
                            colour='viridis', alpha=0.3, rasterise=True):
    """Adds a waterfall plot against frequency with a colour axis.

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
        rasterise : bool, optional
            rasterise plot. Default: True.

    Returns:
        axes
            axes with projected waterfall.
        colorbar
            colourbar for projected data.
    """

    tnames = tp.settings.to_tp()
    axlabels = tp.settings.labels()
    quantity = tnames[quantity] if quantity in tnames else quantity
    projected = tnames[projected] if projected in tnames else projected
    data = tp.data.resolve.resolve(data, [quantity, projected], temperature,
                                                                     direction)
    f = np.ravel(data['frequency'])
    q = np.abs(np.ravel(data[quantity]))
    p = np.abs(np.ravel(data[projected]))

    fs = np.shape(data['frequency'])
    try:
        cmap = mpl.cm.get_cmap(colour)
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            cmap = colour
        elif isinstance(colour, str) and colour == 'skelton':
            from tp import plot
            cmap = plot.colour.skelton()
        else:
            raise Exception('Unrecognised colour argument. '
                            'Expected a colourmap or colourmap name.')

    csort = np.ravel(np.ma.masked_invalid(p).compressed())
    csort = csort[csort.argsort()]
    cmin = csort[int(round(len(csort)/100, 0))]
    extend = False
    if projected == 'occupation' and cmin < 1:
        cmin = 1
        extend = True
    cmax = csort[-1]
    cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    if invert:
        scat = ax.scatter(q, f, s=1, marker='.', c=p, norm=cnorm, cmap=cmap,
                          alpha=alpha, linewidth=0, rasterized=rasterise)
    else:
        scat = ax.scatter(f, q, s=1, marker='.', c=p, norm=cnorm, cmap=cmap,
                   alpha=alpha, linewidth=0, rasterized=rasterise)

    cbar = plt.colorbar(scat, extend='min') if extend else plt.colorbar(scat)
    cbar.set_alpha(1)
    cbar.set_label(axlabels[projected])
    cbar.draw_all()

    if main:
        qsort = np.ma.masked_invalid(q)
        qsort = qsort[qsort.argsort()]
        if invert:
            ax.set_xlabel(axlabels[quantity])
            ax.set_ylabel(axlabels['frequency'])

            ax.set_xlim(qsort[int(round(len(qsort)/100, 0))], qsort[-1])
            ax.set_ylim(0, np.amax(f))

            ax.set_xscale('log')
            ax.yaxis.set_major_locator(tp.settings.locator()['major'])
            ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])

        else:
            ax.set_ylabel(axlabels[quantity])
            ax.set_xlabel(axlabels['frequency'])

            ax.set_ylim(qsort[int(round(len(qsort)/100, 0))], qsort[-1])
            ax.set_xlim(0, np.amax(f))

            ax.set_yscale('log')
            ax.xaxis.set_major_locator(tp.settings.locator()['major'])
            ax.xaxis.set_minor_locator(tp.settings.locator()['minor'])

    return ax, cbar
