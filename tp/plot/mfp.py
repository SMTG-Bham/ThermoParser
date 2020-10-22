"""Functions for plotting against mean free path.

Functions:
    add_cum_kappa:
        cumulative kappa vs mean free path.

    add_markers:
        add markers lines to a linear plot.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tp.settings
import warnings
from tp.data import resolve

def add_cum_kappa(ax, data, kmin=1, temperature=300, direction='avg',
                  xmarkers=None, ymarkers=None, legend='\kappa_l', main=True,
                  scale=False, colour='#000000', fill=False, fillcolour=20,
                  line=True, markerkwargs={}, **kwargs):
    """Cumulates and plots kappa against mean free path.

    Arguments:
        ax : axes
            axes to plot on.
        data : dict
            Phono3py-like data including:
                mode_kappa: array-like
                    frequency and q-point decomposed lattice thermal
                    conductivity.
                mean_free_path : array-like
                    mean free paths.
                temperature : array-like
                    temperature.
        kmin : float or int, optional
            minimum kappa to plot in percent. Default: 1.

        temperature : float, optional
            temperature in K (finds nearest). Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average

        xmarkers : array-like, optional
            mark kappa at certain mean free paths.
        ymarkers : array-like, optional
            mark mean free path at certain kappas.

        main : bool, optional
            set ticks, labels, limits. Default: True.
        scale : bool, optional
            if main, scale to percent. If not main, scale to axis
            limits. Default: False.

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

        markerkwargs : dict, optional
            keyword arguments for the markers, passed to
            matplotlib.pyplot.plot.
            Defaults: {'color':      'black',
                       'rasterized': False}.
        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.fill_between
            if filled or matplotlib.pyplot.plot otherwise.
            Defaults: {'label':      '$\mathregular{\kappa_l}$',
                       'rasterized': False}
    """

    from tp.calculate import cumulate

    # defaults

    defkwargs = {'label':      '$\mathregular{\kappa_l}$',
                 'rasterized': False}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # data formatting and calculation

    data = resolve.resolve(data, ['mode_kappa', 'mean_free_path'],
                           temperature=temperature, direction=direction)
    k = np.ravel(data['mode_kappa'])
    mfp = np.abs(np.ravel(data['mean_free_path']))

    mfp, k = cumulate(mfp, k)
    np.savetxt('cumkappa-mfp-{:.0f}K-{}.dat'.format(temperature, direction),
               np.transpose([mfp, k]), header='mfp(m) k_l(Wm-1K-1)')
    mfp = np.append(mfp, 100*mfp[-1])
    k = np.append(k, k[-1])

    # percentage output

    if scale:
        axscale = [0, 100] if main else None
        k, _ = tp.plot.frequency.scale_to_axes(ax, k, scale=axscale)

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

        ax.fill_between(mfp, k, facecolor=fillcolour, edgecolor=colour,
                        **kwargs)
    else:
        ax.plot(mfp, k, color=colour, **kwargs)

    # axes formatting

    if main:
        mindex = next(x[0] for x in enumerate(k) if x[1] > kmin*k[-2]/100)
        axlabels = tp.settings.labels()
        ax.set_ylabel(axlabels['cumulative_kappa'])
        ax.set_xlabel(axlabels['mean_free_path'])
        ax.set_ylim(0, k[-2])
        ax.set_xlim(mfp[mindex], mfp[-2])
        tp.settings.set_locators(ax, x='log', y='linear')

    if xmarkers is not None or ymarkers is not None:
        add_markers(ax, mfp, k, xmarkers, ymarkers, **markerkwargs)

    return

def add_markers(ax, x, y, xmarkers=None, ymarkers=None, **kwargs):
    """Adds marker lines linking a linear plot to the axes.

    Args:
        ax : axes
            axes to plot on.
        x : array-like
            x data.
        y : array-like
            y data.
        xmarkers : array-like, optional
            where on the x axis to mark.
        ymarkers : array-like, optional
            where on the y axis to mark.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.plot.
            Defaults: {'color':      'black',
                       'linewidth':  axes line width,
                       'rasterized': False}.
    """

    from scipy.interpolate import interp1d

    # defaults

    defkwargs = {'color':      'black',
                 'linewidth':  ax.spines['bottom'].get_linewidth(),
                 'rasterized': False}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # get limits

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    #xticks = list(ax.get_xticks())
    #yticks = list(ax.get_yticks())

    # plot x

    if xmarkers is not None:
        yinter = interp1d(x, y)
        xmarky = yinter(xmarkers)
        for m in range(len(xmarkers)):
            ax.plot([xmarkers[m], xmarkers[m], -1e100],
                    [-1e100,      xmarky[m],   xmarky[m]], **kwargs)
        #xticks = np.append(xticks, xmarkers)
        #yticks = np.append(yticks, xmarky)

    # plot y

    if ymarkers is not None:
        xinter = interp1d(y, x)
        ymarkx = yinter(ymarkers)
        for m in range(len(ymarkers)):
            ax.axvline(ymarkers[m], ymax=ymarkx[m], c='black')
            ax.axhline(ymarkx[m], xmax=ymarkers[m], c='black')
        #yticks = np.append(yticks, ymarkers)
        #xticks = np.append(xticks, ymarkx)

    #ax.set_xticks(xticks)
    #ax.set_yticks(yticks)
    #ax.set_xlim(xlim)
    #ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))

    return
