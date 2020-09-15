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
from tp.data import resolve

def add_cum_kappa(ax, data, kmin=1, temperature=300, direction='avg',
                  xmarkers=None, ymarkers=None, legend='\kappa_l', main=True,
                  scale=False, colour='#000000', fill=False, fillalpha=20,
                  line=True, markercolour='black', markerkwargs={},
                  rasterise=False, **kwargs):
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
        legend : str, optional
            legend entry, accepts maths. Default: \kappa_l.

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
            RGB line colour. Default: #000000.
        fill : bool, optional
            fill below lines. Default: False.
        fillalpha : int, optional
            fill alpha in %. Default: 20.
        line : bool, optional
            plot lines. Default: True.
        markercolour : str, optional
            marker line colours. Default: black.
        markerkwargs : dict, optional
          keyword arguments for the markers, passed to
          matplotlib.pyplot.plot. Default: {}.
        rasterise : bool, optional
            rasterise plot. Default: False.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.fill_between.

    Returns:
        axes
            axes with cumulated kappa.
    """

    from tp.calculate import cumulate

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
        ymin = 0 if main else ax.get_ylim()[0]
        ymax = 100 if main else ax.get_ylim()[1]
        k = np.multiply(k, (ymax-ymin) / k[-1])
        k = np.add(k, ymin)

    # colour
    fillcolour = colour + str(fillalpha) if fill else colour + '00'
    if not line: colour = fillcolour

    # plotting

    ax.fill_between(mfp, k, label='$\mathregular{{{}}}$'.format(legend),
                    facecolor=fillcolour, edgecolor=colour,
                    rasterized=rasterise, **kwargs)

    # axes formatting

    if main:
        mindex = next(x[0] for x in enumerate(k) if x[1] > kmin*k[-2]/100)
        axlabels = tp.settings.labels()
        ax.set_ylabel(axlabels['cumulative_kappa'])
        ax.set_xlabel(axlabels['mean_free_path'])
        ax.set_xscale('log')
        ax.set_ylim(0, k[-2])
        ax.set_xlim(mfp[mindex], mfp[-2])

        ax.xaxis.set_major_locator(tp.settings.locator()['log'])
        ax.yaxis.set_major_locator(tp.settings.locator()['major'])
        ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])

    if xmarkers is not None or ymarkers is not None:
        add_markers(ax, mfp, k, xmarkers, ymarkers, colour=markercolour,
                    **markerkwargs)

    return ax

def add_markers(ax, x, y, xmarkers=None, ymarkers=None, colour='black',
                **kwargs):
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

    Returns:
        axes
            axes with markers.
    """

    from scipy.interpolate import interp1d

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
                    [-1e100,      xmarky[m],   xmarky[m]], c=colour)
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

    return ax
