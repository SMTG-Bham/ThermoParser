"""Functions for plotting against mean free path.

Functions
---------

    add_cum_kappa:
        cumulative kappa vs mean free path.

    add_markers:
        add markers lines to a linear plot.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import tp.settings
import warnings
from tp.data import resolve

warnings.filterwarnings('ignore', module='matplotlib')

def add_cum_kappa(ax, data, kmin=1, temperature=300, direction='avg',
                  xmarkers=None, ymarkers=None, add_xticks=False,
                  add_yticks=False, main=True, scale=False, colour='#000000',
                  fill=False, fillcolour=0.2, line=True, markerkwargs={},
                  **kwargs):
    """Cumulates and plots kappa against mean free path.

    Arguments
    ---------

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
        add_xticks : bool, optional
            add x_ticks for each marker. Doesn't work on log axes.
            Default: False.
        add_yticks : bool, optional
            add y_ticks for each marker. Doesn't work on log axes.
            Default: False.

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
            if a float from 0-1 and colour in #RRGGBB format, sets
            fill colour opacity. Otherwise treats it as a colour.
            Default: 0.2.
        line : bool, optional
            plot lines. Default: True.

        markerkwargs : dict, optional
            keyword arguments for the markers, passed to
            matplotlib.pyplot.plot.
            Defaults:

                color:      black
                rasterized: False

        **kwargs
            keyword arguments passed to matplotlib.pyplot.fill_between
            if filled or matplotlib.pyplot.plot otherwise.
            Defaults:

                label:      $\mathregular{\kappa_l}$
                rasterized: False
    """

    from tp.calculate import cumulate

    # defaults

    defkwargs = {'label':      '$\mathregular{\kappa_l}$',
                 'rasterized': False}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # input checks

    for name, value in zip(['main', 'scale', 'fill', 'line', 'add_xticks',
                                                             'add_yticks'],
                           [ main,   scale,   fill,   line,   add_xticks,
                                                              add_yticks]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)
    assert fill or line, 'fill or line or both must be True.'

    # data formatting and calculation

    data = resolve.resolve(data, ['mode_kappa', 'mean_free_path'],
                           temperature=temperature, direction=direction)
    k = np.ravel(data['mode_kappa'])
    mfp = np.abs(np.ravel(data['mean_free_path']))

    mfp, k = cumulate(mfp, k)
    np.savetxt('cumkappa-mfp-{:.0f}K-{}.dat'.format(
               data['meta']['temperature'], direction),
               np.transpose([mfp, k]), header='mfp(m) k_l(Wm-1K-1)')
    mfp = np.append(mfp, 100*mfp[-1])
    k = np.append(k, k[-1])

    # percentage output

    if scale:
        axscale = [0, 100] if main else None
        k, _ = tp.plot.utilities.scale_to_axis(ax, k, scale=axscale)

    # colour

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
                fillcolour2 = colour
        if not line: colour = fillcolour2

        # plotting

        ax.fill_between(mfp, k, facecolor=fillcolour2, edgecolor=colour,
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
        tp.plot.utilities.set_locators(ax, x='log', y='linear')

    if xmarkers is not None or ymarkers is not None:
        add_markers(ax, mfp, k, xmarkers, ymarkers, add_xticks, add_yticks,
                    **markerkwargs)

    return

def add_markers(ax, x, y, xmarkers=None, ymarkers=None, add_xticks=False,
                add_yticks=False, **kwargs):
    """Adds marker lines linking a linear plot to the axes.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        x : array-like
            x data.
        y : array-like
            y data.

        xmarkers : array-like, int or float, optional
            where on the x axis to mark.
        ymarkers : array-like, int or float, optional
            where on the y axis to mark.
        add_xticks : bool, optional
            add x_ticks for each marker. Doesn't work on log axes.
            Default: False.
        add_yticks : bool, optional
            add y_ticks for each marker. Doesn't work on log axes.
            Default: False.

        **kwargs
            keyword arguments passed to matplotlib.pyplot.plot.
            Defaults:

                color:      black
                linewidth:  axes line width
                rasterized: False
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
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()

    # plot x

    if xmarkers is not None:
        if isinstance(xmarkers, (int, float)): xmarkers = [xmarkers]
        yinter = interp1d(x, y)
        xmarky = yinter(xmarkers)
        for m in range(len(xmarkers)):
            ax.plot([xmarkers[m], xmarkers[m], 0],
                    [0,           xmarky[m],   xmarky[m]], **kwargs)
        xticks = np.append(xticks, xmarkers)
        yticks = np.append(yticks, xmarky)

    # plot y

    if ymarkers is not None:
        if isinstance(ymarkers, (int, float)): ymarkers = [ymarkers]
        xinter = interp1d(y, x)
        ymarkx = xinter(ymarkers)
        for m in range(len(ymarkers)):
            ax.plot([ymarkx[m], ymarkx[m],   0],
                    [0,         ymarkers[m], ymarkers[m]], **kwargs)
        yticks = np.append(yticks, ymarkers)
        xticks = np.append(xticks, ymarkx)

    if add_xticks:
        ax.set_xticks(xticks)
    if add_yticks:
        ax.set_yticks(yticks)

    return
