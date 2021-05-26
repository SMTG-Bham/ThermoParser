"""Function for plotting mean free path on the x-axis.

Functions
---------

    add_cum_kappa:
        cumulative kappa vs mean free path.

    add_markers:
        add markers lines to a linear plot.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import tp.settings
import warnings
import yaml
from tp.data import resolve

warnings.filterwarnings('ignore', module='matplotlib')

try:
    filename = '{}/.config/tprc.yaml'.format(os.path.expanduser("~"))
    with open(filename, 'r') as f:
        conf = yaml.safe_load(f)
except Exception:
    conf = None

def add_cum_kappa(ax, data, kmin=1, temperature=300, direction='avg',
                  label=None, xmarkers=None, ymarkers=None, add_xticks=False,
                  add_yticks=False, main=True, scale=False, colour='#000000',
                  fill=False, fillcolour=0.2, line=True, linestyle='-',
                  marker=None, markerkwargs={}, **kwargs):
    """Cumulates and plots kappa against mean free path.

    Can plot data from multiple data dictionaries and directions.
    Colour, linestyle etc. are looped, so if you have two datasets and
    two directions, but only specify two colours, the first will apply
    to the first direction in both datasets, whereas if you want one for
    the first dataset and one for the second, you would repeat the first
    colour twice and the second twice too.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            (list of sets of) Phono3py data including:

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
            (list of) direction(s) from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.
        label : str, optional
            (list of) legend label(s). Defaults to $\mathregular{\kappa_l}$
            if there is only one line plotted, or direction if there are
            more.

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

        colour : str or list, optional
            (list of) RGB line colour(s). Default: default colour cycle.
        fill : bool, optional
            fill below lines. Default: False.
        fillcolour : int or str or list, optional
            if a float from 0-1 and colour in #RRGGBB format, sets
            fill colour opacity. Otherwise treats it as a colour.
            Default: 0.2.
        line : bool, optional
            plot lines. Default: True.
        linestyle : str or list, optional
            (list of) linestyle(s). Default: "-".
        marker : str or list, optional
            (list of) markers. Default: None.

        markerkwargs : dict, optional
            keyword arguments for the markers, passed to
            matplotlib.pyplot.plot.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                color:      black
                rasterized: False

        kwargs
            keyword arguments passed to matplotlib.pyplot.fill_between
            if filled or matplotlib.pyplot.plot otherwise.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Default:

                rasterized: False

    Returns
    -------

        None
            adds plot directly to ax.
    """

    from tp.calculate import cumulate

    # defaults

    defkwargs = {'rasterized': False}

    if conf is None or 'mfp_cum_kappa_kwargs' not in conf or \
       conf['mfp_cum_kappa_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['mfp_cum_kappa_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['main', 'scale', 'fill', 'line', 'add_xticks',
                                                             'add_yticks'],
                           [ main,   scale,   fill,   line,   add_xticks,
                                                              add_yticks]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)
    assert fill or line, 'fill or line or both must be True.'

    # data formatting and calculation

    if isinstance(direction, str):
        direction = [direction]
    
    if isinstance(data, dict):
        data = [data]
    
    if colour is None:
        colour = plt.rcParams['axes.prop_cycle'].by_key()['color']
    elif isinstance(colour, str):
        colour = [colour]
    
    if fillcolour is None:
        fillcolour = plt.rcParams['axes.prop_cycle'].by_key()['color']
    elif isinstance(fillcolour, (str, float, int)):
        fillcolour = [fillcolour]
    
    if isinstance(linestyle, str):
        linestyle = [linestyle]
    
    if marker is None or isinstance(marker, str):
        marker = [marker]
    
    if label is None:
        if len(data) == 1 and len(direction) == 1:
            label = ['$\mathregular{\kappa_l}$']
        else:
            label = direction
    elif isinstance(label, str):
        label = [label]

    mfpmin, mfpmax, kmax = None, None, None
    xticks, yticks = [], []
    i = 0

    for dat in data:
        for d in direction:
            data2 = resolve.resolve(dat, ['mode_kappa', 'mean_free_path'],
                                    temperature=temperature, direction=d)
            k = np.ravel(data2['mode_kappa'])
            mfp = np.abs(np.ravel(data2['mean_free_path']))

            mfp, k = cumulate(mfp, k)
            np.savetxt('cumkappa-mfp-{:.0f}K-{}.dat'.format(
                       dat['meta']['temperature'], d), np.transpose([mfp, k]),
                       header='mfp(m) k_l(Wm-1K-1)')

            mindex = next(x[0] for x in enumerate(k) if x[1] > kmin*k[-1]/100)
            if mfpmin is None or mfpmin > mfp[mindex]:
                mfpmin = mfp[mindex]
            if mfpmax is None or mfpmax < mfp[-1]:
                mfpmax = mfp[-1]
            if kmax is None or kmax < k[-1]:
                kmax = 100 if main and scale else k[-1]

            mfp = np.append(mfp, 100*mfp[-1])
            k = np.append(k, k[-1])

            # percentage output

            if scale:
                axscale = [0, 100] if main else None
                k = tp.plot.utilities.scale_to_axis(ax, k, scale=axscale)

            colour1 = colour[i % len(colour)]
            fillcolour1 = fillcolour[i % len(fillcolour)]
            linestyle1 = linestyle[i % len(linestyle)]
            marker1 = marker[i % len(marker)]
            label1 = "${}$".format(label[i % len(label)])

            # colour
            # Tries to read the colour as an rgb code, then alpha value.

            if fill:
                try:
                    fillcolour2 = tp.plot.colour.rgb2array(colour1, fillcolour1)
                except Exception:
                    if isinstance(colour1, list) and \
                       isinstance(fillcolour1, (float, int)) and \
                       fillcolour1 >= 0 and fillcolour1 <= 1:
                        fillcolour2 = colour1
                        if len(colour1) == 3:
                            fillcolour2.append(fillcolour1)
                        elif len(colour1) == 4:
                            fillcolour2[3] = fillcolour1
                    else:
                        fillcolour2 = fillcolour1

            # plotting

            if fill and line:
                ax.plot(mfp, k, color=colour1, linestyle=linestyle1,
                        marker=marker1, label=label1, **kwargs)
                ax.fill_between(mfp, k, facecolor=fillcolour2, edgecolor=colour)

            elif fill and not line:
                ax.fill_between(mfp, k, facecolor=fillcolour2, **kwargs)

            else:
                ax.plot(mfp, k, color=colour1, linestyle=linestyle1,
                        marker=marker1, label=label1, **kwargs)

            if xmarkers is not None or ymarkers is not None:
                xtick, ytick = add_markers(ax, mfp, k, xmarkers, ymarkers,
                                           **markerkwargs)
                xticks.append(xtick)
                yticks.append(ytick)

            i += 1

    # axes formatting

    if main:
        axlabels = tp.settings.labels()
        if scale:
            ax.set_ylabel(axlabels['cumulative_percent'])
        else:
            ax.set_ylabel(axlabels['cumulative_kappa'])
        ax.set_xlabel(axlabels['mean_free_path'])
        ax.set_ylim(0, kmax)
        ax.set_xlim(mfpmin, mfpmax)
        tp.plot.utilities.set_locators(ax, x='log', y='linear')

    if add_xticks:
        ax.set_xticks(np.append(ax.get_xticks(), xticks))
    if add_yticks:
        ax.set_yticks(np.append(ax.get_yticks(), yticks))

    return

def add_markers(ax, x, y, xmarkers=None, ymarkers=None, **kwargs):
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

        kwargs
            keyword arguments passed to matplotlib.pyplot.plot.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                color:      black
                linewidth:  axes line width
                rasterized: False

    Returns
    -------

        list
            added x locations.
        list
            added y locations.
    """

    from scipy.interpolate import interp1d

    # defaults

    defkwargs = {'color':      'black',
                 'linewidth':  ax.spines['bottom'].get_linewidth(),
                 'rasterized': False}

    if conf is None or 'marker_kwargs' not in conf or \
       conf['marker_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['marker_kwargs'], **kwargs}

    xticks, yticks = [], []

    # plot x

    if xmarkers is not None:
        if isinstance(xmarkers, (int, float)): xmarkers = [xmarkers]
        yinter = interp1d(x, y)
        xmarky = yinter(xmarkers)
        for m in range(len(xmarkers)):
            ax.plot([xmarkers[m], xmarkers[m], 0],
                    [0,           xmarky[m],   xmarky[m]], **kwargs)
        xticks.append(xmarkers)
        yticks.append(xmarky)

    # plot y

    if ymarkers is not None:
        if isinstance(ymarkers, (int, float)): ymarkers = [ymarkers]
        xinter = interp1d(y, x)
        ymarkx = xinter(ymarkers)
        for m in range(len(ymarkers)):
            ax.plot([ymarkx[m], ymarkx[m],   0],
                    [0,         ymarkers[m], ymarkers[m]], **kwargs)
        yticks.append(ymarkers)
        xticks.append(ymarkx)

    return xticks, yticks
