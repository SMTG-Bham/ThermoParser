"""Heatmap plotters.

Functions:
    add_heatmap:
        heatmap.
    add_ztmap:
        ZT vs temperature and carrier concentration.
    add_kappa_target:
        kappa_l needed to reach a given ZT.
"""

import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tp
import warnings
from matplotlib import cm
from scipy.interpolate import interp1d, interp2d

warnings.filterwarnings('ignore', module='matplotlib')

def add_heatmap(ax, x, y, c, xinterp=None, yinterp=None, kind='linear',
                xscale='linear', yscale='linear', cscale='linear',
                cmin=None, cmax=None, colour='viridis', **kwargs):
    """Adds a heatmap to a set of axes.

    Formats limits, parses extra colourmap options, makes sure data
    isn't obscured.

    Arguments:
        ax : axes
            axes to plot on.
        x : array-like
            x data.
        y : array-like
            y data.
        c : array-like (shape: x, y)
            colour data.

        xinterp : int, optional
            density of interpolation. Default: None.
        yinterp : int, optional
            density of interpolation. Default: None.
        kind : str, optional
            interpolation kind. Default: linear.

        xscale : str, optional
            x scale (linear or log). Default: linear.
        yscale : str, optional
            y scale (linear or log). Default: linear
        cscale : str, optional
            colour scale (linear or log). Default: linear.
        cmin : float, optional
            override colour scale minimum. Default: None.
        cmax : float, optional
            override colour scale maximum. Default: None.

        colour : colourmap or str or array-like, optional
            colourmap or colourmap name; or key colour or min and max
            RGB colours to generate a colour map. Default: viridis.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.pcolormesh.
            Defaults: {'rasterized': False}

    Returns:
        colourbar
            colourbar.
    """

    # defaults

    defkwargs = {'rasterized': False}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # colour

    extend = 'neither'
    if cmin is None:
        cmin = np.amin(c)
    elif cmin > np.amin(c):
        extend = 'min'
    if cmax is None:
        cmax = np.amax(c)
    elif cmax < np.amax(c):
        if extend == 'min':
            extend = 'both'
        else:
            extend = 'max'

    if cscale == 'linear':
        cnorm = colors.Normalize(vmin=cmin, vmax=cmax)
    elif cscale == 'log':
        cnorm = colors.LogNorm(vmin=cmin, vmax=cmax)

    try:
        colours = mpl.cm.get_cmap(colour)
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            colours = colour
        else:
            try:
                colours = tp.plot.colour.elbow(colour)
            except Exception:
                colours = tp.plot.colour.linear(colour[1], colour[0])

    # data interpolation

    if xinterp is not None or yinterp is not None:
        cinterp = interp2d(x, y, c, kind=kind)
        if xinterp is not None:
            if xscale == 'linear': x = np.linspace(x[0], x[-1], xinterp)
            if xscale == 'log': x = np.geomspace(x[0], x[-1], xinterp)
        if yinterp is not None:
            if yscale == 'linear': y = np.linspace(y[0], y[-1], yinterp)
            if yscale == 'log': y = np.geomspace(y[0], y[-1], yinterp)
        c = cinterp(x, y)

    # ensure all data is shown
    if len(x) == len(c[0]):
        x = list(x)
        if xscale == 'linear':
            x.append(2 * x[-1] - x[-2])
        elif xscale =='log':
            x.append(x[-1] ** 2 / x[-2])
    if len(y) == len(c):
        y = list(y)
        if yscale == 'linear':
            y.append(2 * y[-1] - y[-2])
        elif yscale =='log':
            y.append(y[-1] ** 2 / y[-2])

    # plotting

    heat = ax.pcolormesh(x, y, c, cmap=colours, norm=cnorm, **kwargs)
    cbar = plt.colorbar(heat, extend=extend)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    # axes formatting
    tp.settings.set_locators(ax, x=xscale, y=yscale)
    tp.settings.set_locators(cbar.ax, y=cscale)

    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[0], y[-1])

    return cbar

def add_ztmap(ax, data, kdata=None, direction='avg', xinterp=None,
              yinterp=None, kind='linear', cmin=None, cmax=None,
              colour='viridis', **kwargs):
    """Convenience wrapper for plotting ZT heatmaps.

    Calculates ZT, plots and formats labels etc.

    Arguments:
        ax : axes
            axes to plot on.
        data : dict
            dictionary containing temperature and doping, and either zt
            or conductivity, seebeck and electronic_thermal_conductivity.

        kdata : dict, optional
            dictionary containing lattice_thermal_conductivity and
            temperature. Ignored if zt in data, if not ignored and None,
            defaults to 1 at all temperatures.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.

        xinterp : int, optional
            density of interpolation. Default: None.
        yinterp : int, optional
            density of interpolation. Default: None.
        kind : str, optional
            interpolation kind. Default: linear.

        cmin : float, optional
            override colour scale minimum. Default: None.
        cmax : float, optional
            override colour scale maximum. Default: None.
        colour : colourmap or str or array-like, optional
            colourmap or colourmap name; or key colour or min and max
            RGB colours to generate a colour map. Default: viridis.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.pcolormesh.
            Defaults: {'rasterized': False}

    Returns:
        colourbar
            colourbar.
    """

    # data formatting

    if 'zt' in data:
        if np.ndim(data['zt']) == 4:
            data = tp.data.resolve(data, 'zt', direction)
    else:
        data = tp.data.resolve.resolve(data, ['conductivity', 'seebeck',
                                       'electronic_thermal_conductivity'],
                                       direction=direction)

        if kdata is not None:
            kdata = tp.data.resolve.resolve(kdata,
                                            'lattice_thermal_conductivity',
                                            direction=direction)
            kinterp = interp1d(kdata['temperature'],
                               kdata['lattice_thermal_conductivity'],
                               kind='cubic')
            data['lattice_thermal_conductivity'] = kinterp(data['temperature'])
        else: # if no kappa_lat, set to 1
            data['lattice_thermal_conductivity'] = np.ones(
                                                      len(data['temperature']))

        data = tp.calculate.zt_fromdict(data)

    # plotting

    cbar = add_heatmap(ax, data['temperature'], list(np.abs(data['doping'])),
                       np.transpose(data['zt']), xinterp=xinterp,
                       yinterp=yinterp, kind=kind, yscale='log', cmin=cmin,
                       cmax=cmax, colour=colour, **kwargs)

    # axes formatting

    labels = tp.settings.labels()
    ax.set_xlabel(labels['temperature'])
    ax.set_ylabel(labels['doping'])
    cbar.set_label(labels['zt'])

    return cbar

def add_kappa_target(ax, data, zt=2, direction='avg', xinterp=None,
                     yinterp=None, kind='linear', cmin=0, cmax=None,
                     colour='viridis', **kwargs):
    """Plots a heatmap of k_latt required for a target ZT

    Calculates lattice thermal conductivity, plots and formats labels
    etc. May be useful to screen materials to calculate kappa_l for.

    Arguments:
        ax : axes
            axes to plot on.
        data : dict
            dictionary containing temperature and doping, and either zt
            or conductivity, seebeck and electronic_thermal_conductivity.

        zt : float, optional
            target ZT. Default: 2.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.

        xinterp : int, optional
            density of interpolation. Default: None.
        yinterp : int, optional
            density of interpolation. Default: None.
        kind : str, optional
            interpolation kind. Default: linear.

        cmin : float, optional
            override colour scale minimum. Default: 0.
        cmax : float, optional
            override colour scale maximum. Default: None.
        colour : colourmap or str or array-like, optional
            colourmap or colourmap name; or key colour or min and max
            RGB colours to generate a colour map. Default: viridis.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.pcolormesh.
            Defaults: {'rasterized': False}

    Returns:
        colourbar
            colourbar.
    """

    # data formatting

    data = tp.data.resolve.resolve(data, ['conductivity', 'seebeck',
                                   'electronic_thermal_conductivity'],
                                   direction=direction)
    data['zt'] = zt

    data = tp.calculate.kl_fromdict(data)

    # plotting

    cbar = add_heatmap(ax, data['temperature'], list(np.abs(data['doping'])),
                       np.transpose(data['lattice_thermal_conductivity']),
                       xinterp=xinterp, yinterp=yinterp, kind=kind,
                       yscale='log', cmin=0, cmax=None, colour=colour,**kwargs)

    # axes formatting

    labels = tp.settings.labels()
    ax.set_xlabel(labels['temperature'])
    ax.set_ylabel(labels['doping'])
    cbar.set_label(labels['lattice_thermal_conductivity'])

    return cbar
