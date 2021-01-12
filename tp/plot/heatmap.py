"""Heatmap plotters.

Functions
---------

    add_heatmap:
        heatmap.
    add_ztmap:
        ZT vs temperature and carrier concentration.
    add_kappa_target:
        kappa_l needed to reach a given ZT.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import tp
import warnings
from scipy.interpolate import interp1d, interp2d
import matplotlib.ticker as ticker

warnings.filterwarnings('ignore', module='matplotlib')

def add_heatmap(ax, x, y, c, xinterp=None, yinterp=None, kind='linear',
                xscale='linear', yscale='linear', cscale='linear', xmin=None,
                xmax=None, ymin=None, ymax=None, cmin=None, cmax=None,
                colour='viridis', **kwargs):
    """Adds a heatmap to a set of axes.

    Formats limits, parses extra colourmap options, makes sure data
    isn't obscured.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        x : array-like
            x data.
        y : array-like
            y data.
        c : array-like (shape: x, y)
            colour data.

        xinterp : int, optional
            density of interpolation. None turns it off. Default: None.
        yinterp : int, optional
            density of interpolation. None turns it off. Default: None.
        kind : str, optional
            interpolation kind. Default: linear.

        xscale : str, optional
            x scale (linear or log). Default: linear.
        yscale : str, optional
            y scale (linear or log). Default: linear
        cscale : str, optional
            colour scale (linear or log). Default: linear.
        xmin : float, optional
            override x minimum. Default: None.
        xmax : float, optional
            override x maximum. Default: None.
        ymin : float, optional
            override y minimum. Default: None.
        ymax : float, optional
            override y maximum. Default: None.
        cmin : float, optional
            override colour scale minimum. Default: None.
        cmax : float, optional
            override colour scale maximum. Default: None.

        colour : colourmap or str or array-like, optional
            colourmap or colourmap name; or key colour or min and max
            RGB colours to generate a colour map. Default: viridis.

        **kwargs
            keyword arguments passed to matplotlib.pyplot.pcolormesh.
            Defaults:

                rasterized: False

    Returns
    -------

        colourbar
            colourbar.
    """

    # defaults

    defkwargs = {'rasterized': False}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # data trimming

    x = np.array(x)
    y = np.array(y)
    c = np.array(c)
    if xmin is None: xmin = x[0]
    if xmax is None: xmax = x[-1]
    if ymin is None: ymin = y[0]
    if ymax is None: ymax = y[-1]
    xi = np.where((x >= xmin) & (x <= xmax))[0]
    yi = np.where((y >= ymin) & (y <= ymax))[0]
    x = x[xi]
    y = y[yi]

    try:
        c = c[np.ix_(xi, yi)]
    except IndexError:
        c = c[np.ix_(xi[:-1], yi[:-1])]

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
        cnorm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
    elif cscale == 'log':
        cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    try:
        colours = mpl.cm.get_cmap(colour)
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            colours = colour
        else:
            try:
                colours = tp.plot.colour.uniform(colour)
            except Exception:
                colours = tp.plot.colour.linear(colour[1], colour[0])

    # data interpolation

    if xinterp is not None or yinterp is not None:
        cinterp = interp2d(x, y, np.transpose(c), kind=kind)
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

    # axes formatting
    tp.plot.utilities.set_locators(ax, x=xscale, y=yscale)
    tp.plot.utilities.set_locators(cbar.ax, y=cscale)

    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[0], y[-1])

    return cbar

def add_ztmap(ax, data, kdata=None, direction='avg', xinterp=200,
              yinterp=200, kind='linear', xmin=None, xmax=None, ymin=None,
              ymax=None, cmin=None, cmax=None, colour='viridis',
              output='zt.hdf5', **kwargs):
    """Convenience wrapper for plotting ZT heatmaps.

    Calculates ZT and writes to hdf5, plots and formats labels etc.

    Arguments
    ---------

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
            density of interpolation. None turns it off. Default: 200.
        yinterp : int, optional
            density of interpolation. None turns it off. Default: 200.
        kind : str, optional
            interpolation kind. Default: linear.

        xmin : float, optional
            override temperature minimum. Default: None.
        xmax : float, optional
            override temperature maximum. Default: None.
        ymin : float, optional
            override doping minimum. Default: None.
        ymax : float, optional
            override doping maximum. Default: None.
        cmin : float, optional
            override ZT minimum. Default: None.
        cmax : float, optional
            override ZT maximum. Default: None.
        colour : colourmap or str or array-like, optional
            colourmap or colourmap name; or key colour or min and max
            RGB colours to generate a colour map. Default: viridis.

        output : str, optional
            output filename to write to. Only writes if ZT has been
            within this function. Set to None to not write. Can be
            loaded with h5py.File to put back into this function.
            Default: zt.hdf5.

        **kwargs
            keyword arguments passed to matplotlib.pyplot.pcolormesh.
            Defaults:

                rasterized: False

    Returns
    -------

        colourbar
            colourbar.
    """

    import h5py

    # data formatting

    if 'zt' in data:
        if np.ndim(data['zt']) == 4:
            data = tp.data.resolve(data, 'zt', direction)
        else:
            data = dict(data)
    else:
        if np.ndim(data['conductivity']) == 4:
            data = tp.data.resolve.resolve(data, ['conductivity', 'seebeck',
                                           'electronic_thermal_conductivity'],
                                           direction=direction)
        else:
            data = dict(data)

        if 'meta' not in data:
            data['meta'] = {}

        if kdata is not None:
            kdata = tp.data.resolve.resolve(kdata,
                                            'lattice_thermal_conductivity',
                                            direction=direction)

            # scale to smallest data set (don't extrapolate)
            # interpolation takes care of it if data is smaller than kdata
            tindex = np.where((data['temperature'] <= kdata['temperature'][-1])
                            & (data['temperature'] >= kdata['temperature'][0]))
            data['temperature'] = np.array(data['temperature'])[tindex[0]]
            data['conductivity'] = np.array(data['conductivity'])[tindex[0]]
            data['electronic_thermal_conductivity'] = \
                   np.array(data['electronic_thermal_conductivity'])[tindex[0]]
            data['seebeck'] = np.array(data['seebeck'])[tindex[0]]

            kinterp = interp1d(kdata['temperature'],
                               kdata['lattice_thermal_conductivity'],
                               kind='cubic')
            data['lattice_thermal_conductivity'] = kinterp(data['temperature'])
            data['meta']['kappa_source'] = kdata['meta']['kappa_source']
        else: # if no kappa_lat, set to 1
            data['lattice_thermal_conductivity'] = np.ones(
                                                   len(data['temperature']))
            data['meta']['kappa_source'] = 'Set to 1 W m^-1 K^-1'

        data = tp.calculate.zt_fromdict(data)

        if output is not None:
            tp.data.save.hdf5(data, output)

    # plotting

    cbar = add_heatmap(ax, data['temperature'], list(np.abs(data['doping'])),
                       data['zt'], xinterp=xinterp,
                       yinterp=yinterp, kind=kind, yscale='log', xmin=xmin,
                       xmax=xmax, ymin=ymin, ymax=ymax, cmin=cmin, cmax=cmax,
                       colour=colour, **kwargs)

    # axes formatting

    labels = tp.settings.labels()
    ax.set_xlabel(labels['temperature'])
    ax.set_ylabel(labels['doping'])
    cbar.set_label(labels['zt'])

    ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    return cbar

def add_kappa_target(ax, data, zt=2, direction='avg', xinterp=200,
                     yinterp=200, kind='linear', xmin=None, xmax=None,
                     ymin=None, ymax=None, cmin=0, cmax=None, colour='viridis',
                     output='target-kl.hdf5', **kwargs):
    """Plots a heatmap of k_latt required for a target ZT

    Calculates lattice thermal conductivity, plots and formats labels
    etc. May be useful to screen materials to calculate kappa_l for.

    Arguments
    ---------

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
            density of interpolation. None turns it off. Default: 200.
        yinterp : int, optional
            density of interpolation. None turns it off. Default: 200.
        kind : str, optional
            interpolation kind. Default: linear.

        xmin : float, optional
            override temperature minimum. Default: None.
        xmax : float, optional
            override temperature maximum. Default: None.
        ymin : float, optional
            override doping minimum. Default: None.
        ymax : float, optional
            override doping maximum. Default: None.
        cmin : float, optional
            override kappa minimum. Default: 0.
        cmax : float, optional
            override kappa maximum. Default: None.
        colour : colourmap or str or array-like, optional
            colourmap or colourmap name; or key colour or min and max
            RGB colours to generate a colour map. Default: viridis.

        output : str, optional
            output filename to write to. Set to None to not write.
            Default: target-kl.hdf5.

        **kwargs
            keyword arguments passed to matplotlib.pyplot.pcolormesh.
            Defaults:

                rasterized: False

    Returns
    -------

        colourbar
            colourbar.
    """

    # data formatting

    data = tp.data.resolve.resolve(data, ['conductivity', 'seebeck',
                                   'electronic_thermal_conductivity'],
                                   direction=direction)
    data['zt'] = zt

    data = tp.calculate.kl_fromdict(data)

    if output is not None:
        tp.data.save.hdf5(data, output)

    # plotting

    cbar = add_heatmap(ax, data['temperature'], list(np.abs(data['doping'])),
                       data['lattice_thermal_conductivity'],
                       xinterp=xinterp, yinterp=yinterp, kind=kind,
                       yscale='log', xmin=xmin, xmax=xmax, ymin=ymin,
                       ymax=ymax, cmin=cmin, cmax=cmax, colour=colour,**kwargs)

    # axes formatting

    labels = tp.settings.labels()
    ax.set_xlabel(labels['temperature'])
    ax.set_ylabel(labels['doping'])
    cbar.set_label(labels['lattice_thermal_conductivity'])

    return cbar
