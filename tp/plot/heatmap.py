"""Heatmap plotters.

``add_heatmap`` is a base function, which handles axes limits, colourbar
formatting and extra colourmap options compared to ``matplotlib.pyplot.pcolormesh``.
The other functions plot specific quantites using this function.
"""

#Functions
#---------
#
#    add_heatmap:
#        heatmap.
#    add_pfmap:
#        power factor vs temperature and carrier concentration.
#    add_pfdiff:
#        power factor difference vs temperature and carrier concentration.
#    add_ztmap:
#        ZT vs temperature and carrier concentration.
#    add_ztdiff:
#        ZT difference vs temperature and carrier concentration.
#    add_kappa_target:
#        kappa_l needed to reach a given ZT.
#"""

from copy import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import tp
import warnings
import yaml
from scipy.interpolate import interp2d

warnings.filterwarnings('ignore', module='matplotlib')

try:
    filename = '{}/.config/tprc.yaml'.format(os.path.expanduser("~"))
    with open(filename, 'r') as f:
        conf = yaml.safe_load(f)
except yaml.parser.ParserError:
    warnings.warn('Failed to read ~/.config/tprc.yaml')
    conf = None
except FileNotFoundError:
    conf = None

def add_heatmap(ax, x, y, c, xinterp=None, yinterp=None, kind='linear',
                xscale='linear', yscale='linear', cscale='linear', xmin=None,
                xmax=None, ymin=None, ymax=None, cmin=None, cmax=None,
                colour='viridis', undercolour=None, overcolour=None,
                discrete=False, levels=None, contours=None,
                contourcolours='black', contourkwargs=None, **kwargs):
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

        colour : colourmap or str or array-like or dict, optional
            colourmap or colourmap name or highlight colour or
            highlight, min, max colours in that order, or dictionary
            with mid and min and/or max keys. Colour format must be
            hex or rgb (array) or a named colour recognised by
            matplotlib. Default: Blues.
        discrete : bool, optional
            use colour bands rather than a continuously shifting colour
            sequence. Default: False.
        levels : int or array-like, optional
            sets levels for discrete plots. Lists directly specify the
            contour levels while integers specify the maximum-1 number
            of "nice" levels to plot. Default: None.
        undercolour : str or array-like, optional
            colour for values under cmin. Default: None.
        overcolour : str or array-like, optional
            colour for values over cmax. Default: None.

        contours : int or float or array-like, optional
            contours to plot. Default: None.
        contourcolours : string or array-like, optional
            colours of the contours. If fewer are supplied, they
            repeat. Default: black.

        contourkwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.contour.
            Default: None
        kwargs
            keyword arguments passed to matplotlib.pyplot.pcolormesh or
            matplotlib.pyplot.contourf. Defaults are defined below,
            which are overridden by those in ``~/.config/tprc.yaml``,
            both of which are overridden by arguments passed to this
            function. Defaults:

                rasterized: True

    Returns
    -------

        colourbar
            colourbar.
    """

    # defaults

    defkwargs = {'rasterized': True}

    if conf is None or 'heatmap_kwargs' not in conf or \
       conf['heatmap_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['heatmap_kwargs'], **kwargs}

    # data trimming

    x = np.array(x)
    y = np.array(y)
    c = np.array(c)
    if xmin is None: xmin = np.nanmin(x)
    if xmax is None: xmax = np.nanmax(x)
    if ymin is None: ymin = np.nanmin(y)
    if ymax is None: ymax = np.nanmax(y)
    xi = np.where((x >= xmin) & (x <= xmax))[0]
    yi = np.where((y >= ymin) & (y <= ymax))[0]
    y = y[yi]
    x = x[xi]

    try:
        c = c[np.ix_(xi, yi)]
    except IndexError:
        c = c[np.ix_(xi[:-1], yi[:-1])]

    # colour
    # Sets colourbar extension.

    extend = 'neither'
    if cmin is None:
        cmin = np.nanmin(c)
    elif cmin > np.nanmin(c):
        extend = 'min'
    if cmax is None:
        cmax = np.nanmax(c)
    elif cmax < np.nanmax(c):
        if extend == 'min':
            extend = 'both'
        else:
            extend = 'max'

    if cscale == 'linear':
        cnorm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
    elif cscale == 'log':
        cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    # Tries to read as a colourmap or colourmap name, or uses a single
    # #rrggbb colour as the highlight colour for a tp.plot.colour.uniform.

    try:
        try:
            colours = copy(mpl.cm.get_cmap(colour))
        except AttributeError:
            colours = copy(mpl.colormaps[colour])
    except ValueError:
        if isinstance(colour, mpl.colors.ListedColormap):
            colours = colour
        elif isinstance(colour, str):
            colours = tp.plot.colour.uniform(colour)
        elif isinstance(colour, list):
            colours = tp.plot.colour.uniform(*colour)
        elif isinstance(colour, dict):
            colours = tp.plot.colour.uniform(**colour)
        else:
            raise Exception('colour must be a colourmap, colourmap '
                            'name, single highlight colour or '
                            'highlight, min, max colours in that '
                            'order, or a dictionary with mid and min '
                            'and/or max keys. Colour format must be '
                            'hex or rgb (array) or a named colour '
                            'recognised by matplotlib.')
    if undercolour is not None:
        colours.set_under(undercolour)
    if overcolour is not None:
        colours.set_over(overcolour)

    # data interpolation

    if xinterp is not None or yinterp is not None:
        cinterp = interp2d(x, y, np.transpose(c), kind=kind)
        if xinterp is not None:
            if xscale == 'linear': x = np.linspace(x[0], x[-1], xinterp)
            if xscale == 'log': x = np.geomspace(x[0], x[-1], xinterp)
        if yinterp is not None:
            if yscale == 'linear': y = np.linspace(y[0], y[-1], yinterp)
            if yscale == 'log': y = np.geomspace(y[0], y[-1], yinterp)
        x = np.array(x)
        y = np.array(y)
        c = cinterp(x, y)
    else:
        c = np.transpose(c)

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
    if discrete:
        kwargs.pop('rasterized')
        if levels in [None, (), [None], (None,)]:
            levels = 10
        else:
            if isinstance(levels, (list, np.ndarray, tuple)) and len(levels) == 1:
                levels = int(levels[0])
            elif isinstance(levels, float):
                raise TypeError(
                                'levels must be an int or array-like, not float')
        heat = ax.contourf(x[:-1], y[:-1], c, cmap=colours, norm=cnorm,
                           levels=levels, **kwargs)
    else:
        heat = ax.pcolormesh(x, y, c, cmap=colours, norm=cnorm, **kwargs)
    fig = ax.get_figure()
    if 'dos' in fig.__dict__ and fig.__dict__['dos']:
        # place colourbar outside dos
        cbar = plt.colorbar(heat, extend=extend)
    else:
        cbar = plt.colorbar(heat, ax=ax, extend=extend)

    # contours
    if contours not in [None, [None], (), (None,)]:
        if contourkwargs == None:
            contourkwargs = {}
        if not isinstance(contours, (list, np.ndarray)):
            contours = list(contours)
        if not isinstance(contourcolours, (list, np.ndarray)):
            contourcolours = list(contourcolours)
        if len(contours) > 1:
            ctnorm = [(ct-np.amin(contours))/
                      (np.amax(contours)-np.amin(contours))
                       for ct in contours]
        else:
            ctnorm = 0
        cmap = None
        try:
            cmap = contourcolours
            try:
                contourcolours = mpl.cm.get_cmap(contourcolours)(ctnorm)
            except AttributeError:
                contourcolours = mpl.colormaps[contourcolours](ctnorm)
            plt.contour(x[:-1], y[:-1], c, contours, cmap=cmap,
                        **contourkwargs)
        except ValueError:
            cmap = None
            if  isinstance(contourcolours, mpl.colors.ListedColormap):
                cmap = contourcolours
                contourcolours = contourcolours(ctnorm)
                plt.contour(x[:-1], y[:-1], c, contours, cmap=cmap,
                            **contourkwargs)
            elif not isinstance(contourcolours, (list, np.ndarray)):
                contourcolours = [contourcolours]
            while len(contourcolours) < len(contours):
                contourcolours += contourcolours
        if cmap is None:
            plt.contour(x[:-1], y[:-1], c, contours, colors=contourcolours,
                        **contourkwargs)
        for i, c in enumerate(contours):
            cbar.ax.axhline(c, color=contourcolours[i], **contourkwargs)

    # axes formatting
    tp.plot.utilities.set_locators(ax, x=xscale, y=yscale)
    tp.plot.utilities.set_locators(cbar.ax, y=cscale)

    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[0], y[-1])

    return cbar

def add_pfmap(ax, data, direction='avg', xinterp=200, yinterp=200,
              kind='linear', xmin=None, xmax=None, ymin=None, ymax=None,
              cmin=None, cmax=None, colour='viridis', discrete=False,
              levels=None, contours=None, contourcolours='black',
              contourkwargs=None, **kwargs):
    """Convenience wrapper for plotting power factor heatmaps.

    Calculates power factor, plots and formats labels etc.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            dictionary containing temperature and doping, and either
            power_factor or conductivity and seebeck.

        direction : str, optional
            crystal direction, accepts x-z/ a-c or average/ avg.
            Default: average.

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
            RGB colours to generate a colour map. Colour format must be
            hex or rgb (array) or a named colour recognised by
            matplotlib. Default: viridis.
        discrete : bool, optional
            use colour bands rather than a continuously shifting colour
            sequence. Default: False.
        levels : int or array-like, optional
            sets levels for discrete plots. Lists directly specify the
            contour levels while integers specify the maximum-1 number
            of "nice" levels to plot. Default: None.

        contours : int or float or array-like, optional
            PF contours to plot. Default: None.
        contourcolours : string or array-like, optional
            colours of the PF contours. If fewer are supplied, they
            repeat. Default: black.

        contourkwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.contour.
            Default: None
        kwargs
            keyword arguments passed to matplotlib.pyplot.pcolormesh or
            matplotlib.pyplot.contourf. Defaults are defined below,
            which are overridden by those in ``~/.config/tprc.yaml``,
            both of which are overridden by arguments passed to this
            function. Defaults:

                rasterized: True

    Returns
    -------

        colourbar
            colourbar.
    """
    
    # defaults

    defkwargs = {'rasterized': True}

    if conf is None or 'pfmap_kwargs' not in conf or \
       conf['pfmap_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['pfmap_kwargs'], **kwargs}

    equants = ['conductivity', 'seebeck']

    # data formatting

    if 'power_factor' in data:
        data = tp.data.utilities.resolve(data, 'power_factor', direction=direction)
    else:
        data = tp.data.utilities.resolve(data, equants, direction=direction)
        data = tp.calculate.power_factor_fromdict(data, use_tprc=True)

    # plotting

    cbar = add_heatmap(ax, data['temperature'], list(np.abs(data['doping'])),
                       data['power_factor'], xinterp=xinterp, yinterp=yinterp,
                       kind=kind, yscale='log', xmin=xmin, xmax=xmax,
                       ymin=ymin, ymax=ymax, cmin=cmin, cmax=cmax,
                       colour=colour, discrete=discrete, levels=levels,
                       contours=contours, countourcolours=contourcolours,
                       contourkwargs=contourkwargs, **kwargs)

    # axes formatting

    labels = tp.settings.labels()
    ax.set_xlabel(labels['temperature'])
    ax.set_ylabel(labels['doping'])
    cbar.set_label(labels['power_factor'])

    return cbar

def add_pfdiff(ax, data1, data2, direction='avg', xinterp=200, yinterp=200,
               kind='linear', xmin=None, xmax=None, ymin=None, ymax=None,
               cmin=None, cmax=None, colour1='#800080', colour2='#FF8000',
               midcolour='#FFFFFF', label1=None, label2=None, discrete=False,
               levels=None, contours=None, contourcolours='black',
               contourkwargs=None, **kwargs):
    """Plots a difference of two power factors heatmap.

    Calculates power factor, plots and formats labels etc.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data(1,2) : dict
            dictionaries containing temperature and doping, and either
            power_factor or conductivity and seebeck. The result is
            data1 - data2.

        direction : str, optional
            crystal direction, accepts x-z/ a-c or average/ avg.
            Default: average.

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
        colour(1,2) : str, optional
            colours for data(1,2). Colour format must be hex or rgb
            (array) or a named colour recognised by matplotlib.
            Deault: #800080, #FF8000.
        midcolour : str, optional
            colour at 0 difference. Default: #FFFFFF.
        label(1,2) : str, optional
            labels for data(1,2) to go in a legend. Default: None.
        discrete : bool, optional
            use colour bands rather than a continuously shifting colour
            sequence. Default: False.
        levels : int or array-like, optional
            sets levels for discrete plots. Lists directly specify the
            contour levels while integers specify the maximum-1 number
            of "nice" levels to plot. Default: None.

        contours : int or float or array-like, optional
            PF contours to plot. Default: None.
        contourcolours : string or array-like, optional
            colours of the PF contours. If fewer are supplied, they
            repeat. Default: black.

        contourkwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.contour.
            Default: None
        kwargs
            keyword arguments passed to matplotlib.pyplot.pcolormesh or
            matplotlib.pyplot.contourf. Defaults are defined below,
            which are overridden by those in ``~/.config/tprc.yaml``,
            both of which are overridden by arguments passed to this
            function. Defaults:

                rasterized: True

    Returns
    -------

        colourbar
            colourbar.
        list
            legend handles.
        list
            legend labels.
    """

    # defaults

    defkwargs = {'rasterized': True}

    if conf is None or 'pfdiff_kwargs' not in conf or \
       conf['pfdiff_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['pfdiff_kwargs'], **kwargs}

    equants = ['conductivity', 'seebeck']

    # data formatting

    data = [data1, data2]
    for i in [0, 1]:
        data[i]['doping'] = np.abs(data[i]['doping'])
        if 'power_factor' in data[i]:
            data[i] = tp.data.utilities.resolve(data[i], 'power_factor',
                                              direction=direction)
        else:
            data[i] = tp.data.utilities.resolve(data[i], equants,
                                              direction=direction)
            data[i] = tp.calculate.power_factor_fromdict(data[i], use_tprc=True)

    data[0], data[1] = tp.calculate.interpolate(*data, dependent='temperature',
                                                keys1='power_factor',
                                                keys2='power_factor')
    data[0], data[1] = tp.calculate.interpolate(*data, dependent='doping',
                                                keys1='power_factor', axis1=1,
                                                keys2='power_factor', axis2=1)
    diff = np.subtract(data[0]['power_factor'], data[1]['power_factor'])
    dmin, dmax = np.nanmin(diff), np.nanmax(diff)
    mid = dmin / (dmin - dmax)

    # colours

    # legend handles and labels
    h = [mpl.patches.Patch(facecolor=colour1, label=label1),
         mpl.patches.Patch(facecolor=colour2, label=label2)]
    l = [label1, label2]

    # if one is uniformly less than the other, its colour doesn't show
    if dmin > 0:
        colour2 = midcolour
    elif dmax < 0:
        colour1 = midcolour

    colour = tp.plot.colour.elbow(cmin=colour2, cmax=colour1, cmid=midcolour,
                                  midpoint=mid)
                            
    # plotting

    cbar = add_heatmap(ax, data[0]['temperature'], list(data[0]['doping']),
                       diff, xinterp=xinterp, yinterp=yinterp, kind=kind,
                       yscale='log', xmin=xmin, xmax=xmax, ymin=ymin,
                       ymax=ymax, cmin=cmin, cmax=cmax, colour=colour,
                       discrete=discrete, levels=levels, contours=contours,
                       contourcolours=contourcolours,
                       contourkwargs=contourkwargs, **kwargs)

    # axes formatting

    labels = tp.settings.labels()
    ax.set_xlabel(labels['temperature'])
    ax.set_ylabel(labels['doping'])
    cbar.set_label('$\Delta$ '+labels['power_factor'])

    return cbar, h, l

def add_ztmap(ax, data, kdata=1., direction='avg', xinterp=200,
              yinterp=200, kind='linear', xmin=None, xmax=None, ymin=None,
              ymax=None, cmin=None, cmax=None, colour='viridis',
              discrete=False, levels=None, contours=None,
              contourcolours='black', contourkwargs=None, **kwargs):
    """Convenience wrapper for plotting ZT heatmaps.

    Calculates ZT, plots and formats labels etc.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            dictionary containing temperature and doping, and either zt
            or conductivity, seebeck and electronic_thermal_conductivity.

        kdata : dict or int or float, optional
            dictionary containing lattice_thermal_conductivity and
            temperature. Ignored if zt in data. Can specify a constant
            value instead. Default: 1.
        direction : str, optional
            crystal direction, accepts x-z/ a-c or average/ avg.
            Default: average.

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
            colours to generate a colour map. Colour format must be hex
            or rgb (array) or a named colour recognised by matplotlib.
            Default: viridis.
        discrete : bool, optional
            use colour bands rather than a continuously shifting colour
            sequence. Default: False.
        levels : int or array-like, optional
            sets levels for discrete plots. Lists directly specify the
            contour levels while integers specify the maximum-1 number
            of "nice" levels to plot. Default: None.

        contours : int or float or array-like, optional
            ZT contours to plot. Default: None.
        contourcolours : string or array-like, optional
            colours of the ZT contours. If fewer are supplied, they
            repeat. Default: black.

        contourkwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.contour.
            Default: None
        kwargs
            keyword arguments passed to matplotlib.pyplot.pcolormesh or
            matplotlib.pyplot.contourf. Defaults are defined below,
            which are overridden by those in ``~/.config/tprc.yaml``,
            both of which are overridden by arguments passed to this
            function. Defaults:

                rasterized: True

    Returns
    -------

        colourbar
            colourbar.
    """
    
    # defaults

    defkwargs = {'rasterized': True}

    if conf is None or 'ztmap_kwargs' not in conf or \
       conf['ztmap_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['ztmap_kwargs'], **kwargs}

    ltc = 'lattice_thermal_conductivity'
    equants = ['conductivity', 'seebeck', 'electronic_thermal_conductivity']

    # data formatting

    if 'zt' in data:
        data = tp.data.utilities.resolve(data, 'zt', direction=direction)
    else:
        if kdata is None:
            kdata = 1.

        if not isinstance(kdata, (int, float)):
            data, kdata = tp.calculate.interpolate(data, kdata, 'temperature',
                                                   equants, ltc, kind='cubic')
            data[ltc] = kdata[ltc]
            data['meta']['dimensions'][ltc] = kdata['meta']['dimensions'][ltc]
            data['meta']['units'][ltc] = kdata['meta']['units'][ltc]
            data['meta']['kappa_source'] = kdata['meta']['kappa_source']

        else:
            data[ltc] = kdata * np.ones(len(data['temperature']))
            data['meta']['dimensions'][ltc] = ['temperature']
            data['meta']['units'][ltc] = tp.settings.units()[ltc]
            data['meta']['kappa_source'] = 'Set to {} {}'.format(kdata,
                                                   data['meta']['units'][ltc])

        data = tp.data.utilities.resolve(data, [*equants, ltc], direction=direction)
        data = tp.calculate.zt_fromdict(data, use_tprc=True)

    # plotting

    cbar = add_heatmap(ax, data['temperature'], list(np.abs(data['doping'])),
                       data['zt'], xinterp=xinterp, yinterp=yinterp, kind=kind,
                       yscale='log', xmin=xmin, xmax=xmax, ymin=ymin,
                       ymax=ymax, cmin=cmin, cmax=cmax, colour=colour,
                       discrete=discrete, levels=levels, contours=contours,
                       contourcolours=contourcolours,
                       contourkwargs=contourkwargs, **kwargs)

    # axes formatting

    labels = tp.settings.labels()
    ax.set_xlabel(labels['temperature'])
    ax.set_ylabel(labels['doping'])
    cbar.set_label(labels['zt'])

    return cbar

def add_ztdiff(ax, data1, data2, kdata1=1., kdata2=1., direction='avg',
               xinterp=200, yinterp=200, kind='linear', xmin=None, xmax=None,
               ymin=None, ymax=None, cmin=None, cmax=None, colour1='#800080',
               colour2='#FF8000', midcolour='#FFFFFF', label1=None,
               label2=None, discrete=False, levels=None, contours=None,
               contourcolours='black', contourkwargs=None, **kwargs):
    """Plots a difference of two ZTs heatmap.

    Calculates ZT, plots and formats labels etc.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data(1,2) : dict
            dictionaries containing temperature and doping, and either
            zt or conductivity, seebeck and electronic_thermal_conductivity.
            The result is data1 - data2.

        kdata(1,2) : dict or int or float, optional
            dictionaries containing lattice_thermal_conductivity and
            temperature. Ignored if zt in data. Can specify constant
            values instead. Default: 1.
        direction : str, optional
            crystal direction, accepts x-z/ a-c or average/ avg.
            Default: average.

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
        colour(1,2) : str, optional
            colours for data(1,2). Colour format must be hex or rgb
            (array) or a named colour recognised by matplotlib.
            Default: #800080, #FF8000.
        midcolour : str, optional
            colour at 0 difference. Default: #FFFFFF.
        label(1,2) : str, optional
            labels for data(1,2) to go in a legend. Default: None.
        discrete : bool, optional
            use colour bands rather than a continuously shifting colour
            sequence. Default: False.
        levels : int or array-like, optional
            sets levels for discrete plots. Lists directly specify the
            contour levels while integers specify the maximum-1 number
            of "nice" levels to plot. Default: None.

        contours : int or float or array-like, optional
            ZT contours to plot. Default: None.
        contourcolours : string or array-like, optional
            colours of the ZT contours. If fewer are supplied, they
            repeat. Default: black.

        contourkwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.contour.
            Default: None
        kwargs
            keyword arguments passed to matplotlib.pyplot.pcolormesh or
            matplotlib.pyplot.contourf. Defaults are defined below,
            which are overridden by those in ``~/.config/tprc.yaml``,
            both of which are overridden by arguments passed to this
            function. Defaults:

                rasterized: True

    Returns
    -------

        colourbar
            colourbar.
        list
            legend handles.
        list
            legend labels.
    """

    # defaults

    defkwargs = {'rasterized': True}

    if conf is None or 'ztdiff_kwargs' not in conf or \
       conf['ztdiff_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['ztdiff_kwargs'], **kwargs}

    ltc = 'lattice_thermal_conductivity'
    equants = ['conductivity', 'seebeck', 'electronic_thermal_conductivity']

    # data formatting

    data = [data1, data2]
    kdata = [kdata1, kdata2]
    for i in [0, 1]:
        data[i]['doping'] = np.abs(data[i]['doping'])
        if 'zt' in data[i]:
            data[i] = tp.data.utilities.resolve(data[i], 'zt',
                                                direction=direction)
        else:
            if kdata[i] is None:
                kdata[i] = 1.
            if not isinstance(kdata[i], (int, float)):
                data[i], kdata[i] = tp.calculate.interpolate(data[i], kdata[i],
                                                             'temperature',
                                                             equants, ltc,
                                                             kind='cubic')
                data[i][ltc] = kdata[i][ltc]
                data[i]['meta']['dimensions'][ltc] = kdata[i]['meta']['dimensions'][ltc]
                data[i]['meta']['units'][ltc] = kdata[i]['meta']['units'][ltc]
                data[i]['meta']['kappa_source'] = kdata[i]['meta']['kappa_source']

            else:
                data[i][ltc] = kdata[i] * np.ones(len(data[i]['temperature']))
                data[i]['meta']['dimensions'][ltc] = ['temperature']
                data[i]['meta']['units'][ltc] = tp.settings.units()[ltc]
                data[i]['meta']['kappa_source'] = 'Set to {} {}'.format(kdata[i],
                                                   data[i]['meta']['units'][ltc])
                        
            data[i] = tp.data.utilities.resolve(data[i], [*equants, ltc],
                                                direction=direction)
            data[i] = tp.calculate.zt_fromdict(data[i], use_tprc=True)

    data[0], data[1] = tp.calculate.interpolate(*data, dependent='temperature',
                                                keys1='zt', keys2='zt')
    data[0], data[1] = tp.calculate.interpolate(*data, dependent='doping',
                                                keys1='zt', axis1=1,
                                                keys2='zt', axis2=1)
    diff = np.subtract(data[0]['zt'], data[1]['zt'])
    dmin, dmax = np.nanmin(diff), np.nanmax(diff)
    mid = dmin / (dmin - dmax)

    # colours

    # legend handles and labels
    h = [mpl.patches.Patch(facecolor=colour1, label=label1),
         mpl.patches.Patch(facecolor=colour2, label=label2)]
    l = [label1, label2]

    # if one is uniformly less than the other, its colour doesn't show
    if dmin > 0:
        colour2 = midcolour
    elif dmax < 0:
        colour1 = midcolour

    colour = tp.plot.colour.elbow(cmin=colour2, cmax=colour1, cmid=midcolour,
                                  midpoint=mid)
                            
    # plotting

    cbar = add_heatmap(ax, data[0]['temperature'], list(data[0]['doping']),
                       diff, xinterp=xinterp, yinterp=yinterp, kind=kind,
                       yscale='log', xmin=xmin, xmax=xmax, ymin=ymin,
                       ymax=ymax, cmin=cmin, cmax=cmax, colour=colour,
                       discrete=discrete, levels=levels, contours=contours,
                       contourcolours=contourcolours,
                       contourkwargs=contourkwargs, **kwargs)

    # axes formatting

    labels = tp.settings.labels()
    ax.set_xlabel(labels['temperature'])
    ax.set_ylabel(labels['doping'])
    cbar.set_label('$\Delta$ '+labels['zt'])

    return cbar, h, l
    
def add_kappa_target(ax, data, zt=2, direction='avg', xinterp=200,
                     yinterp=200, kind='linear', xmin=None, xmax=None,
                     ymin=None, ymax=None, cmin=0, cmax=None, colour='viridis',
                     negativecolour='grey', discrete=False, levels=None,
                     contours=None, contourcolours='black', contourkwargs=None,
                     **kwargs):
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
            crystal direction, accepts x-z/ a-c or average/ avg.
            Default: average.

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
            RGB colours to generate a colour map. Colour format must be
            hex or rgb (array) or a named colour recognised by
            matplotlib. Default: viridis.
        negativecolour : str or array-like, optional
            colour for values under cmin. Default: grey.
        discrete : bool, optional
            use colour bands rather than a continuously shifting colour
            sequence. Default: False.
        levels : int or array-like, optional
            sets levels for discrete plots. Lists directly specify the
            contour levels while integers specify the maximum-1 number
            of "nice" levels to plot. Default: None.

        contours : int or float or array-like, optional
            kappa contours to plot. Default: None.
        contourcolours : string or array-like, optional
            colours of the kappa contours. If fewer are supplied, they
            repeat. Default: black.

        contourkwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.contour.
            Default: None
        kwargs
            keyword arguments passed to matplotlib.pyplot.pcolormesh or
            matplotlib.pyplot.contourf. Defaults are defined below,
            which are overridden by those in ``~/.config/tprc.yaml``,
            both of which are overridden by arguments passed to this
            function. Defaults:

                rasterized: True

    Returns
    -------

        colourbar
            colourbar.
    """

    # defaults

    defkwargs = {'rasterized': True}

    if conf is None or 'kappa_target_kwargs' not in conf or \
       conf['kappa_target_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['kappa_target_kwargs'], **kwargs}

    ltc = 'lattice_thermal_conductivity'
    equants = ['conductivity', 'seebeck', 'electronic_thermal_conductivity']
    # data formatting

    data['zt'] = zt

    data = tp.calculate.kl_fromdict(data, use_tprc=True)
    data = tp.data.utilities.resolve(data, ltc, direction=direction)

    # plotting

    cbar = add_heatmap(ax, data['temperature'], list(np.abs(data['doping'])),
                       data[ltc], xinterp=xinterp, yinterp=yinterp, kind=kind,
                       yscale='log', xmin=xmin, xmax=xmax, ymin=ymin,
                       ymax=ymax, cmin=cmin, cmax=cmax, colour=colour,
                       undercolour=negativecolour, discrete=discrete,
                       levels=levels, contours=contours,
                       contourcolours=contourcolours,
                       contourkwargs=contourkwargs, **kwargs)

    # axes formatting

    labels = tp.settings.labels()
    ax.set_xlabel(labels['temperature'])
    ax.set_ylabel(labels['doping'])
    cbar.set_label(labels[ltc])

    return cbar
