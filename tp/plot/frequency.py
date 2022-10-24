"""Functions which plot phonon frequency on the x-axis.

Each function takes a set of axes and data dictionary as its main input,
and and an ``invert`` argument to plot sideways by a phonon dispersion,
and a main argument, which determines if axes limits are set and
sometimes how scaling is handled, which helps if plotting multiple
quantities on the same axes.

Functions
---------

    add_dos:
        phonon dos.
    add_cum_kappa:
        cumulative kappa vs frequency.
    add_waterfall:
        scatter plots of various values.
    add_projected_waterfall:
        waterfall, but with a second quantity projected.
    add_density:
        density of phonon modes for a property vs frequency.

    format_waterfall:
        formatting for the waterfall and density plots.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import tp
import warnings
import yaml
from scipy.stats import gaussian_kde

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

def add_dos(ax, data, sigma=None, projected=True, total=False,
            totallabel='Total', main=True, invert=False, scale=False,
            colour='tab10', totalcolour=None, fill=True, fillalpha=0.2,
            line=True, linestyle='-', marker=None, **kwargs):
    """Adds a phonon density of states (DoS) to a set of axes.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            DoS data.
        sigma : float
            Gaussian broadening. 0.2 is a good place to start. Does not
            know if you've already broadened it. Default: None.

        projected : bool, optional
            plot atom-projected DoS. Default: True.
        total : bool, optional
            plot total DoS. Default: False
        totallabel : str, optional
            label for the total line. Other labels are taken directly
            from the input dictionary. Default: Total.

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            plot frequency on y axis. Default: False.
        scale : bool, optional
            if main, scale to percent. If not main, scale to axis
            limits. Default: False.

        colour : dict or list or str or colourmap, optional
            RGB colours per atom as a dictionary or a list in POSCAR
            order, with total as the last colour. Can instead provide a
            colourmap or colourmap name, which don't include a total
            colour. If not projected, can specify a single total colour.
            Total colour is overridden by totalcolour. Default: tab10.
        totalcolour : str, optional
            colour for the total line. Overrides specifying as part of
            colour. Default: black.
        fill : bool, optional
            fill below lines. Default: True.
        fillalpha : float, optional
            fill alpha scaled to 0-1. Default: 0.2.
        line : bool, optional
            plot lines. Default: True.
        linestyle : dict or list or str, optional
            linestyle or dictionary of linestyles per atom and total, or
            list in POSCAR order, with total as the last linestyle.
            Default: -.
        marker : dict or list or str or tuple, optional
            marker or dictionary of markers per atom and total, or list
            in POSCAR order, with total as the last marker.
            Default: None.

        kwargs
            keyword arguments passed to matplotlib.pyplot.plot, unless
            line=False, in which case they are passed to fill_between
            instead. rasterized is always passed to both.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                rasterized: False

    Returns
    -------

        None
            adds plot directly to ax.
    """

    # defaults

    defkwargs = {'rasterized': False}

    if conf is None or 'dos_kwargs' not in conf or conf['dos_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['dos_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['total', 'main', 'invert', 'scale', 'fill', 'line'],
                           [ total,   main,   invert,   scale,   fill,   line]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)
    assert fill or line, 'fill or line or both must be True.'
    assert projected or total, 'projected or total or both must be True.'
    assert isinstance(fillalpha, (float, int)) and fillalpha >= 0 \
                                               and fillalpha <= 1, \
           'fillalpha must be a float/ integer between 0 and 1.'

    # data scaling

    data = dict(data)
    f = data.pop('frequency')
    if 'meta' in data: del data['meta']
    if not total: del data['total']

    # find the true ends of the data and truncate to there
    n = int(np.ceil(0.9*len(f)))
    for key in data:
        while n < len(f) - 1 and data[key][n] != 0:
            n += 1
    n += 1
    f = f[:n]
    for key in data:
        data[key] = data[key][:n]

    n = int(np.ceil(0.1*len(f)))
    for key in data:
        while n > 0 and data[key][n] != 0:
            n -= 1
    f = f[n:]
    for key in data:
        data[key] = data[key][n:]

    fmin, fmax = f[0], f[-1]
    fstep = f[1] - f[0]

    if sigma is not None:
        # extend the data to 3 sigma above and below
        while f[-1] < fmax + 3 * sigma:
            f = np.append(f, f[-1] + fstep)

        while f[0] > fmin - 3 * sigma:
            f = np.insert(f, 0, f[0] - fstep)

        data2 = {}
        for key in data:
            data2[key] = np.zeros(len(f))
            for i, d in enumerate(data[key]):
                data2[key] += d * tp.calculate.gaussian(f, f[i], sigma)
        data = data2

    if scale:
        axscale = [0, 100] if main else None
        axis = 'x' if invert else 'y'
        data = tp.plot.utilities.scale_to_axis(ax, data, scale=axscale,
                                               axis=axis)
    if total:
        totaldata = data.pop('total')

    # colours
    # If projected, tries to read the colour as a colourmap name, then colourmap
    # object, then list of colours, then dictionary, and converts it
    # into a dictionary if necessary. If not projected, colour can be a single
    # colour for the total. The total colour is found first from totalcolour,
    # otherwise it can be specified in the dictionary under total or as the
    # last list item, or defaults to black.

    if not projected:
        if totalcolour is not None:
            colour = {'total': totalcolour}
        elif isinstance(colour, dict) and 'total' in colour:
            pass
        elif isinstance(colour, list):
            if len(colour) == 1:
                colour = {'total': colour[0]}
            elif len(colour) > len(data):
                colour = {'total': colour[len(data)]}
        else:
            colour = {'total': colour}
    else:
        if not isinstance(colour, dict):
            if not isinstance(colour, list):
                try:
                    cmap = plt.cm.get_cmap(colour)(np.linspace(0,1,len(data)))
                except ValueError:
                    cmap = colour(np.linspace(0,1,len(data)))
            else:
                cmap = colour
            colour = {}
            for i, c in enumerate(data):
                colour[c] = cmap[i]
            if len(cmap) > len(data):
                colour['total'] = cmap[len(data)]

        if total and 'total' not in colour:
            if totalcolour is not None:
                colour['total'] = totalcolour
            else:
                colour['total'] = '#000000'

    markers = {}
    if isinstance(marker, (str, tuple)) or marker is None:
        if not projected:
            markers['total'] = marker
        else:
            marker = [marker]
    if isinstance(marker, list):
        while len(marker) < len(data):
            marker.append(marker[-1])
        for i, c in enumerate(data):
            markers[c] = marker[i]
        if len(marker) > len(data):
            markers['total'] = marker[len(data)+1]
    elif isinstance(marker, dict):
        markers = marker
    else:
        raise Exception('marker must be a dict, list, str, tuple or None')
    if total not in markers:
        markers['total'] = None

    linestyles = {}
    if isinstance(linestyle, str):
        if not projected:
            linestyles['total'] = linestyle 
        else:
            linestyle = [linestyle]
    if isinstance(linestyle, list):
        while len(linestyle) < len(data):
            linestyle.append(linestyle[-1])
        for i, c in enumerate(data):
            linestyles[c] = linestyle[i]
        if len(linestyle) > len(data):
            linestyles['total'] = linestyle[len(data)+1]
    if isinstance(linestyle, dict):
        linestyles = linestyle
    if total not in linestyles:
        linestyles['total'] = '-'

    # plotting

    if total:
        if invert:
            if fill and line:
                ax.fill_between(totaldata, f, facecolor=colour['total'],
                                linewidth=0, alpha=fillalpha,
                                rasterized=kwargs['rasterized'])
                ax.plot(totaldata, f, label=totallabel, color=colour['total'],
                        linestyle=linestyles['total'], marker=markers['total'],
                        **kwargs)
            elif fill and not line:
                ax.fill_between(totaldata, f, label=totallabel, linewidth=0,
                                facecolor=colour['total'], alpha=fillalpha,
                                rasterized=kwargs['rasterized'])
            else:
                ax.plot(totaldata, f, label=totallabel, color=colour['total'],
                        linestyle=linestyles['total'], marker=markers['total'],
                        **kwargs)
        else:
            if fill and line:
                ax.fill_between(f, totaldata, facecolor=colour['total'],
                                linewidth=0, alpha=fillalpha,
                                rasterized=kwargs['rasterized'])
                ax.plot(f, totaldata, label=totallabel, color=colour['total'],
                        linestyle=linestyles['total'], marker=markers['total'],
                        **kwargs)
            elif fill and not line:
                ax.fill_between(f, totaldata, label=totallabel,
                                facecolor=colour['total'], linewidth=0,
                                alpha=fillalpha, **kwargs)
            else:
                ax.plot(f, totaldata, label=totallabel, color=colour['total'],
                        linestyle=linestyles['total'], marker=markers['total'],
                        **kwargs)
    if projected:
        for key in data:
            if invert:
                if fill and line:
                    ax.fill_between(data[key], f, facecolor=colour[key],
                                    linewidth=0, alpha=fillalpha)
                    ax.plot(data[key], f, label='${}$'.format(key),
                            color=colour[key], linestyle=linestyles[key],
                            marker=markers[key], **kwargs)
                elif fill and not line:
                    ax.fill_between(data[key], f, label='${}$'.format(key),
                                    facecolor=colour[key], linewidth=0,
                                    alpha=fillalpha)
                else:
                    ax.plot(data[key], f, label='${}$'.format(key),
                            color=colour[key], linestyle=linestyles[key],
                            marker=markers[key], **kwargs)
            else:
                if fill and line:
                    ax.fill_between(f, data[key], facecolor=colour[key],
                                    linewidth=0, alpha=fillalpha)
                    ax.plot(f, data[key], color=colour[key],
                            label='${}$'.format(key), linestyle=linestyles[key],
                            marker=markers[key], **kwargs)
                elif fill and not line:
                    ax.fill_between(f, data[key], label='${}$'.format(key),
                                    facecolor=colour[key], linewidth=0,
                                    alpha=fillalpha)
                else:
                    ax.plot(f, data[key], label='${}$'.format(key),
                            color=colour[key], linestyle=linestyles[key],
                            marker=markers[key], **kwargs)

    # axes formatting

    if sigma is not None:
        fmax += sigma
    else:
        fmax *= 1.01

    if main:
        if invert:
            axlabels = tp.settings.inverted_labels()
            pad = mpl.rcParams['xtick.major.pad']
            ax.set_xlabel(axlabels['dos'], labelpad=pad)
            ax.tick_params(axis='y', labelleft=False)
            ax.set_xlim(left=0)
            ax.set_ylim(top=fmax)
            if round(fmin, 1) >= 0.:
                ax.set_ylim(bottom=0)
            tp.plot.utilities.set_locators(ax, x='null', y='linear')
        else:
            axlabels = tp.settings.labels()
            pad = mpl.rcParams['ytick.major.pad']
            ax.set_xlabel(axlabels['frequency'])
            ax.set_ylabel(axlabels['dos'], labelpad=pad)
            ax.set_ylim(bottom=0)
            ax.set_xlim(right=fmax)
            if round(fmin, 1) >= 0.:
                ax.set_xlim(left=0)
            tp.plot.utilities.set_locators(ax, y='null', x='linear')

    return

def add_cum_kappa(ax, data, temperature=300, direction='avg', label=None,
                  main=True, invert=False, scale=False, colour=None,
                  fill=False, fillcolour=0.2, line=True, linestyle='-',
                  marker=None, verbose=False, **kwargs):
    """Cumulates and plots kappa against frequency.

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
        data : dict or list
            (list of sets of) Phono3py data including:

        temperature : float, optional
            temperature in K (finds nearest). Default: 300.

                mode_kappa: array-like
                    frequency and q-point decomposed lattice thermal
                    conductivity.
                frequency : array-like
                    frequencies.
                temperature : array-like
                    temperature.

        direction : str or list, optional
            (list of) direction(s) from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.
        label : str, optional
            (list of) legend label(s). Defaults to $\mathregular{\kappa_l}$
            if there is only one line plotted, or direction if there are
            more.

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            plot frequency on y axis. Default: False.
        scale : bool, optional
            if main, scale to percent, else scale to axis limits.
            Default: True.

        colour : str or list, optional
            (list of) hex or named line colour(s).
            Default: default colour cycle.
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

        verbose : bool, optional
            Write actual temperature used if applicable.
            Default: False.
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

    # defaults

    defkwargs = {'rasterized': False}

    if conf is None or 'frequency_cum_kappa_kwargs' not in conf or \
       conf['frequency_cum_kappa_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['frequency_cum_kappa_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['main', 'invert', 'scale', 'fill', 'line'],
                           [ main,   invert,   scale,   fill,   line]):
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

    fmax, kmax = None, None
    i = 0

    for dat in data:
        for d in direction:
            data2 = tp.data.utilities.resolve(dat, 'mode_kappa',
                                            temperature=temperature,
                                            direction=d)
            k = np.ravel(data2['mode_kappa'])
            f = np.ravel(data2['frequency'])

            f, k = tp.calculate.cumulate(f, k)

            if fmax is None or fmax < f[-1]:
                fmax = f[-1]
            if kmax is None or kmax < k[-1]:
                kmax = 100 if main and scale else k[-1]

            f = np.append(f, 100*f[-1])
            k = np.append(k, k[-1])

            if scale:
                axscale = [0, 100] if main else None
                axis = 'x' if invert else 'y'
                k = tp.plot.utilities.scale_to_axis(ax, k, scale=axscale, axis=axis)

            colour1 = colour[i % len(colour)]
            fillcolour1 = fillcolour[i % len(fillcolour)]
            linestyle1 = linestyle[i % len(linestyle)]
            marker1 = marker[i % len(marker)]
            label1 = label[i % len(label)]

            # colour
            # Tries to read the colour as an rgb code, then alpha value.

            if fill:
                try:
                    fillcolour2 = mpl.colors.to_rgba(colour1, fillcolour1)
                except ValueError:
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
                if invert:
                    ax.plot(k, f, color=colour1, linestyle=linestyle1,
                            marker=marker1, label=label1, **kwargs)
                    ax.fill_between(k, f, facecolor=fillcolour2, linewidth=0)
                else:
                    ax.plot(f, k, color=colour1, linestyle=linestyle1,
                            marker=marker1, label=label1, **kwargs)
                    ax.fill_between(f, k, facecolor=fillcolour2)

            elif fill and not line:
                    if invert:
                        ax.fill_between(k, f, facecolor=fillcolour2, **kwargs)
                    else:
                        ax.fill_between(f, k, facecolor=fillcolour2, **kwargs)

            else:
                if invert:
                    ax.plot(k, f, color=colour1, linestyle=linestyle1,
                            marker=marker1, label=label1, **kwargs)
                else:
                    ax.plot(f, k, color=colour1, linestyle=linestyle1,
                            marker=marker1, label=label1, **kwargs)

            i += 1

        if verbose:
            print('Using {} {}.'.format(data2['meta']['temperature'],
                                        data2['meta']['units']['temperature']))

    # axes formatting

    if main:
        if invert:
            axlabels = tp.settings.inverted_labels()
            if scale:
                ax.set_xlabel(axlabels['cumulative_percent'])
            else:
                ax.set_xlabel(axlabels['cumulative_kappa'])
            ax.tick_params(axis='y', labelleft=False)
            ax.set_xlim(0, kmax)
            ax.set_ylim(0, fmax)

        else:
            axlabels = tp.settings.labels()
            if scale:
                ax.set_ylabel(axlabels['cumulative_percent'])
            else:
                ax.set_ylabel(axlabels['cumulative_kappa'])
            ax.set_xlabel(axlabels['frequency'])
            ax.set_ylim(0, kmax)
            ax.set_xlim(0, fmax)

        tp.plot.utilities.set_locators(ax, x='linear', y='linear')

    return

def add_waterfall(ax, data, quantity, xquantity='frequency', temperature=300,
                  direction='avg', main=True, invert=False, colour='viridis',
                  verbose=False, **kwargs):
    """Adds a waterfall plot of quantities against frequency.

    Has an option to change the x-quantity.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            data including frequency and quantity.
        quantity : str
            y-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path,
            mode_kappa, occupation or ph_ph_strength.

        xquantity : str, optional
            x-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path,
            mode_kappa, occupation or ph_ph_strength.
            Default: frequency.
        temperature : float, optional
            temperature in K. Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            invert x- and y-axes. Default: False.

        colour : colourmap or str or array-like or dict, optional
            colourmap or colourmap name or list of colours (one for
            each band or one for each point) or two colours for a
            linear colourmap (required hex or rgb or named colours) or
            a dictionary with cmin and cmax keys or a single colour.
            Default: viridis.

        verbose : bool, optional
            Write actual temperature used if applicable.
            Default: False.
        kwargs
            keyword arguments passed to matplotlib.pyplot.scatter.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                alpha:      0.3
                s:          2
                edgecolors: black
                linewidth:  0
                marker:     '.'
                rasterized: True

    Returns
    -------

        None
            adds plot directly to ax.
    """

    # defaults

    defkwargs = {'alpha':      0.3,
                 's':          2,
                 'edgecolors': 'black',
                 'linewidth':  0,
                 'marker':     '.',
                 'rasterized': True}

    if conf is None or 'waterfall_kwargs' not in conf or \
       conf['waterfall_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['waterfall_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['main', 'invert'],
                           [ main,   invert]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)

    # data formatting

    tnames = tp.settings.to_tp()
    if invert: quantity, xquantity = xquantity, quantity
    quantity = tnames[quantity] if quantity in tnames else quantity
    xquantity = tnames[xquantity] if xquantity in tnames else xquantity
    if quantity == 'lattice_thermal_conductivity': quantity = 'mode_kappa'
    if xquantity == 'lattice_thermal_conductivity': xquantity = 'mode_kappa'

    data = tp.data.utilities.resolve(data, [quantity, xquantity],
                                   temperature=temperature,
                                   direction=direction)
    if verbose and 'temperature' in data['meta']:
        print('Using {} {}.'.format(data['meta']['temperature'],
                                    data['meta']['units']['temperature']))
    x = np.ravel(data[xquantity])
    y = np.abs(np.ravel(data[quantity]))

    # colour
    # Tries to read as a colourmap name or colourmap object or list of
    # colours (of varying formats), one per band, and assigns
    # appropriately, or just two colours generate a linear colourmap, or
    # a dictionary with cmin and cmax keys. Otherwise leaves as is,
    # which is appropriate for a single colour or colour per point.

    s = np.shape(data[xquantity])
    try:
        colour = mpl.cm.get_cmap(colour)
        colours = [colour(i) for i in np.linspace(0, 1, s[1])]
        colours = np.tile(colours, (s[0], 1))
    except ValueError:
        if isinstance(colour, mpl.colors.ListedColormap):
            colours = [colour(i) for i in np.linspace(0, 1, s[1])]
            colours = np.tile(colours, (s[0], 1))
        elif isinstance(colour, str) and colour == 'skelton':
            colour = tp.plot.colour.skelton()
            colours = [colour(i) for i in np.linspace(0, 1, s[1])]
            colours = np.tile(colours, (s[0], 1))
        elif isinstance(colour, list) or isinstance(colour, np.ndarray):
            if len(colour) == s[1]:
                if np.ndim(colour) == 1:
                    colours = np.tile(colour, s[0])
                elif np.ndim(colour) == 2:
                    colours = np.tile(colour, (s[0], 1))
            elif len(colour) == 2:
                colour = tp.plot.colour.linear(cmin=colour[0], cmax=colour[1])
                colours = [colour(i) for i in np.linspace(0, 1, s[1])]
                colours = np.tile(colours, (s[0], 1))
        elif isinstance(colour, dict):
                colour = tp.plot.colour.linear(**colour)
                colours = [colour(i) for i in np.linspace(0, 1, s[1])]
                colours = np.tile(colours, (s[0], 1))
        else:
            colours = colour

    # plotting

    ax.scatter(x, y, facecolor=colours, **kwargs)

    # axes formatting

    if main:
        data[xquantity] = x
        data[quantity] = y
        format_waterfall(ax, data, quantity, xquantity)
        if invert:
            axlabels = tp.settings.inverted_labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.tick_params(axis='y', labelleft=False)
        else:
            axlabels = tp.settings.labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.set_ylabel(axlabels[quantity])

    return

def add_projected_waterfall(ax, data, quantity, projected,
                            xquantity='frequency', temperature=300,
                            direction='avg', main=True, invert=False,
                            colour='viridis', cmin=None, cmax=None,
                            cscale=None, unoccupied='grey', verbose=False,
                            **kwargs):
    """Adds a waterfall plot against frequency with a colour axis.

    Has an option to change the x-quantity.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            data including frequency and quantity.
        quantity : str
            y-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path,
            mode_kappa, occupation or ph_ph_strength.
        projected : str
            colour-axis quantity. Accepts frequency, gamma,
            group_velocity, gv_by_gv, heat_capacity, lifetime,
            mean_free_path, mode_kappa, occupation or ph_ph_strength.

        xquantity : str, optional
            x-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path,
            mode_kappa, occupation or ph_ph_strength.
            Default: frequency.
        temperature : float, optional
            temperature in K. Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            invert x- and y-axes. Default: False.

        colour : colourmap or str or array-like or dict, optional
            colourmap or colourmap name or highlight colour or
            highlight, min, max colours in that order, or dictionary
            with cmid and cmin and/or cmax keys. Colour format must be
            hex or rgb (array) or a named colour recognised by
            matplotlib. Default: viridis.
        cmin : float, optional
            colour scale minimum. Default: display 99 % data.
        cmax : float, optional
            colour scale maximum. Default: display 99.9 % data.
        cscale : str optional
            override colour scale (linear/ log). Default: None.
        unoccupied : str or array-like, optional
            if the colour variable is occupation, values below 1 are
            coloured in this colour. If set to None, or cmin is set,
            this feature is turned off. Default: grey.

        verbose : bool, optional
            Write actual temperature used if applicable.
            Default: False.
        kwargs
            keyword arguments passed to matplotlib.pyplot.scatter.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                alpha:      0.3
                s:          2
                edgecolors: black
                linewidth:  0
                marker:     .
                rasterized: True

    Returns
    -------

        colorbar
            colourbar for projected data.
    """

    from copy import copy

    # defaults

    defkwargs = {'alpha':      0.3,
                 's':          2,
                 'edgecolors': 'black',
                 'linewidth':  0,
                 'marker':     '.',
                 'rasterized': True}

    if conf is None or 'projected_waterall_kwargs' not in conf or \
       conf['projected_waterfall_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['projected_waterfall_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['main', 'invert'],
                           [ main,   invert]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)

    # data formatting

    if quantity == 'kappa': quantity = 'mode_kappa'
    if xquantity == 'kappa': xquantity = 'mode_kappa'
    if projected == 'kappa': projected = 'mode_kappa'
    tnames = tp.settings.to_tp()
    axlabels = tp.settings.labels()
    if invert: quantity, xquantity = xquantity, quantity
    quantity = tnames[quantity] if quantity in tnames else quantity
    xquantity = tnames[xquantity] if xquantity in tnames else xquantity
    projected = tnames[projected] if projected in tnames else projected

    data = tp.data.utilities.resolve(data, [quantity, xquantity, projected],
                                   temperature=temperature,
                                   direction=direction)
    if verbose and 'temperature' in data['meta']:
        print('Using {} {}.'.format(data['meta']['temperature'],
                                    data['meta']['units']['temperature']))
    x = np.ravel(data[xquantity])
    y = np.abs(np.ravel(data[quantity]))
    c = np.abs(np.ravel(data[projected]))

    cmap = tp.plot.utilities.parse_colours(colour)
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
        data[xquantity] = x
        data[quantity] = y
        format_waterfall(ax, data, quantity, xquantity)
        if invert:
            axlabels = tp.settings.inverted_labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.tick_params(axis='y', labelleft=False)
        else:
            axlabels = tp.settings.labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.set_ylabel(axlabels[quantity])

    return cbar

def add_density(ax, data, quantity, xquantity='frequency', temperature=300,
                direction='avg', main=True, invert=False, colour='Blues',
                verbose=False, **kwargs):
    """Adds a density plot of quantities against frequency.

    Has an option to change the x-quantity.

    Arguments
    ---------

        ax : axes
            axes to plot on.
        data : dict
            data including frequency and quantity.
        quantity : str
            y-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path,
            mode_kappa, occupation or ph_ph_strength.

        xquantity : str, optional
            x-axis quantity. Accepts frequency, gamma, group_velocity,
            gv_by_gv, heat_capacity, lifetime, mean_free_path,
            mode_kappa, occupation or ph_ph_strength.
            Default: frequency.
        temperature : float, optional
            temperature in K. Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c or
            average/ avg. Default: average.

        main : bool, optional
            set ticks, labels, limits. Default: True.
        invert : bool, optional
            invert x- and y-axes. Default: False.

        colour : colourmap or str or array-like or dict, optional
            colourmap or colourmap name or highlight colour or
            highlight, min, max colours in that order, or dictionary
            with cmid and cmin and/or cmax keys. Colour format must be
            hex or rgb (array) or a named colour recognised by
            matplotlib. Default: Blues.

        verbose : bool, optional
            Write actual temperature used if applicable.
            Default: False.
        kwargs
            keyword arguments passed to matplotlib.pyplot.scatter.
            Defaults are defined below, which are overridden by those in
            ``~/.config/tprc.yaml``, both of which are overridden by
            arguments passed to this function.
            Defaults:

                s:          1  
                rasterized: True

    Returns
    -------

        None
            adds plot directly to ax.
    """

    # defaults

    defkwargs = {'s':      20,
                 'rasterized': True}

    if conf is None or 'density_kwargs' not in conf or \
       conf['density_kwargs'] is None:
        kwargs = {**defkwargs, **kwargs}
    else:
        kwargs = {**defkwargs, **conf['density_kwargs'], **kwargs}

    # input checks

    for name, value in zip(['main', 'invert'],
                           [ main,   invert]):
        assert isinstance(value, bool), '{} must be True or False.'.format(name)

    # data formatting

    if quantity == 'kappa': quantity = 'mode_kappa'
    if xquantity == 'kappa': xquantity = 'mode_kappa'
    tnames = tp.settings.to_tp()
    if invert: quantity, xquantity = xquantity, quantity
    quantity = tnames[quantity] if quantity in tnames else quantity
    xquantity = tnames[xquantity] if xquantity in tnames else xquantity

    data = tp.data.utilities.resolve(data, [quantity, xquantity],
                                   temperature=temperature,
                                   direction=direction)
    if verbose and 'temperature' in data['meta']:
        print('Using {} {}.'.format(data['meta']['temperature'],
                                    data['meta']['units']['temperature']))
    x = np.ravel(data[xquantity])
    y = np.abs(np.ravel(data[quantity]))
    mask = np.ma.masked_invalid(y).mask
    x = np.ma.masked_where(mask, x).compressed()
    y = np.ma.masked_where(mask, y).compressed()
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x_dens, y_dens, z_dens = x[idx], y[idx], z[idx]

    cmap = tp.plot.utilities.parse_colours(colour)

    # plotting

    ax.scatter(x_dens, y_dens, c=z_dens, cmap=cmap, **kwargs)

    # axes formatting

    if main:
        data[xquantity] = x_dens
        data[quantity] = y_dens
        format_waterfall(ax, data, quantity, xquantity)
        if invert:
            axlabels = tp.settings.inverted_labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.tick_params(axis='y', labelleft=False)
        else:
            axlabels = tp.settings.labels()
            ax.set_xlabel(axlabels[xquantity])
            ax.set_ylabel(axlabels[quantity])

    return

def format_waterfall(ax, data, yquantity, xquantity='frequency',
                     temperature=300, direction='avg', invert=False):
    """Formats axes for waterfall plots.

    Arguments
    ---------

        ax : axes
            axes to format.
        data : dict
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
        invert : bool, optional
            invert x- and y-axes. Default: False.

    Returns
    -------

        None
            formats ax directly.
    """


    if invert:
        xquantity, yquantity = yquantity, xquantity
        
    data2 = tp.data.utilities.resolve(data, [xquantity, yquantity],
                                    temperature=temperature,
                                    direction=direction)
    data2 = {'x': np.ravel(data2[xquantity]), 'y': np.ravel(data2[yquantity])}

    limit = {'x': ax.set_xlim,   'y': ax.set_ylim}
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
    ax.scatter(lim['x'], lim['y'], alpha=0)
    tp.plot.utilities.set_locators(ax, x=loc['x'], y=loc['y'])

    return
