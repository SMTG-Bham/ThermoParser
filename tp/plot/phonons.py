"""Tools for dealing with phonon dispersions.

Functions:
    add_dispersion:
        phonon dispersion.
    add_multi:
        phonon dispersions.
    add_alt_dispersion:
        phono3py quantity against high symmetry path.
    add_projected_dispersion:
        phonon dispersion with phono3py quantity on colour axis.
    add_alt_projected_dispersion:
        alt_dispersion + projection.
    add_wideband:
        phonon dispersion broadened according to scattering.

    get_equivalent_qpoint:
        converts phonopy to phono3py qpoints.

    formatting:
        formatting axes.
    tile_properties:
        tiling properties semi-intelligently.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tp
import warnings
warnings.filterwarnings('ignore', module='matplotlib')

def add_dispersion(ax, data, bandmin=None, bandmax=None, main=True, label=None,
                   colour='#800080', linestyle='solid', xmarkkwargs={},
                   **kwargs):
    """Adds a phonon band structure to a set of axes.

    Labels, colours and linestyles can be given one for the whole
    dispersion, or one for each band, with the last entry filling all
    remaining bands.

    Arguments:
        ax : axes
            axes to plot on.
        data : dict
            dispersion data containing:
                x : array-like
                    high-symmetry path.
                frequency : array-like
                    phonon frequencies.
                tick_position : array-like
                    x-tick positions.
                tick_label : array-like
                    x-tick labels.

        bandmin : int, optional
            Zero-indexed minimum band index to plot. Default: 0.
        bandmax : int, optional
            Zero-indexed maximum band index to plot. Default: max index.

        main : bool, optional
            set axis ticks, labels, limits. Default: True.
        label : str, optional
            legend label(s). Default: None

        colour : str or array-like, optional
            line colour(s). Note if retrieved from a colourmap or in
            [r,g,b] notation, the colour should be enclosed in [], e.g.
            colour = plt.get_cmap('winter_r')(linspace(0, 1, len(files)))
            tp.plot.phonons.add_dispersion(..., colour=[colour[0]], ...)
            Default: #800080.
        linestyle : str or array-like, optional
            linestyle(s) ('-', '--', '.-', ':'). Default: solid.

        xmarkkwargs : dict, optional
            keyword arguments for x markers passed to
            matplotlib.pyplot.axvline. Set color to None to turn off.
            Defaults: {'color':      'black',
                       'linewidth':  axis line width,
                       'rasterized': False}
        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.plot.
            Defaults: {'rasterized': False}
    """

    # defaults

    defkwargs = {'rasterized': False}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # data selection

    if bandmin is None:
        bandmin = 0
    else:
        bandmin = np.amax([0, bandmin])
    if bandmax is None:
        bandmax = len(data['frequency'][0])
    else:
        bandmax = np.amin([len(data['frequency'][0]), bandmax])
    f = np.array(data['frequency'])[:,bandmin:bandmax]

    # line appearance

    colour = tile_properties(colour, bandmin, bandmax)
    linestyle = tile_properties(linestyle, bandmin, bandmax)

    # prevents unintentionally repeated legend entries

    label2 = []
    if isinstance(label, str):
        label2.extend(['$\mathregular{{{}}}$'.format(label)])
    else:
        try:
            label = ['$\mathregular{{{}}}$'.format(l) for l in label]
            label2.extend(label)
        except Exception:
            label2.extend([label])
    label2.append(None)
    label2 = tile_properties(label2, bandmin, bandmax)

    # plotting

    for n in range(len(f[0])):
        ax.plot(data['x'], f[:,n], color=colour[n], linestyle=linestyle[n],
                label=label2[n], **kwargs)

    # axes formatting

    if main:
        if round(np.amin(f), 1) == 0:
            ax.set_ylim(bottom=0)
        formatting(ax, data, 'frequency', **xmarkkwargs)

    return

def add_multi(ax, data, bandmin=None, bandmax=None, main=True, label=None,
              colour='winter_r', linestyle='solid', xmarkkwargs={}, **kwargs):
    """Adds multiple phonon band structures to a set of axes.

    Scales the x-scales to match.

    Arguments:
        ax : axes
            axes to plot on.
        data : array-like
            dictionaries of dispersion data containing:
                x : array-like
                    high-symmetry path.
                frequency : array-like
                    phonon frequencies.
                tick_position : array-like
                    x-tick positions.
                tick_label : array-like
                    x-tick labels.

        bandmin : int, optional
            Zero-indexed minimum band index to plot. Default: 0.
        bandmax : int, optional
            Zero-indexed maximum band index to plot. Default: max index.

        main : bool, optional
            set axis ticks, labels, limits. Default: True.
        label : array-like, optional
            legend labels. Default: None

        colour : colourmap or str or array-like, optional
            colourmap or colourmap name or list of colours. Note [r,g,b]
            format colours should be enclosed in and additional [],
            i.e. [[[r,g,b]],...]. Default: winter_r.
        linestyle : str or array-like, optional
            linestyle(s) ('-', '--', '.-', ':'). Default: solid.

        xmarkkwargs : dict, optional
            keyword arguments for x markers passed to
            matplotlib.pyplot.axvline. Set color to None to turn off.
            Defaults: {'color':      'black',
                       'linewidth':  axis line width,
                       'rasterized': False}
        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.plot.
            Defaults: {'rasterized': False}
    """

    # defaults

    defkwargs = {'rasterized': False}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # line appearance

    try:
        colours = mpl.cm.get_cmap(colour)(np.linspace(0, 1, len(data)))
        colours = [[c] for c in colours]
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            colours = [[colour(i)] for i in np.linspace(0, 1, len(data))]
        elif isinstance(colour, str) and colour == 'skelton':
            colour = tp.plot.colour.skelton()
            colours = [[colour(i)] for i in np.linspace(0, 1, len(data))]
        else:
            colours = colour

    if label is None:
        label = np.full(len(data), None)

    if isinstance(linestyle, str) or len(linestyle) == 1:
        linestyle = np.repeat(linestyle, len(data))
    elif len(linestyle) < len(data):
        while len(linestyle) < len(data):
            linestyle.append(linestyle[-1])

    # plotting

    for i in range(len(data)):
        add_dispersion(ax, data[i], bandmin=bandmin, bandmax=bandmax,
                       main=False, label=label[i], colour=colours[i],
                       linestyle=linestyle[i], **kwargs)

    # axes formatting

    if main:
        if bandmin is None:
            bandmin = 0
        else:
            bandmin = np.amax([0, bandmin])
        if bandmax is None:
            bandmax = len(data[0]['frequency'][0])
        else:
            bandmax = np.amin([len(data[0]['frequency'][0]), bandmax])

        f = [d['frequency'] for d in data]
        f = np.array(f)[:,:,bandmin:bandmax]

        if round(np.amin(f), 1) == 0:
            ax.set_ylim(bottom=0)
        formatting(ax, data[0], 'frequency', **xmarkkwargs)

    return

def add_alt_dispersion(ax, data, pdata, quantity, bandmin=None, bandmax=None,
                       temperature=300, direction='avg', label=['Longitudinal',
                       'Transverse_1', 'Transverse_2', 'Optic'],
                       poscar='POSCAR', main=True, log=False,
                       interpolate=10000, smoothing=5, colour=['#44ffff',
                       '#ff8044', '#ff4444', '#00000010'], linestyle='-',
                       workers=32, xmarkkwargs={}, **kwargs):
    """Plots a phono3py quantity on a high-symmetry path.

    Labels, colours and linestyles can be given one for the whole
    dispersion, or one for each band, with the last entry filling all
    remaining bands. Requires a POSCAR.

    Arguments:
        ax : axes
            axes to plot on.
        data : dict
            Phono3py-like data containing:
                qpoint : array-like
                    q-point locations.
                (quantity) : array-like
                    plotted quantity.
                temperature : array-like, optional
                    temperature, if necessary for quantity.
        pdata : dict
            phonon dispersion data containing:
                q-point : array-like
                    qpoint locations.
                x : array-like
                    high-symmetry path.
                tick_position : array-like
                    x-tick positions.
                tick_label : array-like
                    x-tick labels.
        quantity : str
            quantity to plot.

        bandmin : int, optional
            Zero-indexed minimum band index to plot. Default: 0.
        bandmax : int, optional
            Zero-indexed maximum band index to plot. Default: max index.
        temperature : float, optional
            approximate temperature in K (finds closest). Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c,
            average/ avg or normal/ norm. Default: average.
        label : array-like, optional
            labels per line. A single dataset could have a single label,
            or the default labels the lines by type. You'll want to
            change this if a minimum band index is set.
            Default: ['Longitudinal', 'Transverse_1', 'Transverse_2',
            'Optic'].

        poscar : str, optional
            VASP POSCAR filepath. Default: POSCAR.

        main : bool, optional
            set axis ticks, label, limits. Default: True.
        log : bool, optional
            log the y scale. Default: False.
        interpolate : int, optional
            number of points per line. Default: 10,000.
        smoothing : int, optional
            every n points to sample. Default: 5.

        colour : str or array-like, optional
            line colour(s). Note if retrieved from a colourmap or in
            [r,g,b] notation, the colour should be enclosed in [], e.g.
            colour = plt.get_cmap('winter_r')(linspace(0, 1, len(files)))
            tp.plot.phonons.add_dispersion(..., colour=[colour[0]], ...)
            Default: ['#44ffff', '#ff8044', '#ff4444', '#00000010'].
        linestyle : str or array-like, optional
            linestyle(s) ('-', '--', '.-', ':'). Default: solid.

        workers : int, optional
            number of workers for paralellised section. Default: 32.

        xmarkkwargs : dict, optional
            keyword arguments for x markers passed to
            matplotlib.pyplot.axvline. Set color to None to turn off.
            Defaults: {'color':      'black',
                       'linewidth':  axis line width,
                       'rasterized': False}
        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.plot.
            Defaults: {'rasterized': False}
    """

    from functools import partial
    from multiprocessing import Pool
    from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar
    import pymatgen.symmetry.analyzer as pmg
    from scipy.interpolate import interp1d

    # defaults

    defkwargs = {'rasterized': False}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # Phono3py data formatting

    tnames = tp.settings.to_tp()
    if quantity in tnames: quantity = tnames[quantity]
    if quantity == 'kappa': quantity = 'mode_kappa'

    data = tp.data.resolve.resolve(data, quantity, temperature, direction)
    y = data[quantity]
    if bandmin is None:
        bandmin = 0
    else:
        bandmin = np.amax([0, bandmin])
    if bandmax is None:
        bandmax = len(y[0])
    else:
        bandmax = np.amin([len(y[0]), bandmax])
    y = np.array(y)[:,bandmin:bandmax]
    qk = data['qpoint']

    # Phonopy data formatting

    qp = pdata['qpoint']
    x = pdata['x']

    if len(x) > 100:
        xi = [x[i] for i in range(len(x)) if i % smoothing == 0]
        qpi = [qp[i] for i in range(len(qp)) if i % smoothing == 0]
    else:
        xi = x
        qpi = qp
    if xi[-1] != x[-1]: qpi.append(qp[-1]); xi.append(x[-1])

    # map Phono3py to Phonopy

    struct = Poscar.from_file(poscar).structure
    sg = SpacegroupAnalyzer(struct)
    symops = sg.get_point_group_operations(cartesian=False)
    geq = partial(get_equivalent_qpoint, np.array(qk), symops)
    with Pool(processes=workers) as pool:
        min_id = pool.map(geq, qpi)
    y2 = y[min_id, :]

    x, indices = np.unique(xi, return_index=True)
    y2 = np.array(y2)[indices]

    # interpolate

    x2 = np.linspace(min(x), max(x), interpolate)
    yinterp = interp1d(x, y2, kind='cubic', axis=0)
    y2 = np.abs(yinterp(x2))
    ysort = np.ravel(y2)
    ysort = ysort[ysort.argsort()]
    ymin = ysort[int(round(len(ysort)/100, 0))]
    ymax = ysort[-1]

    # line appearance

    colour = tile_properties(colour, bandmin, bandmax)
    linestyle = tile_properties(linestyle, bandmin, bandmax)

    # prevents unintentionally repeated legend entries

    label2 = []
    if isinstance(label, str):
        label2.extend(['$\mathregular{{{}}}$'.format(label)])
    else:
        try:
            label = ['$\mathregular{{{}}}$'.format(l) for l in label]
            label2.extend(label)
        except Exception:
            label2.extend([label])
    label2.append(None)
    label2 = tile_properties(label2, bandmin, bandmax)

    # plotting

    for n in range(len(y2[0])):
        ax.plot(x2, y2[:,n], color=colour[n], linestyle=linestyle[n],
                label=label2[n], **kwargs)

    # axes formatting

    if main:
        ax.set_ylim(ymin, ymax)
        formatting(ax, pdata, quantity, log=log, **xmarkkwargs)

    return

def add_projected_dispersion(ax, data, pdata, quantity, bandmin=None,
                             bandmax=None, temperature=300, direction='avg',
                             poscar='POSCAR', main=True, interpolate=500,
                             colour='viridis_r', cmin=None, cmax=None,
                             cscale=None, unoccupied='grey', workers=32,
                             xmarkkwargs={}, **kwargs):
    """Plots a phonon dispersion with projected colour.

    Plots a phonon dispersion, and projects a quantity onto the colour
    axis. Requires a POSCAR.

    Arguments:
        ax : axes
            axes to plot on.
        data : dict
            Phono3py-like data containing:
                qpoint : array-like
                    q-point locations.
                (quantity) : array-like
                    projected quantity.
                temperature : array-like, optional
                    temperature, if necessary for quantity.
        pdata : dict
            phonon dispersion data containing:
                q-point : array-like
                    qpoint locations.
                x : array-like
                    high-symmetry path.
                frequency : array-like
                    phonon frequencies.
                tick_position : array-like
                    x-tick positions.
                tick_label : array-like
                    x-tick labels.
        quantity : str
            quantity to project.

        bandmin : int, optional
            Zero-indexed minimum band index to plot. Default: 0.
        bandmax : int, optional
            Zero-indexed maximum band index to plot. Default: max index.
        temperature : float, optional
            approximate temperature in K (finds closest). Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c,
            average/ avg or normal/ norm. Default: average.

        poscar : str, optional
            VASP POSCAR filepath. Default: POSCAR.

        main : bool, optional
            set axis ticks, label, limits. Default: True.
        interpolate : int, optional
            number of points per path (e.g. gamma -> X) per line.
            Default: 500.

        colour : colormap or str, optional
            colourmap or colourmap name. Default: viridis_r.
        cmin : float, optional
            colour scale minimum. Default: display 99 % data.
        cmax : float, optional
            colour scale maximum. Default: display 99.9 % data.
        cscale : str, optional
            override colour scale (linear/ log). Default: None.
        unoccupied : str, optional
            if the colour variable is occuption, values below 1 are
            coloured in this colour. Id set to None, or cmin is set,
            this feature is turned off. Default: grey.

        workers : int, optional
            number of workers for paralellised section. Default: 32.

        xmarkkwargs : dict, optional
            keyword arguments for x markers passed to
            matplotlib.pyplot.axvline. Set color to None to turn off.
            Defaults: {'color':      'black',
                       'linewidth':  axis line width,
                       'rasterized': False}
        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.scatter.
            Defaults: {'marker':     '.',
                       'rasterized': True,
                       's':          1}

    Returns:
        colorbar
            colour bar for projected data.
    """

    from functools import partial
    from multiprocessing import Pool
    from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar
    import pymatgen.symmetry.analyzer as pmg
    from copy import copy
    from scipy.interpolate import interp1d

    # defaults

    defkwargs = {'marker':     '.',
                 'rasterized': True,
                 's':          1}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # Phono3py data formatting

    tnames = tp.settings.to_tp()
    quantity = tnames[quantity] if quantity in tnames else quantity
    if quantity == 'kappa': quantity = 'mode_kappa'

    data = tp.data.resolve.resolve(data, quantity, temperature, direction)
    c = data[quantity]
    if bandmin is None:
        bandmin = 0
    else:
        bandmin = np.amax([0, bandmin])
    if bandmax is None:
        bandmax = len(c[0])
    else:
        bandmax = np.amin([len(c[0]), bandmax])
    c = np.array(c)[:, bandmin:bandmax]
    qk = data['qpoint']

    # Phonopy data formatting

    qp = pdata['qpoint']
    x = pdata['x']
    f = np.array(pdata['frequency'])[:, bandmin:bandmax]

    # map Phono3py to Phonopy

    struct = Poscar.from_file(poscar).structure
    sg = SpacegroupAnalyzer(struct)
    symops = sg.get_point_group_operations(cartesian=False)
    geq = partial(get_equivalent_qpoint, np.array(qk), symops)
    with Pool(processes=workers) as pool:
        min_id = pool.map(geq, qp)
    c = c[min_id, :]

    x, indices = np.unique(x, return_index=True)
    f = np.array(f)[indices]
    c = np.array(c)[indices]

    # interpolate

    index = [0, 0]
    x2, f2, c2 = [], [], []
    for d in pdata['tick_position'][1:]:
        index[0] = index[1] + 1
        index[1] = next(i[0] for i in enumerate(x) if i[1] == d)
        xtemp = np.linspace(x[index[0]], x[index[1]], interpolate)
        finterp = interp1d(x[index[0]:index[1]], f[index[0]:index[1]],
                           kind='cubic', axis=0, fill_value='extrapolate')
        x2.append(xtemp)
        f2.append(finterp(xtemp))

        cinterp = interp1d(x[index[0]:index[1]], c[index[0]:index[1]],
                           kind='cubic', axis=0, fill_value='extrapolate')
        c2.append(np.abs(cinterp(xtemp)))

    if isinstance(colour, str):
        cmap = copy(plt.cm.get_cmap(colour))
    else:
        cmap = copy(colour)
    cnorm, extend = tp.plot.utilities.colour_scale(c2, quantity, cmap, cmin,
                                                   cmax, cscale, unoccupied)

    # plotting

    for n in range(len(f2[0][0])):
        for i in range(len(x2)):
            line = ax.scatter(x2[i], np.array(f2[i])[:,n], c=np.array(c2[i])[:,n], cmap=cmap, norm=cnorm,
                              **kwargs)

    # axes formatting

    axlabels = tp.settings.labels()
    cbar = plt.colorbar(line, extend=extend)
    cbar.set_alpha(1)
    cbar.set_label(axlabels[quantity])
    tp.plot.utilities.set_locators(cbar.ax, y=cbar.ax.yaxis.get_scale())
    cbar.draw_all()

    if main:
        if round(np.amin(f), 1) == 0:
            ax.set_ylim(bottom=0)
        formatting(ax, pdata, 'frequency', **xmarkkwargs)

    return cbar

def add_alt_projected_dispersion(ax, data, pdata, quantity, projected,
                                 bandmin=None, bandmax=None, temperature=300,
                                 direction='avg', poscar='POSCAR', main=True,
                                 log=False, interpolate=10000, smoothing=10,
                                 colour='viridis', cmin=None, cmax=None,
                                 cscale=None, unoccupied='grey', workers=32,
                                 xmarkkwargs={}, **kwargs):
    """Plots a phono3py quantity on a high-symmetry path and projection.

    Just because you can, doesn't mean you should. A maxim I may fail to
    live up to, so I leave it to you, dear reader, to decide for
    yourself. Requires a POSCAR.

    Arguments:
        ax : axes
            axes to plot on.
        data : dict
            Phono3py-like data containing:
                qpoint : array-like
                    q-point locations.
                (quantities) : array-like
                    plotted and projected quantities.
                temperature : array-like, optional
                    temperature, if necessary for quantities.
        pdata : dict
            phonon dispersion data containing:
                q-point : array-like
                    qpoint locations.
                x : array-like
                    high-symmetry path.
                tick_position : array-like
                    x-tick positions.
                tick_label : array-like
                    x-tick labels.
        quantity : str
            quantity to plot.
        projected: str
            quantity to project.

        bandmin : int, optional
            Zero-indexed minimum band index to plot. Default: 0.
        bandmax : int, optional
            Zero-indexed maximum band index to plot. Default: max index.
        temperature : float, optional
            approximate temperature in K (finds closest). Default: 300.
        direction : str, optional
            direction from anisotropic data, accepts x-z/ a-c,
            average/ avg or normal/ norm. Default: average.

        poscar : str, optional
            VASP POSCAR filepath. Default: POSCAR.

        main : bool, optional
            set axis ticks, label, limits. Default: True.
        log : bool, optional
            log the y scale. Default: False.
        interpolate : int, optional
            number of points per line. Default: 10,000.
        smoothing : int, optional
            every n points to sample. Default: 10.

        colour : colormap or str, optional
            colourmap or colourmap name. Default: viridis.
        cmin : float, optional
            colour scale minimum. Default: display 99 % data.
        cmax : float, optional
            colour scale maximum. Default: display 99.9 % data.
        cscale : str, optional
            override colour scale (linear/ log). Default: None.
        unoccupied : str, optional
            if the colour variable is occuption, values below 1 are
            coloured in this colour. Id set to None, or cmin is set,
            this feature is turned off. Default: grey.

        workers : int, optional
            number of workers for paralellised section. Default: 32.

        xmarkkwargs : dict, optional
            keyword arguments for x markers passed to
            matplotlib.pyplot.axvline. Set color to None to turn off.
            Defaults: {'color':      'black',
                       'linewidth':  axis line width,
                       'rasterized': False}
        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.scatter.
            Defaults: {'marker':     '.',
                       'rasterized': True,
                       's':          1}

    Returns:
        colorbar
            colour bar for projected data.
    """

    from functools import partial
    from multiprocessing import Pool
    from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar
    import pymatgen.symmetry.analyzer as pmg
    from copy import copy
    from scipy.interpolate import interp1d

    # defaults

    defkwargs = {'marker':     '.',
                 'rasterized': True,
                 's':          1}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # Phono3py data formatting

    tnames = tp.settings.to_tp()
    if quantity in tnames: quantity = tnames[quantity]
    if projected in tnames: projected = tnames[projected]
    if quantity == 'kappa': quantity = 'mode_kappa'
    if projected == 'kappa': projected = 'mode_kappa'
    qs = quantity if quantity == projected else [quantity, projected]

    data = tp.data.resolve.resolve(data, qs, temperature, direction)
    y = data[quantity]
    c = data[projected]
    if bandmin is None:
        bandmin = 0
    else:
        bandmin = np.amax([0, bandmin])
    if bandmax is None:
        bandmax = len(y[0])
    else:
        bandmax = np.amin([len(y[0]), bandmax])
    y = np.array(y)[:, bandmin:bandmax]
    c = np.array(c)[:, bandmin:bandmax]
    qk = data['qpoint']

    # Phonopy data formatting

    qp = pdata['qpoint']
    x = pdata['x']

    # map Phono3py to Phonopy

    struct = Poscar.from_file(poscar).structure
    sg = SpacegroupAnalyzer(struct)
    symops = sg.get_point_group_operations(cartesian=False)
    geq = partial(get_equivalent_qpoint, np.array(qk), symops)
    with Pool(processes=workers) as pool:
        min_id = pool.map(geq, qp)
    y2 = y[min_id, :]
    c2 = c[min_id, :]

    x, indices = np.unique(x, return_index=True)
    y2 = np.array(y2)[indices]
    c2 = np.array(c2)[indices]
    if len(x) > 100:
        xi = [x[i] for i in range(len(x)) if i % smoothing == 0]
        y2 = [y2[i] for i in range(len(y2)) if i % smoothing == 0]
    else:
        xi = x
        y2 = y2
    if xi[-1] != x[-1]: y2.append(y2[-1]); xi.append(x[-1])

    # interpolate

    x2 = np.linspace(min(x), max(x), interpolate)
    yinterp = interp1d(xi, y2, kind='cubic', axis=0)
    y2 = np.abs(yinterp(x2))
    ysort = np.ravel(y2)
    ysort = ysort[ysort.argsort()]
    ymin = ysort[int(round(len(ysort)/100 - 1, 0))]
    ymax = ysort[int(round(len(ysort)*99.9/100 - 1, 0))]

    cinterp = interp1d(x, c2, kind='cubic', axis=0)
    c2 = np.abs(cinterp(x2))
    if isinstance(colour, str):
        cmap = copy(plt.cm.get_cmap(colour))
    else:
        cmap = copy(colour)
    cnorm, extend = tp.plot.utilities.colour_scale(c2, projected, cmap, cmin,
                                                   cmax, cscale, unoccupied)

    # plotting

    for n in range(len(y2[0])):
        line = ax.scatter(x2, y2[:,n], c=c2[:,n], cmap=cmap, norm=cnorm,
                          **kwargs)

    # axes formatting

    axlabels = tp.settings.labels()
    cbar = plt.colorbar(line, extend=extend)
    cbar.set_alpha(1)
    cbar.set_label(axlabels[projected])
    tp.plot.utilities.set_locators(cbar.ax, y=cbar.ax.yaxis.get_scale())
    cbar.draw_all()

    if main:
        ax.set_ylim(ymin, ymax)
        formatting(ax, pdata, quantity, log=log, **xmarkkwargs)

    return cbar

def add_wideband(ax, kdata, pdata, temperature=300, poscar='POSCAR', main=True,
                 smoothing=5, colour='viridis', workers=32, xmarkkwargs={},
                 **kwargs):
    """Plots a phonon dispersion with broadened bands.

    Requires a POSCAR.

    Arguments:

        ax : axes
            axes to plot on.
        kdata : dict
            Phono3py-like data containing:
                qpoint : array-like
                    q-point locations.
                gamma : array-like
                    imaginary component of the self-energy.
                temperature : array-like
                    temperature.
        pdata : dict
            phonon dispersion data containing:
                q-point : array-like
                    qpoint locations.
                x : array-like
                    high-symmetry path.
                frequency : array-like
                    phonon frequencies.
                tick_position : array-like
                    x-tick positions.
                tick_label : array-like
                    x-tick labels.

        temperature : float, optional
            approximate temperature in K (finds closest). Default: 300.

        poscar : str, optional
            VASP POSCAR filepath. Default: POSCAR.

        main : bool, optional
            set axis ticks, label, limits. Default: True.
        smoothing : int, optional
            every n points to sample. Default: 5.

        colour : colormap or str or list, optional
            colourmap or colourmap name or max #RRGGBB colour (fades to
            white) or min and max #RRGGBB colours. Default: viridis.

        workers : int, optional
            number of workers for paralellised section. Default: 32.

        xmarkkwargs : dict, optional
            keyword arguments for x markers passed to
            matplotlib.pyplot.axvline. Set color to None to turn off.
            Defaults: {'color':      None,
                       'linewidth':  axis line width,
                       'rasterized': False}
        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.pcolormesh.
            Defaults: {'rasterized': True,
                       'shading':    'auto'}
    """

    from functools import partial
    from multiprocessing import Pool
    from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar
    import pymatgen.symmetry.analyzer as pmg
    from scipy.interpolate import interp1d
    from tp.calculate import lorentzian

    # defaults

    defkwargs = {'rasterized': True,
                 'shading':    'auto'}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    defxmarkkwargs = {'color': None}
    for key in defxmarkkwargs:
        if key not in xmarkkwargs:
            xmarkkwargs[key] = defxmarkkwargs[key]

    # Phono3py data formatting

    kdata = tp.data.resolve.resolve(kdata, 'gamma', temperature)
    c = np.array(kdata['gamma'])
    qk = kdata['qpoint']

    # Phonopy data formatting

    qp = pdata['qpoint']
    x = pdata['x']

    qpi = [qp[i] for i in range(len(qp)) if i % smoothing == 0]
    xi = [x[i] for i in range(len(x)) if i % smoothing == 0]
    if xi[-1] != x[-1]: qpi.append(qp[-1]); xi.append(x[-1])

    # map Phono3py to Phonopy

    struct = Poscar.from_file(poscar).structure
    sg = SpacegroupAnalyzer(struct)
    symops = sg.get_point_group_operations(cartesian=False)
    geq = partial(get_equivalent_qpoint, np.array(qk), symops)
    with Pool(processes=workers) as pool:
        min_id = pool.map(geq, qpi)
    c2 = c[min_id, :]

    x, indices = np.unique(x, return_index=True)
    f = np.array(pdata['frequency'])[indices]

    # interpolate

    x2 = np.linspace(min(x), max(x), 2500)
    finterp = interp1d(x, f, kind='cubic', axis=0)
    f = finterp(x2)

    cinterp = interp1d(xi, c2, kind='cubic', axis=0)
    c2 = np.abs(cinterp(x2))
    fmax = np.amax(np.add(f, c2))
    fmin = np.amin(np.subtract(f, c2))
    f2 = np.linspace(fmin, fmax, 2500)

    # broadening

    area = np.zeros((len(x2), len(f2)))
    for q in range(len(area)):
        for b in range(len(c2[q])):
            area[q] = np.add(area[q], lorentzian(f2, f[q][b], c2[q][b]))

    cnorm = mpl.colors.LogNorm(vmin=np.amin(area), vmax=np.amax(area))

    # colours

    try:
        colours = mpl.cm.get_cmap(colour)
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            colours = colour
        else:
            try:
                colours = tp.plot.colour.linear(colour)
            except Exception:
                colours = tp.plot.colour.linear(colour[1], colour[0])

    # plotting

    ax.pcolormesh(x2, f2, np.transpose(area), cmap=colours, norm=cnorm,
                  **kwargs)

    # axes formatting

    if main:
        if round(np.amin(f), 1) == 0:
            ax.set_ylim(bottom=0)
        formatting(ax, pdata, 'frequency', **xmarkkwargs)

    return

def get_equivalent_qpoint(qk, symops, qp, tol=1e-2):
    """Finds the closest phono3py qpoint to a phonopy qpoint.

    Arguments:
        qk : array-like
            all qpoints from the phono3py kappa file.
        symmops
            symmetry operations (e.g. from Pymatgen)
        qp : array-like
            single qpoint from the phonon dispersion.

        tol : float, optional
            tolerance. Default: 1e-2.

    Returns:
        int
            nearest qpoint index.
    """

    from pymatgen.util.coord import pbc_diff

    # Equivalent qpoints

    points = np.dot(qp, [m.rotation_matrix for m in symops])

    # Remove duplicates

    rm_list = []
    for i in range(len(points) - 1):
        for j in range(i + 1, len(points)):
            if np.allclose(pbc_diff(points[i], points[j]), [0, 0, 0], tol):
                rm_list.append(i)
                break
    seq = np.delete(points, rm_list, axis=0)

    # Find nearest

    dists = [np.min(np.sum(np.power(k - seq, 2), axis=1)) for k in qk]

    return dists.index(np.min(dists))

def formatting(ax, data, yquantity='frequency', log=False, **kwargs):
    """Formats the axes of phonon plots.

    Arguments:
        ax : axes
            axes to format.
        pdata : dict
            phonon dispersion data containing:
                x : array-like
                    high-symmetry path.
                tick_position : array-like
                    x-tick positions.
                tick_label : array-like
                    x-tick labels.
        yquantity : str, optional
            y variable. Default: frequency.
        log : bool, optional
            log the y scale. Default: False.

        **kwargs : dict, optional
            keyword arguments passed to matplotlib.pyplot.axvline.
            Set color to None to turn off.
            Defaults: {'color':      'black',
                       'linewidth':  axis line width,
                       'rasterized': False}
    """

    # defaults

    defkwargs = {'color':      'black',
                 'linewidth':  ax.spines['bottom'].get_linewidth(),
                 'rasterized': False}
    for key in defkwargs:
        if key not in kwargs:
            kwargs[key] = defkwargs[key]

    # lines

    if kwargs['color'] is not None:
        for d in data['tick_position']:
            ax.axvline(d, **kwargs)

    # formatting

    axlabels = tp.settings.labels()
    ax.set_xlabel(axlabels['wavevector'])
    ax.set_ylabel(axlabels[yquantity])

    ax.set_xticks(data['tick_position'])
    ax.set_xticklabels(data['tick_label'])
    ax.tick_params(axis='x', which='minor', top=False, bottom=False)
    ax.set_xlim(data['x'][0], data['x'][-1])
    scale = 'log' if log else 'linear'
    tp.plot.utilities.set_locators(ax, y=scale)

    return

def tile_properties(properties, bandmin, bandmax):
    """Tiles properties for dispersion and alt_dispersion

    Allows for different colour formats and fills out arrays with the
    last element or selects a subset.

    Arguments:
        property : array-like or str
            array or string to tile.
        bandmin : int
            minimum band.
        bandmax : int
            maximum band.

    Returns:
        tiled
            array.
    """

    bdiff = bandmax - bandmin
    if isinstance(properties, str):
        tiled = np.repeat(properties, bdiff)
    elif len(properties) == 1:
        tiled = np.tile(properties, (bdiff, 1))
    elif len(properties) == bdiff:
        tiled = properties
    elif len(properties) > bdiff:
        if len(properties) >= bandmax:
            tiled = properties[bandmin:bandmax]
        else:
            tiled = properties[:bdiff]
    elif len(properties) < bdiff:
        tiled = list(properties)
        while len(tiled) < bdiff:
            tiled.append(tiled[-1])

    return tiled
