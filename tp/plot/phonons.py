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
        alt_dispersion + projection
    add_wideband:
        phonon dispersion broadened according to scattering.

    get_equivalent_qpoint:
        converts phonopy to phono3py qpoints.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tp.settings
import warnings
warnings.filterwarnings('ignore', module='matplotlib')

def add_dispersion(ax, data, bandrange=None, main=True, label=None,
                   colour='#800080', linestyle='solid', axcolour='black'):
    """Adds a phonon band structure to a set of axes.

    Future: work out a good way to add line style and colour while
    supporting cyclers.

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

        bandrange : array-like, optional
            minumum and maximum band indices to plot. Default: None.

        main : bool, optional
            set axis ticks, labels, limits. Default: True.
        label : str, optional
            legend label. Default: None

        colour : str or array-like, optional
            line colour(s). Note if retrieved from a colourmap or in
            [r,g,b] notation, the colour should be enclosed in [], e.g.
            colour = plt.get_cmap('winter_r')(linspace(0, 1, len(files)))
            tp.plot.phonons.add_dispersion(..., colour=[colour[0]], ...)
            Default: #800080.
        linestyle : str or array-like, optional
            linestyle(s) ('-', '--', '.-', ':'). Default: solid.
        axcolour : str
            axis colour. Default: black.

    Returns:
        axes
            axes with dispersion.
    """

    if bandrange is None:
        bmin = 0; bmax = len(data['frequency'][0])
    else:
        bmin = np.amax([0, bandrange[0]])
        bmax = np.amin([len(data['frequency'][0]), bandrange[1]])
    bdiff = bmax - bmin
    f = np.array(data['frequency'])[:,bmin:bmax]

    # allows single colours/ linestyles or arrays, which can be filled out
    if isinstance(colour, str):
        colours = np.repeat(colour, bdiff)
    elif len(colour) == 1:
        if isinstance(colour, str):
            colours = np.repeat(colour, bdiff)
        else:
            colours = np.tile(colour, (bdiff, 1))
    elif len(colour) == bdiff:
        colours = colour
    elif len(colour) > bdiff:
        colours = colour[bmin:bmax]
    elif len(colour) < bdiff:
        colours = list(colour)
        while len(colours) < bdiff:
            colours.append(colour[-1])

    if isinstance(linestyle, str) or len(linestyle) == 1:
        linestyles = np.repeat(linestyle, bdiff)
    elif len(linestyle) == bdiff:
        linestyles = linestyle
    elif len(linestyle) > bdiff:
        linestyles = linestyle[bmin:bmax]
    elif len(linestyle) < bdiff:
        linestyles = linestyle
        while len(linestyles) < bdiff:
            linestyles.append(linestyle[-1])

    # ensures only one legend entry
    if bdiff == 1:
        ax.plot(data['x'], f, label=label, color=colours[0],
                              linestyle=linestyles[0])
    else:
        ax.plot(data['x'], f[:,0], label=label, color=colours[0],
                                   linestyle=linestyles[0])
        for n in range(bdiff):
            ax.plot(data['x'], f[:,n], color=colours[n],
                                       linestyle=linestyles[n])

    if main:
        axlabels = tp.settings.labels()
        spinewidth = ax.spines['bottom'].get_linewidth()
        for d in data['tick_position']:
            ax.axvline(d, color=axcolour, linewidth=spinewidth)

        ax.set_xlabel(axlabels['wavevector'])
        ax.set_ylabel(axlabels['frequency'])

        ax.set_xticks(data['tick_position'])
        ax.set_xticklabels(data['tick_label'])
        ax.tick_params(axis='x', which='minor', top=False, bottom=False)
        ax.set_xlim(data['x'][0], data['x'][-1])

        ax.yaxis.set_major_locator(tp.settings.locator()['major'])
        ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])
        if round(np.amin(f), 1) == 0:
            ax.set_ylim(bottom=0)

    return ax

def add_multi(ax, data, bandrange=None, main=True, label=None,
              colour='winter_r', linestyle='solid', axcolour='black'):
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

        bandrange : array-like, optional
            minumum and maximum band indices to plot. Default: None.

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
        axcolour : str
            axis colour. Default: black.

    Returns:
        axes
            axes with dispersions.
    """

    try:
        colours = mpl.cm.get_cmap(colour)(np.linspace(0, 1, len(data)))
        colours = [[c] for c in colours]
    except Exception:
        if isinstance(colour, mpl.colors.ListedColormap):
            colours = [[colour(i)] for i in np.linspace(0, 1, len(data))]
        elif isinstance(colour, str) and colour == 'skelton':
            from tp import plot
            colour = plot.colour.skelton()
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

    for i in range(len(data)):
        ax = add_dispersion(ax, data[i], bandrange=bandrange, main=False,
                            label=label[i], colour=colours[i],
                            linestyle=linestyle[i], axcolour=axcolour)

    if main:
        axlabels = tp.settings.labels()
        spinewidth = ax.spines['bottom'].get_linewidth()
        for d in data[0]['tick_position']:
            ax.axvline(d, color=axcolour, linewidth=spinewidth)

        ax.set_xlabel(axlabels['wavevector'])
        ax.set_ylabel(axlabels['frequency'])

        ax.set_xticks(data[0]['tick_position'])
        ax.set_xticklabels(data[0]['tick_label'])
        ax.tick_params(axis='x', which='minor', top=False, bottom=False)
        ax.set_xlim(data[0]['x'][0], data[0]['x'][-1])

        ax.yaxis.set_major_locator(tp.settings.locator()['major'])
        ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])
        f = [d['frequency'] for d in data]
        if round(np.amin(f), 1) == 0:
            ax.set_ylim(bottom=0)

    return ax

def add_alt_dispersion(ax, data, pdata, quantity, bandrange=None,
                       temperature=300, direction='avg',
                       label=['Longitudinal', 'Transverse_1',
                       'Transverse_2', 'Optic'], poscar='POSCAR',
                       main=True, log=False, interpolate=5,
                       colour=['#44ffff', '#ff8044', '#ff4444',
                       '#00000010'], linestyle='-', axcolour='black',
                       rasterise=True, workers=32):
    """Plots a phono3py quantity on a high-symmetry path.

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

        bandrange : array-like, optional
            minumum and maximum band indices to plot. Default: None.
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
            every n points to sample. Default: 5.

        colour : str or array-like, optional
            line colour(s). If too few colours are selected, the last
            one is repeated for the rest of the bands. Note if retrieved
            from a colourmap or in [r,g,b] notation, the colour should
            be enclosed in [], e.g.
            colour = plt.get_cmap('winter_r')(linspace(0, 1, len(files)))
            tp.plot.phonons.add_dispersion(..., colour=[colour[0]], ...)
            Default: ['#44ffff', '#ff8044', '#ff4444', '#00000010'].
        linestyle : str or array-like, optional
            linestyle(s) ('-', '--', '.-', ':'). Default: solid.
        axcolour : str
            axis colour. Default: black.
        rasterise : bool, optional
            rasterise plot. Default: True.

        workers : int, optional
            number of workers for paralellised section. Default: 32.

    Returns:
        axes
            axes with "dispersion".
    """

    from functools import partial
    from multiprocessing import Pool
    from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar
    import pymatgen.symmetry.analyzer as pmg
    from scipy.interpolate import interp1d

    # Phono3py
    tnames = tp.settings.to_tp()
    if quantity in tnames: quantity = tnames[quantity]
    data = tp.data.resolve.resolve(data, quantity, temperature, direction)
    y = data[quantity]
    if bandrange is None:
        bmin = 0; bmax = len(y[0])
    else:
        bmin = np.amax([0, bandrange[0]])
        bmax = np.amin([len(y[0]), bandrange[1]])
    bdiff = bmax - bmin
    y = np.array(y)[:,bmin:bmax]
    qk = data['qpoint']

    # Phonopy
    qp = pdata['qpoint']
    x = pdata['x']

    qpi = [qp[i] for i in range(len(qp)) if i % interpolate == 0]
    xi = [x[i] for i in range(len(x)) if i % interpolate == 0]
    if xi[-1] != x[-1]: qpi.append(qp[-1]); xi.append(x[-1])

    # map Phono3py to Phonopy
    struct = Poscar.from_file(poscar).structure
    sg = SpacegroupAnalyzer(struct)
    symops = sg.get_point_group_operations(cartesian=False)
    geq = partial(get_equivalent_qpoint, np.array(qk), symops)
    with Pool(processes=workers) as pool:
        min_id = pool.map(geq, qpi)
    y2 = y[min_id, :]

    # Interpolate
    x2 = np.linspace(min(xi), max(xi), 2500)
    yinterp = interp1d(xi, y2, kind='cubic', axis=0)
    y2 = np.abs(yinterp(x2))
    ysort = np.ravel(y2)
    ysort = ysort[ysort.argsort()]
    ymin = ysort[int(round(len(ysort)/100, 0))]
    ymax = ysort[-1]

    # allows single colours/ linestyles or arrays, which can be filled out
    if isinstance(colour, str) or len(colour) == 1:
        colours = np.repeat(colour, bdiff)
    elif len(colour) == bdiff:
        colours = colour
    elif len(colour) > bdiff:
        colours = colour[bmin:bmax]
    elif len(colour) < bdiff:
        colours = colour
        while len(colours) < bdiff:
            colours.append(colour[-1])

    if isinstance(linestyle, str) or len(linestyle) == 1:
        linestyles = np.repeat(linestyle, bdiff)
    elif len(linestyle) == bdiff:
        linestyles = linestyle
    elif len(linestyle) > bdiff:
        linestyles = linestyle[bmin:bmax]
    elif len(linestyle) < bdiff:
        linestyles = linestyle
        while len(linestyles) < bdiff:
            linestyles.append(linestyle[-1])

    # ensures only one legend entry
    if bdiff == 1:
        ax.plot(x2, y2, label=label, color=colours[0], linestyle=linestyles[0])
    else:
        for i in range(len(label)):
            ax.plot(x2, y2[:,i], label='$\mathregular{{{}}}$'.format(label[i]),
                    color=colours[i], linestyle=linestyles[i])
        if len(label) < bdiff:
            for n in range(len(label), bdiff):
                ax.plot(x2, y2[:,n], color=colours[n], linestyle=linestyles[n])

    if main:
        axlabels = tp.settings.labels()
        spinewidth = ax.spines['bottom'].get_linewidth()
        for d in pdata['tick_position']:
            ax.axvline(d, color=axcolour, linewidth=spinewidth)

        ax.set_xlabel(axlabels['wavevector'])
        ax.set_ylabel(axlabels[quantity])

        ax.set_xticks(pdata['tick_position'])
        ax.set_xticklabels(pdata['tick_label'])
        ax.tick_params(axis='x', which='minor', top=False, bottom=False)
        ax.set_xlim(x2[0], x2[-1])

        ax.set_ylim(ymin, ymax)
        if log:
            ax.set_yscale('log')
            ax.yaxis.set_major_locator(tp.settings.locator()['log'])
        else:
            ax.yaxis.set_major_locator(tp.settings.locator()['major'])
            ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])

    return ax

def add_projected_dispersion(ax, data, pdata, quantity, bandrange=None,
                             temperature=300, direction='avg', poscar='POSCAR',
                             main=True, interpolate=5, colour='viridis_r',
                             axcolour='black', rasterise=True, workers=32):
    """Plots a phonon dispersion with projected colour.

    Plots a phonon dispersion, and projects a quantity onto the colour
    axis.

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
            quantity to plot.

        bandrange : array-like, optional
            minumum and maximum band indices to plot. Default: None.
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
            every n points to sample. Default: 5.

        colour : colormap or str, optional
            colourmap or colourmap name. Default: viridis_r.
        axcolour : str
            axis colour. Default: black.
        rasterise : bool, optional
            rasterise plot. Default: True.

        workers : int, optional
            number of workers for paralellised section. Default: 32.

    Returns:
        axes
            axes with projected dispersion.
        colorbar
            colour bar for projected data.
    """

    from functools import partial
    from multiprocessing import Pool
    from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar
    import pymatgen.symmetry.analyzer as pmg
    from scipy.interpolate import interp1d

    # Phono3py
    tnames = tp.settings.to_tp()
    quantity = tnames[quantity] if quantity in tnames else quantity
    data = tp.data.resolve.resolve(data, quantity, temperature, direction)
    c = data[quantity]
    if bandrange is None:
        bmin = 0; bmax = len(c[0])
    else:
        bmin = np.amax([0, bandrange[0]])
        bmax = np.amin([len(c[0]), bandrange[1]])
    bdiff = bmax - bmin
    c = np.array(c)[:, bmin:bmax]
    qk = data['qpoint']

    # Phonopy
    qp = pdata['qpoint']
    x = pdata['x']
    f = np.array(pdata['frequency'])[:, bmin:bmax]

    qpi = [qp[i] for i in range(len(qp)) if i % interpolate == 0]
    xi = [x[i] for i in range(len(x)) if i % interpolate == 0]
    if xi[-1] != x[-1]: qpi.append(qp[-1]); xi.append(x[-1])

    # Map Phono3py to Phonopy
    struct = Poscar.from_file(poscar).structure
    sg = SpacegroupAnalyzer(struct)
    symops = sg.get_point_group_operations(cartesian=False)
    geq = partial(get_equivalent_qpoint, np.array(qk), symops)
    with Pool(processes=workers) as pool:
        min_id = pool.map(geq, qpi)
    c2 = c[min_id, :]

    x, indices = np.unique(x, return_index=True)
    f = np.array(f)[indices]

    # Interpolate
    x2 = np.linspace(min(x), max(x), 2500)
    finterp = interp1d(x, f, kind='cubic', axis=0)
    f = finterp(x2)

    cinterp = interp1d(xi, c2, kind='cubic', axis=0)
    c2 = np.abs(cinterp(x2))
    csort = np.ravel(c2)
    csort = csort[csort.argsort()]
    cmin = csort[int(round(len(csort)/100, 0))]
    cmax = csort[-1]
    cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    for n in range(bdiff):
        line = ax.scatter(x2, f[:,n], c=c2[:,n], cmap=colour, norm=cnorm,
                          marker='.', s=1, rasterized=rasterise)

    if main:
        axlabels = tp.settings.labels()
        spinewidth = ax.spines['bottom'].get_linewidth()
        for d in pdata['tick_position']:
            ax.axvline(d, color=axcolour, linewidth=spinewidth)

        cbar = plt.colorbar(line)
        cbar.set_label(axlabels[quantity])

        ax.set_xlabel(axlabels['wavevector'])
        ax.set_ylabel(axlabels['frequency'])

        ax.set_xticks(pdata['tick_position'])
        ax.set_xticklabels(pdata['tick_label'])
        ax.tick_params(axis='x', which='minor', top=False, bottom=False)
        ax.set_xlim(x2[0], x2[-1])

        ax.yaxis.set_major_locator(tp.settings.locator()['major'])
        ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])
        if round(np.amin(f), 1) == 0:
            ax.set_ylim(bottom=0)

    return ax, cbar

def add_alt_projected_dispersion(ax, data, pdata, quantity, projected,
                                 temperature=300, direction='avg',
                                 poscar='POSCAR', main=True, log=False,
                                 interpolate=5, colour='viridis',
                                 axcolour='black', rasterise=True, workers=32):
    """Plots a phono3py quantity on a high-symmetry path and projection.

    Just because you can, doesn't mean you should. A maxim I may fail to
    live up to, so I leave it to you, dear reader, to decide for
    yourself.

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
            every n points to sample. Default: 5.

        colour : colormap or str, optional
            colourmap or colourmap name. Default: viridis.
        axcolour : str
            axis colour. Default: black.
        rasterise : bool, optional
            rasterise plot. Default: True.

        workers : int, optional
            number of workers for paralellised section. Default: 32.

    Returns:
        axes
            axes with projected "dispersion".
        colorbar
            colour bar for projected data.
    """

    from functools import partial
    from multiprocessing import Pool
    from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar
    import pymatgen.symmetry.analyzer as pmg
    from scipy.interpolate import interp1d

    # Phono3py
    tnames = tp.settings.to_tp()
    if quantity in tnames: quantity = tnames[quantity]
    if projected in tnames: projected = tnames[projected]
    qs = quantity if quantity == projected else [quantity, projected]
    data = tp.data.resolve.resolve(data, qs, temperature, direction)
    y = data[quantity]
    c = data[projected]
    qk = data['qpoint']

    # Phonopy
    qp = pdata['qpoint']
    x = pdata['x']

    qpi = [qp[i] for i in range(len(qp)) if i % interpolate == 0]
    xi = [x[i] for i in range(len(x)) if i % interpolate == 0]
    if xi[-1] != x[-1]: qpi.append(qp[-1]); xi.append(x[-1])

    # Map Phono3py to Phonopy
    struct = Poscar.from_file(poscar).structure
    sg = SpacegroupAnalyzer(struct)
    symops = sg.get_point_group_operations(cartesian=False)
    geq = partial(get_equivalent_qpoint, np.array(qk), symops)
    with Pool(processes=workers) as pool:
        min_id = pool.map(geq, qpi)
    y2 = y[min_id, :]
    c2 = c[min_id, :]

    # Interpolate
    x2 = np.linspace(min(xi), max(xi), 2500)
    yinterp = interp1d(xi, y2, kind='cubic', axis=0)
    y2 = np.abs(yinterp(x2))
    ysort = np.ravel(y2)
    ysort = ysort[ysort.argsort()]
    ymin = ysort[int(round(len(ysort)/100, 0))]
    ymax = ysort[-1]

    cinterp = interp1d(xi, c2, kind='cubic', axis=0)
    c2 = np.abs(cinterp(x2))
    csort = np.ravel(c2)
    csort = csort[csort.argsort()]
    cmin = csort[int(round(len(csort)/100, 0))]
    cmax = csort[-1]
    cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    for n in range(len(y2[0])):
        line = ax.scatter(x2, y2[:,n], c=c2[:,n], cmap=colour,
                          norm=cnorm, marker='.', s=1, rasterized=rasterise)

    if main:
        axlabels = tp.settings.labels()
        spinewidth = ax.spines['bottom'].get_linewidth()
        for d in pdata['tick_position']:
            ax.axvline(d, color=axcolour, linewidth=spinewidth)

        cbar = plt.colorbar(line)
        cbar.set_label(axlabels[projected])

        ax.set_xlabel(axlabels['wavevector'])
        ax.set_ylabel(axlabels[quantity])

        ax.set_xticks(pdata['tick_position'])
        ax.set_xticklabels(pdata['tick_label'])
        ax.tick_params(axis='x', which='minor', top=False, bottom=False)
        ax.set_xlim(x2[0], x2[-1])

        ax.set_ylim(ymin, ymax)
        if log:
            ax.set_yscale('log')
            ax.yaxis.set_major_locator(tp.settings.locator()['log'])
        else:
            ax.yaxis.set_major_locator(tp.settings.locator()['major'])
            ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])

    return ax, cbar

def add_wideband(ax, data, pdata, temperature=300, poscar='POSCAR', main=True,
                 interpolate=5, colour='viridis', rasterise=True, workers=32):
    """Plots a phonon dispersion with broadened bands.

    Arguments:
        ax : axes
            axes to plot on.
        data : dict
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
        interpolate : int, optional
            every n points to sample. Default: 5.

        colour : colormap or str, optional
            colourmap or colourmap name. Default: viridis.
        rasterise : bool, optional
            rasterise plot. Default: True.

        workers : int, optional
            number of workers for paralellised section. Default: 32.

    Returns:
        axes
            axes with widebands.
    """

    from functools import partial
    from multiprocessing import Pool
    from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar
    import pymatgen.symmetry.analyzer as pmg
    from scipy.interpolate import interp1d
    from tp.calculate import lorentzian

    # Phono3py
    data = tp.data.resolve.resolve(data, 'gamma', temperature)
    c = np.array(data['gamma'])
    qk = data['qpoint']

    # Phonopy
    qp = pdata['qpoint']
    x = pdata['x']

    qpi = [qp[i] for i in range(len(qp)) if i % interpolate == 0]
    xi = [x[i] for i in range(len(x)) if i % interpolate == 0]
    if xi[-1] != x[-1]: qpi.append(qp[-1]); xi.append(x[-1])

    # Map Phono3py to Phonopy
    struct = Poscar.from_file(poscar).structure
    sg = SpacegroupAnalyzer(struct)
    symops = sg.get_point_group_operations(cartesian=False)
    geq = partial(get_equivalent_qpoint, np.array(qk), symops)
    with Pool(processes=workers) as pool:
        min_id = pool.map(geq, qpi)
    c2 = c[min_id, :]

    x, indices = np.unique(x, return_index=True)
    f = np.array(pdata['frequency'])[indices]

    # Interpolate
    x2 = np.linspace(min(x), max(x), 2500)
    finterp = interp1d(x, f, kind='cubic', axis=0)
    f = finterp(x2)

    cinterp = interp1d(xi, c2, kind='cubic', axis=0)
    c2 = np.abs(cinterp(x2))
    fmax = np.amax(np.add(f, c2))
    fmin = np.amin(np.subtract(f, c2))
    f2 = np.linspace(fmin, fmax, 2500)

    area = np.zeros((len(x2), len(f)))
    for q in range(len(area)):
        for b in range(len(c2[q])):
            area[q] = np.add(area[q], lorentzian(f2, f[q][b], c2[q][b]))

    cnorm = mpl.colors.LogNorm(vmin=np.amin(area), vmax=np.amax(area))

    ax.pcolormesh(x2, f2, np.transpose(area), cmap=colour, norm=cnorm,
                  rasterized=rasterise)

    if main:
        axlabels = tp.settings.labels()
        ax.set_xlabel(axlabels['wavevector'])
        ax.set_ylabel(axlabels['frequency'])

        ax.set_xticks(pdata['tick_position'])
        ax.set_xticklabels(pdata['tick_label'])
        ax.tick_params(axis='x', which='minor', top=False, bottom=False)
        ax.set_xlim(x2[0], x2[-1])

        ax.yaxis.set_major_locator(tp.settings.locator()['major'])
        ax.yaxis.set_minor_locator(tp.settings.locator()['minor'])
        if round(np.amin(f), 1) == 0:
            ax.set_ylim(bottom=0)

    return ax

def get_equivalent_qpoint(qk, symops, qp, tol=1e-2):
    """Finds the closest phono3py qpoint to a phonopy qpoint.

    Arguments:
        qk : array-like
            qpoint from the phono3py kappa file.
        symmops
            symmetry operations (e.g. from Pymatgen)
        qp : array-like
            qpoint from the phonon dispersion.

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
