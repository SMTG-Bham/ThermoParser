"""Provides a command line interface.

Functions
---------

    tp
        gen
            kpar
            kpoints
        get
            amset
            boltztrap
            phono3py
            zt
        run
            boltztrap
        save
            cumkappa
            kappa
            zt
        plot
            avg_rates
            cumkappa
            dos
            kappa
            kappa_target
            phonons
            transport
            waterfall
            wideband
            ztmap
"""

import click
import h5py
import matplotlib as mpl
import numpy as np
import tp
import warnings
from copy import deepcopy
from tp.cli.options import *

@click.group()
def tp_cli():
    """Command line tools for transport properties."""
    return


@tp_cli.group()
def gen():
    """Tools for generating calculation inputs."""
    return


@gen.command()
@kpoints_options
def kpar(kpoints, mesh, poscar):
    """Suggests KPAR values for VASP.

    Reads a KPOINTS file by default, but a different file or mesh can be
    specified.
    """

    if kpoints is not None:
        pass
    elif mesh != ():
        kpoints = mesh
    else:
        kpoints = 'KPOINTS'

    k = tp.setup.vasp.get_kpar(kpoints, poscar=poscar)
    ks = [str(kp) for kp in k]
    click.echo('KPAR = {} or {}'.format(','.join(ks[:-1]), ks[-1]))

    return


@gen.command()
@kpoints_options
@click.option('-z', '--zero-kpoints', '--zero-ibzkpt',
              help='IBZKPT file path. Overrides --zero-mesh.',
              type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('--zero-mesh',
              help='zero-weighted mesh. Overridden by --zero-kpoints.',
              nargs=3,
              type=int)
@click.option('-o', '--output',
              help='Output path.',
              type=click.Path(file_okay=True, dir_okay=False),
              default='KPOINTS',
              show_default=True)
def kpoints(kpoints, mesh, poscar, zero_kpoints, zero_mesh, output):
    """Generates KPOINTS file for VASP with 0-weighted KPOINTS.

    Accepts either meshes, which are not preferred as they may not be
    the same as VASP would choose, or IBZKPT files.
    """

    if kpoints is not None:
        pass
    elif mesh != ():
        kpoints = mesh
    else:
        kpoints = 'KPOINTS'

    if zero_kpoints is None and zero_mesh != ():
        zero_kpoints = zero_mesh

    tp.setup.vasp.get_kpoints(kpoints, zero_kpoints, poscar=poscar,
                              output=output)

    return



@tp_cli.group()
def get():
    """Tools for accessing data."""
    return


@get.command('amset')
@input_argument
@click.option('-q', '--quantity',
              help='Quantity to read.',
              default='conductivity',
              show_default=True)
@doping_type_option
@doping_option
@direction_option
@temperature_option
@click.option('--spin',
              help='Spin direction.',
              type=click.Choice(['up', 'down', 'avg'], case_sensitive=False),
              default='avg',
              show_default=True)
@click.option('-s', '--scattering',
              help='Scattering mechanism.',
              type=click.Choice(['ADP', 'IMP', 'PIE', 'POP', 'total'],
                                case_sensitive=False),
              default='total',
              show_default=True)
@click.option('-b', '--band',
              help='AMSET band index (1-indexed). Only for bands used in '
                   'AMSET postprocessing, not for the original calculation.',
              type=int,
              default=1,
              show_default=True)
@click.option('--kpoint',
              help='k-point coordinates.',
              type=click.FloatRange(0, 1),
              nargs=3,
              default=[0., 0., 0.],
              show_default=True)
@click.option('-p', '--poscar',
              help='POSCAR path. Required for --qpoint.',
              type=click.Path(file_okay=True, dir_okay=False),
              default='POSCAR',
              show_default=True)

def get_amset(filename, quantity, dtype, doping, direction, temperature, spin,
              scattering, band, kpoint, poscar):
    """Prints AMSET values at given conditions.

    Requires an AMSET transport or mesh file.
    """
    # I'm making a few of these for different inputs so it doesn't get
    # too unwieldly checking which source they're from, but they could
    # also be combined.

    if quantity == 'scattering_rates':
        quantity = 'weighted_rates'

    try:
        data = tp.data.load.amset(filename, quantity, doping=dtype)
    except UnicodeDecodeError:
        data = tp.data.load.amset_mesh(filename, quantity, spin=spin,
                                       doping=dtype)

    dims = data['meta']['dimensions'][quantity]
    if 'kpoint' in dims:
        from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
        from pymatgen.io.vasp.inputs import Poscar

        struct = Poscar.from_file(poscar).structure
        sg = SpacegroupAnalyzer(struct)
        symops = sg.get_point_group_operations(cartesian=False)
        ki = tp.plot.phonons.get_equivalent_qpoint(data['kpoint'], symops,
                                                   kpoint)
        kpt = data['kpoint'][ki]
    else:
        kpt = None
    if 'temperature' not in dims:
        temperature = None
    if 'doping' not in dims:
        doping = None

    data = tp.data.resolve.resolve(data, quantity, direction=direction,
                                   temperature=temperature, stype=scattering,
                                   dtype=dtype, doping=doping, kpoint=kpt)
    if 'band' in dims:
        data[quantity] = data[quantity][band-1]
        if len(str(band)) > 1 and str(band)[-2:] in ['11', '12', '13']:
            suffix = 'th'
        elif str(band)[-1] == '1':
            suffix = 'st'
        elif str(band)[-1] == '2':
            suffix = 'nd'
        elif str(band)[-1] == '3':
            suffix = 'rd'
        else:
            suffix = 'th'

    printq = ''
    for i, l in enumerate(quantity):
        if l == '_':
            printq += ' '
        else:
            printq += l
    if printq[-1] == 's':
        printq = printq[:-1]

    print('The {}'.format(printq), end=' ')
    end = ' '
    if 'band' in dims:
        print('of the {}{} band in AMSET'.format(band, suffix), end=end)
    if 'kpoint' in dims:
        k = data['kpoint'][ki]
        print('at ({:.3f}, {:.3f}, {:.3f})'.format(k[0], k[1], k[2]), end=end)
    if 3 in dims:
        print('in the {} direction'.format(direction), end=end)
        dims.remove(3)
        if 3 in dims:
            dims.remove(3)
    if 'stype' in dims:
        print('due to {} scattering'.format(scattering), end=end)
        dims.remove('stype')
    if 'temperature' in dims:
        if len(dims) == 2:
            end = ' and '
        print('at {:d} K'.format(int(data['meta']['temperature'])), end=end)
        dims.remove('temperature')
    if 'doping' in dims:
        print('{:.3e} carriers cm^-3'.format(data['meta']['doping']), end=' ')
    print('is {:.3} {}.'.format(data[quantity], tp.settings.units()[quantity]))

    return

@get.command('boltztrap')
@input_argument
@click.option('-q', '--quantity',
              help='Quantity to read.',
              default='conductivity',
              show_default=True)
@doping_type_option
@doping_option
@direction_option
@temperature_option

def get_boltztrap(filename, quantity, dtype, doping, direction, temperature):
    """Prints BoltzTraP values at given conditions.

    Requires a tp boltztrap.hdf5 file.
    """
    # I'm making a few of these for different inputs so it doesn't get
    # too unwieldly checking which source they're from, but they could
    # also be combined.

    data = tp.data.load.boltztrap(filename, quantity, doping=dtype)

    dims = data['meta']['dimensions'][quantity]
    if 'temperature' not in dims:
        temperature = None
    if 'doping' not in dims:
        doping = None

    data = tp.data.resolve.resolve(data, quantity, direction=direction,
                                   temperature=temperature, dtype=dtype,
                                   doping=doping)
    printq = ''
    for i, l in enumerate(quantity):
        if l == '_':
            printq += ' '
        else:
            printq += l
    if printq[-1] == 's':
        printq = printq[:-1]

    end = ' '
    print('The {}-type {}'.format(data['meta']['doping_type'], printq), end=' ')
    if 3 in dims:
        print('in the {} direction'.format(direction), end=end)
        dims.remove(3)
        if 3 in dims:
            dims.remove(3)
    if 'temperature' in dims:
        if len(dims) == 2:
            end = ' and '
        print('at {:d} K'.format(int(data['meta']['temperature'])), end=end)
        dims.remove('temperature')
    if 'doping' in dims:
        print('{:.3e} carriers cm^-3'.format(data['meta']['doping']), end=' ')
    print('is {:.3} {}.'.format(data[quantity], tp.settings.units()[quantity]))

    return


@get.command('phono3py')
@input_argument
@click.option('-q', '--quantity',
              help='Quantity to read.',
              default='lattice_thermal_conductivity',
              show_default=True)
@direction_option
@temperature_option
@click.option('-b', '--band',
              help='Phonon band index (1-indexed).',
              type=int,
              default=1,
              show_default=True)
@click.option('--qpoint',
              help='q-point coordinates.',
              type=click.FloatRange(0, 1),
              nargs=3,
              default=[0., 0., 0.],
              show_default=True)
@click.option('-p', '--poscar',
              help='POSCAR path. Required for --qpoint.',
              type=click.Path(file_okay=True, dir_okay=False),
              default='POSCAR',
              show_default=True)

def get_phono3py(filename, quantity, direction, temperature, band, qpoint,
                 poscar):
    """Prints Phono3py values at given conditions.

    Requires a Phono3py kappa file.
    """
    # I'm making a few of these for different inputs so it doesn't get
    # too unwieldly checking which source they're from, but they could
    # also be combined.

    data = tp.data.load.phono3py(filename, quantity)

    dims = data['meta']['dimensions'][quantity]
    if 'qpoint' in dims:
        from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
        from pymatgen.io.vasp.inputs import Poscar

        struct = Poscar.from_file(poscar).structure
        sg = SpacegroupAnalyzer(struct)
        symops = sg.get_point_group_operations(cartesian=False)
        qi = tp.plot.phonons.get_equivalent_qpoint(data['qpoint'], symops,
                                                   qpoint)
        qpt = data['qpoint'][qi]
    else:
        qpt = None
    if 'temperature' not in dims:
        temperature = None

    data = tp.data.resolve.resolve(data, quantity, direction=direction,
                                   temperature=temperature, qpoint=qpt)
    if 'band' in dims:
        data[quantity] = data[quantity][band-1]
        if len(str(band)) > 1 and str(band)[-2:] in ['11', '12', '13']:
            suffix = 'th'
        elif str(band)[-1] == '1':
            suffix = 'st'
        elif str(band)[-1] == '2':
            suffix = 'nd'
        elif str(band)[-1] == '3':
            suffix = 'rd'
        else:
            suffix = 'th'

    printq = ''
    for i, l in enumerate(quantity):
        if l == '_':
            printq += ' '
        else:
            printq += l
    if printq[-1] == 's':
        printq = printq[:-1]

    end = ' '
    print('The {}'.format(printq), end=end)
    if 'band' in dims:
        print('of the {}{} band'.format(band, suffix), end=end)
    if 'qpoint' in dims:
        q = data['qpoint'][qi]
        print('at ({:.3f}, {:.3f}, {:.3f})'.format(q[0], q[1], q[2]), end=end)
    if 3 in dims:
        print('in the {} direction'.format(direction), end=end)
    if 'temperature' in dims:
        print('at {:d} K'.format(int(data['meta']['temperature'])), end=end)
    print('is {:.3} {}.'.format(data[quantity], tp.settings.units()[quantity]))

    return


@get.command('zt')
@input_argument
@click.option('-k', '--kappa',
              help='Phono3py kappa-mxxx.hdf5.',
              type=click.Path(exists=True, file_okay=True, dir_okay=False))
@doping_type_option
@doping_option
@direction_option
@temperature_option
def get_zt(filename, kappa, dtype, doping, direction, temperature):
    """Calculates and prints the ZT at given conditions.

    Requires electronic input file, and preferably phononic input, else
    lattice thermal conductivity defaults to 1.
    """

    try:
        edata = tp.data.load.amset(filename)
    except UnicodeDecodeError:
        try:
            edata = tp.data.load.boltztrap(filename, doping=dtype)
        except Exception:
            data = h5py.File(filename, 'r')
            edata = dict(data)
            for key in edata.keys():
                if isinstance(edata[key], dict) and dtype in edata[key]:
                    edata[key] = edata[key][dtype][()]

    equants = ['conductivity', 'seebeck', 'electronic_thermal_conductivity']
    ltc = 'lattice_thermal_conductivity'

    if 'zt' in edata:
        edata = tp.data.resolve.resolve(edata, 'zt', direction=direction,
                                        doping=doping, temperature=temperature)
    elif kappa is not None:
        kdata = tp.data.load.phono3py(kappa)
        edata, kdata = tp.calculate.interpolate(edata, kdata, 'temperature',
                                                equants, ltc, kind='cubic')
        edata[ltc] = kdata[ltc]
        edata['meta']['dimensions'][ltc] = kdata['meta']['dimensions'][ltc]
        edata = tp.calculate.zt_fromdict(edata)
        edata = tp.data.resolve.resolve(edata, 'zt', direction=direction,
                                        doping=doping, temperature=temperature)
    else:
        warnings.warn('Lattice thermal conductivity set to 1. For a more '
                      'accurate calculation, pass a phono3py kappa file to -k.')
        edata[ltc] = np.ones((len(edata['temperature']), 3, 3))

        edata['meta']['dimensions'][ltc] = ['temperature']
        edata = tp.calculate.zt_fromdict(edata)
        edata = tp.data.resolve.resolve(edata, 'zt', direction=direction,
                                        doping=doping, temperature=temperature)

    zt = edata['zt']
    n = edata['meta']['doping']
    t = int(edata['meta']['temperature'])

    print('The ZT in the {} direction at {:d} K and {:.3e} carriers cm^-3 '
            'is {:.3f}'.format(direction, t, n, zt))

    return


@tp_cli.group()
def run():
    """Tools for transport properties postprocessing."""
    return


@run.command()

@click.option('--tmin',
              help='Set minimum temperature in K. Does not speed up '
                   'calculation, just discards data.',
              type=float)
@click.option('--tmax',
              help='Maximum temperature in K.',
              default=1000.,
              show_default=True)
@click.option('--tstep',
              help='Temperature step size in K.',
              default=50.,
              show_default=True)

@click.option('--dmin',
              help='Minimum doping in cm^-3.',
              default=1e18,
              show_default=True)
@click.option('--dmax',
              help='Maximum doping in cm^-3.',
              default=1e21,
              show_default=True)
@click.option('--ndope',
              help='number of doping concetrations to calculate.',
              default=16,
              show_default=True)

@click.option('-r', '--relaxation_time',
              help='Charge carrier relaxation time in s.',
              default=1e-14,
              show_default=True)
@click.option('-l', '--lpfac',
              help='Factor to interpolate the DoS density by.',
              default=10,
              show_default=True)

@click.option('-k', '--kmode',
              help='Kappa calculation method. Options:'
                   'boltzmann: boltztrap method (default); '
                   'wiedemann: Wiedemann-Franz law with constant L; '
                   'snyder:    Wiedemann-Franz law with L dependant on '
                              'Seebeck coefficient.',
              default='boltzmann',
              type=click.Choice(['boltzmann', 'wiedemann', 'snyder'],
                                case_sensitive=False),
              show_default=False)
@click.option('--run/--norun',
              help='Run BoltzTraP.  [default: run]',
              default=True,
              show_default=False)
@click.option('--analyse/--noanalyse',
              help='Analyse BoltzTraP.  [default: analyse]',
              default=True,
              show_default=False)

@click.option('--kpoints',
              help='Path to KPOINTS if there are zero-weighted kpoints.',
              default='KPOINTS',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              show_default=True)
@click.option('-v', '--vasprun',
              help='Path to vasprun.xml.',
              default='vasprun.xml',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              show_default=True)
@click.option('-o', '--output',
              help='Output filename.',
              default='boltztrap.hdf5',
              show_default=True)

def boltztrap(tmin, tmax, tstep, dmin, dmax, ndope, relaxation_time, lpfac,
              kmode, run, analyse, kpoints, vasprun, output):
    """Runs BoltzTraP and sends data to hdf5."""

    doping = np.geomspace(dmin, dmax, ndope)

    tp.data.run.boltztrap(tmax=tmax, tstep=tstep, tmin=tmin, doping=doping,
                          ke_mode=kmode, vasprun=vasprun, kpoints=kpoints,
                          relaxation_time=relaxation_time, lpfac=lpfac,
                          run=run, analyse=analyse, output=output)

    return



@tp_cli.group()
def save():
    """Tools for saving data."""
    return


@save.command('cumkappa')
@input_argument
@click.option('--mfp/--frequency',
              help='x-axis quantity.  [default: frequency]',
              default=False,
              show_default=False)
@direction_option
@temperature_option
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-cumkappa',
              show_default=True)
@click.option('-e', '--extension',
              help='File type.',
              type=click.Choice(['dat', 'csv'], case_sensitive=False),
              multiple=True,
              default=['dat'],
              show_default=True)
def save_cumkappa(filename, mfp, direction, temperature, output, extension):
    """Extracts cumulative kappa from Phono3py hdf5.

    Saves to dat and/ or csv.
    """

    tp.data.save.cumkappa(filename, mfp=mfp, direction=direction,
                          temperature=temperature, output=output,
                          extension=extension)
    for e in extension:
        click.echo('{}.{} written'.format(output, e))

    return


@save.command('kappa')
@input_argument
@directions_option
@click.option('-o', '--output',
              help='Output filenames, sans extension.',
              default='tp-kappa',
              show_default=True)
def save_kappa(filename, direction, output):
    """Extracts kappa from Phono3py hdf5.

    Saves to text file.
    """

    f = tp.data.load.phono3py(filename)

    units = tp.settings.units()
    header = 'T({})'.format(units['temperature'])
    data = [f['temperature']]
    for d in direction:
        aniso = tp.data.resolve.resolve(f, 'lattice_thermal_conductivity',
                                        direction=d)
        header += ' kappa_{}({})'.format(d, units['lattice_thermal_conductivity'])
        data.append(aniso['lattice_thermal_conductivity'])

    np.savetxt('{}.dat'.format(output), np.transpose(data), header=header)
    click.echo('{}.dat written'.format(output))

    return


@save.command('zt')
@input_argument
@click.option('-k', '--kappa',
              help='Phono3py kappa-mxxx.hdf5.',
              type=click.Path(exists=True, file_okay=True, dir_okay=False))
@doping_type_option
@direction_option
@interpolate_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-zt',
              show_default=True)
def save_zt(filename, kappa, dtype, direction, interpolate, kind, output):
    """Extracts ZT from electronic transport data file.

    Saves ZT to hdf5 and highlights to yaml, and prints max to stdout.
    Currently accepts AMSET transport json or BoltzTraP hdf5, along with
    Phono3py hdf5 for lattice thermal conductivity.
    """

    tp.data.save.zt(filename, kappa, direction=direction, doping=dtype,
                    tinterp=interpolate, dinterp=interpolate, kind=kind,
                    output=output)
    click.echo('{0}.yaml and {0}.hdf5 written'.format(output))

    return



@tp_cli.group()
def plot():
    """Tools for plotting."""
    return


@plot.command()
@inputs_argument
@click.option('--total/--nototal',
              help='Plot total scattering rate  [default: total]',
              default=True,
              show_default=False)
@click.option('-x', '--x',
              help='x-axis variable.',
              type=click.Choice(['temperature', 'doping', 'both'],
                   case_sensitive=False),
              default='both',
              show_default=True)
@click.option('--crt',
        help='Constant relaxation time value.  [default: off]',
              type=float)
@doping_option
@temperature_option

@click.option('-c', '--colour',
              help='Colourmap name or list of colours in the order '
                   'plotted. Total and CRT are plotted last.',
              multiple=True,
              default=['tab10'],
              show_default=True)
@line_options

@xy_limit_options
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-avg-rates',
              show_default=True)
@verbose_option

def avg_rates(filenames, total, x, crt, doping, temperature, colour, linestyle,
              marker, xmin, xmax, ymin, ymax, style, large, extension, output,
              verbose):
    """Plots averaged scattering rates.

    Requires AMSET mesh files. Plots scattering rates averaged over
    kpoints and weighted by the derivative of the Fermi-Dirac
    distribution against tempearture or carrier concentration or both.
    If one file is given, it will be used for both plots, or if two are
    specified the one with the most temperatures will be used for the
    temperature plot and the one with the most carrier concentrations
    will be used for the doping plot. x-limits only work for individual
    plots.
    """

    if x == 'both':
        axes = tp.axes.large if large else tp.axes.small
        fig, ax, add_legend = axes.two_h(style)
    else:
        axes = tp.axes.large if large else tp.axes.small
        fig, ax, add_legend = axes.one(style)
    tax = ax[0] if x == 'both' else ax
    dax = ax[1] if x == 'both' else ax

    if len(filenames) == 1:
        data = [tp.data.load.amset_mesh(filenames[0], 'weighted_rates')]
        tindex = 0
        dindex = 0
    elif len(filenames) == 2:
        data = [tp.data.load.amset_mesh(f, 'weighted_rates') for f in filenames]
        tindex = 0 if len(data[0]['temperature']) > len(data[1]['temperature']) else 1
        dindex = 0 if len(data[0]['doping']) > len(data[1]['doping']) else 1
        if verbose:
            print('Using {} for the temperature data.'.format(filenames[tindex]))
            print('Using {} for the doping data.'.format(filenames[dindex]))
    else:
        raise Exception('Please specify at most two filenames')

    nlines = np.amax([len(d['stype']) for d in data])
    if total:
        nlines += 1
    if crt is not None:
        nlines += 1

    linestyle = list(linestyle)
    colour = list(colour)
    marker = list(marker)
    if len(colour) == 1:
        colour = colour[0]

    try:
        colours = mpl.cm.get_cmap(colour)(np.linspace(0, 1, nlines))
        colours = [c for c in colours]
    except ValueError:
        if isinstance(colour, str) and colour == 'skelton':
            colour = tp.plot.colour.skelton()
            colours = [colour(i) for i in np.linspace(0, 1, nlines)]
        else:
            colours = colour

    while len(colours) < nlines:
        colours.append(colours[-1])
    while len(linestyle) < nlines:
        linestyle.append('solid')
    while len(marker) < nlines:
        marker.append(None)

    labels = tp.settings.large_labels() if large else tp.settings.labels()
    if x == 'temperature' or x == 'both':
        tdata = tp.data.resolve.resolve(data[tindex], 'weighted_rates',
                                        doping=doping)
        if verbose:
            print('Using {} {}.'.format(tdata['meta']['doping'],
                                        tdata['meta']['units']['doping']))
        for i, rate in enumerate(data[tindex]['stype']):
            tax.plot(tdata['temperature'], tdata['weighted_rates'][i],
                     color=colours[i], linestyle=linestyle[i],
                     marker=marker[i], label=rate)
        if total:
            i += 1
            totrate = np.sum(tdata['weighted_rates'], axis=0)
            tax.plot(tdata['temperature'], totrate, color=colours[i],
                     linestyle=linestyle[i], marker=marker[i], label='Total')
        if crt is not None:
            i += 1
            crtrate = np.full(len(tdata['temperature']), crt)
            tax.plot(tdata['temperature'], crtrate, color=colours[i],
                     linestyle=linestyle[i], marker=marker[i], label='CRT')
        tax.set_xlabel(labels['temperature'])
        tax.set_ylabel(labels['scattering_rates'])
        tp.plot.utilities.set_locators(tax, x='linear', y='log')

    if x == 'doping' or x == 'both':
        ddata = tp.data.resolve.resolve(data[dindex], 'weighted_rates',
                                        temperature=temperature)
        if verbose:
            print('Using {} {}.'.format(ddata['meta']['temperature'],
                                        ddata['meta']['units']['temperature']))
        for i, rate in enumerate(data[dindex]['stype']):
            dax.plot(np.abs(ddata['doping']), ddata['weighted_rates'][i],
                     color=colours[i], linestyle=linestyle[i],
                     marker=marker[i], label=rate)
        if total:
            i += 1
            totrate = np.sum(ddata['weighted_rates'], axis=0)
            dax.plot(np.abs(ddata['doping']), totrate, color=colours[i],
                     linestyle=linestyle[i], marker=marker[i], label='Total')
        if crt is not None:
            i += 1
            crtrate = np.full(len(ddata['doping']), crt)
            dax.plot(np.abs(ddata['doping']), crtrate, color=colours[i],
                     linestyle=linestyle[i], marker=marker[i], label='CRT')
        dax.set_xlabel(labels['doping'])
        dax.set_ylabel(labels['scattering_rates'])
        tp.plot.utilities.set_locators(dax, x='log', y='log')

    add_legend(title='Rate')

    if x != 'both':
        if xmin is not None:
            if xmax is not None:
                ax.set_xlim(xmin, xmax)
            else:
                ax.set_xlim(left=xmin)
        elif xmax is not None:
            ax.set_xlim(right=xmax)

    if ymin is not None:
        if ymax is not None:
            for a in ax:
                a.set_ylim(ymin, ymax)
        else:
            for a in ax:
                ax.set_ylim(bottom=ymin)
    elif ymax is not None:
        for a in ax:
            ax.set_ylim(top=ymax)

    for ext in extension:
        fig.savefig('{}.{}'.format(output, ext))

    return


@plot.command()
@inputs_argument
@click.option('--mfp/--frequency',
              help='x-axis quantity.  [default: frequency]',
              default=False,
              show_default=False)
@click.option('-p', '--percent/--raw',
              help='Plot kappa in percent of total.  [default: raw]',
              default=False,
              show_default=False)

@directions_option
@temperature_option
@click.option('--xmin', metavar='xmin',
              help='Override minimum x.',
              type=float)
@click.option('--minkappa',
              help='Minimum kappa to plot for --mfp in percent.',
              type=click.FloatRange(0, 100),
              default=1.,
              show_default=True)

@click.option('-c', '--colour',
              help='Colour(s).',
              multiple=True,
              default=None)
@fill_options
@line_options

@xy_limit_options
@legend_options
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-cumkappa',
              show_default=True)
@verbose_option

def cumkappa(filenames, mfp, percent, direction, temperature, minkappa, colour,
             fill, fillalpha, line, linestyle, marker, xmin, xmax, ymin, ymax,
             label, legend_title, location, style, large, extension, output,
             verbose):
    """Plots cumulative kappa against frequency or mean free path.

    Reads Phono3py hdf5. Properties such as colour and linestyle loop,
    so if you have two data files and two directions, only two colours
    need to be specified, one for each direction; however if you want
    one for each datafile, each will need to be repeated twice
    consecutively.
    """

    linestyle = list(linestyle)
    colour = list(colour) if len(colour) > 0 else None
    label = list(label) if len(label) > 0 else None
    marker = list(marker) if len(marker) > 0 else None

    axes = tp.axes.large if large else tp.axes.small
    fig, ax, add_legend = axes.one(style)

    if mfp:
        data = [tp.data.load.phono3py(f, ['mode_kappa', 'mean_free_path']) for f in filenames]
        tp.plot.mfp.add_cum_kappa(ax, data, temperature=temperature,
                                  direction=direction, colour=colour,
                                  fill=fill, fillcolour=fillalpha, line=line,
                                  kmin=minkappa, scale=percent, label=label,
                                  linestyle=linestyle, marker=marker,
                                  verbose=verbose)
    else:
        data = [tp.data.load.phono3py(f, ['mode_kappa', 'frequency']) for f in filenames]
        tp.plot.frequency.add_cum_kappa(ax, data, temperature=temperature,
                                        direction=direction, colour=colour,
                                        fill=fill, fillcolour=fillalpha,
                                        line=line, scale=percent, label=label,
                                        marker=marker, linestyle=linestyle,
                                        verbose=verbose)

    if large:
        if mfp:
            ax.set_xlabel(tp.settings.large_labels()['mean_free_path'])
        if percent:
            ax.set_ylabel(tp.settings.large_labels()['cumulative_percent'])
        else:
            ax.set_ylabel(tp.settings.large_labels()['cumulative_kappa'])

    if label is not None:
        if location is None:
            add_legend(title=legend_title)
        else:
            add_legend(title=legend_title, location=location)

    if xmin is not None:
        if xmax is not None:
            ax.set_xlim(xmin, xmax)
        else:
            ax.set_xlim(left=xmin)
    elif xmax is not None:
        ax.set_xlim(right=xmax)

    if ymin is not None:
        if ymax is not None:
            ax.set_ylim(ymin, ymax)
        else:
            ax.set_ylim(bottom=ymin)
    elif ymax is not None:
        ax.set_ylim(top=ymax)

    for ext in extension:
        fig.savefig('{}.{}'.format(output, ext))

    return


@plot.command()
@input_argument
@dos_options
@click.option('-c', '--colour',
              help='Colour(s) in POSCAR order with total at the end or '
                   'colourmap name. If --notprojected, a single colour '
                   'can be specified. Total colour is overridden by '
                   '--total-colour.',
              multiple=True,
              default=['tab10'],
              show_default=True)
@fill_options
@line_options

@xy_limit_options
@click.option('--legend-title',
              help='Legend title. Accepts maths notation.')
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-dos',
              show_default=True)

def dos(filename, poscar, atoms, projected, total, total_label, total_colour,
        colour, fill, fillalpha, line, linestyle, marker, xmin, xmax, ymin,
        ymax, legend_title, location, style, large, extension, output):
    """Plots a phonon density of states."""

    axes = tp.axes.large if large else tp.axes.small
    fig, ax, add_legend = axes.one(style)

    linestyle = list(linestyle)
    colour = list(colour)
    marker = list(marker)
    if len(colour) == 1:
        colour = colour[0]

    data = tp.data.load.phonopy_dos(filename, poscar, atoms)
    tp.plot.frequency.add_dos(ax, data, projected=projected, total=total,
                              totallabel=total_label, colour=colour,
                              totalcolour=total_colour, fill=fill,
                              fillalpha=fillalpha, line=line,
                              linestyle=linestyle, marker=marker)
    if location is None:
        add_legend(title=legend_title)
    else:
        add_legend(title=legend_title, location=location)

    if xmin is not None:
        if xmax is not None:
            ax.set_xlim(xmin, xmax)
        else:
            ax.set_xlim(left=xmin)
    elif xmax is not None:
        ax.set_xlim(right=xmax)

    if ymin is not None:
        if ymax is not None:
            ax.set_ylim(ymin, ymax)
        else:
            ax.set_ylim(bottom=ymin)
    elif ymax is not None:
        ax.set_ylim(top=ymax)

    for ext in extension:
        fig.savefig('{}.{}'.format(output, ext))

    return


@plot.command()
@click.option('-k', '--kfile',
              help='Thermal data filename(s). Required for a '
                   '--component of lattice or total.',
              multiple=True)
@click.option('-e', '--efile',
              help='Electronic data filename(s). Required for a '
                   '--component of electronic or total.',
              multiple=True)
@click.option('--component',
              help='Thermal conductivity component.',
              default='lattice',
              type=click.Choice(['lattice', 'electronic', 'total'],
                                case_sensitive=False),
              show_default=True)
@directions_option
@click.option('--tmin',
              help='Minimum temperature to plot, by default in K.',
              default=300.,
              show_default=True)
@click.option('--tmax',
              help='Maximum temperature to plot, by default in K.',
              default=np.inf,
              show_default=False)
@doping_type_option
@doping_option

@click.option('-c', '--colour',
              help='Colourmap name or min and max colours or list of '
                   'colours.',
              multiple=True,
              default=['tab10'],
              show_default=True)
@line_options

@xy_limit_options
@auto_legend_options
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-kappa',
              show_default=True)

def kappa(kfile, efile, component, direction, tmin, tmax, dtype, doping,
          colour, linestyle, marker, xmin, xmax, ymin, ymax, label,
          legend_title, legend, location, style, large, extension, output):
    """Plots line graphs of thermal conductivity against temperature.

    Currently not all combinations of inputs work. If multiple --direction
    are specified, only the first file will be read. If --component is
    total and one --direction but multiple sets of files are specified,
    either there must be the same number of each, or only one file of
    one type (electronic or phononic), in which case it is used for all
    instances of the other.
    """
    # Future: If multiple components are specified, only one set of
    # input files are accepted.?

    linestyle = list(linestyle)
    marker = list(marker)
    label = list(label)
    colour = list(colour)

    tc = 'thermal_conductivity'
    etc = 'electronic_thermal_conductivity'
    ltc = 'lattice_thermal_conductivity'

    axes = tp.axes.large if large else tp.axes.small
    fig, ax, add_legend = axes.one(style)

    if component in ['electronic', 'total']:
        if len(efile) != 0:
            edata = []
            for f in efile:
                try:
                    edata.append(tp.data.load.amset(f))
                except UnicodeDecodeError:
                    try:
                        edata.append(tp.data.load.boltztrap(f, doping=dtype))
                    except Exception:
                        data = h5py.File(f, 'r')
                        edata.append(dict(data))
                        for key in edata[-1].keys():
                            if isinstance(edata[-1][key], dict) and \
                               dtype in edata[-1][key]:
                                edata[-1][key] = edata[-1][key][dtype][()]
        else:
            raise Exception('--efile must be specified for a '
                            '--component of electronic or total.')
    if component in ['lattice', 'total']:
        if len(kfile) != 0:
            kdata = []
            for f in kfile:
                kdata.append(tp.data.load.phono3py(f))
        else:
            raise Exception('--kfile must be specified for a '
                            '--component of lattice or total.')

    data = []
    defleg = {'labels': [], 'title': None} # default legend 
    if component == 'total':
        q = tc
        if len(direction) > 1:
            defleg['title'] = 'Direction'
            defleg['labels'] = direction
            for d in direction:
                kdata2 = tp.data.resolve.resolve(kdata[0], ltc, direction=d)
                edata2 = tp.data.resolve.resolve(edata[0], etc, doping=doping,
                                                 direction=d)
                kdata2, edata2 = tp.calculate.interpolate(kdata2, edata2,
                                                          'temperature', ltc,
                                                          etc, kind='cubic')
                data.append({'temperature': kdata2['temperature'],
                             tc:            kdata2[ltc] + edata2[etc]})
        elif len(kdata) == len(edata):
            defleg['title'] = 'Phononic Data'
            defleg['labels'] = kfile
            for i in range(len(kdata)):
                kdata[i] = tp.data.resolve.resolve(kdata[i], ltc,
                                                   direction=direction)
                edata[i] = tp.data.resolve.resolve(edata[i], etc, doping=doping,
                                                   direction=direction)
                kdata[i], edata[i] = tp.calculate.interpolate(kdata[i],
                                                              edata[i],
                                                              'temperature',
                                                              ltc, etc,
                                                              kind='cubic')
                data.append({'temperature': kdata[i]['temperature'],
                             tc:            kdata[i][ltc] + edata[i][etc]})
        elif len(kdata) == 1:
            defleg['title'] = 'Electronic Data'
            defleg['labels'] = efile
            kdata[0] = tp.data.resolve.resolve(kdata[0], ltc,
                                               direction=direction)
            for i in range(len(edata)):
                edata[i] = tp.data.resolve.resolve(edata[i], etc, doping=doping,
                                                   direction=direction)
                kdata2 = kdata[0] # in case of different-sized arrays
                kdata2, edata[i] = tp.calculate.interpolate(kdata2, edata[i],
                                                            'temperature',
                                                            ltc, etc,
                                                            kind='cubic')
                data.append({'temperature': kdata[i]['temperature'],
                             tc:            kdata2[ltc] + edata[i][etc]})
        elif len(edata) == 1:
            defleg['title'] = 'Phononic Data'
            defleg['labels'] = kfile
            edata[0] = tp.data.resolve.resolve(edata[0], etc, doping=doping,
                                               direction=direction)
            for i in range(len(kdata)):
                kdata[i] = tp.data.resolve.resolve(kdata[i], ltc,
                                                   direction=direction)
                edata2 = edata[0]
                kdata[i], edata2 = tp.calculate.interpolate(kdata[i], edata2,
                                                            'temperature',
                                                            ltc, etc,
                                                            kind='cubic')
                data.append({'temperature': kdata[i]['temperature'],
                             tc:            kdata[i][ltc] + edata2[etc]})
    elif component == 'lattice':
        q = ltc
        if len(direction) > 1:
            defleg['title'] = 'Direction'
            defleg['labels'] = direction
            for d in direction:
                kdata2 = tp.data.resolve.resolve(kdata[0], ltc, direction=d)
                data.append({'temperature': kdata2['temperature'],
                             tc:            kdata2[ltc]})
        else:
            defleg['title'] = 'Phononic Data'
            defleg['labels'] = kfile
            for i in range(len(kdata)):
                kdata[i] = tp.data.resolve.resolve(kdata[i], ltc,
                                                   direction=direction)
                data.append({'temperature': kdata[i]['temperature'],
                             tc:            kdata[i][ltc]})
    elif component == 'electronic':
        q = etc
        if len(direction) > 1:
            defleg['title'] = 'Direction'
            defleg['labels'] = direction
            for d in direction:
                edata2 = tp.data.resolve.resolve(edata[0], etc, direction=d)
                data.append({'temperature': edata2['temperature'],
                             tc:            edata2[etc]})
        else:
            defleg['title'] = 'Electronic Data'
            defleg['labels'] = efile
            for i in range(len(edata)):
                edata[i] = tp.data.resolve.resolve(edata[i], etc, doping=doping,
                                                   direction=direction)
                data.append({'temperature': edata[i]['temperature'],
                             tc:            edata[i][etc]})

    if label == []:
        label = defleg['labels']
        if legend_title is None:
            legend_title = defleg['title']

    try:
        colours = mpl.cm.get_cmap(colour[0])(np.linspace(0, 1, len(data)))
        colours = [c for c in colours]
    except ValueError:
        if isinstance(colour[0], str) and colour[0] == 'skelton':
            colour = tp.plot.colour.skelton()
            colours = [colour(i) for i in np.linspace(0, 1, len(data))]
        elif len(colour) == 2 and len(data) != 2:
            colour = tp.plot.colour.linear(*colour)
            colours = [colour(i) for i in np.linspace(0, 1, len(data))]
        else:
            colours = colour

    while len(colours) < len(data):
        colours.append(colours[-1])
    while len(linestyle) < len(data):
        linestyle.append('solid')
    while len(marker) < len(data):
        marker.append(None)
    while len(label) < len(data):
        label.append(None)

    for i in range(len(data)):
        j = np.where((data[i]['temperature'] <= tmax)
                   & (data[i]['temperature'] >= tmin))[0]

        ax.plot(data[i]['temperature'][j], data[i][tc][j],
                label='${}$'.format(label[i]), linestyle=linestyle[i],
                marker=marker[i], c=colours[i])

    if xmin is not None:
        if xmax is not None:
            ax.set_xlim(xmin, xmax)
        else:
            ax.set_xlim(left=xmin)
    elif xmax is not None:
        ax.set_xlim(right=xmax)

    if ymin is not None:
        if ymax is not None:
            ax.set_ylim(ymin, ymax)
        else:
            ax.set_ylim(bottom=ymin)
    elif ymax is not None:
        ax.set_ylim(top=ymax)

    axlabels = tp.settings.large_labels() if large else tp.settings.labels()
    ax.set_xlabel(axlabels['temperature'])
    ax.set_ylabel(axlabels[q])
    tp.plot.utilities.set_locators(ax, 'linear', 'linear')
    if legend:
        if location is None:
            add_legend(title=legend_title)
        else:
            add_legend(title=legend_title, location=location)

    for ext in extension:
        fig.savefig('{}.{}'.format(output, ext))

    return


@plot.command()
@input_argument
@click.option('-z', '--zt',
              help='Target ZT.',
              type=click.FloatRange(0),
              default=2.,
              show_default=True)

@direction_option
@interpolate_options

@click.option('-c', '--colour',
              help='Colourmap name or #rrggbb highlight colour or min '
                   'and max and highlight #rrggbb colours to generate '
                   'a colourmap from.',
              multiple=True,
              default=['viridis'],
              show_default=True)
@click.option('-n', '--negativecolour',
              help='Colour for values below --cmin.',
              default='grey',
              show_default=True)

@xy_limit_options
@c_limit_options
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-kappa-target',
              show_default=True)

def kappa_target(filename, zt, direction, interpolate, kind, colour,
                 negativecolour, xmin, xmax, ymin, ymax, cmin, cmax, style,
                 large, extension, output):
    """Plots lattice thermal conductivity to reach a target ZT.

    Currently accepts AMSET transport json or BoltzTraP hdf5.
    """

    cmin = 0 if cmin is None else cmin
    if len(colour) == 1:
        colour = colour[0]
    else:
        colour = list(colour)

    axes = tp.axes.large if large else tp.axes.small
    try:
        edata = tp.data.load.amset(filename)
    except UnicodeDecodeError:
        edata = tp.data.load.boltztrap(filename)

    fig, ax, _ = axes.one_colourbar(style)

    cbar = tp.plot.heatmap.add_kappa_target(ax, edata, zt=zt,
                                            direction=direction,
                                            xinterp=interpolate,
                                            yinterp=interpolate,
                                            kind=kind, colour=colour,
                                            negativecolour=negativecolour,
                                            xmin=xmin, xmax=xmax, ymin=ymin,
                                            ymax=ymax, cmin=cmin, cmax=cmax)

    if large:
        cbar.set_label(tp.settings.large_labels()['lattice_thermal_conductivity'])

    for ext in extension:
        fig.savefig('{}.{}'.format(output, ext))

    return


@plot.command('phonons')
@inputs_argument
@click.option('--bandmin',
              help='Minimum band index.',
              type=click.IntRange(0))
@click.option('--bandmax',
              help='Maximum band index.',
              type=click.IntRange(0))

@click.option('-c', '--colour',
              help='Colourmap name or min and max colours or list of '
                   'colours.',
              multiple=True,
              default=['winter_r'],
              show_default=True)
@line_options
@click.option('--xmarkcolour',
              help='High-symmetry point line colour.',
              default='black',
              show_default=True)
@click.option('--xmarklinestyle',
              help='High-symmetry point linestyle.',
              default='solid',
              show_default=True)

@click.option('-d', '--dos',
              help='projected_dos.dat or equivalent for optional DoS plot.',
              type=click.Path(exists=True, file_okay=True, dir_okay=False))
@dos_options
@click.option('--doscolour',
              help='Colour(s) in POSCAR order with total at the end or '
                   'colourmap name. If --notprojected, a single colour '
                   'can be specified. Total colour is overridden by '
                   '--total-colour.',
              multiple=True,
              default=['tab10'],
              show_default=True)
@fill_options

@legend_options
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-phonons',
              show_default=True)

def converge_phonons(filenames, bandmin, bandmax, colour, linestyle, marker,
                     xmarkcolour, xmarklinestyle, dos, poscar, atoms,
                     projected, total, total_label, total_colour, doscolour,
                     fill, fillalpha, line, label, legend_title, location,
                     style, large, extension, output):
    """Plots and overlays phonon dispersions.

    Can have optional DoS on the side.
    """

    # Really the DoS bit should be in a chained command, but I'm having
    # trouble getting that working atm.

    linestyle = list(linestyle)
    if len(colour) == 1:
        colour = colour[0]
    else:
        colour = list(colour)
    if len(marker) == 0:
        marker = [None]
    else:
        marker = list(marker)
    if len(label) == 0:
        label = [None]
    else:
        label = list(label)
    if len(doscolour) == 1:
        doscolour = doscolour[0]
    else:
        doscolour = list(doscolour)

    data = [tp.data.load.phonopy_dispersion(f) for f in filenames]

    axes = tp.axes.large if large else tp.axes.small
    if dos is None:
        fig, ax, add_legend = axes.one(style)
    else:
        fig, ax, add_legend = axes.one_dos(style)
        dosax = ax[1]
        ax = ax[0]

    tp.plot.phonons.add_multi(ax, data, colour=colour, linestyle=linestyle,
                              marker=marker, label=label, bandmin=bandmin,
                              bandmax=bandmax,
                              xmarkkwargs={'color':     xmarkcolour,
                                           'linestyle': xmarklinestyle})
    if dos is not None:
        dosdata = tp.data.load.phonopy_dos(dos, poscar, atoms)
        tp.plot.frequency.add_dos(dosax, dosdata, projected=projected,
                                  total=total, totallabel=total_label,
                                  colour=doscolour, totalcolour=total_colour,
                                  fill=fill, fillalpha=fillalpha, line=line,
                                  invert=True)
        dosax.set_ylim(ax.get_ylim())
        tp.plot.utilities.set_locators(dosax, dos=True)

    if label != [None] or dos is not None:
        if location is None:
            add_legend(title=legend_title)
        else:
            add_legend(title=legend_title, location=location)

    for ext in extension:
        fig.savefig('{}.{}'.format(output, ext))

    return


@plot.command()
@inputs_argument
@click.option('-k', '--kfile',
              help='Thermal data filename(s). Required for a --quantity '
                   'of lattice_ or total_thermal_conductivity.',
              multiple=True)
@click.option('-q', '--quantity',
              help='Quantity(/ies) to plot. Max: 4.',
              multiple=True,
              default=['conductivity', 'seebeck',
                       'electronic_thermal_conductivity'],
              show_default=True)
@directions_option
@click.option('--tmin',
              help='Minimum temperature to plot, by default in K.',
              default=300.,
              show_default=True)
@click.option('--tmax',
              help='Maximum temperature to plot, by default in K.',
              default=np.inf,
              show_default=False)
@doping_type_option
@dopings_option

@click.option('-c', '--colour',
              help='Colourmap name or min and max colours or list of '
                   'colours.',
              multiple=True,
              default=['tab10'],
              show_default=True)
@line_options

@xy_limit_options
@auto_legend_options
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-transport',
              show_default=True)

def transport(filenames, kfile, quantity, direction, tmin, tmax, dtype, doping,
              colour, linestyle, marker, xmin, xmax, ymin, ymax, label,
              legend_title, legend, location, style, large, extension, output):
    """Plots line graphs of transport properties against temperature.

    Currently not all combinations of inputs work. The order of
    precedence is lines represent doping > direction > files.
    If using multiple sets of files (so only one doping and direction),
    either there must be the same number of each, or only one file of
    one type, in which case it is used for all instances of the other.

    Changing the y-axis limits is also not currently supported (we
    recommend changing the temperatures plotted instead).
    """

    if len(quantity) < 1 or len(quantity) > 4:
        raise Exception('--quantity must be between 1 and 4 items long.')
    tnames = tp.settings.to_tp()
    quantity = [tnames[q] if q in tnames else q for q in quantity]

    linestyle = list(linestyle)
    marker = list(marker)
    label = list(label)
    colour = list(colour)

    tc = 'thermal_conductivity'
    etc = 'electronic_thermal_conductivity'
    ltc = 'lattice_thermal_conductivity'

    axes = tp.axes.large if large else tp.axes.small
    axf = [axes.one, axes.two_h, axes.three_h, axes.four_square]
    axlabels = tp.settings.large_labels() if large else tp.settings.labels()
    fig, ax, add_legend =  axf[len(quantity) - 1](style)
    if len(quantity) == 4:
        ax = [ax[0][0], ax[0][1], ax[1][0], ax[1][1]]

    edata = []
    for f in filenames:
        try:
            edata.append(tp.data.load.amset(f))
        except UnicodeDecodeError:
            try:
                edata.append(tp.data.load.boltztrap(f, doping=dtype))
            except Exception:
                data = h5py.File(f, 'r')
                edata.append(dict(data))
                for key in edata[-1].keys():
                    if isinstance(edata[-1][key], dict) and \
                       dtype in edata[-1][key]:
                        edata[-1][key] = edata[-1][key][dtype][()]
    if ltc in quantity or tc in quantity:
        if len(kfile) != 0:
            kdata = []
            if len(filenames) != len(kfile) and len(filenames) != 1 and \
               len(kfile) != 1:
                raise Exception('Could not match electronic and phononic data. '
                                'either only one file should be specified for '
                                'one, or both should be the same length.')
            for f in kfile:
                kdata.append(tp.data.load.phono3py(f))
        else:
            raise Exception('--kfile must be specified for a '
                            '--quantity of {} or {}.'.format(ltc, tc))

    for e in edata:
        e['doping'] = np.abs(e['doping'])

    data = [[], [], [], []]
    defleg = {'labels': [], 'title': None} # default legend 
    # The data array could get quite large and memory inefficient, but
    # looping through the plot function at the end is more maintainable
    # compared to putting one in every if statement :shrug:
    lendata = 1
    for i, q in enumerate(quantity):
        if len(doping) > 1: # one line per doping
            lendata = len(doping)
            defleg['title'] = axlabels['doping']
            if q == tc:
                kdata2 = deepcopy(kdata[0])
                kdata2 = tp.data.resolve.resolve(kdata2, ltc,
                                                 direction=direction[0])
                dopelist = []
                for d in doping:
                    edata2 = deepcopy(edata[0])
                    edata2 = tp.data.resolve.resolve(edata2, etc, doping=d,
                                                     direction=direction[0])
                    dopelist.append('{:.2e}'.format(edata2['meta']['doping']))
                    kdata3, edata2 = tp.calculate.interpolate(kdata2, edata2,
                                                              'temperature',
                                                              ltc, etc,
                                                              kind='cubic')
                    data[i].append({'temperature': kdata3['temperature'],
                                    q:             kdata3[ltc] + edata2[etc]})
                defleg['labels'] = dopelist
            elif q in edata[0] and 'temperature' in edata[0]['meta']['dimensions'][q]:
                dopelist = []
                for d in doping:
                    edata2 = deepcopy(edata[0])
                    edata2 = tp.data.resolve.resolve(edata2, q, doping=d,
                                                     direction=direction[0])
                    dopelist.append('{:.2e}'.format(edata2['meta']['doping']))
                    data[i].append({'temperature': edata2['temperature'],
                                    q:             edata2[q]})
                defleg['labels'] = dopelist
            elif q in kdata[0] and 'temperature' in kdata[0]['meta']['dimensions'][q]:
                kdata2 = deepcopy(kdata[0])
                kdata2 = tp.data.resolve.resolve(kdata2, q,
                                                 direction=direction[0])
                data[i].append({'temperature': kdata2['temperature'],
                                q:             kdata2[q]})
        elif len(direction) > 1: # one line per direction
            lendata = len(direction)
            defleg['title'] = 'Direction'
            defleg['labels'] = direction
            if q == tc:
                for d in direction:
                    kdata2 = deepcopy(kdata[0])
                    kdata2 = tp.data.resolve.resolve(kdata2, ltc,
                                                     direction=d)
                    edata2 = deepcopy(edata[0])
                    edata2 = tp.data.resolve.resolve(edata2, etc,
                                                     doping=doping[0],
                                                     direction=d)
                    kdata2, edata2 = tp.calculate.interpolate(kdata2, edata2,
                                                              'temperature',
                                                              ltc, etc,
                                                              kind='cubic')
                    data[i].append({'temperature': kdata2['temperature'],
                                    q:             kdata2[ltc] + edata2[etc]})
            elif q in edata[0] and 'temperature' in edata[0]['meta']['dimensions'][q]:
                for d in direction:
                    edata2 = deepcopy(edata[0])
                    edata2 = tp.data.resolve.resolve(edata2, q,
                                                     doping=doping[0],
                                                     direction=d)
                    data[i].append({'temperature': edata2['temperature'],
                                    q:             edata2[q]})
            elif q in kdata[0] and 'temperature' in kdata[0]['meta']['dimensions'][q]:
                for d in direction:
                    kdata2 = deepcopy(kdata[0])
                    kdata2 = tp.data.resolve.resolve(kdata2, q,
                                                     direction=d)
                    data[i].append({'temperature': kdata2['temperature'],
                                    q:             kdata2[q]})
        else: # one line per file/ one line
            if q == tc:
                lendata = len(kdata)
                if len(kdata) == len(edata):
                    defleg['title'] = 'Electronic Data'
                    defleg['labels'] = filenames
                    for j in range(len(kdata)):
                        kdata2 = deepcopy(kdata[j])
                        edata2 = deepcopy(edata[j])
                        kdata2 = tp.data.resolve.resolve(kdata2, ltc,
                                                         direction=direction[0])
                        edata2 = tp.data.resolve.resolve(edata2, etc,
                                                         doping=doping[0],
                                                         direction=direction[0])
                        kdata2, edata2 = tp.calculate.interpolate(kdata2, edata2,
                                                                  'temperature',
                                                                  ltc, etc,
                                                                  kind='cubic')
                        data[i].append({'temperature': kdata2['temperature'],
                                        q:             kdata2[ltc] + edata2[etc]})
                elif len(kdata) == 1:
                    defleg['title'] = 'Electronic Data'
                    defleg['labels'] = filenames
                    lendata = len(edata)
                    kdata2 = deepcopy(kdata[0])
                    kdata2 = tp.data.resolve.resolve(kdata, ltc,
                                                     direction=direction[0])
                    for j in range(len(edata)):
                        edata2 = deepcopy(edata[j])
                        edata2 = tp.data.resolve.resolve(edata2, etc,
                                                         doping=doping[0],
                                                         direction=direction[0])
                        kdata2, edata2 = tp.calculate.interpolate(kdata2, edata2,
                                                                  'temperature',
                                                                  ltc, etc,
                                                                  kind='cubic')
                        data[i].append({'temperature': kdata2['temperature'],
                                        q:             kdata2[ltc] + edata2[etc]})
                elif len(edata) == 1:
                    defleg['title'] = 'Phononic Data'
                    defleg['labels'] = kfile
                    edata2 = deepcopy(edata[0])
                    edata2 = tp.data.resolve.resolve(edata2, etc,
                                                     doping=doping[0],
                                                     direction=direction[0])
                    for j in range(len(kdata)):
                        kdata2 = deepcopy(kdata[j])
                        kdata2 = tp.data.resolve.resolve(kdata2, ltc,
                                                         direction=direction[0])
                        kdata2, edata2 = tp.calculate.interpolate(kdata2, edata2,
                                                                  'temperature',
                                                                  ltc, etc,
                                                                  kind='cubic')
                        data[i].append({'temperature': kdata2['temperature'],
                                        q:             kdata2[ltc] + edata2[etc]})
            elif q in edata[0] and 'temperature' in edata[0]['meta']['dimensions'][q]:
                if len(edata) > 1:
                    lendata = len(edata)
                    defleg['title'] = 'Electronic Data'
                    defleg['labels'] = filenames
                elif defleg['title'] is None:
                    defleg['title'] = 'Electronic Data'
                    defleg['labels'] = filenames
                for j in range(len(edata)):
                    edata2 = deepcopy(edata[j])
                    edata2 = tp.data.resolve.resolve(edata2, q,
                                                     doping=doping[0],
                                                     direction=direction[0])
                    data[i].append({'temperature': edata2['temperature'],
                                    q:             edata2[q]})
            elif q in kdata[0] and 'temperature' in kdata[0]['meta']['dimensions'][q]:
                if len(kdata) > 1:
                    lendata = len(kdata)
                    defleg['title'] = 'Phononic Data'
                    defleg['labels'] = kfile
                elif defleg['title'] is None:
                    defleg['title'] = 'Phononic Data'
                    defleg['labels'] = kfile
                for j in range(len(kdata)):
                    kdata2 = deepcopy(kdata[j])
                    kdata2 = tp.data.resolve.resolve(kdata2, q,
                                                     direction=direction[0])
                    data[i].append({'temperature': kdata2['temperature'],
                                    q:             kdata2[q]})

    try:
        colours = mpl.cm.get_cmap(colour[0])(np.linspace(0, 1, lendata))
        colours = [c for c in colours]
    except ValueError:
        if isinstance(colour[0], str) and colour[0] == 'skelton':
            colour = tp.plot.colour.skelton()
            colours = [colour(i) for i in np.linspace(0, 1, lendata)]
        elif len(colour) == 2 and len(data) != 2:
            colour = tp.plot.colour.linear(*colour)
            colours = [colour(i) for i in np.linspace(0, 1, lendata)]
        else:
            colours = colour

    if label == []:
        label = defleg['labels']
        if legend_title is None:
            legend_title = defleg['title']

    while len(colours) < lendata:
        colours.append(colours[-1])
    while len(linestyle) < lendata:
        linestyle.append('solid')
    while len(marker) < lendata:
        marker.append(None)
    while len(label) < lendata:
        label.append(None)

    for i, d in enumerate(data):
        for j, d2 in enumerate(d):
            d2['temperature'] = np.array(d2['temperature'])
            k = np.where((d2['temperature'] <= tmax)
                       & (d2['temperature'] >= tmin))[0]

            ax[i].plot(d2['temperature'][k], d2[quantity[i]][k],
                       linestyle=linestyle[j], marker=marker[j], c=colours[j],
                       label="${}$".format(label[j]))

    for i, q in enumerate(quantity):
        ax[i].set_xlabel(axlabels['temperature'])
        ax[i].set_ylabel(axlabels[q])
        if len(doping) > 1 and q in ['conductivity', etc]:
            tp.plot.utilities.set_locators(ax[i], 'linear', 'log')
        else:
            tp.plot.utilities.set_locators(ax[i], 'linear', 'linear')
    if legend:
        if location is None:
            add_legend(title=legend_title)
        else:
            add_legend(title=legend_title, location=location)

    for a in ax:
        if xmin is not None:
            if xmax is not None:
                a.set_xlim(xmin, xmax)
            else:
                a.set_xlim(left=xmin)
        elif xmax is not None:
            a.set_xlim(right=xmax)

    #if ymin is not None:
    #    if ymax is not None:
    #        ax.set_ylim(ymin, ymax)
    #    else:
    #        ax.set_ylim(bottom=ymin)
    #elif ymax is not None:
    #    ax.set_ylim(top=ymax)

    for ext in extension:
        fig.savefig('{}.{}'.format(output, ext))

    return


@plot.command()
@input_argument
@click.option('-y',
              help='y-axis quantity. Options include frequency, kappa, '
                   'group_velocity, lifetime, mean_free_path, '
                   'occupation and ph_ph_strength.',
              default='kappa',
              show_default=True)
@click.option('-x',
              help='x-axis quantity. Options include frequency, kappa, '
                   'group_velocity, lifetime, mean_free_path, '
                   'occupation and ph_ph_strength.',
              default='frequency',
              show_default=True)
@click.option('-p', '--projected',
              help='Optional colour quantity. Options include density, '
                   'which colours by point density, or the usual '
                   'frequency, kappa, group_velocity, lifetime, '
                   'mean_free_path, occupation and ph_ph_strength.')

@direction_option
@temperature_option

@click.option('-c', '--colour',
              help='Colourmap name or colour or list of colours (one '
                   'per band). If --projected, only colourmap name.',
              multiple=True,
              default=['viridis'],
              show_default=True)
@click.option('-a', '--alpha',
              help='Marker opacity.',
              type=click.FloatRange(0., 1.),
              default=0.3,
              show_default=True)
@click.option('-l', '--linewidth',
              help='Marker edge linewidth.',
              type=click.FloatRange(0),
              default=0.,
              show_default=True)
@click.option('--edgecolour',
              help='Marker edge colour.',
              default='black',
              show_default=True)
@click.option('-m', '--marker',
              help='Marker shape.',
              default='.',
              show_default=True)
@click.option('--markersize',
              help='Marker size.',
              type=click.FloatRange(0),
              default=1.,
              show_default=True)

@click.option('--xscale',
              help='Override x-scale.',
              type=click.Choice(['linear', 'log'], case_sensitive=False))
@click.option('--yscale',
              help='Override y-scale.',
              type=click.Choice(['linear', 'log'], case_sensitive=False))
@click.option('--cscale',
              help='Override colour-scale if projected.',
              type=click.Choice(['linear', 'log'], case_sensitive=False))
@xy_limit_options
@c_limit_options

@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-waterfall',
              show_default=True)
@verbose_option

def waterfall(filename, y, x, projected, direction, temperature, colour, alpha,
              linewidth, edgecolour, marker, markersize, xscale, yscale,
              cscale, xmin, xmax, ymin, ymax, cmin, cmax, style, large,
              extension, output, verbose):
    """Plots 3rd-order phonon properties per band per q-point."""

    axes = tp.axes.large if large else tp.axes.small
    if len(colour) == 1:
        colour = colour[0]
    else:
        colour = list(colour)

    if x == 'kappa': x = 'mode_kappa'
    if y == 'kappa': y = 'mode_kappa'
    quantities = [x, y]
    if projected == 'density':
        pass
    elif projected is not None:
        if projected == 'kappa': projected = 'mode_kappa'
        quantities.append(projected)
    data = tp.data.load.phono3py(filename, quantities)

    if projected is None or projected == 'density':
        fig, ax, _ = axes.one(style)
    else:
        fig, ax, _ = axes.one_colourbar(style)

    if projected == 'density':
        cbar = tp.plot.frequency.add_density(ax, data, y, xquantity=x,
                                             temperature=temperature,
                                             direction=direction,
                                             colour=colour, alpha=alpha,
                                             linewidth=linewidth,
                                             edgecolors=edgecolour,
                                             marker=marker, s=markersize,
                                             verbose=verbose)
    elif projected is not None:
        cbar = tp.plot.frequency.add_projected_waterfall(ax, data, y,
                  projected, xquantity=x, temperature=temperature,
                  direction=direction, colour=colour, alpha=alpha,
                  linewidth=linewidth, edgecolors=edgecolour, marker=marker,
                  s=markersize, cmin=cmin, cmax=cmax, cscale=cscale,
                  verbose=verbose)
    else:
        tp.plot.frequency.add_waterfall(ax, data, y, xquantity=x,
                                        temperature=temperature,
                                        direction=direction, colour=colour,
                                        alpha=alpha, linewidth=linewidth,
                                        edgecolors=edgecolour, marker=marker,
                                        s=markersize, verbose=verbose)

    tp.plot.utilities.set_locators(ax, xscale, yscale)

    if large:
        ax.set_xlabel(tp.settings.large_labels()[x])
        ax.set_ylabel(tp.settings.large_labels()[y])
        if projected is not None:
            cbar.set_label(tp.settings.large_labels()[projected])

    if xmin is not None:
        if xmax is not None:
            ax.set_xlim(xmin, xmax)
        else:
            ax.set_xlim(left=xmin)
    elif xmax is not None:
        ax.set_xlim(right=xmax)

    if ymin is not None:
        if ymax is not None:
            ax.set_ylim(ymin, ymax)
        else:
            ax.set_ylim(bottom=ymin)
    elif ymax is not None:
        ax.set_ylim(top=ymax)

    for ext in extension:
        fig.savefig('{}.{}'.format(output, ext))

    return


@plot.command()
@click.argument('phonons',
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('kappa',
                type=click.Path(exists=True, file_okay=True, dir_okay=False))

@temperature_option
@click.option('-p', '--poscar',
              help='POSCAR path.',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              default='POSCAR',
              show_default=True)

@click.option('-c', '--colour',
              help='Colourmap name or max colour (fades to white) or '
                   'min and max colours to generate a colourmap from.',
              multiple=True,
              default=['viridis'],
              show_default=True)
@click.option('--smoothing',
              help='Every n points to sample.',
              default=5,
              show_default=True)

@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-wideband',
              show_default=True)
@verbose_option

def wideband(phonons, kappa, temperature, poscar, colour, smoothing, style,
             large, extension, output, verbose):
    """Plots a broadened phonon dispersion."""

    axes = tp.axes.large if large else tp.axes.small
    if len(colour) == 1:
        colour = colour[0]
    else:
        colour = list(colour)

    pdata = tp.data.load.phonopy_dispersion(phonons)
    kdata = tp.data.load.phono3py(kappa, 'wideband')

    fig, ax, _ = axes.one(style)

    tp.plot.phonons.add_wideband(ax, kdata, pdata, temperature=temperature,
                                 poscar=poscar, smoothing=smoothing,
                                 colour=colour, verbose=verbose)

    for ext in extension:
        fig.savefig('{}.{}'.format(output, ext))

    return


@plot.command()
@input_argument
@click.option('-k', '--kappa',
              help='Phono3py kappa hdf5. Ignored if ZT is in file. If '
                   'otherwise unspecified, set to 1 (W m-1 K-1).',
              type=click.Path(file_okay=True, dir_okay=False))

@direction_option
@doping_type_option
@interpolate_options

@click.option('-c', '--colour',
              help='Colourmap name or #rrggbb highlight colour or min '
                   'and max and highlight #rrggbb colours to generate '
                   'a colourmap from.',
              multiple=True,
              default=['viridis'],
              show_default=True)

@xy_limit_options
@c_limit_options
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-ztmap',
              show_default=True)

def ztmap(filename, kappa, direction, dtype, interpolate, kind, colour, xmin,
          xmax, ymin, ymax, cmin, cmax, style, large, extension, output):
    """Plots ZT against temperature and carrier concentration."""

    axes = tp.axes.large if large else tp.axes.small
    if len(colour) == 1:
        colour = colour[0]

    try:
        edata = tp.data.load.amset(filename)
    except UnicodeDecodeError:
        try:
            edata = tp.data.load.boltztrap(filename, doping=dtype)
        except Exception:
            data = h5py.File(filename, 'r')
            edata = dict(data)
            for key in edata.keys():
                if isinstance(edata[key], dict) and dtype in edata[key]:
                    edata[key] = edata[key][dtype][()]

    if kappa is not None:
        kdata = tp.data.load.phono3py(kappa)
    else:
        kdata = None

    fig, ax, _ = axes.one_colourbar(style)

    tp.plot.heatmap.add_ztmap(ax, edata, kdata=kdata, direction=direction,
                              xinterp=interpolate, yinterp=interpolate,
                              kind=kind, colour=colour, xmin=xmin, xmax=xmax,
                              ymin=ymin, ymax=ymax, cmin=cmin, cmax=cmax)

    for ext in extension:
        fig.savefig('{}.{}'.format(output, ext))

    return



# Welcome to the crystal ball :ghost:
#@plot.group(chain=True)
#def frequency():
#    """Chainable tools for plotting."""
#    return
