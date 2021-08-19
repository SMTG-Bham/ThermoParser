"""Provides a command line interface."""

import click
import matplotlib.pyplot as plt
import numpy as np
import tp
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
    """Run BoltzTraP and sends data to hdf5."""

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
def save_cumkappa(filename, mfp, direction, temperature, output):
    """Extracts cumulative kappa from Phono3py hdf5.

    Saves to text file.
    """

    tp.data.save.cumkappa(filename, mfp=mfp, direction=direction,
                          temperature=temperature, output=output)
    click.echo('{}.dat written'.format(output))

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
@line_fill_options

@xy_limit_options
@legend_options
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-cumkappa',
              show_default=True)

def cumkappa(filenames, mfp, percent, direction, temperature, minkappa, colour,
             fill, fillalpha, line, linestyle, marker, xmin, xmax, ymin, ymax,
             label, legend_title, style, large, extension, output):
    """Plot cumulative kappa against frequency or mean free path.

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

    axes = tp.axes.one_large if large else tp.axes.one
    if label is None:
        fig, ax = axes.plain(style)
    else:
        fig, ax, add_legend = axes.medium_legend(style)

    if mfp:
        data = [tp.data.load.phono3py(f, ['mode_kappa', 'mean_free_path']) for f in filenames]
        tp.plot.mfp.add_cum_kappa(ax, data, temperature=temperature,
                                  direction=direction, colour=colour,
                                  fill=fill, fillcolour=fillalpha, line=line,
                                  kmin=minkappa, scale=percent, label=label,
                                  linestyle=linestyle, marker=marker)
    else:
        data = [tp.data.load.phono3py(f, ['mode_kappa', 'frequency']) for f in filenames]
        tp.plot.frequency.add_cum_kappa(ax, data, temperature=temperature,
                                        direction=direction, colour=colour,
                                        fill=fill, fillcolour=fillalpha,
                                        line=line, scale=percent, label=label,
                                        marker=marker, linestyle=linestyle)

    if large:
        if mfp:
            ax.set_xlabel(tp.settings.large_labels()['mean_free_path'])
        if percent:
            ax.set_ylabel(tp.settings.large_labels()['cumulative_percent'])
        else:
            ax.set_ylabel(tp.settings.large_labels()['cumulative_kappa'])

    if label is not None:
        add_legend(title="${}$".format(legend_title))

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
        plt.savefig('{}.{}'.format(output, ext))

    return

@plot.command()
@input_argument
@click.option('-p', '--poscar',
              help='POSCAR path. Ignored if --atoms specified.',
              type=click.Path(file_okay=True, dir_okay=False),
              default='POSCAR',
              show_default=True)
@click.option('-a', '--atoms',
              help='Atoms in POSCAR order. Repeated names have their '
                   'contributions summed, or different names can be '
                   'used to separate environments. E.g. "Ba 1 Sn 1 O 3", '
                   '"Ba Sn O O O" and "Ba Sn O 3" are all valid and '
                   'equivalent. Overrides --poscar.')
@click.option('--projected/--notprojected',
              help='Plot atom-projected DoS.  [default: projected]',
              default=True,
              show_default=False)
@click.option('-t', '--total/--nototal',
              help='Plot total DoS.  [default: nototal]',
              default=False,
              show_default=False)
@click.option('--total-label',
              help='Label for the total line.',
              default='Total',
              show_default=True)

@click.option('-c', '--colour',
              help='Colour(s) in POSCAR order with total at the end or '
                   'colourmap name. If --notprojected, a single colour '
                   'can be specified. Total colour is overridden by '
                   '--total-colour.',
              multiple=True,
              default=['tab10'],
              show_default=True)
@click.option('--total-colour',
              help='Colour for the total line. Overrides --colour.')
@line_fill_options

@xy_limit_options
@click.option('--legend-title',
              help='Legend title. Accepts maths notation.')
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-dos',
              show_default=True)

def dos(filename, poscar, atoms, projected, total, total_label, colour,
        total_colour, fill, fillalpha, line, linestyle, marker, xmin, xmax,
        ymin, ymax, legend_title, style, large, extension, output):
    """Plot a phonon density of states."""

    axes = tp.axes.one_large if large else tp.axes.one
    fig, ax, add_legend = axes.medium_legend(style)

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
    if legend_title is not None:
        legend_title = "${}$".format(legend_title)
    add_legend(title=legend_title)

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
        plt.savefig('{}.{}'.format(output, ext))

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

@xyc_limit_options
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

    axes = tp.axes.one_large if large else tp.axes.one
    try:
        edata = tp.data.load.amset(filename)
    except UnicodeDecodeError:
        edata = tp.data.load.boltztrap(filename)

    fig, ax = axes.colourbar(style)

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
        plt.savefig('{}.{}'.format(output, ext))

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
              help='Colourmap name or colours.',
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

@legend_options
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-phonons',
              show_default=True)

def converge_phonons(filenames, bandmin, bandmax, colour, linestyle, marker,
                     xmarkcolour, xmarklinestyle, label, legend_title, style,
                     large, extension, output):
    """Plots and overlays phonon dispersions."""

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

    data = [tp.data.load.phonopy_dispersion(f) for f in filenames]

    axes = tp.axes.one_large if large else tp.axes.one
    if label != [None]:
        fig, ax, add_legend = axes.medium_legend(style)
    else:
        fig, ax = axes.plain(style)

    tp.plot.phonons.add_multi(ax, data, colour=colour, linestyle=linestyle,
                              marker=marker, label=label, bandmin=bandmin,
                              bandmax=bandmax,
                              xmarkkwargs={'color':     xmarkcolour,
                                           'linestyle': xmarklinestyle})

    if label != [None]:
        add_legend(title=legend_title)

    for ext in extension:
        plt.savefig('{}.{}'.format(output, ext))

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
@xyc_limit_options

@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-waterfall',
              show_default=True)

def waterfall(filename, y, x, projected, direction, temperature, colour, alpha,
              linewidth, edgecolour, marker, markersize, xscale, yscale,
              cscale, xmin, xmax, ymin, ymax, cmin, cmax, style, large,
              extension, output):
    """Plot 3rd-order phonon properties per band per q-point."""

    axes = tp.axes.one_large if large else tp.axes.one
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

    fig, ax = axes.plain(style)

    if projected == 'density':
        cbar = tp.plot.frequency.add_density(ax, data, y, xquantity=x,
                                             temperature=temperature,
                                             direction=direction,
                                             colour=colour, alpha=alpha,
                                             linewidth=linewidth,
                                             edgecolors=edgecolour,
                                             marker=marker, s=markersize)
    elif projected is not None:
        cbar = tp.plot.frequency.add_projected_waterfall(ax, data, y,
                  projected, xquantity=x, temperature=temperature,
                  direction=direction, colour=colour, alpha=alpha,
                  linewidth=linewidth, edgecolors=edgecolour, marker=marker,
                  s=markersize, cmin=cmin, cmax=cmax, cscale=cscale)
    else:
        tp.plot.frequency.add_waterfall(ax, data, y, xquantity=x,
                                        temperature=temperature,
                                        direction=direction, colour=colour,
                                        alpha=alpha, linewidth=linewidth,
                                        edgecolors=edgecolour, marker=marker,
                                        s=markersize)

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
        plt.savefig('{}.{}'.format(output, ext))

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

def wideband(phonons, kappa, temperature, poscar, colour, smoothing, style,
             large, extension, output):
    """Plots a broadened phonon dispersion."""

    axes = tp.axes.one_large if large else tp.axes.one
    if len(colour) == 1:
        colour = colour[0]
    else:
        colour = list(colour)

    pdata = tp.data.load.phonopy_dispersion(phonons)
    kdata = tp.data.load.phono3py(kappa, 'wideband')

    fig, ax = axes.plain(style)

    tp.plot.phonons.add_wideband(ax, kdata, pdata, temperature=temperature,
                                 poscar=poscar, smoothing=smoothing,
                                 colour=colour)

    for ext in extension:
        plt.savefig('{}.{}'.format(output, ext))

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

@xyc_limit_options
@plot_io_options
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-ztmap',
              show_default=True)

def ztmap(filename, kappa, direction, dtype, interpolate, kind, colour, xmin,
          xmax, ymin, ymax, cmin, cmax, style, large, extension, output):
    """Plots ZT against temperature and carrier concentration."""

    axes = tp.axes.one_large if large else tp.axes.one
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

    fig, ax = axes.colourbar(style)

    tp.plot.heatmap.add_ztmap(ax, edata, kdata=kdata, direction=direction,
                              xinterp=interpolate, yinterp=interpolate,
                              kind=kind, colour=colour, xmin=xmin, xmax=xmax,
                              ymin=ymin, ymax=ymax, cmin=cmin, cmax=cmax)

    for ext in extension:
        plt.savefig('{}.{}'.format(output, ext))

    return



#@plot.group(chain=True)
#def frequency():
#    """Chainable tools for plotting."""
#    return
