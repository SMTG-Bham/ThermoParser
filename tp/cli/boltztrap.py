"""Provides a command line BoltzTraP runner."""

import click
import numpy as np
from tp.data.run import boltztrap

@click.command(help='Runs boltzTraP and outputs to hdf5.')
@click.option('--dmin',
              help='Minimum doping in cm^-3.',
              default=1e18,
              show_default=True)
@click.option('--dmax',
              help='Maximum doping in cm^-3.',
              default=1e21,
              show_default=True)
@click.option('-k', '--kmode',
              help='Kappa calculation method. Options:'
                   'boltzmann: boltztrap method (default); '
                   'wiedemann: Wiedemann-Franz law with constant L; '
                   'snyder:    Wiedemann-Franz law with L dependant on '
                              'Seebeck coefficient.')
              default='boltzmann',
              type=click.Choice(['boltzmann', 'wiedemann', 'snyder'],
                                case_sensitive=False),
              show_default=False)
@click.option('--kpoints',
              help='Path to KPOINTS if there are zero-weighted kpoints.',
              default='KPOINTS',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              show_default=True)
@click.option('-l', '--lpfac',
              help='Factor to interpolate the DoS density by.',
              default=10,
              show_default=True)
@click.option('--ndope',
              help='number of doping concetrations to calculate.',
              default=16,
              show_default=True
@click.option('--analyse/--no_analyse',
              help='Analyse BoltzTraP.',
              default=True,
              show_default=True)
@click.option('--run/--no_run',
              help='Run BoltzTraP.',
              default=True,
              show_default=True)
@click.option('-o', '--output',
              help='Output filename.',
              default='boltztrap.hdf5',
              show_default=True)
@click.option('-r', '--relaxation_time',
              help='Charge carrier relaxation time in s.',
              default=1e-14,
              show_default=True)
@click.option('--tmin',
              help='Set minimum temperature in K. Does not speed up '
                   'calculation, just discards data.',
              type=float)
@click.option('--tmax',
              help='Maximum temperature in K.')
              default=1000.,
              show_default=True)
@click.option('--tstep',
              help='Temperature step size in K.',
              default=50.,
              show_default=True)
@click.option('-v', '--vasprun',
              help='Path to vasprun.xml.',
              default='vasprun.xml',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              show_default=True)
def cli():
    """CLI to run BoltzTraP."""

    doping = np.geomspace(dmin, dmax, ndope)

    boltztrap(tmax=tmax, tstep=tstep, tmin=tmin, doping=doping, ke_mode=kmode,
              vasprun=vasprun, kpoints=kpoints, relaxation_time=relaxation_time,
              lpfac=lpfac, run=run, analyse=analyse, output=output)

    return
