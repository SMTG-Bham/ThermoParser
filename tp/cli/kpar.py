"""Provides CLI for generating KPAR."""

import click
import tp

@click.command(help='Suggests KPAR values for VASP. Reads a KPOINTS '
                    'file by default, but a different file or mesh can '
                    'be specified.')
@click.option('-k', '--kpoints', '-i', '--ibzkpt',
              help='KPOINTS/IBZKPT file. Overrides --mesh.',
              type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-m', '--mesh',
              help='k-point mesh. Overrides --kpoints.',
              nargs=3,
              type=int)
@click.option('-p', '--poscar',
              help='POSCAR path.',
              default='POSCAR',
              show_default=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False))
def cli():
    """CLI for generating KPAR."""

    if kpoints is not None:
        pass
    elif mesh is not None:
        kpoints = mesh
    else:
        kpoints = 'KPOINTS'

    k = tp.setup.vasp.get_kpar(kpoints, poscar=poscar)
    ks = [str(kp) for kp in k]
    click.echo('KPAR = {} or {}'.format(','.join(ks[:-1]), ks[-1]))

    return
