"""Provides a CLI for converging phonons visually."""

import click
import matplotlib.pyplot as plt
import tp

@click.command(help='Plots phonon dispersions for comparison.')
@click.argument('files',
                help='band.yamls from Phonopy or sumo.',
                type=click.Path(exists=True, file_okay=True, dir_okay=False),
                nargs=-1)
@click.option('--bandmin',
              help='Minimum band index.')
              type=int)
@click.option('--bandmax',
              help='Maximum band index.',
              type=int)
@click.option('-c', '--colour',
              help='Colourmap name or colours.'
              multiple=True,
              default=['winter_r'],
              show_default=True)
@click.option('-e', '--extension',
              help='Output extension(s).',
              multiple=True,
              default=['pdf'],
              show_default=True)
@click.option('-l', '--labels',
              multiple=True
              help='legend labels, one for each dispersion.')
@click.option('--large/--small',
              help='Axes size.',
              default=False,
              show_default=True)
@click.option('--linestyle',
              help='Linestyle(s).'
              multiple=True,
              default=['solid'],
              show_default=True)
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-converge-phonons',
              show_default=True)
@click.option('-s', '--style',
              help='Style sheet to overlay. Later ones override '
                   'earlier ones.',
              multiple=True,
              default=[],
              show_default=False)
@click.option('-t', '--title',
              help='Legend title.')
def cli():
    """CLI for plotting phonon convergence."""

    data = [tp.data.load.phonopy_dispersion(f) for f in files]

    axes = tp.axes.one_large if large else tp.axes.one
    if labels is not None:
        fig, ax, add_legend = axes.medium_legend(style)
    else:
        fig, ax = axes.plain(style)

    tp.plot.phonons.add_multi(ax, data, colour=colour, linestyle=linestyle,
                              label=labels, bandmin=bandmin, bandmax=bandmax)

    if labels is not None:
        add_legend(title=title)

    for ext in extension:
        plt.savefig('{}.{}'.format(output, ext))

    return
