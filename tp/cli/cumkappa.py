"""Provides a CLI for plotting cumulative kappa"""

import click
import matplotlib.pyplot as plt
import tp

@click.command(help='Plot cumulative kappa against frequency or mean '
                    'free path. Properties such as colour and '
                    'linestyle loop, so if you have two data files and '
                    'two directions, only two colours need to be '
                    'specified, one for each direction; however if you '
                    'want one for each datafile, each will need to be '
                    'repeated twice consecutively.')
@click.argument('file',
                help='Phono3py kappa-mxxx.hdf5(s).')
                nargs=-1)
@click.option('-c', '--colour',
              help='Colour(s).',
              multiple=True)
@click.option('-d', '--direction',
              help='Direction(s) for anisotropic data.',
              multiple=True,
              type=click.Choice(['a','b','c', 'x','y','z', 'average','avg'],
                                case_sensitive=False)
              default=['avg'],
              show_default=True)
@click.option('-e', '--extension',
              help='Output extension(s).',
              multiple=True,
              default=['pdf'],
              show_default=True)
@click.option('-f', '--fill/--nofill',
              help='Fill under line',
              default=True,
              show_default=True)
@click.option('--fillalpha',
              help='Fill opacity (0-1). Only works if --colour is #RRGGBB.')
              default=0.2,
              show_default=True)
@click.option('-l', '--label',
              help='Legend label(s). Accepts maths notation.')
@click.option('--large/--small',
              help='Axes size.',
              default=True,
              show_default=True)
@click.option('--legend_title',
              help='Legend title. Accepts maths.')
@click.option('--linestyle',
              help='linestyle(s).',
              multiple=Truem
              default=['solid'],
              show_default=True)
@click.option('--mfp/--frequency',
              help='x-axis quantity.',
              default=True,
              show_default=True)
@click.option('-m', '--marker',
              help='Marker(s).',
              multiple=True)
@click.option('--minkappa',
              help='Minimum kappa to plot for --mfp in percent.',
              default=1.,
              show_default=True)
@click.option('--line/--noline',
              help='Plot line.',
              default=True,
              show_default=True)
@click.option('-o', '--output',
              help='Output filename, sans extension.',
              default='tp-cumkappa',
              show_default=True)
@click.option('-p', '--percent/--raw',
              help='Plot kappa in percent of total.',
              default=False,
              show_default=True)
@click.option('-s', '--style',
              help='Style sheet to overlay. Later ones override '
                   'earlier ones.',
              multiple=True,
              default=[],
              show_default=False)
@click.option('-t', '--temperature',
              help='Temperature in K.',
              default=300.,
              show_default=True)
@click.option('--xmin', metavar='xmin',
              help='Override minimum x.',
              type=float)
def cli():
    """CLI for plotting cumulative kappa."""
    axes = tp.axes.one_large if large else tp.axes.one
    if label is None:
        fig, ax = axes.plain(style)
    else:
        fig, ax, add_legend = axes.medium_legend(style)

    if mfp:
        data = [tp.data.load.phono3py(f, ['mode_kappa', 'mean_free_path']) for f in file]
        tp.plot.mfp.add_cum_kappa(ax, data, temperature=temperature,
                                  direction=direction, colour=colour,
                                  fill=fill, fillcolour=fillalpha, line=line,
                                  kmin=minkappa, scale=percent, label=label,
                                  linestyle=linestyle, marker=marker)
    else:
        data = [tp.data.load.phono3py(f, ['mode_kappa', 'frequency']) for f in file]
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

    for ext in extension:
        plt.savefig('{}.{}'.format(output, ext))

    return
