"""Provides groups of options to streamline the CLI

Functions
---------

    direction_option
    directions_option
    doping_type_option
    dopings_option
    input_argument
    inputs_argument
    interpolate_options
    kpoints_options
    legend_options
    line_options
    line_fill_options
    plot_io_options
    temperature_option
    xy_limit_options
    xyc_limit_options
"""

import click

def direction_option(f):
    """Option for singular option for anisotropic data."""

    f = click.option('-d', '--direction',
                     help='Direction(s) for anisotropic data.',
                     type=click.Choice(['a', 'b', 'c',
                                        'x', 'y', 'z',
                                        'average', 'avg',
                                        'normal', 'norm'],
                                       case_sensitive=False),
                     default='avg',
                     show_default=True)(f)

    return f

def directions_option(f):
    """Option for multiple directions for anisotropic data."""

    f = click.option('-d', '--direction',
                     help='Direction(s) for anisotropic data.',
                     multiple=True,
                     type=click.Choice(['a', 'b', 'c',
                                        'x', 'y', 'z',
                                        'average', 'avg',
                                        'normal', 'norm'],
                                       case_sensitive=False),
                     default=['avg'],
                     show_default=True)(f)

    return f

def doping_type_option(f):
    """Option for doping type."""

    f = click.option('-t', '--type', 'dtype',
                     help='Type of doping.',
                     type=click.Choice(['n', 'p']),
                     default='n',
                     show_default=True)(f)

    return f

def dopings_option(f):
    """Option for doping concentrations."""

    f = click.option('--n', '--concentration', 'doping',
                     help='Doping concentration (will be rounded).',
                     multiple=True,
                     default=[1.e19],
                     show_default=True)(f)

    return f

def input_argument(f):
    """Option for an input file."""

    f = click.argument('filename',
                       type=click.Path(exists=True, file_okay=True,
                                       dir_okay=False))(f)

    return f

def inputs_argument(f):
    """Option for an input file."""

    f = click.argument('filenames',
                       type=click.Path(exists=True, file_okay=True,
                                       dir_okay=False),
                       nargs=-1)(f)

    return f

def interpolate_options(f):
    """Options for interpolation."""

    f = click.option('-i', '--interpolate',
                     help='Number of points to interpolate to on each axis.',
                     type=click.IntRange(1),
                     default=200,
                     show_default=True)(f)
    f = click.option('--kind',
                     help='Interpolation kind.',
                     default='linear',
                     show_default=True)(f)

    return f

def kpoints_options(f):
    """Group of options for handling KPOINTS files."""

    f = click.option('-k', '--kpoints', '-i', '--ibzkpt',
                     help='KPOINTS/IBZKPT file. Overrides --mesh.',
                     type=click.Path(exists=True,
                                     file_okay=True, dir_okay=False))(f)
    f = click.option('-m', '--mesh',
                     help='k-point mesh. Overridden by --kpoints.',
                     nargs=3,
                     type=int)(f)
    f = click.option('-p', '--poscar',
                     help='POSCAR path.',
                     type=click.Path(file_okay=True, dir_okay=False),
                     default='POSCAR',
                     show_default=True)(f)

    return f

def legend_options(f):
    """Group of options for plot legends."""

    f = click.option('-l', '--label',
                     help='Legend label(s). Accepts maths notation.',
                     multiple=True,
                     default=None)(f)
    f = click.option('--legend_title',
                     help='Legend title. Accepts maths notation.')(f)

    return f

def line_options(f):
    """Group of options for line plots"""

    f = click.option('--linestyle',
                     help='linestyle(s).',
                     multiple=True,
                     default=['solid'],
                     show_default=True)(f)
    f = click.option('-m', '--marker',
                    help='Marker(s).',
                    multiple=True,
                    default=[None])(f)

    return f

def line_fill_options(f):
    """Group of options for fillable line plots"""

    f = click.option('-f', '--fill/--nofill',
                     help='Fill under line.  [default: fill]',
                     default=True,
                     show_default=False)(f)
    f = click.option('--fillalpha',
                     help='Fill opacity (0-1). Only works if --colour is #RRGGBB.',
                     type=click.FloatRange(0, 1),
                     default=0.2,
                     show_default=True)(f)

    f = click.option('--line/--noline',
                     help='Plot line.  [default: line]',
                     default=True,
                     show_default=False)(f)
    f = click.option('--linestyle',
                     help='linestyle(s).',
                     multiple=True,
                     default=['solid'],
                     show_default=True)(f)
    f = click.option('-m', '--marker',
                    help='Marker(s).',
                    multiple=True,
                    default=[None])(f)

    return f

def plot_io_options(f):
    """Group of options for plot file I/O."""

    f = click.option('-s', '--style',
                     help='Style sheet to overlay. Later ones override '
                          'earlier ones.',
                     multiple=True,
                     default=[],
                     show_default=False)(f)
    f = click.option('--large/--small',
                     help='Axes size.  [default: small]',
                     default=False,
                     show_default=False)(f)
    f = click.option('-e', '--extension',
                     help='Output extension(s).',
                     multiple=True,
                     default=['pdf'],
                     show_default=True)(f)

    return f

def temperature_option(f):
    """Option for temperature selection."""

    f = click.option('-t', '--temperature',
                     help='Temperature (by default in K).',
                     type=click.FloatRange(0),
                     default=300.,
                     show_default=True)(f)

    return f

def xy_limit_options(f):
    """Options for x and y axes limits."""

    f = click.option('--xmin',
                     help='Override minimum x.',
                     type=float)(f)
    f = click.option('--xmax',
                     help='Override maximum x.',
                     type=float)(f)
    f = click.option('--ymin',
                     help='Override minimum y.',
                     type=float)(f)
    f = click.option('--ymax',
                     help='Override maximum y.',
                     type=float)(f)

    return f

def xyc_limit_options(f):
    """Options for x, y and colour axes limits."""

    f = click.option('--xmin',
                     help='Override minimum x.',
                     type=float)(f)
    f = click.option('--xmax',
                     help='Override maximum x.',
                     type=float)(f)
    f = click.option('--ymin',
                     help='Override minimum y.',
                     type=float)(f)
    f = click.option('--ymax',
                     help='Override maximum y.',
                     type=float)(f)
    f = click.option('--cmin',
                     help='Override minimum colour value.',
                     type=float)(f)
    f = click.option('--cmax',
                     help='Override maximum colour value.',
                     type=float)(f)

    return f
