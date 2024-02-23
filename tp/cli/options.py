"""Provides groups of options to streamline the CLI

Functions
---------

    direction_function
    
    doping_type_option
    
    doping_function
    
    dos_function
    
    heatmap_options
    
    input_argument
    
    inputs_function
    
    interpolate_options
    
    kpoints_options
    
    legend_options
    
    line_options
    
    fill_options
    
    plot_io_function
    
    temperature_option
    
    verbose_option

    axes_limit_function
"""

import click

def direction_function(multiple=False):
    default = ['avg'] if multiple else 'avg'
    def direction_option(f):
        """Option for anisotropic data."""
    
        f = click.option('-d', '--direction',
                         help='Direction(s) for anisotropic data.',
                         type=click.Choice(['a', 'b', 'c',
                                            'x', 'y', 'z',
                                            'average', 'avg',
                                            'normal', 'norm'],
                                           case_sensitive=False),
                         multiple=multiple,
                         default=default,
                         show_default=True)(f)
    
        return f
    return direction_option

def doping_type_option(f):
    """Option for doping type."""

    f = click.option('-t', '--type', 'dtype',
                     help='Type of doping.',
                     type=click.Choice(['n', 'p']),
                     default='n',
                     show_default=True)(f)

    return f

def doping_function(multiple=False):
    default = [1e19] if multiple else 1e19
    def doping_option(f):
        """Option for doping concentration."""
    
        f = click.option('-n', '--concentration', 'doping',
                         help='Doping concentration(s) (will be rounded).',
                         multiple=multiple,
                         default=default,
                         type=float,
                         show_default=True)(f)
    
        return f
    return doping_option

def heatmap_options(f):
    """Options for heatmaps."""

    f = click.option('--discrete/--continuous',
                     help='Discretise colourmap.  [default: continuous]',
                     default=False,
                     show_default=False)(f)
    f = click.option('-l', '--levels',
                     help='Levels for discrete plots. Lists directly '
                          'specify the contour levels, while integers '
                          'specify the maximum-1 number of "nice" '
                          'levels to plot.',
                     multiple=True,
                     default=None)(f)
    f = click.option('--contours',
                     help='Contours to plot.',
                     multiple=True,
                     type=float,
                     default=None)(f)
    f = click.option('--contourcolours',
                     help='contour colours',
                     multiple=True,
                     default='black')(f)
    return f

def dos_function(dosargs=['-c', '--colour']):
    if isinstance(dosargs, str):
        dosargs = [dosargs]
    def dos_options(f):
        """Options for DoS plots."""
    
        f = click.option('-p', '--poscar',
                         help='POSCAR path. Ignored if --atoms specified.',
                         type=click.Path(file_okay=True, dir_okay=False),
                         default='POSCAR',
                         show_default=True)(f)
        f = click.option('--atoms',
                         help='Atoms in POSCAR order. Repeated names have '
                              'their contributions summed, or different names '
                              'can be used to separate environments. E.g. '
                              '"Ba 1 Sn 1 O 3", "Ba Sn O O O" and "Ba Sn O 3" '
                              'are all valid and equivalent. '
                              'Overrides --poscar.')(f)
        f = click.option('--sigma',
                         help='Standard deviation of Gaussian broadening. 0.2 '
                              'is a good place to start. Does not know if '
                              'you\'ve already broadened it. Off by default.',
                        type=float,
                        default=None)(f)
        f = click.option('--projected/--notprojected',
                         help='Plot atom-projected DoS.  [default: projected]',
                         default=True,
                         show_default=False)(f)
        f = click.option(*dosargs,
                         help='Colour(s) in POSCAR order with total at the '
                              'end or colourmap name. If --notprojected, a '
                              'single colour can be specified. Total colour '
                              'is overridden by --total-colour.',
                         multiple=True,
                         default=['tab10'],
                         show_default=True)(f)
        f = click.option('-t', '--total/--nototal',
                         help='Plot total DoS.  [default: nototal]',
                         default=False,
                         show_default=False)(f)
        f = click.option('--total-label',
                         help='Label for the total line.',
                         default='Total',
                         show_default=True)(f)
        f = click.option('--total-colour',
                         help='Colour for the total line. Overrides --colour.')(f)
    
        return f
    return dos_options

def inputs_function(name='filenames', nargs=-1):
    def inputs_argument(f):
        """Option for input files."""
        f = click.argument(name,
                           type=click.Path(exists=True, file_okay=True,
                                           dir_okay=False),
                           nargs=nargs)(f)

        return f
    return inputs_argument

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

def legend_function(toggle=True, label=True):
    def legend_options(f):
        """Group of options for plot legends."""
    
        if label:
            f = click.option('-l', '--label',
                             help='Legend label(s).',
                             multiple=True,
                             type=str,
                             default=None)(f)
        f = click.option('--legend_title',
                         help='Legend title.')(f)
        if toggle:
            f = click.option('--legend/--nolegend',
                             help='Show legend.  [default: legend]',
                             default=True,
                             show_default=False)(f)
        f = click.option('--location',
                         help='Legend location. Accepts strings such as '
                              'above and below, as well as 1-indexed '
                              'ordinals corresponding to inside those axes',
                         type=str)(f)
    
        return f
    return legend_options

def line_options(f):
    """Group of options for line plots"""

    f = click.option('--linestyle',
                     help='linestyle(s).',
                     multiple=True,
                     default=['solid'],
                     type=str,
                     show_default=True)(f)
    f = click.option('-m', '--marker',
                    help='Marker(s).',
                    multiple=True,
                    default=[None],
                    type=str)(f)

    return f

def fill_options(f):
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

    return f

def plot_io_function(name):
    def plot_io_options(f):
        """Group of options for plot file I/O."""
    
        f = click.option('-s', '--style',
                         help='Style sheet to overlay. Later ones override '
                              'earlier ones.',
                         multiple=True,
                         default=[],
                         type=str,
                         show_default=False)(f)
        f = click.option('--large/--small',
                         help='Axes size.  [default: small]',
                         default=False,
                         show_default=False)(f)
        f = click.option('--save/--nosave',
                         help='Write to file.  [default: save]',
                         default=True,
                         show_default=False)(f)
        f = click.option('--show/--noshow',
                         help='Show plot.  [default: noshow]',
                         default=False,
                         show_default=False)(f)
        f = click.option('--extension',
                         help='Output extension(s). Requires --save.',
                         multiple=True,
                         default=['pdf'],
                         type=str,
                         show_default=True)(f)
        f = click.option('-o', '--output',
                         help='Output filename, sans extension.',
                         default=name,
                         show_default=True)(f)

        return f
    return plot_io_options

def temperature_option(f):
    """Option for temperature selection."""

    f = click.option('-t', '--temperature',
                     help='Temperature (by default in K).',
                     type=click.FloatRange(0),
                     default=300.,
                     show_default=True)(f)

    return f

def verbose_option(f):
    """Option for verbose output."""

    f = click.option('--verbose/--notverbose',
                     help='Output plot conditions.  [default: verbose]',
                     default=True,
                     show_default=False)(f)

    return f

def axes_limit_function(multiple=False, c=False):
    def axes_limit_options(f):
        """Options for axes limits."""
    
        f = click.option('--xmin',
                         help='Override minimum x.',
                         multiple=multiple,
                         type=float)(f)
        f = click.option('--xmax',
                         help='Override maximum x.',
                         multiple=multiple,
                         type=float)(f)
        f = click.option('--ymin',
                         help='Override minimum y.',
                         multiple=multiple,
                         type=float)(f)
        f = click.option('--ymax',
                         help='Override maximum y.',
                         multiple=multiple,
                         type=float)(f)
        if c:
            f = click.option('--cmin',
                             help='Override minimum colour value.',
                             multiple=multiple,
                             type=float)(f)
            f = click.option('--cmax',
                             help='Override maximum colour value.',
                             multiple=multiple,
                             type=float)(f)
        return f
    return axes_limit_options
