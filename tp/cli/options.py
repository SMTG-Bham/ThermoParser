"""Provides groups of options to streamline the CLI

Functions
---------

    direction_function:
        function for picking the --direction (-d).
    doping_type_option:
        function for picking the doping --type (-t).
    doping_function:
        function for picking the doping --concentration (-n).
    dos_function:
        function for setting DoS formatting.
    inputs_function:
        function for picking the input file(s) argument.
    interpolate_options:
        function for setting interpolation options.
    kpoints_options:
        function for handling KPOINTS files.
    legend_options:
        function for formatting legends.
    line_options:
        function for formatting line plots.
    fill_options:
        function for formatting fillable line plots.
    plot_io_function:
        function for formatting plot outputs.
    temperature_option:
        function for picking the --temperature (-t).
    verbose_option:
        function for increasing the verbosity.
    axes_limit_function
        function for setting the axis limits.
"""

import click

def direction_function(multiple=False):
    """Function to create direction options.
        
    Arguments
    ---------

        multiple : bool, optional
            whether to allow multiple directions. Default: False.

    Returns
    -------
        decorator
            direction option decorator.
    """

    default = ['avg'] if multiple else 'avg'
    def direction_option(f):
        """Option for anisotropic data.
        
        Options
        -------

            --direction, -d : str or array-like, optional
                direction. Default: avg.

        Returns
        -------

            decorator
                direction option decorator.
        """
    
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
    """Option for doping type.
        
        Options
        -------

            --type : str, optional
                doping type. Default: n.

        Returns
        -------

            decorator
                doping type option decorator.
        """

    f = click.option('--type', 'dtype',
                     help='Type of doping.',
                     type=click.Choice(['n', 'p']),
                     default='n',
                     show_default=True)(f)

    return f

def doping_function(multiple=False):
    """Function to create doping options.
        
    Arguments
    ---------

        multiple : bool, optional
            whether to allow multiple concentrations. Default: False.

    Returns
    -------
        decorator
            carrier concentration option decorator.
    """

    default = [1e19] if multiple else 1e19
    def doping_option(f):
        """Option for doping concentration.
        
        Options
        -------

            --concentration, -n : int or float or array-like, optional
                carrier concentration. Default: 1e19.

        Returns
        -------

            decorator
                carrier concentration option decorator.
        """
    
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
    """Function for creating DoS options.

    Arguments
    ---------

        dosargs : array-like, optional
            names for the colour argument. Default: ['-c', '--colour'].

    Returns
    -------
        decorator
            DoS options decorator.
    """

    if isinstance(dosargs, str):
        dosargs = [dosargs]
    def dos_options(f):
        """Options for DoS plots.
        
        Options
        -------

            --poscar, -p : str, optional
                POSCAR path. Ignored if atoms specified. Default: POSCAR.
            --atoms : str, optional
                atoms in POSCAR order. Overrides poscar. Default: None.
            --sigma : float, optional
                standard deviation for gaussian broadening.
                Recommended: 0.2. Default: None.
            --projected/--notprojected : bool, optional
                plot atom-projected DoS. Default: --projected.
            --colour : str or array-like, optional
                colours in atom order followed by total or colourmap.
                Default: tab10.
            --total, -t/--nototal : bool, optional
                plot the total DoS. Default: --nototal.
            --total_label : str, optional
                total label for the legend. Default: total.
            --total_colour : str, optional
                total line colour. Overrides --colour.

        Returns
        -------

            decorator
                DoS options decorator.
        """
    
        f = click.option('-p', '--poscar',
                         help='POSCAR path. Ignored if --atoms specified.',
                         type=click.Path(file_okay=True, dir_okay=False,
                                         exists=False),
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
    """Function for creating input arguments.

    Arguments
    ---------

        name : str, optional
            filename argument. Default: 'filenames'.
        nargs : int, optional
            number of input files allowed. Default: -1.

    Returns
    -------
        decorator
            input argument decorator.
    """

    def inputs_argument(f):
        """Option for input files.
        
        Options
        -------

            name : str or array-like
                input files.

        Returns
        -------

            decorator
                input argument decorator.
        """

        f = click.argument(name,
                           type=click.Path(exists=True, file_okay=True,
                                           dir_okay=False),
                           nargs=nargs)(f)

        return f
    return inputs_argument

def interpolate_options(f):
    """Options for interpolation.
        
        Options
        -------

            --interpolate, -i : int, optional
                number of points to interpolate. Default: 200.
            --kind : str, optional
                interpolation kind. Default: linear.

        Returns
        -------

            decorator
                interpolation options decorator.
        """

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
    """Group of options for handling KPOINTS files.
    
        Options
        -------

            --kpoints, -k, --ibzkbt, -i : str, optional
                KPOINTS/IBZKPT file. Overrides --mesh.
            --mesh, -m : 3x int, optional
                k-point mesh. Overridden by --kpoints.
            --poscar, -p : str, optional
                POSCAR. Required for --mesh.

        Returns
        -------

            decorator
                kpoints options decorator.
        """

    f = click.option('-k', '--kpoints', '-i', '--ibzkpt',
                     help='KPOINTS/IBZKPT file. Overrides --mesh.',
                     type=click.Path(exists=True,
                                     file_okay=True, dir_okay=False))(f)
    f = click.option('-m', '--mesh',
                     help='k-point mesh. Overridden by --kpoints.',
                     nargs=3,
                     type=int)(f)
    f = click.option('-p', '--poscar',
                     help='POSCAR path. Required for --mesh.',
                     type=click.Path(file_okay=True, dir_okay=False),
                     default='POSCAR',
                     show_default=True)(f)

    return f

def legend_function(toggle=True, label=True):
    """Function for creating legend options.

    Arguments
    ---------

        toggle : bool, optional
            include option to remove legend. Default: True.
        label : bool, optional
            include option to set labels manually. Default: True.

    Returns
    -------
        decorator
            legend options decorator.
    """

    def legend_options(f):
        """Group of options for plot legends.
    
        Options
        -------

            --label, -l : str, optional
                legend label(s).
            --legend_title : str, optional
                legend title.
            --legend/--nolegend : bool, optional
                show legend. Default: --legend.
            --location : str or int, optional
                legend location.

        Returns
        -------

            decorator
                legend options decorator.
        """
    
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
    """Group of options for line plots
    
        Options
        -------

            --linestyle : str, optional
                linestyle(s). Default: solid.
            --marker : str, optional
                marker(s). Default: None.

        Returns
        -------

            decorator
                line options decorator.
        """

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
    """Group of options for fillable line plots
    
        Options
        -------

            --fill, -f/--nofill : bool, optional
                fill under line. Default: --fill.
            --fillalpha : str, optional
                fill opacity (0-1). Only works for #RRGGBB colours.
                Default: 0.2.
            --line/--noline : bool, optional
                show line. Default: --line.

        Returns
        -------

            decorator
                fill options decorator.
        """

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
    """Function for creating plot I/O options.

    Arguments
    ---------

        name : str
            output filename.

    Returns
    -------
        decorator
            plot I/O options decorator.
    """

    def plot_io_options(f):
        """Group of options for plot file I/O.
    
        Options
        -------

            --style, -s : str, optional
                style sheet to overlay. Later ones override earlier
                ones.
            --large/--small : str, optional
                axes size. Default: small.
            --save/--nosave : bool, optional
                write to file. Default: --save.
            --show/--noshow : bool, optional
                show plot. Default: --noshow.
            --extension : str, optional
                output extensions(s). Default: pdf.
            --output : str, optional.
                output filenames, sans extension.

        Returns
        -------

            decorator
                plot I/O options decorator.
        """
    
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
    """Option for temperature selection.
        
        Options
        -------

            --temperature, -t : str, optional
                temperature. Default: 300.

        Returns
        -------

            decorator
                temperature option decorator.
        """

    f = click.option('-t', '--temperature',
                     help='Temperature (by default in K).',
                     type=click.FloatRange(0),
                     default=300.,
                     show_default=True)(f)

    return f

def verbose_option(f):
    """Option for output verbosity.
        
        Options
        -------

            --verbose/--notverbose : bool, optional
                output plot conditions. Default: --verbose.

        Returns
        -------

            decorator
                verbosity option decorator.
        """

    f = click.option('--verbose/--notverbose',
                     help='Output plot conditions.  [default: verbose]',
                     default=True,
                     show_default=False)(f)

    return f

def axes_limit_function(multiple=False, c=False):
    """Function for creating axes limit options.

    Arguments
    ---------

        multiple : bool, optional
            allow multiple limits. Default: False.
        c : bool, optional
            include colour limits. Default: False.

    Returns
    -------
        decorator
            axes limit options decorator.
    """

    def axes_limit_options(f):
        """Options for axes limits.
        
        Options
        -------

            --xmin : float, optional
                override minimum x.
            --xmax : float, optional
                override maximum x.
            --ymin : float, optional
                override minimum y.
            --ymax : float, optional
                override maximum y.
            --cmin : float, optional
                override minimum colour.
            --cmax : float, optional
                override maximum colour.

        Returns
        -------

            decorator
                axes limits options decorator.
        """
    
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
