"""Colour scheme and colourmap generators.

Functions:
    linear:
        linear between two colours.
    elbow:
        attempt at bigradient linear for higher contrast. WIP.
    highlight:
        takes an existing map and highlights specific entries.
    skelton:
        rainbowy discreet colourmap. May be broken.

    hsb2rgb:
        colour converter.
"""

from matplotlib.colors import ListedColormap
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.interpolate import interp1d

def linear(cmax, cmin='#ffffff', alpha=1., density=512):
    """Generates single-gradient colour maps.

    Arguments:
        cmax : str
            colour at maximum.

        cmin : str, optional
            colour at minimum. Default: #ffffff.
        alpha : float, optional
            colour alpha (from 0-1). Default: 1.0.

        density : int
            number of colours to output. Default: 512.

    Returns:
        colourmap
            colourmap.
    """
    colours = np.full((density, 4), alpha)
    for n in range(0,3):
        cmin2 = int(cmin[2*n+1:2*n+3], 16)
        cmax2 = int(cmax[2*n+1:2*n+3], 16)
        colours[:,n] = np.linspace(cmin2/256, cmax2/256, density)

    return ListedColormap(colours)

def elbow(cmid, cmin='#ffffff', cmax='#000000',
               midpoint=0.5, alpha=1., density=512, kind='linear'):
    """Attmept at higher contrast colour maps. Requires refinement.

    Arguments:
        cmid : str
            colour at midpoint.

        cmin : str, optional.
            colour at minimum. Default: #ffffff.
        cmax : str, optional
            colour at maximum. Default: #000000.
        midpoint : float, optional
            midpoint position (from 0-1). Default: 0.5.
        alpha : float, optional
            colour alpha (from 0-1). Default: 1.0.

        density : int, optional
            number of colours to output. Default: 512.
        kind : str, optional
            interpolation scheme (see scipy.interpolate.interp1d).
            Default: linear.

    Returns:
        colormap
            colourmap.
    """

    colours = []
    for n in range(3):
        cmin2 = int(cmin[2*n+1:2*n+3], 16)
        cmid2 = int(cmid[2*n+1:2*n+3], 16)
        cmax2 = int(cmax[2*n+1:2*n+3], 16)
        x = [0, midpoint, 1]
        y = [float(cmin2)/255, float(cmid2)/255, float(cmax2)/255]
        x2 = np.linspace(0, 1, density)
        c = interp1d(x, y, kind)
        colour = c(x2)
        colours.append(list(colour))
    colours.append(list(np.full((density), alpha)))
    colours = np.array(colours).reshape(4, density)

    return ListedColormap(colours.transpose())

def highlight(cmap, colour, position=[0], density=512):
    """Highlights values in a colourmap.

    Arguments:
        cmap : colormap
            colourmap to edit
        colour : str or array-like
            colours.
        position : int or array-like, optional
            position of the colours within density. Default: 0.
        density : int, optional
            number of colours. Default: 512

    Returns:
        colormap
            highlighted colourmap.
    """

    if type(colour) == str: colour = [colour]
    if type(position) == int: position = [position]
    colours = list(cmap(np.linspace(0, 1, density)))
    for n in range(len(position)):
        colours[n] = colour[n]

    return ListedColormap(colours)


def skelton(density=512, alpha=1.):
    """Generates Jonathan Skelton's rainbowy colourmap.

    Arguments:
        density : int, optional
            number of colours. Default: 512.
        alpha : float, optional
            alpha (from 0-1). Default: 1.

    Returns:
        colourmap
            colourmap.
    """

    increment = 150. / (density - 1)
    c = [hsb2rgb((240. + i * increment) % 360., 1., 1., alpha) for i in range(density)]

    return ListedColormap(c)

def hsb2rgb(h, s, b, alpha=1):
    """Converts hsb to rgba colours.

    Arguments:
        h : float
            hue.
        s : float
            saturation.
        b : float
            brightness.

        alpha : float, optional
            colour alpha (from 0-1). Default: 1.

    Returns:
        list
            red, green, blue, alpha.
    """

    import math

    tempC = s * b
    tempMin = b - tempC

    tempHPrime = h / 60.0
    tempX = tempC * (1.0 - math.fabs((tempHPrime % 2.0) - 1.0))

    r, g, b = 0.0, 0.0, 0.0

    if tempHPrime < 1.0:
        r = tempC
        g = tempX
        b = 0
    elif tempHPrime < 2.0:
        r = tempX
        g = tempC
        b = 0
    elif tempHPrime < 3.0:
        r = 0
        g = tempC
        b = tempX
    elif tempHPrime < 4.0:
        r = 0
        g = tempX
        b = tempC
    elif tempHPrime < 5.0:
        r = tempX
        g = 0
        b = tempC
    else:
        r = tempC
        g = 0
        b = tempX

    return [r + tempMin, g + tempMin, b + tempMin, alpha]

#def watercolour(ax, data, command, temperature=300, direction='avg'):
#    """Basic string parser to generate colours for waterfall plots.
#
#    Its a mess, better use add_waterfall or add_projected_waterfall as
#    approriate.
#
#    TODO: highlight, occupation, labels.
#
#    Arguments:
#        ax : axes
#            axes to add labels to the legend of
#        data : dict
#            dictionary containing waterfall data.
#        command : str
#            what you want to do. Commands:
#
#        temperature : float, optional
#            temperature in K. Default: 300.
#        direction : str, optional
#            direction. Accepts a-c/ x-y, average/ avg or normal/ norm.
#            Default: average.
#
#    Returns:
#        axes
#            axes with new legend labels.
#        np.ndarray
#            colours.
#    """
#
#    import matplotlib as mpl
#    import tp
#    from tp.data import aniso
#
#    command = command.split()
#    ti = np.abs(np.subtract(data['temperature'], temperature)).argmin()
#    fs = np.shape(data['frequency'])
#    colours = []
#    cbars = []
#    labels = []
#
#    # lists about phon3py quantities
#    twod = ['frequency', 'gamma', 'group_velocity', 'gv_by_gv', 'heat_capacity',
#            'lifetime', 'mean_free_path', 'mode_kappa']
#    hast = ['gamma', 'heat_capacity', 'mode_kappa']
#    iso = ['group_velocity', 'gv_by_gv', 'mode_kappa']
#    log = ['gamma', 'group_velocity', 'gv_by_gv', 'heat_capacity', 'lifetime',
#           'mean_free_path', 'mode_kappa']
#
#    # lists about colours
#    generators = ['elbow', 'linear']
#
#    i = 0
#    phase = None
#    while i < len(command):
#        print('START ', command[i])
#        if command[i] in twod: # colour by 2d quantity
#            print('TWOD ', command[i])
#
#            # Sort out anisotropy etc.
#
#            quantity = np.abs(data[command[i]])
#            if command[i] in hast:
#                quantity = quantity[ti]
#            if command[i] in iso:
#                quantity = aniso.three(quantity, direction)
#
#            # rescale to log axis
#
#            if command[i] in log:
#                qsort = np.ma.masked_invalid(np.ravel(quantity)).compressed()
#                qsort = qsort[qsort.argsort()]
#                cmin = qsort[int(round(len(qsort)/100, 0))]
#                cmax = qsort[-1]
#                quantity = np.log10(quantity)
#                notinf = np.amin(np.ma.masked_invalid(quantity).compressed())
#                for j in np.argwhere(np.isinf(quantity)):
#                    quantity[j[0]][j[1]] = notinf
#                quantity = np.subtract(quantity, np.log10(cmin))
#                quantity = np.divide(quantity, np.log10(cmax))
#            else:
#                cmin = np.amin(quantity)
#                cmax = np.amax(quantity)
#                quantity = np.subtract(quantity, cmin)
#                quantity = np.divide(quantity, cmax)
#
#            i += 1
#            phase = 'cmap'
#
#        if command[i] == 'band': # colour by band
#            phase = 'cmap'
#            i += 1
#
#        if phase == 'cmap': # parsing colourmaps
#            print('CMAP ', command[i])
#            j = 0
#            try:
#                cmap = mpl.cm.get_cmap(command[i])
#            except Exception:
#                if command[i] == 'skelton':
#                    cmap = tp.plot.colour.skelton()
#                elif command[i] in generators:
#                    stop = False
#                    while not stop:
#                        if len(command) > i+j+1 and command[i+j+1].startswith('#'):
#                            print('COLOUR ', command[i + j + 1])
#                            j += 1
#                        else:
#                            stop = True
#                    if command[i] == 'linear':
#                        if j == 1:
#                            cmap = tp.plot.colour.linear(command[i+1])
#                        elif j == 2:
#                            cmap = tp.plot.colour.linear(command[i+2],
#                                                         command[i+1])
#                        else:
#                            raise Exception('Expected 1 or 2 #RGB colours '
#                                            'after linear.')
#                    elif command[i] == 'elbow':
#                        if j == 1:
#                            cmap = tp.plot.colour.elbow(command[i+1])
#                        elif j == 2:
#                            cmap = tp.plot.colour.elbow(command[i+1],
#                                                        cmax=command[i+2])
#                        elif j == 3:
#                            cmap = tp.plot.colour.elbow(command[i+2],
#                                                        cmin=command[i+1],
#                                                        cmax=command[i+3])
#                        else:
#                            raise Exception('Expected 1, 2 or 3 #RGB colours '
#                                            'after elbow.')
#                elif command[i-1] == 'band':
#                    stop = False
#                    while not stop:
#                        if len(command) > i+j and command[i+j] != 'stop':
#                            print('COLOUR ', command[i + j])
#                            j += 1
#                        else:
#                            stop = True
#                    colours == [command[i+k] for k in range(j)]
#                    if j == fs[1]:
#                        pass
#                    elif j < fs:
#                        colours == [command[i+k] for k in range(j)]
#                        while j < fs[1]:
#                            colours.append(command[i+j])
#                    elif j > fs[1]:
#                        colours = colours[:fs[1]]
#                else:
#                    raise Exception('Expected a colourmap name or colourmap '
#                                    'generator name, got {}'.format(command[i]))
#            if command[i-1] in twod:
#                colours = [list(map(cmap, q)) for q in quantity]
#                cs = np.shape(colours)
#                colours = np.reshape(colours, (cs[0]*cs[1], cs[2]))
#                cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)
#                if command[i-1] in log:
#                    cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap,
#                                                              norm=cnorm))
#                    cbar.ax.yaxis.set_major_locator(ticker.LogLocator(subs=[0,1]))
#                else:
#                    cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap,
#                                                              norm=cnorm))
#                    cbar.ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
#                    cbar.ax.yaxis.set_minor_locator(ticker.AutoMinorLoactor(2))
#            elif command[i-1] == 'band':
#                colours = [cmap(i) for i in np.linspace(0, 1, fs[1])]
#                colours = np.tile(colours, (fs[0], 1))
#            phase = None
#            i += j + 1
#
#        if phase == colour:
#
#    return ax, colours

