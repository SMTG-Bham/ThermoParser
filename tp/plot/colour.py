"""Colour scheme and colourmap generators.

Functions:
    linear:
        linear between two colours.
    elbow:
        attempt at bigradient linear for higher contrast. WIP.
    highlight:
        takes an existing map and highlights specific entries.
    skelton:
        rainbowy discreet colourmap.

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
