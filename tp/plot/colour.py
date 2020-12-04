"""Colour scheme and colourmap generators.

Functions:
    linear:
        linear between two colours.
    uniform:
        attempt at three-colour linear for higher contrast. WIP.
    elbow:
        Like uniform, except one can chose the midpoint location.
    highlight:
        takes an existing map and highlights specific entries.
    skelton:
        rainbowy discreet colourmap.

    hsb2rgb:
        colour converter.
    rgb2array:
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

    cmin2 = rgb2array(cmin, alpha)
    cmax2 = rgb2array(cmax, alpha)
    colours = np.linspace(cmin2, cmax2, density)

    return ListedColormap(colours)

def uniform(cmid, cmin='#ffffff', cmax='#000000', alpha=1.,
          density=512):
    """Attmept at higher contrast colour maps. Requires refinement

    Adjusts mid colour position relative to similarity to endpoint
    colours.

    Arguments:
        cmid : str
            colour at midpoint.

        cmin : str, optional.
            colour at minimum. Default: #ffffff.
        cmax : str, optional
            colour at maximum. Default: #000000.
        alpha : float, optional
            colour alpha (from 0-1). Default: 1.0.

        density : int, optional
            number of colours to output. Default: 512.

    Returns:
        colormap
            colourmap.
    """

    cmin2 = np.array(rgb2array(cmin, alpha))
    cmid2 = np.array(rgb2array(cmid, alpha))
    cmax2 = np.array(rgb2array(cmax, alpha))
    cnorm = (cmid2[:3] - cmin2[:3]) / (cmax2[:3] - cmin2[:3])
    # pythagoras
    midpoint = np.sqrt(((1-cnorm[0])**2 + (1-cnorm[1])**2 + (1-cnorm[2])**2)/3)
    x = [0, midpoint, 1]
    y = [cmin2, cmid2, cmax2]
    colours = []
    x2 = np.linspace(0, 1, density)
    c = interp1d(x, y, 'linear', axis=0)
    colour = c(x2)

    return ListedColormap(colour)

def elbow(cmid, cmin='#ffffff', cmax='#000000', midpoint=0.7, alpha=1.,
          density=512):
    """Attmept at higher contrast colour maps. Requires refinement.

    Arguments:
        cmid : str
            colour at midpoint.

        cmin : str, optional.
            colour at minimum. Default: #ffffff.
        cmax : str, optional
            colour at maximum. Default: #000000.
        midpoint : float, optional
            midpoint position (from 0-1). Default: 0.7.
        alpha : float, optional
            colour alpha (from 0-1). Default: 1.0.

        density : int, optional
            number of colours to output. Default: 512.

    Returns:
        colormap
            colourmap.
    """

    x = [0, midpoint, 1]
    cmin2 = rgb2array(cmin, alpha)
    cmid2 = rgb2array(cmid, alpha)
    cmax2 = rgb2array(cmax, alpha)
    y = [cmin2, cmid2, cmax2]
    colours = []
    x2 = np.linspace(0, 1, density)
    c = interp1d(x, y, 'linear', axis=0)
    colour = c(x2)

    return ListedColormap(colour)

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
    """Converts hsb to an rgba colour array.

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
            rgba.
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

def rgb2array(colour, alpha=1.):
    """Converts #RRGGBB string to rgba array.

    Arguments:
        colour : str
            #RRGGBB colour string

        alpha : float, optional
            colour alpha (from 0-1). Default: 1.0.

    Returns:
        list
            rgba.
    """

    colour2 = [1, 1, 1, alpha]
    colour = colour.strip('#')
    for n in range(0,3):
        colour2[n] = int(colour[2*n:2*n+2], 16) / 255

    return colour2
