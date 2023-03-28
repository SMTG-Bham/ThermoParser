"""Colour scheme and colourmap generators.

Functions
---------

    linear:
        linear between two colours.
    uniform:
        bigradient colourmap for higher contrast.
    elbow:
        like uniform, except one can chose the midpoint location.
    highlight:
        takes an existing map and highlights specific entries.
    skelton:
        rainbowy discreet colourmap.


    hsb2rgb:
        colour converter.
"""

from matplotlib.colors import ListedColormap, to_rgba
import numpy as np
from scipy.interpolate import interp1d

def linear(cmax, cmin='white', alpha=1., density=512):
    """Generates single-gradient colour maps.

    Accepts named, hex and RGB (values 0-1) formats.

    Arguments
    ---------

        cmax : str
            colour at maximum.

        cmin : str, optional
            colour at minimum. Default: white.
        alpha : float, optional
            colour alpha (from 0-1). Default: 1.0.

        density : int
            number of colours to output. Default: 512.

    Returns
    -------

        colourmap
            colourmap.
    """
    cmax = to_rgba(cmax, alpha)
    cmin = to_rgba(cmin, alpha)
    colours = np.abs(np.linspace(cmin, cmax, density))

    return ListedColormap(colours)

def uniform(cmid, cmin='white', cmax='#333333', alpha=1.,
            density=512):
    """Generates bigradient colourmaps.

    Adjusts mid colour position to keep the overall gradient even.
    Accepts named, hex and RGB (values 0-1) formats.

    Arguments
    ---------

        cmid : str
            colour at midpoint.

        cmin : str, optional.
            colour at minimum. Default: white.
        cmax : str, optional
            colour at maximum. Default: #333333.
        alpha : float, optional
            colour alpha (from 0-1). Default: 1.0.

        density : int, optional
            number of colours to output. Default: 512.

    Returns
    -------

        colormap
            colourmap.
    """

    cmax = np.array(to_rgba(cmax, alpha))
    cmin = np.array(to_rgba(cmin, alpha))
    cmid = np.array(to_rgba(cmid, alpha))
    cnorm = (cmid[:3] - cmin[:3]) / (cmax[:3] - cmin[:3])
    # pythagoras
    midpoint = np.sqrt(cnorm[0]**2 + cnorm[1]**2 + cnorm[2]**2)/3
    x = [0, midpoint, 1]
    y = [cmin, cmid, cmax]
    x2 = np.linspace(0, 1, density)
    c = interp1d(x, y, 'linear', axis=0)
    colour = np.abs(c(x2))

    return ListedColormap(colour)

def elbow(cmid, cmin='white', cmax='black', midpoint=0.7, alpha=1.,
          density=512):
    """Generates bigradient colourmaps.

    Allows for full customisation of colours and midpoint location.
    Accepts named, hex and RGB (values 0-1) formats.

    Arguments
    ---------

        cmid : str
            colour at midpoint.

        cmin : str, optional.
            colour at minimum. Default: white.
        cmax : str, optional
            colour at maximum. Default: black.
        midpoint : float, optional
            midpoint position (from 0-1). Default: 0.7.
        alpha : float, optional
            colour alpha (from 0-1). Default: 1.0.

        density : int, optional
            number of colours to output. Default: 512.

    Returns
    -------

        colormap
            colourmap.
    """

    x = [0, midpoint, 1]
    cmin = to_rgba(cmin, alpha)
    cmid = to_rgba(cmid, alpha)
    cmax = to_rgba(cmax, alpha)
    y = [cmin, cmid, cmax]
    x2 = np.linspace(0, 1, density)
    c = interp1d(x, y, 'linear', axis=0)
    colour = np.abs(c(x2))

    return ListedColormap(colour)

def highlight(cmap, colour, position=[0], density=512):
    """Highlights values in a colourmap.

    Arguments
    ---------

        cmap : colormap
            colourmap to edit
        colour : str or array-like
            colours.
        position : int or array-like, optional
            position of the colours within density. Default: 0.
        density : int, optional
            number of colours. Default: 512

    Returns
    -------

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

    Arguments
    ---------

        density : int, optional
            number of colours. Default: 512.
        alpha : float, optional
            alpha (from 0-1). Default: 1.

    Returns
    -------

        colourmap
            colourmap.
    """

    increment = 150. / (density - 1)
    c = [hsb2rgb((240. + i * increment) % 360., 1., 1., alpha) for i in range(density)]

    return ListedColormap(c)

def hsb2rgb(h, s, b, alpha=1):
    """Converts hsb to an rgba colour array.

    Arguments
    ---------

        h : float
            hue.
        s : float
            saturation.
        b : float
            brightness.

        alpha : float, optional
            colour alpha (from 0-1). Default: 1.

    Returns
    -------

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
