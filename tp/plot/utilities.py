"""Utilities to aid the plotting scripts

Functions:
    colour_scale:
        Sorts colour limits and colourbar format.
"""

import matplotlib as mpl
import numpy as np

def colour_scale(c, name, cmap, cmin=None, cmax=None, cscale=None,
                 unoccupied='grey'):
    """Formats the colour scale for phono3py quantities.

    Arguments:
        c : array-like
            colour data.
        name : str
            name of colour variable.
        cmap : colormap
            colourmap.

        cmin : float, optional
            colour scale minimum. Default: display 99 % data.
        cmax : float, optional
            colour scale maximum. Default: display 99.9 % data.
        cscale : str, optional
            override colour scale (linear/ log). Default: None.
        unoccupied : str, optional
            if the colour variable is occupation, values below 1 are
            coloured in this colour. If set to None, or cmin is set,
            this feature is turned off. Default: grey.

    Returns:
        cnorm
            colour normalisation object.
        extend
            which direction to extend the colourbar in.
    """

    extend = [False, False]
    if name in ['frequency', 'heat_capacity']:
        if cmin is None:
            cmin = np.amin(c)
        elif cmin > np.amin(c):
            extend[0] = True
        if cmax is None:
            cmax = np.amax(c)
        elif cmax < np.amax(c):
            extend[1] = True
        cnorm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)

    elif name == 'occupation':
        csort = np.ravel(np.ma.masked_invalid(np.ma.masked_equal(c, 0)).compressed())
        csort = csort[csort.argsort()]
        clim = csort[int(round(len(csort)*99.9/100 - 1, 0))]
        if cmin is None:
            cmin = np.amin(c)
            if unoccupied is not None and cmin < 1:
                cmin = 1
                extend[0] = True
                cmap.set_under(unoccupied)
        elif cmin > np.amin(c):
            extend[0] = True
        if cmax is None:
            cmax = clim
            extend[1] = True
        elif cmax < np.amax(c):
            extend[1] = True
        if cmax/cmin < 10:
            cnorm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
        else:
            cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    else:
        csort = np.ravel(np.ma.masked_invalid(np.ma.masked_equal(c, 0)).compressed())
        csort = csort[csort.argsort()]
        clim = [csort[int(round(len(csort)/100 - 1, 0))],
                csort[int(round(len(csort)*99.9/100 - 1, 0))]]
        if cmin is None:
            cmin = clim[0]
        elif cmin > clim[0]:
            extend[0] = True
        if cmax is None:
            cmax = clim[1]
        elif cmax < clim[1]:
            extend[1] = True
        cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    if cscale is not None:
        if cscale == 'linear':
            cnorm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
        elif cscale == 'log':
            cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)

    if extend[0] and extend[1]:
        extend = 'both'
    elif extend[0] and not extend[1]:
        extend = 'min'
    elif not extend[0] and extend[1]:
        extend = 'max'
    elif not extend[0] and not extend[1]:
        extend = 'neither'

    return cnorm, extend
