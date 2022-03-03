"""Functions for dealing with legends.

Functions
---------

    consolidate
        combine legends.
    add_add_legend
        add add_legend function.
"""

import numpy as np
import tp

def consolidate(axes):
    """Combine legends.

    Also removes duplicates.

    Arguments
    ---------

        axes : array-like
            axes to combine the legends for.

    Returns
    -------

        array-like
            legend handles.
        array-like
            legend labels.
    """

    axes = np.ravel(axes)
    handles = []
    labels = []

    for ax in axes:
        if ax is not None:
            h, l = ax.get_legend_handles_labels()
            for i in range(len(l)-1, -1, -1):
                if l[i] in labels:
                    del h[i]; del l[i]
            handles.extend(h)
            labels.extend(l)

    return handles, labels

def add_add_legend(ax, locations, names, defloc):
    """Adds an add_legend function.

    Arguments
    ---------

        ax : list or axes
            axes to add add_legend to.
        locations : list of dicts
            outside-axes legend arguments (loc and bbox_to_anchor) and
            axes ordinal (e.g. 'axes': 0, if there is more than one) in
            a flat array in the order of the flattened axes array by row.
        names : list of lists
            names of locations in an flat array in the order of the
            flattened axes array by row, with additional locations in
            order at the end. Ordinal integers added automatically.
        defloc : str
            default location name.

    Returns
    -------

        function
            add_legend.
    """

    ax = np.ravel(ax)
    for i in range(len(ax)):
        names[i].append(str(i+1))
    namestr = ' or '.join(['/ '.join(n) for n in names])
    exceptstr = 'location must be {}.'.format(namestr)
    for l in locations:
        if 'ncol' not in l:
            l['ncol'] = 1

    if len(ax) > 1:
        def add_legend(location=defloc, custom=False, *args, **kwargs):

            if isinstance(location, (int, float)):
                location = str(location)
            if 'ncol' in kwargs:
                for l in locations:
                    l['ncol'] = kwargs['ncol']
                    del kwargs['ncol']

            fin = False
            if custom:
                for i, a in enumerate(ax):
                    if a is not None and location in names[i]:
                        legend = a.legend(loc='best', *args, **kwargs)
                        fin = True
                        break
                if not fin:
                    for j, l in enumerate(locations):
                        if location in names[i+j+1]:
                            legend = ax[l['axes']].legend(loc=l['loc'],
                                            bbox_to_anchor=l['bbox_to_anchor'],
                                            ncol=l['ncol'], *args, **kwargs)
                            fin = True
                            break
                if not fin:
                    raise Exception(exceptstr)
            else:
                handles, labels = tp.axes.legend.consolidate(ax)
                for i, a in enumerate(ax):
                    if a is not None and location in names[i]:
                        legend = a.legend(loc='best', handles=handles,
                                          labels=labels, *args, **kwargs)
                        fin = True
                        break
                if not fin:
                    for j, l in enumerate(locations):
                        if location in names[i+j+1]:
                            legend = ax[l['axes']].legend(loc=l['loc'],
                                            bbox_to_anchor=l['bbox_to_anchor'],
                                            handles=handles, labels=labels,
                                            ncol=l['ncol'], *args, **kwargs)
                            fin = True
                            break
                if not fin:
                    raise Exception(exceptstr)

            return legend

        customstr = \
        """custom : bool, optional
                enable manual editing of handles and labels arguments.
                Default: False.
            """

    else:
        def add_legend(location=defloc, *args, **kwargs):

            if isinstance(location, (int, float)):
                location = str(location)
            if 'ncol' in kwargs:
                for l in locations:
                    l['ncol'] = kwargs['ncol']
                    del kwargs['ncol']

            fin = False
            for i, a in enumerate(ax):
                if location in names[i]:
                    legend = a.legend(loc='best', *args, **kwargs)
                    fin = True
                    break
            if not fin:
                for j, l in enumerate(locations):
                    if location in names[i+j+1]:
                        legend = ax[0].legend(loc=l['loc'],
                                        bbox_to_anchor=l['bbox_to_anchor'],
                                        ncol=l['ncol'], *args, **kwargs)
                        fin = True
                        break
            if not fin:
                raise Exception(exceptstr)

            return legend

        customstr = ''

        # docstrings don't accept str.format the normal way :(
    add_legend.__doc__ ="""Adds a pre-positioned legend.

        Accepts all normal plt.legend inputs (title etc.) except loc
        and bbox_to_anchor.

        Arguments
        ---------

            location : str, optional
                legend location. Accepts {}. Default: {}.
            {}*args, **kwargs
                passed to ax.legend.

        Returns
        -------

            legend
                legend.
        """.format(namestr, defloc, customstr)

    return add_legend
