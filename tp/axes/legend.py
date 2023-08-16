"""Functions for dealing with legends.

Functions
---------

    consolidate
        combine legends.
    add_add_legend
        add add_legend function.
    alphabetise
        alphabetises or enumerates axes.
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
            # docstring is at the bottom

            if isinstance(location, (int, float)):
                location = str(location)
            incol = 1
            if 'ncol' in kwargs:
                for l in locations:
                    l['ncol'] = kwargs['ncol']
                incol = kwargs['ncol']
                del kwargs['ncol']

            fin = False
            if custom:
                for i, a in enumerate(ax):
                    if a is not None and location in names[i]:
                        legend = a.legend(loc='best', ncol=incol,
                                          *args, **kwargs)
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
                                          labels=labels, ncol=incol,
                                          *args, **kwargs)
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
            # docstring is at the bottom

            if isinstance(location, (int, float)):
                location = str(location)
            incol = 1
            if 'ncol' in kwargs:
                for l in locations:
                    l['ncol'] = kwargs['ncol']
                incol = kwargs['ncol']
                del kwargs['ncol']

            fin = False
            for i, a in enumerate(ax):
                if location in names[i]:
                    legend = a.legend(loc='best', ncol=incol, *args, **kwargs)
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

def alphabetise(ax, labels=None, preset='latin', prefix='', suffix='',
                x=0., y=1.01, label_dos=True):
    """Enumerates or alphabetises plot axes

    Can manually define, or some presets are available.

    Arguments
    ---------

        ax : axes or list
            axes to enumerate.

        labels : str or list, optional
            manually defined labels. Overrides preset.
        preset : str
            preset label sequence. Options:

                latin
                    a, b, c... (default)
                Latin
                    A, B, C...
                arabic
                    1, 2, 3...
                roman
                    i, ii, iii...
                Roman
                    I, II, III...
                greek
                    \\alpha, \\beta, \gamma...
                Greek
                    \Alpha, \Beta, \Gamma...

        prefix : str, optional
            prefix to all labels, e.g. "(". Default: None.
        suffix : str, optional
            suffix to all labels, e.g. ")". Default: None.

        x : float, optional
            x-position, where the axis is a scale from 0-1. Default: 0.
        y : float, optional
            y-position, where the axis is a scale from 0-1. Default: 1.01.

        label_dos : bool, optional
            label DoS axes. Only works with tp axes. Default: True.
    """

    if isinstance(ax, list):
        ax = np.ravel(ax)
    else:
        ax = list(ax)

    fig = ax[0].get_figure()
    if 'dos' in fig.__dict__ and fig.__dict__['dos'] and not label_dos:
        ax = ax[:-1]
    ax = list(ax[np.where(ax != None)])

    presets = {'latin':  list('abcdefghijklmnopqrstuvwxyz'),
               'Latin':  list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'),
               'arabic': np.array(range(len(ax))) + 1,
               'roman':  'i ii iii iv v vi vii viii ix x '
                         'xi xii xiii xiv xv xvi xvii xviii xix xx'.split(),
               'Roman':  'I II III IV V VI VII VIII IX X '
                         'XI XII XIII XIV XV XVI XVII XVIII XIX XX'.split(),
               'greek':  ['$\\alpha$', '$\\beta$', '$\gamma$', '$\delta$',
                          '$\\varepsilon$', '$\zeta$', '$\eta$', '$\\theta$',
                          '$\iota$', '$\kappa$', '$\lambda$', '$\mu$',
                          '$\\nu$', '$\\xi$', 'o', '$\pi$',
                          '$\\rho$', '$\sigma$', '$\\tau$', '$\\upsilon$',
                          '$\phi$', '$\chi$', '$\psi$', '$\omega$'],
               'Greek':  ['A', 'B', '$\Gamma$', '$\Delta$',
                          'E', 'Z', 'H', '$\Theta$',
                          'I', 'K', '$\Lambda$', 'M',
                          'N', '$\Xi$', 'O', '$\Pi$',
                          'P', '$\Sigma$', 'T', '$\\Upsilon$',
                          '$\Phi$', 'X', '$\Psi$', '$\Omega$']}

    if labels is None:
        labels = presets[preset][:len(ax)]
    if len(labels) < len(ax):
        ax = ax[:len(labels)]

    for i, a in enumerate(ax):
        a.text(x, y, prefix + str(labels[i]) + suffix, ha='left',
               transform=a.transAxes)

    return
