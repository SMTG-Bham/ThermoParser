"""Functions for dealing with legends.

Functions
---------

    consolidate
        combine legends.
"""

import numpy as np

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
