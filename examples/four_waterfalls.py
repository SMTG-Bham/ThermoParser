#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import tp

kappafile = 'kappa-m505028.hdf5'
direction = [['avg', 'norm'], ['norm', 'norm']]
temperature = 300
quantities = 'mode_kappa mfp lifetime gv waterfall'

# Axes
fig, ax = tp.plot.axes.four_square(['tp', 'dark_background'])

# Load
data = tp.data.load.phono3py(kappafile, quantities=quantities)

colour = [['', ''], ['', '']]
colour[0][0] = '#800080'
colour[0][1] = 'viridis'
colour[1][1] = tp.plot.colour.linear('#FF00FF', '#00FF00')
colour[1][0] = np.repeat('#555555', 15)
colour[1][0][8] = 'red'
plots = [['mode_kappa', 'group_velocity'],
         ['lifetime', 'mean_free_path']]

# Add
for (i,j), a in np.ndenumerate(ax):
    ax[i][j] = tp.plot.frequency.add_waterfall(a, data, plots[i][j],
                                               main=True, colour=colour[i][j],
                                               temperature=temperature,
                                               direction=direction[i][j])

# Save
plt.savefig('four-waterfall.pdf')
