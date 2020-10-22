#!/usr/bin/env python3

import tp
from matplotlib import pyplot as plt

f = '../data/amset_data_85x85x47.json'
target = 2
direction = 'x'
quantities = 'seebeck conductivity ke'

colour = tp.plot.colour.elbow('#ff8080')
colour = tp.plot.colour.highlight(colour, 'grey')
interp=200

# Axes

fig, ax = tp.plot.axes.one_colourbar()

# Load

data = tp.data.load.amset(f, quantities=quantities)

# Add

cbar = tp.plot.heatmap.add_kappa_target(ax, data, zt=target, colour=colour,
                                        direction=direction, xinterp=interp,
                                        yinterp=interp)

# Save

plt.savefig('target_kl.pdf')
