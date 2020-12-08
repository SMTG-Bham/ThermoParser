#!/usr/bin/env python3

import tp
from matplotlib import pyplot as plt

f = 'data/zno/boltztrap.hdf5'
target = 2
direction = 'x'

colour = tp.plot.colour.uniform('#008080')
colour = tp.plot.colour.highlight(colour, 'grey')

# Axes

fig, ax = tp.axes.one.colourbar()

# Load

data = tp.data.load.boltztrap(f)

# Add

tp.plot.heatmap.add_kappa_target(ax, data, zt=target, colour=colour,
                                 direction=direction)

# Save

plt.savefig('target-kl.pdf')
