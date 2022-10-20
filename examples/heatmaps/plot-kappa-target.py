#!/usr/bin/env python3

import tp

f = '../data/zno/boltztrap.hdf5'
target = 2
direction = 'x'

colour = tp.plot.colour.uniform('#008080')

# Axes
fig, ax, _ = tp.axes.small.one_colourbar()

# Load
data = tp.data.load.boltztrap(f)

# Add
tp.plot.heatmap.add_kappa_target(ax, data, zt=target, colour=colour,
                                 direction=direction)

# Save
fig.savefig('kappa-target.pdf')
fig.savefig('kappa-target.png')
