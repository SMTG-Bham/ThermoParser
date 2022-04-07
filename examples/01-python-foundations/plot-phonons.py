#!/usr/bin/env python3

import tp

# Axes
fig, ax, _ = tp.axes.small.one()

# Load
dispersion = tp.data.load.phonopy_dispersion('../data/zno/band.yaml')

# Plot
tp.plot.phonons.add_dispersion(ax, dispersion)

# Save
fig.savefig('phonons.png')
