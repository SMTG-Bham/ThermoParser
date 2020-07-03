#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = 'band.yaml'
kappafile = 'kappa-m505028.hdf5'
temperature = 300
colour = tp.plot.colour.linear('#ff0000', '#000000')

# Axes
fig, ax = tp.plot.axes.one(['pretty2', 'dark_background'])

# Load
data = tp.data.load.phono3py(kappafile, quantities=['gamma', 'frequency',
                                                    'temperature', 'qpoint'])
pdata = tp.data.load.phonopy_dispersion(phile)

# Add
ax = tp.plot.phonons.add_wideband(ax, data, pdata, temperature=temperature,
                                  colour=colour)

# Save
plt.savefig('wideband.pdf')
