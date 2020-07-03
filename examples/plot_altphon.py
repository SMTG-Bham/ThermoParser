#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = 'band.yaml'
kappafile = 'kappa-m505028.hdf5'
direction = 'norm'
temperature = 300
colour = ['#44ffff', '#ff8044', '#ff4444', '#00000010']

# Axes
fig, ax = tp.plot.axes.one()

# Load
data = tp.data.load.phono3py(kappafile, temperature=temperature,
                             quantities=['group_velocity', 'frequency',
                                         'temperature', 'qpoint'])
pdata = tp.data.load.phonopy_dispersion(phile)

# Add
ax = tp.plot.phonons.add_alt_dispersion(ax, data, pdata, 'group_velocity',
                                        direction=direction, colour=colour)

# Save
plt.savefig('altphon.pdf')
