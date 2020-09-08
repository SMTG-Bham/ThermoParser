#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = '../data/band.yaml'
kappafile = '../data/kappa-m505028.hdf5'
poscar = '../data/POSCAR'
temperature = 300

colour = tp.plot.colour.linear('#ff0000', '#000000')

# Axes

fig, ax = tp.plot.axes.one(['tp', 'dark_background'])

# Load

data = tp.data.load.phono3py(kappafile, quantities='wideband')
pdata = tp.data.load.phonopy_dispersion(phile)

# Add

ax = tp.plot.phonons.add_wideband(ax, data, pdata, temperature=temperature,
                                  colour=colour, poscar=poscar)

# Save

plt.savefig('wideband.pdf')
