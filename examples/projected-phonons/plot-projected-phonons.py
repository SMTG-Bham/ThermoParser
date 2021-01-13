#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = '../data/zno/band.yaml'
kappafile = '../data/zno/kappa-m404021.hdf5'
poscar = '../data/zno/POSCAR'

projected = 'lifetime'
quantities = ['frequency', projected, 'dispersion']
temperature = 300
colour = 'viridis'

# Axes

fig, ax = tp.axes.one.colourbar()

# Load

data = tp.data.load.phono3py(kappafile, quantities=quantities)
pdata = tp.data.load.phonopy_dispersion(phile)

# Add

tp.plot.phonons.add_projected_dispersion(ax, data, pdata, projected,
                                         temperature=temperature,
                                         colour=colour, poscar=poscar)

# Save

plt.savefig('prophon.pdf')
plt.savefig('prophon.png')
