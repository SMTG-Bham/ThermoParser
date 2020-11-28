#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = '../data/band.yaml'
kappafile = '../data/kappa-m505028.hdf5'
poscar = '../data/POSCAR'
direction = 'norm'
temperature = 300
projected = 'group_velocity'
quantities = ['frequency', 'dispersion', projected]

# Axes

fig, ax, add_legend = tp.axes.one.wide_large_legend()

# Load

data = tp.data.load.phono3py(kappafile, quantities=quantities)
pdata = tp.data.load.phonopy_dispersion(phile)

# Add

tp.plot.phonons.add_alt_dispersion(ax, data, pdata, 'group_velocity',
                                   direction=direction,
                                   temperature=temperature, poscar=poscar)
add_legend()

# Save

plt.savefig('altphon.pdf')
