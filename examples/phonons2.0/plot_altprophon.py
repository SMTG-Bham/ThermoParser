#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import tp

phile = '../data/band.yaml'
kappafile = '../data/kappa-m505028.hdf5'
poscar = '../data/POSCAR'
dispersion = 'group_velocity'
projected = 'mode_kappa'
direction = 'norm'
temperature = 300
quantities = [dispersion, projected, 'dispersion']

# Axes
fig, ax = tp.plot.axes.one_colourbar()

# Load
data = tp.data.load.phono3py(kappafile, quantities=quantities)
pdata = tp.data.load.phonopy_dispersion(phile)

# Add
tp.plot.phonons.add_alt_projected_dispersion(ax, data, pdata, dispersion,
                                             projected, direction=direction,
                                             temperature=temperature,
                                             poscar=poscar)

# Save
plt.savefig('altprophon.pdf')
