#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import tp

phile = 'band.yaml'
kappafile = 'kappa-m505028.hdf5'
dispersion = 'group_velocity'
projected = 'occupation'
direction = 'norm'
temperature = 50
quantities = [dispersion, projected, 'dispersion']

# Axes
fig, ax = tp.plot.axes.one_colourbar()

# Load
data = tp.data.load.phono3py(kappafile, quantities=quantities)
pdata = tp.data.load.phonopy_dispersion(phile)

# Add
ax = tp.plot.phonons.add_alt_projected_dispersion(ax, data, pdata, dispersion,
                                                  projected,
                                                  direction=direction,
                                                  temperature=temperature)

# Save
plt.savefig('altprophon.pdf')
