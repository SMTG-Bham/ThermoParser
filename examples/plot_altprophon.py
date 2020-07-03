#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import tp

phile = 'band.yaml'
kappafile = 'kappa-m505028.hdf5'
direction = 'norm'
temperature = 50

# Axes
fig, ax = tp.plot.axes.one_colourbar()

# Load
data = tp.data.load.phono3py(kappafile, quantities=['occupation', 'mode_kappa',
                                                    'temperature', 'qpoint'])
pdata = tp.data.load.phonopy_dispersion(phile)

# Add
ax = tp.plot.phonons.add_alt_projected_dispersion(ax, data, pdata,
                                                  'mode_kappa', 'occupation',
                                                  direction=direction,
                                                  temperature=temperature)

# Save
plt.savefig('altprophon.pdf')
