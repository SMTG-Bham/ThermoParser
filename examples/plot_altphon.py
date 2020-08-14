#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = 'band.yaml'
kappafile = 'kappa-m505028.hdf5'
direction = 'norm'
temperature = 300
quantities = 'frequency gv dispersion'

# Axes
fig, ax = tp.plot.axes.one_wide_large_legend()

# Load
data = tp.data.load.phono3py(kappafile, quantities=quantities)
pdata = tp.data.load.phonopy_dispersion(phile)

# Add
ax = tp.plot.phonons.add_alt_dispersion(ax, data, pdata, 'group_velocity',
                                        direction=direction,
                                        temperature=temperature)

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Save
plt.savefig('altphon.pdf')
