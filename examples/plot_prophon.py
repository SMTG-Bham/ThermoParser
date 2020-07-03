#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = 'band.yaml'
kappafile = 'kappa-m505028.hdf5'
projected = 'mean_free_path'
waterfall = 'mode_kappa'
direction = 'norm'
temperature = 300
colour = 'viridis_r'

# Axes
fig, ax = tp.plot.axes.one_dos_colourbar()

# Load
data = tp.data.load.phono3py(kappafile, temperature=temperature,
                             direction=direction,
                             quantities=['mean_free_path', 'frequency',
                                         'temperature', 'qpoint', 'mode_kappa'])
pdata = tp.data.load.phonopy_dispersion(phile)

# Add
ax[0], cbar = tp.plot.phonons.add_projected_dispersion(ax[0], data, pdata,
                                                       projected, colour=colour)
ax[1], cbar = tp.plot.frequency.add_projected_waterfall(ax[1], data, waterfall,
                                                        projected, colour=colour,
                                                        invert=True)
cbar.remove()

ax[1].set_xlabel('$\mathrm{\kappa_l}$')
ax[1].set_ylim(ax[0].get_ylim())
ax[1].set_yticks([])
ax[1].set_yticklabels([])
ax[1].set_ylabel('')

plt.savefig('prophon.pdf')
