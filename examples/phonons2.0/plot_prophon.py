#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = '../data/band.yaml'
kappafile = '../data/kappa-m505028.hdf5'
poscar = '../data/POSCAR'
projected = 'mean_free_path'
waterfall = 'mode_kappa'
direction = 'norm'
temperature = 300
quantities=[projected, waterfall, 'dispersion', 'waterfall']

colour = 'viridis_r'

# Axes

fig, ax = tp.plot.axes.one_dos_colourbar()

# Load

data = tp.data.load.phono3py(kappafile, quantities=quantities)
pdata = tp.data.load.phonopy_dispersion(phile)

# Add

cbar = tp.plot.phonons.add_projected_dispersion(ax[0], data, pdata, projected,
                                                colour=colour, poscar=poscar,
                                                temperature=temperature,
                                                direction=direction)
cbar.remove()
cbar = tp.plot.frequency.add_projected_waterfall(ax[1], data, waterfall,
                                                 projected, colour=colour,
                                                 temperature=temperature,
                                                 direction=direction,
                                                 invert=True)

ax[1].set_xlabel('$\mathrm{\kappa_l}\ (W\ m^{-1}\ K^{-1})$')
ax[1].set_ylim(ax[0].get_ylim())
tp.settings.set_locators(ax[1], dos=True, x='log')

# Save

plt.savefig('prophon.pdf')
