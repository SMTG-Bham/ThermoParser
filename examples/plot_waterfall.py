#!/usr/bin/env python3

from matplotlib import cm
import matplotlib.pyplot as plt
import tp

kappafile = 'kappa-m505028.hdf5'
dosfile = 'projected_dos.dat'
direction = 'avg'
temperature = 20
waterfall = 'mode_kappa'
projected = 'occupation'
atoms = ['Sb', 2, 'Mg', 3]
colours = {'Sb': '#00ff00',
           'Mg': '#800080'}
colour = tp.plot.colour.highlight(cm.get_cmap('viridis'), 'grey')
quantities = ['waterfall', waterfall, projected]

# Axes
fig, ax = tp.plot.axes.one_colourbar_small_legend()

# Load
data = tp.data.load.phono3py(kappafile, quantities=quantities)
dos = tp.data.load.phonopy_dos(dosfile, atoms)

# Add
ax, cbar = tp.plot.frequency.add_projected_waterfall(ax, data, waterfall,
                                                     projected, main=True,
                                                     colour=colour,
                                                     temperature=temperature,
                                                     direction=direction)
ax = tp.plot.frequency.add_dos(ax, dos, colours)

ax.legend(loc="center left", bbox_to_anchor=(1.25, 0.5))

# Save
plt.savefig('waterfall-kappa.pdf')
