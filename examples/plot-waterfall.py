#!/usr/bin/env python3

from matplotlib import cm
import matplotlib.pyplot as plt
import tp

kappafile = 'data/basno3/kappa-m363636.hdf5'
dosfile = 'data/basno3/projected_dos.dat'
poscar = 'data/basno3/POSCAR'
direction = 'avg'
temperature = 300
waterfall = 'mean_free_path'
projected = 'mode_kappa'
quantities = ['waterfall', waterfall, projected]

colours = {'Ba': '#800080',
           'Sn': '#ff8080',
           'O':  '#00ffff'}
colour = cm.get_cmap('viridis')

# Axes

fig, ax, add_legend = tp.axes.one.colourbar_small_legend()

# Load

data = tp.data.load.phono3py(kappafile, quantities=quantities)
dos = tp.data.load.phonopy_dos(dosfile, poscar=poscar)

# Add

"""
First the axes limits are set for the waterfall plot, so the dos can
be scaled to the limits. Once the dos is plotted, the waterfall is
plotted on top, so the dos colour doesn't obscure or distort the
waterfall colour.
"""
tp.plot.frequency.format_waterfall(ax, data, waterfall, direction=direction,
                                   temperature=temperature)
tp.plot.frequency.add_dos(ax, dos, colour=colours, scale=True, main=False)
cbar = tp.plot.frequency.add_projected_waterfall(ax, data, waterfall,
                                                 projected, main=True,
                                                 colour=colour,
                                                 temperature=temperature,
                                                 direction=direction)
add_legend()

# Save

plt.savefig('waterfall.pdf')
