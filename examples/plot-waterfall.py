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

"""
Plots with legends return an add_legend function, which will place the
legend nicely and still accepts arguments like title. This doesn't stop
you using plt.legend instead (see line 51).
"""
fig, ax, add_legend = tp.axes.one_large.colourbar_small_legend()

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
                                                 projected, colour=colour,
                                                 temperature=temperature,
                                                 direction=direction)
add_legend()

"""
Some of the longer names don't fit on presentation style plots, but it's
easy to switch.
"""
labels = tp.settings.short_labels()
cbar.set_label(labels[projected])

# Save

plt.savefig('waterfall.pdf')
plt.savefig('waterfall.png')
