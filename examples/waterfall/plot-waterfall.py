#!/usr/bin/env python3

from matplotlib import cm
import matplotlib.pyplot as plt
from os import path
import tp

kappafile = '../data/basno3/kappa-m363636.hdf5'
dosfile = '../data/basno3/projected_dos.dat'
poscar = '../data/basno3/POSCAR'
if not path.isfile(kappafile) or (path.getsize(kappafile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh')

direction = 'avg'
temperature = 300
waterfall = 'mean_free_path'
projected = 'mode_kappa'
quantities = ['waterfall', waterfall, projected]

colours = {'Ba': '#ffcf06',
           'Sn': '#59c605',
           'O':  '#00b1f7'}
colour = cm.get_cmap('viridis')

# Axes

"""
Plots with legends return an add_legend function, which will place the
legend nicely and still accepts arguments like title. This doesn't stop
you using plt.legend instead (see line 50).
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
tp.plot.frequency.add_dos(ax, dos, colour=colours, scale=True, main=False, alpha=0.6)
cbar = tp.plot.frequency.add_projected_waterfall(ax, data, waterfall,
                                                 projected, colour=colour,
                                                 temperature=temperature,
                                                 direction=direction, s=50)
add_legend()

"""
Some of the longer names don't fit on presentation style plots, but it's
easy to switch. Here large_labels have been used (for large axes), which
points to mid-length labels in the settings file by default, but
short_labels could also be used, which uses a symbol (e.g. kappa). The
default can be changed in the settings file.
"""
labels = tp.settings.large_labels()
cbar.set_label(labels[projected])

# Save

plt.savefig('waterfall.pdf')
plt.savefig('waterfall.png')
