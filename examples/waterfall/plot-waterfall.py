#!/usr/bin/env python3

from os import path
import tp

kfile = '../data/basno3/kappa-m363636.hdf5'
dfile = '../data/basno3/projected_dos.dat'
poscar = '../data/basno3/POSCAR'
if not path.isfile(kfile) or (path.getsize(kfile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh')

direction = 'avg'
temperature = 300
waterfall = 'mean_free_path'
projected = 'mode_kappa'
quantities = ['waterfall', waterfall, projected]

colours = {'Ba': '#ffcf06',
           'Sn': '#59c605',
           'O':  '#00b1f7'}
colour = 'viridis'

# Axes
fig, ax, add_legend = tp.axes.small.one_colourbar()

# Load
data = tp.data.load.phono3py(kfile, quantities=quantities)
dos = tp.data.load.phonopy_dos(dfile, poscar=poscar)

# Add
tp.plot.frequency.format_waterfall(ax, data, waterfall, direction=direction,
                                   temperature=temperature)
tp.plot.frequency.add_dos(ax, dos, colour=colours, scale=True, main=False,
                          alpha=0.6, line=False)
cbar = tp.plot.frequency.add_projected_waterfall(ax, data, waterfall,
                                                 projected, colour=colour,
                                                 temperature=temperature,
                                                 direction=direction)
add_legend()

# Save

fig.savefig('waterfall.pdf')
fig.savefig('waterfall.png')
