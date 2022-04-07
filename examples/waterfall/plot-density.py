#!/usr/bin/env python3

from os import path
import tp

kfile = '../data/basno3/kappa-m363636.hdf5'
if not path.isfile(kfile) or (path.getsize(kfile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh')

direction = 'avg'
temperature = 300
waterfall = 'mean_free_path'
quantities = ['waterfall', waterfall]
# waterfall is an alias for frequency

colours = {'Ba': '#ffcf06',
           'Sn': '#59c605',
           'O':  '#00b1f7'}
colour = 'Blues'

# Axes
fig, ax, add_legend = tp.axes.small.one()

# Load
data = tp.data.load.phono3py(kfile, quantities=quantities)

# Add
tp.plot.frequency.add_density(ax, data, waterfall, colour=colour,
                              temperature=temperature, direction=direction)

# Save

fig.savefig('density.pdf')
fig.savefig('density.png')
