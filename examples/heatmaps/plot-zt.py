#!/usr/bin/env python3

from os import path
import tp

bfile = '../data/zno/boltztrap.hdf5'
kfile = '../data/zno/kappa-m404021.hdf5'
if not path.isfile(kfile) or (path.getsize(kfile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh in the folder above.')

direction = 'x'
colour = '#800080'

# Axes
fig, ax, add_legend = tp.axes.small.one_colourbar()

# Load
adata = tp.data.load.boltztrap(bfile)
kdata = tp.data.load.phono3py(kfile)

# Add
tp.plot.heatmap.add_ztmap(ax, adata, kdata, direction=direction, colour=colour,
                          discrete=True)

# Save
fig.savefig('ztmap.pdf')
fig.savefig('ztmap.png')
