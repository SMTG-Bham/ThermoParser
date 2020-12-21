#!/usr/bin/env python3

import tp
from matplotlib import pyplot as plt

bfile = 'data/zno/boltztrap.hdf5'
kfile = 'data/zno/kappa-m404021.hdf5'
direction = 'x'
colour = '#800080'

# Axes

fig, ax = tp.axes.one_large.colourbar()

# Load

adata = tp.data.load.boltztrap(bfile)
kdata = tp.data.load.phono3py(kfile)

# Add

tp.plot.heatmap.add_ztmap(ax, adata, kdata, direction=direction, colour=colour)

# Save

plt.savefig('ztmap.pdf')
plt.savefig('ztmap.png')
