#!/usr/bin/env python3

import tp
from matplotlib import pyplot as plt

bfile = '../data/boltztrap.hdf5'
bquants = 'seebeck conductivity ke'
kfile = '../data/kappa-m505028.hdf5'
kquants = 'kappa'
direction = 'x'

colour = '#ff8080'

# Axes

fig, ax = tp.plot.axes.one_colourbar()

# Load

bdata = tp.data.load.boltztrap(bfile, quantities=bquants)
kdata = tp.data.load.phono3py(kfile, quantities=kquants)

# Add

cbar = tp.plot.heatmap.add_ztmap(ax, bdata, kdata, direction=direction,
                                 colour=colour, xinterp=200, yinterp=200)

# Save

plt.savefig('ztmap-b.pdf')
