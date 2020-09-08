#!/usr/bin/env python3

import tp
from matplotlib import pyplot as plt

afile = '../data/amset_data_85x85x47.json'
aquants = 'seebeck conductivity ke'
kfile = '../data/kappa-m505028.hdf5'
kquants = 'kappa'
direction = 'x'

colour = '#ff8080'

# Axes

fig, ax = tp.plot.axes.one_colourbar()

# Load

adata = tp.data.load.amset(afile, quantities=aquants)
kdata = tp.data.load.phono3py(kfile, quantities=kquants)

# Add

ax, cbar = tp.plot.heatmap.add_ztmap(ax, adata, kdata, direction=direction,
                                     colour=colour, xinterp=200, yinterp=200)

# Save

plt.savefig('ztmap-a.pdf')
