#!/usr/bin/env python3

import tp
from matplotlib import pyplot as plt

afile = 'amset_data_85x85x47.json'
aquants = ['temperature', 'doping', 'seebeck', 'conductivity',
           'electronic_thermal_conductivity']
kfile = 'kappa-m505028.hdf5'
kquants = ['kappa', 'temperature']


# Axes

fig, ax = tp.plot.axes.one_colourbar()

# Load

data = tp.data.load.amset(afile, quantities=aquants)
kdata = tp.data.load.phono3py(kfile, quantities=kquants)

# Add

ax, cbar = tp.plot.heatmap.add_ztmap(ax, data, kdata, colour='#ff8080',
                                     direction='x', xinterp=200, yinterp=200)

# Save

plt.savefig('ztmap.pdf')
