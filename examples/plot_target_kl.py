#!/usr/bin/env python3

import tp
from matplotlib import pyplot as plt

f = 'amset_data_85x85x47.json'
quants = ['temperature', 'doping', 'seebeck', 'conductivity',
          'electronic_thermal_conductivity']
colour = tp.plot.colour.elbow('#ff8080')
colour = tp.plot.colour.highlight(colour, 'grey')

# Axes

fig, ax = tp.plot.axes.one_colourbar()

# Load

data = tp.data.load.amset(f, quantities=quants)

# Add

ax, cbar = tp.plot.heatmap.add_kappa_target(ax, data, colour=colour,
                                            direction='x', xinterp=200,
                                            yinterp=200)

# Save

plt.savefig('target_kl.pdf')
