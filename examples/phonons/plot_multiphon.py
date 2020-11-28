#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

files = ['../data/band.yaml', '../data/band.yaml']

label = ['1\ x\ 1\ x\ 1', '2\ x\ 2\ x\ 2']
legend_title = 'Supercell Size'
colour = 'winter_r'
linestyle = ['-', ':']

# Axes

fig, ax, add_legend = tp.axes.one.medium_legend()

# Load

data = [tp.data.load.phonopy_dispersion(f) for f in files]

# Add

tp.plot.phonons.add_multi(ax, data, colour=colour, label=label,
                          linestyle=linestyle)
add_legend(title=legend_title)

# Save

plt.savefig('multiphon.pdf')
