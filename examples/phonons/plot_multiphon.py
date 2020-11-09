#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

files = ['../data/band.yaml', '../data/band.yaml']

label = ['1\ x\ 1\ x\ 1', '2\ x\ 2\ x\ 2']
legend_title = 'Supercell Size'
colour = 'winter_r'
linestyle = ['-', ':']

# Axes

fig, ax = tp.plot.axes.one_medium_legend()

# Load

data = [tp.data.load.phonopy_dispersion(f) for f in files]

# Add

tp.plot.phonons.add_multi(ax, data, colour=colour, label=label,
                          linestyle=linestyle)

# Save

ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), title=legend_title)
plt.savefig('multiphon.pdf')
