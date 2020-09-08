#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = '../data/band.yaml'
dosfile = '../data/projected_dos.dat'
atoms = ['Sb', 2, 'Mg', 3]

colour = '#ff8000'
colours = {'Sb': '#00ff00',
           'Mg': '#ff00ff'}

# Axes

fig, ax = tp.plot.axes.one_dos_small_legend()

# Load

dispersion = tp.data.load.phonopy_dispersion(phile)
dos = tp.data.load.phonopy_dos(dosfile, atoms)

# Add

ax[0] = tp.plot.phonons.add_dispersion(ax[0], dispersion, colour=colour)
ax[1] = tp.plot.frequency.add_dos(ax[1], dos, colours, invert=True)

ax[1].set_xlim(left=0)
ax[1].set_ylim(ax[0].get_ylim())
ax[1].set_xticks([])
ax[1].set_xticklabels([])
ax[1].set_yticks([])
ax[1].set_yticklabels([])

# Save

ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.savefig('phonons.pdf')
