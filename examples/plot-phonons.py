#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = 'data/zno/band.yaml'
dosfile = 'data/zno/projected_dos.dat'
poscar = 'data/zno/POSCAR'

colour = '#ff8000'
colours = {'Zn': '#ffff00',
           'O':  '#00ffff'}

# Axes

fig, ax, add_legend = tp.axes.one.dos_small_legend()

# Load

dispersion = tp.data.load.phonopy_dispersion(phile)
dos = tp.data.load.phonopy_dos(dosfile, poscar=poscar)

# Add

tp.plot.phonons.add_dispersion(ax[0], dispersion, colour=colour)
tp.plot.frequency.add_dos(ax[1], dos, colour=colours, invert=True, line=True)

ax[1].set_ylim(ax[0].get_ylim())
tp.plot.utilities.set_locators(ax[1], dos=True)
add_legend()

# Save

plt.savefig('phonons.pdf')
