#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = '../data/band.yaml'
dosfile = '../data/projected_dos.dat'
poscar = '../data/POSCAR'

colour = '#ff8000'
colours = {'Sb': '#00ff00',
           'Mg': '#ff00ff'}

# Axes

fig, ax = tp.plot.axes.one_dos_small_legend()

# Load

dispersion = tp.data.load.phonopy_dispersion(phile)
dos = tp.data.load.phonopy_dos(dosfile, poscar=poscar)

# Add

tp.plot.phonons.add_dispersion(ax[0], dispersion, colour=colour)
tp.plot.frequency.add_dos(ax[1], dos, colours, invert=True, line=True)

ax[1].set_ylim(ax[0].get_ylim())
tp.plot.utilities.set_locators(ax[1], dos=True)

# Save

ax[1].legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.savefig('phonons.pdf')
