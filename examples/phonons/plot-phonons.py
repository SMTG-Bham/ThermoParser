#!/usr/bin/env python3

import tp

phile = '../data/zno/band.yaml'
dosfile = '../data/zno/projected_dos.dat'
poscar = '../data/zno/POSCAR'

colour = '#f0901f'
colours = {'Zn': '#d46ef9',
           'O':  '#7b8eff'}

# Axes
fig, ax, add_legend = tp.axes.large.one_dos()

# Load
dispersion = tp.data.load.phonopy_dispersion(phile)
dos = tp.data.load.phonopy_dos(dosfile, poscar=poscar)

# Add
tp.plot.phonons.add_dispersion(ax[0], dispersion, colour=colour)
tp.plot.frequency.add_dos(ax[1], dos, colour=colours, invert=True, line=True)

ax[1].set_ylim(ax[0].get_ylim())
add_legend()

# Save
fig.savefig('phonons.pdf')
fig.savefig('phonons.png')
