#!/usr/bin/env python3

import tp

pfile = '../data/zno/band.yaml'
dfile = '../data/zno/projected_dos.dat'
poscar = '../data/zno/POSCAR'

colour = '#f0901f'
colours = {'Zn': '#d46ef9',
           'O':  '#7b8eff'}
location = 2

# Axes
fig, ax, add_legend = tp.axes.small.one_dos()

# Load
dispersion = tp.data.load.phonopy_dispersion(pfile)
dos = tp.data.load.phonopy_dos(dfile, poscar=poscar)

# Plot
tp.plot.phonons.add_dispersion(ax[0], dispersion, colour=colour)
tp.plot.frequency.add_dos(ax[1], dos, colour=colours, invert=True, line=True)

# Formatting
ax[1].set_ylim(ax[0].get_ylim())
add_legend(title='ZnO')

# Save
fig.savefig('tutorial-03.png')
