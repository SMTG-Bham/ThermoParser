#!/usr/bin/env python3

import tp

dosfile = '../data/basno3/projected_dos.dat'
poscar = '../data/basno3/POSCAR'

colours = {'Ba':  '#ff00ff',
           'Sn':  '#00ffff',
           'O':   '#ff0000',
           'O_2': '#ff8000'}

fig, ax, add_legend = tp.axes.small.two_h()

# Unsmeared, read from POSCAR with total
dos1 = tp.data.load.phonopy_dos(dosfile, poscar=poscar)
tp.plot.frequency.add_dos(ax[0], dos1, total=True, colour=colours)

# Smeared (sigma=0.2), custom atoms with no total
dos2 = tp.data.load.phonopy_dos(dosfile, atoms='Ba Sn O O_2 2')
tp.plot.frequency.add_dos(ax[1], dos2, sigma=0.2, colour=colours)

add_legend(title='$BaSnO_3$', location=2, ncol=1)

fig.savefig('dos.png')