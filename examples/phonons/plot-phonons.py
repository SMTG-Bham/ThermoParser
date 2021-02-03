#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

phile = '../data/zno/band.yaml'
dosfile = '../data/zno/projected_dos.dat'
poscar = '../data/zno/POSCAR'

colour = '#f0901f'
colours = {'Zn': '#d46ef9',
           'O':  '#7b8eff'}

# Axes

"""
Plots with legends return an add_legend function, which will place the
legend nicely and still accepts arguments like title. This doesn't stop
you using plt.legend instead (see line 41).
"""
fig, ax, add_legend = tp.axes.one_large.dos_small_legend()

# Load

dispersion = tp.data.load.phonopy_dispersion(phile)
dos = tp.data.load.phonopy_dos(dosfile, poscar=poscar)

# Add

tp.plot.phonons.add_dispersion(ax[0], dispersion, colour=colour)
tp.plot.frequency.add_dos(ax[1], dos, colour=colours, invert=True, line=True)

"""
DoS-style plots need to have their y-axis rescaled to match the main
plot (line 39). They can also have their labels removed (line 40). If
something has a meaningful x-axis quantity, ticks can be reimplemented
by adding xscale='linear' or xscale='log' to line 40.
"""
ax[1].set_ylim(ax[0].get_ylim())
add_legend()

# Save

plt.savefig('phonons.pdf')
plt.savefig('phonons.png')
