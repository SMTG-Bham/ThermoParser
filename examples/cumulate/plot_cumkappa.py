#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

kappafile = '../data/kappa-m505028.hdf5'
dosfile = '../data/projected_dos.dat'
poscar = '../data/POSCAR'
direction = 'avg'
quantities = 'mode_kappa frequency'

colours = {'Sb': '#00ff00',
           'Mg': '#800080'}

# Axes

fig, ax = tp.plot.axes.one_small_legend()

# Load

data = tp.data.load.phono3py(kappafile, quantities=quantities)
dos = tp.data.load.phonopy_dos(dosfile, poscar=poscar)

# Add

tp.plot.frequency.add_cum_kappa(ax, data, direction=direction, main=True)
tp.plot.frequency.add_dos(ax, dos, colours, main=False, scale=True)

# Save

ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.savefig('cumkappa-{}.pdf'.format(direction))
