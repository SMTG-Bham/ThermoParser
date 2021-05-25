#!/usr/bin/env python3

import matplotlib.pyplot as plt
from os import path
import tp

kappafile = '../data/zno/kappa-m404021.hdf5'
temperature = 300
direction = ['x', 'z']
quantities = 'frequency mode_kappa mfp'
if not path.isfile(kappafile) or (path.getsize(kappafile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh')

dosfile = '../data/zno/projected_dos.dat'
poscar = '../data/zno/POSCAR'

colour = ['#59c605', '#ffcf06']
colours = {'Zn': '#d46ef9',
           'O':  '#7b8eff'}
# Axes

fig, ax, add_legend = tp.axes.two.h_small_legend()

# Load

data = tp.data.load.phono3py(kappafile, quantities=quantities)
dos = tp.data.load.phonopy_dos(dosfile, poscar=poscar)

# Add

tp.plot.frequency.add_cum_kappa(ax[0], data, temperature=temperature,
                                direction=direction, colour=colour)
tp.plot.mfp.add_cum_kappa(ax[1], data, temperature=temperature,
                          direction=direction, colour=colour, xmarkers=1e-7)
tp.plot.frequency.add_dos(ax[0], dos, colour=colours, scale=True, main=False)

add_legend()

# Save

plt.savefig('cumkappa.pdf')
plt.savefig('cumkappa.png')
