#!/usr/bin/env python3

import tp

kappafile = '../data/zno/kappa-m404021.hdf5'
temperature = 300
direction = ['x', 'z']
quantities = 'frequency mode_kappa mfp'
# Note for cumkappa and waterfall plots, mode_kappa and not kappa is required

# You can ignore this section
from os import path
if not path.isfile(kappafile) or (path.getsize(kappafile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh in the folder above.')
# End ignore

dosfile = '../data/zno/projected_dos.dat'
poscar = '../data/zno/POSCAR'

colour = ['#59c605', '#ffcf06']
colours = {'Zn': '#d46ef9',
           'O':  '#7b8eff'}
# Axes
fig, ax, add_legend = tp.axes.small.two_h()

# Load
data = tp.data.load.phono3py(kappafile, quantities=quantities)
dos = tp.data.load.phonopy_dos(dosfile, poscar=poscar)

# Add
tp.plot.frequency.add_cum_kappa(ax[0], data, temperature=temperature,
                                direction=direction, colour=colour)
tp.plot.mfp.add_cum_kappa(ax[1], data, temperature=temperature, scale=True,
                          direction=direction, colour=colour, xmarkers=1e-7)
tp.plot.frequency.add_dos(ax[0], dos, colour=colours, scale=True, main=False)

add_legend(location=2, ncol=2)

# Save
fig.savefig('cumkappa.pdf')
fig.savefig('cumkappa.png')
