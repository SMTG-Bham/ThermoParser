#!/usr/bin/env python3

from os import path
import tp

pfile = '../data/zno/band.yaml'
kfile = '../data/zno/kappa-m404021.hdf5'
poscar = '../data/zno/POSCAR'
if not path.isfile(kfile) or (path.getsize(kfile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh')

projected = 'lifetime'
quantities = ['frequency', projected, 'dispersion']
temperature = 300
colour = 'viridis'

# Axes
fig, ax, add_legend = tp.axes.small.one_colourbar()

# Load
data = tp.data.load.phono3py(kfile, quantities=quantities)
pdata = tp.data.load.phonopy_dispersion(pfile)

# Add
tp.plot.phonons.add_projected_dispersion(ax, data, pdata, projected,
                                         temperature=temperature,
                                         colour=colour, poscar=poscar)

# Save

fig.savefig('prophon.pdf')
fig.savefig('prophon.png')
