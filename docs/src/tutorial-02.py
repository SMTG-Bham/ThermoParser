#!/usr/bin/env python3

import tp

# Variables
pfile = '../data/zno/band.yaml'
kfile = '../data/zno/kappa-m404021.hdf5'
poscar = '../data/zno/POSCAR'
temperature = 300

colour = ['#000000', '#ff0000']

# You can ignore down to line 20!
from os import path
if not path.isfile(kfile) or (path.getsize(kfile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh')
# Stop ignoring!

# Axes
fig, ax, _ = tp.axes.large.one('dark_background')

# Load
kdata = tp.data.load.phono3py(kfile, quantities='wideband')
pdata = tp.data.load.phonopy_dispersion(pfile)

# Plot
tp.plot.phonons.add_wideband(ax, kdata, pdata, temperature=temperature,
                             colour=colour, poscar=poscar)

# Save
fig.savefig('tutorial-02.png')
