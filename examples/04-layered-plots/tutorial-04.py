#!/usr/bin/env python3

import tp

scs = '222 333 444 555'.split()
pfiles = ['../data/basno3/band-{}.yaml'.format(s) for s in scs]
kfile = '../data/basno3/kappa-m363636.hdf5'
dfile = '../data/basno3/projected_dos.dat'
poscar = '../data/basno3/POSCAR'

direction = 'avg'
temperature = 300
waterfall = 'mean_free_path'
quantities = ['waterfall', waterfall]

colour = 'winter_r'
colours = {'Ba': '#ffcf06',
           'Sn': '#59c605',
           'O':  '#00b1f7'}
cmap = 'viridis'

# You can ignore down to line 23!
from os import path
if not path.isfile(kfile) or (path.getsize(kfile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh in the folder above.')
# Stop ignoring!

# Axes
fig, ax, add_legend = tp.axes.small.two_h()

# Load
dispersions = [tp.data.load.phonopy_dispersion(f) for f in pfiles]
kappa = tp.data.load.phono3py(kfile, quantities=quantities)
dos = tp.data.load.phonopy_dos(dfile, poscar=poscar)

# Plot
tp.plot.phonons.add_multi(ax[0], dispersions, colour=colour, label=scs)
tp.plot.frequency.format_waterfall(ax[1], kappa, waterfall, direction=direction,
                                   temperature=temperature, invert=True)
tp.plot.frequency.add_dos(ax[1], dos, colour=colours, scale=True, main=False,
                          alpha=0.6, line=False, invert=True)
tp.plot.frequency.add_waterfall(ax[1], kappa, waterfall, colour=cmap,
                                direction=direction, temperature=temperature,
                                invert=True)

# Formatting

tp.plot.utilities.set_locators(ax[1], x='log', y='linear')
axlabels = tp.settings.labels()
ax[1].set_xlabel(axlabels['mean_free_path'])
ax[1].set_ylabel(axlabels['frequency'])
add_legend()

# Save
fig.savefig('tutorial-04.png')
