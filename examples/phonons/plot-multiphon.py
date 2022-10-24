#!/usr/bin/env python3

import tp

scs = '222 333 444 555'.split()
files = ['../data/basno3/band-{}.yaml'.format(s) for s in scs]
label = [' x '.join(list(s)) for s in scs]
legend_title = 'Supercell Size'

# Axes
fig, ax, add_legend = tp.axes.large.one()

# Load
data = [tp.data.load.phonopy_dispersion(f) for f in files]

# Add
tp.plot.phonons.add_multi(ax, data, label=label)
add_legend(title=legend_title)

# Save
fig.savefig('multiphon.pdf')
fig.savefig('multiphon.png')
