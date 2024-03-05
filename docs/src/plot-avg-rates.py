#!/usr/bin/env python3

import tp
import numpy as np

f = '../data/basno3/mesh_75x75x75.h5'
q = 'weighted_rates' # 'weighted_mfp' #
temperature = 1000
doping = 1e19
colour = {'IMP':   'red',
          'POP':   'blue',
          'Total': 'magenta'}

# Example only shenanegans
from os import path
if not path.isfile(f) or path.getsize(f) < 1024*1024*10:
    raise Exception('File not found, please use get-data.sh in the folder above.')
# End of example only shenanegans

# Axes
fig, ax, add_legend = tp.axes.large.two_h()

# Load
data = tp.data.load.amset_mesh(f, q)
data['doping'] = np.abs(data['doping'])
tdata = tp.data.utilities.resolve(data, q, doping=doping)
ddata = tp.data.utilities.resolve(data, q, temperature=temperature)

# Add
for i, rate in enumerate(data['stype']):
    ax[0].plot(tdata['temperature'], tdata[q][i], label=rate, color=colour[rate])
    ax[1].plot(ddata['doping'], ddata[q][i], label=rate, color=colour[rate])

axlabels = tp.settings.large_labels()
for a in ax:
    a.set_ylabel(axlabels[q])
ax[0].set_xlabel(axlabels['temperature'])
ax[1].set_xlabel(axlabels['doping'])
tp.plot.utilities.set_locators(ax[0], x='linear', y='log')
tp.plot.utilities.set_locators(ax[1], x='log', y='log')

add_legend(title='Rate', location=2)

# Save
fig.savefig('avg-rates.pdf')
fig.savefig('avg-rates.png')
