#!/usr/bin/env python3

import tp
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt

afile = '../data/basno3/transport_75x75x75.json'
kfile = '../data/basno3/kappa-m363636.hdf5'

doping = [-1e18, -1e19, -1e20, -1e21]
direction = 'avg'
quantities = ['conductivity', 'seebeck']
scale = ['log', 'linear', 'linear']

from os import path
if not path.isfile(kfile) or (path.getsize(kfile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh')

# Axes
plt.style.use('tp')
fig = plt.figure(figsize=(27/2.54, 8.3/2.54))
grid = GridSpec(1, 14)
ax = [fig.add_subplot(grid[0, :4]),
      fig.add_subplot(grid[0, 5:9]),
      fig.add_subplot(grid[0, 10:])]
plt.subplots_adjust(left=0.06, right=0.98,
                    bottom=0.12, top=0.95)

# Load
adata = tp.data.load.amset(afile)
kdata = tp.data.load.phono3py(kfile)

# Plot
for i, q in enumerate(quantities):
    for d in doping:
        data = tp.data.resolve.resolve(adata, q, direction=direction, doping=d)
        ax[i].plot(data['temperature'], data[q],
                   label="{:.2e}".format(data['meta']['doping']))

q = 'lattice_thermal_conductivity'
data = tp.data.resolve.resolve(kdata, q, direction=direction)
ax[2].plot(data['temperature'], data[q])

# Formatting
axlabels = tp.settings.labels()
for i, q in enumerate([*quantities, q]):
    ax[i].set_xlabel(axlabels['temperature'])
    ax[i].set_ylabel(axlabels[q])
    tp.plot.utilities.set_locators(ax[i], x='linear', y=scale[i])

handles, labels = tp.axes.legend.consolidate(ax)
ax[2].legend(loc='best', title=axlabels['doping'], handles=handles,
             labels=labels)
tp.axes.legend.alphabetise(ax, preset='roman', suffix=')', x=-0.12)

# Save
fig.savefig('transport.png')
