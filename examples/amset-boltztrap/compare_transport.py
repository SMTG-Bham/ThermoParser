#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tp

afile = '../data/amset_data_85x85x47.json'
bfile = '../data/boltztrap.hdf5'
qs = ['conductivity', 'seebeck', 'electronic_thermal_conductivity']
scale = ['log', 'linear', 'log']
direction = 'avg'

leg = ['10$\mathrm{{^{{{}}}}}$'.format(i) for i in range(18, 22)]
acmap = plt.get_cmap('winter')(np.linspace(0, 1, 4))
bcmap = plt.get_cmap('autumn')(np.linspace(0, 1, 4))

# Axes

fig, ax, add_legend = tp.axes.three.h_top_legend()

# Load

adata = tp.data.load.amset(afile, qs)
adata = tp.data.resolve.resolve(adata, qs, direction=direction)

bdata = tp.data.load.boltztrap(bfile, qs)
bdata = tp.data.resolve.resolve(bdata, qs, direction=direction)

bdi = [np.where(np.round(bdata['doping']) == np.round(1e18))[0],
       np.where(np.round(bdata['doping']) == np.round(1e19))[0],
       np.where(np.round(bdata['doping']) == np.round(1e20))[0],
       np.where(np.round(bdata['doping']) == np.round(1e21))[0]]

# Add

axlabels = tp.settings.labels()

for q in range(len(qs)):
    for i in range(4):
        ax[q].plot(adata['temperature'], adata[qs[q]][:,i+2],c=acmap[i],
                   label='A-{}'.format(leg[i]))
        ax[q].plot(bdata['temperature'], bdata[qs[q]][:,bdi[i]], c=bcmap[i],
                   label='B-{}'.format(leg[i]))
    ax[q].set_xlabel(axlabels['temperature'])
    ax[q].set_ylabel(axlabels[qs[q]])

for a in range(len(ax)):
    tp.plot.utilities.set_locators(ax[a], x='linear', y=scale[a])

# Save

add_legend(title=axlabels['doping'])
plt.savefig('compare-transport.pdf')
