#!/usr/bin/env python3


import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tp

afile = 'amset_data_85x85x47.json'
bfile = 'boltztrap.hdf5'
run_boltztrap = True
direction = 'avg'
acmap = plt.get_cmap('winter')(np.linspace(0, 1, 4))
bcmap = plt.get_cmap('autumn')(np.linspace(0, 1, 4))
leg = ['10$\mathrm{{^{{{}}}}}$'.format(i) for i in range(18, 22)]
quantities = 'conductivity seebeck'# ke'
aquantities = 'conductivity seebeck ke'

if run_boltztrap:
    tp.data.run.boltztrap(zero_weighted=False, run=True)

# Axes

fig, ax = tp.plot.axes.three_h_top_legend()

# Load

adata = tp.data.load.amset(afile, aquantities)
adata = tp.data.resolve.resolve(adata, aquantities, direction=direction)
bdata = tp.data.load.boltztrap(bfile, quantities)
bdata = tp.data.resolve.resolve(bdata, quantities, direction=direction)
adi = [np.where(np.round(adata['doping']) == np.round(-1e18))[0],
       np.where(np.round(adata['doping']) == np.round(-1e19))[0],
       np.where(np.round(adata['doping']) == np.round(-1e20))[0],
       np.where(np.round(adata['doping']) == np.round(-1e21))[0]]
bdi = [np.where(np.round(bdata['doping']) == np.round(1e18))[0],
       np.where(np.round(bdata['doping']) == np.round(1e19))[0],
       np.where(np.round(bdata['doping']) == np.round(1e20))[0],
       np.where(np.round(bdata['doping']) == np.round(1e21))[0]]
print(bdata['doping'][:])
# Add

for i in range(4):
    ax[0].plot(adata['temperature'], adata['conductivity'][:,i+2],
               c=acmap[i], label='A-{}'.format(leg[i]))
    ax[1].plot(adata['temperature'], adata['seebeck'][:,i+2],
               c=acmap[i], label='A-{}'.format(leg[i]))
    ax[2].plot(adata['temperature'], adata['electronic_thermal_conductivity'][:,i+2],
               c=acmap[i], label='A-{}'.format(leg[i]))
for i in range(4):
    ax[0].plot(bdata['temperature'], bdata['conductivity'][:,bdi[i]],
               c=bcmap[i], label='B-{}'.format(leg[i]))
    ax[1].plot(bdata['temperature'], bdata['seebeck'][:,bdi[i]],
               c=bcmap[i], label='B-{}'.format(leg[i]))
    #    ax[2].plot(bdata['temperature'], bdata['electronic_thermal_conductivity'][:,bdi[i]],
    #           c=bcmap[i], label='B-{}'.format(leg[i]))

ax[0].set_yscale('log')
ax[2].set_yscale('log')
axlabels = tp.settings.labels()
ax[0].set_xlabel(axlabels['temperature'])
ax[1].set_xlabel(axlabels['temperature'])
ax[2].set_xlabel(axlabels['temperature'])
ax[0].set_ylabel(axlabels['conductivity'])
ax[1].set_ylabel(axlabels['seebeck'])
ax[2].set_ylabel(axlabels['electronic_thermal_conductivity'])

ax[0].xaxis.set_major_locator(ticker.MaxNLocator(4))
ax[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax[0].yaxis.set_major_locator(ticker.LogLocator())

ax[1].xaxis.set_major_locator(ticker.MaxNLocator(4))
ax[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax[1].yaxis.set_major_locator(ticker.MaxNLocator(4))
ax[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())

ax[2].xaxis.set_major_locator(ticker.MaxNLocator(4))
ax[2].xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax[2].yaxis.set_major_locator(ticker.LogLocator())

# Save

legend = ax[1].legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=8,
                      title=axlabels['doping'])
plt.savefig('compare-transport.pdf')
