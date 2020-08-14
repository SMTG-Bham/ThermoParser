#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tp

afile = 'amset_data_85x85x47.json'
direction = 'avg'
cmap = plt.get_cmap('winter')(np.linspace(0, 1, 6))
leg = ['10$\mathrm{{^{{{}}}}}$'.format(i) for i in range(16, 22)]
quantities = 'conductivity seebeck ke'

# Axes

fig, ax = tp.plot.axes.three_h_top_legend()

# Load

data = tp.data.load.amset(afile)
data = tp.data.resolve.resolve(data, quantities, direction=direction)

# Add

for i, d in enumerate(data['doping']):
    ax[0].plot(data['temperature'], data['conductivity'][i], c=cmap[i], label=leg[i])
    ax[1].plot(data['temperature'], data['seebeck'][i], c=cmap[i], label=leg[i])
    ax[2].plot(data['temperature'], data['electronic_thermal_conductivity'][i],
                                                       c=cmap[i], label=leg[i])

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

legend = ax[1].legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=6,
                      title=axlabels['doping'])
plt.savefig('transport.pdf')
