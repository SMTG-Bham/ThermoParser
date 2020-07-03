#!/usr/bin/env python3

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import tp

kappafile = 'kappa-m505028.hdf5'
temperature = 10

data = tp.data.load.phono3py(kappafile, quantities=['frequency', 'occupation',
                                                    'temperature'])
occ = data['occupation'][np.where(data['temperature'][:] == temperature)[0][0]]
occ = np.ravel(occ)
f = np.ravel(data['frequency'])
tocc = [len(np.where(data['occupation'][t] < 1)[0]) for t in range(len(data['temperature']))]

fig, ax = tp.plot.axes.two_h()

ax[0].scatter(f[15:], occ[15:], rasterized=True)
ax[1].scatter(data['temperature'], tocc, rasterized=True)
ax[1].set_yscale('log')

plt.savefig('occupation.pdf')
