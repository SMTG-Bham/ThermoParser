#!/usr/bin/env python3

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as con
import tp

kappafile = 'kappa-m505028.hdf5'
temperature = 10

hbar = con.physical_constants['Planck constant over 2 pi in eV s'][0]
kb = con.physical_constants['Boltzmann constant in eV/K'][0]

data = tp.data.load.phono3py(kappafile, quantities=['frequency', 'occupation',
                                                    'temperature'])
o1 = np.divide(np.multiply(data['frequency'], 1e12 * hbar),
                           kb * temperature)
o2 = np.expm1(o1)
o3 = np.reciprocal(o2)
o4 = data['occupation'][np.where(data['temperature'][:] == temperature)[0][0]]
f = data['frequency']

fig, ax = tp.plot.axes.four_square()

ax[0][0].scatter(f, o1, rasterized=True)
ax[0][1].scatter(f, o2, rasterized=True)
ax[1][0].scatter(f, o3, rasterized=True)
ax[1][1].scatter(f, o4, rasterized=True)

o3sort = np.ravel(o3)
o3sort = o3sort[o3sort.argsort()]
ax[1][0].set_ylim(o3sort[int(round(len(o3sort)/100, 0))],
                  o3sort[int(round(len(o3sort)*99/100, 0))])
o4sort = np.ravel(o4)
o4sort = o4sort[o4sort.argsort()]
ax[1][1].set_ylim(1/(o3sort[int(round(len(o3sort)/100, 0))]),
                  1/(o3sort[int(round(len(o3sort)*99/100, 0))]))

plt.savefig('occupation-test.pdf')
