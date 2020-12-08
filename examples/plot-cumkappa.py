#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

kappafile = 'data/zno/kappa-m404021.hdf5'
temperature = 300
direction = ['x', 'z']
quantities = 'frequency mode_kappa mfp'

main = [False, True]
colour = ['#ff8000', '#0000ff']

# Axes

fig, ax, add_legend = tp.axes.two.h_small_legend()

# Load

data = tp.data.load.phono3py(kappafile, quantities=quantities)

# Add

for i in [0, 1]:
    tp.plot.frequency.add_cum_kappa(ax[0], data, temperature=temperature,
                                    direction=direction[i], colour=colour[i],
                                    main=main[i], label=direction[i])
    tp.plot.mfp.add_cum_kappa(ax[1], data, temperature=temperature,
                              direction=direction[i], colour=colour[i],
                              xmarkers=1e-6, main=main[i], label=direction[i])
add_legend(title='Direction')

# Save

plt.savefig('cumkappa.pdf')
