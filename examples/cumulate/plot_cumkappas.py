#!/usr/bin/env python3

import matplotlib.pyplot as plt
import tp

kappafile = '../data/kappa-m505028.hdf5'
temperature = 300
direction = ['x', 'z']
quantities = 'frequency mode_kappa mfp'

main = [True, False]
colour = ['#ff8000', '#0000ff']

# Axes

fig, ax = tp.plot.axes.two_h_small_legend()

# Load

data = tp.data.load.phono3py(kappafile, quantities=quantities)

# Add

for i in [0, 1]:
    tp.plot.frequency.add_cum_kappa(ax[0], data, temperature=temperature,
                                    direction=direction[i], colour=colour[i],
                                    main=main[i], **{'label': direction[i]})
    tp.plot.mfp.add_cum_kappa(ax[1], data, temperature=temperature,
                              direction=direction[i], colour=colour[i],
                              xmarkers=[2e-8], main=main[i],
                              **{'label': direction[i]})

# Save

ax[1].legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Direction')
plt.savefig('cumkappas.pdf')
