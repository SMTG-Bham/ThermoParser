#!/usr/bin/env python3

from os import path
import tp

efiles = ['../data/basno3/transport_75x75x75.json',
          '../data/basno3/boltztrap.hdf5']
kfile = '../data/basno3/kappa-m363636.hdf5'
label = ['MRTA', 'CRTA']
if not path.isfile(kfile) or (path.getsize(kfile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh in the folder above.')

# Axes
fig, ax, add_legend = tp.axes.small.one_colourbar()

# Load
adata = [tp.data.load.amset(efiles[0]),
         tp.data.load.boltztrap(efiles[1])]
kdata = tp.data.load.phono3py(kfile)

# Add
_, h, l = tp.plot.heatmap.add_ztdiff(ax, *adata, kdata, kdata, label1=label[0],
                                     label2=label[1])

add_legend(handles=h, labels=l)

# Save
fig.savefig('ztdiff.pdf')
fig.savefig('ztdiff.png')
