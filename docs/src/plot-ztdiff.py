#!/usr/bin/env python3

from os import path
import tp

bfiles = ['../data/zno/boltztrap.hdf5',
          '../data/basno3/boltztrap.hdf5']
kfiles = ['../data/zno/kappa-m404021.hdf5',
          '../data/basno3/kappa-m363636.hdf5']
label = ['ZnO', 'BaSnO$_3$']
for f in kfiles:
    if not path.isfile(f) or (path.getsize(f) < 1024*1024*100):
        raise Exception('File not found, please use get-data.sh in the folder above.')

# Axes
fig, ax, add_legend = tp.axes.small.one_colourbar()

# Load
adata = [tp.data.load.boltztrap(b) for b in bfiles]
kdata = [tp.data.load.phono3py(k) for k in kfiles]

# Add
_, h, l = tp.plot.heatmap.add_ztdiff(ax, *adata, *kdata, label1=label[0],
                                     label2=label[1])

add_legend(handles=h, labels=l)

# Save
fig.savefig('ztdiff.pdf')
fig.savefig('ztdiff.png')
