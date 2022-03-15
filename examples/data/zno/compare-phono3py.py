#!/usr/bin/env python3

import tp

pfile = 'band.yaml'
kfile = 'kappa-m404021.hdf5'
poscar = 'POSCAR'

pplabel = 'Phonopy'
p3plabel = 'Phono3py'

fig, ax, add_legend = tp.axes.small.one()

pp = tp.data.load.phonopy_dispersion(pfile)
p3p = tp.data.load.phono3py(kfile, 'frequency')

tp.plot.phonons.add_dispersion(ax, pp, colour='grey', label=pplabel)
tp.plot.phonons.add_alt_dispersion(ax, p3p, pp, 'frequency', poscar=poscar,
                                   label=p3plabel, colour='purple')

add_legend(title='ZnO')
fig.savefig('compare-phono3py.pdf')
