#!/usr/bin/env python3

from os import path
import tp

phile = '../data/zno/band.yaml'
kappafile = '../data/zno/kappa-m404021.hdf5'
poscar = '../data/zno/POSCAR'
temperature = 300
if not path.isfile(kappafile) or (path.getsize(kappafile) < 1024*1024*100):
    raise Exception('File not found, please use get-data.sh in the folder above.')

colour = ['#000000', '#ff0000']

# running this section as a function is particularly important for mac users
# due to its use of multiprocessing (this is the case for all projected phonon
# plots including alt_phonons and wideband)
def main():
    # Axes
    fig, ax, add_legend = tp.axes.large.one('dark_background')
    
    # Load
    data = tp.data.load.phono3py(kappafile, quantities='wideband')
    pdata = tp.data.load.phonopy_dispersion(phile)
    
    # Add
    tp.plot.phonons.add_wideband(ax, data, pdata, temperature=temperature,
                                 colour=colour, poscar=poscar)
    
    # Save
    fig.savefig('wideband.pdf')
    fig.savefig('wideband.png')

if __name__ == "__main__":
    main()
