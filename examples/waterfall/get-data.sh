#!/bin/bash

f='../data/basno3/kappa-m363636.hdf5'
repo='https://zenodo.org/record/4428179/files/kappa-m363636.hdf5?download=1'

if [[ ! -f $f ]]; then
    wget $repo -O $f
else
    echo Already downloaded!
fi
