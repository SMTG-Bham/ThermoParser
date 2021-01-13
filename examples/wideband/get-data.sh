#!/bin/bash

f='../data/zno/kappa-m404021.hdf5'
repo='https://zenodo.org/record/4428179/files/kappa-m404021.hdf5?download=1'

if [[ ! -f $f ]]; then
    wget $repo -O $f
else
    echo Already downloaded!
fi
