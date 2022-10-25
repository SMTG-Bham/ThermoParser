#!/bin/bash

f1='../data/basno3/kappa-m363636.hdf5'
f2='../data/zno/kappa-m404021.hdf5'
repo1='https://zenodo.org/record/4428179/files/kappa-m363636.hdf5?download=1'
repo2='https://zenodo.org/record/4428179/files/kappa-m404021.hdf5?download=1'

if [[ ! -f $f1 ]]; then
    wget $repo1 -O $f1
else
    echo Already downloaded!
fi

if [[ ! -f $f2 ]]; then
    wget $repo2 -O $f2
else
    echo Already downloaded!
fi
