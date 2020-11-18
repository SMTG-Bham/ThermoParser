## `aniso.py` and `resolve.py`

`aniso.py` module contains functions for selecting or averaging
directions in anisotropic data.
`resolve.py` selects which `aniso.py` function to apply to which
quantities, and also resolves by temperature where appropriate.

## `load.py`

Loads data from standard output files from the various codes into
dictionary form.
It also applys a consisitent naming convention, data structure and unit
convention across codes, with help of `tp.settings`; and attaches
metadata including the data source and units.
`tp.load.phono3py` includes derived quantities (lifetime, mean free path
and occupation), which can be written to hdf5.

## `run.py`

Currently contains a BoltzTraP runner. Temperamental.

## `save.py`

Currently contains a function to save nested dictionaries to hdf5.
