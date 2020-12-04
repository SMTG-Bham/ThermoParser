-------------------------------
``aniso.py`` and ``resolve.py``
-------------------------------

``aniso.py`` module contains functions for selecting or averaging
directions in anisotropic data. ``resolve.py`` selects which
``aniso.py`` function to apply to which quantities, and also resolves by
temperature where appropriate.

-----------
``load.py``
-----------

Loads data from standard output files from the various codes into
dictionary form. It also applys a consisitent naming convention, data
structure and unit convention across codes, with help of
``tp.settings``; and attaches metadata including the data source and
units.

* `AMSET`_

    * ``tp.load.amset`` loads transport properties from ``json``
    * ``tp.load.amset_mesh`` loads scattering properties from ``hdf5``

* BoltzTraP:

    * ``tp.load.boltztrap`` reads the ``hdf5`` outputted by
      ``tp.run.boltztrap``, although this can show questionable
      reliability

* `Phono3py`_

    * ``tp.load.phono3py`` includes derived quantities (lifetime, mean
      free path and occupation), which can be written to ``hdf5``

* `Phonopy`_

    * ``tp.load.phonopy_dispersion`` reads both Phonopy and `sumo`_
      ``disp.yaml``\ s
    * ``tp.load.phonopy_dos`` reads ``projected_dos.dat``

.. _AMSET: <https://hackingmaterials.lbl.gov/amset/>
.. _Phono3py: <https://phonopy.github.io/phono3py/hdf5_howto.html#kappa-hdf5-file>
.. _Phonopy: <https://phonopy.github.io/phonopy/>
.. _sumo: <https://github.com/SMTG-UCL/sumo>

----------
``run.py``
----------

Currently contains a BoltzTraP runner. Temperamental.

-----------
``save.py``
-----------

Currently contains a function to save nested dictionaries to ``hdf5``.
