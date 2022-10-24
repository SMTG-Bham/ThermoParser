Data handling modules.

----------------
``utilities.py``
----------------

``utilities.py`` uses the metadata stored by ThermoPlotter to
manipulate data.

* ``tp.data.utilities.resolve`` breaks down data by dependent variable,
  such as by temperature or direction.
* ``tp.data.utilities.merge`` merges data dictionaries along a common
  axis. This is particularly useful for codes that have heavy memory
  requirements, such as AMSET, so you can run them on the first half of
  temperatures and the second, and then stitch them together here to
  get denser data along a large temperature range. Only merges along
  one direction at a time, so a two-dimensional merge would require
  multiple merges: one along the horizontal axis for each set of data
  along the vertical axis, and then a final merge along the vertical
  axis, or vice versa.

-----------
``load.py``
-----------

Loads data from standard output files from the various codes into
dictionary form. It also applys a consisitent naming convention, data
structure and unit convention across codes, with help of
``tp.settings``; and attaches metadata including the data source and
units.

* `AMSET`_

  * ``tp.data.load.amset`` loads transport properties from ``json``.
  * ``tp.data.load.amset_mesh`` loads scattering properties from ``hdf5``,
    including derived properties such as weighted rates.

* BoltzTraP:

  * ``tp.data.load.boltztrap`` reads the ``hdf5`` outputted by
    ``tp.run.boltztrap``.

* `Phono3py`_

  * ``tp.data.load.phono3py`` includes derived quantities (lifetime, mean
    free path and occupation), which can be written to ``hdf5``.

* `Phonopy`_

  * ``tp.data.load.phonopy_dispersion`` reads both Phonopy and `sumo`_
    ``disp.yaml`` files (we recommend the latter).
  * ``tp.data.load.phonopy_dos`` reads ``projected_dos.dat``

.. _AMSET: <https://hackingmaterials.lbl.gov/amset/>
.. _Phono3py: <https://phonopy.github.io/phono3py/hdf5_howto.html#kappa-hdf5-file>
.. _Phonopy: <https://phonopy.github.io/phonopy/>
.. _sumo: <https://github.com/SMTG-UCL/sumo>

----------
``run.py``
----------

Currently contains a BoltzTraP runner. Faster and less error-prone than
the Pymatgen version, and writes to hdf5 for standardisation purposes.

-----------
``save.py``
-----------

Saves data to files.

* ``tp.data.save.phono3py``: Save calculated properties to hdf5.
* ``tp.data.save.zt``: Save ZT to hdf5 and highlights to yaml. CLI:
  ``tp save zt``.
* ``tp.data.save.kappa_target``: Save maximimum lattice thermal
  conductivity required to reach a target ZT to hdf5.
* ``tp.data.save.cumkappa``: Save cumulative lattice thermal
  conductivity to plain text. CLI: ``tp save cumkappa``.
* ``tp .data.save.hdf5``: Saves nested dictionaries to hdf5.

