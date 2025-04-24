----------------------------------
Tutorial-0X: Miscellaneous Assists
----------------------------------

Some of the calculations required to obtain the data to run e.g. AMSET
and Phono3py can be very expensive. To reduce that burden, there are a
few basic functions to make calculations more efficient. The first is
specific to VASP, but the rest are more general.

Zero-Weighted k-points
----------------------

AMSET and BoltzTraP require dense k-point grids to get accurate results.
Not all k-points are created equal, however, so we have provided a tool
to combine two KPOINTS files, the converged KPOINTS file which must be
weighted, and a second file of less equal k-points which can be zero-
weighted, to increase the population without costing so much as their
bretheren.

.. code-block:: python

    tp.setup.get_kpoint('weighted_KPOINTS', 'unweighted_KPOINTS')

.. code-block::

    tp gen kpoints -k weighted_KPOINTS -z unweighted_KPOINTS

When setting KPAR, unweighted k-points should not be considered. Our
KPAR generator suggests suitable KPAR values, ignoring the zero-
weighted k-points.

.. code-block:: python

    tp.setup.get_kpar('KPOINTS')

.. code-block::

    tp gen kpar -k KPOINTS

Phonopy Config Files
--------------------

If you generate phonopy data from the command line, configuration files are
useful to save time regenerating and record inputs. ThermoParser will generate
these for you based on a ``POSCAR`` and your inputs.

.. code-block:: python

    tp.setup.get_band_conf('supercell size')
    tp.setup.get_dos_conf('supercell size')

.. code-block::

    tp gen band-conf 'supercell size'
    tp gen dos-conf 'supercell size'

The required argument ``'dim'``, i.e. the supercell size, can be a string or an
array (the latter in python only), and accepts 1x1, 3x1, 3x3 and 6x1 arrays,
that is to say ``2``, ``'2 2 2'``, ``'2 0 0  0 2 0  0 0 2'`` and
``'2 2 2 0 0 0'`` all give the same result.

Target Lattice Thermal Conducitivity
------------------------------------

The kappa-target plot shows what lattice thermal conductivity would be
required to achieve a specified ZT. If it's too low, you may not want
to bother with the expensive third-order+ phonon caculations!

Merge
-----

``tp.data.utilities.merge`` uses the tp metadata to combine
multiple data dictionaries, so one can obtain denser data for memory-
intensive calculations (such as AMSET) by running multiple times and
merging the data dictionaries before plotting.
