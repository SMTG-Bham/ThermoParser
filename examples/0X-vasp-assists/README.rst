-------------------------
Tutorial-0X: VASP Assists
-------------------------

Some of the calculations required to obtain the data to run e.g. AMSET
and Phono3py can be very expensive. To reduce that burden, there are a
few basic functions to make calculations more efficient, at least in
VASP.

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

    tp gen kpar-k KPOINTS

Target Lattice Thermal Conducitivity
------------------------------------

While not related to VASP, here seems a good point to mention the
kappa-target plot, which shows what lattice thermal conductivity would
be required to achieve a specified ZT. If its too low, you may not want
to bother with the expensive third-order+ phonon caculations!

Merge
-----

Bonus 2: ``tp.data.utilities.merge`` uses the tp metadata to combine
multiple data dictionaries, so one can obtain denser data for memory-
intensive calculations (such as AMSET) by running multiple times and
merging the data dictionaries before plotting.
