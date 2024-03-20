-------
Phonons
-------

.. image:: ../../examples/phonons/phonons.png
   :alt: Phonon dispersion and density of states.

This shows a phonon dispersion and density of states (DoS). This is the
only combined plot currenty possible at the command line, with:

.. literalinclude:: ../../examples/phonons/plot-phonons.sh
   :language: bash

and in python:

.. literalinclude:: ../../examples/phonons/plot-phonons.py
   :language: python
   :linenos:
   :emphasize-lines: 22,24

All plot-types in ``tp.plot.frequency`` have an invert argument to plot
them side-on by a phonon dispersion, which also shortens their x-axis
labels and removes their y-axis labels and tick labels (line 22).
Currently the axes of the two plots have to be aligned manually, but
this is simple (line 24).

.. image:: ../../examples/phonons/multiphon.png
   :alt: Phonon dispersions for different supercell sizes.

This shows phonons dispersions for various supercell sizes for
convergence:

.. literalinclude:: ../../examples/phonons/plot-multiphon.sh
   :language: bash

and in python:

.. literalinclude:: ../../examples/phonons/plot-multiphon.py
   :language: python
   :linenos:
   :emphasize-lines: 17,18

At the command line, the same command is used for both, but in python
a separate command, ``tp.plot.phonons.add_multi`` is used. While this
could be done by looping over ``tp.plot.phonons.add_dispersion``, it
is bundled into a one-liner for convenience. As phonon dispersions do
not always have the same x-scale (at least in phonopy), this rescales
the dispersions so they always match, which is also useful to compare
different materials (so long as they are at least in the same space
group and are preferably closely related), or materials under expansion
or compression. This also demonstrates the large axes style, more
appropriate for presentations or posters than the default style, which
is better for papers; and the add_legend function, which adds a
pre-positioned legend, and accepts all the other usual arguments
including ``title`` (line 18).
