-------------
Average Rates
-------------

.. image:: avg-rates.png
   :alt: A plot of weighted average scattering rates against temperature and carrier concnetration.

This shows the average scattering rates against temperature, averaged
across k-points and weighted by the derivative of the Fermi-Dirac
distribution, as they are when calculating the conductivity, thereby
giving a representative image of the effect the scattering processes
play in the material. This can be plotted at the command line with:

.. code-block::

   tp plot avg-rates ../data/basno3/mesh_75x75x75.h5 -t 1000 -n -1e19 --location 2 --large -c red -c blue -c magenta

If rates are included in the file but are very low, you may want to
``--exclude`` them; ``--exclude PIE`` may be particularly popular.
While ThermoParser does not currently have python functions to plot
line graphs, it does calculate the weighted rates and the ancillary
functions should make this relatively straightforward:

.. literalinclude:: plot-avg-rates.py
   :language: python
   :linenos:
   :emphasize-lines: 7,24,26,27

The weighted rates can be calculated directly via the
``tp.data.load.amset_mesh`` function, by specifying ``weighted_rates``
as the ``quantity`` (lines 7, 24). Selecting which data to show is then
as simple as using the resolve function (lines 26/27). To plot the
graph matplotlib's normal ``plot`` function is used, and the axes,
labels, ticks and legend are all sorted out by tp functions (lines 21,
34-38. 39/40 and 42 respectively, see also `Tutorial-05`_).

.. _Tutorial-05: https://smtg-ucl.github.io/ThermoParser/tutorial-05.html

This graph can instead show weighted mean free paths, by adding the
``--mfp`` tag at the command line or changing ``q`` to equal
``weighted_mfp`` in python.
