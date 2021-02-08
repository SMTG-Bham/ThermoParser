Please see the `examples`_.

.. _examples: https://github.com/smtg-ucl/ThermoPlotter/tree/master/examples

All plotting functions can read in defaults from ``~/.config/tprc.yaml``
to enable consistent customisation. A template ``tprc.yaml`` is
available in the top directory.

-------------
``colour.py``
-------------

Contains colourmap generators:

* ``linear``: a linear colourmap between a max and min colour.
* ``uniform`` and ``elbow``: colourmaps where the highlight colour is
  part way through the map. ``uniform`` calculates the midpoint relative
  to the end points for a more uniform plot, while in ``elbow`` the
  position can be chosen.
* ``highlight``: adds a highlight colour to an existing colourmap.
  may be useful in ``tp.plot.frequency.add_waterfall`` or
  ``tp.plot.phonons``.
* ``skelton``: Jonathan Skelton's rainbowy colourmap from his waterfall
  plots. He has lots of useful stuff in `Phono3py-Power-Tools`_.

.. _Phono3py-Power-Tools: https://github.com/skelton-group/Phono3py-Power-Tools

----------------
``frequency.py``
----------------

Plots frequency on the x-axis (mostly).
Written in such a way as to allow plotting several quantities on the
same axes, as information can often be gleaned from the interaction of
several of these plots.

.. tip::

    Use the ``scale`` argument in ``add_dos`` and ``add_cum_kappa`` to
    scale them to the axes scale and limits when plotting with other
    quantities.

.. tip::

    Use the ``invert`` argument to invert the axes to plot alongside a
    phonon dispersion from ``phonons.py``.

* ``add_dos``: adds a phonon density of states (DoS).
* ``add_cum_kappa``: adds a cumulative lattice thermal conductivity vs
  frequency plot.
* ``add_waterfall``: adds a scatter plot of a third-order phonon
  property against frequency, broken down per mode per q-point. Also
  supports alternate x-axis quantities.
* ``add_projected_waterfall``: like ``add_waterfall``, but with a second
  (or third) third-order phonon quantity projected on the colour axis.
* ``format_waterfall``: formats the axes for waterfall plots.

.. tip::
    If you run ``format_waterfall`` then ``add_dos`` and then
    ``add_waterfall``, the waterfall will be overlayed on top of the DoS
    so the waterfall colours remain clear.

--------------
``heatmap.py``
--------------

Plots heatmaps.

* ``add_heatmap``: The base function, which wraps around
  ``matplotlib.pyplot.pcolormesh``, but with extra functionality
  including avoiding data being accidentally cropped, setting the axes
  limits within the function so colourmaps can be renormalised,
  automated colourbar extension and integration of ``tp.plot.colours``
  colourmaps so a single or pair of ``#RRGGBB`` colours can be used as
  the colour argument.
* ``add_ztmap``: Reads ZT data, or reads raw electronic transport data
  and lattice thermal conductivity data and calculates ZT (in which case
  the data is written to ``hdf5``), and plots against temperature and
  carrier concentration.
* ``add_kappa_target``: reads electronic transport data and plots the
  lattice thermal conductivity required to reach a specified ZT.
  Envisioned as a screening step before expensive phonon calculations.

--------------
``phonons.py``
--------------

Plots phonon dispersions and variations thereof. For those that require
third-order phonon quantities, the keyword ``dispersion`` can be used
instead of ``qpoint`` in ``tp.data.load``. Some functions are there for
completeness and should probably not be used. Heed Uncle Ben.

* ``add_dispersion``: The classic. Phonon frequencies along a
  high-symmetry path.
* ``add_multi``: Loops over ``add_dispersion``. Great for convergence,
  rescales x-data to match each other.
* ``add_alt_dispersion``: Replaces frequency with a third-order phonon
  quantity.
* ``add_projected_dispersion``: projects a third-order phonon quantity
  onto the colour axis of a standard dispersion.
* ``add_alt_projected_dispersion``: Third-order phonon quantities on
  the y and colour axes.
* ``add_wideband``: Frequencies are lorentzian-broadened with imaginary
  self-energy to visualise scattering.\ :sup:`1` Has its own keyword
  in ``tp.data.load``: ``wideband``.

----------------
``utilities.py``
----------------

Functions which assist in plotting.

* ``colour_scale``: Attempts to optimise the colour scale for
  third-order phonon properties. Optimised for larger q-point meshes.
* ``scale_to_axis``: Scales data to an axis.
* ``set_locators``: A one-liner for setting up locators for both axes,
  with a ``dos`` argument to get rid of labels from an inverted
  ``tp.plot.frequency`` plot alongside a ``tp.plot.phonons`` plot.

---------
Reference
---------

.. code-block:: bibtex

@article{Maradudin1962,
  title={Scattering of neutrons by an anharmonic crystal},
  author={Maradudin, AA and Fein, AE},
  journal={Physical Review},
  volume={128},
  number={6},
  pages={2589},
  year={1962},
  publisher={APS}
}
