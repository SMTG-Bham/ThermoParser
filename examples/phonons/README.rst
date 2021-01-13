.. image:: phonons.png
   :alt: Phonon dispersion and density of states.

This shows a phonon dispersion and density of states (DoS). All plot-
types in ``tp.plot.frequency`` have an invert argument to plot them
side-on by a phonon dispersion, which also shortens their x-axis labels
and removes their y-axis labels and tick labels. Currently the axes of
the two plots have to be aligned manually, but this is simple (see
script).

This is also available as a command line script ``tp-phonons``.

This also demonstrates the large axes style, more appropriate for
presentations or posters than the default style, which is better for
papers; and also the add_legend function, which adds a pre-positioned
legend.
