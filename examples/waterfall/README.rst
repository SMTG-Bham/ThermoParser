.. image:: waterfall.png
   :alt: Waterfall plot of mean free path against frequency with lattice thermal conductivity projected.

This shows a plot of mean free path against frequency with lattice
thermal conductivity projected on the colour axis, as well as a density
of states (DoS).This example highlights a complication: in order for the
waterfall to be on top of the DoS, so as not to obscure the colour, an
additional command, ``format_waterfall``, must be used so the DoS can be
scaled correctly (line 43).
