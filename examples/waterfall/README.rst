.. image:: waterfall.png
   :alt: Waterfall plot of mean free path against frequency with lattice thermal conductivity projected.

This shows a plot of mean free path against frequency with lattice
thermal conductivity projected on the colour axis, as well as a density
of states (DoS). It shows the connections between phonon scattering and the
constituent elements, for example there seems to be a big drop in phonon mean
free path in the frequencies with high amounts of Ba character.
This example highlights a complication: in order for the waterfall to be on top
of the DoS, so as not to obscure the colour, an additional command,
``format_waterfall``, must be used so the DoS can be scaled correctly (line 43).
