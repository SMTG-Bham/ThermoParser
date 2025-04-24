"""Plotting tools.

The main tools are called ``add_`` something, which take a set of axes
and data dictionaries as inputs, amoung other things; these are
accompanied by ancillary scripts including formatting scripts and a
module of colourmap generators.
"""

#Modules
#-------
#
#    colour
#        colourmap generators and colour converters.
#    frequency
#        plots frequency on the x-axis.
#    heatmap
#        generic heatmap plotter and specific wrappers.
#    mfp
#        cumulative lattice thermal conductivitity vs. mean free path.
#    phonons
#        plots q-points on the x axis.
#    utilities
#        ancillary scripts.
#"""

from . import colour, frequency, heatmap, mfp, phonons, utilities
