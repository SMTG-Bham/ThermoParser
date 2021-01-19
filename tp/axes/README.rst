Please see the `examples`_.

.. _examples: https://github.com/smtg-ucl/ThermoPlotter/tree/master/examples

Pre-sized axes for plotting convenience. Many users may want to ignore
this module.

The names are descriptive and split into modules by number of primary
sets of axes (DoS axes are not counted in this number). The default axes
roughly follow the Nature guidelines, while the ``_large`` axes look
better in posters and presentations. Where there are multiple sets of
axes (including DoS axes), they are returned in an array of the shape
they appear in, i.e. ``four.square`` returns a figure and a 2x2 array of
axes. Where there is legend space, a function is returned which will
generate and position the legend for you, and still takes the normal
arguments (title etc.).
