Please see the `examples`_.

.. _examples: <https://github.com/kbspooner/ThermoPlotter/tree/master/examples>

Pre-sized axes for plotting convenience. Many users may want to ignore
this module.

The names are descriptive and split into modules by number of primary
sets of axes (DoS axes are not counted in this number). Where there are
multiple sets of axes (including DoS axes), they are returned in an
array of the shape they appear in, i.e. ``four.square`` returns a figure
and a 2x2 array of axes. Where there is legend space, a function is
returned which will position the legend for you, and still takes the
normal title etc. arguments.
