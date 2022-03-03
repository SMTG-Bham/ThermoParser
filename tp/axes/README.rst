Please see the `examples`_.

.. _examples: https://github.com/smtg-ucl/ThermoPlotter/tree/master/examples

Pre-sized axes for plotting convenience. Many users may want to ignore
this module.

The names are descriptive and split into modules by size. ``small`` has
axes sized for papers and the like, and roughly follows the Nature
guidelines, while ``large`` axes are higher quality and have relatively
much larger text, which work well for presentations and posters but can
take much longer to render.

Each function returns a figure, a set of axes or a list of sets of axes
if there are more than one, and a function to add a legend. The
function takes all the normal ax.legend options, excluding ``loc`` and
``bbox_to_anchor``, as well as ``location``, which determines the
legend position (e.g. ``above``), and if there are multiple axes,
``custom``, which allows ``handles`` and ``labels`` to be set manually.
Arrays of axes are in the shape of the axes, i.e. ``four_square``
returns a 2x2 array of axes.

``legend`` contains helper functions for the other modules, but may
also be useful to people using their own axes layouts.
``legend.consolidate``, for example, takes the legend entries from a
list of axes and combines them, removing duplicates.