-------------------------------
Tutorial-01: Python Foundations
-------------------------------

.. image:: tutorial-01.png
   :alt: Phonon dispersion of ZnO.

The most basic plots in ThermoPlotter consist of four commands:

 1. `Axes`_
 2. `Load`_
 3. `Plot`_
 4. Save

ThermoPlotter offers ways to do the first three (the last is trivial),
but each part can be replaced by your own code if you prefer.

The CLI version of this code is:

.. code-block:: bash

   tp plot phonons ../data/zno/band.yaml

And the python version is:

.. literalinclude:: tutorial-01.py
   :language: python
   :linenos:
   :emphasize-lines: 6,9,12,15
   
Axes (line 6)
-------------

ThermoPlotter offers a number of pre-sized axes in the ``axes`` module.
Within, there is a ``large`` module, which provides presentation-
oriented figures, but these take relatively long to load so here we use
``small``, which is better for papers. The functions consist of the
number of primary axes, and a description where necessary. In this
case, we use the most basic ``one``.

Every function returns three things:

 1. figure: The whole plot area, used in step 4.
 2. axes: Where the data goes, used in step 3.
 3. add_legend: Adds a legend. We don't use it here, so we've thrown it
    away by assigning it to ``_``. Giving it any name and not using it
    has the same effect. Legends will be discussed in `Tutorial-03`_.

This command can easily be replaced with ``matplotlib.pyplot`` commands
such as ``figure`` and ``subfigs``.

Load (line 9)
-------------

ThermoPlotter contains several data loading functions for different
inputs in the ``data.load`` module. At their most basic, these take a
file to read from, and return a dictionary of the data.

While you can load your own data, ThermoPlotter takes several steps to
ensure consistent format and units between codes and also with some
old versions of these codes; and with the plotting functions. It also
consistently applies custom unit conversions and provides metadata
including array shapes, units and data sources, so care should be
taken if this step is done manually.

Plot (line 12)
--------------

ThermoPlotter contains a number of plotting functions in the ``plot``
module. At their most basic, they take a set of axes to plot on, and a
data dictionary to read from. All plot functions also take a number of
customisation option, including all the ``kwargs`` from the underlying
``matplotlib`` function.

This stage can be replaced with ordinary ``matplotlib`` functions such
as ``ax.plot``, which can be assisted by a number of ancillary
ThermoPlotter functions discussed in `Tutorial-05`_.

Save (line 15)
--------------

Simply ``figure.savefig(name.extension)``.

More Help
---------

ThermoPlotter has extensive docs. Docstrings are available throughout
the code itself, if you are using IPython or Jupyter and you can access
them by typing the command name followed by two question marks, and
IDEs normally have options too. If you are on the master branch,
they can also be found `here <https://smtg-ucl.github.io/ThermoPlotter/>`_.

.. _Tutorial-03: https://smtg-ucl.github.io/ThermoPlotter/tutorial-03.html
.. _Tutorial-05: https://smtg-ucl.github.io/ThermoPlotter/tutorial-05.html
.. _Axes: https://smtg-ucl.github.io/ThermoPlotter/tp.axes.html
.. _Load: https://smtg-ucl.github.io/ThermoPlotter/tp.data.html#module-tp.data.load
.. _Plot: https://smtg-ucl.github.io/ThermoPlotter/tp.plot.html