-------------------------------
Tutorial-01: Python Foundations
-------------------------------

Here we show a phonon dispersion, calculated via the supercell method
in Phonopy, which shows the vibrational frequencies of ZnO along a
high symmetry path in reciprocal space. The gradient of the bands is
the group velocity of the phonons, an important factor in the lattice
thermal conductivity; and if there are negative-frequency (imaginary)
modes, it indicates kinetic instability. More information can also be
gleaned, which is covered in later tutorials.

.. image:: tutorial-01.png
   :alt: Phonon dispersion of ZnO.

The most basic plots in ThermoParser consist of four commands:

 1. `Axes`_
 2. `Load`_
 3. `Plot`_
 4. Save

ThermoParser offers ways to do the first three (the last is trivial),
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

ThermoParser offers a number of pre-sized axes in the ``axes`` module.
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

ThermoParser contains several data loading functions for different
inputs in the ``data.load`` module. At their most basic, these take a
file to read from, and return a dictionary of the data.

While you can load your own data, ThermoParser takes several steps to
ensure consistent format and units between codes and also with some
old versions of these codes; and with the plotting functions. It also
consistently applies custom unit conversions and provides metadata
including array shapes, units and data sources, so care should be
taken if this step is done manually.

Plot (line 12)
--------------

ThermoParser contains a number of plotting functions in the ``plot``
module. At their most basic, they take a set of axes to plot on, and a
data dictionary to read from. All plot functions also take a number of
customisation option, including all the ``kwargs`` from the underlying
``matplotlib`` function.

This stage can be replaced with ordinary ``matplotlib`` functions such
as ``ax.plot``, which can be assisted by a number of ancillary
ThermoParser functions discussed in `Tutorial-05`_.

Save (line 15)
--------------

Simply ``figure.savefig(name.extension)``.

More Help
---------

ThermoParser has extensive docs. Docstrings are available throughout
the code itself, if you are using IPython or Jupyter and you can access
them by typing the command name followed by two question marks, and
IDEs normally have options too. If you are on the master branch,
they can also be found `here <https://smtg-bham.github.io/ThermoParser/>`_.

.. _Tutorial-03: https://smtg-bham.github.io/ThermoParser/tutorial-03.html
.. _Tutorial-05: https://smtg-bham.github.io/ThermoParser/tutorial-05.html
.. _Axes: https://smtg-bham.github.io/ThermoParser/tp.axes.html
.. _Load: https://smtg-bham.github.io/ThermoParser/tp.data.html#module-tp.data.load
.. _Plot: https://smtg-bham.github.io/ThermoParser/tp.plot.html