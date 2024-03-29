----------------
``calculate.py``
----------------

This module contains functions for calculating some properties used in
ThermoParser.

---------------
``settings.py``
---------------

``settings.py`` contains the default style sheet, tick locators
and axis labels, as well as a means to convert the units on loading 
It can load in custom defaults from ``~/.config/tprc.yaml``, to enable
you to plot in your own style.

--------
``axes``
--------

Contains pre-sized axes. The ``small`` roughly follow Nature
guidelines, while the ``large`` axes look better for posters and
presentations.

--------
``data``
--------

Contains modules that run, load, parse and save data.

--------
``plot``
--------

Contains modules that plot data, as well as colourmap generators and
functions which aid or enhance plotting. Plotting functions can read in
defaults from ``~/.config/tprc.yaml``.

---------
``setup``
---------

Functions to aid with setting up calculations. Currently geared towards
electronic property calculations using VASP.

-------
``cli``
-------

Provides a command-line interface.