.. image:: figures/tp-logo-header.png
    :target: https://smtg-ucl.github.io/ThermoPlotter/
    :align: center

.. code-block::

     ________
    ///// \\\\
    \________/_______________________________________________________________
    |_____                            : ___                                  \
    | |   |                           :|   \ |           |     |             /
    | |   |__   __  |___  |_____   __ :|___/ |    ___  __|__ __|__  __  |__  \
    | |   |  | /  \ |   \ |  |  | /  \:|     |   /   \   |     |   /  \ |  \ /
    | |   |  | |__/ |     |  |  | |  |:|     |   |   |   |     |   |__/ |    \_____
    | |   |  | \__  |     |  |  | \__/:|     \__ \___/   \__   \__ \__  |    :1.0.0\
    |_________________________________:______________________________________:_____/


.. image:: https://travis-ci.com/SMTG-UCL/ThermoPlotter.svg?branch=master
    :target: https://travis-ci.com/SMTG-UCL/ThermoPlotter

ThermoPlotter is a toolkit used to simplify the production of
high-quality plots from the outputs of specialised analytical codes. It
is focused on computational materials science and particularly
thermoelectrics materials. Traditionally, the steps required to
transform raw data, produce appropriate plots and adjust their
appearance are arduous and often result in long, unweildy python
scripts. ThermoPlotter is built on top of `matplotlib`_ and greatly
simplifies this process. It can be used to build short,
easy-to-customise plotting scripts and there are also several basic
command-line interface options.

Click on the image to go to the `gallery`_!

.. image:: figures/wideband.png
   :alt: A phonon dispersion where widened bands show phonon scattering
   :target: https://smtg-ucl.github.io/ThermoPlotter/gallery.html

.. _gallery: https://smtg-ucl.github.io/ThermoPlotter/gallery.html

Installation
------------

ThermoPlotter can easily be installed with git and pip:

.. code-block:: bash

    git clone git@github.com:SMTG-UCL/ThermoPlotter.git
    cd ThermoPlotter
    python3 -m pip install --user .

After installing, you may want to copy ``ThermoPlotter/tprc.yaml`` to
``~/.config/tprc.yaml``, if you want to set your own default axis
labels, unit conversions, default style sheets (two are provided),
other aesthetic alterations and more!

Usage
-----

ThermoPlotter is designed to have four main stages:

#. *Axes*:
     Pick an axis layout from ``tp.axes``.
#. *Load*:
     Use the functions is ``tp.data.load`` to load the relevant data.
#. *Add*:
     Use functions in ``tp.plot`` modules to add graphs to the axes.
#. *Save*:
     Use ``plt.savefig`` or equivalent to produce the figure.

As ThermoPlotter is dependent on matplotlib, each stage can be
substituted with bespoke code, e.g. using ``matplotlib.pyplot.subplots``
or ``matplotlib.axes.Axes.scatter``. These can still be
supplemented with ThermoPlotter helper functions, such as default labels 
which the user can set in ``tp.settings`` 
or colourmap generators in ``tp.plot.colour``.

The best way to get a feel for ThermoPlotter is to see it in action:
Take a look at our  `examples <https://github.com/smtg-ucl/ThermoPlotter/tree/master/examples>`_ scripts.

Currently supported codes are:

* Phononic properties:

  * `Phonopy <https://phonopy.github.io/phonopy/>`_
  * `Phono3py <http://phonopy.github.io/phono3py/>`_

* Electronic properties:

  * `AMSET <https://hackingmaterials.lbl.gov/amset/>`_
  * `BoltzTraP <https://www.imc.tuwien.ac.at/forschungsbereich_theoretische_chemie/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap/>`_

Current plotting modes are split into four areas.

* ``tp.plot.phonons`` contains plots along a high-symmetry path,
  including phonon dispersions and plots which project other quantities
  onto these paths in various ways.
* ``tp.plot.frequency`` plots frequency on the x-axis, including density
  of states (DoS), cumulative kappa, "waterfall" and density plots.
  Each function has a ``main`` argument, which can be useful when
  plotting multiple quantities on the same set of axes; and an
  ``invert`` argument, which swaps the x and y axes to let you plot
  DoS-style next to a ``tp.plot.phonons`` plot.
* ``tp.plot.mfp`` contains a cumulative kappa against mean free path
  plot.
* ``tp.plot.heatmap`` contains a heatmap plotter, and wrappers which
  format appropriately for ZT against temperature and doping
  concentration; and one which plots the lattice thermal conductivity
  required to reach a target ZT, again against temperature and doping.

A set of example scripts is provided in the ``tp/examples`` folder, and
there is `documentation`_.

Contributing
------------

We welcome any contributions, whether they be a feature request or a new
piece of code (or anything else). Adding options is inteded to be
straightforward; the modularity of the code means that each step is mostly 
independent of the others. 

Bugs and feature requests can be submitted to the `issue tracker`_,
while contributions can be made using the `fork and pull`_ approach.
Contributions should include comprehensive docstrings, and where
appropriate `examples`_, further `documentation`_ and `tests`_ are greatly
appreciated. Documentation uses the `sphinx`_ package, and can be built from the docs
directory with ``sphinx-build -b html src/ .``.

.. _issue tracker: https://github.com/smtg-ucl/ThermoPlotter/issues
.. _fork and pull: https://guides.github.com/activities/forking
.. _examples: https://github.com/smtg-ucl/ThermoPlotter/tree/master/examples
.. _documentation: https://smtg-ucl.github.io/ThermoPlotter/
.. _tests: https://github.com/smtg-ucl/ThermoPlotter/tree/master/tests

Testing
-------

Tests use the `unittest`_ package, and can be run from the test directory
with ``python3 -m unittest``.

.. _unittest: https://docs.python.org/3/library/unittest.html

Contributors
------------

Many thanks to all those who contributed code or ideas to ThermoPlotter!
Roughly chronologically, they are so far:

* Kieran B. Spooner
* Maud Einhorn
* David O. Scanlon
* Daniel W. Davies
* Bonan Zhu
* Sean R. Kavanagh
* Warda Rahim

License
-------

ThermoPlotter is licensed under the GNU Affero General Public License v3
(AGPLv3).

Requirements
------------

ThermoPlotter uses the following open-source packages:

* `h5py <http://docs.h5py.org/>`_
* `json <https://docs.python.org/3/library/json.html>`_
* `matplotlib <https://matplotlib.org>`_
* `numpy <https://numpy.org>`_
* `pymatgen <https://pymatgen.org>`_
* `scipy <https://www.scipy.org>`_
* `sphinx <https://www.sphinx-doc.org>`_
* `yaml <https://pyyaml.org/>`_
