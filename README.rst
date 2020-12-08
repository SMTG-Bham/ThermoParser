.. code-block::

     ________
    ///// \\\\
    \________/_______________________________________________________________
    |_____                            : ___                                  \
    | |   |                           :|   \ |           |     |             /
    | |   |__   __  |___  |_____   __ :|___/ |    ___  __|__ __|__  __  |__  \
    | |   |  | /  \ |   \ |  |  | /  \:|     |   /   \   |     |   /  \ |  \ /
    | |   |  | |__/ |     |  |  | |  |:|     |   |   |   |     |   |__/ |    \_____
    | |   |  | \__  |     |  |  | \__/:|     \__ \___/   \__   \__ \__  |    :0.3.1\
    |_________________________________:______________________________________:_____/


ThermoPlotter is a toolkit for quickly, easily and prettily plotting the
outputs of specialised analytical codes. It is focused on computational
materials science and particularly thermoelectrics materials. It
essentially wraps around `matplotlib`_ functions, and is intended to be
used as a python package, to build easily customisable plotting scripts.
There are also several basic command-line scripts.

------------
Installation
------------

ThermoPlotter can easily be installed with git and pip:

.. code-block:: bash

    git clone https://github.com/kbspooner/ThermoPlotter
    cd ThermoPlotter
    python3 -m pip install --user -e .

.. tip::

    The ``-e`` or editable option is encouraged so you can add your own
    defaults in ``tp.settings`` and elsewhere.

As it is designed to be customisable, when upgrading you may want to
preserve your changes:

.. code-block:: bash

    git stash
    git pull
    python3 -m pip install --user -e .
    git stash apply

.. tip::

    You only need to repeat the pip stage if you want to update the command-
    line script list.

-----
Usage
-----

ThermoPlotter is designed to have four main stages:

#. *Axes*:
   Pick an axis layout from ``tp.axes``.
#. *Load*:
     Use the functions is ``tp.data.load`` to load the relevant data.
#. *Add*:
     Use functions in modules in ``tp.plot`` to add graphs to the axes.
#. *Save*:
     Use ``plt.savefig`` or equivalent to produce the figure.

As ThermoPlotter is simply a scripting code, each option can be
substituted with bespoke code (i.e. using ``matplotlib.pyplot.subplots``
or ``matplotlib.axes.Axes.scatter``); but these can still be
supplemented with helper functions, such as default labels in
``tp.settings`` or colourmap generators in ``tp.plot.colour``.

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
  of states (DoS), cumulative kappa and "waterfall" plots.
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

A full set of example scripts is provided in the ``tp/examples`` folder.

------------
Contributing
------------

We welcome any contributions, whether they be a feature request or a new
piece of code (or anything else). Adding options is inteded to be
straightforward, as each step is mostly independant of the others, so
only one new function should be required. We would of course be happy to
discuss, if desired.

Bugs and feature requests can be submitted to the `issue tracker`_,
while contributions can be made using the `fork and pull`_ approach.
Contributions should have comprehensive docstrings, and where
appropriate `examples`_, documentation and `tests`_ are greatly
appreciated.

.. _issue tracker: https://github.com/kbspooner/ThermoPlotter/issues
.. _fork and pull: https://guides.github.com/activities/forking
.. _examples: https://github.com/kbspooner/ThermoPlotter/tree/master/examples
.. _tests: https://github.com/kbspooner/ThermoPlotter/tree/master/tests

-------
License
-------

ThermoPlotter is licensed under the GNU Affero General Public License v3
(AGPLv3).

------------
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
