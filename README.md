```
 ________
///// \\\\
\________/_______________________________________________________________
|_____                            : ___                                  \
| |   |                           :|   \ |           |     |             /
| |   |__   __  |___  |_____   __ :|___/ |    ___  __|__ __|__  __  |__  \
| |   |  | /  \ |   \ |  |  | /  \:|     |   /   \   |     |   /  \ |  \ /
| |   |  | |__/ |     |  |  | |  |:|     |   |   |   |     |   |__/ |    \_____
| |   |  | \__  |     |  |  | \__/:|     \__ \___/   \__   \__ \__  |    :0.3.0\
|_________________________________:______________________________________:_____/
```

ThermoPlotter is a toolkit for quickly, easily and prettily plotting the
outputs of specialised analytical codes.
It is focused on computational materials science and particularly
thermoelectrics materials.
It essentially wraps around [Matplotlib](https://matplotlib.org/)
functions, and is intended to be used as a python package, to build
easily customisable plotting scripts.
There are also several basic command-line scripts.

## Installation

ThermoPlotter can easily be installed with git and pip:

```bash
git clone https://github.com/kbspooner/ThermoPlotter
cd ThermoPlotter
python3 -m pip install --user -e .
```

As it is designed to be customisable, when upgrading you may want to
preserve your changes:

```bash
git stash
git pull
python3 -m pip install --user -e .
git stash apply
```

You only need to repeat the pip stage if you want to update the command-
line script list.

## Usage

ThermoPlotter is designed to have four main stages:

  1. *Axes*:
     Pick an axis layout from `tp.axes`.
  2. *Load*:
     Use the functions is `tp.data.load` to load the relevant data.
  3. *Add*:
     Use functions in modules in `tp.plot` to add graphs to the axes.
  4. *Save*:
     Use `plt.savefig` or equivalent to produce the figure.

As ThermoPlotter is simply a scripting code, each option can be
substituted with bespoke code (i.e. using `plt.subplots` or 
`ax.scatter`); but these can still be supplemented with helper
functions, such as default labels in `tp.settings` or colourmap
generators in `tp.plot.colour`.

Currently supported codes are:

Phononic properties:
* [Phonopy](https://phonopy.github.io/phonopy/)
* [Phono3py](http://phonopy.github.io/phono3py/)

Electronic properties:
* [AMSET](https://hackingmaterials.lbl.gov/amset/)
* [BoltzTraP](https://www.imc.tuwien.ac.at/forschungsbereich_theoretische_chemie/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap/)

Current plotting modes are split into four areas.

* `tp.plot.phonons` contains plots along a high-symmetry path, including
phonon dispersions and plots which project other quantities onto these
paths in various ways.
* `tp.plot.frequency` plots frequency on the x-axis, including density
of states (DoS), cumulative kappa and "waterfall" plots.
Each function has a `main` argument, which can be useful when plotting
multiple quantities on the same set of axes; and an `invert` argument,
which swaps the x and y axes to let you plot DoS-style next to a
`tp.plot.phonons` plot.
* `tp.plot.mfp` contains a cumulative kappa against mean free path plot.
* `tp.plot.heatmap` contains a heatmap plotter, and wrappers which
format appropriately for ZT against temperature and doping
concentration; and one which plots the lattice thermal conductivity
required to reach a target ZT, again against temperature and doping.

A full set of example scripts is provided in the `tp/examples` folder.

## Contributing

We welcome any contributions, whether they be a feature request or a new
piece of code (or anything else).
Adding options is inteded to be straightforward, for example adding a
new axes layout should be as simple as copying another function in
`tp.axes` and altering the numbers; and adding support for a new code
should only require adding a function in the `tp.data.load` module.
We would of course be happy to discuss, if desired.

Bugs and feature requests can be submitted to the [issue tracker](https://github.com/kbspooner/ThermoPlotter/issues),
while contributions can be made using the [fork and pull](https://guides.github.com/activities/forking)
approach.
Contributions should have comprehensive docstrings, and where
appropriate [examples](https://github.com/kbspooner/ThermoPlotter/tree/master/examples),
documentation and [tests](https://github.com/kbspooner/ThermoPlotter/tree/master/tests)
are greatly appreciated.

## License

ThermoPlotter is licensed under the GNU Affero General Public License v3
(AGPLv3).

## Requirements

ThermoPlotter requires the following open-source packages:

* [h5py](http://docs.h5py.org/en/stable/)
* [json](https://docs.python.org/3/library/json.html)
* [matplotlib](https://matplotlib.org/)
* [numpy](https://numpy.org/)
* [pymatgen](https://pymatgen.org/)
* [scipy](https://www.scipy.org/)
* [yaml](https://pyyaml.org/wiki/PyYAMLDocumentation)
