Command-line scripts for rapid analysis.

ThermoPlotter uses `click`_ to organise the command line interface (cli).
This allows a more navigable, flexible and maintainable cli, but is slightly
different to the commonly-used ``argparse``. There are two main
practical differences for the user:

1. The modular structure: the overarching command is ``tp``, running
   this lists the sub commands, likewise with each of those (e.g. ``tp
   plot``). The full command currently has three parts, e.g. ``tp plot
   waterfall``, after which you add the arguments. ``--help`` grants
   information at every step.
2. Multiple arguments must use the ``--thing`` every time, i.e. in
   ``argparse`` you may type ``-d a b c``, but in ``click`` you must
   type ``-d a -d b -d c``.

.. _click: https://click.palletsprojects.com/en/8.0.x/

----------
``tp gen``
----------

Generates inputs for VASP.

- ``tp gen kpar``: generates suggestions for KPAR in VASP. Which you
  choose should depend on factors such as computer architecture.
- ``tp gen kpoints``: Generates a VASP KPOINTS file with a mix of
  weighted and zero-weighted k-points, which are useful for generating
  a dense but cheap mesh for AMSET, amoung other things. Be sure to do
  proper convergence testing!

----------
``tp get``
----------

Gets specific data points from a data file. Rounds dependent variables
to the nearest datum in the file and doesn't interpolate, but will tell
you the exact conditions it used (as opposed to the ones you inputted).
Based around the ``tp.data.resolve.resolve`` function.

- ``tp get amset``
- ``tp get boltztrap``
- ``tp get phonopy``
- ``tp get zt``

-----------
``tp plot``
-----------

Plotting tools.

- ``tp plot avg-rates``: amset scattering rates against doping
   concentration and/ or temperature.
- ``tp plot cumkappa``: cumulative lattice thermal conductivity against
  frequency or mean free path. Python: ``tp.plot.frequency.add_cumkappa``
  or ``tp.plot.mfp.add_cumkappa``.
- ``tp plot dos``: phonon density of states. Python:
  ``tp.plot.frequency.add_dos``.
- ``tp plot kappa``: line graphs of thermal conductivity against
  temperature. Can plot a lines for different directions, or different
  files. Auto-generates legend.
- ``tp plot kappa-target``: heatmap of the highest lattice thermal
  conductivity possible to reach a target ZT against temperature and
  carrier concentration. Useful to help decide if a material is worth
  doing expensive third-order phonon calculations on. Python:
  ``tp.plot.heatmap.add_kappa_target``.
- ``tp plot phonons``: plots and overlays phonon dispersions. Useful for
  convergence. Python: ``tp.plot.phonons.add_dispersion`` or
  ``tp.plot.phonons.add_multi``.
- ``tp plot transport``: up to four electronic transport properties (or
  lattice thermal conductivity) against temperature. Multiple lines can
  represent doping concentration, direction or files. Auto-generates
  legend.
- ``tp plot waterfall``: third-order phonon properties per band and
  qpoint. Can have an optional quantitt projected on the colour axis,
  or be converted to a density map. Python:
  ``tp.plot.frequency.add_waterfall`` or
  ``tp.plot.frequency.add_projected_waterfall`` or
  ``tp.plot.frequency.add_density``.
- ``tp plot wideband``: broadened phonon dispersion, indicating
  scattering. Python: ``tp.plot.phonons.add_wideband``.
- ``tp plot ztmap``: heatmap of zt against temperature and carrier
  concentration. Python: ``tp.plot.heatmap.add_ztmap``.

----------
``tp run``
----------

Runs postprocessing codes efficiently.

- ``tp run boltztrap``: Runs BoltzTraP from VASP outputs and writes to
  hdf5. Significantly faster and more robust than Pytmatgen. Still
  requires the BoltzTraP ``x_trans``. Python: ``tp.data.run.boltztrap``.

-----------
``tp save``
-----------

Saves data to files.

- ``tp save cumkappa``: saves cumulative lattice thermal conductivity to
   a plain text file. Python: ``tp.data.save.cumkappa``.
- ``tp save kappa``: saves lattice thermal conductivity to a plain text
  file.
- ``tp save zt``: saves ZT to a hdf5 file, as well as the max lattice
  thermal conductivity at each temperature and the associated carrier
  concentration to a human-readable yaml, and prints the max ZT across
  all temperatures to the terminal. Python: ``tp.data.save.zt``.
