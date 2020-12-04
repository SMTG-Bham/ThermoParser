Command-line scripts for rapid plotting.

--------------
``tp-phonons``
--------------

Plots a phonon dispersion, and optionally a phonon dos. Wraps around
``tp.plot.phonons.add_dispersion`` and ``tp.plot.frequency.add_dos``.

-----------------------
``tp-converge-phonons``
-----------------------

Plots multiple phonon dispersions, good for convergence. Wraps around
``tp.plot.phonons.add_multi``.

---------------
``tp-wideband``
---------------

Plots a phonon dispersion with bands widened due to scattering. Wraps
around ``tp.plot.phonons.add_wideband``.

----------------
``tp-waterfall``
----------------

Plots thermal transport properties against each other, optionally with a
third projected on the colour axis. Wraps around
``tp.plot.frequency.add_waterfall`` and
``tp.plot.frequency.add_projected_waterfall``.

---------------
``tp-cumkappa``
---------------

Plots cumulative lattice thermal conductivity against frequency or mean
free path. Wraps around ``tp.plot.frequency.add_cum_kappa`` and
``tp.plot.mfp.add_cum_kappa``.

------------
``tp-ztmap``
------------

Plots a heatmap of ZT against temperature and carrier concentration.
Wraps around ``tp.plot.heatmap.add_ztmap``.

-------------------
``tp-kappa-target``
-------------------

Plots a heatmap of lattice thermal conductivity required to reach a
minimum ZT against temperature and carrier concentration. Wraps around
``tp.plot.heatmap.add_kappa_target``.
