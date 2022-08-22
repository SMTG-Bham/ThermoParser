.. image:: cumkappa.png
   :alt: Cumulative lattice thermal conductivity against frequency and mean free path.

This shows two plots of cumulative lattice thermal conductivity. The
first one is against frequency, and particularly in conjunction with a
density of states (DoS) or a phonon band structure, this can show the
relative contributions of the constituent atoms to the lattice thermal
conductivity. The second is against mean free path, and can be used to
approximate the effect of nanostructuring, for which reason the ability
to add markers has been included, e.g. here you can see nanostructuring
to 10 nanometers may reduce the lattice thermal conductivity by nearly
50 %.

This also demonstrates use of the ``add_legend`` function supplied with
``tp.axes`` functions, which combines the legends of all plots and
places itself in one of several pre-programmed position, which you can
select with the ``location`` argument. Numbers will place the legend
in one of the axes with ``loc='best'``, and there are also several
descriptive positions, such as ``above`` and ``right``. It also accepts
most ``ax.legend`` arguements such as ``title`` and ``ncol``. If you
want to use custom handles and labels in a multi-axes figure (including
DoS axes), you must specify ``custom=True``.
