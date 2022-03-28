The simplest way to use ThermoPlotter is via the command-line interface
(CLI). This is easily navigable due to its modular structure. To get
started, simply type ``tp`` into the command line, and a list of
modules will be returned along with short descriptions. Pick one and
type the full thing in, e.g. ``tp plot``, to get a similar list of the
available functions. Then type in all three parts and ``--help`` to get
a detailed description of how to use the function, e.g.
``tp plot avg-rates --help``. Note ``-h`` is *not* a substitute for
``--help``.

Typically, you will need one or more data files, which are inputted
directly after the command, e.g. ``tp plot avg-rates mesh.h5``. Then
there is a list of optional inputs for customisation. Options with a
slash, e.g. ``--total / --nototal`` are toggles, which can be appended
by themselves, while other options normally require something else
written afterwards, e.g. ``-n 1e19``. If there is a list in square
brackets, the appended text must be from that list. Many options can
take multiple values, in which case the tag will need to be repeated:
``-m '*' -m 'o'``.

Later examples will mention if the CLI can be used instead, but will
focus on the python interface.