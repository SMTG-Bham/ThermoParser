----------------------------------
Tutorial-06: Package Customisation
----------------------------------

First we'll look at the ``tprc.yaml``, which can customise how
ThermoPlotter behaves, and then how you might expand it.

tprc.yaml
---------

``tprc.yaml`` is a configuration file which allows you to customise
many aspects of ThermoPlotter. A template is provided in the top
directory, but to use it it must be saved as ``~/.config/tprc.yaml``.
The file is basically a dictionary. If you prefer your own style sheet,
you can set it as default by setting ``style: yourstyle``, and the same
for ``large_style``, which will need a larger font size. To change the
number of ticks on an axis, you can change the ``locator``, e.g. to get
a maximum of five major ticks, with 10 minor ticks each you would do:

.. code-block:: yaml

    locator:
        major: 5
        minor: 10

The quantity aliases, conversions, units  and kwargs are all nested
dictionaries like the locators. The aliases are listed as:

.. code-block:: yaml

    to_xxx:
        alias: xxx_name

If you add unit conversions, remember to update the units and labels
too! An example has been provided for converting from
S m<super>-1</super> to S cm<super>-1</super>. There are six labels
dictionaries. ``long_``, ``medium_`` and ``short_labels`` contain the
actual labels of those lengths, while ``labels``, ``inverted_labels``
and ``large_labels`` point to which of those you would like to default
to for small, inverted and large axes, respectively. For example, if
you prefer to put a inverted DoSs in a separate axes rather than a DoS
axes, you may want to set ``inverted_labels: long``. Finally, there are
the default ``kwargs`` for each function, which are passed to the
matplotlib plotting function, such as ``plt.plot``. These override
defaults set in the ThermoPlotter plotting function, but are overridden
by arguments specified by the user.

Extending ThermoPlotter
-----------------------

Due to its modular nature, extending ThermoPlotter should be relatively
easy. So long as it reads or writes data in the standard ThermoPlotter
format, a function can be added to a module without having to worry
about conflicts with the rest of the code. The exception are the
``settings.py`` and ``tprc.yaml``, which may need extra labels, kwargs,
etc. if a new function is added elsewhere. We welcome contributions!
