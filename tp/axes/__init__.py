"""Provides pre-sized axes.

The default style is based somewhat on the Nature guidelines (we can all
dream!), where the axes.one figure has a width of 8.6 cm. All axes are 
square, and subsequent figures maintain the height of axes.one per plot
(single height plots have a height of 8.3 cm whilst double height plots
have a height of 16.6 cm, which more effectively maintains relative
scale, as legends etc. are normally added on the side). As per the
guidelines, font size is 8 pt.
The large style uses larger axes and relatively even larger fonts, which
is better suited for presentations, posters and non-columnated text.
legend is for legend handling.

Modules
-------

    one
    two
    three
    four

    one_large
    two_large
    three_large
    four_large

    legend
"""

from . import one, two, three, four, one_large, two_large, three_large, four_large, legend
