"""Provides pre-sized axes.

The name convention is tp.axes.[size].[number](_[description]), where
decription starts with a shape if the number is greater than one, and
lists the additional things (dos, colourbar) from left to right.
The small style is based somewhat on the Nature guidelines (we can all
dream!), where the axes.small.one figure has a width of 8.6 cm. All
axes are square, and subsequent figures maintain the height of axes.one
per plot (single height plots have a height of 8.3 cm whilst double
height plots have a height of 16.6 cm, which more effectively maintains
relative scale, as colourbars etc. are normally added on the side). As
per the guidelines, font size is 8 pt. The large style uses larger axes
and relatively even larger fonts, which is better suited for
presentations, posters and non-columnated text. legend is for legend
handling, which is primarily used internally but may be useful if you
make your own axes from scratch.

Modules
-------

    small
    large
    legend
"""

from . import small, large, legend
