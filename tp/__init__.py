"""A package for postprocessings transport calculations."""

#Modules
#-------
#
#    calculate
#        tools for calculating derived properties and interpolating
#    settings
#        default values and metadata
#"""

__version__ = '3.3.0' # Remember to update the main README header too!

def docstring_replace(**kwargs):
    def d(f):
        f.__doc__ = f.__doc__.format(**kwargs)
        return f
    return d

from . import settings # required for other imports, prevents circularity
from . import axes, calculate, data, plot, setup, cli
