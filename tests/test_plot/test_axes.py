"""Tests tp.plot.axes functions

Failures here can have knock on effects on tp.plot modules.
"""

import unittest
from matplotlib import pyplot as plt
from tp.plot import axes

class AxesTestCase(unittest.TestCase):
    def test_fig(self):
        for name, func in axes.__dict__.items():
            if name in axes.__dir__():
                fig, ax = func()
                self.assertIsInstance(fig.get_dpi(), float,
                                      '{} failed.'.format(name))
                plt.close()

    def test_ax(self):
        for name, func in axes.__dict__.items():
            if name in axes.__dir__():
                fig, ax = func()
                while isinstance(ax, list):
                    ax = ax[0]
                self.assertIsInstance(ax.get_xscale(), str,
                                      '{} failed.'.format(name))
                plt.close()

if __name__ == '__main__':
    unittest.main()
