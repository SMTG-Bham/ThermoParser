"""Tests the tp.plot.mfp.add_cum_kappa function

Failures in tp.plot.axes and tp.calculate can have knock-on effects here.
"""

import unittest
from tp.plot import mfp
from tp.plot.axes import one_small_legend
import matplotlib.pyplot as plt

class BaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'mean_free_path': [[[[[0, 1, 2, 3]]]]],
                 'mode_kappa':     [[[[[0, 2, 11, 100]]]]],
                 'temperature':    [11],
                 'meta':           {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ax = mfp.add_cum_kappa(cls.ax, cls.d, direction='x',
                                   main=True, scale=False)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (1, 3))

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), (0, 113))

class KMinTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'mean_free_path': [[[[[0, 1, 2, 3]]]]],
                 'mode_kappa':     [[[[[0, 2, 11, 100]]]]],
                 'temperature':    [11],
                 'meta':           {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ax = mfp.add_cum_kappa(cls.ax, cls.d, direction='x', kmin=10,
                                   main=True, scale=False)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (2, 3))

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), (0, 113))

class ScaleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'mean_free_path': [[[[[0, 2, 4]]]]],
                 'mode_kappa':     [[[[[0, 1, 2]]]]],
                 'temperature':    [11],
                 'meta':           {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ylim = (1, 2)
        cls.ax.set_ylim(cls.ylim)
        cls.ax = mfp.add_cum_kappa(cls.ax, cls.d, direction='x',
                                         main=False, scale=True)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), self.ylim)

if __name__ == '__main__':
    unittest.main()
